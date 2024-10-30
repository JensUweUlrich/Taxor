// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#include <lemon/list_graph.h> /// Must be first include.

#include <iostream>
#include <filesystem>
#include <regex>
#include <type_traits>

#include <chopper/count/read_data_file.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/count/check_filenames.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/count/output.hpp>
#include <chopper/layout/configuration.hpp>
#include <chopper/layout/data_store.hpp>
#include <chopper/layout/filenames_data_input.hpp>
#include <chopper/layout/hierarchical_binning.hpp>
#include <chopper/layout/output.hpp>

#include <build/create_ixfs_from_chopper_pack.hpp>
#include <build/build_data.hpp>
#include <build/strong_types.hpp>
#include <build/build_arguments.hpp>
#include <build/dna4_traits.hpp>
#include <build/adjust_seed.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <syncmer.hpp>

#include <parse_ncbi_taxonomy.hpp>

#include "taxor_build.hpp"
#include "index.hpp"
#include "store_index.hpp"
#include "taxor_build_configuration.hpp"

namespace taxor::build
{

using sequence_file_t = seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::seq>>;

void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::build::configuration & config)
{
    parser.info.version = "0.2.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Creates an HIXF index of a given set of fasta files";

    parser.info.description.emplace_back("Creates an HIXF index using either k-mers or syncmers");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.input_file_name,
                      '\0', "input-file", "tab-separated-value file containing taxonomy information and reference file names",
                      seqan3::option_spec::required);

    parser.add_option(config.input_sequence_folder, '\0', "input-sequence-dir", "directory containing the fasta reference files");

    parser.add_option(config.output_file_name, '\0', "output-filename", "A file name for the resulting index.");

    parser.add_option(config.kmer_size, '\0', "kmer-size", "size of kmers used for index construction",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(64)});

    parser.add_option(config.syncmer_size, '\0', "syncmer-size", "size of syncmer used for index construction",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(26)});

    parser.add_option(config.window_size, '\0', "window-size", "window size of minimizer scheme used for index construction",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(96)});

    parser.add_option(config.scaling, '\0', "scaling", "factor for scaling down syncmer/minimizer sketches",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(10), static_cast<size_t>(1000)});

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

    parser.add_flag(config.use_syncmer,'\0', "use-syncmer", "enable using syncmers for smaller index size");

    parser.add_flag(config.output_verbose_statistics,
                    '\0', "output-verbose-statistics",
                    "Enable verbose statistics to be printed",
                    seqan3::option_spec::hidden);

    parser.add_flag(config.debug,
                    '\0', "debug",
                    "Enables debug output in layout file.",
                    seqan3::option_spec::hidden);
}

std::vector<std::string> str_split(std::string &str, char delimiter)
{

    std::stringstream str_stream(str);
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(str_stream, segment, delimiter))
    {
        seglist.push_back(segment);
    }

    return std::move(seglist);
}

void sanity_checks(taxor::build::configuration & config)
{
    if (config.use_syncmer)
    {
        if (config.kmer_size > 30)
        {
            throw seqan3::argument_parser_error{"The chosen k-mer size is too large for the syncmer scheme. Please choose a k-mer size <= 30 or use the minimizer scheme"};
        }
    }

    // enable using seveal input file
    config.input_files = str_split(config.input_file_name, ',');

    // check whether all input files exist
    for (std::string & f : config.input_files)
    {   
        std::filesystem::path filepath{f};
        if (!std::filesystem::exists(filepath))
            throw seqan3::argument_parser_error{"Please check the given input file(s). \nThe following input file does not exist: " + f};
    }

    // check whether taxonomy files have correct input format
    for (std::string & f : config.input_files)
    {
        try
        {
            std::vector<taxonomy::Species> orgs = taxonomy::parse_refseq_taxonomy_file(f);
        }
        catch (std::out_of_range const &e)
        {
            throw seqan3::argument_parser_error{"Error parsing the taxonomy file: " + f};    
        }
    }

    // enable using several input directories
    config.input_folders = str_split(config.input_sequence_folder, ',');
    
    // check whether all input files exist
    for (std::string & f : config.input_folders)
    {   
        std::filesystem::path filepath{f};
        if (!std::filesystem::exists(filepath))
            throw seqan3::argument_parser_error{"Please check the given input folder(s). \nThe following input folder does not exist: " + f};
    }


}

size_t determine_best_number_of_technical_bins(chopper::layout::data_store & data, chopper::layout::configuration & config)
{
    std::stringstream * const output_buffer_original = data.output_buffer;
    std::stringstream * const header_buffer_original = data.header_buffer;

    std::set<size_t> potential_t_max = [&] ()
    {
        std::set<size_t> result;

        for (size_t t_max = 64; t_max <= config.tmax; t_max *= 2)
            result.insert(t_max);

        // Additionally, add the t_max that is closest to the sqrt() of the number of
        // user bins, as it is expected to evenly spread bins and may perform well.
        size_t const user_bin_count{std::ranges::size(data.kmer_counts)};
        size_t const sqrt_t_max{chopper::next_multiple_of_64(std::ceil(std::sqrt(user_bin_count)))};
        result.insert(sqrt_t_max);

        return result;
    }();

    // with -determine-best-tmax the algorithm is executed multiple times and result with the minimum
    // expected query costs are written to the standard output

    double best_expected_HIXF_query_cost{std::numeric_limits<double>::infinity()};
    size_t best_t_max{};
    size_t max_hixf_id{};
    size_t t_max_64_memory{};

    for (size_t const t_max : potential_t_max)
    {
        std::stringstream output_buffer_tmp;
        std::stringstream header_buffer_tmp;
        config.tmax = t_max;                               // overwrite tmax
        data.output_buffer = &output_buffer_tmp;           // overwrite buffer
        data.header_buffer = &header_buffer_tmp;           // overwrite buffer
        data.previous = chopper::layout::previous_level{}; // reset previous IBF, s.t. data refers to top level IBF

        chopper::layout::hibf_statistics global_stats{config, data.fp_correction, data.kmer_counts};
        data.stats = &global_stats.top_level_ibf;

        // execute the actual algorithm
        size_t const max_hixf_id_tmp = chopper::layout::hierarchical_binning{data, config}.execute();
        global_stats.finalize();

        global_stats.print_summary(t_max_64_memory, config.output_verbose_statistics);

        // Use result if better than previous one.
        if (global_stats.expected_HIBF_query_cost < best_expected_HIXF_query_cost)
        {
            *output_buffer_original = std::move(output_buffer_tmp);
            *header_buffer_original = std::move(header_buffer_tmp);
            max_hixf_id = max_hixf_id_tmp;
            best_t_max = t_max;
            best_expected_HIXF_query_cost = global_stats.expected_HIBF_query_cost;
        }
        else if (!config.force_all_binnings)
        {
            break;
        }
    }

    std::cout << "# Best t_max (regarding expected query runtime): " << best_t_max << '\n';
    config.tmax = best_t_max;
    return max_hixf_id;
}

template < bool RECURSIVE > std::vector<std::filesystem::path> file_list( std::filesystem::path dir, std::regex pattern )
{
    std::vector<std::filesystem::path> result ;

    using iterator = std::conditional< RECURSIVE, 
                                       std::filesystem::recursive_directory_iterator, std::filesystem::directory_iterator >::type ;

    const iterator end ;
    for( iterator iter { dir } ; iter != end ; ++iter )
    {
        const std::string fname = iter->path().filename().string() ;
        if( std::filesystem::is_regular_file(*iter) && std::regex_match( fname, pattern ) ) result.push_back( *iter ) ;
    }
    
    return result ;
}

inline auto create_filename_clusters(taxor::build::configuration const taxor_config, 
                                     std::vector<taxonomy::Species> &orgs,
                                     std::map<std::string, uint64_t> &user_bin_map)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters;

    for (uint64_t org_index = 0; org_index < orgs.size(); ++org_index) //auto& species : orgs)
    {
        std::string reg = "^" + orgs[org_index].file_stem + "[\\_a-z]*\\.[a-z\\.]+";
        std::regex reg1(reg);
        bool found = false;
        std::string filepath{};
        for (std::string folder : taxor_config.input_folders)
        {

            std::filesystem::path fpath{folder};
            std::vector<std::filesystem::path> files = file_list<false>(fpath, reg1);
            
            if (files.size() == 1)
            {
                filepath = files[0].string();
                found = true;
                break;
            }

            if (files.size() > 1)
                throw std::logic_error{"More than one file was found for " + orgs[org_index].accession_id + " in " + folder};
        }

        if (!found)
            throw std::logic_error{"Could not find a genome file for " + orgs[org_index].accession_id};

        filename_clusters[orgs.at(org_index).accession_id].push_back(filepath);
        user_bin_map.emplace(std::make_pair(filepath, org_index));
    }

    return filename_clusters;
}

inline void count_minimizers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                           chopper::count::configuration const & count_config,
                           taxor::build::configuration const taxor_config)
{
    using traits_type = seqan3::sequence_file_input_default_traits_dna;
    using sequence_file_t = seqan3::sequence_file_input<traits_type, seqan3::fields<seqan3::field::seq>>;

    std::ofstream fout{count_config.count_filename};

    if (!fout.good())
        throw std::runtime_error{"Could not open file" + count_config.count_filename.string() + " for reading."};

    // create the hll dir if it doesn't already exist
    if (!count_config.disable_sketch_output)
        std::filesystem::create_directory(count_config.sketch_directory);

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    seqan3::shape const shape = seqan3::ungapped{taxor_config.kmer_size};
    auto minimizer_view = seqan3::views::minimiser_hash(shape,
                                                        seqan3::window_size{taxor_config.window_size},
                                                        seqan3::seed{hixf::adjust_seed(shape.count())});

    #pragma omp parallel for schedule(static) num_threads(count_config.threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        chopper::sketch::hyperloglog sketch(count_config.sketch_bits);

        std::vector<seqan3::dna5_vector> refs{};
        for (auto const & filename : cluster_vector[i].second)
        {
            for (auto && [seq] : sequence_file_t{filename})
            {
			    for (auto hash : seq | minimizer_view)
                {
                    if (taxor_config.scaling > 1)
                    {
                        uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                        if (double(v) <= double(UINT64_MAX) / double(taxor_config.scaling))
                        {
                            sketch.add(hash);
                        }
                    }
                    else
                    {
                        sketch.add(hash);
                    }
                }
            }
		}
        //process_sequence_files(cluster_vector[i].second, config, sketch);

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t const weight = sketch.estimate();

        #pragma omp critical
        chopper::count::write_count_file_line(cluster_vector[i], weight, fout);

        if (!count_config.disable_sketch_output)
            chopper::count::write_sketch_file(cluster_vector[i], sketch, count_config);
    }
}



inline void count_syncmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                           chopper::count::configuration const & count_config,
                           taxor::build::configuration const taxor_config)
{
    using traits_type = seqan3::sequence_file_input_default_traits_dna;
    using sequence_file_t = seqan3::sequence_file_input<traits_type, seqan3::fields<seqan3::field::seq>>;

    size_t t_syncmer = ceil((taxor_config.kmer_size - taxor_config.syncmer_size + 1) / 2);

    std::ofstream fout{count_config.count_filename};

    if (!fout.good())
        throw std::runtime_error{"Could not open file" + count_config.count_filename.string() + " for reading."};

    // create the hll dir if it doesn't already exist
    if (!count_config.disable_sketch_output)
        std::filesystem::create_directory(count_config.sketch_directory);

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector{};
    for (auto const & cluster : filename_clusters)
        cluster_vector.emplace_back(cluster.first, cluster.second);

    #pragma omp parallel for schedule(static) num_threads(count_config.threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        chopper::sketch::hyperloglog sketch(count_config.sketch_bits);

        std::vector<seqan3::dna5_vector> refs{};
        for (auto const & filename : cluster_vector[i].second)
        {
            for (auto && [seq] : sequence_file_t{filename})
            {
			    ankerl::unordered_dense::set<size_t> syncmer_hashes = hashing::seq_to_syncmers(taxor_config.kmer_size, seq, taxor_config.syncmer_size, t_syncmer);   
                for (auto &hash : syncmer_hashes)
                {
                    if (taxor_config.scaling > 1)
                    {
                        uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                        if (double(v) <= double(UINT64_MAX) / double(taxor_config.scaling))
                        {
                            sketch.add(hash);
                        }
                    }
                    else
                    {
                        sketch.add(hash);
                    }
                }
            }
		}
        //process_sequence_files(cluster_vector[i].second, config, sketch);

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t const weight = sketch.estimate();

        #pragma omp critical
        chopper::count::write_count_file_line(cluster_vector[i], weight, fout);

        if (!count_config.disable_sketch_output)
            chopper::count::write_sketch_file(cluster_vector[i], sketch, count_config);
    }
}


void create_layout(taxor::build::configuration const taxor_config, 
                   std::vector<taxonomy::Species> &orgs,
                   std::map<std::string, uint64_t> &user_bin_map)
{
	chopper::count::configuration count_config{};
	//config.data_file = taxor_config.input_file_name;
	count_config.k = taxor_config.kmer_size;
	count_config.threads = taxor_config.threads;
	count_config.output_prefix = "chopper";
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters {};
    try
    {
        filename_clusters = create_filename_clusters(taxor_config, orgs, user_bin_map);
        chopper::detail::apply_prefix(count_config.output_prefix, count_config.count_filename, count_config.sketch_directory);
        chopper::count::check_filenames(filename_clusters, count_config);
    }
    catch(const std::exception& e)
    {
        throw;
    }
    
    
    return;
	
    if (taxor_config.use_syncmer)
    {
        count_syncmers(filename_clusters, count_config, taxor_config);
    }
    else if (taxor_config.window_size > taxor_config.kmer_size)
    {
        count_minimizers(filename_clusters, count_config, taxor_config);
    }
    else
        chopper::count::count_kmers(filename_clusters, count_config);

	chopper::layout::configuration layout_config;
    chopper::layout::data_store data;

	layout_config.input_prefix = count_config.output_prefix;
	
	chopper::detail::apply_prefix(layout_config.input_prefix, layout_config.count_filename, layout_config.sketch_directory);
    // Read in the data file containing file paths, kmer counts and additional information.
    chopper::layout::read_filename_data_file(data, layout_config); 

    layout_config.rearrange_user_bins = true;
	layout_config.determine_best_tmax = true;
	layout_config.estimate_union = true;
	layout_config.tmax = 4096;
	layout_config.threads = taxor_config.threads;

	std::stringstream output_buffer;
    std::stringstream header_buffer;

    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;

    size_t max_hixf_id;
    max_hixf_id = determine_best_number_of_technical_bins(data, layout_config);
    
    std::cout << "write Layout header" << std::endl;

    // brief Write the output to the layout file.
    std::ofstream fout{layout_config.output_filename};
    chopper::layout::write_layout_header_to(layout_config, max_hixf_id, header_buffer.str(), fout);
    fout << output_buffer.str();

}

void build_hixf(taxor::build::configuration const config, 
                std::vector<taxonomy::Species> &orgs,
                std::map<std::string, uint64_t> &user_bin_map)
{
  
	hixf::build_arguments args{};
	args.bin_file = std::filesystem::path{"binning.out"};
	args.out_path = config.output_file_name;
    args.kmer_size = config.kmer_size;
    args.window_size = config.window_size;
    args.shape = seqan3::shape{seqan3::ungapped{args.kmer_size}};
    args.syncmer_size = config.syncmer_size;
    args.threads = config.threads;
	args.compute_syncmer = config.use_syncmer;
    args.scaling = config.scaling;
    if (config.use_syncmer)
        args.t_syncmer = ceil((args.kmer_size - args.syncmer_size + 1) / 2);
	
    hixf::build_data data{};
   
    hixf::create_ixfs_from_chopper_pack(data, args);
    
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hixf.user_bins.num_user_bins(); ++i)
    {
        bin_path.push_back(std::vector<std::string>{data.hixf.user_bins.filename_of_user_bin(i)});
        uint64_t org_index = user_bin_map[data.hixf.user_bins.filename_of_user_bin(i)];
        orgs.at(org_index).user_bin = i;
        orgs.at(org_index).seq_len = 0;
        for (auto && [seq] : sequence_file_t{data.hixf.user_bins.filename_of_user_bin(i)})
	    {
            orgs.at(org_index).seq_len += seq.size();
        }
    }
    
    taxor_index<hixf::hierarchical_interleaved_xor_filter<uint8_t>> index{hixf::window{args.window_size},
                                                                                args.shape,
                                                                                args.kmer_size,
                                                                                args.syncmer_size,
                                                                                args.t_syncmer,
                                                                                args.parts,
                                                                                args.compute_syncmer,
                                                                                args.scaling,
                                                                                args.compressed,
                                                                                bin_path,
                                                                                orgs,
                                                                                std::move(data.hixf)};
    
    store_index(args.out_path, index);
}

int execute(seqan3::argument_parser & parser)
{
    taxor::build::configuration config;
    //chopper::layout::data_store data;

    set_up_subparser_layout(parser, config);

    try
    {
        parser.parse();
        std::cout << "checking input ... " << std::flush;
        sanity_checks(config);
        std::cout << "done!" << std::endl;

    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[TAXOR BUILD ERROR] " << ext.what() << '\n';
        return -1;
    }


    // parse taxonomy input files
    std::cout << "parsing taxonomy input files ... " << std::flush;
    std::vector<taxonomy::Species> orgs {};

    for (std::string &f : config.input_files)
    {
        std::vector<taxonomy::Species> org_tmp = taxonomy::parse_refseq_taxonomy_file(f);
        orgs.insert(orgs.end(), org_tmp.begin(), org_tmp.end());
    }
    std::cout << "done!" << std::endl;

    std::cout << "creating HIXF layout ... " << std::flush;
    // map filename to index of species in orgs vector 
    std::map<std::string, uint64_t> user_bin_map{};
    try
    {
        create_layout(config, orgs, user_bin_map);
    }
    catch (std::exception const &e)
    {
        std::cerr << "[TAXOR BUILD ERROR] " << e.what() << '\n';
        return -1;
    }
    
    std::cout << "done!" << std::endl;
    std::cout << "building HIXF index ... " << std::flush;
    build_hixf(config, orgs, user_bin_map);
    std::cout << "done!" << std::endl;

    return 0;
}

} // namespace chopper::layout
