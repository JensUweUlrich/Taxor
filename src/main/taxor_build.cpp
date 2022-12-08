#include <lemon/list_graph.h> /// Must be first include.

#include <iostream>

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

#include <syncmer.hpp>

#include <parse_ncbi_taxonomy.hpp>

#include "taxor_build.hpp"
#include "index.hpp"
#include "store_index.hpp"
#include "taxor_build_configuration.hpp"

namespace taxor::build
{

void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::build::configuration & config)
{
    parser.info.version = "0.1.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Creates and HIXF index of a given set of fasta files";

    parser.info.description.emplace_back("Creates an HIXF index using either kmers, syncmers or FracMinHash sketches");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.input_file_name,
                      '\0', "input-file", "",
 /*                     "Provide the prefix you used for the output prefix in the chopper count --output-prefix option. "
                      "If you have different means of estimating the k-mer counts of your input data, make sure that a "
                      "file [INPUT-PREFIX].count exists. It needs to be tab-separated and consist of two columns: "
                      "\"[filepath] [tab] [weight/count]\".",*/
                      seqan3::option_spec::required);
    parser.add_list_item("", "Example file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "GCF_000839185.1 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/185/GCF_000839185.1_ViralProj14174 Cowpox virus    k__Viruses;p__Nucleocytoviricota;c__Pokkesviricetes;o__Chitovirales;f__Poxviridae;g__Orthopoxvirus;s__Cowpox virus");
    parser.add_list_item("", "GCF_000860085.1 https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/860/085/GCF_000860085.1_ViralProj15241 Vaccinia virus  k__Viruses;p__Nucleocytoviricota;c__Pokkesviricetes;o__Chitovirales;f__Poxviridae;g__Orthopoxvirus;s__Vaccinia virus");
    parser.add_list_item("", "```");

    parser.add_option(config.input_sequence_folder, '\0', "input-sequence_dir", "directory containing the fasta reference files");

    parser.add_option(config.output_file_name, '\0', "output-filename", "A file name for the resulting index.");

    parser.add_option(config.kmer_size, '\0', "kmer-size", "size of kmers used for index construction",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

    parser.add_option(config.syncmer_size, '\0', "syncmer-size", "size of syncmer used for index construction",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(30)});

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use. Currently, only merging of sketches is parallelized, so if option "
                      "--rearrange-user-bins is not set, --threads will have no effect.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

    parser.add_flag(config.use_syncmer,'\0', "use-syncmer", "enable using syncmers for smaller index size");

    parser.add_flag(config.output_verbose_statistics,
                    '\0', "output-verbose-statistics",
                    "Enable verbose statistics to be "
                    "printed to std::cout. If the flag --determine-best-tmax is not set, this flag is ignored "
                    "and has no effect.",
                    seqan3::option_spec::hidden);

    parser.add_flag(config.debug,
                    '\0', "debug",
                    "Enables debug output in layout file.",
                    seqan3::option_spec::hidden);
}
/*
void sanity_checks(layout::data_store const & data, chopper::layout::configuration & config)
{
    if (config.rearrange_user_bins)
        config.estimate_union = true;

    if (config.estimate_union &&
        (!std::filesystem::exists(config.sketch_directory) || std::filesystem::is_empty(config.sketch_directory)))
    {
        throw seqan3::argument_parser_error{"The directory " + config.sketch_directory.string() + " must be present "
                                            "and not empty in order to enable --estimate-union or "
                                            "--rearrange-user-bins (created with chopper count)."};
    }

    if (data.filenames.empty())
        throw seqan3::argument_parser_error{seqan3::detail::to_string("The file ", config.count_filename.string(),
                                                                      " appears to be empty.")};

    if (config.aggregate_by_column != -1 && data.extra_information[0].empty())
    {
        throw seqan3::argument_parser_error{"Aggregate Error: You want to aggregate by something but your "
                                            "file does not contain any extra information columns."};
    }

    // note that config.aggregate_by_column cannot be 0 or 1 because the parser check the valid range [2, max]
    assert(config.aggregate_by_column == -1 || config.aggregate_by_column > 1);
    if ((config.aggregate_by_column - 2/*extrainfo starts at 2*/  /* ) > static_cast<int>(data.extra_information[0].size()))
    {
        throw seqan3::argument_parser_error{"Aggregate Error: You want to aggregate by a column index that "
                                            "is larger than the number of extra information columns."};
    }
}

*/


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

inline auto create_filename_clusters(taxor::build::configuration const taxor_config, 
                                     std::vector<taxonomy::Species> &orgs,
                                     std::map<std::string, uint64_t> &user_bin_map)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters;

    for (uint64_t org_index = 0; org_index < orgs.size(); ++org_index) //auto& species : orgs)
    {
        std::string filepath = taxor_config.input_sequence_folder + "/" + orgs[org_index].file_stem + "_genomic.fna.gz";
        filename_clusters[orgs.at(org_index).accession_id].push_back(filepath);
        user_bin_map.emplace(std::make_pair(filepath, org_index));
    }

    return filename_clusters;
}

inline void count_syncmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                           chopper::count::configuration const & count_config,
                           taxor::build::configuration const taxor_config)
{
    using traits_type = seqan3::sequence_file_input_default_traits_dna;
    using sequence_file_t = seqan3::sequence_file_input<traits_type, seqan3::fields<seqan3::field::seq>>;

    uint8_t t_syncmer = ceil((taxor_config.kmer_size - taxor_config.syncmer_size) / 2) + 1;

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
			    std::vector<uint64_t> syncmer_hashes = hashing::seq_to_syncmers(taxor_config.kmer_size, seq, taxor_config.syncmer_size, t_syncmer);
                for (auto &hash : syncmer_hashes)
                    sketch.add(hash);
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
    auto filename_clusters = create_filename_clusters(taxor_config, orgs, user_bin_map);
	//auto filename_clusters = chopper::count::read_data_file(config);

	chopper::detail::apply_prefix(count_config.output_prefix, count_config.count_filename, count_config.sketch_directory);
    chopper::count::check_filenames(filename_clusters, count_config);
    if (taxor_config.use_syncmer)
        count_syncmers(filename_clusters, count_config, taxor_config);
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
    args.syncmer_size = config.syncmer_size;
    args.threads = config.threads;
	args.compute_syncmer = config.use_syncmer;
    if (config.use_syncmer)
        args.t_syncmer = ceil((args.kmer_size - args.syncmer_size) / 2) + 1;
	
    hixf::build_data data{};
   
    hixf::create_ixfs_from_chopper_pack(data, args);
    
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hixf.user_bins.num_user_bins(); ++i)
    {
        bin_path.push_back(std::vector<std::string>{data.hixf.user_bins.filename_of_user_bin(i)});
        orgs.at(user_bin_map[data.hixf.user_bins.filename_of_user_bin(i)]).user_bin = i;
    }

   taxor_index<hixf::hierarchical_interleaved_xor_filter<uint8_t>> index{hixf::window{args.window_size},
                                                                                args.shape,
                                                                                args.kmer_size,
                                                                                args.syncmer_size,
                                                                                args.t_syncmer,
                                                                                args.parts,
                                                                                args.compute_syncmer,
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

        // TODO: sanity check of parameters

    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[TAXOR BUILD ERROR] " << ext.what() << '\n';
        return -1;
    }

    std::vector<taxonomy::Species> orgs = taxonomy::parse_refseq_taxonomy_file(config.input_file_name);
    // map filename to index of species in orgs vector 
    std::map<std::string, uint64_t> user_bin_map{};
    create_layout(config, orgs, user_bin_map);

    build_hixf(config, orgs, user_bin_map);


    return 0;
}

} // namespace chopper::layout
