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
#include <cmath>

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

/*!\file taxor_build.cpp
 * \brief Implements the `taxor build` subcommand: parses a taxonomy TSV plus one or
 * more reference genome directories (searched recursively), computes a hierarchical
 * binning layout for the reference genomes via the third-party `chopper` library, and
 * builds the resulting HIXF (Hierarchical Interleaved XOR Filter) index. The actual
 * filter construction from the chopper layout/pack files is delegated to hixf/build/
 * (see hixf::create_ixfs_from_chopper_pack()); this file is responsible for CLI
 * option handling, input validation, taxonomy/genome-file matching, and orchestrating
 * the counting -> layout -> index-build pipeline.
 */

namespace taxor::build
{

using sequence_file_t = seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::seq>>;

/*!\brief Register all `taxor build` CLI options/flags on the sub-parser and bind
 * them to the corresponding fields of \p config.
 * \param parser The `build` sub-parser to configure.
 * \param config The configuration struct whose fields are bound as option targets;
 * populated once `parser.parse()` is called by the caller.
 */
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

/*!\brief Split \p str into substrings on every occurrence of \p delimiter.
 * \param str The string to split (not modified).
 * \param delimiter The single character to split on.
 * \return The substrings between delimiters, in order; empty segments (e.g. from
 * consecutive delimiters) are included as empty strings.
 */
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

/*!\brief Validate and finalize the raw CLI configuration before the build pipeline runs.
 *
 * Rejects a syncmer k-mer size above 30, splits the comma-separated
 * `input_file_name`/`input_sequence_folder` strings into config.input_files /
 * config.input_folders, checks that every resulting path exists, and eagerly
 * parses each taxonomy file to fail fast on a malformed input file rather than
 * partway through the (much more expensive) layout/build steps.
 * \param config The configuration to validate; config.input_files and
 * config.input_folders are populated as a side effect.
 */
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
        catch (std::exception const &e)
        {
            throw seqan3::argument_parser_error{"Error parsing the taxonomy file: " + f + " (" + e.what() + ")"};
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

/*!\brief Try several candidate technical-bin counts (t_max) for the HIXF layout and
 * keep the one with the lowest expected query cost.
 *
 * Candidates are the powers of two from 64 up to config.tmax, plus the multiple of 64
 * closest to sqrt(number of user bins) (expected to spread bins evenly). For each
 * candidate, chopper::layout::hierarchical_binning is executed from scratch (into a
 * temporary output/header buffer pair so a worse candidate doesn't clobber a better
 * one already found) and its expected HIBF query cost is compared to the best found so
 * far. Candidates are tried in increasing t_max order and the search stops at the first
 * non-improving candidate, unless config.force_all_binnings requests exhaustively
 * trying every candidate. On return, config.tmax is overwritten with the best t_max
 * found, and data.output_buffer / data.header_buffer contain the layout for that t_max.
 * \param data Layout data store; its output/header buffers and `previous` level are
 * mutated across trials and end up holding the winning layout's contents.
 * \param config Layout configuration; config.tmax is overwritten with the best t_max found.
 * \return The maximum HIXF id assigned by the winning hierarchical_binning run.
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

        // chopper::layout::hibf_statistics (a vendored, third-party class - not something we can
        // edit directly, and any local patch to the fetched copy under build/ would be lost on a
        // clean rebuild) estimates memory assuming an Interleaved BLOOM Filter: its private
        // compute_bin_size() uses the standard Bloom filter bit-size formula
        // -n*k / ln(1 - p^(1/k)), driven by config.num_hash_functions (k) and
        // config.false_positive_rate (p). Taxor actually builds a Hierarchical Interleaved XOR
        // Filter, whose size does not depend on a target false positive rate or a number of hash
        // functions at all: an XOR filter needs a fixed ~1.23x overhead over the number of stored
        // elements (see e.g. XorFilter::XorFilter in xorfilter.hpp), at a constant number of
        // fingerprint bits per element (8, since every hierarchical_interleaved_xor_filter Taxor
        // builds uses uint8_t fingerprints - see hixf::hierarchical_interleaved_xor_filter<uint8_t>
        // in index.hpp/taxor_build.cpp). compute_bin_size() is linear in the number of elements
        // stored, so the ratio of XOR-filter bits/element to Bloom-filter bits/element rescales
        // chopper's reported (Bloom filter) byte count into an exact - not approximate - HIXF byte
        // count, without needing to reimplement chopper's private per-level bin-size aggregation.
        double const bloom_bits_per_element =
            -static_cast<double>(config.num_hash_functions) /
            std::log(1.0 - std::exp(std::log(config.false_positive_rate) / static_cast<double>(config.num_hash_functions)));
        double constexpr xor_filter_overhead_factor = 1.23; // standard XOR filter fingerprint-array overhead (arrayLength = 32 + 1.23*size)
        double constexpr xor_fingerprint_bits = 8.0;        // Taxor's HIXF always uses uint8_t fingerprints
        double const hibf_to_hixf_size_factor = (xor_filter_overhead_factor * xor_fingerprint_bits) / bloom_bits_per_element;

        size_t const hixf_size_bytes =
            static_cast<size_t>(std::ceil(static_cast<double>(global_stats.total_hibf_size_in_byte()) * hibf_to_hixf_size_factor));
        std::cout << "# Estimated size for the actual Hierarchical Interleaved XOR Filter (HIXF): "
                  << chopper::byte_size_to_formatted_str(hixf_size_bytes) << '\n';

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

/*!\brief List regular files found in \p input_folders, keyed by an accession id
 * extracted from each file's name.
 *
 * For each regular file, the accession id is reconstructed from the first two
 * underscore-separated parts of the filename stem (e.g. "GCF_000123456.1_..."
 * yields the accession "GCF_000123456.1"); files whose stem has no underscore
 * are silently skipped. Only the `RECURSIVE == true` instantiation is used
 * anywhere in this codebase (via file_list<true>() in create_filename_clusters()),
 * to walk genome_updater-style nested per-species/per-assembly directory layouts;
 * the `RECURSIVE == false` instantiation has no call site.
 * \tparam RECURSIVE If true, each folder is walked recursively
 * (std::filesystem::recursive_directory_iterator); if false, only its immediate
 * entries are listed (std::filesystem::directory_iterator).
 * \param input_folders Directories to search.
 * \return Map from accession id to the matched file's path. If the same accession
 * id occurs in more than one file, only the first one encountered is kept
 * (std::map::emplace does not overwrite an existing key).
 */
template < bool RECURSIVE > std::map<std::string, std::filesystem::path> file_list( std::vector<std::string> input_folders)
{
    std::map<std::string, std::filesystem::path> result ;

    using iterator = std::conditional< RECURSIVE, 
                                       std::filesystem::recursive_directory_iterator, std::filesystem::directory_iterator >::type ;
    for (std::string folder : input_folders)
    {
        std::filesystem::path fpath{folder};
   
        const iterator end ;
        for( iterator iter {folder} ; iter != end ; ++iter )
        {
            const std::string fname = iter->path().filename().string() ;
            if( std::filesystem::is_regular_file(*iter)) 
            {
                std::string filename = (*iter).path().stem().string();
                std::vector<std::string> file_parts = str_split(filename, '_');
                if (file_parts.size() > 1)
                {
                    std::string accession = file_parts[0] + "_" + file_parts[1];
                    result.emplace(std::make_pair(accession, *iter)) ;
                }
            }
        }   
    }
    
    return result ;
}

/*!\brief Build a one-genome-file-per-cluster mapping for every taxon in \p orgs,
 * and record each file's index into \p orgs for later lookup.
 *
 * Every species is placed in its own single-file cluster (keyed by accession id);
 * the actual genome file is looked up by accession id in the map built by
 * file_list<true>(). \p user_bin_map is populated as a side effect so that later
 * pipeline stages (see build_hixf()) can map a resolved HIXF user-bin filename
 * back to the corresponding entry in \p orgs.
 * \param taxor_config Build configuration; only input_folders is used here.
 * \param orgs Parsed taxonomy entries; not modified, but indexed to build \p user_bin_map.
 * \param user_bin_map Output parameter: filled with (genome file path -> index into orgs).
 * \return Map from accession id to a single-element vector containing that
 * genome's file path (the vector form matches chopper's expected cluster shape).
 */
inline auto create_filename_clusters(taxor::build::configuration const taxor_config,
                                     std::vector<taxonomy::Species> &orgs,
                                     std::map<std::string, uint64_t> &user_bin_map)
{
    robin_hood::unordered_map<std::string, std::vector<std::string>> filename_clusters;
    // recursive: --input-sequence-dir may contain the genome files nested in
    // subdirectories rather than directly inside it (e.g. genome_updater's
    // per-species/per-assembly directory layout)
    std::map<std::string, std::filesystem::path> files = file_list<true>(taxor_config.input_folders);
    for (uint64_t org_index = 0; org_index < orgs.size(); ++org_index) //auto& species : orgs)
    {
        // reg/reg1/found below are computed but never used: the actual lookup a few
        // lines down is a direct accession_id map lookup, not a regex match against
        // file_stem. Left in place as-is (not part of this documentation pass' scope).
        std::string reg = "^" + orgs[org_index].file_stem + "[\\_a-z]*\\.[a-z\\.]+";
        std::regex reg1(reg);
        bool found = false;
        std::string filepath{};
        if (files.contains(orgs[org_index].accession_id))
        {
            filepath = files.at(orgs[org_index].accession_id).string();
        }
        else
        {
            throw std::logic_error{"Could not find a genome file for " + orgs[org_index].accession_id};
        }
        filename_clusters[orgs.at(org_index).accession_id].push_back(filepath);
        user_bin_map.emplace(std::make_pair(filepath, org_index));
    }

    return filename_clusters;
}

/*!\brief Compute per-cluster HyperLogLog sketches and k-mer-count estimates using the
 * (open) minimizer scheme, and write them to \p count_config's count/sketch files.
 *
 * Runs one iteration per cluster in an OpenMP parallel-for loop (each iteration builds
 * its own local sketch, so no locking is needed there); only the shared count/sketch
 * file writes are wrapped in an `omp critical` section, since std::ofstream is not
 * safe for concurrent writes from multiple threads. If taxor_config.scaling > 1, hashes
 * are downsampled by re-hashing each minimizer hash with wyhash and only keeping it
 * when the re-hash falls below UINT64_MAX / scaling, approximating a 1/scaling
 * subsample while keeping the same hashes selected across all files that share a hash
 * value (this is what lets scaled sketches from different genomes stay comparable).
 * \param filename_clusters One or more genome files per cluster key.
 * \param count_config chopper counting configuration (output paths, thread count, HLL
 * sketch precision, etc.).
 * \param taxor_config Build configuration; kmer_size, window_size and scaling are used here.
 */
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
                        // downsample: keep hash iff its wyhash falls in the bottom 1/scaling
                        // of the uint64 range. Re-hashing avoids bias from the minimizer
                        // hash's own value distribution and keeps the *same* hashes selected
                        // across every genome, since the threshold test is a pure function of v.
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



/*!\brief Compute per-cluster HyperLogLog sketches and k-mer-count estimates using the
 * syncmer scheme, and write them to \p count_config's count/sketch files.
 *
 * Analogous to count_minimizers(), but hashes sequences with
 * hashing::seq_to_syncmers() (parameterized by the closed-syncmer offset
 * `t_syncmer`, derived from kmer_size/syncmer_size) instead of a minimizer view.
 * The same wyhash-based downsampling by taxor_config.scaling and OpenMP
 * parallel-for-with-critical-section pattern as count_minimizers() applies here.
 * \param filename_clusters One or more genome files per cluster key.
 * \param count_config chopper counting configuration (output paths, thread count, HLL
 * sketch precision, etc.).
 * \param taxor_config Build configuration; kmer_size, syncmer_size and scaling are used here.
 */
inline void count_syncmers(robin_hood::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                           chopper::count::configuration const & count_config,
                           taxor::build::configuration const taxor_config)
{
    using traits_type = seqan3::sequence_file_input_default_traits_dna;
    using sequence_file_t = seqan3::sequence_file_input<traits_type, seqan3::fields<seqan3::field::seq>>;

    size_t t_syncmer = ceil((taxor_config.kmer_size - taxor_config.syncmer_size + 1) / 2.0);

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
                        // see count_minimizers() for why this downsampling scheme is used
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


/*!\brief Drive the chopper counting -> layout pipeline to produce the HIXF binning
 * layout file consumed later by build_hixf().
 *
 * Builds per-species file clusters (create_filename_clusters()), then counts
 * (approximate) k-mer/minimizer/syncmer set sizes per cluster - using
 * count_syncmers() if taxor_config.use_syncmer is set, count_minimizers() if a
 * minimizer window larger than the k-mer size is configured, or chopper's own
 * exact count_kmers() otherwise - before feeding the resulting counts into
 * chopper::layout::hierarchical_binning via determine_best_number_of_technical_bins(),
 * and finally writing the layout (header + body) to layout_config.output_filename.
 * \param taxor_config Build configuration (k-mer/window/syncmer sizes, scaling, thread count, ...).
 * \param orgs Parsed taxonomy entries; passed through to create_filename_clusters().
 * \param user_bin_map Output parameter populated by create_filename_clusters(), mapping
 * each genome file path to its index in \p orgs.
 */
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
    
    std::cout << "write Layout header" << std::endl << std::flush;

    // brief Write the output to the layout file.
    std::ofstream fout{layout_config.output_filename};
    chopper::layout::write_layout_header_to(layout_config, max_hixf_id, header_buffer.str(), fout);
    fout << output_buffer.str();

}

/*!\brief Build the HIXF index from the chopper layout written by create_layout(),
 * reattach species metadata to each resulting user bin, and store the finished
 * index to config.output_file_name.
 *
 * Translates \p config into a hixf::build_arguments struct and delegates the actual
 * filter construction to hixf::create_ixfs_from_chopper_pack() (see hixf/build/).
 * For every user bin produced, the owning genome file is matched back to its entry
 * in \p orgs via \p user_bin_map (see the comment at the `.at()` call below for why
 * a mismatch is treated as a hard error rather than silently defaulting), the
 * species' `user_bin` id and total sequence length are filled in, and a taxor_index
 * wrapping the built hixf::hierarchical_interleaved_xor_filter is written to disk
 * via store_index().
 * \param config Build configuration (kmer/window/syncmer sizes, scaling, thread count,
 * output path, use_syncmer flag, ...).
 * \param orgs Parsed taxonomy entries; user_bin and seq_len are filled in as a side effect.
 * \param user_bin_map Maps a genome file path (as produced by create_filename_clusters())
 * to its index in \p orgs.
 */
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
        args.t_syncmer = ceil((args.kmer_size - args.syncmer_size + 1) / 2.0);
	
    hixf::build_data data{};
   
    hixf::create_ixfs_from_chopper_pack(data, args);
    
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hixf.user_bins.num_user_bins(); ++i)
    {
        bin_path.push_back(std::vector<std::string>{data.hixf.user_bins.filename_of_user_bin(i)});
        // .at() (not operator[]): if this filename doesn't exactly match a key inserted
        // into user_bin_map (e.g. a formatting mismatch surviving the chopper count/layout/
        // pack round-trip), operator[] would silently default org_index to 0, clobbering
        // species #0's user_bin assignment and leaving the *actual* owner of this bin
        // without one - the corruption only surfaces later as an unrelated "map::at" crash
        // during search, on a different query hit. Fail loudly here instead, at build time.
        uint64_t org_index;
        try
        {
            org_index = user_bin_map.at(data.hixf.user_bins.filename_of_user_bin(i));
        }
        catch (std::out_of_range const &)
        {
            throw std::runtime_error{"Could not match built HIXF user bin filename '" +
                                      data.hixf.user_bins.filename_of_user_bin(i) +
                                      "' to any input genome file. The index would be built with corrupted "
                                      "species metadata; aborting instead of silently mis-assigning it."};
        }
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

/*!\brief Entry point for `taxor build`, invoked from main.cpp's dispatch in main().
 *
 * Registers and parses the build CLI options, runs sanity_checks(), parses every
 * taxonomy input file into \p orgs, computes the HIXF layout (create_layout()),
 * and builds+stores the HIXF index (build_hixf()). Each stage is wrapped in its
 * own try/catch so a failure is reported with a `[TAXOR BUILD ERROR]`-prefixed
 * message and a non-zero exit code without a stack of unrelated errors.
 * \param parser The `build` sub-parser, already selected by main.cpp's top-level parser.
 * \return 0 on success, -1 if CLI parsing/validation, taxonomy parsing, layout
 * creation, or index building fails.
 */
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

    try
    {
        for (std::string &f : config.input_files)
        {
            std::vector<taxonomy::Species> org_tmp = taxonomy::parse_refseq_taxonomy_file(f);
            orgs.insert(orgs.end(), org_tmp.begin(), org_tmp.end());
        }
    }
    catch (std::exception const &e)
    {
        std::cerr << "[TAXOR BUILD ERROR] " << e.what() << '\n';
        return -1;
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
    try
    {
        build_hixf(config, orgs, user_bin_map);
    }
    catch (std::exception const &e)
    {
        std::cerr << "[TAXOR BUILD ERROR] " << e.what() << '\n';
        return -1;
    }
    std::cout << "done!" << std::endl;

    return 0;
}

} // namespace chopper::layout
