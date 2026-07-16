// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#include <seqan3/core/debug_stream.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/fracmin_hash.hpp>

#include <search/search_arguments.hpp>
#include <search/sync_out.hpp>
#include <search/threshold.hpp>
#include <search/do_parallel.hpp>

#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>

#include <syncmer.hpp>

#include "index.hpp"
#include "load_index.hpp"
#include "taxor_search.hpp"
#include "taxor_search_configuration.hpp"

/*!\file taxor_search.cpp
 * \brief Implements the `taxor search` subcommand: loads a prebuilt HIXF index, computes
 * k-mer/syncmer hashes for each query read, queries the index's membership_agent for
 * matching reference user bins, applies a statistical (or fixed) threshold to decide which
 * hits are significant, and writes a TSV report of read-to-reference matches.
 */

namespace taxor::search
{

/*!\brief Registers all `taxor search` command-line options on the given argument_parser.
 * \param parser The seqan3::argument_parser to configure.
 * \param config The configuration struct whose members are bound to the individual options.
 */
void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::search::configuration & config)
{
    parser.info.version = "0.2.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Queries files of DNA sequences against a list of HIXF index files";

    parser.info.description.emplace_back("Query sequences against the taxor HIXF index structure");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.index_file,
                      '\0', "index-file", "taxor index file containing HIXF index and reference sequence information",
                      seqan3::option_spec::required);

    parser.add_option(config.query_file, '\0', "query-file", "file containing sequences to query against the index");

    parser.add_option(config.report_file, '\0', "output-file", "A file name for the resulting output.");

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

    parser.add_option(config.threshold,
                    '\0', "percentage",
                    "If set, this threshold is used instead of the k-mer/syncmer models.",
                    seqan3::option_spec::standard,
                    seqan3::arithmetic_range_validator{static_cast<double>(0.0), static_cast<double>(1.0)});

    parser.add_option(config.error_rate,
                    '\0', "error-rate",
                    "Expected error rate of reads that will be queried",
                    seqan3::option_spec::standard,
                    seqan3::arithmetic_range_validator{static_cast<double>(0.0), static_cast<double>(1.0)});

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

/*!\brief Splits a string into a list of substrings on a delimiter character.
 * \param str The string to split (unused after the call; not modified).
 * \param delimiter The character to split on.
 * \return The list of substrings between occurrences of delimiter.
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

/*!\brief Validates the parsed configuration and splits comma-separated index/query file lists.
 *
 * Splits config.index_file and config.query_file (each possibly a comma-separated list) into
 * config.index_file_list and config.query_file_list, verifies that every referenced file exists,
 * and - when more than one index file is given - loads each index just to compare its
 * k-mer/syncmer/window/scaling parameters, throwing if the index files are not all built with
 * the same scheme (since results from mismatched indexes could not be combined meaningfully).
 * \param config The configuration to validate and populate the *_list members of.
 */
void sanity_checks(taxor::search::configuration & config)
{
    // enable using seveal index files
    config.index_file_list = str_split(config.index_file, ',');

    // check whether all index files exist
    // check whether all index files use the same k-mer scheme
    bool compute_syncmer{false};
    uint16_t scaling{1u};
    uint32_t window_size{1u};
    uint8_t kmer_size{1u};
    uint8_t syncmer_size{1u};
    uint8_t t_syncmer{1u};
    for (std::string & f : config.index_file_list)
    {   
        std::filesystem::path filepath{f};
        if (!std::filesystem::exists(filepath))
            throw seqan3::argument_parser_error{"Please check the given index file(s). \nThe following index file does not exist: " + f};
        
        if (config.index_file_list.size() > 1)
        {
            auto index = taxor_index<hixf_t>{};
            double index_io_time{0.0};
            load_index(index, filepath, index_io_time);
            if (kmer_size == 1)
            {
                kmer_size = index.kmer_size();
                window_size = index.window_size();
                scaling = index.scaling();
                syncmer_size = index.syncmer_size();
                t_syncmer = index.t_syncmer();
                compute_syncmer = index.use_syncmer();
                continue;
            }
            
            if (kmer_size != index.kmer_size() || window_size != index.window_size() || scaling != index.scaling() ||
                syncmer_size != index.syncmer_size() || t_syncmer != index.t_syncmer() || compute_syncmer != index.use_syncmer())
                throw seqan3::argument_parser_error{"At least two index files have been created with different kmer selection schemes.\n Please provide only index files using the same kmer-/syncmer-/window-size!"};
            
        }
    }

    // enable using several input directories
    config.query_file_list = str_split(config.query_file, ',');
    
    // check whether all input files exist
    for (std::string & f : config.query_file_list)
    {   
        std::filesystem::path filepath{f};
        if (!std::filesystem::exists(filepath))
            throw seqan3::argument_parser_error{"Please check the given input query files. \nThe following query file does not exist: " + f};
    }


}

/*!\brief Searches all reads in one query file against one HIXF index and appends matches to outstrm.
 *
 * Loads the index (asynchronously, overlapped with reading the first chunk of query reads),
 * then processes the query file in chunks of 1024 records. For each chunk, hashes are computed
 * per read (syncmers or minimiser k-mers, depending on the index) and the chunk is searched in
 * parallel via hixf::do_parallel; each worker thread queries the index's membership_agent and
 * writes matching reference lines for its share of the reads to the shared outstrm.
 * \param arguments Search parameters (index/query file paths, thread count, threshold settings, ...);
 *                   populated further from the loaded index's k-mer/syncmer scheme once available.
 * \param index The (not-yet-loaded) taxor_index to load and search against; taken by rvalue reference
 *              since this function loads it in place.
 * \param outstrm Shared, mutex-guarded output stream that results are appended to.
 */
void search_single(hixf::search_arguments & arguments, taxor_index<hixf_t> && index, hixf::sync_out &outstrm)
{
    double s_factor = 10.0; // currently unused: only referenced from the commented-out scaling lambda below
    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};
    hixf::threshold::threshold thresholder;
    // map hixf user bin to species list index position
    std::map<size_t, size_t> user_bin_index{};
    // Loading the (potentially large) index from disk is independent of reading query records,
    // so it is kicked off asynchronously here and only waited on just before the first chunk is
    // searched (see index_load_checked below), overlapping index I/O with query file I/O.
    auto cereal_worker = [&]()
    {
        load_index(index, arguments, index_io_time);
        arguments.compute_syncmer = index.use_syncmer();
        arguments.shape_size = index.kmer_size();
        arguments.window_size = index.window_size();
        arguments.scaling = index.scaling();
        arguments.shape = seqan3::shape{seqan3::ungapped{arguments.shape_size}};
        hixf::threshold_parameters param = arguments.make_threshold_parameters();
        thresholder = hixf::threshold::threshold{param};
        for (size_t i = 0; i < index.species().size(); ++i)
        {
            user_bin_index.emplace(std::make_pair(index.species().at(i).user_bin, i));
            //if (index.species().at(i).taxid.compare("160232") == 0)
            //            std::cerr << index.species().at(user_bin_index.at(count.first)).organism_name << "\t" << count.first << "\t" << count.second << std::endl;

        }
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);
    seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

/*    hixf::sync_out synced_out{arguments.out_file};

    {
        synced_out << "#QUERY_NAME\tACCESSION\tREFERENCE_NAME\tTAXID\tREF_LEN\tQUERY_LEN\tQHASH_COUNT\tQHASH_MATCH\tTAX_STR\tTAX_ID_STR\n";
    }
*/
    
    std::mutex count_mutex;
    double mean_sum{0.0}; // currently unused: never assigned to or read outside its declaration
    uint16_t reads{0};
    // Some already-built indexes may have a user bin with no matching species metadata
    // (a data-integrity bug in the build step, now fixed - see taxor_build.cpp). Rather
    // than aborting the whole search on the first affected match, skip just that match and
    // warn once (shared/thread-safe across worker threads) so the user knows to rebuild.
    std::atomic_flag warned_missing_user_bin = ATOMIC_FLAG_INIT;
    /*!\brief Per-chunk worker: hashes and searches records [start, end) of the current chunk.
     *
     * Invoked once per thread by hixf::do_parallel, each call handling a disjoint sub-range of
     * `records`. Thread-safety within this worker relies on a few things actually present in the
     * code: (1) each invocation constructs its own local `counter` (membership_agent), so worker
     * threads never share one membership_agent's internal state; (2) `hashes` and `result_string`
     * are local to the lambda invocation (effectively per-thread buffers reused across the loop
     * body); (3) the shared `reads` counter is only mutated under count_mutex; (4) the one-time
     * "missing user bin" warning is guarded by the std::atomic_flag warned_missing_user_bin so it
     * is printed at most once even though multiple threads may hit the condition concurrently;
     * (5) all writes to the shared `outstrm` go through hixf::sync_out::write, which internally
     * locks its own mutex around the underlying std::ofstream.
     * \param start Index of the first record (in `records`) to process.
     * \param end Index one past the last record (in `records`) to process.
     */
    auto worker = [&](size_t const start, size_t const end)
    {
        // Constructed fresh per worker invocation (i.e. per thread, per chunk) rather than
        // shared, so concurrent threads never contend on or corrupt one membership_agent's
        // internal state.
        auto counter = [&index]()
        {
            return index.ixf().membership_agent();
        }();
        std::string result_string{};
        std::vector<uint64_t> hashes;
        
        /*auto scaling = [s_factor](uint64_t i) { uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(i);
                                                return double(v) <= double(UINT64_MAX) / s_factor; 
                                              } ;
        */

        auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{hixf::adjust_seed(arguments.shape.count())});      

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            result_string.clear();

            
            // Compute the read's hash set, either as syncmers or as classic minimiser k-mers,
            // depending on how the index was built (arguments.compute_syncmer, set from the
            // loaded index). When arguments.scaling > 1, apply FracMinHash-style subsampling:
            // a hash is kept only if wyhash(hash) falls below UINT64_MAX/scaling, i.e. roughly
            // a 1/scaling fraction of hashes are retained (must mirror the same scaling scheme
            // used when the index was built, or matches would be missed/skewed).
            if (arguments.compute_syncmer)
            {
                seqan3::dna5_vector dna5_vector{seq.begin(), seq.end()};
                ankerl::unordered_dense::set<size_t> tmp = hashing::seq_to_syncmers(index.kmer_size(),dna5_vector, index.syncmer_size(), index.t_syncmer());
                if (arguments.scaling > 1)
                {
                    for (auto &hash : tmp)
                    {
                        uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                        if (double(v) <= double(UINT64_MAX) / double(arguments.scaling))
                        {
                            hashes.push_back(hash);
                        }
                    }
                }
                else
                {
			        hashes.assign(std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }

            }
            else
            {
                for (auto hash :  seq | hash_adaptor)
                {
                    if (arguments.scaling > 1)
                    {
                        uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                        if (double(v) <= double(UINT64_MAX) / double(arguments.scaling))
                        {
                            hashes.push_back(hash);
                        }
                    }
                    else
                    {
                        hashes.push_back(hash);
                    }
                }
                //auto minimiser_view = seq | hash_adaptor | std::views::common;
                //hashes.assign(minimiser_view.begin(), minimiser_view.end());

            }
            size_t const hash_count{hashes.size()};
            size_t fp_correction = hash_count * 0.003; // currently unused: computed but never read afterward
            // thresholder.get() derives the minimum number of matching hashes required for a hit
            // to count as significant, from the observed hash_count and the read's approximate
            // per-hash "coverage" (hash_count relative to the number of possible k-mer positions
            // in the read), combined with the configured error_rate - see hixf/search/threshold.hpp.
            size_t threshold = thresholder.get(hash_count, (double)hash_count / ((double)seq.size() - (double)index.kmer_size() + 1.0));

            auto & result = counter.bulk_contains(hashes, threshold); // Results contains user bin IDs
            hashes.clear();
            // write one line per reference match 
            if (result.empty())
            {
                result_string += id + '\t';
                result_string += "-\t-\t-\t-\t";
                result_string += std::to_string(seq.size()) + "\n";
            }
            else{
                // Report near-best hits rather than only the single best one: find the highest
                // hash-match count among all candidate user bins for this read, then (below)
                // keep every user bin within 80% of that best count.
                uint64_t max_count = 0;
                for (auto && count : result)
                {
                    if (count.second > max_count)
                        max_count = count.second;
                }

                bool wrote_any = false;
                for (auto && count : result)
                {
                    // filter out counts that have less than max_count*0.8 matching hashes
                    if (static_cast<double>(count.second) < static_cast<double>(max_count) * 0.8)
                        continue;

                    auto ub_it = user_bin_index.find(count.first);
                    if (ub_it == user_bin_index.end())
                    {
                        // Known data-integrity issue in some already-built indexes (fixed at
                        // the source in taxor_build.cpp): a user bin with no matching species
                        // metadata. Skip this one match instead of aborting the whole search.
                        if (!warned_missing_user_bin.test_and_set())
                            std::cerr << "[TAXOR SEARCH WARNING] Index contains a user bin (id "
                                      << count.first << ") with no matching species metadata; "
                                      << "skipping affected match(es). Consider rebuilding the index.\n";
                        continue;
                    }

                    auto const & species = index.species().at(ub_it->second);
                    result_string += id + '\t';
                    result_string += species.accession_id;
                    result_string += '\t';
                    result_string += species.organism_name;
                    result_string += '\t';
                    result_string += species.taxid;
                    result_string += '\t';
                    result_string += std::to_string(species.seq_len);
                    result_string += '\t';
                    result_string += std::to_string(seq.size());
                    result_string += '\t';
                    result_string += std::to_string(hash_count);
                    result_string += '\t';
                    result_string += std::to_string(count.second);
                    result_string += '\t';
                    result_string += species.taxnames_string;
                    result_string += '\t';
                    result_string += species.taxid_string;
                    result_string += '\n';
                    wrote_any = true;
                }

                // if every match for this read happened to hit a corrupted/unmapped user
                // bin, still emit a line for it (as "no match") rather than silently
                // dropping the read from the output entirely
                if (!wrote_any)
                {
                    result_string += id + '\t';
                    result_string += "-\t-\t-\t-\t";
                    result_string += std::to_string(seq.size()) + "\n";
                }
            }
            count_mutex.lock();
            reads++;
            count_mutex.unlock();
            outstrm.write(result_string);
        }
    };
  
    bool index_load_checked = false;
    for (auto && chunked_records : fin | seqan3::views::chunk(1024))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        // future::wait() never propagates exceptions raised in the async task (unlike get()),
        // so a failed/corrupt index load would otherwise pass silently and the search below
        // would run against a partially-initialized index. future::get() may only be called
        // once (it consumes the future's shared state), so only do this on the first chunk -
        // once index loading has completed and been checked, it stays completed for the rest
        // of the (possibly many) chunks in this loop.
        if (!index_load_checked)
        {
            cereal_handle.get();
            index_load_checked = true;
        }

        hixf::do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed << std::setprecision(2) << index_io_time << '\t' << reads_io_time << '\t'
                    << compute_time;
    }
    
}

/*!\brief Runs search_single for every combination of query file and index file in config.
 *
 * Opens the shared TSV output file once (writing the header line), then loops over every
 * query file and, for each one, every index file, running a fresh search against each. Note
 * that this means each query file is searched independently against every index file (results
 * are not merged/deduplicated across indexes).
 * \param config The validated search configuration (see sanity_checks), taken by value.
 */
void search_hixf(taxor::search::configuration const config)
{
    hixf::sync_out synced_out{config.report_file};
    synced_out << "#QUERY_NAME\tACCESSION\tREFERENCE_NAME\tTAXID\tREF_LEN\tQUERY_LEN\tQHASH_COUNT\tQHASH_MATCH\tTAX_STR\tTAX_ID_STR\n";
    for (std::string query : config.query_file_list)
    {
        for (std::string hixf_file : config.index_file_list)
        {
            hixf::search_arguments search_args{};
            search_args.index_file = hixf_file; 
            search_args.query_file = query;
    //	    search_args.out_file = config.report_file;
            search_args.threshold = config.threshold;
            search_args.threads = config.threads;
            search_args.seq_error_rate = config.error_rate;
        
            auto index = taxor_index<hixf_t>{};
            search_single(search_args, std::move(index), synced_out);
        }
    }
}



/*!\brief Parses command-line arguments and runs the `taxor search` subcommand.
 *
 * Registers the subcommand's options, parses and validates them, then runs the search over
 * all configured query/index file combinations. Argument-parsing/validation errors and any
 * exception raised during the search itself are each caught separately and reported to stderr.
 * \param parser The seqan3::argument_parser instance used to parse the subcommand's options.
 * \return 0 on success, -1 if argument parsing/sanity checks fail or the search itself throws.
 */
int execute(seqan3::argument_parser & parser)
{
    taxor::search::configuration config;
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
        std::cerr << "[TAXOR SEARCH ERROR] " << ext.what() << '\n';
        return -1;
    }

    try
    {
        search_hixf(config);
    }
    catch (std::exception const & e)
    {
        std::cerr << "[TAXOR SEARCH ERROR] " << e.what() << '\n';
        return -1;
    }

    return 0;
}
}