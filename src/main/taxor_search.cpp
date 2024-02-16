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

namespace taxor::search
{

void set_up_subparser_layout(seqan3::argument_parser & parser, taxor::search::configuration & config)
{
    parser.info.version = "0.1.0";
    parser.info.author = "Jens-Uwe Ulrich";
    parser.info.email = "jens-uwe.ulrich@hpi.de";
    parser.info.short_description = "Queries a file of DNA sequences against an HIXF index";

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

void search_single(hixf::search_arguments & arguments, taxor_index<hixf_t> && index)
{
    double s_factor = 10.0;
    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};
    hixf::threshold::threshold thresholder;
    // map hixf user bin to species list index position
    std::map<size_t, size_t> user_bin_index{};
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
            user_bin_index.emplace(std::make_pair(index.species().at(i).user_bin, i));
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);
    seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    hixf::sync_out synced_out{arguments.out_file};

    {
        synced_out << "#QUERY_NAME\tACCESSION\tREFERENCE_NAME\tTAXID\tREF_LEN\tQUERY_LEN\tQHASH_COUNT\tQHASH_MATCH\tTAX_STR\tTAX_ID_STR\n";
    }

    
    std::mutex count_mutex;
    double mean_sum{0.0};
    uint16_t reads{0};
    auto worker = [&](size_t const start, size_t const end)
    {
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
            size_t fp_correction = hash_count * 0.003;
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
                uint64_t max_count = 0;
                for (auto && count : result)
                {
                    if (count.second > max_count)
                        max_count = count.second;
                }

                for (auto && count : result)
                {
                    // filter out counts that have less than max_count*0.8 matching hashes
                    if (static_cast<double>(count.second) < static_cast<double>(max_count) * 0.8)
                        continue;
                    result_string += id + '\t';
                    result_string += index.species().at(user_bin_index[count.first]).accession_id;
                    result_string += '\t';
                    result_string += index.species().at(user_bin_index[count.first]).organism_name;
                    result_string += '\t';
                    result_string += index.species().at(user_bin_index[count.first]).taxid;
                    result_string += '\t';
                    result_string += std::to_string(index.species().at(user_bin_index[count.first]).seq_len);
                    result_string += '\t';
                    result_string += std::to_string(seq.size());
                    result_string += '\t';
                    result_string += std::to_string(hash_count);
                    result_string += '\t';
                    result_string += std::to_string(count.second);
                    result_string += '\t';
                    result_string += index.species().at(user_bin_index[count.first]).taxnames_string;
                    result_string += '\t';
                    result_string += index.species().at(user_bin_index[count.first]).taxid_string;
                    result_string += '\n';
                }
            }
            count_mutex.lock();
            reads++;
            count_mutex.unlock();
            synced_out.write(result_string);
        }
    };
  
    for (auto && chunked_records : fin | seqan3::views::chunk(1024))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

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

void search_hixf(taxor::search::configuration const config)
{
    hixf::search_arguments search_args{};
	search_args.index_file = config.index_file; 
	search_args.query_file = config.query_file;
	search_args.out_file = config.report_file;
    search_args.threshold = config.threshold;
    search_args.threads = config.threads;
    search_args.seq_error_rate = config.error_rate;
    
	auto index = taxor_index<hixf_t>{};
    search_single(search_args, std::move(index));
}



int execute(seqan3::argument_parser & parser)
{
    taxor::search::configuration config;
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

    search_hixf(config);

    return 0;
}
}