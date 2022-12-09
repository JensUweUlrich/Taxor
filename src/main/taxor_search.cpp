#include <seqan3/core/debug_stream.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>

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
 /*                     "Provide the prefix you used for the output prefix in the chopper count --output-prefix option. "
                      "If you have different means of estimating the k-mer counts of your input data, make sure that a "
                      "file [INPUT-PREFIX].count exists. It needs to be tab-separated and consist of two columns: "
                      "\"[filepath] [tab] [weight/count]\".",*/
                      seqan3::option_spec::required);

    parser.add_option(config.query_file, '\0', "query-file", "file containing sequences to query against the index");

    parser.add_option(config.report_file, '\0', "output-file", "A file name for the resulting output.");

    parser.add_option(config.threads,
                      '\0', "threads",
                      "The number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{static_cast<size_t>(1), static_cast<size_t>(32)});

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
    //constexpr bool is_ibf = std::same_as<index_t, raptor_index<index_structure::ibf>>
    //                     || std::same_as<index_t, raptor_index<index_structure::ibf_compressed>>;

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
        arguments.shape = index.shape();
        arguments.shape_size = index.kmer_size();
        // TODO: thresholding should be set based on used pattern
        hixf::threshold_parameters param = arguments.make_threshold_parameters();
        if (arguments.compute_syncmer)
            param.fracminhash = true;
        param.seq_error_rate = 0.05;
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
        size_t position{};
        std::string line{};
        for (auto const & file_list : arguments.bin_path)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            for (auto const & filename : file_list)
            {
                line += filename;
                line += ',';
            }
            line.back() = '\n';
            synced_out << line;
            ++position;
        }
        synced_out << "#QUERY_NAME\tREFERENCE_NAME\n";
    }

    

    auto worker = [&](size_t const start, size_t const end)
    {
        auto counter = [&index]()
        {
            return index.ixf().membership_agent();
        }();
        std::string result_string{};
        std::vector<uint64_t> hashes;

        // TODO: choose between minimizer and syncmers
        auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{hixf::adjust_seed(arguments.shape_weight)});

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            if (arguments.compute_syncmer)
            {
                seqan3::dna5_vector dna5_vector{seq.begin(), seq.end()};
                std::vector<uint64_t> strobe_hashes = hashing::seq_to_syncmers(index.kmer_size(),dna5_vector, index.syncmer_size(), index.t_syncmer());
                hashes.assign(std::make_move_iterator(strobe_hashes.begin()), std::make_move_iterator(strobe_hashes.end()));
            }
            else
            {
                auto minimiser_view = seq | hash_adaptor | std::views::common;
                hashes.assign(minimiser_view.begin(), minimiser_view.end());
            }
            size_t const hash_count{hashes.size()};
            size_t threshold = thresholder.get(hash_count, (double)hash_count / ((double)seq.size() - (double)index.kmer_size() + 1.0));
            std::cout << "Threshold: " << threshold << std::endl;
            std::cout << "Minimizer count: " << hash_count << std::endl;

            auto & result = counter.bulk_contains(hashes, threshold); // Results contains user bin IDs
            // write one line per reference match 
            if (result.empty())
            {
                result_string += id + '\t';
                result_string += "-\n";
            }
            else{
                for (auto && count : result)
                {
                    result_string += id + '\t';
                    result_string += index.species().at(user_bin_index[count]).organism_name;
                    result_string += '\n';
                }
            }
            /*
            if (auto & last_char = result_string.back(); last_char == ',')
                last_char = '\n';
            else
                result_string += '\n';
            */
            synced_out.write(result_string);
            break;
        }
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        hixf::do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    // GCOVR_EXCL_START
    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed << std::setprecision(2) << index_io_time << '\t' << reads_io_time << '\t'
                    << compute_time;
    }
    // GCOVR_EXCL_STOP
}

void search_hixf(taxor::search::configuration const config)
{
    hixf::search_arguments search_args{};
	search_args.index_file = config.index_file; //std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor_kmer.hixf"};
	search_args.query_file = config.query_file;//std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/files.renamed/GCF_000839085.1_genomic.fna.gz"};
	search_args.out_file = config.report_file;//std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor.hixf_search.out"};
    search_args.threads = config.threads;
    // TODO: should be set after loading the index
	//search_args.compute_syncmer = false;
	//search_args.threshold = 0.2;
    
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