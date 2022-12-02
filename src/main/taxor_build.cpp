#include <iostream>

#include <chopper/count/read_data_file.hpp>
#include <chopper/detail_apply_prefix.hpp>
#include <chopper/count/check_filenames.hpp>
#include <chopper/count/count_kmers.hpp>
#include <chopper/layout/configuration.hpp>
#include <chopper/layout/data_store.hpp>
#include <chopper/layout/filenames_data_input.hpp>
#include <chopper/layout/hierarchical_binning.hpp>
#include <chopper/layout/output.hpp>

#include <build/chopper_build.hpp>
#include <build/build_arguments.hpp>

#include "taxor_build.hpp"
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
    parser.add_list_item("", "Example count file:");
    parser.add_list_item("", "```");
    parser.add_list_item("", "/absolute/path/to/file1.fasta     500");
    parser.add_list_item("", "/absolute/path/to/file2.fa.gz     600");
    parser.add_list_item("", "```");

    
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


void create_layout(taxor::build::configuration const taxor_config)
{
	chopper::count::configuration config{};
	config.data_file = taxor_config.input_file_name;
	config.k = taxor_config.kmer_size;
	config.threads = taxor_config.threads;
	config.output_prefix = "chopper";
	/*auto filename_clusters = chopper::count::read_data_file(config);

	chopper::detail::apply_prefix(config.output_prefix, config.count_filename, config.sketch_directory);
    chopper::count::check_filenames(filename_clusters, config);
    chopper::count::count_kmers(filename_clusters, config);
*/
	chopper::layout::configuration layout_config;
    chopper::layout::data_store data;

	layout_config.input_prefix = config.output_prefix;
	
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

void build_hixf(taxor::build::configuration const config)
{

	hixf::build_arguments args{};
	args.bin_file = std::filesystem::path{"binning.out"};
	args.out_path = config.output_file_name;
    args.kmer_size = config.kmer_size;
    args.syncmer_size = config.syncmer_size;
    args.threads = config.threads;
	args.compute_syncmer = config.use_syncmer;
	hixf::chopper_build(args);
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

    create_layout(config);

    build_hixf(config);


    return 0;
}

} // namespace chopper::layout
