
#include <search/search_arguments.hpp>
#include <search/search_hixf.hpp>

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

    parser.info.description.emplace_back("Creates an HIXF index using either kmers, syncmers or FracMinHash sketches");

    parser.add_subsection("Main options:");
    // -----------------------------------------------------------------------------------------------------------------
    parser.add_option(config.index_file,
                      '\0', "index-file", "",
 /*                     "Provide the prefix you used for the output prefix in the chopper count --output-prefix option. "
                      "If you have different means of estimating the k-mer counts of your input data, make sure that a "
                      "file [INPUT-PREFIX].count exists. It needs to be tab-separated and consist of two columns: "
                      "\"[filepath] [tab] [weight/count]\".",*/
                      seqan3::option_spec::required);

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

void search_hixf()
{
    hixf::search_arguments search_args{};
	search_args.index_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor_kmer.hixf"};
	search_args.query_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/files.renamed/GCF_000839085.1_genomic.fna.gz"};
	search_args.out_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor.hixf_search.out"};
	search_args.compute_syncmer = false;
	search_args.threshold = 0.2;
	hixf::search_hixf(search_args);
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


    return 0;
}
}