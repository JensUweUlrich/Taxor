#pragma once

#include <string>

/*!\file taxor_search_configuration.hpp
 * \brief Holds the configuration struct populated from command-line options for `taxor search`.
 */

namespace taxor::search
{

//!\brief Configuration options for the `taxor search` subcommand, filled in by the argument parser.
struct configuration
{
    std::string index_file{}; //!< Raw (possibly comma-separated) --index-file argument as given on the command line.
    std::vector<std::string> index_file_list{}; //!< index_file split on ',' into individual index file paths.
    std::string query_file{}; //!< Raw (possibly comma-separated) --query-file argument as given on the command line.
    std::vector<std::string> query_file_list{}; //!< query_file split on ',' into individual query file paths.
    std::string report_file{}; //!< Output TSV file path for the search results.
    double threshold{-1.0}; //!< Fixed match-percentage threshold (0..1); if unset (<0), the statistical k-mer/syncmer model is used instead.
    double error_rate{0.04}; //!< Expected per-read sequencing error rate, used by the statistical threshold model.
    uint8_t threads{1u}; //!< Number of worker threads used to process query reads in parallel.
    bool output_verbose_statistics{false}; //!< Enables verbose statistics output (currently unused by taxor_search.cpp).
    bool debug{false}; //!< Enables debug output (currently unused by taxor_search.cpp).
};

}