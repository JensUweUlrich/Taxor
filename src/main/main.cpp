

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <filesystem>

#include <sys/resource.h>
#include <sys/time.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/dream_index/interleaved_binary_fuse_filter.hpp>


#include <cereal/archives/binary.hpp>

#include <zlib.h>
#include "StopClock.hpp"

#include "taxor_build.hpp"
#include "taxor_search.hpp"
#include "taxor_profile.hpp"


using namespace seqan3::literals;

/*!\file main.cpp
 * \brief Top-level entry point of the taxor executable: parses the global CLI
 * (which of the `build`, `search`, `profile` subcommands to run) and dispatches
 * to the corresponding subcommand's `execute()`, then reports CPU time and peak
 * memory usage.
 */

/*!\brief Return the process' cumulative CPU time (user + system) in seconds since start.
 * \return CPU time in seconds, as a floating point value with microsecond resolution.
 */
double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

/*!\brief Return the process' peak resident set size (maximum memory used so far).
 * \return Peak RSS in bytes.
 */
long getPeakRSS(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss * 1024;
}

/*!\brief Program entry point: parses the top-level `taxor` argument parser and
 * dispatches to the `build`, `search`, or `profile` subcommand, then prints
 * CPU time and peak RSS to stdout.
 * \param argc Argument count, forwarded to seqan3::argument_parser.
 * \param argv Argument vector, forwarded to seqan3::argument_parser.
 * \return The error code returned by the executed subcommand (0 on success,
 * -1 on a CLI parsing error, or the subcommand's own error code).
 */
int main(int argc, char const **argv)
{

	seqan3::argument_parser top_level_parser{"taxor", argc, argv, seqan3::update_notifications::off,
                                             {"build", "search", "profile"}};
    top_level_parser.info.version = "0.2.0";

    try
    {
        top_level_parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[TAXOR ERROR] " << ext.what() << '\n'; // customize your error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

    int error_code{};

    // dispatch to the subcommand selected on the command line; seqan3::argument_parser
    // already restricts app_name to one of "build"/"search"/"profile", so exactly one
    // branch below is taken
    if (sub_parser.info.app_name == std::string_view{"taxor-build"})
        error_code = taxor::build::execute(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"taxor-search"})
        error_code = taxor::search::execute(sub_parser);
	else if (sub_parser.info.app_name == std::string_view{"taxor-profile"})
        // taxor_profile.hpp is documented separately; only invoked here.
        error_code = taxor::profile::execute(sub_parser);

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	//std::cout << "Real time : " << ReadBouncerTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

    return error_code;

}

