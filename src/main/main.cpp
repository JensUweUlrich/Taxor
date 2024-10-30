

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

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long getPeakRSS(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_maxrss * 1024;
}

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

    if (sub_parser.info.app_name == std::string_view{"taxor-build"})
        error_code = taxor::build::execute(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"taxor-search"})
        error_code = taxor::search::execute(sub_parser);
	else if (sub_parser.info.app_name == std::string_view{"taxor-profile"})
        error_code = taxor::profile::execute(sub_parser);

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	//std::cout << "Real time : " << ReadBouncerTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

    return error_code;

}

