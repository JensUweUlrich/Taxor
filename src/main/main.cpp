

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

#include <syncmer.hpp>
//#include "build.hpp"

#include "taxor_build.hpp"
#include "taxor_search.hpp"

#include <search/search_arguments.hpp>
#include <search/search_hixf.hpp>
#include <build/temp_hash_file.hpp>

#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>

#include <robin_hood.h>

using namespace seqan3::literals;

//std::string gtdb_root = "D:\\gtdb_genomes_reps_r202";

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
                                             {"build", "search"}};
    top_level_parser.info.version = "1.0.0";

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

	size_t peakSize = getPeakRSS();
	int peakSizeMByte = (int)(peakSize / (1024 * 1024));

	//std::cout << "Real time : " << ReadBouncerTime.elapsed() << " sec" << std::endl;
	std::cout << "CPU time  : " << cputime() << " sec" << std::endl;
	std::cout << "Peak RSS  : " << peakSizeMByte << " MByte" << std::endl;

    return error_code;

	// Print debuging of bulk count method to fid error when building ixf
	/*using sequence_file_t = seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::seq>>;
	hixf::build_arguments args{};
	std::vector<std::vector<size_t>> bins{};

	seqan3::interleaved_xor_filter<> ixf{64, 7000000};
	
	bool success = false;
	while (!success)
	{
		size_t bin_idx{0};
		for (int64_t pos = 1160; pos < 1170;++pos)
		{
			std::vector<size_t> hashes{};
			hixf::read_from_temp_hash_file(pos, hashes);
			if (hashes.empty())
				continue;
			//bins.emplace_back(hashes);
			success = iff.add_bin_elements(bin_idx, hashes);
			std::cout << bin_idx << "\t" << hashes.size() << std::endl;
			bin_idx++;
			if (!success)
				break;
		}
		if (!success)
		{
			iff.clear();
            iff.set_seed();
			std::cout << "Reset seed after bin " << bin_idx << std::endl;
		}
	}

	robin_hood::unordered_flat_set<size_t> read_hashes{};
	for (auto && [seq] : sequence_file_t{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/files.renamed/GCF_000839085.1_genomic.fna.gz"})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(args.shape,
                                                         seqan3::window_size{args.window_size},
                                                         seqan3::seed{hixf::adjust_seed(args.shape.count())}))
                    read_hashes.insert(hash);
	std::vector<size_t> c{};
    std::ranges::copy(read_hashes, std::back_inserter(c));
	typedef seqan3::interleaved_binary_fuse_filter<>::counting_agent_type<u_int64_t> TIFFAgent;
	TIFFAgent ixf_count_agent = iff.counting_agent< uint64_t >();

	auto result = ixf_count_agent.bulk_count(c);
    seqan3::debug_stream << "Root result: " << result << "\n";
	
	*/

	
	//create_layout();

	/*hixf::search_arguments search_args{};
	search_args.index_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor_kmer.hixf"};
	search_args.query_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/files.renamed/GCF_000839085.1_genomic.fna.gz"};
	search_args.out_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor.hixf_search.out"};
	search_args.compute_syncmer = false;
	search_args.threshold = 0.2;
	hixf::search_hixf(search_args);
*/

	return 1;
	
/*
	
	
//	multi_interleaved_xor_filter mixf = load_multi_ixf_index("GTDB_r202.ixf");
//	std::cout << "Filter loaded successfully" << std::endl;
	
	//std::vector<std::vector<Species*>> bins2vec(mixf.count_single_filter());

	//bins_to_species(mixf, bins2vec);
	

	double ratio_sum = 0.0;
	uint64_t nr_reads = 0;
	std::vector<double> sensitivity(mixf.species_vector.size());
	std::fill(sensitivity.begin(), sensitivity.end(), 0.0);
	std::vector<uint64_t> classified(mixf.species_vector.size());
	std::fill(classified.begin(), classified.end(), 0);
	int count = 0;
	for ( const auto& [seq, header, qual] : fin)
	{

		if (seq.size() < 1000)
			continue;

		// rhv holds hash vectors for forward and reverse complement
		std::vector<size_t> rhv = hashing::seq_to_syncmers(kmer_size, seq, sync_size, t_syncmer);
		uint64_t strobe_nr = rhv.size();
//		TInterval ci = calculateCI(0.1, kmer_size, strobe_nr, 0.95);
//		uint64_t threshold = strobe_nr - ci.second;
		std::vector<uint64_t> max_found(mixf.species_vector.size());
		std::fill(max_found.begin(), max_found.end(), 0);
//		std::cout << rhv[0].first.size() << "\t" << seq.size() << std::endl;
		
			std::vector<TIXFAgent::counting_vector> count_vectors = mixf.bulk_count(rhv);

			std::vector<uint64_t> result(mixf.species_vector.size());
			for (uint64_t i = 0; i < result.size(); ++i)
			{
				result[i] = 0;
				if (mixf.species_vector[i].filter_index == UINT16_MAX)
					continue;
				for (uint64_t spec_bin = mixf.species_vector[i].first_bin; spec_bin <= mixf.species_vector[i].last_bin; ++spec_bin)
				{
					result[i] += count_vectors[mixf.species_vector[i].filter_index][spec_bin];
					if (i == 355)
					{
						std::cout << spec_bin << "\t" << count_vectors[mixf.species_vector[i].filter_index][spec_bin] << "\t" << strobe_nr << std::endl;
					}
				}
			}
			
			// find maximum matches between forward and reverse complement matches
			for (int i = 0; i < max_found.size(); ++i)
			{
				if (result[i] > max_found[i])
					max_found[i] = result[i];
			}
		

		//break;

		double max_ratio = 0.0;

		std::vector<std::pair<uint64_t, double>> potential_indexes{};
		for (int i = 0; i < max_found.size(); ++i)
		{
			
			double ratio = (double) max_found[i] / (double) strobe_nr;
			sensitivity[i] += ratio;

			// Why does 0.05 work so well?
			if ( ratio > 0.1)
			{

				if (ratio > max_ratio)
					max_ratio = ratio;

				potential_indexes.emplace_back(std::make_pair(i, ratio));
			}
		}

		int spec_count = 0;
		for (const auto & p : potential_indexes)
		{
			if (p.second >= max_ratio)
			{
				classified[p.first] += 1;
//				if (++spec_count > 1)
//					std::cout << header << std::endl;

			}
		}
		
		nr_reads++;
//		if (nr_reads == 10)
//			break;
	}

	for (int i=0; i < sensitivity.size(); ++i)
	{
		sensitivity[i] /= (double) nr_reads;
	}


	for (uint64_t i = 0; i < mixf.species_vector.size(); ++i)
	{
		
		
		if (classified[i] > 1)
		{
			std::cout << i << "\t" << mixf.species_vector[i].name << "\t" << classified[i] << "\t" << sensitivity[i] << std::endl;
		}
	}

//	seqan3::debug_stream << sensitivity << "\n";
//	seqan3::debug_stream << classified << "\n";
	seqan3::debug_stream << nr_reads << "\n";

	return 0;
	*/
}

