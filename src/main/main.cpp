

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <filesystem>


#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <cereal/archives/binary.hpp>

#include <zlib.h>
#include "StopClock.hpp"
#include "interleaved_binary_fuse_filter.hpp"

#include "index.hpp"
#include "build.hpp"

#include <chopper_build.hpp>
#include <build_arguments.hpp>

using namespace seqan3::literals;

//std::string gtdb_root = "D:\\gtdb_genomes_reps_r202";



void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines) 
{
	std::ifstream ifs(fname);
	if (ifs.fail()) {
		std::cerr << "error" << std::endl;
		return;
	}
	std::string line;
	while (getline(ifs, line)) {
		std::stringstream ss(line);
		std::vector<std::string> item;
		std::string tmp;
		while (getline(ss, tmp, '\t')) {
			item.push_back(tmp);
		}
		lines.push_back(item);
	}
}

inline uint64_t get_bin_size( double false_positive, uint16_t hash_functions, uint64_t max_hashes )
{
    return std::ceil( max_hashes
                      * ( -hash_functions / std::log( 1 - std::exp( std::log( false_positive ) / hash_functions ) ) ) );
}

 /*void getGenomeRepsByGenus(std::string& filename, TGenVec& genvec)
 {
	 std::vector<std::vector<std::string> > gtdb_reps_ncbi{};
	 read_tsv(filename, gtdb_reps_ncbi);
	 TGenMap genera{};
	 TGenMap::iterator gen_it;
	 for (std::vector<std::string> line : gtdb_reps_ncbi)
	 {
		 if (line[0].compare("organism_name") == 0)
			 continue;

		 Species sp{};
		 sp.accession_id = line[0];
		 sp.organism_name = line[1];
		 sp.taxonomy = line[3];
		 std::vector<std::string> item;
		 std::string tmp;
		 std::stringstream taxonomy(line[3]);
		 std::string genus{};
		 while (getline(taxonomy, tmp, ';'))
		 {
			 if (tmp.substr(0, 1).compare("g") == 0)
			 {
				 genus = tmp.substr(3);
			 }
			 else if (tmp.substr(0, 1).compare("s") == 0)
			 {
				 sp.name = tmp.substr(3);
			 }
		 }
		 gen_it = genera.find(genus);
		 if (gen_it != genera.end())
		 {
			 (*gen_it).second.emplace_back(sp);
		 }
		 else
		 {
			 std::vector<Species> sp_vec{};
			 sp_vec.emplace_back(sp);
			 genera.insert(std::make_pair(genus, std::move(sp_vec)));
		 }
	 }

	genvec.insert(genvec.end(), genera.begin(), genera.end());

 }
*/
 void load_species_from_file(std::string& filename, std::vector<Species> & species)
 {
	std::vector<std::vector<std::string> > gtdb_reps_ncbi{};
	read_tsv(filename, gtdb_reps_ncbi);
	uint64_t species_counter = 0;

	for (std::vector<std::string> line : gtdb_reps_ncbi)
	{
		if (line[0].compare("organism_name") == 0)
			continue;

		Species sp{};
		sp.accession_id = line[0];
		sp.organism_name = line[1];
		sp.taxonomy = line[3];
		 
		std::string tmp;
		std::stringstream taxonomy(line[3]);
		std::string genus{};
		while (getline(taxonomy, tmp, ';'))
		{
			if (tmp.substr(0, 1).compare("s") == 0)
			{
				sp.name = tmp.substr(3);
			}
		}
		 
		species.emplace_back(sp);

		if (++species_counter >= breakpoint)
			break;
	 }
 }

uint64_t construct_and_query_ixf(size_t bins, size_t elements_per_bin)
{
	StopClock ixf_construct{};
	ixf_construct.start();
	//ulrich::interleaved_xor_filter<> ixf(elems);
	seqan3::interleaved_xor_filter<> ixf(bins, elements_per_bin);
	std::vector<uint64_t> elems{};
	while (true)
	{
		bool success = true;
		for (int e = 0; e < bins ; ++e)
		{
			std::vector<uint64_t> tmp{};
			for (uint64_t i = 0; i < elements_per_bin; ++i)
			{
				uint64_t key = (e*elements_per_bin) + i;
				tmp.emplace_back(key);
			}
			success = ixf.add_bin_elements(e, tmp);
			if (!success)
			{
				ixf.clear();
				ixf.set_seed();
				std::cout << e << std::endl;
				break;
			}
			if (e == 2)
				elems=std::move(tmp);
		}

		if (success)
			break;
	}
	ixf_construct.stop();
	std::cout << "Built IXF in " << ixf_construct.elapsed() << " seconds" << std::endl;
	std::cout << ixf.bin_count() << std::endl;
	std::cout << (double) ixf.bit_size() / (double) ixf.bin_count() << std::endl;
	StopClock ixf_query{};
	ixf_query.start();
	typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;
/*
	auto agent = ixf.membership_agent();
    auto & result = agent.bulk_contains(5768566);
*/    //seqan3::debug_stream << result << '\n';
	TIXFAgent ixf_count_agent = ixf.counting_agent< uint64_t >();
	//auto result = count_agent.bulk_count(readHs);
	auto ixf_result = ixf_count_agent.bulk_count(elems);
	ixf_query.stop();
//	seqan3::debug_stream << ixf_result << "\n";
	double fpr {0.0};
	for (int i =0; i < ixf_result.size(); ++i)
	{
		if (i == 2)
		{
			std::cout << ixf_result[i] << "/" << elems.size() << std::endl;
			return ixf_result[i];
			continue;
		}
		fpr += (double) ixf_result[i] / (double) elements_per_bin;
	}
	fpr /= (double) (bins - 1);
	std::cout << "FPR of the IXF: " << fpr << std::endl;
	double mbytes = (double) ixf.bit_size() / (double) 8388608;
	std::cout << "Size of the IXF: " << mbytes << " MBytes" << std::endl;
	std::cout << "Queried " << elems.size() <<" keys in IXF in " << ixf_query.elapsed() << " seconds" << std::endl;

	StopClock ixf_store{};
	ixf_store.start();
	std::string filter_file = "test.ixf";
	
	std::ofstream               os( filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( ixf );  
	ixf_store.stop();
	std::cout << "Stored IXF in " << ixf_store.elapsed() << " seconds" << std::endl;

}

void construct_and_query_iff(size_t bins, size_t elements_per_bin)
{
	
	//ulrich::interleaved_xor_filter<> ixf(elems);
//	ulrich::interleaved_xor_filter<> ixf(bins, elements_per_bin);
	std::vector<std::vector<uint64_t>> elems{};
	
	for (int e = 0; e < bins ; ++e)
	{
		std::vector<uint64_t> tmp{};
		for (uint64_t i = 0; i < elements_per_bin; ++i)
		{
			uint64_t key = (e*elements_per_bin) + i;
			tmp.emplace_back(key);
		}
		elems.emplace_back(std::move(tmp));
	}
	StopClock iff_construct{};
	iff_construct.start();
	ulrich2::interleaved_binary_fuse_filter<> iff(elems);

/*	while (true)
	{
		bool success = true;
		for (int e = 0; e < bins ; ++e)
		{
			std::vector<uint64_t> tmp{};
			for (uint64_t i = 0; i < elements_per_bin; ++i)
			{
				uint64_t key = (e*elements_per_bin) + i;
				tmp.emplace_back(key);
			}
			success = ixf.add_bin_elements(e, tmp);
			if (!success)
			{
				ixf.clear();
				ixf.set_seed();
				break;
			}
			if (e == 2)
				elems=std::move(tmp);
		}

		if (success)
			break;
	}
	*/
	iff_construct.stop();
	std::cout << "Built IFF in " << iff_construct.elapsed() << " seconds" << std::endl;

	StopClock iff_query{};
	iff_query.start();
	typedef ulrich2::interleaved_binary_fuse_filter<>::counting_agent_type< uint64_t > TIFFAgent;
	TIFFAgent iff_count_agent = iff.counting_agent< uint64_t >();
	//auto result = count_agent.bulk_count(readHs);
	auto iff_result = iff_count_agent.bulk_count(elems[2]);
	iff_query.stop();
	//seqan3::debug_stream << iff_result << "\n";

	double fpr {0.0};
	for (int i =0; i < iff_result.size(); ++i)
	{
		if (i == 2)
			continue;
		fpr += (double) iff_result[i] / (double) elements_per_bin;
	}
	fpr /= (double) (bins - 1);
	std::cout << "FPR of the IFF: " << fpr << std::endl;
	double mbytes = (double) iff.bit_size() / (double) 8388608;
	std::cout << "Size of the IFF: " << mbytes << " MBytes" << std::endl;
	std::cout << "Queried " << elems[2].size() <<" keys in IFF in " << iff_query.elapsed() << " seconds" << std::endl;


	StopClock iff_store{};
	iff_store.start();
	std::string filter_file = "test.iff";
	
	std::ofstream               os( filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( iff );  
	iff_store.stop();
	std::cout << "Stored IFF in " << iff_store.elapsed() << " seconds" << std::endl;

}

void construct_and_query_ibf(size_t bins, size_t elements_per_bin)
{
	
	StopClock ibf_construct{};
	ibf_construct.start();
	uint64_t BinSizeBits = get_bin_size(0.0039, 3, elements_per_bin);
	std::cout << "Size of 1 Bin in MByte: " << (double) BinSizeBits / (double) 8388608 << std::endl;
	seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{bins}, seqan3::bin_size{BinSizeBits}, seqan3::hash_function_count{3u}};
	std::vector<uint64_t> elems{};
	for (uint16_t e = 0; e < bins; ++e)
	{
		for (uint64_t i = 0; i < elements_per_bin; ++i)
		{
			uint64_t key = (e*elements_per_bin) + i;
			ibf.emplace(key, seqan3::bin_index{e});
			if (e == 2)
				elems.emplace_back(key);

		}
	}
	ibf_construct.stop();
	std::cout << "Built IBF in " << ibf_construct.elapsed() << " seconds" << std::endl;

	StopClock ibf_query{};
	ibf_query.start();
	typedef seqan3::interleaved_bloom_filter<>::counting_agent_type< uint64_t > TIBFAgent;
	TIBFAgent ibf_count_agent = ibf.counting_agent< uint64_t >();
	auto ibf_result = ibf_count_agent.bulk_count(elems);
	ibf_query.stop();
	//seqan3::debug_stream << ibf_result << '\n';
	double mbytes = (double) ibf.bit_size() / (double) 8388608;
	std::cout << "Size of the IBF: " << mbytes << " MBytes" << std::endl;
	std::cout << "Queried " << elems.size() <<" keys in IBF in " << ibf_query.elapsed() << " seconds" << std::endl;

	StopClock ibf_store{};
	ibf_store.start();
	std::string filter_file = "test.ibf";
	
	std::ofstream               os( filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( ibf );  
	ibf_store.stop();
	std::cout << "Stored IBF in " << ibf_store.elapsed() << " seconds" << std::endl;

}

void load_and_query_ixf(const std::string filter_file, size_t elements, size_t batches)
{
	seqan3::interleaved_xor_filter<> ixf{};
	std::ifstream              is( filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( ixf );

	std::cout << "IXF successfully loaded" << std::endl;
//	auto agent = ixf.membership_agent();
//    auto & result = agent.bulk_contains(2876);
    //seqan3::debug_stream << result << '\n';
	typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;
	TIXFAgent ixf_count_agent = ixf.counting_agent< uint64_t >();
	StopClock ixf_query{};
	ixf_query.start();
	for (uint64_t b = 0; b < batches; ++b)
	{
//		std::cout << "Batch: " << b+1 << std::endl;
		std::vector<uint64_t> tmp{};
		for (uint64_t key = 0; key < elements; ++key)
		{
			tmp.emplace_back(key);
		}
		auto ixf_result = ixf_count_agent.bulk_count(tmp);
	}
	
	ixf_query.stop();
	std::cout << "Queried " << batches * elements <<" keys in IBF in " << ixf_query.elapsed() << " seconds" << std::endl;

}

void bins_to_species(multi_interleaved_xor_filter& mixf, std::vector<std::vector<Species*>>& bins2vec)
{
	for (int i = 0 ; i < mixf.count_single_filter(); ++i)
	{
		std::vector<Species*> tmp(mixf.get(i).bin_count());
		bins2vec[i] = std::move(tmp);
	}

	for (auto & species : mixf.species_vector)
	{
		if (species.filter_index == UINT16_MAX)
			continue;
		auto &filter = bins2vec[species.filter_index];
		for (uint64_t bin = species.first_bin; bin <= species.last_bin; ++bin)
		{
			filter[bin] = &species;
		}

	}
}



int main(int argc, char const **argv)
{

	hixf::build_arguments args{};
	args.bin_file = std::filesystem::path{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/binning.out"};
	args.out_path = "/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/raptor.hixf";
	args.compute_syncmer = false;
	hixf::chopper_build(args);
	
	return 1;
	uint64_t q_p = pow (2, 8) - 1;
	int kmer_size = 16;
	int sync_size = 4;
	int low_sync_offset = 2;
	int up_sync_offset = 12;
	int max_dist = 255;
	int w_min = kmer_size/(kmer_size-sync_size+1) + low_sync_offset > 1 ? kmer_size/(kmer_size-sync_size+1) + low_sync_offset : 1;
    int w_max = kmer_size/(kmer_size-sync_size+1) + up_sync_offset;
    int t_syncmer = (kmer_size-sync_size)/2 + 1;

	//for (int i = 4001 ; i < 5000 ; ++i)
	//{
	//	if(construct_and_query_ixf(4500, 100000) < 100000)
	//construct_and_query_ixf(4363, 100000);
	//construct_and_query_ixf(4364, 110000);
	//		std::cout << "Hier : " << i << std::endl;
	//}

	

	std::string gtdb_root = "/media/jens/INTENSO/gtdb_genomes_reps_r202";
	std::string filename = gtdb_root + "/gtdb_reps_accessions_NCBI_name.tsv";
	std::filesystem::path bin_meta_file = "GTDB_r202_bin_meta.tsv";
	
	//TGenVec genera{};
	std::vector<Species> species;
	load_species_from_file(filename, std::ref(species));
	std::cout << "species loaded" << std::endl;
	//uint64_t max_hashes = 0;
	
	uint64_t min_syncmers = UINT64_MAX;
	uint64_t max_syncmers = 0;
	uint64_t all_syncmers = 0;
	

	//std::vector<std::vector<uint64_t>> genome_hashes{};

	// compute number of bins based on maximum number of elements per bin
	
	multi_interleaved_xor_filter mixf = build_ixf_index(species, "GTDB_r202.ixf", bin_meta_file, gtdb_root, kmer_size, sync_size, t_syncmer);	
	std::cout << "Filter built successfully " << std::endl;

/*
	std::filesystem::path genomes_path = argv[1];
	std::vector<std::vector<uint64_t>> genome_hashes{};
	for (const auto & entry : std::filesystem::directory_iterator(argv[1]))
	{
		std::vector<Seqs> ref_seqs{};
		parse_ref_seqs(ref_seqs, entry.path());
		uint64_t kmers = 0;
		uint64_t strobes = 0;
		std::set<uint64_t> randstrobes{};
		for (const auto & ref_s : ref_seqs)
		{
			std::vector<uint64_t> strobe_hashes = seq_to_randstrobes(kmer_size, w_min, w_max, ref_s.seq, sync_size, t_syncmer, q_p, max_dist);
			randstrobes.insert(strobe_hashes.begin(), strobe_hashes.end());
			strobes += strobe_hashes.size();
			kmers += ref_s.seq.size() - kmer_size + 1;
		}
		float ratio = (float) strobes / (float) kmers;
		std::cout << entry.path() << std::endl;
		std::cout << "kmers: " << kmers << "\t randstrobes: " << strobes << "\t ratio: " << ratio <<std::endl;
		genome_hashes.emplace_back(std::move(std::vector<uint64_t>(randstrobes.begin(), randstrobes.end())));
	}

	seqan3::interleaved_xor_filter<> ixf(genome_hashes);
	double mbytes = (double) ixf.bit_size() / (double) 8388608;
	std::cout << "Size of the IXF: " << mbytes << " MBytes" << std::endl;
*/	
	std::filesystem::path read_file = argv[1];
	seqan3::sequence_file_input fin{read_file};
	
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
}

