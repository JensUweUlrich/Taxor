

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

using namespace seqan3::literals;

//std::string gtdb_root = "D:\\gtdb_genomes_reps_r202";



 void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines) {
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
/*
 bool gzipInflate(const std::string& compressedBytes, std::string& uncompressedBytes) {
	 if (compressedBytes.size() == 0) {
		 uncompressedBytes = compressedBytes;
		 return true;
	 }

	 uncompressedBytes.clear();

	 unsigned full_length = compressedBytes.size();
	 unsigned half_length = compressedBytes.size() / 2;

	 unsigned uncompLength = full_length;
	 char* uncomp = (char*)calloc(sizeof(char), uncompLength);

	 z_stream strm;
	 strm.next_in = (Bytef*)compressedBytes.c_str();
	 strm.avail_in = compressedBytes.size();
	 strm.total_out = 0;
	 strm.zalloc = Z_NULL;
	 strm.zfree = Z_NULL;

	 bool done = false;

	 if (inflateInit2(&strm, (16 + MAX_WBITS)) != Z_OK) {
		 free(uncomp);
		 return false;
	 }

	 while (!done) {
		 // If our output buffer is too small  
		 if (strm.total_out >= uncompLength) {
			 // Increase size of output buffer  
			 char* uncomp2 = (char*)calloc(sizeof(char), uncompLength + half_length);
			 memcpy(uncomp2, uncomp, uncompLength);
			 uncompLength += half_length;
			 free(uncomp);
			 uncomp = uncomp2;
		 }

		 strm.next_out = (Bytef*)(uncomp + strm.total_out);
		 strm.avail_out = uncompLength - strm.total_out;

		 // Inflate another chunk.  
		 int err = inflate(&strm, Z_SYNC_FLUSH);
		 if (err == Z_STREAM_END) done = true;
		 else if (err != Z_OK) {
			 break;
		 }
	 }

	 if (inflateEnd(&strm) != Z_OK) {
		 free(uncomp);
		 return false;
	 }

	 for (size_t i = 0; i < strm.total_out; ++i) {
		 uncompressedBytes += uncomp[i];
	 }
	 free(uncomp);
	 return true;
 }
*/
 /* Reads a file into memory. */
 /*
 bool loadBinaryFile(const std::string& filename, std::string& contents) {
	 // Open the gzip file in binary mode  
	 FILE* f = fopen(filename.c_str(), "rb");
	 if (f == NULL)
		 return false;

	 // Clear existing bytes in output vector  
	 contents.clear();

	 // Read all the bytes in the file  
	 int c = fgetc(f);
	 while (c != EOF) {
		 contents += (char)c;
		 c = fgetc(f);
	 }
	 fclose(f);

	 return true;
 }
*/
 /** Compress a STL string using zlib with given compression level and return
 * the binary data. */
 /*
 std::string compress_string(const std::string& str, int compressionlevel = Z_BEST_COMPRESSION)
 {
	 z_stream zs;                        // z_stream is zlib's control structure
	 memset(&zs, 0, sizeof(zs));

	 if (deflateInit(&zs, compressionlevel) != Z_OK)
		 throw(std::runtime_error("deflateInit failed while compressing."));

	 zs.next_in = (Bytef*)str.data();
	 zs.avail_in = str.size();           // set the z_stream's input

	 int ret;
	 char outbuffer[32768];
	 std::string outstring;

	 // retrieve the compressed bytes blockwise
	 do {
		 zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
		 zs.avail_out = sizeof(outbuffer);

		 ret = deflate(&zs, Z_FINISH);

		 if (outstring.size() < zs.total_out) {
			 // append the block to the output string
			 outstring.append(outbuffer, zs.total_out - outstring.size());
		 }
	 } while (ret == Z_OK);

	 deflateEnd(&zs);

	 if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
		 std::ostringstream oss;
		 oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
		 throw(std::runtime_error(oss.str()));
	 }

	 return outstring;
 }
*/

 void getGenomeRepsByGenus(std::string& filename, std::map<std::string, std::vector<species>>& genera)
 {
	 std::vector<std::vector<std::string> > gtdb_reps_ncbi{};
	 read_tsv(filename, gtdb_reps_ncbi);
	 std::map<std::string, std::vector<species>>::iterator gen_it;
	 for (std::vector<std::string> line : gtdb_reps_ncbi)
	 {
		 if (line[0].compare("organism_name") == 0)
			 continue;

		 species sp{};
		 sp.accession_id = line[0];
		 sp.organism_name = line[1];
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
			 std::vector<species> sp_vec{};
			 sp_vec.emplace_back(sp);
			 genera.insert(std::make_pair(genus, std::move(sp_vec)));
		 }
	 }
 }

 

void construct_and_query_ixf(size_t bins, size_t elements_per_bin)
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
			continue;
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



int main(int argc, char const **argv)
{
	uint64_t q_p = pow (2, 8) - 1;
	int kmer_size = 16;
	int sync_size = 4;
	int low_sync_offset = 2;
	int up_sync_offset = 12;
	int max_dist = 255;
	int w_min = kmer_size/(kmer_size-sync_size+1) + low_sync_offset > 1 ? kmer_size/(kmer_size-sync_size+1) + low_sync_offset : 1;
    int w_max = kmer_size/(kmer_size-sync_size+1) + up_sync_offset;
    int t_syncmer = (kmer_size-sync_size)/2 + 1;



	std::string gtdb_root = "/media/jens/INTENSO/gtdb_genomes_reps_r202";
	std::string filename = gtdb_root + "/gtdb_reps_accessions_NCBI_name.tsv";
	std::filesystem::path bin_meta_file = "GTDB_r202_bin_meta.tsv";
	
	TGenMap genera{};
	getGenomeRepsByGenus(filename, std::ref(genera));

	//uint64_t max_hashes = 0;
	
	uint64_t min_syncmers = UINT64_MAX;
	uint64_t max_syncmers = 0;
	uint64_t all_syncmers = 0;
	

	//std::vector<std::vector<uint64_t>> genome_hashes{};

	// compute number of bins based on maximum number of elements per bin
	
	build_ixf_index(genera, "GTDB_r202.ixf", bin_meta_file, gtdb_root, kmer_size, sync_size, t_syncmer);	
	

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
	
	std::filesystem::path read_file = argv[2];
	seqan3::sequence_file_input fin{read_file};
	typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;

	double ratio_sum = 0.0;
	uint64_t nr_reads = 0;
	std::vector<double> sensitivity(ixf.bin_count());
	std::fill(sensitivity.begin(), sensitivity.end(), 0.0);
	std::vector<uint64_t> classified(ixf.bin_count());
	std::fill(classified.begin(), classified.end(), 0);
	int count = 0;
	for ( const auto& [seq, header, qual] : fin)
	{

		if (seq.size() < 1000)
			continue;

		read_hash_vector rhv = seq_to_randstrobes_read(kmer_size, w_min, w_max, seq, sync_size, t_syncmer, q_p, max_dist);
		uint64_t strobe_nr = rhv[0].first.size();
//		TInterval ci = calculateCI(0.1, kmer_size, strobe_nr, 0.95);
//		uint64_t threshold = strobe_nr - ci.second;
		std::vector<uint64_t> max_found(ixf.bin_count());
		std::fill(max_found.begin(), max_found.end(), 0);
		for (const auto & p : rhv)
		{
			TIXFAgent ixf_count_agent = ixf.counting_agent< uint64_t >();
			auto ixf_result = ixf_count_agent.bulk_count(p.first);
			for (int i = 0; i < max_found.size(); ++i)
			{
				if (ixf_result[i] > max_found[i])
					max_found[i] = ixf_result[i];
			}
		}
		double max_ratio = 0.0;
		uint64_t max_index = UINT64_MAX;
		std::vector<std::pair<uint64_t, double>> potential_indexes{};
		for (int i = 0; i < max_found.size(); ++i)
		{
			double ratio = (double) max_found[i] / (double) strobe_nr;
			sensitivity[i] += ratio;

			// Why does 0.05 work so well?
			if ( ratio > 0.05)
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
				if (++spec_count > 1)
					std::cout << header << std::endl;
			}
		}
		
		nr_reads++;
	}

	for (int i=0; i < sensitivity.size(); ++i)
	{
		sensitivity[i] /= (double) nr_reads;
	}

	seqan3::debug_stream << sensitivity << "\n";
	seqan3::debug_stream << classified << "\n";
	seqan3::debug_stream << nr_reads << "\n";
*/
	return 0;
}

