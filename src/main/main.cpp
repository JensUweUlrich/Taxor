

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include <cereal/archives/binary.hpp>

#include <zlib.h>
#include "StopClock.hpp"
#include "interleaved_binary_fuse_filter.hpp"

using namespace seqan3::literals;

//std::string gtdb_root = "D:\\gtdb_genomes_reps_r202";

struct species
{
	std::string name;
	std::string organism_name;
	std::string accession_id;
};

struct Seqs {
	
	std::string       seqid;
	seqan3::dna4_vector seq;
	
};

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
}

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

 /**
   * split underlying sequence based on stretches of N's
   * @seq      : reference sequence as a string
   * @seqlen   : length of the sequence
   * @return   : vector of subsequences, that results from splitting the original sequence
   */
 std::vector<std::string> cutOutNNNs(std::string& seq, uint64_t seqlen)
 {
	 std::vector<std::string> splittedStrings;
	 size_t start = 0;
	 size_t end = 0;

	 while ((start = seq.find_first_not_of("N", end)) != std::string::npos)
	 {
		 end = seq.find("N", start);
		 if (end > seqlen)
		 {
			 std::string s = seq.substr(start, seqlen - start - 1);
			 splittedStrings.push_back(s);
			 break;
		 }
		 std::string s = seq.substr(start, end - start);

		 splittedStrings.push_back(s);
	 }
	 return splittedStrings;
 }


 /**
	   read reference sequences from files and store them in reference sequence queue
	   @queue_refs:    thread safe queue that stores reference sequences to put in the ibf
	   @config:        settings used to build the IBF (reference file path, number of references, kmer size, fragment length)
	   @stats:         statistics (like runtime for specific tasks) for building the ibf
	   @throws:        FileParserException
   */

 void parse_ref_seqs(std::vector< Seqs >& queue_refs, const std::filesystem::path& reference_file)
 {
	 //seqan::SeqFileIn seqFileIn;

	 seqan3::sequence_file_input fin{reference_file};


	 for ( auto& [seq, header, qual] : fin)
	 {
		 const std::string seqid = get_seqid( header );

		 int counter = 1;
		 auto charstring = seq | seqan3::views::to_char;
		 std::string s = std::string(charstring.begin(), charstring.end());

	 	 // remove all Ns from the sequence
		 std::stringstream buf;
		 for (std::string subseq : cutOutNNNs(s, s.size()))
		 {
			 buf << subseq;
		 }
		 s =  buf.str();
		 auto r = s | seqan3::views::char_to<seqan3::dna4>;
		 seqan3::dna4_vector newseq(r.begin(), r.end());

		 queue_refs.emplace_back(Seqs{ seqid, newseq });
					
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
	/*
	std::string filename = gtdb_root + "\\gtdb_reps_accessions_NCBI_name.tsv";
	std::map<std::string, std::vector<species>> genera{};
	getGenomeRepsByGenus(filename, std::ref(genera));

	for (std::pair<std::string, std::vector<species>> genus : genera)
	{
		auto filepath = std::filesystem::path(gtdb_root) / "genera" / genus.first;
		filepath.replace_extension("fna.gz");
		FILE* outfile = fopen(seqan::toCString(filepath.string()), "wb");
		std::cout << filepath.string() << std::endl;
		std::string outy{};
		for (species sp : genus.second)
		{
			auto inpath = std::filesystem::path(gtdb_root) / "gtdb_genomes_reps_r202";
			if (sp.accession_id.substr(0, 3).compare("GCA") == 0)
			{
				inpath /= "GCA";
			}
			else
			{
				inpath /= "GCF";
			}
			std::string infilename = sp.accession_id + "_genomic.fna.gz";
			inpath = inpath / sp.accession_id.substr(4, 3) / sp.accession_id.substr(7, 3) / sp.accession_id.substr(10, 3) / infilename;
			
			std::string fileData;
			if (!loadBinaryFile(inpath.string(), fileData)) {
				printf("Error loading input file.");
				return 0;
			}

			// Uncompress the file data  
			std::string data;
			if (!gzipInflate(fileData, data)) {
				printf("Error decompressing file.");
				return 0;
			}


			std::stringstream fdata(data);
			std::string line;
			while (getline(fdata, line)) {
				if (line.substr(0, 1).compare(">") == 0)
				{
					std::cout << line << std::endl;
				}
			}
			*/

	std::vector<std::vector<uint64_t>> elems{};
/*	std::vector<Seqs> ref_seqs{};
	parse_ref_seqs(ref_seqs, std::filesystem::path(argv[1]));
	auto refHashes = seqan3::views::kmer_hash(ref_seqs[0].seq, seqan3::ungapped{15});
	std::set<uint64_t> kmerHashSet(refHashes.begin(), refHashes.end());
	std::vector<uint64_t> v1 = std::vector(kmerHashSet.begin(), kmerHashSet.end());
//	std::vector<uint64_t> v1{127, 2876, 2875, 457, 458, 875, 7654, 7566, 5876};
	elems.emplace_back(std::move(v1));

	ref_seqs.clear();
	kmerHashSet.clear();
	parse_ref_seqs(ref_seqs, std::filesystem::path(argv[2]));
	refHashes = seqan3::views::kmer_hash(ref_seqs[0].seq, seqan3::ungapped{15});
	kmerHashSet = std::set<uint64_t>(refHashes.begin(), refHashes.end());
	std::vector<uint64_t> v2 = std::vector(kmerHashSet.begin(), kmerHashSet.end());
//	std::vector<uint64_t> v2{7465, 2876, 8576, 856, 6586, 5876, 7654, 5785};
	elems.emplace_back(std::move(v2));
*/
	size_t nr_elems = atoi(argv[1]);
	size_t nr_bins = atoi(argv[2]);

//	construct_and_query_iff(nr_bins, nr_elems);
	construct_and_query_ibf(nr_bins, nr_elems);
	construct_and_query_ixf(nr_bins, nr_elems);

//	load_and_query_ixf("test_mum3.ixf", nr_elems, nr_bins);

/*	std::vector<uint64_t> v1{127, 2876, 2875, 457, 458, 875, 7654, 7566, 5876};
	elems.emplace_back(std::move(v1));
	std::vector<uint64_t> v2{7465, 2876, 8576, 856, 6586, 5876, 7654, 5785};
	elems.emplace_back(std::move(v2));
*/

//	std::vector<uint64_t> readHs{7465, 2876, 8576, 856, 6586, 5876, 7654};
	
	
	/*for (int e = 0; e < nr_bins ; ++e)
	{
		std::vector<uint64_t> tmp{};
		for (uint64_t i = 0; i < nr_elems; ++i)
		{
			uint64_t key = (e*nr_elems) + i;
			tmp.emplace_back(key);
		}
		elems.emplace_back(std::move(tmp));
	}
	*/
	
/*	ulrich2::interleaved_binary_fuse_filter<> iff(elems);
	auto agent = iff.membership_agent();
    auto & result = agent.bulk_contains(2876);
    seqan3::debug_stream << result << '\n';
*/
//	ulrich::interleaved_xor_filter<> ixf(elems);
/*
	typedef ulrich2::interleaved_binary_fuse_filter<>::counting_agent_type< uint64_t > TIFFAgent;
	TIFFAgent iff_count_agent = iff.counting_agent< uint64_t >();
	auto result = iff_count_agent.bulk_count(readHs);
	seqan3::debug_stream << result << '\n';
*/
	
	
	


/*
	seqan3::interleaved_bloom_filter ibf2{};
	std::ifstream              is( filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( ibf2 );

	TAgent count_agent2 = ibf2.counting_agent< uint64_t >();
	auto result2 = count_agent2.bulk_count(elems[2]);
	seqan3::debug_stream << result2 << '\n';
*/

/*	ulrich::interleaved_xor_filter<> ixf2{};
	std::ifstream              is( filter_file, std::ios::binary );
    cereal::BinaryInputArchive archive( is );
    archive( ixf2 );

	TAgent count_agent2 = ixf2.counting_agent< uint64_t >();
	auto result2 = count_agent2.bulk_count(elems[2]);
	seqan3::debug_stream << result2 << '\n';
*/
	return 0;
}

