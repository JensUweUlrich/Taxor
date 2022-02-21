

#ifndef build_hpp
#define build_hpp

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

struct species
{
	std::string name;
	std::string organism_name;
	std::string accession_id;
};

typedef std::map<std::string, std::vector<species>> TGenMap;


struct Seqs {
	
	std::string       seqid;
	seqan3::dna5_vector seq;
	
};

typedef std::tuple<std::string, uint64_t, uint64_t> TSpecBin;

inline std::string get_seqid( std::string header )
{
    return header.substr( 0, header.find( ' ' ) );
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
		 auto r = s | seqan3::views::char_to<seqan3::dna5>;
		 seqan3::dna5_vector newseq(r.begin(), r.end());

		 queue_refs.emplace_back(Seqs{ seqid, newseq });
					
	 }
			 
 }

uint64_t compute_bin_number(TGenMap& genera, const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer)
{
	uint64_t bin_nr = 0;
	uint64_t species_counter = 0;
	for (std::pair<std::string, std::vector<species>> genus : genera)
	{
		if (genus.first.compare("Pseudomonas") != 0)
			continue;
//		auto filepath = std::filesystem::path(gtdb_root) / "genera" / genus.first;
//		filepath.replace_extension("fna.gz");
//		FILE* outfile = fopen(seqan::toCString(filepath.string()), "wb");
//		std::cout << filepath.string() << std::endl;
		std::string outy{};
		for (species sp : genus.second)
		{
//			std::cout << sp.name << std::endl;
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

			std::vector<Seqs> ref_seqs{};
			parse_ref_seqs(ref_seqs, inpath);
			
			
			std::set<uint64_t> randstrobes{};
			for (const auto & ref_s : ref_seqs)
			{
				std::vector<uint64_t> strobe_hashes = seq_to_syncmers(kmer_size, ref_s.seq, sync_size, t_syncmer);
				randstrobes.insert(strobe_hashes.begin(), strobe_hashes.end());
			}
			
//			genome_hashes.emplace_back(std::move(std::vector<uint64_t>(randstrobes.begin(), randstrobes.end())));

			bin_nr += (randstrobes.size() / 50000) + 1;

			if (++species_counter % 1000 == 0)
			{
				std::cout << species_counter << std::endl;
			}
			
		}
	}
	return bin_nr;
}

void write_meta_file(std::vector<TSpecBin>& bin_meta, const std::filesystem::path& file)
{
    std::ofstream outfile(file.string());
    for (const auto & sp_bin : bin_meta)
	{
		outfile << std::get<0>(sp_bin) << "\t" << std::get<1>(sp_bin) << "\t" << std::get<2>(sp_bin) << std::endl;
	}
    outfile.close();
}

void build_ixf_index(TGenMap& genera, const std::string& filter_file, const std::filesystem::path& bin_meta_file,
                     const std::string& gtdb_root, int kmer_size, int sync_size, int t_syncmer)
{
    uint64_t bin_nr = compute_bin_number(std::ref(genera), gtdb_root,kmer_size, sync_size, t_syncmer);

	std::cout << bin_nr << std::endl;

	seqan3::interleaved_xor_filter<> ixf(bin_nr, 50000);
	bool success = true;
	bin_nr = 0;

	std::vector<TSpecBin> bin_meta{};
	while (true)
	{
		uint64_t species_counter = 0;
		
		for (std::pair<std::string, std::vector<species>> genus : genera)
		{
		if (genus.first.compare("Pseudomonas") != 0)
			continue;

			std::string outy{};
			for (species sp : genus.second)
			{
				std::cout << sp.name << std::endl;
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

				std::vector<Seqs> ref_seqs{};
				parse_ref_seqs(ref_seqs, inpath);
			
			
				std::set<uint64_t> randstrobes{};
				for (const auto & ref_s : ref_seqs)
				{
					std::vector<uint64_t> strobe_hashes = seq_to_syncmers(kmer_size, ref_s.seq, sync_size, t_syncmer);
					randstrobes.insert(strobe_hashes.begin(), strobe_hashes.end());
				}

				std::vector<uint64_t> syncmers(randstrobes.begin(), randstrobes.end());
				uint64_t start_bin = bin_nr;
				for (int j = 0; j < (randstrobes.size() / 50000) + 1 ; j++)
				{
					std::vector<uint64_t>::iterator first = syncmers.begin() + (j*50000);
					std::vector<uint64_t>::iterator last = syncmers.begin() + ((j+1)*50000);
					if (j == randstrobes.size() / 50000)
						last = syncmers.end();
					std::vector<uint64_t> subvec = std::vector<uint64_t>(first, last);
					success = ixf.add_bin_elements(bin_nr++, subvec);
//					std::cout << sp.name << "\t" << syncmers.size() << "\t" << bin_nr-1 << "/" << j << "\t" << subvec.size() << "\t" << success << std::endl;
					if (!success)
						break;
				}
				bin_meta.emplace_back(std::move(std::make_tuple(sp.name ,start_bin, bin_nr-1)));
				if (++species_counter % 1000 == 0)
				{
					std::cout << species_counter << std::endl;
				}
					
				if (!success)
					break;
			}

			if (!success)
			{
				ixf.clear();
				ixf.set_seed();
				bin_nr = 0;
				species_counter = 0;
				bin_meta.clear();
				break;
			}
		}

		if (success)
			break;
				
	}

	write_meta_file(bin_meta, bin_meta_file);
    std::ofstream os( filter_file, std::ios::binary );
    cereal::BinaryOutputArchive archive( os );
    archive( ixf ); 
}

#endif