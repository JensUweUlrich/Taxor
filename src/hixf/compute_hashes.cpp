
#include <seqan3/search/views/minimiser_hash.hpp>
#include <syncmer.hpp>

#include "adjust_seed.hpp"
#include "compute_hashes.hpp"
#include "dna4_traits.hpp"

namespace hixf
{

using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

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
 void parse_ref_seqs(std::vector< seqan3::dna5_vector >& refs, const std::filesystem::path& reference_file)
 {

	for (auto && [seq] : sequence_file_t{reference_file})
	 {

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

		 refs.emplace_back(newseq);
					
	 }
			 
 }

void compute_hashes(robin_hood::unordered_flat_set<size_t> & hashes,
                   build_arguments const & arguments,
                   chopper_pack_record const & record)
{
    if (arguments.compute_syncmer) 
    {
        for (auto const & filename : record.filenames)
        {
            std::vector<seqan3::dna5_vector> refs{};
            parse_ref_seqs(refs, filename);
            for (const auto & seq : refs)
	    	{
			    std::vector<uint64_t> strobe_hashes = hashing::seq_to_syncmers(arguments.kmer_size, seq, arguments.syncmer_size, arguments.t_syncmer);
                hashes.insert(std::make_move_iterator(strobe_hashes.begin()), std::make_move_iterator(strobe_hashes.end()));
		    }
        }
    }
    else if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : record.filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                hashes.insert(minimiser_value);
        }
    }
    else
    {
        //using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        for (auto const & filename : record.filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape.count())}))
                    hashes.insert(hash);
    }
}

}
