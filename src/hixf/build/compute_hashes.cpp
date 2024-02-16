
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/core/debug_stream.hpp>
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

void compute_hashes(ankerl::unordered_dense::set<size_t> & hashes,
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
                //std::cout << "Before compute syncmers" << std::endl << std::flush;
                ankerl::unordered_dense::set<size_t> tmp = hashing::seq_to_syncmers(arguments.kmer_size, seq, arguments.syncmer_size, arguments.t_syncmer);
                if (arguments.scaling > 1)
                {
                    for (auto &hash : tmp)
                    {
                        uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                        if (double(v) <= double(UINT64_MAX) / double(arguments.scaling))
                        {
                            hashes.insert(hash);
                        }
                    }
                }
                else
                {
			        hashes.insert(std::make_move_iterator(tmp.begin()), std::make_move_iterator(tmp.end()));
                }
                //std::cout << "After compute syncmers" << std::endl << std::flush;
                //hashes.insert(hashing::seq_to_syncmers(arguments.kmer_size, seq, arguments.syncmer_size, arguments.t_syncmer));
                
		    }
            
        }
    }
    else
    {
        //using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        for (auto const & filename : record.filenames)
        {
            for (auto && [seq] : sequence_file_t{filename})
            {
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape.count())}))
                    {
                        if (arguments.scaling > 1)
                        {
                            uint64_t v = ankerl::unordered_dense::detail::wyhash::hash(hash);
                            if (double(v) <= double(UINT64_MAX) / double(arguments.scaling))
                            {
                                hashes.insert(hash);
                            }
                        }
                        else
                        {
                            hashes.insert(hash);
                        }
                    }
            }
        }
    }
}

}
