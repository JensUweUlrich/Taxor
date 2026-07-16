
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <syncmer.hpp>

#include "adjust_seed.hpp"
#include "compute_hashes.hpp"
#include "dna4_traits.hpp"

namespace hixf
{

/*!\file compute_hashes.cpp
 * \brief Implements hixf::compute_hashes: reading reference sequences and turning them into k-mer/syncmer hashes.
 */

using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

/*!\brief Splits a sequence into substrings at stretches of 'N' characters, discarding the 'N' runs themselves.
 * \param seq The reference sequence as a string.
 * \param seqlen The length of `seq`; used to bound the final chunk when the last non-N stretch is not
 *               terminated by another 'N'.
 * \return A vector of the substrings found between (and excluding) runs of 'N'.
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
			 std::string s = seq.substr(start, seqlen - start);
			 splittedStrings.push_back(s);
			 break;
		 }
		 std::string s = seq.substr(start, end - start);

		 splittedStrings.push_back(s);
	 }
	 return splittedStrings;
 }

/*!\brief Reads all sequences from a reference file and appends the cleaned sequences to `refs`.
 * \param refs Output vector that the parsed sequences are appended to.
 * \param reference_file Path to the (optionally compressed) sequence file to read.
 * \details
 *
 * Each record is converted to a plain string, any stretches of 'N' are removed via cutOutNNNs(), and the
 * remaining pieces are concatenated back together before being stored as a seqan3::dna5_vector.
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

/*!\brief Computes the set of distinct hash values (syncmers or minimisers) for all sequences of a user bin.
 * \param hashes Output set; computed hashes are inserted here (not cleared beforehand).
 * \param arguments The build arguments; `compute_syncmer` selects the hashing strategy, `scaling` optionally
 *                  subsamples the resulting hash set.
 * \param record The chopper pack record whose reference sequence file(s) are read and hashed.
 * \details
 *
 * If `arguments.compute_syncmer` is set, sequences are read, cleaned of 'N' runs, and hashed via
 * hashing::seq_to_syncmers(); otherwise a plain/gapped minimiser hash (seqan3::views::minimiser_hash) is used.
 * If `arguments.scaling > 1`, only a `1/scaling` fraction of hashes are kept: each hash value is re-hashed with
 * wyhash and only kept if the re-hash falls below `UINT64_MAX / scaling`. This is a FracMinHash-style
 * downsampling scheme that keeps (approximately) the same subset of hash values across different sequences,
 * which is required for the downsampled comparisons to be meaningful.
 */
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
