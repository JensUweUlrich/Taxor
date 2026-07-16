
#ifndef index_hpp
#define index_hpp

#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
//#include "xxhash.h"
#include <inttypes.h>
#include <seqan3/alphabet/nucleotide/all.hpp>

/*!\file syncmer.hpp
 * \brief Declares the public interface for computing "open syncmers" from a nucleotide sequence.
 *
 * Open syncmers are a k-mer sampling (sketching) scheme: a k-mer is selected ("is a syncmer") if the
 * lexicographically/numerically smallest of its overlapping s-mers (s < k) occurs at a fixed offset \c t
 * within the k-mer. This yields a deterministic, locally-defined subset of k-mers that is used downstream
 * (instead of all k-mers) to build and query the Hierarchical Interleaved XOR Filter index, reducing the
 * number of hashes that need to be stored/compared while preserving good coverage properties.
 */
namespace hashing
{
//uint64_t hash(std::string kmer);

//robin_hood::unordered_flat_set<uint64_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t);
/*!\brief Computes the set of canonical open-syncmer hash values for a nucleotide sequence.
 * \param k   Length of the k-mers considered for syncmer selection.
 * \param seq The input nucleotide sequence (seqan3::dna5_vector).
 * \param s   Length of the internal s-mers (s < k) used to decide whether a k-mer is a syncmer.
 * \param t   Offset (0-based position) within the k-mer at which the minimal s-mer must start for the
 *            k-mer to be selected as an open syncmer.
 * \return An unordered set of hash values (size_t), one for each canonical k-mer selected as a syncmer.
 */
ankerl::unordered_dense::set<size_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t);

}


#endif /* index_hpp */



