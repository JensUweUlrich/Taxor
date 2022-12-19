
#ifndef index_hpp
#define index_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <deque>
#include <tuple>
#include "robin_hood.h"
#include "xxhash.h"
#include <inttypes.h>
#include <seqan3/alphabet/nucleotide/all.hpp>

namespace hashing
{
//uint64_t hash(std::string kmer);

robin_hood::unordered_flat_set<uint64_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t);

}                                        


#endif /* index_hpp */



