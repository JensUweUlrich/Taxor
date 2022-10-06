//
//  index.hpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

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

uint64_t hash(std::string kmer);
static inline uint64_t hash64(uint64_t key, uint64_t mask);


typedef std::vector< std::tuple<uint64_t, unsigned int, int >> mers_vector;
//typedef std::vector< std::tuple<uint64_t, unsigned int, unsigned int>> mers_vector_reduced;
typedef robin_hood::unordered_map< uint64_t, std::tuple<unsigned int, unsigned int >> kmer_lookup;

typedef std::vector< std::pair<std::vector<uint64_t>, bool>> read_hash_vector;

//static inline void make_string_to_hashvalues(std::string &seq, std::vector<uint64_t> &string_hashes, int k, uint64_t kmask);
static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, 
                                   unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  
                                   unsigned int w_start, unsigned int w_end, uint64_t q);


//mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index);
//std::vector<uint64_t> seq_to_randstrobes(int k, int w_min, int w_max, const seqan3::dna5_vector &seq, 
//                                         int s, int t, uint64_t q, int max_dist);

std::vector<uint64_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t);
                                         
read_hash_vector seq_to_randstrobes_read(int k, int w_min, int w_max, const seqan3::dna5_vector &seq, 
                                         int s, int t, uint64_t q, int max_dist);


void process_flat_vector(mers_vector &flat_vector, uint64_t &unique_elements);
unsigned int index_vector(mers_vector  &mers_vector, kmer_lookup &mers_index, float f);

struct hit {
    int query_s;
    int query_e;
    int ref_s;
    int ref_e;
    bool is_rc = false;
};

struct nam {
    int nam_id;
    int query_s;
    int query_e;
    int query_prev_hit_startpos;
    int ref_s;
    int ref_e;
    int ref_prev_hit_startpos;
    int n_hits = 0;
    int ref_id;
    float score;
//    unsigned int previous_query_start;
//    unsigned int previous_ref_start;
    bool is_rc = false;
};


#endif /* index_hpp */


