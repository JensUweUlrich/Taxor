// Copyright (c) 2021- Kristoffer Sahlin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// modified by Jens-Uwe Ulrich

#include "syncmer.hpp"
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>
#include <climits>
#include <inttypes.h>


/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:
/*
uint64_t hash(std::string kmer)
{
    unsigned long hash = 5381;
    int c;
    for (std::string::size_type i=0; i< kmer.length(); i++) {
        c = kmer[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
/*    }
    return hash;
}
*/

namespace hashing
{

static unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table


static inline uint64_t syncmer_kmer_hash(uint64_t packed) {
    //return robin_hood::hash_int(packed);
    return ankerl::unordered_dense::detail::wyhash::hash(packed);
    //return XXH64(&packed, sizeof(uint64_t), 0);
}


static inline void make_string_to_hashvalues_open_syncmers_canonical(const seqan3::dna5_vector &seq, 
                                                                     ankerl::unordered_dense::set<size_t> &string_hashes, 
                                                                     const size_t k, 
                                                                     const size_t s, 
                                                                     const size_t t) {
    std::vector<unsigned int> pos_to_seq_coordinate;
    const uint64_t kmask = (1ULL << 2*k) - 1;
    const uint64_t smask = (1ULL << 2*s) - 1;
    const uint64_t kshift = (k - 1) * 2;
    const uint64_t sshift = (s - 1) * 2;
    std::deque<uint64_t> qs;  // s-mer hashes
    uint64_t qs_min_val = UINT64_MAX;
    size_t qs_min_pos = -1;

    size_t l = 0;
    uint64_t xk[] = {0, 0};
    uint64_t xs[] = {0, 0};
    for (int i = 0; i < seq.size(); i++) {
        // tranlate A,C,T,G,a,c,t,g in 0,1,2,3
        int c = seq_nt4_table[(uint8_t) seq[i].to_char()];
        if (c < 4) { // not an "N" base
            xk[0] = (xk[0] << 2 | c) & kmask;                  // forward strand
            xk[1] = xk[1] >> 2 | (uint64_t)(3 - c) << kshift;  // reverse strand
            xs[0] = (xs[0] << 2 | c) & smask;                  // forward strand
            xs[1] = xs[1] >> 2 | (uint64_t)(3 - c) << sshift;  // reverse strand
            if (++l < s) {
                continue;
            }
            // we find an s-mer
            uint64_t ys = std::min(xs[0], xs[1]);
            uint64_t hash_s = ys;
            qs.push_back(hash_s);
            // not enough hashes in the queue, yet
            if (qs.size() < k - s + 1) {
                continue;
            }
            if (qs.size() == k - s + 1) { // We are at the last s-mer within the first k-mer, need to decide if we add it
                for (size_t j = 0; j < qs.size(); j++) {
                    if (qs[j] < qs_min_val) {
                        qs_min_val = qs[j];
                        qs_min_pos = i - k + j + 1;
                    }
                }
            }
            else {
                // update queue and current minimum and position
                qs.pop_front();

                if (qs_min_pos == i - k) { // we popped the previous minimizer, find new brute force
                    qs_min_val = UINT64_MAX;
                    qs_min_pos = i - s + 1;
                    for (int j = qs.size() - 1; j >= 0; j--) { //Iterate in reverse to choose the rightmost minimizer in a window
                        if (qs[j] < qs_min_val) {
                            qs_min_val = qs[j];
                            qs_min_pos = i - k + j + 1;
                        }
                    }
                } else if (hash_s < qs_min_val) { // the new value added to queue is the new minimum
                    qs_min_val = hash_s;
                    qs_min_pos = i - s + 1;
                }
            }
            if (qs_min_pos == i - k + t) { // occurs at t:th position in k-mer
                
                uint64_t yk = std::min(xk[0], xk[1]);
                string_hashes.insert(syncmer_kmer_hash(yk));
            }
        } else {
            // if there is an "N", restart
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
            qs.clear();
        }
    }
}

ankerl::unordered_dense::set<size_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t)
{
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    ankerl::unordered_dense::set<size_t> string_hashes{};

    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, k, s, t);
    return std::move(string_hashes);
    
}



}