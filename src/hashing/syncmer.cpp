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
//#include <inttypes.h>

/*!\file syncmer.cpp
 * \brief Implements the "open syncmer" k-mer sampling scheme.
 *
 * A k-mer is emitted as a (canonical) syncmer if the minimal of its overlapping s-mers (s < k) starts at a
 * fixed offset \c t within the k-mer. The implementation packs nucleotides into 2-bit-per-base integers and
 * maintains both the forward-strand and reverse-complement encoding side by side so that the canonical
 * (strand-independent) representation of each k-mer/s-mer can be used, which is required for symmetric
 * matching between reference and query sequences regardless of which strand they were sequenced from.
 */

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

/*!\brief Lookup table mapping an ASCII nucleotide character to a 2-bit code (A=0, C=1, T=2, G=3).
 *
 * Any character that is not one of A/C/G/T/a/c/g/t (including 'N'/'n' and other ambiguity codes) maps to 4,
 * which callers use as a sentinel to detect and skip/reset on non-ACGT bases.
 */
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


/*!\brief Hashes a 2-bit-packed k-mer integer into a final hash value.
 * \param packed The 2-bit-per-base packed representation of a k-mer (or s-mer).
 * \return A 64-bit hash of \p packed, computed with the wyhash implementation from ankerl::unordered_dense.
 */
static inline uint64_t syncmer_kmer_hash(uint64_t packed) {
    //return robin_hood::hash_int(packed);
    return ankerl::unordered_dense::detail::wyhash::hash(packed);
    //return XXH64(&packed, sizeof(uint64_t), 0);
}


/*!\brief Scans a sequence and inserts the hash of every canonical open syncmer k-mer into \p string_hashes.
 * \param seq           The input nucleotide sequence.
 * \param string_hashes Output set that hashes of selected syncmer k-mers are inserted into.
 * \param k             Length of the k-mers considered for syncmer selection.
 * \param s             Length of the internal s-mers (s < k) used to determine the syncmer offset.
 * \param t             Required 0-based offset of the minimal s-mer within the k-mer for the k-mer to be
 *                       selected as an open syncmer (i.e. a k-mer is emitted iff its minimal s-mer starts
 *                       exactly at position \c t within the k-mer).
 *
 * The sequence is scanned once, maintaining 2-bit-packed forward- and reverse-complement encodings of the
 * current k-mer and s-mer (\c xk / \c xs) so that canonical (strand-independent) hash values can be derived
 * without re-reading previously seen bases. A sliding window (\c qs, a deque of the last k-s+1 s-mer hashes
 * that overlap the current k-mer) is maintained together with the position and value of its current minimum
 * (\c qs_min_val / \c qs_min_pos), so the minimal s-mer of each k-mer window can be tracked incrementally
 * rather than recomputed from scratch for every new k-mer.
 */
static inline void make_string_to_hashvalues_open_syncmers_canonical(const seqan3::dna5_vector &seq,
                                                                     ankerl::unordered_dense::set<size_t> &string_hashes,
                                                                     const size_t k,
                                                                     const size_t s,
                                                                     const size_t t) {
    std::vector<unsigned int> pos_to_seq_coordinate;
    // masks/shifts to maintain 2-bit-packed forward (left-shift) and reverse-complement (right-shift)
    // encodings of the current k-mer/s-mer window as a fixed-width sliding integer
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

                // "i - k" is the sequence position of the s-mer that just fell out of the window;
                // if that was where the current minimum lived, it is no longer valid and the whole
                // window must be rescanned (this is the O(k) worst case of an otherwise O(1) update)
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
            // the k-mer ending at position i is a syncmer iff its minimal s-mer starts exactly at
            // offset t within it; only then is the (canonical) k-mer hash emitted
            if (qs_min_pos == i - k + t) { // occurs at t:th position in k-mer

                uint64_t yk = std::min(xk[0], xk[1]);
                string_hashes.insert(syncmer_kmer_hash(yk));
            }
        } else {
            // if there is an "N", restart
            // an ambiguous base invalidates any k-mer/s-mer spanning it, so the packed encodings, the
            // valid-base run length l, and the sliding-window state are all reset; scanning resumes as
            // if starting a fresh sequence right after this position
            qs_min_val = UINT64_MAX;
            qs_min_pos = -1;
            l = xs[0] = xs[1] = xk[0] = xk[1] = 0;
            qs.clear();
        }
    }
}

/*!\brief Computes the set of canonical open-syncmer hash values for a nucleotide sequence.
 * \param k   Length of the k-mers considered for syncmer selection.
 * \param seq The input nucleotide sequence.
 * \param s   Length of the internal s-mers (s < k) used to decide whether a k-mer is a syncmer.
 * \param t   Required offset of the minimal s-mer within the k-mer for it to be selected.
 * \return The set of hash values of all canonical open-syncmer k-mers found in \p seq.
 */
ankerl::unordered_dense::set<size_t> seq_to_syncmers(int k, const seqan3::dna5_vector &seq, int s, int t)
{
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    ankerl::unordered_dense::set<size_t> string_hashes{};

    make_string_to_hashvalues_open_syncmers_canonical(seq, string_hashes, k, s, t);
    return std::move(string_hashes);

}



}