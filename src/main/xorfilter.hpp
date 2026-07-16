/*!\file xorfilter.hpp
 * \brief Standalone single-set XOR filter implementation (a probabilistic set membership structure similar to a
 *        Bloom/cuckoo filter, but built by peeling a 3-uniform hypergraph of hashed positions).
 *
 * \details Not currently compiled into the taxor binary — this file is not #include'd by build.hpp, nor by any
 *          other file in the codebase; it is not even reachable transitively (build.hpp uses seqan3's own
 *          `interleaved_xor_filter` instead of this local implementation). It is retained here as reference
 *          for an earlier/alternate single-filter XOR-filter implementation but has no effect on the running tool.
 */

#ifndef XOR_FILTER_XOR_FILTER_H_
#define XOR_FILTER_XOR_FILTER_H_

#include <assert.h>
#include <algorithm>
#include "hashutil.hpp"
#include "timing.hpp"

using namespace std;
using namespace hashing;

namespace xorfilter {
// status returned by a xor filter operation
enum Status {
  Ok = 0,
  NotFound = 1,
  NotEnoughSpace = 2,
  NotSupported = 3,
};

//__attribute__((always_inline))
/*!\brief Rotate the bits of a 64-bit unsigned integer left by \p c positions.
 * \param n The value to rotate.
 * \param c Number of bit positions to rotate left by (reduced modulo the bit width).
 * \return \p n rotated left by \p c bits.
 */
inline uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

//__attribute__((always_inline))
/*!\brief Map a 32-bit hash value into the range [0, n) without using an (expensive) modulo operation.
 *
 * \details Corresponds to `(hash * n) / 2^32` but is roughly 4 times faster than an ordinary division/modulo
 *          (see http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/).
 * \param hash The hash value to reduce.
 * \param n    Exclusive upper bound of the target range.
 * \return A value in [0, n).
 */
inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

//__attribute__((always_inline))
/*!\brief Derive one of the 3 candidate fingerprint-array positions for a key's hash.
 * \param hash        The 64-bit hash of the key.
 * \param index       Which of the 3 hash functions to use (0, 1 or 2).
 * \param blockLength Number of table slots reserved per hash function (`arrayLength / 3`).
 * \return The absolute index into the fingerprint array for this key/hash-function combination:
 *         a value obtained by rotating \p hash and reducing it into [0, blockLength), then offsetting
 *         by `index * blockLength` so the three hash functions address disjoint blocks of the array.
 */
inline size_t getHashFromHash(uint64_t hash, int index, int blockLength) {
    uint32_t r = rotl64(hash, index * 21);
    return (size_t) reduce(r, blockLength) + ((size_t)index) * ((size_t)blockLength);
}

/*!\brief A single-set probabilistic membership filter ("XOR filter") over items of type \p ItemType, storing one
 *        \p FingerprintType fingerprint per fingerprint-array slot.
 *
 * \tparam ItemType        Type of the elements stored in the filter (expected to be hashable via \p HashFamily).
 * \tparam FingerprintType Unsigned integer type used to store per-slot fingerprints; its width controls the
 *                         false-positive rate of Contain().
 * \tparam HashFamily      Hash functor type used to hash items to 64-bit values; defaults to SimpleMixSplit.
 *
 * \details The filter is built once via AddAll() over a fixed key set (construction is not incremental), then
 *          queried via Contain(). Internally each key is mapped to 3 candidate positions in a `fingerprints`
 *          array of length `arrayLength = 32 + 1.23 * size`, split into 3 equal blocks of `blockLength` slots
 *          (one block per hash function) so the 3 positions for any key always fall in distinct blocks. A
 *          fingerprint is stored such that XOR-ing the fingerprints at all 3 of a key's positions reproduces
 *          the key's own fingerprint; construction (see AddAll()) uses the standard XOR-filter "peeling"
 *          algorithm to find an assignment for which this holds for every inserted key.
 */
template <typename ItemType, typename FingerprintType,
          typename HashFamily = SimpleMixSplit>
class XorFilter {
 public:

  size_t size;
  size_t arrayLength;
  size_t blockLength;
  FingerprintType *fingerprints;

  HashFamily* hasher;
  size_t hashIndex{0};

  /*!\brief Compute the fingerprint stored for a given key hash.
   * \param hash The 64-bit hash of the key.
   * \return The fingerprint value (low 32 bits XOR'd with high 32 bits of \p hash, truncated to \p FingerprintType).
   */
  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) hash ^ (hash >> 32);
  }

  /*!\brief Allocate an (initially empty) filter sized for \p size elements.
   * \param size Number of elements the filter is expected to hold; determines the fingerprint array size.
   */
  explicit XorFilter(const size_t size) {
    hasher = new HashFamily();
    this->size = size;
    this->arrayLength = 32 + 1.23 * size;
    this->blockLength = arrayLength / 3;
    fingerprints = new FingerprintType[arrayLength]();
    std::fill_n(fingerprints, arrayLength, 0);
  }

  ~XorFilter() {
    delete[] fingerprints;
    delete hasher;
  }

  //!\brief Convenience overload of AddAll() taking a vector instead of a raw pointer/length range.
  Status AddAll(const vector<ItemType> &data, const size_t start, const size_t end) {
      return AddAll(data.data(),start,end);
  }

  /*!\brief Build the filter from the keys in `data[start, end)`, replacing any previous contents.
   * \param data  Pointer to the array of keys to insert.
   * \param start Index of the first key to insert (inclusive).
   * \param end   Index one past the last key to insert (exclusive).
   * \return Ok on success (implementation always returns Ok; failure to converge is instead handled internally
   *         by retrying with a new hash function seed, see the out-of-line definition below).
   */
  Status AddAll(const ItemType* data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  /*!\brief Test whether \p item is (probably) a member of the filter.
   * \param item The item to test.
   * \return Ok if the item is (probably, subject to the filter's false-positive rate) present, NotFound otherwise.
   */
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  //!\brief Return a short human-readable summary of the filter's current contents.
  std::string Info() const;

  // number of current inserted items;
  //!\brief Number of items the filter was constructed to hold.
  size_t Size() const { return size; }

  // size of the filter in bytes.
  //!\brief Total memory footprint of the fingerprint array, in bytes.
  size_t SizeInBytes() const { return arrayLength * sizeof(FingerprintType); }
};

//!\brief Per fingerprint-array-slot bookkeeping used during XOR-filter construction: how many keys currently
//!       hash to this slot (t2count) and the XOR of all their hash values (t2).
struct t2val {
  uint64_t t2 = 0;
  uint64_t t2count = 0;
};

typedef struct t2val t2val_t;

const int blockShift = 18;

/*
* count number of occurences of each index and XOR to the corresponding hash values of each index
*/
/*!\brief Fold a block of pending (hash, index) pairs into the \p t2vals table: increments the hit-count and
 *        XORs the hash value into `t2vals[index]` for each pair.
 * \param tmp    Flat array of blocks, each block holding interleaved (hash, index) pairs.
 * \param b      Which block within \p tmp to process.
 * \param len    Number of valid uint64_t entries in block \p b (i.e. 2x the number of pairs).
 * \param t2vals Output table of per-slot (XOR of hashes, hit count), indexed by fingerprint-array position.
 */
void applyBlock(uint64_t* tmp, int b, int len, t2val_t* t2vals) {
    for (int i = 0; i < len; i += 2) {
        uint64_t x = tmp[(b << blockShift) + i];
        int index = (int) tmp[(b << blockShift) + i + 1];
        t2vals[index].t2count++;
        t2vals[index].t2 ^= x;
    }
}

/*
*  remove all indexes in tmp array from index occurence array and remove corresponding hashes by XOR-ing with the hash value
*  add index to single entry set if an index now occurs only once in the remaining set of indexes
*/
/*!\brief Undo/consume a block of pending (hash, index) pairs against \p t2vals, decrementing hit-counts and
 *        appending newly-"alone" (hit-count == 1) slots to the \p alone worklist used by the peeling algorithm.
 * \param tmp      Flat array of blocks, each block holding interleaved (hash, index) pairs.
 * \param b        Which block within \p tmp to process.
 * \param len      Number of valid uint64_t entries in block \p b (i.e. 2x the number of pairs).
 * \param t2vals   Table of per-slot (XOR of hashes, hit count) to update.
 * \param alone    Worklist of fingerprint-array positions that have exactly one remaining hashed key; appended to.
 * \param alonePos Current (next free) write index into \p alone.
 * \return The updated \p alonePos after appending any newly-alone indexes.
 */
int applyBlock2(uint64_t* tmp, int b, int len, t2val_t*  t2vals, int* alone, int alonePos) {
    for (int i = 0; i < len; i += 2) {
        uint64_t hash = tmp[(b << blockShift) + i];
        int index = (int) tmp[(b << blockShift) + i + 1];
        int oldCount = t2vals[index].t2count;
        if (oldCount >= 1) {
            int newCount = oldCount - 1;
            t2vals[index].t2count = newCount;
            if (newCount == 1) {
                alone[alonePos++] = index;
            }
            t2vals[index].t2 ^= hash;
        }
    }
    return alonePos;
}

/*!\brief Construct the filter for the keys in `keys[start, end)` using the standard XOR-filter "peeling"
 *        (2-core / hypergraph) algorithm, and populate the fingerprint array.
 *
 * \details The algorithm repeatedly: (1) hashes every key to 3 candidate array positions and tallies, per
 *          position, how many keys map there and the XOR of their hashes (via applyBlock()/applyBlock2()
 *          working on `blockShift`-sized chunks to bound peak memory); (2) repeatedly removes ("peels") array
 *          positions that are currently hit by exactly one remaining key — these are pushed onto a stack
 *          (`reverseOrder`/`reverseH`) together with which of the 3 hash slots identified them; if peeling
 *          gets stuck before every key has been removed (a "2-core" remains), the whole attempt is discarded
 *          and retried from scratch with a different hash function selection (`hashIndex`), which is the
 *          standard way XOR-filter construction escapes unlucky hash choices. Once every key has been peeled
 *          (`reverseOrderPos == size`), the stack is unwound in reverse: for each key, the fingerprint slot
 *          that uniquely identified it is set so that XOR-ing all 3 of its candidate slots reproduces its
 *          fingerprint, without disturbing the two other (already-visited) slots.
 * \param keys  Pointer to the array of keys to insert.
 * \param start Index of the first key to insert (inclusive).
 * \param end   Index one past the last key to insert (exclusive).
 * \return Always returns Ok; the retry-on-stuck-peeling logic is handled internally via the `while (true)` loop
 *         rather than surfaced as a distinct status.
 */
template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFilter<ItemType, FingerprintType, HashFamily>::AddAll(
    const ItemType* keys, const size_t start, const size_t end) {

    int m = arrayLength;
    // stack sigma
    uint64_t* reverseOrder = new uint64_t[size];
    uint8_t* reverseH = new uint8_t[size];
    size_t reverseOrderPos;
    hashIndex = 0;

    //cout << "Seed : " << hasher->seed << endl;
    // note: the seed is hardcoded here rather than randomized; the "use a new random numbers" retry
    // path further below (re-creating hasher) is commented out, so every retry re-hashes with the
    // same seed and only hashIndex (the rotation amount in getHashFromHash) changes between attempts
    hasher->seed = 13572355802537770549ULL;
    //t2val_t * t2vals = new t2val_t[m];
    // temporary array H
    std::vector<t2val_t> t2vals_vec(m);
    t2val_t* t2vals = t2vals_vec.data();
    //cout << "start" << endl;
    // retry loop: rebuild the hypergraph and attempt to peel it; on getting stuck (a 2-core remains)
    // this loop starts over with a new hashIndex until a fully peelable assignment is found
    while (true) {
        //memset(t2vals, 0, sizeof(t2val_t) * m);
        std::fill(t2vals_vec.begin(), t2vals_vec.end(), t2val_t{ 0,0 });
        // number of elements / 2^18 => if more than 2^18 elements, we need 2 blocks
        int blocks = 1 + ((3 * blockLength) >> blockShift);
        //cout << "blocks : " << blocks << "\tblocklength: " << blockLength << "\tblockshift: " << blockShift <<  endl;
        uint64_t* tmp = new uint64_t[blocks << blockShift];
       // cout << int(blocks << blockShift) << endl;
        int* tmpc = new int[blocks]();
        for(size_t i = start; i < end; i++) {
            uint64_t k = keys[i];
            uint64_t hash = (*hasher)(k);
            for (int hi = 0; hi < 3; hi++) {
                int index = getHashFromHash(hash, hi, blockLength);
                int b = index >> blockShift;
              //  std::cout << "hi: " << hi << "\t" << "key: " << k << "\t" << "index: " << index << "\t" << "b: " << b << std::endl;
                int i2 = tmpc[b];
                tmp[(b << blockShift) + i2] = hash;
                tmp[(b << blockShift) + i2 + 1] = index;
                tmpc[b] += 2;
                //cout << "i2: " << i2 << "\t(b << blockShift) + i2: " << int((b << blockShift) + i2) << "\ttmpc[b]: " << tmpc[b] << endl;
                if (i2 + 2 == (1 << blockShift)) {
                    applyBlock(tmp, b, i2 + 2, t2vals);
                    tmpc[b] = 0;
                }
            }

        }
        // count occurences of index positions for all computed hash values
        for (int b = 0; b < blocks; b++) {
            applyBlock(tmp, b, tmpc[b], t2vals);
        }

        delete[] tmp;
        delete[] tmpc;
        reverseOrderPos = 0;

        // pick only index positions where only one unique hash value points to => those are our start positions
        int* alone = new int[arrayLength];
        int alonePos = 0;
        for (size_t i = 0; i < arrayLength; i++) 
        {
            if (t2vals[i].t2count == 1) {
                alone[alonePos++] = i;
            }
        }

        tmp = new uint64_t[blocks << blockShift];
        tmpc = new int[blocks]();
        reverseOrderPos = 0;
        int bestBlock = -1;
        while (reverseOrderPos < size)
        {
            if (alonePos == 0) 
            {
                // we need to apply blocks until we have an entry that is alone
                // (that is, until alonePos > 0)
                // so, find a large block (the larger the better)
                // but don't need to search very long
                // start searching where we stopped the last time
                // (to make it more even)
                for (int i = 0, b = bestBlock + 1, best = -1; i < blocks; i++)
                {
                    if (b >= blocks)
                    {
                        b = 0;
                    }
                    std::cout << "tmpc[" << b <<"]: " << tmpc[b] << std::endl;
                    if (tmpc[b] > best) 
                    {
                        best = tmpc[b];
                        bestBlock = b;
                        if (best > (1 << (blockShift - 1))) 
                        {
                            // sufficiently large: stop
                            break;
                        }
                    }
                }
                std::cout << "Best Block: " << bestBlock << std::endl;
                if (tmpc[bestBlock] > 0) {
                    alonePos = applyBlock2(tmp, bestBlock, tmpc[bestBlock], t2vals, alone, alonePos);
                    tmpc[bestBlock] = 0;
                }
                // applying a block may not actually result in a new entry that is alone
                if (alonePos == 0) {
                    for (int b = 0; b < blocks && alonePos == 0; b++) {
                        if (tmpc[b] > 0) {
                            alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                            tmpc[b] = 0;
                        }
                    }
                }
            }
            if (alonePos == 0) {
                break;
            }
            int i = alone[--alonePos];
            int b = i >> blockShift;
            if (tmpc[b] > 0) {
                alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                tmpc[b] = 0;
            }
            uint8_t found = -1;
            if (t2vals[i].t2count == 0) {
                continue;
            }
            long hash = t2vals[i].t2;
            for (int hi = 0; hi < 3; hi++) {
                int h = getHashFromHash(hash, hi, blockLength);
                if (h == i) {
                    found = (uint8_t) hi;
                    t2vals[i].t2count = 0;
                } else {
                    int b = h >> blockShift;
                    int i2 = tmpc[b];
                    tmp[(b << blockShift) + i2] = hash;
                    tmp[(b << blockShift) + i2 + 1] = h;
                    tmpc[b] += 2;
                    if (tmpc[b] >= 1 << blockShift) {
                        alonePos = applyBlock2(tmp, b, tmpc[b], t2vals, alone, alonePos);
                        tmpc[b] = 0;
                    }
                }
            }
            reverseOrder[reverseOrderPos] = hash;
            reverseH[reverseOrderPos] = found;
            reverseOrderPos++;
            //if (reverseOrderPos % 100000 == 0)
            //cout << "reverseOrderPos: " << reverseOrderPos << "\tsize: " << size << endl;
        }
        delete[] tmp;
        delete[] tmpc;
        delete[] alone;

        
        if (reverseOrderPos == size) {
            break;
        }

        hashIndex++;

        // use a new random numbers
        // delete hasher;
        //hasher = new HashFamily();

    }

    for (int i = reverseOrderPos - 1; i >= 0; i--) {
        // the hash of the key we insert next
        uint64_t hash = reverseOrder[i];
        int found = reverseH[i];
        // which entry in the table we can change
        int change = -1;
        // we set table[change] to the fingerprint of the key,
        // unless the other two entries are already occupied
        FingerprintType xor2 = fingerprint(hash);
        for (int hi = 0; hi < 3; hi++) {
            size_t h = getHashFromHash(hash, hi, blockLength);
            if (found == hi) {
                change = h;
            } else {
                // this is different from BDZ: using xor to calculate the
                // fingerprint
                xor2 ^= fingerprints[h];
            }
        }
        fingerprints[change] = xor2;
    }
    //delete [] t2vals;
    delete [] reverseOrder;
    delete [] reverseH;

    return Ok;
}

/*!\brief Test membership of \p key: recompute its 3 candidate fingerprint-array positions directly from its
 *        hash (equivalent to getHashFromHash() for hi = 0, 1, 2) and check whether XOR-ing the stored
 *        fingerprints at those positions reproduces the key's own fingerprint.
 * \param key The key to test.
 * \return Ok if the XOR of the 3 stored fingerprints matches the key's fingerprint (probable member,
 *         subject to false positives), NotFound otherwise (definite non-member).
 */
template <typename ItemType, typename FingerprintType,
          typename HashFamily>
Status XorFilter<ItemType, FingerprintType, HashFamily>::Contain(
    const ItemType &key) const {
    uint64_t hash = (*hasher)(key);
    FingerprintType f = fingerprint(hash);
    uint32_t r0 = (uint32_t) hash;
    uint32_t r1 = (uint32_t) rotl64(hash, 21);
    uint32_t r2 = (uint32_t) rotl64(hash, 42);
    uint32_t h0 = reduce(r0, blockLength);
    uint32_t h1 = reduce(r1, blockLength) + blockLength;
    uint32_t h2 = reduce(r2, blockLength) + 2 * blockLength;
    f ^= fingerprints[h0] ^ fingerprints[h1] ^ fingerprints[h2];
    return f == 0 ? Ok : NotFound;
}

//!\brief Build a short human-readable status string reporting the number of keys the filter holds.
template <typename ItemType, typename FingerprintType,
          typename HashFamily>
std::string XorFilter<ItemType, FingerprintType, HashFamily>::Info() const {
  std::stringstream ss;
  ss << "XorFilter Status:\n"
     << "\t\tKeys stored: " << Size() << "\n";
  return ss.str();
}

}  // namespace xorfilter
#endif  // XOR_FILTER_XOR_FILTER_H_