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
inline uint64_t rotl64(uint64_t n, unsigned int c) {
    // assumes width is a power of 2
    const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
    // assert ( (c<=mask) &&"rotate by type width or more");
    c &= mask;
    return (n << c) | ( n >> ((-c) & mask));
}

//__attribute__((always_inline))
/*
*  corresponds to (hash * n) / 2^32 
*  but is 4 times faster than ordinary division
*  more efficient use of modulo operation
*/
inline uint32_t reduce(uint32_t hash, uint32_t n) {
    // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
    return (uint32_t) (((uint64_t) hash * n) >> 32);
}

//__attribute__((always_inline))
inline size_t getHashFromHash(uint64_t hash, int index, int blockLength) {
    uint32_t r = rotl64(hash, index * 21);
    return (size_t) reduce(r, blockLength) + ((size_t)index) * ((size_t)blockLength);
}

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

  inline FingerprintType fingerprint(const uint64_t hash) const {
    return (FingerprintType) hash ^ (hash >> 32);
  }

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

  Status AddAll(const vector<ItemType> &data, const size_t start, const size_t end) {
      return AddAll(data.data(),start,end);
  }

  Status AddAll(const ItemType* data, const size_t start, const size_t end);

  // Report if the item is inserted, with false positive rate.
  Status Contain(const ItemType &item) const;

  /* methods for providing stats  */
  // summary infomation
  std::string Info() const;

  // number of current inserted items;
  size_t Size() const { return size; }

  // size of the filter in bytes.
  size_t SizeInBytes() const { return arrayLength * sizeof(FingerprintType); }
};

struct t2val {
  uint64_t t2 = 0;
  uint64_t t2count = 0;
};

typedef struct t2val t2val_t;

const int blockShift = 18;

/*
* count number of occurences of each index and XOR to the corresponding hash values of each index
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
    
    //t2val_t * t2vals = new t2val_t[m];
    // temporary array H
    std::vector<t2val_t> t2vals_vec(m);
    t2val_t* t2vals = t2vals_vec.data();
    //cout << "start" << endl;
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

        // pick only index positions where only one unique has value points to => those are our start positions
        int* alone = new int[arrayLength];
        int alonePos = 0;
        for (size_t i = 0; i < arrayLength; i++) 
        {
            if (t2vals[i].t2count == 1) {
                alone[alonePos++] = i;
            }
        }
        cout << "AlonePos: " << alonePos << endl;
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
            cout << "reverseOrderPos: " << reverseOrderPos << "\tsize: " << size << endl;
        }
        delete[] tmp;
        delete[] tmpc;
        delete[] alone;

        
        if (reverseOrderPos == size) {
            break;
        }

        hashIndex++;

        // use a new random numbers
        delete hasher;
        hasher = new HashFamily();

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
