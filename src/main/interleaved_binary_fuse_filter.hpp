/*!\file
 * \author Jens-Uwe Ulrich <jens-uwe.ulrich AT hpi.de>
 * \brief Provides seqan3::interleaved_xor_filter.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/bit>

#include <sdsl/bit_vectors.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/detail/strong_type.hpp>

#include <bitset>

//#define XXH_INLINE_ALL

//#include "xxhash.h"

namespace ulrich2
{

/*!\brief The IXF binning directory. A data structure that efficiently answers set-membership queries for multiple bins.
 * \ingroup search_dream_index
 * \implements seqan3::cerealisable
 *
 * \details
 *
 * ### Binning Directory
 *
 * A binning directory is a data structure that can be used to determine set membership for elements.
 * For example, a common use case is dividing a database into a fixed number (e.g. 1024) bins by some means
 * of clustering (e.g. taxonomic binning or k-mer similarity clustering for genomic sequences).
 * For a query, the binning directory can now answer in which bins the query (probably) occurs.
 * In SeqAn we provide the Interleaved Bloom Filter (IBF) that can answer these queries efficiently.
 *
 * ### Interleaved XOR Filter (IXF)
 *
 * TODO: add description
 *
 * The Interleaved XOR Filter now applies the concept of a XOR Filter to multiple sets and provides a *global*
 * data structure to determine set membership of a query in `b` data sets/bins.
 * Conceptually, a XOR Filter is created for each bin using the same fixed length and fixed hash functions for each
 * filter. The resulting `b` XOR Filters are then interleaved such that the `i`'th l-bits if each XOR Filter are
 * adjacent to each other, with l being 2, 4, 8, 16, 32 or 64. For l=2: 
 * ```
 * Bloom Filter 0       Bloom Filter 1      Bloom Filter 2      Bloom Filter 3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * Where `x.y` denotes the `y`'th bit of the `x`'th XOR Filter.
 * ```
 * Interleaved XOR Filter
 * |0.0|0.1|1.0|1.1|2.0|2.1|3.0|3.1|0.2|0.3|1.2|1.3|2.2|2.3|3.2|3.3|
 * ```
 * A query can now be searched in all `b` bins by computing the `h` hash functions, retrieving the `h` sub-bitvectors of
 * length `b * l` starting at the positions indicated by the hash functions. The bitwise XOR of these sub-bitvectors yields
 * the binningvector, a bitvector of length `b * l` where the `i*l`'th bits indicate set membership in the `i`'th bin iff all 
 * `l` bits of `i` are 0.
 *
 * ### Querying
 * To query the Interleaved XOR Filter for a value, call seqan3::interleaved_xor_filter::membership_agent() and use
 * the returned seqan3::interleaved_xor_filter::membership_agent.
 *
 * To count the occurrences of a range of values in the Interleaved XOR Filter, call
 * seqan3::interleaved_xor_filter::counting_agent() and use
 * the returned seqan3::interleaved_xor_filter::counting_agent_type.
 *
 *
 * ### Thread safety
 *
 * The Interleaved XOR Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time). Furthermore, XOR Filters are immutable data structures
 * by itself.
 *
 * 
 */
template <typename FingerprintType = uint8_t>
//!\cond
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
//!\endcond
class interleaved_binary_fuse_filter
{
private:
    
    using data_type = sdsl::int_vector<>;

    //!\brief The number of bins specified by the user.
    size_t bins{};
    //!\brief The number of bins stored in the IXF (next multiple of 64 of `bins`).
    size_t technical_bins{};
    //!\brief The size of each bin in bits.
    size_t bin_size_{};
    //!\brief number of elements that can be stored in each bin.
    size_t max_bin_elements{};
    
    //!\brief The number of 64-bit integers needed to store `bins` many bits (e.g. `bins = 50` -> `bin_words = 1`).
    size_t bin_words{};
    //!\brief The length of each xor_filter block (filter_size/3)
    size_t segment_length{};
    size_t segment_length_mask{};
    size_t segment_count_length{};
    size_t segment_count{};
    //!\brief The int vector of fingerprints
    data_type data{};
    //!\brief number of bits used to store one hashed item in the XOR filter
    size_t ftype{};
    //!\brief seed for multiplicative hashing
    size_t seed{13572355802537770549ULL};

    size_t bins_per_batch{};

    uint64_t seed_set{0};

    inline constexpr uint64_t murmur64(uint64_t h) const
    {
        h += seed;
        h ^= h >> 33;
        h *= UINT64_C(0xff51afd7ed558ccd);
        h ^= h >> 33;
        h *= UINT64_C(0xc4ceb9fe1a85ec53);
        h ^= h >> 33;
        return h;
    }

    inline constexpr uint64_t rotl64(uint64_t n, unsigned int c) const
    {
        // assumes width is a power of 2
        const unsigned int mask = (CHAR_BIT * sizeof(n) - 1);
        // assert ( (c<=mask) &&"rotate by type width or more");
        c &= mask;
        return (n << c) | ( n >> ((-c) & mask));
    }

    inline FingerprintType fingerprint(const uint64_t hash) const 
    {
        return (FingerprintType)hash;
    }

    /*
    *  corresponds to (hash * n) / 2^32 
    *  but is 4 times faster than ordinary division
    *  more efficient use of modulo operation
    */
    inline constexpr uint32_t reduce(uint32_t hash, uint32_t n) const
    {
        // http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
        return (uint32_t) (((uint64_t) hash * n) >> 32);
    }


    inline __attribute__((always_inline)) size_t getHashFromHash(uint64_t hash, int index) 
    {
        __uint128_t x = (__uint128_t)hash * (__uint128_t)segment_count_length;
        uint64_t h = (uint64_t)(x >> 64);
        h += index * segment_length;
        // keep the lower 36 bits
        uint64_t hh = hash & ((1UL << 36) - 1);
        // index 0: right shift by 36; index 1: right shift by 18; index 2: no shift
        h ^= (size_t)((hh >> (36 - 18 * index)) & segment_length_mask);
        return h;
    }


    inline __attribute__((always_inline)) uint8_t mod3(uint8_t x) 
    {
        if (x > 2) 
        {
            x -= 3;
        }
        return x;
    }

    /*
    * count number of occurences of each index and XOR to the corresponding hash values of each index
    */
/*    void applyBlock(uint64_t* tmp, int b, int len, t2val_t* t2vals)
    {
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
/*    int applyBlock2(uint64_t* tmp, int b, int len, t2val_t*  t2vals, sdsl::int_vector<>& alone, int alonePos) 
    {
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
*/
    int find_alone_positions(std::vector<size_t>& elements, 
                             sdsl::int_vector<>& t2count, 
                             sdsl::int_vector<>& t2hash,
                             sdsl::int_vector<>& alone_positions,
                             sdsl::int_vector<>& rev_order)
    {
        rev_order[elements.size()] = 1;

        int block_bits = 1;
        while((size_t(1) << block_bits) < this->segment_count) 
        { 
            block_bits++; 
        }
        size_t block = size_t(1) << block_bits;
        sdsl::int_vector<> start_pos = sdsl::int_vector(block,0,64);

        for(uint32_t i = 0; i < uint32_t(1) << block_bits; i++) 
        { 
            start_pos[i] = i * this->max_bin_elements / block; 
        }

        for (size_t k : elements) 
        {
            uint64_t hash = murmur64(k);
//            uint64_t hash = XXH3_64bits_withSeed((char*) &k, sizeof(k), seed);
            size_t segment_index = hash >> (64 - block_bits);
            // We only overwrite when the hash was zero. Zero hash values
            // may be misplaced (unlikely).
            while(rev_order[start_pos[segment_index]] != 0) 
            {
                segment_index++;
                segment_index &= (size_t(1) << block_bits) - 1;
            }
            rev_order[start_pos[segment_index]] = hash;
            start_pos[segment_index] += 1;
        }
        uint8_t count_mask = 0;
        for (size_t i = 0; i < elements.size(); i++) 
        {
            uint64_t hash = rev_order[i];
            for (int hi = 0; hi < 3; hi++) 
            {
                int index = getHashFromHash(hash, hi);
                t2count[index] += 4;
                t2count[index] = t2count[index] ^ hi;
                t2hash[index] = t2hash[index] ^ hash;
                count_mask |= t2count[index];
            }
        }
        
        if (count_mask >= 0x80) 
        {
            // we have a possible counter overflow
            // this branch is never taken except if there is a problem in the hash code
            // in which case construction fails
            return 0;
        }

        size_t alonePos = 0;
        for (size_t i = 0; i < this->bin_size_; i++) 
        {
            alone_positions[alonePos] = i;
            int inc = (t2count[i] >> 2) == 1 ? 1 : 0;
            alonePos += inc;
        }
        return alonePos;
    }

    size_t fill_stack(sdsl::int_vector<>& reverse_order, 
                      sdsl::int_vector<>& reverse_h, 
                      sdsl::int_vector<>& t2count, 
                      sdsl::int_vector<>& t2hash, 
                      sdsl::int_vector<>& alone_positions,
                      int alone_pos)
    {
        size_t reverse_order_pos = 0;
        size_t h012[5];
        while (alone_pos > 0) 
        {
            alone_pos--;
            size_t index = alone_positions[alone_pos];
            if ((t2count[index] >> 2) == 1) 
            {

                // It is still there!
                uint64_t hash = t2hash[index];
                int found = t2count[index] & 3;
        
                reverse_h[reverse_order_pos] = found;
                reverse_order[reverse_order_pos] = hash;
        
                h012[0] = getHashFromHash(hash, 0);
                h012[1] = getHashFromHash(hash, 1);
                h012[2] = getHashFromHash(hash, 2);

                size_t index3 = h012[mod3(found + 1)];
                alone_positions[alone_pos] = index3;
                alone_pos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
                t2count[index3] -= 4;
                t2count[index3] = t2count[index3] ^ mod3(found + 1);
                t2hash[index3] = t2hash[index3] ^ hash;

                index3 = h012[mod3(found + 2)];
                alone_positions[alone_pos] = index3;
                alone_pos += ((t2count[index3] >> 2) == 2 ? 1 : 0);
                t2count[index3] -= 4;
                t2count[index3] = t2count[index3] ^ mod3(found + 2);
                t2hash[index3] = t2hash[index3] ^ hash;

                reverse_order_pos++;
            }
        }
        return reverse_order_pos;
    }

    void fill_filter(sdsl::int_vector<>& reverse_order, sdsl::int_vector<>& reverse_h, uint bin)
    {
        // the array h0, h1, h2, h0, h1, h2
        size_t h012[5];
        
        for (int i = reverse_order.size() - 1; i >= 0; i--) 
        {
            // the hash of the key we insert next
            uint64_t hash = reverse_order[i];
            int found = reverse_h[i];

            FingerprintType xor2 = fingerprint(hash);
            h012[0] = (this->bins * getHashFromHash(hash, 0)) + bin;
            h012[1] = (this->bins * getHashFromHash(hash, 1)) + bin;
            h012[2] = (this->bins * getHashFromHash(hash, 2)) + bin;
            h012[3] = h012[0];
            h012[4] = h012[1];
            data[h012[found]] = xor2 ^ data[h012[found + 1]] ^ data[h012[found + 2]];
            
        }
    }
    

    void add_elements(std::vector<std::vector<size_t>>& elements)
    {
        // stack sigma
        // order in which elements will be inserted into their corresponding bins
        std::vector<sdsl::int_vector<>> reverse_orders;
        // order in which hash seeds are used for element insertion
        std::vector<sdsl::int_vector<>> reverse_hs;
        size_t reverse_order_pos;
        //hashIndex = 0;
        // repeat until all xor filters can be build with same seed

        while (true)
        {
            // sets the same seed for all xor filters
            int i = 0;
            bool success = true;
            reverse_hs.clear();
            reverse_order_pos = 0;
            reverse_orders.clear();

            for (std::vector<size_t> vec : elements)
            {
                sdsl::int_vector<> t2count = sdsl::int_vector(this->bin_size_,0,8);
                sdsl::int_vector<> t2hash = sdsl::int_vector(this->bin_size_,0,64);
                sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
                sdsl::int_vector<> rev_order_i = sdsl::int_vector(vec.size(),0,64);
                sdsl::int_vector<> reverse_hi = sdsl::int_vector(vec.size(),0,8);
                int alone_position_nr = find_alone_positions(vec, t2count, t2hash, alone_positions, rev_order_i);
                
//                std::cout << alone_position_nr << std::endl;

                
                reverse_order_pos = fill_stack(std::ref(rev_order_i), std::ref(reverse_hi), t2count, t2hash, alone_positions, alone_position_nr);
                reverse_orders.emplace_back(std::move(rev_order_i));
                reverse_hs.emplace_back(std::move(reverse_hi));
                if (reverse_order_pos != vec.size())
                {
                    success = false;
                    set_seed();
                    t2count.clear();
                    t2hash.clear();
                    alone_positions.clear();
                    rev_order_i.clear();
                    reverse_hi.clear();
                    break;
                }

                i++;
            }
           
            if (success)
                break;

        }

        for (size_t i = 0; i < elements.size(); ++i)
        {
            fill_filter(reverse_orders[i], reverse_hs[i], i);
        }
        
//        std::cout << "Changed Seed " << seed_set << " times!" << std::endl;
    }

public:

    class membership_agent; // documented upon definition below
    template <std::integral value_t>
    class counting_agent_type; // documented upon definition below

    /*!\name Constructors, destructor and assignment
     * \{
     */
    interleaved_binary_fuse_filter() = default; //!< Defaulted.
    interleaved_binary_fuse_filter(interleaved_binary_fuse_filter const &) = default; //!< Defaulted.
    interleaved_binary_fuse_filter & operator=(interleaved_binary_fuse_filter const &) = default; //!< Defaulted.
    interleaved_binary_fuse_filter(interleaved_binary_fuse_filter &&) = default; //!< Defaulted.
    interleaved_binary_fuse_filter & operator=(interleaved_binary_fuse_filter &&) = default; //!< Defaulted.
    ~interleaved_binary_fuse_filter() = default; //!< Defaulted.

    /*!\brief Construct an uncompressed Interleaved XOR Filter.
     * \param bins_ The number of bins.
     * \param size  maximum number of elements to store in a single filter
     *
     * \attention This constructor can only be used to construct **uncompressed** Interleaved XOR Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/interleaved_xor_filter_constructor.cpp
     */
    interleaved_binary_fuse_filter(std::vector<std::vector<size_t>>& elements)
    {
        bins = elements.size();
        for (std::vector<size_t> v : elements)
        {
            if (v.size() > max_bin_elements)
                max_bin_elements = v.size();
        }
        
        this->segment_length = 1L << (int)floor(log(max_bin_elements) / log(3.33) + 2.25);;
        // the current implementation hardcodes a 18-bit limit to
        // to the segment length.
        if (this->segment_length > (1 << 18)) {
            this->segment_length = (1 << 18);
        }
        double size_factor = fmax(1.125, 0.875 + 0.25 * log(1000000) / log(max_bin_elements));
        size_t capacity = max_bin_elements * size_factor;
        size_t segment_count_ = (capacity + segment_length - 1) / segment_length - 2;
        this->bin_size_ = (segment_count_ + 2) * segment_length;
        this->segment_length_mask = this->segment_length - 1;
        this->segment_count = (this->bin_size_ + this->segment_length - 1) / this->segment_length;
        this->segment_count = this->segment_count <= 2 ? 1 : this->segment_count - 2;
        this->bin_size_ = (this->segment_count + 2) * this->segment_length;
        this->segment_count_length = this->segment_count * this->segment_length;
        //fingerprints = new FingerprintType[arrayLength]();

        //std::fill_n(fingerprints, arrayLength, 0);
        
        
        ftype = CHAR_BIT * sizeof(FingerprintType);
        std::cout << bin_size_ << "\t" << ftype << std::endl;
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (bin_size_ == 0)
            throw std::logic_error{"The size of a bin must be > 0."};


        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
            technical_bins  = bin_words << 3; // = bin_words * 8
        }
        else
        {
            bin_words = (bins + 15) >> 4; // = ceil(bins/8)
            technical_bins  = bin_words << 4; // = bin_words * 16
        }

        data = sdsl::int_vector<>(bin_size_ * bins, 0, ftype);
        //data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
    
        add_elements(elements);

        
    }

/*
    interleaved_xor_filter(size_t bins_, size_t max_bin_elements)
    {
        this->bins = bins_;
        std::cout << "Max Elements: " << max_bin_elements << std::endl;
        bin_size_ = 32 + 1.23 * max_bin_elements;
        block_length = bin_size_ / 3;
        ftype = CHAR_BIT * sizeof(FingerprintType);
        std::cout << bin_size_ << "\t" << ftype << std::endl;
        bins_per_batch = 64/ftype;

        if (bins == 0)
            throw std::logic_error{"The number of bins must be > 0."};
        if (bin_size_ == 0)
            throw std::logic_error{"The size of a bin must be > 0."};

        hash_shift = std::countl_zero(bin_size_);

        if (ftype == 8)
        {
            bin_words = (bins + 7) >> 3; // = ceil(bins/8)
            technical_bins  = bin_words << 3; // = bin_words * 8
        }
        else
        {
            bin_words = (bins + 15) >> 4; // = ceil(bins/8)
            technical_bins  = bin_words << 4; // = bin_words * 16
        }
        data = sdsl::int_vector<>(bins * bin_size_, 0, ftype);
    }


    bool add_bin_elements(size_t bin, std::vector<size_t>& elements)
    {
       
        size_t reverse_order_pos = 0;
        std::vector<t2val_t> t2vals_vec(bin_size_);
        sdsl::int_vector<> alone_positions = sdsl::int_vector(bin_size_);
        int alone_position_nr = find_alone_positions(elements, t2vals_vec, alone_positions);
        sdsl::int_vector<> rev_order = sdsl::int_vector(elements.size(),0,64);
        sdsl::int_vector<> reverse_h = sdsl::int_vector(elements.size(),0,8);
        reverse_order_pos = fill_stack(std::ref(rev_order), std::ref(reverse_h), t2vals_vec, alone_positions, elements.size(), alone_position_nr);
        if (reverse_order_pos != elements.size())
        {
            return false;
        }

        fill_filter(rev_order, reverse_h, bin);

        return true;
    }
*/   
    void clear()
    {
        for (size_t idx = 0; idx < data.size(); ++idx)
            data[idx] = 0;
    }


    void set_seed()
    {
        ::std::random_device random;
        seed = random();
        seed <<= 32;
        seed |= random();
        seed_set++;
        std::cout << "Seed change number: " << seed_set << std::endl;
    }

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns a seqan3::interleaved_xor_filter::membership_agent to be used for lookup.
     * `seqan3::interleaved_xor_filter::membership_agent`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_construction.cpp
     * \sa seqan3::interleaved_bloom_filter::membership_agent::bulk_contains
     */
    membership_agent membership_agent() const
    {
        return typename interleaved_binary_fuse_filter<FingerprintType>::membership_agent{*this};
    }

    /*!\brief Returns a seqan3::interleaved_xor_filter::counting_agent_type to be used for counting.
     * `seqan3::interleaved_xor_filter::counting_agent_type`s constructed for this Interleaved XOR Filter.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent_construction.cpp
     * \sa seqan3::interleaved_bloom_filter::counting_agent_type::bulk_count
     */
    template <typename value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
    //!\}

        /*!\brief Returns the number of bins that the Interleaved XOR Filter manages.
     * \returns The number of bins.
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }

    /*!\brief Returns the size of a single bin that the Interleaved XOR Filter manages.
     * \returns The size in bits of a single bin.
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return data.bit_size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Test for equality.
     * \param[in] lhs A `seqan3::interleaved_xor_filter`.
     * \param[in] rhs `seqan3::interleaved_xor_filter` to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    
    friend bool operator==(interleaved_binary_fuse_filter const & lhs, interleaved_binary_fuse_filter const & rhs) noexcept
    {
        return std::tie(lhs.bins, lhs.bins_per_batch, lhs.bin_size_, lhs.ftype, lhs.bin_words, lhs.segment_count_length,
                        lhs.segment_length, lhs.segment_length_mask, lhs.seed, lhs.data) ==
               std::tie(rhs.bins, rhs.bins_per_batch, rhs.bin_size_, rhs.ftype, rhs.bin_words, rhs.segment_count_length,
                        rhs.segment_length, rhs.segment_length_mask, rhs.seed, rhs.data);
    }

    /*!\brief Test for inequality.
     * \param[in] lhs A `seqan3::interleaved_xor_filter`.
     * \param[in] rhs `seqan3::interleaved_xor_filter` to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(interleaved_binary_fuse_filter const & lhs, interleaved_binary_fuse_filter const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(bins);
        archive(bin_size_);
        archive(bin_words);
        archive(ftype);
        if (ftype != (CHAR_BIT * sizeof(FingerprintType)))
        {
            throw std::logic_error{"The interleaved XOR filter was built with a fingerprint size of " + std::to_string(ftype) +
                                   " but it is being read into an interleaved XOR filter with fingerprint of size " +
                                   std::to_string(CHAR_BIT * sizeof(FingerprintType)) + "."};
        }
        archive( bins_per_batch);
        archive(segment_count_length);
        archive(segment_length);
        archive(segment_length_mask); 
        archive(seed);
        archive(data);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the seqan3::interleaved_xor_filter.
 * 
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/membership_agent_construction.cpp
 */
template <typename FingerprintType>
class interleaved_binary_fuse_filter<FingerprintType>::membership_agent
{
private:
    //!\brief The type of the augmented seqan3::interleaved_xor_filter.
    using iff_t = interleaved_binary_fuse_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_xor_filter.
    iff_t const * iff_ptr{nullptr};

public:
    class binning_bitvector;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default; //!< Defaulted.
    membership_agent(membership_agent const &) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default; //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default; //!< Defaulted.
    ~membership_agent() = default; //!< Defaulted.

    /*!\brief Construct a membership_agent from a seqan3::interleaved_xor_filter.
     * \private
     * \param ibf The seqan3::interleaved_xor_filter.
     */
    explicit membership_agent(iff_t const & iff) :
        iff_ptr(std::addressof(iff)), result_buffer(iff.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    binning_bitvector result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of a given value.
     * \param[in] value The raw value to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_xor_filter::membership_agent for each thread.
     */

     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) & noexcept
    {
        //std::cout << "Start bulk contains " << std::endl;
        assert(iff_ptr != nullptr);
        assert(result_buffer.size() == iff_ptr->bin_count());

        size_t bins = iff_ptr->bin_count();

        uint8_t ftype = iff_ptr->ftype;
        uint8_t bins_per_batch = iff_ptr->bins_per_batch;
        size_t segment_count_length = iff_ptr->segment_count_length;
        size_t segment_length = iff_ptr->segment_length;
        size_t segment_length_mask = iff_ptr->segment_length_mask;

        uint64_t hash = iff_ptr->murmur64(value);
//        uint64_t hash = XXH3_64bits_withSeed((char*) &value, sizeof(value), iff_ptr->seed);
        FingerprintType f = iff_ptr->fingerprint(hash);
        //std::cout << "Fingerprint for " << value << " : " << f << std::endl;
        __uint128_t x = (__uint128_t)hash * (__uint128_t) segment_count_length;
        uint32_t h0 = (uint64_t)(x >> 64);
        uint32_t h1 = h0 + segment_length;
        uint32_t h2 = h1 + segment_length;
        uint64_t hh = hash;
        h1 ^= (size_t)((hh >> 18) & segment_length_mask);
        h2 ^= (size_t)((hh) & segment_length_mask);

        h0 = h0*bins*ftype;
        h1 = h1*bins*ftype;
        h2 = h2*bins*ftype;

        // for bins/(64/ftype) : 

        // concatenate (64/ftype) the fingerprint
        
        uint64_t fc64 = f;
        for (uint8_t b = 1; b < bins_per_batch ; ++b)
        {
            fc64 = (fc64 << ftype) | f;
        }
        //std::cout << std::bitset<64>(fc64) << std::endl;


        for (size_t batch = 0; batch < iff_ptr->bin_words; ++batch)
        {
            size_t batch_start = batch * 64;
            uint64_t v = fc64 ^ iff_ptr->data.get_int(h0 + batch_start, 64) 
                              ^ iff_ptr->data.get_int(h1 + batch_start, 64) 
                              ^ iff_ptr->data.get_int(h2 + batch_start, 64);
//            std::cout << std::bitset<64>(v) << std::endl;
            size_t used_bins = batch * bins_per_batch;
            uint8_t bits = 0;
            for (size_t bin = 0; bin < bins_per_batch; ++bin)
            {
                if (used_bins + bin == result_buffer.size())
                    break;

                uint64_t tmp = v << ((bins_per_batch - (bin+1)) * ftype ); //v << (bin*ftype);
                //uint8_t bits = (uint8_t) (tmp >> (64-ftype));
                uint8_t tmpb = std::bitset<8>(tmp >> (64-ftype)).none() << bin;
                bits |= tmpb;
                //std::cout << std::bitset<8>(bits) << std::endl;
                   
            }
            //std::cout << std::bitset<8>(bits) << std::endl;
            result_buffer.data.set_int(used_bins, bits, ftype);
        }

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    //!\}
     [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) && noexcept = delete;
};

//!\brief A bitvector representing the result of a call to `bulk_contains` of the seqan3::interleaved_bloom_filter.
template <typename FingerprintType>
class interleaved_binary_fuse_filter<FingerprintType>::membership_agent::binning_bitvector
{
private:
    //!\brief The underlying datatype to use.
    using data_type = sdsl::bit_vector;
    //!\brief The bitvector.
    data_type data{};

    friend class membership_agent;

    template <std::integral value_t>
    friend class counting_vector;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_bitvector() = default; //!< Defaulted.
    binning_bitvector(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector const &) = default; //!< Defaulted.
    binning_bitvector(binning_bitvector &&) = default; //!< Defaulted.
    binning_bitvector & operator=(binning_bitvector &&) = default; //!< Defaulted.
    ~binning_bitvector() = default; //!< Defaulted.

    //!\brief Construct with given size.
    explicit binning_bitvector(size_t const size) :
        data(size)
    {}
    //!\}

    //!\brief Returns the number of elements.
    size_t size() const noexcept
    {
        return data.size();
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the first element of the container.
    auto begin() noexcept
    {
        return data.begin();
    }

    //!\copydoc begin()
    auto begin() const noexcept
    {
        return data.begin();
    }

    //!\brief Returns an iterator to the element following the last element of the container.
    auto end() noexcept
    {
        return data.end();
    }

    //!\copydoc end()
    auto end() const noexcept
    {
        return data.end();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Test for equality.
    friend bool operator==(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return lhs.data == rhs.data;
    }

    //!\brief Test for inequality.
    friend bool operator!=(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
     //!\brief Return the i-th element.
    auto operator[](size_t const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    //!\copydoc operator[]()
    auto operator[](size_t const i) const noexcept
    {
        assert(i < size());
        return data[i];
    }

    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}
};

/*!\brief A data structure that behaves like a std::vector and can be used to consolidate the results of multiple calls
 *        to seqan3::interleaved_xor_filter::membership_agent::bulk_contains.
 * \ingroup search_dream_index
 * \tparam value_t The type of the count. Must model std::integral.
 *
 * \details
 *
 * When using the seqan3::interleaved_xor_filter::membership_agent::bulk_contains operation, a common use case is to
 * add up, for example, the results for all k-mers in a query. This yields, for each bin, the number of k-mers of a
 * query that are in the respective bin. Such information can be used to apply further filtering or abundance estimation
 * based on the k-mer counts.
 *
 * The seqan3::counting_vector offers an easy way to add up the individual
 * seqan3::interleaved_xor_filter::membership_agent::binning_bitvector by offering an `+=` operator.
 *
 * The `value_t` template parameter should be chosen in a way that no overflow occurs if all calls to `bulk_contains`
 * return a hit for a specific bin. For example, `uint8_t` will suffice when processing short Illumina reads, whereas
 * long reads will require at least `uint32_t`.
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_vector.cpp
 */

template<std::integral value_t>
class counting_vector : public std::vector<value_t>
{
private:
    //!\brief The base type.
    using base_t = std::vector<value_t>;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_vector() = default; //!< Defaulted.
    counting_vector(counting_vector const &) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector const &) = default; //!< Defaulted.
    counting_vector(counting_vector &&) = default; //!< Defaulted.
    counting_vector & operator=(counting_vector &&) = default; //!< Defaulted.
    ~counting_vector() = default; //!< Defaulted.

    using base_t::base_t;
    //typename ixf_t::membership_agent membership_agent;
    //!\}

    /*!\brief Bin-wise adds the bits of a seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \tparam rhs_t The type of the right-hand side.
     *         Must be seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \param rhs The seqan3::interleaved_xor_filter::membership_agent::binning_bitvector.
     * \attention The counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    template <typename rhs_t>
//    //!\cond
//        requires std::same_as<rhs_t,
//                              interleaved_xor_filter<FingerprintType, data_layout_mode>::membership_agent::binning_bitvector>
//                 ||
//                 std::same_as<rhs_t,
//                              interleaved_xor_filter<FingerprintType, data_layout::compressed>::membership_agent::binning_bitvector>
//    //!\endcond
    counting_vector & operator+=(rhs_t const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        // Each iteration can handle 64 bits, so we need to iterate `((rhs.size() + 63) >> 6` many times
        for (size_t batch = 0, bin = 0; batch < ((rhs.size() + 63) >> 6); bin = 64 * ++batch)
        {
            size_t tmp = rhs.data.get_int(batch * 64); // get 64 bits starting at position `batch * 64`
            if (tmp ^ (1ULL<<63)) // This is a special case, because we would shift by 64 (UB) in the while loop.
            {
                while (tmp > 0)
                {
                    // Jump to the next 1 and increment the corresponding vector entry.
                    uint8_t step = std::countr_zero(tmp);
                    bin += step++;
                    tmp >>= step;
                    ++(*this)[bin++];
                }
            }
            else
            {
                ++(*this)[bin + 63];
            }
        }
        return *this;
    }

    /*!\brief Bin-wise addition of two `seqan3::counting_vector`s.
     * \param rhs The other seqan3::counting_vector.
     * \attention The seqan3::counting_vector must be at least as big as `rhs`.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
        assert(this->size() >= rhs.size()); // The counting vector may be bigger than what we need.

        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }
};

/*!\brief Manages counting ranges of values for the seqan3::interleaved_xor_filter.
 * \attention Calling seqan3::interleaved_xor_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <typename FingerprintType>
template <std::integral value_t>
class interleaved_binary_fuse_filter<FingerprintType>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::interleaved_xor_filter.
    using iff_t = interleaved_binary_fuse_filter<FingerprintType>;

    //!\brief A pointer to the augmented seqan3::interleaved_xor_filter.
    iff_t const * iff_ptr{nullptr};

    //!\brief Store a seqan3::interleaved_xor_filter::membership_agent to call `bulk_contains`.
    typename iff_t::membership_agent membership_agent;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    explicit counting_agent_type(iff_t const & iff) :
        iff_ptr(std::addressof(iff)), membership_agent(iff), result_buffer(iff.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * seqan3::interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector<value_t> const & bulk_count(value_range_t && values) & noexcept
    {
        assert(iff_ptr != nullptr);
        assert(result_buffer.size() == iff_ptr->bin_count());

        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, 0);

        for (auto && value : values)
            result_buffer += membership_agent.bulk_contains(value);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector<value_t> const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}

};

} // namespace seqan3
