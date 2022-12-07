#pragma once

#include <ranges>

#include <seqan3/search/dream_index/interleaved_binary_fuse_filter.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace hixf
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [seqan3::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [seqan3::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [seqan3::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [seqan3::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [seqan3::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
 * , or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
 * another IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used
 * to further distinguish between merged bins.
 *
 * In this example, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single technical
 * bin, and another IBF was added for the merged bin.
 * \image html hibf.svg
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
//template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
template <typename FingerprintType = uint8_t>
//!\cond
        requires std::same_as<FingerprintType,uint8_t> || std::same_as<FingerprintType,uint16_t>
//!\endcond
class hierarchical_interleaved_binary_fuse_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;


    // Forward declaration
    template <std::integral value_t>
    class counting_agent_type;

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    //static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of an individual Bloom filter.
    using ixf_t = seqan3::interleaved_binary_fuse_filter<FingerprintType>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_binary_fuse_filter() = default;                                              //!< Defaulted.
    hierarchical_interleaved_binary_fuse_filter(hierarchical_interleaved_binary_fuse_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_binary_fuse_filter &
    operator=(hierarchical_interleaved_binary_fuse_filter const &) = default;                        //!< Defaulted.
    hierarchical_interleaved_binary_fuse_filter(hierarchical_interleaved_binary_fuse_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_binary_fuse_filter &
    operator=(hierarchical_interleaved_binary_fuse_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_binary_fuse_filter() = default;            //!< Defaulted.

    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector<ixf_t> ixf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ixf_id;

    //!\brief The underlying user bins.
    user_bins user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    //template <std::integral value_t = uint8_t>
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_binary_fuse_filter::membership_agent{*this};
    }
   

    /*!\brief Returns a counting_agent_type to be used for counting.
     * \tparam value_t The type to use for the counters; must model std::integral.
     */
    template <std::integral value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }


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
        archive(ixf_vector);
        archive(next_ixf_id);
        archive(user_bins);
    }
    //!\endcond
};

/*!\brief Bookkeeping for user and technical bins.
 */
//template <seqan3::data_layout data_layout_mode>
template <typename FingerprintType>
class hierarchical_interleaved_binary_fuse_filter<FingerprintType>::user_bins
{
private:
    //!\brief Contains filenames of all user bins.
    std::vector<std::string> user_bin_filenames;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_filename_position[i][b]`.
     * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
     * Otherwise, the returned value `j` can be used to access the corresponding filename `user_bin_filenames[j]`.
     */
    std::vector<std::vector<int64_t>> ixf_bin_to_filename_position{};

public:
/*
    user_bins() = default;                                              //!< Defaulted.
    user_bins(user_bins const &) = default; //!< Defaulted.
    user_bins &operator=(user_bins const &) = default;                        //!< Defaulted.
    user_bins(user_bins &&) = default; //!< Defaulted.
    user_bins &operator=(user_bins &&) = default; //!< Defaulted.
    ~user_bins() = default;           
*/
    //!\brief Returns the number of managed user bins.
    size_t num_user_bins() const noexcept
    {
        return user_bin_filenames.size();
    }

    //!\brief Changes the number of managed IBFs.
    void set_ixf_count(size_t const size)
    {
        ixf_bin_to_filename_position.resize(size);
    }

    //!\brief Changes the number of managed user bins.
    void set_user_bin_count(size_t const size)
    {
        user_bin_filenames.resize(size);
    }

    //!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
    std::vector<int64_t> & bin_indices_of_ixf(size_t const idx)
    {
        return ixf_bin_to_filename_position[idx];
    }

    //!\brief Returns the filename of the `idx`th user bin.
    std::string & filename_of_user_bin(size_t const idx)
    {
        return user_bin_filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return user_bin_filenames[ixf_bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`th IBF.
     *        An empty string is returned for merged bins.
     */
    auto operator[](size_t const ixf_idx) const
    {
        return ixf_bin_to_filename_position[ixf_idx]
             | std::views::transform(
                   [this](int64_t i)
                   {
                       if (i == -1)
                           return std::string{};
                       else
                           return user_bin_filenames[i];
                   });
    }

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index(size_t const ixf_idx, size_t const bin_idx) const
    {
        return ixf_bin_to_filename_position[ixf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
    template <typename stream_t>
    void write_filenames(stream_t & out_stream) const
    {
        size_t position{};
        std::string line{};
        for (auto const & filename : user_bin_filenames)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        archive(user_bin_filenames);
        archive(ixf_bin_to_filename_position);
    }
    //!\endcond
};

/*!\brief Manages membership queries for the hibf::hierarchical_interleaved_bloom_filter.
 * \see hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
 * \details
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 */
// TODO: value_t as template?
template <typename FingerprintType>
class hierarchical_interleaved_binary_fuse_filter<FingerprintType>::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hixf_t = hierarchical_interleaved_binary_fuse_filter<FingerprintType>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hixf_t const * const hixf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void bulk_contains_impl(value_range_t && values, int64_t const ixf_idx, size_t const threshold)
    {
        auto agent = hixf_ptr->ixf_vector[ixf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);
        uint16_t sum{};
        uint16_t max_sum{0};
       
        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hixf_ptr->user_bins.filename_index(ixf_idx, bin);

            //if (sum >= 50)
                

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_contains_impl(values, hixf_ptr->next_ixf_id[ixf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hixf_ptr->user_bins.filename_index(ixf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer.emplace_back(current_filename_index);
                sum = 0u;
            }

            if (sum > max_sum)
                max_sum = sum;
        }
        
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default;                                     //!< Defaulted.
    membership_agent(membership_agent const &) = default;             //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default;                  //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default;      //!< Defaulted.
    ~membership_agent() = default;                                    //!< Defaulted.

    /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent(hixf_t const & hixf) : hixf_ptr(std::addressof(hixf))
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector<int64_t> result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a user bin if there are at least this many hits.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
    {
        assert(hixf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        bulk_contains_impl(values, 0, threshold);

        std::ranges::sort(result_buffer); // TODO: necessary?

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                             size_t const threshold) && noexcept = delete;
    //!\}
};


/*!\brief Manages counting ranges of values for the hibf::hierarchical_interleaved_bloom_filter.
 */
//template <seqan3::data_layout data_layout_mode>
template <typename FingerprintType>
template <std::integral value_t>
class hierarchical_interleaved_binary_fuse_filter<FingerprintType>::counting_agent_type
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hixf_t = hierarchical_interleaved_binary_fuse_filter<FingerprintType>;
    typedef seqan3::interleaved_binary_fuse_filter<FingerprintType>::counting_agent_type TIXFAgent;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hixf_t const * const hixf_ptr{nullptr};

    //!\brief Helper for recursive bulk counting.
    template <std::ranges::forward_range value_range_t>
    void bulk_count_impl(value_range_t && values, int64_t const ixf_idx, size_t const threshold)
    {
        auto agent = hixf_ptr->ixf_vector[ixf_idx].template counting_agent<value_t>();
        auto & result = agent.bulk_count(values);

        value_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];
            auto const current_filename_index = hixf_ptr->user_bins.filename_index(ixf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_count_impl(values, hixf_ptr->next_ixf_id[ixf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hixf_ptr->user_bins.filename_index(ixf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer[current_filename_index] = sum;
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default;                                        //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default;             //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default;                  //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default;      //!< Defaulted.
    ~counting_agent_type() = default;                                       //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit counting_agent_type(hixf_t const & hixf) :
        hixf_ptr(std::addressof(hixf)),
        result_buffer(hixf_ptr->user_bins.num_user_bins())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    TIXFAgent::counting_vector result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::forward_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     * \param[in] threshold Do not recurse into merged bins with less than this many hits. Default: 1.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] TIXFAgent::counting_vector const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) & noexcept
    {
        assert(hixf_ptr != nullptr);
        assert(threshold > 0u);
        assert(result_buffer.size() == hixf_ptr->user_bins.num_user_bins());

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, static_cast<value_t>(0u));

        bulk_count_impl(values, 0, threshold);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] TIXFAgent::counting_vector const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) && noexcept = delete;
    //!\}
};

} // namespace hixf
