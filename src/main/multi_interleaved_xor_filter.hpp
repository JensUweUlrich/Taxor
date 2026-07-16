/*!\file multi_interleaved_xor_filter.hpp
 * \brief Container aggregating multiple seqan3::interleaved_xor_filter instances (each capped in bin count)
 *        plus the Species metadata they index, used by the alternate build.hpp indexing pipeline.
 *
 * \details Not currently compiled into the taxor binary — this file is only #include'd by build.hpp, which is
 *          itself not #include'd by anything reachable from main.cpp as of this writing, so this class is
 *          effectively dead code too.
 */

#pragma once


#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>
#include <seqan3/core/concept/cereal.hpp>

#include "Species.hpp"

//!\brief Counting agent type used to bulk-count query hashes against a single seqan3::interleaved_xor_filter.
typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;

/*!\brief Holds a collection of independent seqan3::interleaved_xor_filter instances (each one holding up to
 *        \c max_bins_per_filter bins, see build.hpp) together with the Species metadata describing which
 *        (filter, bin range) each species was assigned to.
 *
 * \details This class exists because a single interleaved_xor_filter cannot grow indefinitely (bin count is
 *          fixed at construction and very large bin counts become impractical), so the alternate build.hpp
 *          pipeline shards species across several separate filters and this class bundles them together for
 *          querying and (de)serialisation as one unit.
 */
class multi_interleaved_xor_filter
{
private:

    std::vector<seqan3::interleaved_xor_filter<>> multi_filter{};


public:
    /*!\name Constructors, destructor and assignment
     */
    multi_interleaved_xor_filter() = default; //!< Defaulted.
    multi_interleaved_xor_filter(multi_interleaved_xor_filter const &) = default; //!< Defaulted.
    multi_interleaved_xor_filter & operator=(multi_interleaved_xor_filter const &) = default; //!< Defaulted.
    multi_interleaved_xor_filter(multi_interleaved_xor_filter &&) = default; //!< Defaulted.
    multi_interleaved_xor_filter & operator=(multi_interleaved_xor_filter &&) = default; //!< Defaulted.
    ~multi_interleaved_xor_filter() = default; //!< Defaulted.

    //!\brief Species metadata, indexed in the same order as filter/bin assignment; see Species.hpp for the
    //!       filter_index/first_bin/last_bin members that tie a species to its location within \c multi_filter.
    std::vector<Species> species_vector{};

    //!\brief Append an already-constructed filter to the end of the collection.
    //!\param ixf The filter to move in; left in a moved-from state afterwards.
    //!\note No call sites of this overload were found anywhere in the source tree; build.hpp uses the
    //!      3-argument overload below instead.
    void add_filter(seqan3::interleaved_xor_filter<>& ixf)
    {
        multi_filter.emplace_back(std::move(ixf));
    }

    //!\brief Insert an already-constructed filter at a specific position in the collection.
    //!\param ixf   The filter to move in; left in a moved-from state afterwards.
    //!\param index Position to insert at.
    //!\note No call sites of this overload were found anywhere in the source tree.
    void add_filter(seqan3::interleaved_xor_filter<>& ixf, const uint16_t index)
    {
        multi_filter.insert(multi_filter.begin() + index,  std::move(ixf));
    }

    //!\brief Construct a new, empty seqan3::interleaved_xor_filter and insert it at position \p index.
    //!\param index            Position to insert the new filter at.
    //!\param bins             Number of bins the new filter should have.
    //!\param elements_per_bin Expected maximum number of elements per bin, used to size the filter.
    void add_filter(const uint16_t index, uint64_t bins, uint64_t elements_per_bin)
    {
        multi_filter.insert(multi_filter.begin() + index,  std::move(seqan3::interleaved_xor_filter<>(bins, elements_per_bin)));
    }

    //!\brief Access the filter at \p index by reference (e.g. to populate its bins).
    //!\param index Position of the filter to access.
    //!\return Reference to the requested seqan3::interleaved_xor_filter.
    seqan3::interleaved_xor_filter<>& get(uint16_t index)
    {
        return std::ref(multi_filter[index]);
    }

    //!\brief Number of individual interleaved_xor_filter instances held by this container.
    uint16_t count_single_filter()
    {
        return multi_filter.size();
    }

    /*!\brief Count, per underlying filter, how many of \p read_hashes hit each bin.
     * \param read_hashes Query k-mer/syncmer hashes to count against every bin of every filter.
     * \return One counting_vector per underlying filter (same order as the filters), each holding one count
     *         per bin of that filter.
     */
    std::vector<TIXFAgent::counting_vector> bulk_count(const std::vector<uint64_t>& read_hashes )
    {
        std::vector<TIXFAgent::counting_vector> count_vectors;
        for (auto & ixf : multi_filter)
        {
            TIXFAgent ixf_count_agent = ixf.counting_agent< uint64_t >();
			auto ixf_result = ixf_count_agent.bulk_count(read_hashes);
            count_vectors.push_back(std::move(ixf_result));
        }
        return std::move(count_vectors);
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
        archive(multi_filter);
        archive(species_vector);
    }
};