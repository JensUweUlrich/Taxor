
#pragma once


#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>
#include <seqan3/core/concept/cereal.hpp>

#include "Species.hpp"

typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;

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

    std::vector<Species> species_vector{};

    void add_filter(seqan3::interleaved_xor_filter<>& ixf)
    {
        multi_filter.emplace_back(std::move(ixf));
    }

    void add_filter(seqan3::interleaved_xor_filter<>& ixf, const uint16_t index)
    {
        multi_filter.insert(multi_filter.begin() + index,  std::move(ixf));
    }

    void add_filter(const uint16_t index, uint64_t bins, uint64_t elements_per_bin)
    {
        multi_filter.insert(multi_filter.begin() + index,  std::move(seqan3::interleaved_xor_filter<>(bins, elements_per_bin)));
    }

    seqan3::interleaved_xor_filter<>& get(uint16_t index)
    {
        return std::ref(multi_filter[index]);
    }

    uint16_t count_single_filter()
    {
        return multi_filter.size();
    }

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