
#pragma once

#include <string>
#include <seqan3/core/concept/cereal.hpp>

class Species
{

public:
	std::string name;
	std::string organism_name;
	std::string accession_id;
	std::string taxonomy;
    uint16_t filter_index = UINT16_MAX;
    uint64_t first_bin = UINT64_MAX;
    uint64_t last_bin = UINT64_MAX;


    /*!\name Constructors, destructor and assignment
     */
    Species() = default; //!< Defaulted.
    Species(Species const &) = default; //!< Defaulted.
    Species & operator=(Species const &) = default; //!< Defaulted.
    Species(Species &&) = default; //!< Defaulted.
    Species & operator=(Species &&) = default; //!< Defaulted.
    ~Species() = default; //!< Defaulted.
	
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
        archive(name);
        archive(organism_name);
        archive(accession_id);
        archive(taxonomy);
        archive(filter_index);
        archive(first_bin);
        archive(last_bin);
    }

};