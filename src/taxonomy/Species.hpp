
#pragma once

#include <string>
#include <seqan3/core/concept/cereal.hpp>

namespace taxor::taxonomy
{

class Species
{

public:
	std::string organism_name;
	std::string accession_id; // ncbi accession
    std::string taxid;
    std::string file_stem;
    std::string taxnames_string;
    std::string taxid_string;
    uint16_t user_bin;
    uint64_t seq_len;


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
        archive(organism_name);
        archive(accession_id);
        archive(taxid);
        archive(taxnames_string);
        archive(taxid_string);
        archive(user_bin);
        archive(seq_len);
    }

};
}