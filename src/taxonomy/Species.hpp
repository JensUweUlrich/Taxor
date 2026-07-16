
#pragma once

#include <string>
#include <seqan3/core/concept/cereal.hpp>

/*!\file Species.hpp
 * \brief Defines taxor::taxonomy::Species, the per-reference-genome metadata
 *        record produced by the taxonomy parser and stored inside a built
 *        `taxor build` index.
 */

namespace taxor::taxonomy
{

/*!\brief Per-reference-genome metadata, one instance per row of the taxonomy
 *        input file consumed by `taxor build`.
 *
 * Instances are created by parse_refseq_taxonomy_file() and are serialised as
 * part of the built index (see CEREAL_SERIALIZE_FUNCTION_NAME() below), so
 * `taxor search`/`taxor profile` can recover the taxonomic lineage of a
 * reference genome that a query matched against.
 */
class Species
{

public:
	//!\brief Organism name as given in the taxonomy input file (may be empty).
	std::string organism_name;
	std::string accession_id; // ncbi accession
	//!\brief NCBI taxid of the reference genome, as a string.
    std::string taxid;
    //!\brief Final path component (file name) of the reference FASTA file, derived
    //!       from the taxonomy input file's path column.
    std::string file_stem;
    //!\brief GTDB-style lineage name string, e.g. "k__...;p__...;...;s__...".
    std::string taxnames_string;
    //!\brief GTDB-style lineage taxid string, aligned position-for-position with
    //!       taxnames_string (e.g. "1;2;...;1234").
    std::string taxid_string;
    //!\brief Index of the user bin this Species was placed into within the index.
    uint64_t user_bin;
    //!\brief Total sequence length of the reference genome.
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