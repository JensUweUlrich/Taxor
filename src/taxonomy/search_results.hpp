#pragma once

#include <string>

/*!\file search_results.hpp
 * \brief Defines taxor::taxonomy::Search_Result, the record type used to parse
 *        and hold one read-to-reference candidate match from a `taxor search`
 *        TSV output line, as consumed by `taxor profile`.
 */

namespace taxor::taxonomy
{

/*!\brief One read-to-reference candidate match, as parsed from a single line of
 *        `taxor search` TSV output.
 *
 * \note Sentinel convention: accession_id == "-" (with tax_id left unpopulated)
 *       marks a read that had no match, i.e. "unclassified". This convention is
 *       relied upon throughout the profiling code that consumes Search_Result
 *       (see profile_output.hpp's write_biobox_binning_file()), so it must be
 *       preserved by anything that produces or forwards Search_Result values.
 */
struct Search_Result
{
    //!\brief Identifier of the query read this result belongs to.
    std::string read_id;
    //!\brief Accession id of the matched reference genome, or "-" if unclassified.
    std::string accession_id;
    //!\brief NCBI taxid of the matched reference genome (unpopulated when unclassified).
    std::string tax_id;
    //!\brief Length of the matched reference genome/sequence.
    uint64_t ref_len;
    //!\brief Length of the query read.
    uint64_t query_len;
    //!\brief Total number of k-mer/syncmer hashes considered for the query.
    uint64_t query_hash_count;
    //!\brief Number of those hashes that matched the reference.
    uint64_t query_hash_match;

};

}