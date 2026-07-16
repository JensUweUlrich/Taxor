#pragma once

#include <seqan3/argument_parser/all.hpp>

/*!\file taxor_profile.hpp
 * \brief Public interface of the `taxor profile` subcommand; see taxor_profile.cpp for the implementation.
 */

namespace taxor::profile
{

    /*!\brief Per-reference bookkeeping used by filter_ref_associations to detect near-identical
     * strain references that "explain" each other's reads (see taxor_profile.cpp).
     */
    struct Ref_Map_Info
    {
        //!\brief Number of reads that map uniquely (with no other candidate) to this reference.
        uint64_t unique_assign_reads = 0;
        //!\brief Total number of reads (unique + ambiguous) that map to this reference.
        uint64_t all_assigned_reads = 0;
        //!\brief Other accession ids this reference co-occurs with in ambiguous reads, and how many reads they share.
        std::map<std::string, size_t> associated_species{};
    };

    //!\brief Entry point for the `taxor profile` subcommand (see taxor_profile.cpp for details).
    int execute(seqan3::argument_parser & parser);

}
