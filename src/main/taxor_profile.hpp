#pragma once

#include <seqan3/argument_parser/all.hpp>

namespace taxor::profile
{

    struct Ref_Map_Info
    {
        uint64_t unique_assign_reads = 0;
        uint64_t all_assigned_reads = 0;
        // <taxid, nr_shared_read_mappings>
        std::map<std::string, size_t> associated_species{};
    };

    int execute(seqan3::argument_parser & parser);

}