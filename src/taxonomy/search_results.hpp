#pragma once

#include <string>

namespace taxor::taxonomy
{

struct Search_Result
{
    std::string read_id;
    std::string accession_id;
    std::string tax_id;
    uint64_t ref_len;
    uint64_t query_len;
    uint64_t query_hash_count;
    uint64_t query_hash_match;

};

}