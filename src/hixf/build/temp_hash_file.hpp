
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_data.hpp"

namespace hixf
{

void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes);

void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes);
}