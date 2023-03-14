
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_data.hpp"

namespace hixf
{

void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes);

void create_temp_hash_file(size_t const ixf_pos, size_t const bin_index, ankerl::unordered_dense::set<size_t> &node_hashes);

void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files);

void read_from_temp_hash_file(size_t const ixf_position,
                              uint16_t const bin_index,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files);
}