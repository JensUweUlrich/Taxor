#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_bins(ankerl::unordered_dense::set<size_t> const & hashes,
                      std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins,
                     size_t const number_of_bins,
                     size_t const bin_index);

void insert_into_bins(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins);

}