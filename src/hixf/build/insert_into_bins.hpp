#pragma once

#include "robin_hood.h"
#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_bins(robin_hood::unordered_flat_set<size_t> const & hashes,
                      std::vector<robin_hood::unordered_flat_set<size_t>> & ixf_bins,
                     size_t const number_of_bins,
                     size_t const bin_index);

void insert_into_bins(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     std::vector<robin_hood::unordered_flat_set<size_t>> & ixf_bins);

}