
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ixf(ankerl::unordered_dense::set<size_t> & parent_kmers,
                     ankerl::unordered_dense::set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root);

void insert_into_ixf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf);

} // namespace raptor::hibf
