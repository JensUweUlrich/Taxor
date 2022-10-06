// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include "robin_hood.h"

#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ixf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root);

void insert_into_ixf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf);

} // namespace raptor::hibf
