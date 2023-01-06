// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
seqan3::interleaved_xor_filter<> construct_ixf(ankerl::unordered_dense::set<size_t> & parent_kmers,
                                                 ankerl::unordered_dense::set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data & data,
                                                 build_arguments const & arguments,
                                                 bool is_root);

seqan3::interleaved_xor_filter<> construct_ixf(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes,
                                               size_t const & current_node_ixf_position);

seqan3::interleaved_xor_filter<> construct_ixf(std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes);


} // namespace raptor::hibf
