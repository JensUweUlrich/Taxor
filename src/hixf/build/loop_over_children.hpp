/ --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void loop_over_children(std::vector<ankerl::unordered_dense::set<size_t>> & parent_hashes,
                        std::vector<int64_t> & ixf_positions,
                        lemon::ListDigraph::Node & current_node,
                        build_data & data,
                        build_arguments const & arguments,
                        bool is_root,
                        bool is_second);

} 
