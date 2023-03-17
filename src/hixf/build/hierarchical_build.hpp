
#pragma once

#include <ankerl/unordered_dense.h>

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(ankerl::unordered_dense::set<size_t> &parent_hashes,
                          lemon::ListDigraph::Node & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root,
                          bool is_second,
                          bool is_third);

} 
