
#pragma once

#include "robin_hood.h"

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> parent_hashes,
                          lemon::ListDigraph::Node const & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root,
                          bool parent_is_root);

} 
