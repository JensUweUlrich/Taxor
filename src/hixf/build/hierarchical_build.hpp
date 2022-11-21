
#pragma once

#include "robin_hood.h"

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(lemon::ListDigraph::Node const & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root);

} 
