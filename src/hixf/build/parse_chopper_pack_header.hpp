
#pragma once

#include <iosfwd>

#include "node_data.hpp"

namespace hixf
{

size_t parse_chopper_pack_header(lemon::ListDigraph & ixf_graph,
                                 lemon::ListDigraph::NodeMap<node_data> & node_map,
                                 std::istream & chopper_pack_file);

} 
