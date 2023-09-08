/ --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <iosfwd>

#include "node_data.hpp"

namespace hixf
{

size_t parse_chopper_pack_header(lemon::ListDigraph & ixf_graph,
                                 lemon::ListDigraph::NodeMap<node_data> & node_map,
                                 std::istream & chopper_pack_file);

} 
