// --------------------------------------------------------------------------------------------------
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

/*!\file parse_chopper_pack_header.hpp
 * \brief Provides hixf::parse_chopper_pack_header.
 */

/*!\brief Parses the header of a chopper pack file and builds the IXF tree structure (nodes and edges) from it.
 * \param ixf_graph The (empty) graph to populate; one node is added per IXF (root plus one per merged bin).
 * \param node_map The node property map to populate with a hixf::node_data entry per node.
 * \param chopper_pack_file The already-opened chopper pack file stream, positioned at the start of the file;
 *                          on return, the stream is positioned at the first non-header (`#FILES` marker) line.
 * \return The number of merged-bin header records found (i.e. the number of non-root IXF nodes created).
 * \details
 *
 * The header consists of one line for the root/high-level IXF (prefixed with hixf::hixf_prefix) followed by
 * zero or more lines for merged bins (prefixed with hixf::merged_bin_prefix), each carrying a dot/underscore
 * separated chain of technical bin indices that locate it within the hierarchy plus its `max_bin_id`. Because
 * a merged bin's line may appear before or after its own parent's line in the file, all merged-bin header
 * records are first collected and then sorted by the length of their bin-index chain (i.e. by tree depth)
 * before being inserted, so that each record's parent node is guaranteed to already exist in `ixf_graph` when
 * the record is processed. For each new node, if it occupies the largest bin index (`max_bin_index`) of its
 * parent, the parent's `favourite_child` is updated to point at it.
 */
size_t parse_chopper_pack_header(lemon::ListDigraph & ixf_graph,
                                 lemon::ListDigraph::NodeMap<node_data> & node_map,
                                 std::istream & chopper_pack_file);

}
