// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <ankerl/unordered_dense.h>

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

/*!\file hierarchical_build.hpp
 * \brief Provides hixf::hierarchical_build, the recursive core of the HIXF construction algorithm.
 */

/*!\brief Recursively builds the interleaved XOR filter for `current_node` and all of its descendants.
 * \param parent_hashes Output set that this node's hashes are inserted into, so the caller (parent node) can
 *                       include them in its own IXF; unused when `is_root` or `is_second` is true (those levels
 *                       instead persist their hashes via temp files, see hixf::create_temp_hash_file).
 * \param current_node The IXF-tree node to build (and recurse from) in this call.
 * \param data The shared build state; the finished filter for this node is written into `data.hixf`.
 * \param arguments The build arguments.
 * \param is_root Whether `current_node` is the root (top-level) IXF node.
 * \param is_second Whether `current_node` is one level below the root.
 * \param is_third Whether `current_node` is two levels below the root.
 * \return The index assigned to `current_node`'s filter within `data.hixf.ixf_vector`.
 * \details
 *
 * For each node, this function first recurses into all child nodes (merged bins that themselves hold a
 * lower-level IXF, via hixf::loop_over_children — in parallel at the root, sequentially further down), then
 * computes/hashes the sequences of any of this node's own (non-merged) user bins (hixf::compute_hashes),
 * and finally assembles this node's own interleaved XOR filter from all of that (hixf::construct_ixf).
 * The `is_root`/`is_second`/`is_third` flags select between three different memory-management strategies for
 * intermediate hash sets — kept fully in memory near the root, spilled to temp files one level down to bound
 * peak memory, and read back per-child further down the tree — since the whole hash-set tree does not fit in
 * memory at once for large inputs.
 */
//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(ankerl::unordered_dense::set<size_t> &parent_hashes,
                          lemon::ListDigraph::Node & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root,
                          bool is_second,
                          bool is_third);

}
