// --------------------------------------------------------------------------------------------------
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

/*!\file loop_over_children.hpp
 * \brief Provides hixf::loop_over_children, which recursively builds all merged-bin child IXFs of a node.
 */

/*!\brief Recursively builds the IXF for every merged-bin child of `current_node`, in parallel when building the
 *        root level.
 * \param parent_hashes Per-bin hash sets of `current_node`; hashes returned from non-root, non-second-level
 *                       children are inserted into the bin matching that child's `parent_bin_index`.
 * \param ixf_positions Per-bin IXF/temp-file positions of `current_node`; entries for merged-bin children are
 *                       overwritten with the index/position returned from hixf::hierarchical_build for that child.
 * \param current_node The IXF-tree node whose children should be built.
 * \param data The shared build state.
 * \param arguments The build arguments (used for thread count when building the root level).
 * \param is_root Whether `current_node` is the root IXF node.
 * \param is_second Whether `current_node` is one level below the root.
 * \details
 *
 * Does nothing if `current_node` has no children (i.e. it is a leaf of the IXF tree). Otherwise, each child is
 * built via a recursive call to hixf::hierarchical_build. At the root level, children are shuffled and built
 * concurrently across `arguments.threads` threads (shuffling makes it less likely that threads contend for the
 * same underlying resources at the same time); at all other levels, children are built with a single thread,
 * since parallelism is already exploited across the top-level children.
 *
 * ### Thread safety
 *
 * Concurrent child builds share `current_node`'s data via `local_ixf_mutex` (bucketed by
 * `parent_bin_index / 64`, since technical bins map onto 64-bit words in the underlying filter) and a single
 * `max_bin_hashes_mutex` guarding the running maximum bin size.
 */
//template <seqan3::data_layout data_layout_mode>
void loop_over_children(std::vector<ankerl::unordered_dense::set<size_t>> & parent_hashes,
                        std::vector<int64_t> & ixf_positions,
                        lemon::ListDigraph::Node & current_node,
                        build_data & data,
                        build_arguments const & arguments,
                        bool is_root,
                        bool is_second);

}
