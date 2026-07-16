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

/*!\file construct_ixf.hpp
 * \brief Provides the hixf::construct_ixf overloads used to build an individual level's IXF.
 */

/*!\brief Builds a single-level interleaved XOR filter directly from an in-memory hash set, splitting it evenly
 *        across `number_of_bins` technical bins.
 * \param parent_kmers Output set that hashes are additionally inserted into when `is_root` is false, so that the
 *                      parent level's IXF also contains them.
 * \param kmers The hashes to insert into the new filter.
 * \param number_of_bins The number of technical bins to split `kmers` across.
 * \param node The IXF-tree node these hashes belong to.
 * \param data The shared build state (used for false-positive correction and node metadata).
 * \param arguments The build arguments.
 * \param is_root Whether `node` is the root of the IXF tree; if false, hashes are propagated to `parent_kmers`.
 * \return The newly constructed seqan3::interleaved_xor_filter.
 * \deprecated No current call sites in this codebase; superseded by the two overloads below, which build the
 *             filter from data already spilled to per-bin temp files.
 */
seqan3::interleaved_xor_filter<> construct_ixf(ankerl::unordered_dense::set<size_t> & parent_kmers,
                                                 ankerl::unordered_dense::set<size_t> & kmers,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data & data,
                                                 build_arguments const & arguments,
                                                 bool is_root);

/*!\brief Builds the interleaved XOR filter for `current_node`, combining hashes from child IXFs (read back from
 *        temporary hash files) with hashes of any bins computed directly at this level.
 * \param data The shared build state; used to look up node/child metadata and hash files.
 * \param current_node The IXF-tree node whose filter is being built.
 * \param ixf_positions For each technical bin of `current_node`, the temp-file index that holds its hashes
 *                       (bins with no dedicated child IXF use `current_node_ixf_position`).
 * \param is_second Whether `current_node` is one level below the root; if true, all inserted hashes are also
 *                  collected and written back out as a single per-node temp hash file for a later re-read
 *                  (see hixf::create_temp_hash_file), instead of being merged straight into the parent's memory.
 * \param current_node_ixf_position The temp-file index under which hashes computed directly at this level
 *                                  (not from a child IXF) were stored.
 * \return The newly constructed seqan3::interleaved_xor_filter.
 * \details
 *
 * Because seqan3::interleaved_xor_filter uses a randomized construction (XOR filters can fail to build for a
 * given seed/element set and must then be retried with a new seed), this function retries the whole bin-by-bin
 * insertion with a freshly reseeded filter (`ixf.set_seed()`) whenever `add_bin_elements` reports failure for
 * any bin, until all bins are inserted successfully.
 */
seqan3::interleaved_xor_filter<> construct_ixf(build_data & data,
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_position);

/*!\brief Builds an interleaved XOR filter directly from a vector of already-partitioned per-bin hash sets.
 * \param node_hashes One hash set per technical bin, in bin order.
 * \return The newly constructed seqan3::interleaved_xor_filter, sized and filled from `node_hashes`.
 */
seqan3::interleaved_xor_filter<> construct_ixf(std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes);


} // namespace raptor::hibf
