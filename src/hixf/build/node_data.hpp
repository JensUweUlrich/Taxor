// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <lemon/list_graph.h>

#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file node_data.hpp
 * \brief Provides hixf::node_data.
 */

/*!\brief Per-node metadata attached to each node of the IXF tree (`build_data::ixf_graph` /
 *        `build_data::node_map`) while parsing the chopper pack file and building the filter.
 * \details
 *
 * One node_data exists per IXF in the hierarchy (the root/high-level IXF and one per merged bin that owns its
 * own lower-level IXF). It is first populated while parsing the chopper pack header and lines
 * (hixf::parse_chopper_pack_header, hixf::read_chopper_pack_file) and then consulted/updated while the actual
 * filter is built (hixf::hierarchical_build, hixf::construct_ixf).
 */
struct node_data // rename:ibf_data? or ibf_node_data
{
    size_t parent_bin_index{};    //!< The technical bin index of this node within its parent's IXF (0 for the root).
    size_t max_bin_index{};       //!< The largest technical bin index used by this node's IXF (from the header's `max_bin_id`).
    size_t number_of_technical_bins{}; //!< The total number of technical bins in this node's IXF.
    size_t max_bin_hashes{};      //!< The largest number of hashes stored in any single technical bin of this node's IXF.
    size_t number_of_hashes{};    //!< The total number of distinct hashes stored across this node's whole IXF.
    lemon::ListDigraph::Node favourite_child{lemon::INVALID}; //!< The child node occupying the bin at #max_bin_index, if any.
    std::vector<chopper_pack_record> remaining_records{}; // non-merged bins (either split or single)

    //!\brief Compares parent_bin_index, max_bin_index, number_of_technical_bins, favourite_child, and
    //!       remaining_records for equality.
    bool operator==(node_data const & rhs) const
    {
        bool res = std::tie(parent_bin_index, max_bin_index, number_of_technical_bins, favourite_child)
                == std::tie(rhs.parent_bin_index, rhs.max_bin_index, rhs.number_of_technical_bins, rhs.favourite_child);

        if (remaining_records.size() != rhs.remaining_records.size())
            return false;

        for (size_t i = 0; i < remaining_records.size(); ++i)
            res &= (remaining_records[i] == rhs.remaining_records[i]);

        return res;
    }

    //!\brief The negation of operator==.
    bool operator!=(node_data const & rhs) const
    {
        return !(*this == rhs);
    }
};

} 
