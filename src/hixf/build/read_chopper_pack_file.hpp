// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------
#pragma once

#include <string>

#include "build_data.hpp"

namespace hixf
{

/*!\file read_chopper_pack_file.hpp
 * \brief Provides hixf::read_chopper_pack_file.
 */

/*!\brief Reads and fully parses a chopper pack file into `data`'s IXF tree.
 * \param data The build state to populate; on return, `data.ixf_graph`/`data.node_map` describe the full IXF
 *             tree (technical bin counts included) and `data` has been resized (see build_data::resize) to
 *             hold all resulting IXFs and user bins.
 * \param chopper_pack_filename Path to the chopper pack file to read.
 * \throws std::logic_error if the file cannot be opened for reading.
 * \details
 *
 * First delegates to hixf::parse_chopper_pack_header to build the tree of IXF nodes from the file's header.
 * Then parses each remaining data line (hixf::parse_chopper_pack_line) and, for each resulting record, walks
 * down the tree following the record's chain of bin indices to find which node it belongs to, growing that
 * node's `number_of_technical_bins` as needed and appending the record to that node's `remaining_records`.
 */
//template <seqan3::data_layout data_layout_mode>
void read_chopper_pack_file(build_data & data, std::string const & chopper_pack_filename);

}
