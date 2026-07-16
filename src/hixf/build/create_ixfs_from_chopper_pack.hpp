// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

/*!\file create_ixfs_from_chopper_pack.hpp
 * \brief Provides hixf::create_ixfs_from_chopper_pack, the top-level entry point for building a HIXF.
 */

/*!\brief Builds a complete hierarchical_interleaved_xor_filter from a chopper pack layout file.
 * \param data The (empty) build state to populate; on return, `data.hixf` holds the finished filter.
 * \param arguments The build arguments, including the path to the chopper pack file (`arguments.bin_file`).
 * \details
 *
 * This is the top-level driver called by the `taxor build` subcommand. It reads and parses the chopper pack
 * file into the IXF node tree (hixf::read_chopper_pack_file), then recursively builds the filter for every
 * node starting at the root (hixf::hierarchical_build).
 */
//template <seqan3::data_layout data_layout_mode>
void create_ixfs_from_chopper_pack(build_data & data, build_arguments const & arguments);

}
