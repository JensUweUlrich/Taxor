// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------
#pragma once

#include "build_data.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file update_user_bins.hpp
 * \brief Provides hixf::update_user_bins.
 */

/*!\brief Registers a chopper pack record as a new user bin and records which technical bins it occupies.
 * \param data The shared build state; a new user bin index is reserved from `data` and its filename(s) stored
 *             in `data.hixf.user_bins`.
 * \param filename_indices Per-technical-bin user-bin index array for the current IXF node; the entries
 *                          corresponding to `record`'s technical bins are set to the newly reserved user bin
 *                          index.
 * \param record The chopper pack record describing the user bin (its reference filenames and the technical
 *               bins it occupies at the current hierarchy level).
 * \details
 *
 * A user bin may consist of several reference sequence files (when `chopper` grouped them together upstream);
 * all filenames are joined into a single `;`-separated string, matching the format used when writing/parsing
 * chopper pack lines (see hixf::parse_chopper_pack_line).
 */
//template <seqan3::data_layout data_layout_mode>
void update_user_bins(build_data & data,
                      std::vector<int64_t> & filename_indices,
                      chopper_pack_record const & record);

}
