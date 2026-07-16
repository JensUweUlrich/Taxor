// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>

#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file parse_chopper_pack_line.hpp
 * \brief Provides hixf::parse_chopper_pack_line.
 */

/*!\brief Parses a single tab-separated data line of a chopper pack file into a hixf::chopper_pack_record.
 * \param current_line One non-header line of the chopper pack file.
 * \return The parsed record, with `filenames`, `bin_indices`, and `number_of_bins` populated
 *         (`estimated_sizes` is left empty; this format does not carry per-record size estimates).
 * \details
 *
 * The expected line format is `filename1;filename2;...\tbin_idx_0;bin_idx_1;...\tnum_bins_0;num_bins_1;...`,
 * i.e. a `;`-separated list of filenames (one user bin, possibly composed of several sequence files), followed
 * by tab-separated, `;`-separated chains of technical bin indices and technical bin counts describing this
 * record's placement at each level of the IXF hierarchy. Any filename that is a symlink is resolved to its
 * target path before being stored.
 */
chopper_pack_record parse_chopper_pack_line(std::string const & current_line);

}
