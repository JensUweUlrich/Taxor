// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <tuple>
#include <vector>

namespace hixf
{

/*!\file chopper_pack_record.hpp
 * \brief Provides hixf::chopper_pack_record.
 */

/*!\brief The parsed representation of a single non-header line of a chopper pack file.
 * \details
 *
 * A chopper pack file line describes where in the IXF hierarchy one user bin (one or more reference sequence
 * files that were grouped together upstream by `chopper`) is placed. Since a user bin's data can be split
 * across several technical bins (to balance load), and its placement can be nested inside merged bins at
 * several IXF levels, most fields are parallel vectors: entry `i` in each vector describes the technical bin
 * assignment at hierarchy level `i`, with the last entry (`.back()`) describing the level at which the record
 * ultimately lives (see hixf::parse_chopper_pack_line).
 */
struct chopper_pack_record
{
    std::vector<std::string> filenames{};        //!< The reference sequence file(s) that make up this user bin.
    std::vector<size_t> bin_indices{};            //!< Per-level starting technical bin index for this record.
    std::vector<size_t> number_of_bins{};         //!< Per-level number of technical bins this record occupies (>1 means split).
    std::vector<size_t> estimated_sizes{};        //!< Per-level estimated sizes (currently unused/unpopulated by the parser).

    //!\brief Compares all member vectors for equality.
    bool operator==(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes)
            == std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }

    //!\brief The negation of operator==.
    bool operator!=(chopper_pack_record const & other) const
    {
        return std::tie(filenames, bin_indices, number_of_bins, estimated_sizes)
            != std::tie(other.filenames, other.bin_indices, other.number_of_bins, other.estimated_sizes);
    }
};

}
