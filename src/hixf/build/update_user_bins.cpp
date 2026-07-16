// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> // Must be first include.

#include "update_user_bins.hpp"

namespace hixf
{

/*!\file update_user_bins.cpp
 * \brief Implements hixf::update_user_bins.
 */

//!\brief See update_user_bins.hpp for documentation.
//template <seqan3::data_layout data_layout_mode>
void update_user_bins(build_data & data,
                      std::vector<int64_t> & filename_indices,
                      chopper_pack_record const & record)
{
    size_t const idx = data.request_user_bin_idx();

    std::string & user_bin_filenames = data.hixf.user_bins.filename_of_user_bin(idx);
    for (auto const & filename : record.filenames)
    {
        user_bin_filenames += filename;
        user_bin_filenames += ';';
    }
    assert(!user_bin_filenames.empty());
    // Every filename was appended with a trailing ';'; drop the last one so the stored string uses ';' purely
    // as a separator (matching the format parse_chopper_pack_line expects when reading it back).
    user_bin_filenames.pop_back();

    std::fill_n(filename_indices.begin() + record.bin_indices.back(), record.number_of_bins.back(), idx);
}

} 
