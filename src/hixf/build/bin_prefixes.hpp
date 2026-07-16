// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <string_view>

namespace hixf
{

/*!\file bin_prefixes.hpp
 * \brief Provides the string literals used to tag header lines in a chopper pack file.
 * \details
 *
 * The chopper pack file format prefixes header comment lines with one of these tokens to indicate
 * whether the line describes the top-level (high-level) IXF or a merged bin at a deeper level. These
 * constants are used by hixf::parse_chopper_pack_header to recognise and split up the header lines.
 */

//!\brief Marks the header line describing the root/high-level IXF.
constexpr std::string_view hixf_prefix{"HIGH_LEVEL_IBF"};
//!\brief Marks a header line describing a merged bin (i.e. a bin that owns its own lower-level IXF).
constexpr std::string_view merged_bin_prefix{"MERGED_BIN"};
//!\brief Marks a bin that was split across multiple technical bins (currently unused by the header parser).
constexpr std::string_view split_bin_prefix{"SPLIT_BIN"};
//!\brief Cached length of #merged_bin_prefix, used to locate the bin-index substring in header lines.
constexpr size_t merged_bin_prefix_length{merged_bin_prefix.size()};

}
