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

constexpr std::string_view hixf_prefix{"HIGH_LEVEL_IBF"};
constexpr std::string_view merged_bin_prefix{"MERGED_BIN"};
constexpr std::string_view split_bin_prefix{"SPLIT_BIN"};
constexpr size_t merged_bin_prefix_length{merged_bin_prefix.size()};

} 
