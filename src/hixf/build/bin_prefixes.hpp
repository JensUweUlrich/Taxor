
#pragma once

#include <string_view>

namespace hixf
{

constexpr std::string_view hixf_prefix{"HIGH_LEVEL_IBF"};
constexpr std::string_view merged_bin_prefix{"MERGED_BIN"};
constexpr std::string_view split_bin_prefix{"SPLIT_BIN"};
constexpr size_t merged_bin_prefix_length{merged_bin_prefix.size()};

} 
