
#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include "strong_types.hpp"

namespace hixf
{

struct build_arguments
{
    // Related to k-mers
    uint8_t kmer_size{20u};
    uint32_t window_size{kmer_size};
    window window_size_strong{kmer_size};
    std::string shape_string{};
    seqan3::shape shape{seqan3::ungapped{kmer_size}};
    bool compute_syncmer{true};
    bool compute_minimiser{false};
    bool disable_cutoffs{false};
    uint8_t syncmer_size{10u};
    uint8_t t_syncmer{6u};
    uint16_t scaling{1u};

    // Related to IXF
    std::filesystem::path out_path{"./"};
    std::string size{"1k"};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    uint8_t parts{1u};
    double fpr{0.05};
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path bin_file{};
    uint8_t threads{1u};
    bool is_socks{false};
    bool is_hixf{false};
    bool is_minimiser{false};
    bool is_syncmer{false};
};

}
