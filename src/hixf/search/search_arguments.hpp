
#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include "threshold_parameters.hpp"

namespace hixf
{

// For costum default message in argparser
struct pattern_size
{
    uint64_t v{};
};

struct search_arguments
{
    // Related to k-mers
    uint32_t window_size{20u};
    seqan3::shape shape{seqan3::ungapped{20u}};
    uint8_t shape_size{shape.size()};
    uint8_t shape_weight{shape.count()};
    uint8_t threads{1u};
    uint8_t parts{1u};
    bool compute_syncmer{true};

    // Related to thresholding
    double tau{0.9999};
    double threshold{-1.0};
    double p_max{0.15};
    double fpr{0.05};
    double seq_error_rate{0.04};
    uint64_t pattern_size{};
    hixf::pattern_size pattern_size_strong{};
    uint8_t errors{0};

    // Related to IBF
    std::filesystem::path index_file{};
    bool compressed{false};

    // General arguments
    std::vector<std::vector<std::string>> bin_path{};
    std::filesystem::path query_file{};
    std::filesystem::path out_file{"search.out"};
    bool write_time{false};
    bool is_socks{false};
    bool is_hixf{true};
    bool cache_thresholds{false};

    hixf::threshold_parameters make_threshold_parameters() noexcept
    {
        return hixf::threshold_parameters{window_size,
                shape,
                shape_size,
                pattern_size,
                errors,
                threshold,
                p_max,
                fpr,
                tau,
                seq_error_rate,
                cache_thresholds,
                index_file.parent_path()};
    }
};

} 
