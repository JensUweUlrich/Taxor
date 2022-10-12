
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

    // Related to thresholding
    double tau{0.9999};
    double threshold{std::numeric_limits<double>::quiet_NaN()};
    double p_max{0.15};
    double fpr{0.05};
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
    bool is_hibf{false};
    bool cache_thresholds{false};

    hixf::threshold_parameters make_threshold_parameters() const noexcept
    {
        return {.window_size{window_size},
                .shape{shape},
                .pattern_size{pattern_size},
                .errors{errors},
                .percentage{threshold},
                .p_max{p_max},
                .fpr{fpr},
                .tau{tau},
                .cache_thresholds{cache_thresholds},
                .output_directory{index_file.parent_path()}};
    }
};

} 
