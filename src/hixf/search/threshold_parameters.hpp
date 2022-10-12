
#pragma once

#include <filesystem>

#include <seqan3/search/kmer_index/shape.hpp>

namespace hixf
{

struct threshold_parameters
{
    // Basic.
    uint32_t window_size{};
    seqan3::shape shape{};
    uint64_t pattern_size{};

    // Threshold.
    uint8_t errors{};                                            // threshold_kinds::(probabilistic|lemma)
    double percentage{std::numeric_limits<double>::quiet_NaN()}; // threshold_kinds::percentage
    double p_max{};                                              // threshold_kinds::probabilistic
    double fpr{};                                                // threshold_kinds::probabilistic
    double tau{};                                                // threshold_kinds::probabilistic

    // Cache results.
    bool cache_thresholds{};
    std::filesystem::path output_directory{};
};

} 
