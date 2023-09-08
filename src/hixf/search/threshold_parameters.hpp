// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


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
    uint8_t kmer_size{};
    uint64_t pattern_size{};

    // Threshold.
    uint8_t errors{};                                            // threshold_kinds::(probabilistic|lemma)
    double percentage{-1.0}; // threshold_kinds::percentage
    double p_max{};                                              // threshold_kinds::probabilistic
    double fpr{};                                                // threshold_kinds::probabilistic
    double tau{};                                                // threshold_kinds::probabilistic
    double seq_error_rate{};                                     // threshold_kinds::confidence_interval
    bool fracminhash{};
    bool use_syncmer{};

    // Cache results.
    bool cache_thresholds{};
    std::filesystem::path output_directory{};
};

} 
