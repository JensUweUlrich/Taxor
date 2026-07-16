// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#pragma once

#include <filesystem>

#include <seqan3/search/kmer_index/shape.hpp>

/*!\file threshold_parameters.hpp
 * \brief Declares the parameter struct used to configure hixf::threshold::threshold, i.e. which statistical
 *        model is used to compute the match-count threshold for a search run and with which parameters.
 */
namespace hixf
{

/*!\brief Plain aggregate of parameters that determine which threshold model hixf::threshold::threshold
 *        selects and how it computes the match-count threshold.
 *
 * Populated from search_arguments::make_threshold_parameters(). Fields are grouped by which threshold_kind
 * they are relevant for (see hixf::threshold::threshold::threshold_kinds); unused fields for a given kind
 * are simply ignored.
 */
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
