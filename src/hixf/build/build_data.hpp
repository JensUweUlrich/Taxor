// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <atomic>
#include <seqan3/std/new>

#include "node_data.hpp"
#include "hierarchical_interleaved_xor_filter.hpp"

namespace hixf
{

/*!\file build_data.hpp
 * \brief Provides hixf::build_data.
 */

/*!\brief Mutable, shared state used while building a hierarchical_interleaved_xor_filter from a chopper pack file.
 * \details
 *
 * A single build_data instance is threaded through the whole build process (hixf::read_chopper_pack_file,
 * hixf::hierarchical_build, hixf::loop_over_children, ...). It owns the tree of IXF nodes (#ixf_graph and
 * #node_map, built from the chopper pack header), the resulting hierarchical_interleaved_xor_filter (#hixf)
 * that is filled in as the tree is processed, and atomic counters (#ixf_number, #user_bin_number) used to
 * hand out unique IXF/user-bin indices to worker threads running in parallel over sibling subtrees.
 */
struct build_data
{
    //!\brief Atomically-issued counter used to assign each IXF node a unique index into hixf::ixf_vector.
    //!\details Cache-line aligned (padded to hardware_destructive_interference_size) together with
    //!         #user_bin_number to avoid false sharing between threads that increment the two counters
    //!         concurrently while building sibling subtrees.
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> ixf_number{};
    //!\brief Atomically-issued counter used to assign each user bin a unique index.
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> user_bin_number{};

    size_t number_of_user_bins{}; //!< Total number of user bins found in the chopper pack file.
    size_t number_of_ixfs{};      //!< Total number of IXF nodes (root + merged-bin subtrees) to be built.

    lemon::ListDigraph ixf_graph{}; //!< The tree of IXF nodes; edges point from a merged bin to its child IXF.
    lemon::ListDigraph::NodeMap<node_data> node_map{ixf_graph}; //!< Per-node metadata, keyed by #ixf_graph node.

    hierarchical_interleaved_xor_filter<uint8_t> hixf{}; //!< The filter being constructed.
    std::vector<double> fp_correction{}; //!< Per-merged-bin-count false-positive correction factors.

    //!\brief Atomically reserves and returns the next free IXF index.
    //!\return A unique index into hixf::ixf_vector, safe to call concurrently from multiple threads.
    size_t request_ixf_idx()
    {
        return std::atomic_fetch_add(&ixf_number, 1u);
    }

    //!\brief Atomically reserves and returns the next free user bin index.
    //!\return A unique index into hixf::user_bins, safe to call concurrently from multiple threads.
    size_t request_user_bin_idx()
    {
        return std::atomic_fetch_add(&user_bin_number, 1u);
    }

    //!\brief Resizes #hixf's internal vectors to hold #number_of_ixfs IXFs and #number_of_user_bins user bins.
    //!\details Must be called once #number_of_ixfs and #number_of_user_bins are known (i.e. after parsing the
    //!         chopper pack header) and before any IXF/user-bin index is dereferenced.
    void resize()
    {
        hixf.ixf_vector.resize(number_of_ixfs);
        hixf.user_bins.set_ixf_count(number_of_ixfs);
        hixf.user_bins.set_user_bin_count(number_of_user_bins);
        hixf.next_ixf_id.resize(number_of_ixfs);
    }

    /*!\brief Precomputes false-positive-rate correction factors for merged bins containing up to `tmax` user bins.
     * \param tmax The maximum number of user bins that may be merged into a single technical bin.
     * \param hash The number of hash functions used by the filter.
     * \param fpr The target false-positive rate of an individual technical bin.
     * \details
     *
     * When multiple user bins are merged into one technical bin, a positive hit on that bin does not by itself
     * tell us which of the merged user bins actually matched; the effective false-positive rate against any one
     * of them is higher than for a bin holding a single user bin. #fp_correction[i] stores the multiplicative
     * factor by which a per-bin hit threshold should be scaled up when i user bins share a technical bin, so
     * that comparisons against a threshold remain calibrated to the same target FPR regardless of merge depth.
     */
    void compute_fp_correction(size_t const tmax, size_t const hash, double const fpr)
    {
        fp_correction.resize(tmax + 1, 1.0);

        double const denominator = std::log(1 - std::exp(std::log(fpr) / hash));

        for (size_t i = 2; i <= tmax; ++i)
        {
            double const tmp = 1.0 - std::pow(1 - fpr, static_cast<double>(i));
            fp_correction[i] = std::log(1 - std::exp(std::log(tmp) / hash)) / denominator;
            assert(fp_correction[i] >= 1.0);
        }
    }
};

}
