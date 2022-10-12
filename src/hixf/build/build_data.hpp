#pragma once

#include <atomic>
#include <seqan3/std/new>

#include "node_data.hpp"
#include "hierarchical_interleaved_xor_filter.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
struct build_data
{
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> ixf_number{};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_t> user_bin_number{};

    size_t number_of_user_bins{};
    size_t number_of_ixfs{};

    lemon::ListDigraph ixf_graph{};
    lemon::ListDigraph::NodeMap<node_data> node_map{ixf_graph};

    hierarchical_interleaved_xor_filter<uint8_t> hixf{};
    std::vector<double> fp_correction{};

    size_t request_ixf_idx()
    {
        return std::atomic_fetch_add(&ixf_number, 1u);
    }

    size_t request_user_bin_idx()
    {
        return std::atomic_fetch_add(&user_bin_number, 1u);
    }

    void resize()
    {
        hixf.ixf_vector.resize(number_of_ixfs);
        hixf.user_bins.set_ixf_count(number_of_ixfs);
        hixf.user_bins.set_user_bin_count(number_of_user_bins);
        hixf.next_ixf_id.resize(number_of_ixfs);
    }

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
