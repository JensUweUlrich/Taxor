// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>

#include "hierarchical_build.hpp"
#include "insert_into_bins.hpp"
#include "loop_over_children.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void loop_over_children(std::vector<ankerl::unordered_dense::set<size_t>> & parent_hashes,
                        std::vector<int64_t> & ixf_positions,
                        lemon::ListDigraph::Node &current_node,
                        build_data & data,
                        build_arguments const & arguments,
                        bool is_root,
                        bool is_second)
{
    auto & current_node_data = data.node_map[current_node];
    std::vector<lemon::ListDigraph::Node> children{};

    for (lemon::ListDigraph::OutArcIt arc_it(data.ixf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
        children.emplace_back(data.ixf_graph.target(arc_it));

    if (children.empty())
        return;

    size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
    std::vector<std::mutex> local_ixf_mutex(number_of_mutex);
    bool is_third{false};
    if (is_second)
        is_third = true;
    
    auto worker = [&](auto && index, auto &&)
    {
        auto & child = children[index];
        ankerl::unordered_dense::set<size_t> hashes{};
        size_t const ixf_pos = hierarchical_build(hashes, child, data, arguments, false, is_root, is_third); // gets back 793 for low level IXF
        auto parent_bin_index = data.node_map[child].parent_bin_index;
        {
            size_t const mutex_id{parent_bin_index / 64};
            std::lock_guard<std::mutex> guard{local_ixf_mutex[mutex_id]};
            ixf_positions[parent_bin_index] = ixf_pos;
            // insert into parent_hashes if current node is not root
            // parent_hashes of level under root need to be stored on disk
            // to reduce peak memory
            if (!is_root && !is_second)
                insert_into_bins(hashes, parent_hashes, 1, parent_bin_index);
            
            // number of hashes stored in child IXF equals number of hashes in one corresponding bin
            // of current node => maximum bin size for current node's IXF needs to be known before
            // constructing the IXF at the current level
            auto &child_node_data = data.node_map[child];
            if (child_node_data.number_of_hashes > current_node_data.max_bin_hashes)
                current_node_data.max_bin_hashes = child_node_data.number_of_hashes;

        }
        
    };

    size_t number_of_threads{};
    auto indices_view = std::views::iota(0u, children.size()) | std::views::common;
    std::vector<size_t> indices{indices_view.begin(), indices_view.end()};

    if (is_root)
    {
        // Shuffle indices: More likely to not block each other. Optimal: Interleave
        std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
        number_of_threads = arguments.threads;
    }
    else
    {
        number_of_threads = 1u;
    }

    seqan3::detail::execution_handler_parallel executioner{number_of_threads};
    executioner.bulk_execute(std::move(worker), std::move(indices), []() {});
    executioner.wait();
    // insert all parent hashes before leaving the method
}

/*
template void loop_over_children<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                    seqan3::interleaved_bloom_filter<> &,
                                                                    std::vector<int64_t> &,
                                                                    lemon::ListDigraph::Node const &,
                                                                    build_data<seqan3::data_layout::uncompressed> &,
                                                                    build_arguments const &,
                                                                    bool);

template void loop_over_children<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                  seqan3::interleaved_bloom_filter<> &,
                                                                  std::vector<int64_t> &,
                                                                  lemon::ListDigraph::Node const &,
                                                                  build_data<seqan3::data_layout::compressed> &,
                                                                  build_arguments const &,
                                                                  bool);
*/
} 
