// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#include <lemon/list_graph.h> /// Must be first include.

#include "compute_hashes.hpp"
#include "construct_ixf.hpp"
#include "hierarchical_build.hpp"
//#include "initialise_max_bin_hashes.hpp"
#include "insert_into_ixf.hpp"
#include "insert_into_bins.hpp"
#include "loop_over_children.hpp"
#include "update_user_bins.hpp"
#include "temp_hash_file.hpp"

namespace hixf
{
std::set<int64_t> investigate{};
//robin_hood::unordered_flat_set<size_t> test_hashes{};
ankerl::unordered_dense::set<size_t> test_hashes{};
//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(ankerl::unordered_dense::set<size_t> &parent_hashes,
                          lemon::ListDigraph::Node & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root,
                          bool is_second,
                          bool is_third)
{
    auto & current_node_data = data.node_map[current_node];
    
    size_t const ixf_pos{data.request_ixf_idx()};

    std::vector<int64_t> ixf_positions(current_node_data.number_of_technical_bins, ixf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    std::vector<ankerl::unordered_dense::set<size_t>> node_hashes{};
    // initialize hash sets of IXF bins
    for (int i = 0; i < current_node_data.number_of_technical_bins; ++i)
    {
        ankerl::unordered_dense::set<size_t> bin_data{};
        node_hashes.emplace_back(bin_data);
    }

    // parse all other children (merged bins) of the current ixf
    // does nothing on lowest level
    loop_over_children(node_hashes, ixf_positions, current_node, data, arguments, is_root, is_second);

    for (size_t i = 0; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        // debugging only
        /*for (auto const & filename : record.filenames)
        {
            if (filename.compare("files/GCF_000147815.2_ASM14781v3_genomic.fna.gz") == 0)
            {
                seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
                investigate.emplace(ixf_pos);
                compute_hashes(test_hashes, arguments, record);
            }
        }
        */
        // add single bin
        if ((is_root || is_second) && record.number_of_bins.back() == 1) // no splitting needed
        {
            //insert_into_bins(arguments, record, node_hashes);
            // compute hashes and create temp file
            auto const bin_index = static_cast<size_t>(record.bin_indices.back());
            ankerl::unordered_dense::set<size_t> hashes{};
            compute_hashes(hashes, arguments, record);

            create_temp_hash_file(ixf_pos, bin_index ,hashes);

            // update maximum bin size if this bin is largest
            if (hashes.size() > current_node_data.max_bin_hashes)
                current_node_data.max_bin_hashes = hashes.size();

            hashes.clear();
        }
        else
        {
            ankerl::unordered_dense::set<size_t> hashes{};
            compute_hashes(hashes, arguments, record);
            if (is_root || is_second)
            {
                size_t const chunk_size = hashes.size() / record.number_of_bins.back() + 1;
                size_t chunk_number{};
                // this is a split bin 
                // create one temp file per splitted bin
                for (auto chunk : hashes | seqan3::views::chunk(chunk_size))
                {
                    assert(chunk_number < record.number_of_bins.back());
                    size_t bin_idx = record.bin_indices.back() + chunk_number;
                    ++chunk_number;
                    ankerl::unordered_dense::set<size_t> tmp_hashes{};
                    for (size_t const value : chunk)
                    {
                        tmp_hashes.insert(value);
                    }
                    // update maximum bin size if this split bin is largest
                    if (tmp_hashes.size() > current_node_data.max_bin_hashes)
                        current_node_data.max_bin_hashes = tmp_hashes.size();

                    create_temp_hash_file(ixf_pos, bin_idx ,tmp_hashes);
                }
            }
            else
            {
                // only insert into hashes and parent_hashes and not into IXF directly
                insert_into_bins(hashes, node_hashes, record.number_of_bins.back(), record.bin_indices.back());
            }
            hashes.clear();
        }

        update_user_bins(data, filename_indices, record);
    }

    // fo debugging only
    /*for (int64_t p : ixf_positions)
    {
        if (investigate.contains(p))
        {
            //seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
            investigate.emplace(ixf_pos);
        }
    }
    */
    // store hashes in parent_hash_set 
    // only wite hashes to file if parent is root IXF
    bool low_mem = false;

    if (is_root || is_second)
    {
        
        auto && ixf = construct_ixf(data, current_node, ixf_positions, is_second, ixf_pos);

        // only for debugging
        /*std::vector<size_t> c{};
        std::ranges::copy(test_hashes, std::back_inserter(c));
        typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;  
        for (uint64_t p : investigate)
        {
            if (p == 0)
                continue;
            TIXFAgent ixf_count_agent = data.hixf.ixf_vector[p].counting_agent< uint64_t >();
		    auto result = ixf_count_agent.bulk_count(c);
            seqan3::debug_stream << "Index " << p << ": " << result << "\n";

            
        }
        */
        
        data.hixf.ixf_vector[ixf_pos] = std::move(ixf);
        data.hixf.next_ixf_id[ixf_pos] = std::move(ixf_positions);
        data.hixf.user_bins.bin_indices_of_ixf(ixf_pos) = std::move(filename_indices);

        // only for debugging
        /*TIXFAgent ixf_count_agent = data.hixf.ixf_vector[ixf_pos].counting_agent< uint64_t >();
		auto result = ixf_count_agent.bulk_count(c);
        seqan3::debug_stream << "Root Index " << ixf_pos << ": " << result << "\n";
        */
    }
    else
    {
        // for level below root, we store hashes of the IXF in a temp file
        // reduces peak memory
        if (is_third)
        //if (low_mem)
        {
            ankerl::unordered_dense::set<size_t> hashset{};
            for (auto hash_bin : node_hashes)
                for (size_t hash : hash_bin)
                    hashset.insert(hash);
            
            current_node_data.number_of_hashes = hashset.size();

            create_temp_hash_file(ixf_pos, hashset);
            hashset.clear();
        }
        else
        {
            for (auto &hash_bin : node_hashes)
            {
                for (size_t hash : hash_bin)
                    parent_hashes.insert(hash);
            }
        }

        if (low_mem)
        {
            for (auto & hashset : node_hashes)
            {
            //std::cout << hashset.size() << std::endl << std::flush;
                if (hashset.size() > current_node_data.max_bin_hashes)
                    current_node_data.max_bin_hashes = hashset.size();
            }
            //std::cout << current_node_data.max_bin_hashes << std::endl << std::flush;
            auto && ixf = construct_ixf(data, current_node, ixf_positions, is_second, ixf_pos);
            data.hixf.ixf_vector[ixf_pos] = std::move(ixf);
            data.hixf.next_ixf_id[ixf_pos] = std::move(ixf_positions);
            data.hixf.user_bins.bin_indices_of_ixf(ixf_pos) = std::move(filename_indices);
        }
        else
        {

            // insert all hashes of all technical bins into newly created IXF
            auto && ixf = construct_ixf(node_hashes);
            data.hixf.ixf_vector[ixf_pos] = std::move(ixf);
            data.hixf.next_ixf_id[ixf_pos] = std::move(ixf_positions);
            data.hixf.user_bins.bin_indices_of_ixf(ixf_pos) = std::move(filename_indices);

            /*std::vector<size_t> c{};
            std::ranges::copy(test_hashes, std::back_inserter(c));
            typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;  
            for (uint64_t p : investigate)
            {
                if (p == 0)
                    continue;
                TIXFAgent ixf_count_agent = data.hixf.ixf_vector[p].counting_agent< uint64_t >();
		        auto result = ixf_count_agent.bulk_count(c);
                seqan3::debug_stream << "Index " << p << ": " << result << "\n";
            }
            */
        }
        
    }

    node_hashes.clear();

    return ixf_pos;
}

}
