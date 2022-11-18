
#include <lemon/list_graph.h> /// Must be first include.

#include "compute_hashes.hpp"
#include "construct_ixf.hpp"
#include "hierarchical_build.hpp"
#include "initialise_max_bin_hashes.hpp"
#include "insert_into_ixf.hpp"
#include "insert_into_bins.hpp"
#include "loop_over_children.hpp"
#include "update_user_bins.hpp"

namespace hixf
{
std::set<int64_t> investigate{};
robin_hood::unordered_flat_set<size_t> test_hashes{};
//template <seqan3::data_layout data_layout_mode>
size_t hierarchical_build(robin_hood::unordered_flat_set<size_t> & parent_hashes,
                          lemon::ListDigraph::Node const & current_node,
                          build_data & data,
                          build_arguments const & arguments,
                          bool is_root)
{
    auto & current_node_data = data.node_map[current_node];
    
    size_t const ixf_pos{data.request_ixf_idx()};

    std::vector<int64_t> ixf_positions(current_node_data.number_of_technical_bins, ixf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    std::vector<robin_hood::unordered_flat_set<size_t>> node_hashes{};
    // initialize hash sets of IXF bins
    for (int i = 0; i < current_node_data.number_of_technical_bins; ++i)
    {
        robin_hood::unordered_flat_set<size_t> bin_data{};
        node_hashes.emplace_back(bin_data);
    }

    // initialize lower level IXF
    // deprecated since wefirst create lower level IXFs and return hashes to higher levels before creating
    // higher level IXF
    /*
    size_t const max_bin_tbs =
        initialise_max_bin_hashes(hashes, ixf_positions, filename_indices, current_node, data, arguments);
    auto && ixf = construct_ixf(parent_hashes, hashes, max_bin_tbs, current_node, data, arguments, is_root);
    */
    //hashes.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ixf
    // does nothing on lowest level
    loop_over_children(node_hashes, ixf_positions, current_node, data, arguments, is_root);
    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    //size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    if (is_root)
        std::cout << "Number of hashes in Bin 0: " <<node_hashes[0].size() << std::endl;

    for (size_t i = 0; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];
        
        for (auto const & filename : record.filenames)
        {
            if (filename.compare("files.renamed/GCF_000839085.1_genomic.fna.gz") == 0)
            {
                seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
                investigate.emplace(ixf_pos);
                compute_hashes(test_hashes, arguments, record);
            }
        }

        if (is_root && record.number_of_bins.back() == 1) // no splitting needed
        {
            insert_into_bins(arguments, record, node_hashes);
        }
        else
        {
            robin_hood::unordered_flat_set<size_t> hashes{};
            compute_hashes(hashes, arguments, record);
            // only insert into hashes and parent_hashes and not into IXF directly
            insert_into_bins(hashes, node_hashes, record.number_of_bins.back(), record.bin_indices.back());
            for (size_t h : hashes)
                parent_hashes.insert(h);
            //parent_hashes.emplace(hashes);
            hashes.clear();
        }

        update_user_bins(data, filename_indices, record);
    }

    for (int64_t p : ixf_positions)
    {
        if (investigate.contains(p))
        {
            seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
            investigate.emplace(ixf_pos);
        }
    }

    // insert all hashes of all technical bins into newly created IXF
    auto && ixf = construct_ixf(node_hashes);

    if (is_root)
    {
        std::vector<size_t> c{};
        std::ranges::copy(test_hashes, std::back_inserter(c));
        typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;
        for (uint64_t p : investigate)
        {
            TIXFAgent ixf_count_agent = data.hixf.ixf_vector[p].counting_agent< uint64_t >();
		    auto result = ixf_count_agent.bulk_count(c);
            seqan3::debug_stream << "Index " << p << ": " << result << "\n";
        }
        TIXFAgent ixf_count_agent = ixf.counting_agent< uint64_t >();
		auto result = ixf_count_agent.bulk_count(c);
        seqan3::debug_stream << "Root result: " << result << "\n";
    }

    data.hixf.ixf_vector[ixf_pos] = std::move(ixf);
    data.hixf.next_ixf_id[ixf_pos] = std::move(ixf_positions);
    data.hixf.user_bins.bin_indices_of_ixf(ixf_pos) = std::move(filename_indices);

    return ixf_pos;
}

/*
template size_t hierarchical_build<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                      lemon::ListDigraph::Node const &,
                                                                      build_data<seqan3::data_layout::uncompressed> &,
                                                                      build_arguments const &,
                                                                      bool);

template size_t hierarchical_build<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                    lemon::ListDigraph::Node const &,
                                                                    build_data<seqan3::data_layout::compressed> &,
                                                                    build_arguments const &,
                                                                    bool);
*/
}
