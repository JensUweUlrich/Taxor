
#include <lemon/list_graph.h> /// Must be first include.

#include "compute_hashes.hpp"
#include "construct_ixf.hpp"
#include "hierarchical_build.hpp"
#include "initialise_max_bin_hashes.hpp"
#include "insert_into_ixf.hpp"
#include "loop_over_children.hpp"
#include "update_user_bins.hpp"

namespace hixf
{
std::set<int64_t> investigate{};
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
    robin_hood::unordered_flat_set<size_t> hashes{};

    // initialize lower level IXF
    size_t const max_bin_tbs =
        initialise_max_bin_hashes(hashes, ixf_positions, filename_indices, current_node, data, arguments);
    auto && ixf = construct_ixf(parent_hashes, hashes, max_bin_tbs, current_node, data, arguments, is_root);
    hashes.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    loop_over_children(parent_hashes, ixf, ixf_positions, current_node, data, arguments, is_root);
    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];
        
        for (auto const & filename : record.filenames)
        {
            if (filename.compare("files.renamed/GCF_000839085.1_genomic.fna.gz") == 0)
            {
                seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
                investigate.emplace(ixf_pos);
            }
        }

        if (is_root && record.number_of_bins.back() == 1) // no splitting needed
        {
            insert_into_ixf(arguments, record, ixf);
        }
        else
        {
            compute_hashes(hashes, arguments, record);
            insert_into_ixf(parent_hashes, hashes, record.number_of_bins.back(), record.bin_indices.back(), ixf, is_root);
        }

        update_user_bins(data, filename_indices, record);
        hashes.clear();
    }

    for (int64_t p : ixf_positions)
    {
        if (investigate.contains(p))
        {
            seqan3::debug_stream << "IXF vector index: " << ixf_pos << "\n";
            investigate.emplace(ixf_pos);
        }
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
