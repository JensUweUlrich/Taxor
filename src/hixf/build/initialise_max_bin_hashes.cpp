
#include <lemon/list_graph.h> /// Must be first include.

#include "compute_hashes.hpp"
#include "hierarchical_build.hpp"
#include "initialise_max_bin_hashes.hpp"
#include "update_user_bins.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
size_t initialise_max_bin_hashes(robin_hood::unordered_flat_set<size_t> & hashes,
                                std::vector<int64_t> & ixf_positions,
                                std::vector<int64_t> & filename_indices,
                                lemon::ListDigraph::Node const & node,
                                build_data & data,
                                build_arguments const & arguments)
{
    auto & node_data = data.node_map[node];

    if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
    {
        // recursively initialize favourite child first
        ixf_positions[node_data.max_bin_index] =
            hierarchical_build(hashes, node_data.favourite_child, data, arguments, false);
        return 1;
    }
    else // max bin is not a merged bin
    {
        // we assume that the max record is at the beginning of the list of remaining records.
        auto const & record = node_data.remaining_records[0];
        compute_hashes(hashes, arguments, record);
        update_user_bins(data, filename_indices, record);

        return record.number_of_bins.back();
    }
}
/*
template size_t
initialise_max_bin_kmers<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                            std::vector<int64_t> &,
                                                            std::vector<int64_t> &,
                                                            lemon::ListDigraph::Node const &,
                                                            build_data<seqan3::data_layout::uncompressed> &,
                                                            build_arguments const &);

template size_t initialise_max_bin_kmers<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                                                          std::vector<int64_t> &,
                                                                          std::vector<int64_t> &,
                                                                          lemon::ListDigraph::Node const &,
                                                                          build_data<seqan3::data_layout::compressed> &,
                                                                          build_arguments const &);
*/
} 
