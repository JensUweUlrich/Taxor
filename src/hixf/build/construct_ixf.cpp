
#include <lemon/list_graph.h> /// Must be first include.

#include "bin_size_in_bits.hpp"
#include "construct_ixf.hpp"
#include "insert_into_ixf.hpp"
#include "temp_hash_file.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
// @deprecated
seqan3::interleaved_xor_filter<> construct_ixf(robin_hood::unordered_flat_set<size_t> & parent_hashes,
                                                 robin_hood::unordered_flat_set<size_t> & hashes,
                                                 size_t const number_of_bins,
                                                 lemon::ListDigraph::Node const & node,
                                                 build_data & data,
                                                 build_arguments const & arguments,
                                                 bool is_root)
{
    auto & node_data = data.node_map[node];

    size_t const hashes_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(hashes.size()) / number_of_bins))};
    //double const bin_bits{static_cast<double>(bin_size_in_bits(arguments, kmers_per_bin))};
    //seqan3::bin_size const bin_size{static_cast<size_t>(std::ceil(bin_bits * data.fp_correction[number_of_bins]))};
    //seqan3::bin_count const bin_count{node_data.number_of_technical_bins};
    seqan3::interleaved_xor_filter<> ixf{node_data.number_of_technical_bins, hashes_per_bin};

    insert_into_ixf(parent_hashes, hashes, number_of_bins, node_data.max_bin_index, ixf, is_root);

    return ixf;
}

seqan3::interleaved_xor_filter<> construct_ixf(std::vector<robin_hood::unordered_flat_set<size_t>> &node_hashes)
{
    std::vector<std::vector<size_t>> tmp{};
    for (auto hash_bin : node_hashes)
    {
        std::vector<size_t> c{};
        std::ranges::copy(hash_bin, std::back_inserter(c));
        tmp.emplace_back(c);
    }
    seqan3::interleaved_xor_filter<> ixf{tmp};

    return ixf;
}


seqan3::interleaved_xor_filter<> construct_ixf(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               std::vector<robin_hood::unordered_flat_set<size_t>> &node_hashes,
                                               size_t const & current_node_ixf_pos)
{
    auto &current_node_data = data.node_map[current_node];
    // create empty IXF based on number of technical bins and max number of hashes per bin
    seqan3::interleaved_xor_filter<> ixf{current_node_data.number_of_technical_bins, current_node_data.max_bin_hashes};
    // first iterate over all child IXFs 
    
    bool success{false};

    while (!success)
    {
        // open file to store all hashes of current node IXF
        std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(current_node_ixf_pos) + ".tmp";
        auto tmp_file = std::filesystem::temp_directory_path() / ixf_tmp_name;
        std::ofstream tmp_stream{tmp_file};
        //std::cout << "create IXF number " << current_node_ixf_pos << std::endl;
        success = true;
        for (lemon::ListDigraph::OutArcIt arc_it(data.ixf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
        {
            auto child = data.ixf_graph.target(arc_it);
            auto& child_node_data = data.node_map[child];
            size_t child_ixf_pos = ixf_positions[child_node_data.parent_bin_index];
            std::vector<size_t> hashes{};
            // read in hashes of child IXF and add all hashes to corresponding bin of current node IXF
            read_from_temp_hash_file(ixf_positions[data.node_map[child].parent_bin_index], hashes);
            success = ixf.add_bin_elements(data.node_map[child].parent_bin_index, hashes);
            if(!success)
                break;
            
            for (size_t h : hashes)
                tmp_stream << h << " ";

        }
        // reset seed if adding bin to IXF was not successful
        if (!success)
        {   
            ixf.clear();
            ixf.set_seed();
            tmp_stream.close();
            continue;
        }

        // iterate over new hashes
        // add hashes of bins for newly computed hashes on that level
        uint16_t bin_idx{0};
        for (auto hash_bin : node_hashes)
        {
            if (hash_bin.size() == 0)
            {
                bin_idx++;
                continue;
            }
            std::vector<size_t> c{};
            std::ranges::copy(hash_bin, std::back_inserter(c));
            success = ixf.add_bin_elements(bin_idx, c);
            if(!success)
                break;

            bin_idx++;

            for (size_t h : c)
                tmp_stream << h << " ";
        }

        tmp_stream.close();
        // reset seed if adding bin to IXF was not successful
        if (!success)
        {   
            ixf.clear();
            ixf.set_seed();
            continue;
        }
        //tmp_files.emplace_back(tmp_file);
    }
    

    return ixf;
}

/*template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::uncompressed>(robin_hood::unordered_flat_set<size_t> &,
                                                 robin_hood::unordered_flat_set<size_t> &,
                                                 size_t const,
                                                 lemon::ListDigraph::Node const &,
                                                 build_data<seqan3::data_layout::uncompressed> &,
                                                 build_arguments const &,
                                                 bool);

template seqan3::interleaved_bloom_filter<>
construct_ibf<seqan3::data_layout::compressed>(robin_hood::unordered_flat_set<size_t> &,
                                               robin_hood::unordered_flat_set<size_t> &,
                                               size_t const,
                                               lemon::ListDigraph::Node const &,
                                               build_data<seqan3::data_layout::compressed> &,
                                               build_arguments const &,
                                               bool);
*/
} 
