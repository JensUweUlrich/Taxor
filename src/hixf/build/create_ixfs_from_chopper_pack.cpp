
#include <lemon/list_graph.h> /// Must be first include.

#include "create_ixfs_from_chopper_pack.hpp"
#include "hierarchical_build.hpp"
#include "read_chopper_pack_file.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void create_ixfs_from_chopper_pack(build_data& data, build_arguments const & arguments)
{
    read_chopper_pack_file(data, arguments.bin_file);
    lemon::ListDigraph::Node root = data.ixf_graph.nodeFromId(0); // root node = high level IBF node
    robin_hood::unordered_flat_set<size_t> root_hashes{};

    size_t const t_max{data.node_map[root].number_of_technical_bins};

    hierarchical_build(root_hashes, root, data, arguments, true);
}

/*
template void
create_ibfs_from_chopper_pack<seqan3::data_layout::uncompressed>(build_data<seqan3::data_layout::uncompressed> &,
                                                                 build_arguments const &);
template void
create_ibfs_from_chopper_pack<seqan3::data_layout::compressed>(build_data<seqan3::data_layout::compressed> &,
                                                               build_arguments const & arguments);
*/
} 
