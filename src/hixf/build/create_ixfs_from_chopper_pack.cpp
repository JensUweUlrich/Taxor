
#include <lemon/list_graph.h> /// Must be first include.

#include "create_ixfs_from_chopper_pack.hpp"
#include "hierarchical_build.hpp"
#include "read_chopper_pack_file.hpp"
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void create_ixfs_from_chopper_pack(build_data& data, build_arguments const & arguments)
{
    read_chopper_pack_file(data, arguments.bin_file);
    lemon::ListDigraph::Node root = data.ixf_graph.nodeFromId(0); // root node = high level IXF node
    robin_hood::unordered_flat_set<size_t> root_hashes{};

    size_t const t_max{data.node_map[root].number_of_technical_bins};

    hierarchical_build(root_hashes, root, data, arguments, true, false);
    typedef seqan3::interleaved_xor_filter<>::counting_agent_type< uint64_t > TIXFAgent;
    TIXFAgent ixf_count_agent = data.hixf.ixf_vector[0].counting_agent< uint64_t >();
    robin_hood::unordered_flat_set<size_t> read_hashes{};
    using sequence_file_t = seqan3::sequence_file_input<hixf::dna4_traits, seqan3::fields<seqan3::field::seq>>;
	for (auto && [seq] : sequence_file_t{"/media/jens/INTENSO/refseq-viral/2022-03-23_22-07-02/files.renamed/GCF_000839085.1_genomic.fna.gz"})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{hixf::adjust_seed(arguments.shape.count())}))
                    read_hashes.insert(hash);
	std::vector<size_t> c{};
    std::ranges::copy(read_hashes, std::back_inserter(c));

	auto result = ixf_count_agent.bulk_count(c);
    seqan3::debug_stream << "Root result: " << result << "\n";
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