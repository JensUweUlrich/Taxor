#include <lemon/list_graph.h> /// Must be first include.

#include "chopper_build.hpp"
#include "create_ixfs_from_chopper_pack.hpp"
#include "store_index.hpp"
#include "index.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void chopper_build(build_arguments const & arguments)
{
    build_data data{};

    create_ixfs_from_chopper_pack(data, arguments);
    std::cout << data.number_of_ixfs << std::endl;
	std::cout << data.number_of_user_bins << std::endl;
    
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hixf.user_bins.num_user_bins(); ++i)
        bin_path.push_back(std::vector<std::string>{data.hixf.user_bins.filename_of_user_bin(i)});

    raptor_index<hierarchical_interleaved_xor_filter<uint8_t>> index{window{arguments.window_size},
                                                                                arguments.shape,
                                                                                arguments.kmer_size,
                                                                                arguments.syncmer_size,
                                                                                arguments.t_syncmer,
                                                                                arguments.parts,
                                                                                arguments.compressed,
                                                                                bin_path,
                                                                                std::move(data.hixf)};
    
    store_index(arguments.out_path, index, arguments);
}

//template void chopper_build<seqan3::data_layout::uncompressed>(build_arguments const &);

//template void chopper_build<seqan3::data_layout::compressed>(build_arguments const &);

} 
