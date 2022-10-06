// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/list_graph.h> /// Must be first include.

#include "chopper_build.hpp"
#include "create_ibfs_from_chopper_pack.hpp"
#include "store_index.hpp"
#include "index.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void chopper_build(build_arguments const & arguments)
{
    build_data data{};

    create_ibfs_from_chopper_pack(data, arguments);

    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hixf.user_bins.num_user_bins(); ++i)
        bin_path.push_back(std::vector<std::string>{data.hixf.user_bins.filename_of_user_bin(i)});

    raptor_index<hierarchical_interleaved_bloom_filter<data_layout_mode>> index{window{arguments.window_size},
                                                                                arguments.shape,
                                                                                arguments.parts,
                                                                                arguments.compressed,
                                                                                bin_path,
                                                                                std::move(data.hibf)};

    store_index(arguments.out_path, index, arguments);
}

//template void chopper_build<seqan3::data_layout::uncompressed>(build_arguments const &);

//template void chopper_build<seqan3::data_layout::compressed>(build_arguments const &);

} // namespace raptor::hibf
