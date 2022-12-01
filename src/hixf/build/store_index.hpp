
#pragma once

#include <filesystem>

//#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>
#include <seqan3/search/dream_index/interleaved_binary_fuse_filter.hpp>

#include "index.hpp"
#include "strong_types.hpp"

namespace hixf
{

template <typename data_t, typename arguments_t>
static inline void
store_index(std::filesystem::path const & path, raptor_index<data_t> const & index, arguments_t const &)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

//template <seqan3::data_layout layout, typename arguments_t>
template <typename arguments_t>
static inline void store_index(std::filesystem::path const & path,
//                               seqan3::interleaved_xor_filter<> && ixf,
                               seqan3::interleaved_binary_fuse_filter<> && ixf,
                               arguments_t const & arguments)
{
//    raptor_index<seqan3::interleaved_xor_filter<>> index{window{arguments.window_size},
    raptor_index<seqan3::interleaved_binary_fuse_filter<>> index{window{arguments.window_size},
                                                                 arguments.shape,
                                                                 arguments.parts,
                                                                 arguments.compressed,
                                                                 arguments.bin_path,
                                                                 std::move(ixf)};

    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

}
