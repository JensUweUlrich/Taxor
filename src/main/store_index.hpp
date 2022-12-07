
#pragma once

#include <filesystem>

#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include "index.hpp"
#include <build/strong_types.hpp>

namespace taxor
{

template <typename data_t, typename arguments_t>
static inline void
store_index(std::filesystem::path const & path, taxor_index<data_t> const & index, arguments_t const &)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

//template <seqan3::data_layout layout, typename arguments_t>
template <typename arguments_t>
static inline void store_index(std::filesystem::path const & path,
                               seqan3::interleaved_xor_filter<> && ixf,
                               arguments_t const & arguments)
{
    taxor_index<seqan3::interleaved_xor_filter<>> index{hixf::window{arguments.window_size},
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
