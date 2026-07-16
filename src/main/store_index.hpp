// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>

#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include "index.hpp"
#include <build/strong_types.hpp>

/*!\file store_index.hpp
 * \brief Thin wrappers around cereal binary-archive serialisation of a taxor::taxor_index to disk.
 */

namespace taxor
{

/*!\brief Serialise a taxor_index to a binary cereal archive at the given path.
 *
 * Opens \p path as a binary ofstream, wraps it in a cereal::BinaryOutputArchive and writes
 * \p index into it (invoking taxor_index::CEREAL_SERIALIZE_FUNCTION_NAME() internally). This is
 * the overload used by taxor_build.cpp to persist a finished index.
 * \tparam data_t Filter data type of \p index (see taxor_index).
 * \param path Destination path the index is written to.
 * \param index The taxor_index to serialise.
 */
template <typename data_t>
static inline void
store_index(std::filesystem::path const & path, taxor_index<data_t> const & index)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

//template <seqan3::data_layout layout, typename arguments_t>
/*!\brief Build a taxor_index<seqan3::interleaved_xor_filter<>> from a filter and build arguments, then serialise it.
 *
 * \note No call sites for this overload were found elsewhere in the codebase as of this writing.
 * It also constructs `taxor_index<seqan3::interleaved_xor_filter<>>` with only 6 arguments
 * (window, shape, parts, compressed, bin_path, ixf), whereas the corresponding taxor_index
 * constructor in index.hpp currently expects 12 (it additionally takes kmer_size, syncmer_size,
 * t_syncmer, use_syncmer, scaling and species); as a template, this is only checked by the
 * compiler if/when it is instantiated.
 * \tparam arguments_t Type of \p arguments (expected to expose window_size, shape, parts,
 *                     compressed and bin_path members).
 * \param path Destination path the index is written to.
 * \param ixf The filter data to move into the newly constructed taxor_index.
 * \param arguments Build arguments the taxor_index is constructed from.
 */
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
