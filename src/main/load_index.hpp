
#pragma once

#include <chrono>

#include <search/search_arguments.hpp>
#include "index.hpp"

/*!\file load_index.hpp
 * \brief Thin wrappers around cereal binary-archive deserialisation of a taxor::taxor_index from disk.
 */

namespace taxor
{

/*!\brief Load one part of a partitioned index, given search arguments and the part number.
 *
 * Appends "_<part>" to the base index file path taken from \p arguments.index_file and delegates
 * to the std::filesystem::path overload below.
 *
 * \note No call sites for this overload were found elsewhere in the codebase as of this writing;
 * taxor_search.cpp only uses the two overloads below.
 * \tparam index_t Type of the taxor_index (or compatible type) to load into.
 * \param[out] index The index object to deserialise into.
 * \param arguments Search arguments; only arguments.index_file is used here.
 * \param part Index of the partition to load; appended as a "_<part>" suffix to the index file name.
 * \param[in,out] index_io_time Accumulator (in seconds) that the time spent deserialising is added to.
 */
template <typename index_t>
void load_index(index_t & index, hixf::search_arguments const & arguments, size_t const part, double & index_io_time)
{
    std::filesystem::path index_file{arguments.index_file};
    index_file += "_" + std::to_string(part);

    load_index(index, index_file, index_io_time);
}

/*!\brief Load a (non-partitioned) index given search arguments.
 * \tparam index_t Type of the taxor_index (or compatible type) to load into.
 * \param[out] index The index object to deserialise into.
 * \param arguments Search arguments; only arguments.index_file is used here.
 * \param[in,out] index_io_time Accumulator (in seconds) that the time spent deserialising is added to.
 */
template <typename index_t>
void load_index(index_t & index, hixf::search_arguments const & arguments, double & index_io_time)
{
    load_index(index, arguments.index_file, index_io_time);
}

/*!\brief Deserialise a taxor_index from a binary cereal archive at the given path.
 *
 * Opens \p path as a binary ifstream, wraps it in a cereal::BinaryInputArchive and reads \p index
 * from it (invoking taxor_index::CEREAL_SERIALIZE_FUNCTION_NAME() internally). The wall-clock time
 * taken by the deserialisation call is measured and added to \p index_io_time.
 * \tparam index_t Type of the taxor_index (or compatible type) to load into.
 * \param[out] index The index object to deserialise into.
 * \param path Path of the index file to read.
 * \param[in,out] index_io_time Accumulator (in seconds) that the time spent deserialising is added to.
 */
template <typename index_t>
void load_index(index_t & index, std::filesystem::path const & path, double & index_io_time)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(index);
    auto end = std::chrono::high_resolution_clock::now();

    index_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

}
