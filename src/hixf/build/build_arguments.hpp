
#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

#include "strong_types.hpp"

namespace hixf
{

/*!\file build_arguments.hpp
 * \brief Provides hixf::build_arguments.
 */

/*!\brief Plain aggregate of all user-configurable settings for building a Hierarchical Interleaved XOR Filter.
 * \details
 *
 * Instances of this struct are populated from command line arguments (see the `taxor build` subcommand) and
 * are then threaded through essentially every function in this module (hixf::create_ixfs_from_chopper_pack,
 * hixf::hierarchical_build, hixf::compute_hashes, ...) to control k-mer/syncmer computation and IXF sizing.
 */
struct build_arguments
{
    // Related to k-mers
    uint8_t kmer_size{20u};             //!< The k-mer size used for (s)k-mer/minimiser hashing.
    uint32_t window_size{kmer_size};    //!< The minimiser window size (only used when `compute_minimiser` is set).
    window window_size_strong{kmer_size}; //!< Strongly-typed variant of #window_size.
    std::string shape_string{};         //!< User-provided textual representation of the (possibly gapped) shape.
    seqan3::shape shape{seqan3::ungapped{kmer_size}}; //!< The shape (contiguous or gapped) used for hashing.
    bool compute_syncmer{true};         //!< If true, use closed syncmers instead of plain k-mers/minimisers.
    bool compute_minimiser{false};      //!< If true, use minimisers instead of syncmers.
    bool disable_cutoffs{false};        //!< If true, disable any minimum-size/coverage cutoffs applied upstream.
    uint8_t syncmer_size{10u};          //!< The s-mer size used when computing syncmers.
    uint8_t t_syncmer{6u};              //!< Syncmer parameter `t`, passed through to hashing::seq_to_syncmers.
    uint16_t scaling{1u};               //!< Downsampling factor; only ~1/`scaling` of hashes are kept (FracMinHash-style).

    // Related to IXF
    std::filesystem::path out_path{"./"}; //!< Output directory/path for the resulting index.
    std::string size{"1k"};             //!< Human-readable target size string (currently informational).
    uint64_t bins{64};                  //!< Default number of technical bins.
    uint64_t bits{4096};                //!< Default number of bits per technical bin.
    uint64_t hash{2};                   //!< Number of hash functions used by the underlying filter.
    uint8_t parts{1u};                  //!< Number of parts the IXF is split into (currently unused by this module).
    double fpr{0.05};                   //!< Target false-positive rate used to size technical bins.
    bool compressed{false};             //!< Whether the resulting filter should be compressed.

    // General arguments
    std::vector<std::vector<std::string>> bin_path{}; //!< Per-bin lists of input reference filenames.
    std::filesystem::path bin_file{};   //!< Path to the chopper pack layout file describing the bin hierarchy.
    uint8_t threads{1u};                //!< Number of threads to use while building.
    bool is_socks{false};               //!< Set when building in "socks" (legacy/alternative) mode.
    bool is_hixf{false};                //!< Set when building a Hierarchical Interleaved XOR Filter.
    bool is_minimiser{false};           //!< Set when the input already consists of precomputed minimisers.
    bool is_syncmer{false};             //!< Set when the input already consists of precomputed syncmers.
};

}
