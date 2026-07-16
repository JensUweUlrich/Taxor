
#pragma once

//#include <sharg/exceptions.hpp>

#include <stdexcept>

#include <build/build_arguments.hpp>
#include <build/hierarchical_interleaved_xor_filter.hpp>
#include <build/strong_types.hpp>

#include <Species.hpp>

/*!\file index.hpp
 * \brief Defines taxor::taxor_index, the in-memory representation of a (built or loaded)
 *        Hierarchical Interleaved XOR Filter (HIXF) based taxonomic index, plus its cereal-based
 *        binary (de)serialisation support.
 */

namespace taxor
{

using ixf_t = seqan3::interleaved_xor_filter<uint8_t>;
//using ibf_compressed = seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>;
using hixf_t = hixf::hierarchical_interleaved_xor_filter<uint8_t>;
//using hibf_compressed = hierarchical_interleaved_bloom_filter<seqan3::data_layout::compressed>;

/*template <typename return_t, typename input_t>
concept compressible_from = (std::same_as<return_t, ibf_compressed> && std::same_as<input_t, ibf>)
                         || (std::same_as<return_t, hibf_compressed> && std::same_as<input_t, hibf>);
*/

/*!\brief In-memory representation of a taxor index.
 *
 * A taxor_index bundles the parameters an index was built with (window size, k-mer/syncmer
 * shape and sizes, scaling factor, number of parts, whether the filter is stored compressed),
 * the per-user-bin input file paths (bin_path_) and taxonomic metadata (species_), together with
 * the underlying filter data structure (ixf_, typically a
 * hixf::hierarchical_interleaved_xor_filter) that membership queries are actually run against.
 *
 * Instances are (de)serialised to/from disk in a single cereal binary archive via
 * CEREAL_SERIALIZE_FUNCTION_NAME(); load_parameters() offers a cheaper path that reads back only
 * the build parameters and species metadata without touching the (potentially large) filter data.
 *
 * \tparam data_t Type of the underlying filter data structure, e.g. seqan3::interleaved_xor_filter<uint8_t>
 *                (single-level index) or hixf_t == hixf::hierarchical_interleaved_xor_filter<uint8_t>
 *                (hierarchical index). Defaults to ixf_t.
 */
template <typename data_t = ixf_t>
class taxor_index
{
private:
    // Grants taxor_index<other_data_t> access to this instance's private members, needed by the
    // converting constructors below that build a taxor_index<data_t> out of a taxor_index<other_data_t>.
    template <typename friend_data_t>
    friend class taxor_index;

    uint64_t window_size_{}; //!< Minimizer/syncmer window size used when building the index.
    seqan3::shape shape_{}; //!< Shape (gapped/ungapped k-mer mask) used when building the index.
    uint8_t parts_{}; //!< Number of parts the index is split into.
    bool compressed_{}; //!< Whether the underlying filter is stored in compressed layout.
    bool use_syncmer_{}; //!< Whether syncmers (true) or plain minimizers (false) were used.
    std::vector<std::vector<std::string>> bin_path_{}; //!< Per user bin: list of input file paths that were indexed into it.
    std::vector<taxonomy::Species> species_{}; //!< Per user bin: taxonomic metadata, aligned with bin_path_.
    data_t ixf_{}; //!< The underlying (hierarchical) interleaved XOR filter holding the index data.
    uint8_t syncmer_size_{}; //!< Size of the syncmers used when building the index (if use_syncmer_ is set).
    uint8_t t_syncmer_{}; //!< The t parameter of the syncmer scheme (position of the required minimal s-mer within the k-mer window).
    uint8_t kmer_size_{}; //!< Size of the k-mers used when building the index.
    uint16_t scaling_{}; //!< Scaling factor applied during index construction.


public:
    //static constexpr seqan3::data_layout data_layout_mode = data_t::data_layout_mode;
    static constexpr uint32_t version{1u};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor.
    taxor_index() = default;
    //!\brief Copy constructor.
    taxor_index(taxor_index const &) = default;
    //!\brief Move constructor.
    taxor_index(taxor_index &&) = default;
    //!\brief Copy assignment operator.
    taxor_index & operator=(taxor_index const &) = default;
    //!\brief Move assignment operator.
    taxor_index & operator=(taxor_index &&) = default;
    //!\brief Destructor.
    ~taxor_index() = default;
    //!\}

    /*!\brief Construct a taxor_index from explicit build parameters and an already-built filter.
     * \param window_size Minimizer/syncmer window size used when building the filter.
     * \param shape Shape (gapped/ungapped k-mer mask) used when building the filter.
     * \param kmer_size Size of the k-mers used when building the filter.
     * \param syncmer_size Size of the syncmers used when building the filter (if use_syncmer is set).
     * \param t_syncmer The t parameter of the syncmer scheme.
     * \param parts Number of parts the index is split into.
     * \param use_syncmer Whether syncmers (true) or plain minimizers (false) were used to build the filter.
     * \param scaling Scaling factor applied during index construction.
     * \param compressed Whether the underlying filter is stored in compressed layout.
     * \param bin_path Per user bin: list of input file paths that were indexed into it.
     * \param species Per user bin: taxonomic metadata, aligned with bin_path.
     * \param ixf The already-built filter data structure; ownership is taken via move.
     */
    explicit taxor_index(hixf::window const window_size,
                         seqan3::shape const shape,
                         uint8_t kmer_size,
                         uint8_t syncmer_size,
                         uint8_t t_syncmer,
                         uint8_t const parts,
                         bool const use_syncmer,
                         uint16_t scaling,
                         bool const compressed,
                         std::vector<std::vector<std::string>> const & bin_path,
                         std::vector<taxonomy::Species> &species,
                         data_t && ixf) :
        window_size_{window_size.v},
        shape_{shape},
        kmer_size_{kmer_size},
        syncmer_size_{syncmer_size},
        t_syncmer_{t_syncmer},
        parts_{parts},
        use_syncmer_{use_syncmer},
        scaling_{scaling},
        compressed_{compressed},
        bin_path_{bin_path},
        species_{species},
        ixf_{std::move(ixf)}
    {}

    /*!\brief Construct an (as yet empty) taxor_index by copying the build parameters out of a hixf::build_arguments object.
     *
     * The underlying filter (ixf_) and species metadata (species_) are left default-constructed;
     * they are expected to be filled in later, in the course of building the index.
     * \param arguments Parsed build arguments describing how the index is to be built.
     */
    explicit taxor_index(hixf::build_arguments const & arguments) :
        window_size_{arguments.window_size},
        shape_{arguments.shape},
        kmer_size_{arguments.kmer_size},
        syncmer_size_{arguments.syncmer_size},
        t_syncmer_{arguments.t_syncmer},
        parts_{arguments.parts},
        use_syncmer_{arguments.compute_syncmer},
        scaling_{arguments.scaling},
        compressed_{arguments.compressed},
        bin_path_{arguments.bin_path},
        species_{},
        ixf_{}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
    }

    /*!\brief Converting copy constructor: builds a taxor_index around a taxor_index<other_data_t>'s data.
     *
     * Copies all build parameters and species metadata from \p other and reconstructs the filter
     * as a data_t from other's filter (data_t{other.ixf_}).
     * \tparam other_data_t Filter data type of \p other.
     * \param other The taxor_index to convert from.
     */
    template <typename other_data_t>
    explicit taxor_index(taxor_index<other_data_t> const & other)
    {
        //static_assert(index_structure::compressible_from<data_t, other_data_t>);
        window_size_ = other.window_size_;
        shape_ = other.shape_;
        kmer_size_ = other.kmer_size_;
        syncmer_size_ = other.syncmer_size_;
        t_syncmer_ = other.t_syncmer_;
        parts_ = other.parts_;
        use_syncmer_ = other.use_syncmer_;
        scaling_ = other.scaling_;
        // Converting to a different data_t is only meaningful for producing a compressed index,
        // so compressed_ is unconditionally set to true here regardless of other.compressed_.
        compressed_ = true;
        bin_path_ = other.bin_path_;
        species_ = other.species_;
        ixf_ = data_t{other.ixf_};
    }

    /*!\brief Converting move constructor: builds a taxor_index by moving out of a taxor_index<other_data_t>.
     *
     * Moves all build parameters and species metadata out of \p other and reconstructs the filter
     * as a data_t from other's (moved-from) filter.
     * \tparam other_data_t Filter data type of \p other.
     * \param other The taxor_index to convert-move from; left in a moved-from state afterwards.
     */
    template <typename other_data_t>
    explicit taxor_index(taxor_index<other_data_t> && other)
    {
        //static_assert(index_structure::compressible_from<data_t, other_data_t>);
        window_size_ = std::move(other.window_size_);
        shape_ = std::move(other.shape_);
        kmer_size_ = std::move(other.kmer_size_);
        syncmer_size_ = std::move(other.syncmer_size_);
        t_syncmer_ = std::move(other.t_syncmer_);
        parts_ = std::move(other.parts_);
        use_syncmer_ = std::move(other.use_syncmer_);
        scaling_ = std::move(other.scaling_);
        // See the copy-converting constructor above: compressed_ is always true after a conversion.
        compressed_ = true;
        bin_path_ = std::move(other.bin_path_);
        species_ = std::move(other.species_);
        ixf_ = std::move(data_t{std::move(other.ixf_)});
    }

    //!\brief Returns the window size used to build the index.
    uint64_t window_size() const
    {
        return window_size_;
    }

    //!\brief Returns the shape used to build the index.
    seqan3::shape shape() const
    {
        return shape_;
    }

    //!\brief Returns the k-mer size used to build the index.
    uint8_t kmer_size() const
    {
        return kmer_size_;
    }

    //!\brief Returns the syncmer size used to build the index (meaningful only if use_syncmer() is true).
    uint8_t syncmer_size() const
    {
        return syncmer_size_;
    }

    //!\brief Returns the t parameter of the syncmer scheme used to build the index.
    uint8_t t_syncmer() const
    {
        return t_syncmer_;
    }

    //!\brief Returns the number of parts the index is split into.
    uint8_t parts() const
    {
        return parts_;
    }

    //!\brief Returns whether syncmers (true) or plain minimizers (false) were used to build the index.
    bool use_syncmer() const
    {
        return use_syncmer_;
    }

    //!\brief Returns the scaling factor applied during index construction.
    uint16_t scaling() const
    {
        return scaling_;
    }

    //!\brief Returns whether the underlying filter is stored in compressed layout.
    bool compressed() const
    {
        return compressed_;
    }

    //!\brief Returns, per user bin, the list of input file paths that were indexed into it.
    std::vector<std::vector<std::string>> const & bin_path() const
    {
        return bin_path_;
    }

    //!\brief Returns the per-user-bin taxonomic metadata, aligned with bin_path().
    std::vector<taxonomy::Species> const & species() const
    {
        return species_;
    }

    //!\brief Returns a mutable reference to the underlying filter data structure.
    data_t & ixf()
    {
        return ixf_;
    }

    //!\brief Returns a const reference to the underlying filter data structure.
    data_t const & ixf() const
    {
        return ixf_;
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        uint32_t parsed_version{taxor_index<>::version};
        archive(parsed_version);
        if (parsed_version == taxor_index<>::version)
        {
            try
            {
                archive(window_size_);
                archive(shape_);
                archive(kmer_size_);
                archive(syncmer_size_);
                archive(t_syncmer_);
                archive(parts_);
                archive(use_syncmer_);
                archive(scaling_);
                archive(compressed_);
                /*if ((data_layout_mode == seqan3::data_layout::compressed && !compressed_)
                    || (data_layout_mode == seqan3::data_layout::uncompressed && compressed_))
                {
                    throw sharg::parser_error{"Data layouts of serialised and specified index differ."};
                }*/
                archive(bin_path_);
                archive(species_);
                archive(ixf_);
            }
            catch (std::exception const & e)
            {
                throw std::runtime_error{"Cannot read index: " + std::string{e.what()}};
            }
        }
        else
        {
            throw std::runtime_error{"Unsupported index version. Check taxor upgrade."}; // GCOVR_EXCL_LINE
        }
    }

    /*!\brief Serialisation support function. Do not load the actual data.
     *
     * Reads back only the build parameters and species metadata (not the filter itself, ixf_),
     * which is considerably cheaper than a full load and is used where just the index's
     * k-mer/syncmer scheme and taxonomy need to be inspected without loading the whole filter.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_input_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        uint32_t parsed_version{};
        archive(parsed_version);
        if (parsed_version == version)
        {
            try
            {
                archive(window_size_);
                archive(shape_);
                archive(kmer_size_);
                archive(syncmer_size_);
                archive(t_syncmer_);
                archive(parts_);
                archive(use_syncmer_);
                archive(scaling_);
                archive(compressed_);
                archive(bin_path_);
                archive(species_);
            }
            // GCOVR_EXCL_START
            catch (std::exception const & e)
            {
                throw std::runtime_error{"Cannot read index: " + std::string{e.what()}};
            }
            // GCOVR_EXCL_STOP
        }
        else
        {
            throw std::runtime_error{"Unsupported index version. Check taxor upgrade."}; // GCOVR_EXCL_LINE
        }
    }
    //!\endcond
};

}
