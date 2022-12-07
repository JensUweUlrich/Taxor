
#pragma once

//#include <sharg/exceptions.hpp>

#include <build/build_arguments.hpp>
#include <build/hierarchical_interleaved_xor_filter.hpp>
#include <build/strong_types.hpp>

#include <Species.hpp>

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

template <typename data_t = ixf_t>
class taxor_index
{
private:
    template <typename friend_data_t>
    friend class taxor_index;

    uint64_t window_size_{};
    seqan3::shape shape_{};
    uint8_t parts_{};
    bool compressed_{};
    bool use_syncmer_{};
    std::vector<std::vector<std::string>> bin_path_{};
    std::vector<taxonomy::Species> species_{};
    data_t ixf_{};
    uint8_t syncmer_size_{};
    uint8_t t_syncmer_{};
    uint8_t kmer_size_{};


public:
    //static constexpr seqan3::data_layout data_layout_mode = data_t::data_layout_mode;
    static constexpr uint32_t version{1u};

    taxor_index() = default;
    taxor_index(taxor_index const &) = default;
    taxor_index(taxor_index &&) = default;
    taxor_index & operator=(taxor_index const &) = default;
    taxor_index & operator=(taxor_index &&) = default;
    ~taxor_index() = default;

    explicit taxor_index(hixf::window const window_size,
                         seqan3::shape const shape,
                         u_int8_t kmer_size,
                         uint8_t syncmer_size,
                         uint8_t t_syncmer,
                         uint8_t const parts,
                         bool const use_syncmer,
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
        compressed_{compressed},
        bin_path_{bin_path},
        species_{species},
        ixf_{std::move(ixf)}
    {}

    explicit taxor_index(hixf::build_arguments const & arguments) :
        window_size_{arguments.window_size},
        shape_{arguments.shape},
        kmer_size_{arguments.kmer_size},
        syncmer_size_{arguments.syncmer_size},
        t_syncmer_{arguments.t_syncmer},
        parts_{arguments.parts},
        use_syncmer_{arguments.compute_syncmer},
        compressed_{arguments.compressed},
        bin_path_{arguments.bin_path},
        species_{},
        ixf_{}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
    }

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
        compressed_ = true;
        bin_path_ = other.bin_path_;
        species_ = other.species_;
        ixf_ = data_t{other.ibf_};
    }

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
        compressed_ = true;
        bin_path_ = std::move(other.bin_path_);
        species_ = std::move(other.species_);
        ixf_ = std::move(data_t{std::move(other.ixf_)});
    }

    uint64_t window_size() const
    {
        return window_size_;
    }

    seqan3::shape shape() const
    {
        return shape_;
    }

    uint8_t kmer_size() const
    {
        return kmer_size_;
    }

    uint8_t syncmer_size() const
    {
        return syncmer_size_;
    }

    uint8_t t_syncmer() const
    {
        return t_syncmer_;
    }

    uint8_t parts() const
    {
        return parts_;
    }

    bool use_syncmer() const
    {
        return use_syncmer_;
    }

    bool compressed() const
    {
        return compressed_;
    }

    std::vector<std::vector<std::string>> const & bin_path() const
    {
        return bin_path_;
    }

    std::vector<taxonomy::Species> const & species() const
    {
        return species_;
    }

    data_t & ixf()
    {
        return ixf_;
    }

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
                //throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
            }
        }
        else
        {
            //throw sharg::parser_error{"Unsupported index version. Check taxor upgrade."}; // GCOVR_EXCL_LINE
        }
    }

    /* \brief Serialisation support function. Do not load the actual data.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_input_archive.
     * \param[in] archive The archive being serialised from/to.
     * \param[in] version Index version.
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
                archive(compressed_);
                archive(bin_path_);
                archive(species_);
            }
            // GCOVR_EXCL_START
            catch (std::exception const & e)
            {
                //throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
            }
            // GCOVR_EXCL_STOP
        }
        else
        {
            //throw sharg::parser_error{"Unsupported index version. Check taxor upgrade."}; // GCOVR_EXCL_LINE
        }
    }
    //!\endcond
};

}
