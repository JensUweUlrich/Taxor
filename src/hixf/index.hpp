
#pragma once

//#include <sharg/exceptions.hpp>

#include "build_arguments.hpp"
#include "hierarchical_interleaved_xor_filter.hpp"
#include "strong_types.hpp"

namespace hixf
{

namespace index_structure
{

using ixf = seqan3::interleaved_xor_filter<uint8_t>;
//using ibf_compressed = seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>;
using hixf = hierarchical_interleaved_xor_filter<uint8_t>;
//using hibf_compressed = hierarchical_interleaved_bloom_filter<seqan3::data_layout::compressed>;

/*template <typename return_t, typename input_t>
concept compressible_from = (std::same_as<return_t, ibf_compressed> && std::same_as<input_t, ibf>)
                         || (std::same_as<return_t, hibf_compressed> && std::same_as<input_t, hibf>);
*/
} // namespace index_structure

template <typename data_t = index_structure::ixf>
class raptor_index
{
private:
    template <typename friend_data_t>
    friend class raptor_index;

    uint64_t window_size_{};
    seqan3::shape shape_{};
    uint8_t parts_{};
    bool compressed_{};
    std::vector<std::vector<std::string>> bin_path_{};
    data_t ixf_{};

public:
    //static constexpr seqan3::data_layout data_layout_mode = data_t::data_layout_mode;
    static constexpr uint32_t version{1u};

    raptor_index() = default;
    raptor_index(raptor_index const &) = default;
    raptor_index(raptor_index &&) = default;
    raptor_index & operator=(raptor_index const &) = default;
    raptor_index & operator=(raptor_index &&) = default;
    ~raptor_index() = default;

    explicit raptor_index(window const window_size,
                          seqan3::shape const shape,
                          uint8_t const parts,
                          bool const compressed,
                          std::vector<std::vector<std::string>> const & bin_path,
                          data_t && ixf) :
        window_size_{window_size.v},
        shape_{shape},
        parts_{parts},
        compressed_{compressed},
        bin_path_{bin_path},
        ixf_{std::move(ixf)}
    {}

    explicit raptor_index(build_arguments const & arguments) :
        window_size_{arguments.window_size},
        shape_{arguments.shape},
        parts_{arguments.parts},
        compressed_{arguments.compressed},
        bin_path_{arguments.bin_path},
        ixf_{}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
    }

    template <typename other_data_t>
    explicit raptor_index(raptor_index<other_data_t> const & other)
    {
        //static_assert(index_structure::compressible_from<data_t, other_data_t>);
        window_size_ = other.window_size_;
        shape_ = other.shape_;
        parts_ = other.parts_;
        compressed_ = true;
        bin_path_ = other.bin_path_;
        ixf_ = data_t{other.ibf_};
    }

    template <typename other_data_t>
    explicit raptor_index(raptor_index<other_data_t> && other)
    {
        //static_assert(index_structure::compressible_from<data_t, other_data_t>);
        window_size_ = std::move(other.window_size_);
        shape_ = std::move(other.shape_);
        parts_ = std::move(other.parts_);
        compressed_ = true;
        bin_path_ = std::move(other.bin_path_);
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

    uint8_t parts() const
    {
        return parts_;
    }

    bool compressed() const
    {
        return compressed_;
    }

    std::vector<std::vector<std::string>> const & bin_path() const
    {
        return bin_path_;
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
        uint32_t parsed_version{raptor_index<>::version};
        archive(parsed_version);
        if (parsed_version == raptor_index<>::version)
        {
            try
            {
                archive(window_size_);
                archive(shape_);
                archive(parts_);
                archive(compressed_);
                /*if ((data_layout_mode == seqan3::data_layout::compressed && !compressed_)
                    || (data_layout_mode == seqan3::data_layout::uncompressed && compressed_))
                {
                    throw sharg::parser_error{"Data layouts of serialised and specified index differ."};
                }*/
                archive(bin_path_);
                archive(ixf_);
            }
            catch (std::exception const & e)
            {
                //throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
            }
        }
        else
        {
            //throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
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
                archive(parts_);
                archive(compressed_);
                archive(bin_path_);
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
            //throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
        }
    }
    //!\endcond
};

} 
