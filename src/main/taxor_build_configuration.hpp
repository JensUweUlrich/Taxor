#pragma once

#include <filesystem>

/*!\file taxor_build_configuration.hpp
 * \brief Defines taxor::build::configuration, the struct populated from the
 * `taxor build` CLI options and threaded through taxor_build.cpp's build pipeline.
 */

namespace taxor::build
{

/*!\brief Configuration for the `taxor build` subcommand, filled in by
 * taxor::build::set_up_subparser_layout() from CLI options and consumed
 * throughout taxor_build.cpp.
 */
struct configuration
{
    std::string input_file_name{}; // provided by user, comma-separated list; split into input_files by sanity_checks()
    std::vector<std::string> input_files{}; // derived from input_file_name (one entry per taxonomy TSV file)
    std::string input_sequence_folder{}; // provided by user, comma-separated list; split into input_folders by sanity_checks()
    std::vector<std::string> input_folders{}; // derived from input_sequence_folder (one entry per reference genome directory)
    std::string output_file_name{}; // provided by user
    int threads{1u};
    int kmer_size{20u};
    int window_size{20u};
    int syncmer_size{10u};
    int scaling{1u}; // downsampling factor: roughly a 1/scaling fraction of hashes are kept when > 1
    bool output_verbose_statistics{false};
    bool debug{false};
    bool use_syncmer{false}; // if true, use the syncmer scheme instead of the (open) minimizer scheme
/*
private:
    friend class cereal::access;

    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        uint32_t const version{1};

        archive(CEREAL_NVP(version));
        archive(CEREAL_NVP(input_prefix));
        archive(CEREAL_NVP(count_filename));
        archive(CEREAL_NVP(sketch_directory));
        archive(CEREAL_NVP(output_filename));
        archive(CEREAL_NVP(tmax));
        // archive(CEREAL_NVP(aggregate_by_column));
        archive(CEREAL_NVP(num_hash_functions));
        archive(CEREAL_NVP(false_positive_rate));
        archive(CEREAL_NVP(alpha));
        archive(CEREAL_NVP(max_rearrangement_ratio));
        archive(CEREAL_NVP(threads));
        archive(CEREAL_NVP(estimate_union));
        archive(CEREAL_NVP(rearrange_user_bins));
        archive(CEREAL_NVP(determine_best_tmax));
        archive(CEREAL_NVP(force_all_binnings));
    }
    */
};

} // namespace taxor::build
