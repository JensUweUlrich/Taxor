#pragma once

#include <filesystem>

namespace taxor::build
{

struct configuration
{
    std::string input_file_name{}; // provided by user
    std::string input_sequence_folder{};
    std::string output_file_name{}; // provided by user
    uint8_t threads{1u};
    uint8_t kmer_size{20u};
    uint8_t syncmer_size{10u};
    bool output_verbose_statistics{false};
    bool debug{false};
    bool use_syncmer{false};
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
