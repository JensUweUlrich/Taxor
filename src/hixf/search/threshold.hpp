
#pragma once

#include "threshold_parameters.hpp"
#include "kmer_model.hpp"
#include "fracminhash_model.hpp"
#include "syncmer_model.hpp"

namespace hixf::threshold
{

class threshold
{
public:
    threshold() = default;
    threshold(threshold const &) = default;
    threshold & operator=(threshold const &) = default;
    threshold(threshold &&) = default;
    threshold & operator=(threshold &&) = default;
    ~threshold() = default;

    threshold(threshold_parameters & arguments)
    {
        kmer_size = arguments.kmer_size;
        error_rate = arguments.seq_error_rate;
        size_t kmers_per_window = arguments.window_size - kmer_size + 1;
        //std::cout << arguments.percentage << std::endl << std::flush;
        if (arguments.percentage > 0.0 && arguments.percentage <= 1.0)
        {
            threshold_kind = threshold_kinds::percentage;
            threshold_percentage = arguments.percentage;
            std::cout << "use percentage-model\t" << arguments.percentage << std::endl << std::flush;
        }
        else if (arguments.use_syncmer)
        {
            threshold_kind = threshold_kinds::syncmer_model;
            std::cout << "use syncmer-model" << std::endl << std::flush;
        }
        else if (kmers_per_window == 1 && !arguments.fracminhash)
        {
            std::cout << "use kmer-model" << std::endl << std::flush;
            threshold_kind = threshold_kinds::kmer_model;
        }
        else
        {
            threshold_kind = threshold_kinds::fracminhash;
            std::cout << "use frac minhash" << std::endl << std::flush;
        }
    }

    size_t get(size_t minimiser_count, double scaling_factor) noexcept
    {
        size_t fp_correction = minimiser_count * 0.0039;
        //std::cout << error_rate << std::endl << std::flush;
        switch (threshold_kind)
        {
            case threshold_kinds::syncmer_model:
            {
                double syncmer_match_ratio = hixf::threshold::get_min_syncmer_match_ratio(kmer_size, error_rate);
                return static_cast<size_t>(minimiser_count * syncmer_match_ratio);
            }
            case threshold_kinds::kmer_model:
            {
                hixf::threshold::TInterval ci = hixf::threshold::calculate_nmut_kmer_CI(error_rate, (size_t) kmer_size, minimiser_count, 0.95);
                return minimiser_count - ci.second - fp_correction;
            }
            case threshold_kinds::fracminhash:
            {
                std::pair<double, double> cont_dist_ci = hixf::threshold::calculate_containment_index_CI(error_rate, 
                                                                                                         kmer_size, 
                                                                                                         minimiser_count, 
                                                                                                         scaling_factor, 
                                                                                                         0.95);
                return static_cast<size_t>(cont_dist_ci.first * minimiser_count) - fp_correction;
            }
            default:
            {
                return static_cast<size_t>(minimiser_count * threshold_percentage);
            }
        }
    }

private:
    
    // TODO: add specific syncmer model based on empirical data
    //       minimum number of OCS found for each error rate and each kmer
    //       0.90 <= r <=0.99 and 16 <= k <= 32
    enum class threshold_kinds
    {
        fracminhash,
        percentage,
        kmer_model,
        syncmer_model,
    };

    threshold_kinds threshold_kind{threshold_kinds::percentage};
    uint8_t kmer_size;
    double threshold_percentage{};
    double error_rate{};

    

    
};

} 
