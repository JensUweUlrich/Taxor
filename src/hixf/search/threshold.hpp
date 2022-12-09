
#pragma once

#include "threshold_parameters.hpp"
#include "kmer_model.hpp"
#include "fracminhash_model.hpp"

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
        size_t kmers_per_window = arguments.window_size - kmer_size + 1;
        if (arguments.percentage > 0.0 && arguments.percentage <= 1.0)
        {
            threshold_kind = threshold_kinds::percentage;
            threshold_percentage = arguments.percentage;
        }
        else if (kmers_per_window == 1 && !arguments.fracminhash)
        {
            error_rate = arguments.seq_error_rate;
            threshold_kind = threshold_kinds::kmer_model;
        }
        else
        {
            threshold_kind = threshold_kinds::fracminhash;
            error_rate = arguments.seq_error_rate;
        }
    }

    size_t get(size_t minimiser_count, double scaling_factor) noexcept
    {
        switch (threshold_kind)
        {
            case threshold_kinds::kmer_model:
            {
                std::cout << "kmer_model" << std::endl;
                hixf::threshold::TInterval ci = hixf::threshold::calculate_nmut_kmer_CI(error_rate, (size_t) kmer_size, minimiser_count, 0.95);
                std::cout << ci.first << "\t" << ci.second << std::endl;
                return minimiser_count - ci.second;
            }
            case threshold_kinds::fracminhash:
            {
                std::pair<double, double> cont_dist_ci = hixf::threshold::calculate_containment_index_CI(error_rate, 
                                                                                                         kmer_size, 
                                                                                                         minimiser_count, 
                                                                                                         scaling_factor, 
                                                                                                         0.95);
                std::cout << "C_low  : " << cont_dist_ci.first << std::endl;
                std::cout << "C_high : " << cont_dist_ci.second << std::endl;
                return cont_dist_ci.first * minimiser_count;
            }
            default:
            {
                return static_cast<size_t>(minimiser_count * threshold_percentage);
            }
        }
    }

private:
    
    enum class threshold_kinds
    {
        fracminhash,
        percentage,
        kmer_model,
    };

    threshold_kinds threshold_kind{threshold_kinds::percentage};
    uint8_t kmer_size;
    double threshold_percentage{};
    double error_rate{};

    

    
};

} 
