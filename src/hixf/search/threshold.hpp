
#pragma once

#include "threshold_parameters.hpp"

namespace hixf
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
        else if (kmers_per_window == 1)
        {
            error_rate = arguments.seq_error_rate;
            threshold_kind = threshold_kinds::confidence_interval;
        }
        else
        {
            threshold_kind = threshold_kinds::probabilistic;
            size_t const kmers_per_pattern = arguments.pattern_size - kmer_size + 1;
            minimal_number_of_minimizers = kmers_per_pattern / kmers_per_window;
            maximal_number_of_minimizers = arguments.pattern_size - arguments.window_size + 1;
        }
    }

    size_t get(size_t minimiser_count) noexcept
    {
        switch (threshold_kind)
        {
            case threshold_kinds::confidence_interval:
            {
                TInterval ci = calculateCI(error_rate, kmer_size, minimiser_count, 0.95);
                return minimiser_count - ci.second;
            }
            case threshold_kinds::percentage:
                return static_cast<size_t>(minimiser_count * threshold_percentage);
            case threshold_kinds::probabilistic:
            {
                size_t const index = std::clamp(minimiser_count, minimal_number_of_minimizers, maximal_number_of_minimizers)
                                - minimal_number_of_minimizers;
                return precomp_thresholds[index] + precomp_correction[index];
            }
            default:
            {
                
            }
        }
    }

private:
    typedef std::pair<size_t, size_t> TInterval;
    enum class threshold_kinds
    {
        probabilistic,
        percentage,
        confidence_interval,
    };

    threshold_kinds threshold_kind{threshold_kinds::confidence_interval};
    std::vector<size_t> precomp_correction{};
    std::vector<size_t> precomp_thresholds{};
    uint8_t kmer_size;
    size_t minimal_number_of_minimizers{};
    size_t maximal_number_of_minimizers{};
    //size_t query_length{};
    double threshold_percentage{};
    double error_rate{};

    /**
    *   calculate die confidence interval for the number of errorneous kmers based on read length, kmer size and error rate
    *   based on "statistics of kmers from a sequence undergoing a simple mutation process without spurious matches" 
    *   by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
    *   @r          : assumed sequencing error rate 
    *   @kmer_size  : size of kmers used
    *   @kmer_count : number of kmers of a given query sequence
    *   @confidence : significance level, e.g. 0.95 for a 95% confidence interval
    *   @return     : pair of integer values representing the boundaries of the confidence interval 
    */
    TInterval calculateCI(double r, uint8_t kmer_size, size_t kmer_count, double confidence)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        // expected number of errorneous/mutated kmers in sequence of length readlen
        //double Nmut = L * q; //@warning: unused variable
        // compute variance
        double varN = (double)kmer_count * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q);
        double alpha = 1 - confidence;
        
        double z = NormalCDFInverse(1.0 - alpha / 2.0);
        size_t low = static_cast<size_t>(floor(kmer_count * q - z * sqrt(varN)));
        size_t high = static_cast<size_t>(ceil(kmer_count * q + z * sqrt(varN)));
        TInterval ci{ low , high };
        return ci;
    }

    /**
    *   Abramowitz-Stegun-Approximation for the inverse normal CDF
    */
    inline double RationalApproximation(double t)
    {
        // Abramowitz and Stegun formula 26.2.23.
        // The absolute value of the error should be less than 4.5 e-4.
        double c[] = { 2.515517, 0.802853, 0.010328 };
        double d[] = { 1.432788, 0.189269, 0.001308 };

        return t - ((c[2] * t + c[1]) * t + c[0]) /
            (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
    }

    /**
    *   approximates the value of the inverse normal cumulative distribution function
    *   @p      : probability, has to be between 0 and 1
    *   @return : z score
    */
    inline double NormalCDFInverse(double p)
    {
        
        if (p <= 0.0 || p >= 1.0)
        {
            std::stringstream os;
            os << "Invalid input argument (" << p
                << "); must be larger than 0 but less than 1.";
            throw std::invalid_argument(os.str());
        }

        // See article above for explanation of this section.
        if (p < 0.5)
        {
            // F^-1(p) = - G^-1(p)
            return -RationalApproximation(sqrt(-2.0 * log(p)));
        }
        else
        {
            // F^-1(p) = G^-1(1-p)
            //std::cout << "RationalApproximation(sqrt(-2.0 * log(1.0 - p))) is: " << RationalApproximation(sqrt(-2.0 * log(1.0 - p))) << std::endl;
            return RationalApproximation(sqrt(-2.0 * log(1.0 - p)));
            
        }
    }
};

} 
