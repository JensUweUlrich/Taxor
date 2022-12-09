

#pragma once

#include "kmer_model.hpp"
#include <cstdlib>

namespace hixf::threshold
{

    /**
    *   calculates the expected containment index based on read length, kmer size and error rate
    *   based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across 
    *   a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
    *   @r double          : assumed sequencing error rate 
    *   @kmer_size size_t  : size of kmers used
    *   @kmer_count size_t : number of kmers of a given query sequence
    *   @return double     : expected number of errornous kmers
    */
    double expected_containment_index(double r, size_t kmer_size);


    /**
    *   calculates the variance of the containment index based on read length, kmer size and error rate
    *   based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across 
    *   a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
    *   @r double          : assumed sequencing error rate 
    *   @kmer_size size_t  : size of kmers used
    *   @kmer_count size_t : number of kmers of a given query sequence
    *   @return double     : variance in the number of errornous kmers
    */
    double variance_containment_index(double r, size_t kmer_size, size_t kmer_count, double scaling_factor);


    /**
    *   calculates the confidence interval of the containment index based on read length, kmer size and error rate
    *   based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across 
    *   a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
    *   @r          : assumed sequencing error rate 
    *   @kmer_size  : size of kmers used
    *   @kmer_count : number of kmers of a given query sequence
    *   @confidence : significance level, e.g. 0.95 for a 95% confidence interval
    *   @return     : pair of integer values representing the boundaries of the confidence interval 
    */
    std::pair<double, double> calculate_containment_index_CI(double r, 
                                                             size_t kmer_size, 
                                                             size_t kmer_count, 
                                                             double scaling_factor, 
                                                             double confidence);
}
