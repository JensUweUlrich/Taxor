

#pragma once

#include <utility>
#include <cstdlib>

namespace hixf::threshold
{

    typedef std::pair<size_t, size_t> TInterval;

    /**
    *   calculates the expected number of errorneous kmers based on read length, kmer size and error rate
    *   based on "statistics of kmers from a sequence undergoing a simple mutation process without spurious matches" 
    *   by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
    *   @r double          : assumed sequencing error rate 
    *   @kmer_size size_t  : size of kmers used
    *   @kmer_count size_t : number of kmers of a given query sequence
    *   @return double     : expected number of errornous kmers
    */
    double expected_nmut_kmer(double r, size_t kmer_size, size_t kmer_count);

    double expected_nmut_kmer_squared(double r, size_t kmer_size, size_t kmer_count);

    /**
    *   calculates the variance of errorneous kmers based on read length, kmer size and error rate
    *   based on "statistics of kmers from a sequence undergoing a simple mutation process without spurious matches" 
    *   by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
    *   @r double          : assumed sequencing error rate 
    *   @kmer_size size_t  : size of kmers used
    *   @kmer_count size_t : number of kmers of a given query sequence
    *   @return double     : variance in the number of errornous kmers
    */
    double variance_nmut_kmer(double r, size_t kmer_size, size_t kmer_count);


    /**
    *   calculates the confidence interval for the number of errorneous kmers based on read length, kmer size and error rate
    *   based on "statistics of kmers from a sequence undergoing a simple mutation process without spurious matches" 
    *   by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
    *   @r          : assumed sequencing error rate 
    *   @kmer_size  : size of kmers used
    *   @kmer_count : number of kmers of a given query sequence
    *   @confidence : significance level, e.g. 0.95 for a 95% confidence interval
    *   @return     : pair of integer values representing the boundaries of the confidence interval 
    */
    TInterval calculate_nmut_kmer_CI(double r, size_t kmer_size, size_t kmer_count, double confidence);
}
