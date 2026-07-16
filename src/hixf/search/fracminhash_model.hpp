

#pragma once

#include "kmer_model.hpp"
#include <cstdlib>

/*!\file fracminhash_model.hpp
 * \brief Declares the FracMinHash containment-index statistical model.
 *
 * FracMinHash is a scaled MinHash sketching scheme: instead of storing every k-mer hash, only hashes below
 * a fraction 1/scaling_factor of the hash space are kept. This file provides the model used to derive a
 * confidence interval for the *containment index* (the expected fraction of a query's sketch that is
 * shared with a reference) under an assumed per-base sequencing error rate, following Hera, Pierce-Ward and
 * Koslicki, "Debiasing FracMinHash and deriving confidence intervals for mutation rates across a wide range
 * of evolutionary distances". The resulting interval is used by hixf::threshold::threshold to derive the
 * minimum number of matching hashes required to call a query sequence a hit.
 */
namespace hixf::threshold
{

    /*!\brief Calculates the expected containment index based on kmer size and assumed sequencing error rate.
     *
     * Based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across
     * a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
     * \param r         Assumed sequencing/mutation error rate (probability a given base differs).
     * \param kmer_size Size of kmers used.
     * \return Expected containment index, i.e. (1 - r)^kmer_size, the probability an individual k-mer survives
     *         unmutated.
     */
    double expected_containment_index(double r, size_t kmer_size);


    /*!\brief Calculates the variance of the containment index based on kmer size, kmer count and error rate.
     *
     * Based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across
     * a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
     * \param r              Assumed sequencing/mutation error rate.
     * \param kmer_size      Size of kmers used.
     * \param kmer_count     Number of kmers of a given query sequence.
     * \param scaling_factor FracMinHash scaling factor (fraction of hash space retained in the sketch).
     * \return Variance of the containment index estimator.
     */
    double variance_containment_index(double r, size_t kmer_size, size_t kmer_count, double scaling_factor);


    /*!\brief Calculates the confidence interval of the containment index based on kmer size, kmer count and error rate.
     *
     * Based on "Debiasing FracMinHash and deriving confidence intervals for mutation rates across
     * a wide range of evolutionary distances" by Hera, M.R., Pierce-Ward, T. and Koslicki, D.
     * \param r              Assumed sequencing/mutation error rate.
     * \param kmer_size      Size of kmers used.
     * \param kmer_count     Number of kmers of a given query sequence.
     * \param scaling_factor FracMinHash scaling factor (fraction of hash space retained in the sketch).
     * \param confidence     Significance level, e.g. 0.95 for a 95% confidence interval.
     * \return Pair (low, high) giving the lower and upper bounds of the containment index confidence interval.
     */
    std::pair<double, double> calculate_containment_index_CI(double r,
                                                             size_t kmer_size,
                                                             size_t kmer_count,
                                                             double scaling_factor,
                                                             double confidence);
}
