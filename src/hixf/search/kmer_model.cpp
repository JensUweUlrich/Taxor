
#include "kmer_model.hpp"
#include "gaussian_inverse.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

/*!\file kmer_model.cpp
 * \brief Implements the k-mer-count statistical model (see kmer_model.hpp) used to derive confidence
 *        intervals for the number of erroneous/mutated k-mers under an assumed sequencing error rate.
 */
namespace hixf::threshold
{

    /*!\brief Calculates the confidence interval for the number of erroneous kmers.
     *
     * Uses a normal approximation around the expected number of mutated kmers, with variance derived from
     * the mutation-process model of Blanca et al. (see file docs). \c q is the per-kmer probability that a
     * kmer contains at least one mutated base, given a per-base error rate \p r.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \param confidence Significance level, e.g. 0.95 for a 95% confidence interval.
     * \return Pair (low, high) giving the lower and upper bounds (rounded down/up and clamped at 0) of the
     *         confidence interval for the number of erroneous kmers.
     */
    TInterval calculate_nmut_kmer_CI(double r, size_t kmer_size, size_t kmer_count, double confidence)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        // expected number of errorneous/mutated kmers in sequence of length readlen
        //double Nmut = L * q; //@warning: unused variable
        // compute variance
        double varN = (double)kmer_count * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q);
        double alpha = 1 - confidence;
        //std::cout << varN << "\t" << q << std::endl;
        double z = NormalCDFInverse(1.0 - alpha / 2.0);
        size_t low = static_cast<size_t>(std::max(0.0, floor(kmer_count * q - z * sqrt(varN))));
        size_t high = static_cast<size_t>(std::max(0.0, ceil(kmer_count * q + z * sqrt(varN))));
        TInterval ci = std::make_pair(low , high );
        return ci;
    }


    /*!\brief Calculates the expected number of erroneous kmers.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return kmer_count * q, where q = 1 - (1-r)^kmer_size is the per-kmer probability of containing a
     *         mutated base.
     */
    double expected_nmut_kmer(double r, size_t kmer_size, size_t kmer_count)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        return kmer_count * q;
    }

    /*!\brief Calculates the variance of the number of erroneous kmers.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return Variance in the number of erroneous kmers, following the mutation-process model of Blanca et al.
     */
    double variance_nmut_kmer(double r, size_t kmer_size, size_t kmer_count)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        double varN = (double)kmer_count * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q);
        return varN;
    }

    /*!\brief Calculates the second moment E[N^2] of the number of erroneous kmers N.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return E[N]^2 + Var[N], i.e. the second raw moment of the number of erroneous kmers.
     */
    double expected_nmut_kmer_squared(double r, size_t kmer_size, size_t kmer_count)
    {
        return pow(expected_nmut_kmer(r, kmer_size, kmer_count), 2) + variance_nmut_kmer(r, kmer_size, kmer_count);
    }
}