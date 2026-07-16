

#pragma once

#include <utility>
#include <cstdlib>

/*!\file kmer_model.hpp
 * \brief Declares the k-mer-count statistical model used to compute confidence intervals for the number of
 *        "mutated" (i.e. non-matching) k-mers of a query sequence under an assumed sequencing error rate.
 *
 * Based on "Statistics of kmers from a sequence undergoing a simple mutation process without spurious
 * matches" by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P. This model is used by
 * hixf::threshold::threshold when every window of the query yields exactly one k-mer (i.e. no minimiser
 * subsampling is applied), so the full k-mer-based mutation statistics apply directly.
 */
namespace hixf::threshold
{

    typedef std::pair<size_t, size_t> TInterval;

    /*!\brief Calculates the expected number of erroneous kmers based on kmer size, kmer count and error rate.
     *
     * Based on "Statistics of kmers from a sequence undergoing a simple mutation process without spurious
     * matches" by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
     * \param r          Assumed sequencing/mutation error rate (per-base probability of a mismatch).
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return Expected number of erroneous (mutated) kmers.
     */
    double expected_nmut_kmer(double r, size_t kmer_size, size_t kmer_count);

    /*!\brief Calculates the second moment (square of expectation plus variance) of the number of erroneous kmers.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return E[N^2] for the number of mutated kmers N, i.e. E[N]^2 + Var[N].
     */
    double expected_nmut_kmer_squared(double r, size_t kmer_size, size_t kmer_count);

    /*!\brief Calculates the variance of erroneous kmers based on kmer size, kmer count and error rate.
     *
     * Based on "Statistics of kmers from a sequence undergoing a simple mutation process without spurious
     * matches" by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \return Variance in the number of erroneous kmers.
     */
    double variance_nmut_kmer(double r, size_t kmer_size, size_t kmer_count);


    /*!\brief Calculates the confidence interval for the number of erroneous kmers based on kmer size, kmer
     *        count and error rate.
     *
     * Based on "Statistics of kmers from a sequence undergoing a simple mutation process without spurious
     * matches" by Blanca, A., Harris, R., Koslicki, D. and Medvedev, P.
     * \param r          Assumed sequencing/mutation error rate.
     * \param kmer_size  Size of kmers used.
     * \param kmer_count Number of kmers of a given query sequence.
     * \param confidence Significance level, e.g. 0.95 for a 95% confidence interval.
     * \return Pair (low, high) giving the integer-rounded lower and upper bounds of the confidence interval
     *         for the number of erroneous kmers.
     */
    TInterval calculate_nmut_kmer_CI(double r, size_t kmer_size, size_t kmer_count, double confidence);
}
