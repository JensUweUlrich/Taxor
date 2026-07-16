
#include "fracminhash_model.hpp"
#include <cmath>
#include "gaussian_inverse.hpp"

/*!\file fracminhash_model.cpp
 * \brief Implements the FracMinHash containment-index confidence-interval model (see fracminhash_model.hpp).
 */
namespace hixf::threshold
{

    /*!\brief Calculates the expected containment index based on kmer size and assumed sequencing error rate.
     * \param r         Assumed sequencing/mutation error rate.
     * \param kmer_size Size of kmers used.
     * \return (1 - r)^kmer_size, the probability that a k-mer is unaffected by errors/mutations.
     */
    double expected_containment_index(double r, size_t kmer_size)
    {
        return pow((1.0 - r), kmer_size);
    }

    /*!\brief Calculates the variance of the containment index estimator.
     * \param r              Assumed sequencing/mutation error rate.
     * \param kmer_size      Size of kmers used.
     * \param kmer_count     Number of kmers of a given query sequence.
     * \param scaling_factor FracMinHash scaling factor (fraction of hash space retained in the sketch).
     * \return Variance of the containment index.
     */
    double variance_containment_index(double r, size_t kmer_size, size_t kmer_count, double scaling_factor)
    {
        double term3 = variance_nmut_kmer(r, kmer_size, kmer_count) / pow(kmer_count, 2);
        double term2 = kmer_count * expected_nmut_kmer(r, kmer_size, kmer_count) - expected_nmut_kmer_squared(r, kmer_size, kmer_count);
        double denominator = scaling_factor * pow(kmer_count, 3) * pow(1.0 - pow(1.0 - scaling_factor, kmer_count),2);
        double term1 = (1.0 - scaling_factor) / denominator;
        return term1 * term2 + term3;
    }

    /*!\brief Calculates the confidence interval of the containment index.
     *
     * Uses a normal approximation: the interval is the expected containment index plus/minus a
     * z-score (obtained via the inverse normal CDF) scaled by the standard deviation of the
     * containment index estimator.
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
                                                             double confidence)
    {
        double z_alpha = NormalCDFInverse(1.0 - (1.0 - confidence) / 2.0);
        double clow = expected_containment_index(r, kmer_size) - z_alpha * sqrt(variance_containment_index(r, kmer_size, kmer_count, scaling_factor));
        double chigh = expected_containment_index(r, kmer_size) + z_alpha * sqrt(variance_containment_index(r, kmer_size, kmer_count, scaling_factor));
        return std::make_pair(clow, chigh);
    }

} // namespace hixf::threshold
