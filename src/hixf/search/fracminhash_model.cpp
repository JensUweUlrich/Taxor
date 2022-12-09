
#include "fracminhash_model.hpp"
#include <cmath>
#include "gaussian_inverse.hpp"

namespace hixf::threshold
{
    
    double expected_containment_index(double r, size_t kmer_size)
    {
        return pow((1.0 - r), kmer_size);
    }

    double variance_containment_index(double r, size_t kmer_size, size_t kmer_count, double scaling_factor)
    {
        double term3 = variance_nmut_kmer(r, kmer_size, kmer_count) / pow(kmer_count, 2);
        double term2 = kmer_count * expected_nmut_kmer(r, kmer_size, kmer_count) - expected_nmut_kmer_squared(r, kmer_size, kmer_count);
        double denominator = scaling_factor * pow(kmer_count, 3) * pow(1.0 - pow(1.0 - scaling_factor, kmer_count),2);
        double term1 = (1.0 - scaling_factor) / denominator;
        return term1 * term2 + term3;
    }

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
