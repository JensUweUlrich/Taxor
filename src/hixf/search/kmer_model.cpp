
#include "kmer_model.hpp"
#include "gaussian_inverse.hpp"
#include <cmath>
#include <iostream>

namespace hixf::threshold
{
   
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
        std::cout << varN << "\t" << q << std::endl;
        double z = NormalCDFInverse(1.0 - alpha / 2.0);
        size_t low = static_cast<size_t>(floor(kmer_count * q - z * sqrt(varN)));
        size_t high = static_cast<size_t>(ceil(kmer_count * q + z * sqrt(varN)));
        TInterval ci = std::make_pair(low , high );
        return ci;
    }


    double expected_nmut_kmer(double r, size_t kmer_size, size_t kmer_count)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        return kmer_count * q;
    }

    double variance_nmut_kmer(double r, size_t kmer_size, size_t kmer_count)
    {
        double q = 1.0 - pow(1.0 - r, kmer_size);
        double varN = (double)kmer_count * (1.0 - q) * (q * (2.0 * (double)kmer_size + (2.0 / r) - 1.0) - 2.0 * (double)kmer_size)
                        + (double)kmer_size * ((double)kmer_size - 1.0) * pow((1.0 - q), 2.0)
                        + (2.0 * (1.0 - q) / (pow(r, 2.0))) * ((1.0 + ((double)kmer_size - 1.0) * (1.0 - q)) * r - q);
        return varN;
    }

    double expected_nmut_kmer_squared(double r, size_t kmer_size, size_t kmer_count)
    {
        return pow(expected_nmut_kmer(r, kmer_size, kmer_count), 2) + variance_nmut_kmer(r, kmer_size, kmer_count);
    }
}