
#include <cmath>

#include "bin_size_in_bits.hpp"

/*!\file bin_size_in_bits.cpp
 * \brief Implements hixf::bin_size_in_bits.
 */

namespace hixf
{
size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored)
{
    // Standard Bloom filter sizing formula, solved for the number of bits m given the number of
    // elements n, the number of hash functions k, and the target false-positive rate fpr:
    //   m = ceil( -(n * k) / ln(1 - fpr^(1/k)) )
    double const numerator{-static_cast<double>(number_of_kmers_to_be_stored * arguments.hash)};
    double const denominator{std::log(1 - std::exp(std::log(arguments.fpr) / arguments.hash))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

}
