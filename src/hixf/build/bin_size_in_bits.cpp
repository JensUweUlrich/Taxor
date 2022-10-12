
#include <cmath>

#include "bin_size_in_bits.hpp"

namespace hixf
{
size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored)
{
    double const numerator{-static_cast<double>(number_of_kmers_to_be_stored * arguments.hash)};
    double const denominator{std::log(1 - std::exp(std::log(arguments.fpr) / arguments.hash))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

} 
