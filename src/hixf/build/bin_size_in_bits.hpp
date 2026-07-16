
#pragma once

#include "build_arguments.hpp"

namespace hixf
{

/*!\file bin_size_in_bits.hpp
 * \brief Provides hixf::bin_size_in_bits.
 */

/*!\brief Computes the required size (in bits) of a single technical bin's Bloom filter for a target false-positive
 *        rate.
 * \param arguments The build arguments; used for the number of hash functions (`hash`) and the target
 *                  false-positive rate (`fpr`).
 * \param number_of_kmers_to_be_stored The number of distinct k-mers/hashes to be stored in the bin.
 * \return The minimal number of bits needed so that the Bloom filter achieves the configured false-positive rate.
 */
size_t bin_size_in_bits(build_arguments const & arguments, size_t const number_of_kmers_to_be_stored);

}
