
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file compute_hashes.hpp
 * \brief Provides hixf::compute_hashes.
 */

/*!\brief Computes the set of distinct hash values (syncmers or minimisers) for all sequences of a user bin.
 * \param kmers The output set that computed hashes are inserted into (not cleared; existing entries are kept).
 * \param arguments The build arguments; selects syncmer vs. minimiser hashing and the associated parameters.
 * \param record The chopper pack record whose reference sequence file(s) (record.filenames) are read and hashed.
 */
void compute_hashes(ankerl::unordered_dense::set<size_t> & kmers,
                   build_arguments const & arguments,
                   chopper_pack_record const & record);

}