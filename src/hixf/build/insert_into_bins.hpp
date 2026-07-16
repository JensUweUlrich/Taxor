#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file insert_into_bins.hpp
 * \brief Provides the hixf::insert_into_bins overloads for distributing hashes across in-memory technical bins.
 */

/*!\brief Distributes `hashes` evenly across `number_of_bins` consecutive technical bins, starting at `bin_index`.
 * \param hashes The hashes to distribute.
 * \param ixf_bins The per-bin hash sets to insert into (indexed by absolute technical bin index).
 * \param number_of_bins The number of consecutive bins (starting at `bin_index`) to spread `hashes` across.
 * \param bin_index The first technical bin index to insert into.
 * \details
 *
 * Automatically performs naive splitting if `number_of_bins > 1`: `hashes` is chunked into `number_of_bins`
 * roughly equal-sized pieces (via seqan3::views::chunk) and chunk `i` is inserted into bin `bin_index + i`.
 */
void insert_into_bins(ankerl::unordered_dense::set<size_t> const & hashes,
                      std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins,
                     size_t const number_of_bins,
                     size_t const bin_index);

/*!\brief Computes the hashes for a chopper pack record and inserts them into its (single) technical bin.
 * \param arguments The build arguments used for hashing.
 * \param record The chopper pack record to hash; must occupy exactly one technical bin (no splitting).
 * \param ixf_bins The per-bin hash sets to insert into.
 * \deprecated No current call sites in this codebase (the only call site is commented out in
 *             hierarchical_build.cpp); superseded by computing hashes explicitly and calling the other overload.
 */
void insert_into_bins(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins);

}