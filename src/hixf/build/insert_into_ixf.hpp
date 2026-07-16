
#pragma once

//#include "robin_hood.h"
#define __STDC_FORMAT_MACROS
#include <ankerl/unordered_dense.h>
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

/*!\file insert_into_ixf.hpp
 * \brief Provides the hixf::insert_into_ixf overloads for inserting hashes directly into an interleaved XOR
 *        filter's technical bins.
 */

/*!\brief Distributes `hashes` evenly across `number_of_bins` consecutive technical bins of `ixf`, starting at
 *        `bin_index`, inserting each chunk directly into the filter.
 * \param parent_hashes Output set that inserted hashes are additionally collected into when `is_root` is false
 *                       (so the caller can propagate them to the parent level's filter).
 * \param hashes The hashes to distribute and insert.
 * \param number_of_bins The number of consecutive bins (starting at `bin_index`) to spread `hashes` across.
 * \param bin_index The first technical bin index to insert into.
 * \param ixf The filter to insert into.
 * \param is_root Whether the caller is building the root filter; if false, inserted hashes are also recorded in
 *                `parent_hashes`.
 * \details
 *
 * Automatically performs naive splitting if `number_of_bins > 1`. Note: inserting into bins of an
 * seqan3::interleaved_xor_filter that has already had other bins built can be problematic/unsupported by the
 * underlying filter (see the inline warning printed on failure in the implementation).
 * \deprecated No current call sites outside of the deprecated hixf::construct_ixf overload; the current build
 *             path instead accumulates hashes per bin and constructs the filter in one shot
 *             (see hixf::construct_ixf(build_data&, ...)).
 */
void insert_into_ixf(ankerl::unordered_dense::set<size_t> & parent_kmers,
                     ankerl::unordered_dense::set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root);

/*!\brief Computes the hashes for a chopper pack record and inserts them into its (single) technical bin of `ixf`.
 * \param arguments The build arguments used for hashing.
 * \param record The chopper pack record to hash; must occupy exactly one technical bin (no splitting).
 * \param ixf The filter to insert into.
 * \deprecated No current call sites in this codebase.
 */
void insert_into_ixf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf);

} // namespace raptor::hibf
