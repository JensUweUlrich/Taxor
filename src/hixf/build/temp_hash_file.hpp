
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_data.hpp"

namespace hixf
{

/*!\file temp_hash_file.hpp
 * \brief Provides functions to spill per-node/per-bin hash sets to temporary binary files on disk and read
 *        them back, used to bound peak memory usage while building a HIXF.
 * \details
 *
 * Files are written under a `hixf_tmp` directory created in the current working directory, named either
 * `interleavedXOR_<ixf_pos>.tmp` (whole-node hash sets) or `interleavedXOR_<ixf_pos>_<bin_index>.tmp`
 * (single-bin hash sets), and store raw, unsorted `size_t` hash values back-to-back in binary form.
 */

/*!\brief Writes the combined hash set of an entire IXF node to a temporary file.
 * \param ixf_pos The IXF node's index; used to name the temp file.
 * \param node_hashes The hashes to write.
 */
void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes);

/*!\brief Writes the hash set of a single technical bin to a temporary file.
 * \param ixf_pos The IXF node's index the bin belongs to.
 * \param bin_index The technical bin's index within that node.
 * \param node_hashes The hashes to write.
 */
void create_temp_hash_file(size_t const ixf_pos, size_t const bin_index, ankerl::unordered_dense::set<size_t> &node_hashes);

/*!\brief Reads back a whole-node hash set previously written by the single-argument create_temp_hash_file().
 * \param ixf_position The IXF node's index whose temp file should be read.
 * \param node_hashes Output vector; hashes are appended here (not cleared beforehand). Left unchanged if the
 *                     file does not exist.
 * \param tmp_files Output set that the file's path is added to on a successful read, so callers can later
 *                  delete all temp files they consumed.
 */
void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files);

/*!\brief Reads back a single technical bin's hash set previously written by the two-argument
 *        create_temp_hash_file() overload.
 * \param ixf_position The IXF node's index the bin belongs to.
 * \param bin_index The technical bin's index within that node.
 * \param node_hashes Output vector; hashes are appended here (not cleared beforehand). Left unchanged if the
 *                     file does not exist.
 * \param tmp_files Output set that the file's path is added to on a successful read, so callers can later
 *                  delete all temp files they consumed.
 */
void read_from_temp_hash_file(size_t const ixf_position,
                              uint16_t const bin_index,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files);
}