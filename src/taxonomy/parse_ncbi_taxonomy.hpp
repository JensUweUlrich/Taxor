#pragma once

#include "Species.hpp"

/*!\file parse_ncbi_taxonomy.hpp
 * \brief Declares the parser that turns an NCBI/RefSeq-style taxonomy input
 *        file (as consumed by `taxor build`) into a list of taxor::taxonomy::Species.
 */

namespace taxor::taxonomy
{

    /*!\brief Parse a tab-separated NCBI/RefSeq taxonomy file into a list of Species.
     * \param filepath Path to the tab-separated taxonomy input file. Expected columns
     *                 (0-indexed): 0 accession id, 1 taxid, 2 reference file/ftp path,
     *                 3 organism name (optional), 4 GTDB-style lineage name string
     *                 (optional, e.g. "k__...;p__...;...;s__..."), 5 matching lineage
     *                 taxid string (optional).
     * \return A vector of Species, one per input row that has a usable file path;
     *         rows whose path column contains no path component (e.g. ftp_path "na"
     *         for withdrawn/suppressed assemblies) are skipped.
     *
     * See parse_refseq_taxonomy_file() in the corresponding .cpp file for details on
     * how the reference file stem is derived and which rows get skipped.
     */
    std::vector<Species> parse_refseq_taxonomy_file(std::string const filepath);

}