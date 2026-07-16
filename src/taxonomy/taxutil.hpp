#pragma once

#include <string>
#include <vector>

/*!\file taxutil.hpp
 * \brief Declares a small generic tab-separated-values reading utility used by
 *        the taxonomy parser (see parse_ncbi_taxonomy.cpp).
 */

namespace taxor::taxonomy
{

    /*!\brief Read a tab-separated file into a vector of tokenised lines.
     * \param fname Path of the file to read.
     * \param[out] lines Filled with one entry per input line, each entry being the
     *                   line split on tab ('\t') characters. Left untouched if the
     *                   file cannot be opened.
     */
    void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines);
}