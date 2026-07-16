#include <fstream>
#include <iostream>
#include <sstream>
#include "taxutil.hpp"

/*!\file taxutil.cpp
 * \brief Implements a small generic tab-separated-values reading utility used
 *        by the taxonomy parser.
 */

namespace taxor::taxonomy
{

    /*!\brief Read a tab-separated file into a vector of tokenised lines.
     * \param fname Path of the file to read.
     * \param[out] lines Filled with one entry per input line, each entry being the
     *                   line split on tab ('\t') characters. Left untouched if the
     *                   file cannot be opened.
     */
    void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines)
    {
        std::ifstream ifs(fname);
        if (ifs.fail()) {
            // Note: only a generic message is printed; the offending fname is not
            // included in the diagnostic.
            std::cerr << "error" << std::endl;
            return;
        }
        std::string line;
        while (getline(ifs, line)) {
            std::stringstream ss(line);
            std::vector<std::string> item;
            std::string tmp;
            while (getline(ss, tmp, '\t')) {
                item.push_back(tmp);
            }
            lines.push_back(item);
        }
    }
}