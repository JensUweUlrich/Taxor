#include <fstream>
#include <iostream>
#include <sstream>
#include "taxutil.hpp"

namespace taxor::taxonomy
{

    void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines) 
    {
        std::ifstream ifs(fname);
        if (ifs.fail()) {
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