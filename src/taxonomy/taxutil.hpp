#pragma once

#include <string>
#include <vector>

namespace taxor::taxonomy
{

    void read_tsv(std::string fname, std::vector<std::vector<std::string> > &lines);
}