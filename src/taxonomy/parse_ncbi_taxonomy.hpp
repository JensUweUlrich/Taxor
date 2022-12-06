#pragma once

#include "Species.hpp"

namespace taxor::taxonomy
{
    

    std::vector<Species> parse_refseq_taxonomy_file(std::string const filepath);

}