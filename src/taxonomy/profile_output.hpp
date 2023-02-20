#pragma once

#include <string>

namespace taxor::taxonomy
{

struct Profile_Output
{
    std::string rank;
    std::string taxid;
    std::string taxid_string;
    std::string taxname_string;
    double percentage;

};

}