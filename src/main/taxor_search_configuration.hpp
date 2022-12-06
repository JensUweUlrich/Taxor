#pragma once

#include <string>

namespace taxor::search
{

struct configuration
{
    std::string index_file{};
    std::string query_file{};
    std::string report_file{};
    double threshold{0.2};
};

}