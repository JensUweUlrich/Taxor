#pragma once

#include <string>

namespace taxor::profile
{

struct configuration
{
    std::string search_file{};
    std::string taxonomy_file{};
    std::string report_file{};
    //double threshold{-1.0};
    //double error_rate{0.04};
    uint8_t threads{1u};
    bool output_verbose_statistics{false};
    bool debug{false};
};

}