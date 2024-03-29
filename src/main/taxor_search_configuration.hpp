#pragma once

#include <string>

namespace taxor::search
{

struct configuration
{
    std::string index_file{};
    std::string query_file{};
    std::string report_file{};
    double threshold{-1.0};
    double error_rate{0.04};
    uint8_t threads{1u};
    bool output_verbose_statistics{false};
    bool debug{false};
};

}