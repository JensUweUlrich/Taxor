#pragma once

#include <string>

namespace taxor::search
{

struct configuration
{
    std::string index_file{};
    std::vector<std::string> index_file_list{};
    std::string query_file{};
    std::vector<std::string> query_file_list{};
    std::string report_file{};
    double threshold{-1.0};
    double error_rate{0.04};
    uint8_t threads{1u};
    bool output_verbose_statistics{false};
    bool debug{false};
};

}