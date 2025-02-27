#pragma once

#include <string>

namespace taxor::profile
{

struct configuration
{
    std::string search_file{};
    std::string binning_file{};
    std::string report_file{};
    std::string sequence_abundance_file{};
    std::string sample_id{};
    double threshold{0.001};
    //double error_rate{0.04};
    uint8_t threads{1u};
    bool output_verbose_statistics{false};
    bool debug{false};
    uint16_t em_steps{100};
};

}