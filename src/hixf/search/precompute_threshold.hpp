
#pragma once

#include "threshold_parameters.hpp"

namespace hixf
{

[[nodiscard]] std::vector<size_t> precompute_threshold(threshold_parameters const & arguments);

} 
