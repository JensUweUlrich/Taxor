
#pragma once

#include "threshold_parameters.hpp"

namespace hixf
{

[[nodiscard]] std::vector<size_t> precompute_correction(threshold_parameters const & arguments);

} 
