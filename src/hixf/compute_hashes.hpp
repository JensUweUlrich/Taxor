
#pragma once

#include "robin_hood.h"

#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

void compute_hashes(robin_hood::unordered_flat_set<size_t> & kmers,
                   build_arguments const & arguments,
                   chopper_pack_record const & record);

} 