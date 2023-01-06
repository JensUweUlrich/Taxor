
#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "chopper_pack_record.hpp"

namespace hixf
{

void compute_hashes(ankerl::unordered_dense::set<size_t> & kmers,
                   build_arguments const & arguments,
                   chopper_pack_record const & record);

} 