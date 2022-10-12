#pragma once

#include <string>

#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void read_chopper_pack_file(build_data & data, std::string const & chopper_pack_filename);

} 
