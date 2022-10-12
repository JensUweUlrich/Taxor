
#pragma once

#include <string>

#include "chopper_pack_record.hpp"

namespace hixf
{

chopper_pack_record parse_chopper_pack_line(std::string const & current_line);

} 
