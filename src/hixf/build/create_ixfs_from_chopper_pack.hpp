/ --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void create_ixfs_from_chopper_pack(build_data & data, build_arguments const & arguments);

}
