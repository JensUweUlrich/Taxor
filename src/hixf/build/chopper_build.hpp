// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include "build_arguments.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void chopper_build(build_arguments const & arguments);

} // namespace raptor::hibf
