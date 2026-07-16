// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/io/sequence_file/input.hpp>

namespace hixf
{

/*!\file dna4_traits.hpp
 * \brief Provides hixf::dna4_traits.
 */

/*!\brief seqan3::sequence_file_input traits that read sequences using the 4-letter (A/C/G/T) DNA alphabet.
 * \details
 *
 * The default seqan3 traits use the 5-letter seqan3::dna5 alphabet (which also allows 'N'). Reading with
 * seqan3::dna4 instead is used where ambiguous bases are not needed and a more compact in-memory
 * representation is preferred (e.g. minimiser hashing in hixf::compute_hashes).
 */
struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    //!\brief The alphabet used to represent sequence characters; overrides the seqan3::dna5 default.
    using sequence_alphabet = seqan3::dna4;
};

}
