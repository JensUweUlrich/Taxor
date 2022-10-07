
#include <seqan3/search/views/minimiser_hash.hpp>

#include "adjust_seed.hpp"
#include "insert_into_ixf.hpp"
#include "dna4_traits.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ixf(robin_hood::unordered_flat_set<size_t> & parent_kmers,
                     robin_hood::unordered_flat_set<size_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};
    // TODO : give all hashed kmers/syncmers/minimizer to given xor filter
/*
    for (auto chunk : kmers | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan3::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
        {
            ixf.emplace(value, bin_idx);
            if (!is_root)
                parent_kmers.insert(value);
        }
    }
    */
}

void insert_into_ibf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf)
{
    //auto const bin_index = seqan3::bin_index{static_cast<size_t>(record.bin_indices.back())};
    // TODO: same here as above
/*    if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : record.filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                ibf.emplace(minimiser_value, bin_index);
        }
    }
    else
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;

        auto hash_view = seqan3::views::minimiser_hash(arguments.shape,
                                                       seqan3::window_size{arguments.window_size},
                                                       seqan3::seed{adjust_seed(arguments.shape.count())});

        for (auto const & filename : record.filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash : seq | hash_view)
                    ibf.emplace(hash, bin_index);
    }
    */
}

} 
