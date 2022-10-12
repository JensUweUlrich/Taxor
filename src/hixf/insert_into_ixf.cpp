
#include <seqan3/search/views/minimiser_hash.hpp>

#include "adjust_seed.hpp"
#include "insert_into_ixf.hpp"
#include "dna4_traits.hpp"
#include "compute_hashes.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ixf(robin_hood::unordered_flat_set<size_t> & parent_hashes,
                     robin_hood::unordered_flat_set<size_t> const & hashes,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root)
{
    size_t const chunk_size = hashes.size() / number_of_bins + 1;
    size_t chunk_number{};
    // TODO : give all hashed kmers/syncmers/minimizer to given xor filter

    for (auto chunk : hashes | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        size_t bin_idx = bin_index + chunk_number;
        ++chunk_number;
        std::vector<size_t> c{};
        std::cout << "insert_into_ixf 30" << std::endl;
        std::ranges::copy(chunk, std::back_inserter(c));
        ixf.add_bin_elements(bin_idx, c);
        std::cout << "insert_into_ixf 33" << std::endl;
        for (size_t const value : chunk)
        {
        //    ixf.emplace(value, bin_idx);
            if (!is_root)
                parent_hashes.insert(value);
        }
        std::cout << "insert_into_ixf 40" << std::endl;
    }
    
}

void insert_into_ixf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf)
{
    auto const bin_index = static_cast<size_t>(record.bin_indices.back());
    robin_hood::unordered_flat_set<size_t> hashes{};
    std::cout << "before compute_hashes" << std::endl;
    compute_hashes(hashes, arguments, record);
    std::vector<size_t> h{hashes.begin(), hashes.end()};
    ixf.add_bin_elements(bin_index, h);
    std::cout << "after compute_hashes" << std::endl;

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
