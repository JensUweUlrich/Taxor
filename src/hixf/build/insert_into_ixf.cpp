
//#include <seqan3/search/views/minimiser_hash.hpp>

//#include "adjust_seed.hpp"
#include "insert_into_ixf.hpp"
//#include "dna4_traits.hpp"
#include "compute_hashes.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ixf(ankerl::unordered_dense::set<size_t> & parent_hashes,
                     ankerl::unordered_dense::set<size_t> const & hashes,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan3::interleaved_xor_filter<> & ixf,
                     bool is_root)
{
    size_t const chunk_size = hashes.size() / number_of_bins + 1;
    size_t chunk_number{};
   
    for (auto chunk : hashes | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        size_t bin_idx = bin_index + chunk_number;
        ++chunk_number;
        std::vector<size_t> c{};
        std::ranges::copy(chunk, std::back_inserter(c));
        // adds hashes to already built bins => problematic
        bool success = ixf.add_bin_elements(bin_idx, c);
        if (success)
            std::cout << "Building IXF not successful" << std::endl;
        for (size_t const value : chunk)
        {
        //    ixf.emplace(value, bin_idx);
            if (!is_root)
                parent_hashes.insert(value);
        }
    }
    
}


void insert_into_ixf(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     seqan3::interleaved_xor_filter<> & ixf)
{
    auto const bin_index = static_cast<size_t>(record.bin_indices.back());
    ankerl::unordered_dense::set<size_t> hashes{};
    compute_hashes(hashes, arguments, record);
    std::vector<size_t> h{hashes.begin(), hashes.end()};
    ixf.add_bin_elements(bin_index, h);
    
}

} 
