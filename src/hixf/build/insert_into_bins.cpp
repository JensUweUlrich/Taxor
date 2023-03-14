#include "insert_into_bins.hpp"
#include <seqan3/utility/views/chunk.hpp>
#include "compute_hashes.hpp"

namespace hixf
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_bins(ankerl::unordered_dense::set<size_t> const & hashes,
                      std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins,
                     size_t const number_of_bins,
                     size_t const bin_index)
{
    size_t const chunk_size = hashes.size() / number_of_bins + 1;
    size_t chunk_number{};
   
    for (auto chunk : hashes | seqan3::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        size_t bin_idx = bin_index + chunk_number;
        ++chunk_number;
        
        for (auto value : chunk)
        {
        //    ixf.emplace(value, bin_idx);
            ixf_bins.at(bin_idx).insert(value);
        }
    }

}

void insert_into_bins(build_arguments const & arguments,
                     chopper_pack_record const & record,
                     std::vector<ankerl::unordered_dense::set<size_t>> & ixf_bins)
{
    auto const bin_index = static_cast<size_t>(record.bin_indices.back());
    ankerl::unordered_dense::set<size_t> hashes{};
    compute_hashes(hashes, arguments, record);
    for (size_t const value : hashes)
    {
        ixf_bins.at(bin_index).insert(value);
    }
    
}

}