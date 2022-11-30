
#include <lemon/list_graph.h> // Must be first include.

#include "update_user_bins.hpp"

namespace hixf
{

//template <seqan3::data_layout data_layout_mode>
void update_user_bins(build_data & data,
                      std::vector<int64_t> & filename_indices,
                      chopper_pack_record const & record)
{
    size_t const idx = data.request_user_bin_idx();

    std::string & user_bin_filenames = data.hixf.user_bins.filename_of_user_bin(idx);
    for (auto const & filename : record.filenames)
    {
        user_bin_filenames += filename;
        user_bin_filenames += ';';
    }
    assert(!user_bin_filenames.empty());
    user_bin_filenames.pop_back();

    std::fill_n(filename_indices.begin() + record.bin_indices.back(), record.number_of_bins.back(), idx);
}

/*
template void update_user_bins<seqan3::data_layout::uncompressed>(build_data<seqan3::data_layout::uncompressed> &,
                                                                  std::vector<int64_t> &,
                                                                  chopper_pack_record const &);

template void update_user_bins<seqan3::data_layout::compressed>(build_data<seqan3::data_layout::compressed> &,
                                                                std::vector<int64_t> &,
                                                                chopper_pack_record const &);
*/
} 