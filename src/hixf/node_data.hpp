
#pragma once

#include <vector>

#include <lemon/list_graph.h>

#include "chopper_pack_record.hpp"

namespace hixf
{

struct node_data // rename:ibf_data? or ibf_node_data
{
    size_t parent_bin_index{};
    size_t max_bin_index{};
    size_t number_of_technical_bins{};
    lemon::ListDigraph::Node favourite_child{lemon::INVALID};
    std::vector<chopper_pack_record> remaining_records{}; // non-merged bins (either split or single)

    bool operator==(node_data const & rhs) const
    {
        bool res = std::tie(parent_bin_index, max_bin_index, number_of_technical_bins, favourite_child)
                == std::tie(rhs.parent_bin_index, rhs.max_bin_index, rhs.number_of_technical_bins, rhs.favourite_child);

        if (remaining_records.size() != rhs.remaining_records.size())
            return false;

        for (size_t i = 0; i < remaining_records.size(); ++i)
            res &= (remaining_records[i] == rhs.remaining_records[i]);

        return res;
    }

    bool operator!=(node_data const & rhs) const
    {
        return !(*this == rhs);
    }
};

} 
