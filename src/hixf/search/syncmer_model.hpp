#pragma once

#include <utility>
#include <cstdlib>
#include <cassert>
#include <cmath>

namespace hixf::threshold
{
    // Matrix that stores minimal matching ratios for each combination of kmer size and read accuracy
    // rows correspond to read accuracies 80% <= x <= 100%
    // columns represent kmer sizes 10,12,14,...,30
    double matching_ratios[21][11] = {
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
                                    };

    double get_min_syncmer_match_ratio(size_t kmer_size, double error_rate)
    {
        assert(kmer_size % 2 == 0);
        assert(kmer_size >= 10u);
        assert(kmer_size <= 30u);
        assert(error_rate >= 0);
        assert(error_rate <= 0.2);

        size_t accuracy = static_cast<size_t>(floor((1 - error_rate) * 100));
        size_t row_index = accuracy - 80;
        size_t col_index = kmer_size - 10 - ((kmer_size-10)/2);
        return matching_ratios[row_index][col_index];
    }

}