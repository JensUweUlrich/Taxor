
#include <seqan3/search/views/minimiser_hash.hpp>

#include "adjust_seed.hpp"
#include "compute_kmers.hpp"
#include "dna4_traits.hpp"

namespace hixf
{

void compute_kmers(robin_hood::unordered_flat_set<size_t> & kmers,
                   build_arguments const & arguments,
                   chopper_pack_record const & record)
{
    if (arguments.is_minimiser)
    {
        uint64_t minimiser_value{};
        for (auto const & filename : record.filenames)
        {
            std::ifstream infile{filename, std::ios::binary};

            while (infile.read(reinterpret_cast<char *>(&minimiser_value), sizeof(minimiser_value)))
                kmers.insert(minimiser_value);
        }
    }
    else
    {
        using sequence_file_t = seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::seq>>;
        for (auto const & filename : record.filenames)
            for (auto && [seq] : sequence_file_t{filename})
                for (auto hash :
                     seq
                         | seqan3::views::minimiser_hash(arguments.shape,
                                                         seqan3::window_size{arguments.window_size},
                                                         seqan3::seed{adjust_seed(arguments.shape.count())}))
                    kmers.insert(hash);
    }
}

}
