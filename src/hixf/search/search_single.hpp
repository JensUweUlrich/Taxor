
#pragma once

#include <seqan3/search/views/minimiser_hash.hpp>

#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>
#include "do_parallel.hpp"
#include "load_index.hpp"
#include "sync_out.hpp"
#include "threshold.hpp"

namespace hixf
{

//template <typename index_t>
void search_single(search_arguments const & arguments, raptor_index<index_structure::hixf> && index)
{
    //constexpr bool is_ibf = std::same_as<index_t, raptor_index<index_structure::ibf>>
    //                     || std::same_as<index_t, raptor_index<index_structure::ibf_compressed>>;

    double index_io_time{0.0};
    double reads_io_time{0.0};
    double compute_time{0.0};

    auto cereal_worker = [&]()
    {
        load_index(index, arguments, index_io_time);
    };
    auto cereal_handle = std::async(std::launch::async, cereal_worker);

    seqan3::sequence_file_input<dna4_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{
        arguments.query_file};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    sync_out synced_out{arguments.out_file};

    {
        size_t position{};
        std::string line{};
        for (auto const & file_list : arguments.bin_path)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            for (auto const & filename : file_list)
            {
                line += filename;
                line += ',';
            }
            line.back() = '\n';
            synced_out << line;
            ++position;
        }
        synced_out << "#QUERY_NAME\tUSER_BINS\n";
    }

    hixf::threshold const thresholder{arguments.make_threshold_parameters()};

    auto worker = [&](size_t const start, size_t const end)
    {
        auto counter = [&index]()
        {
            /*if constexpr (is_ibf)
                return index.ibf().template counting_agent<uint16_t>();
            else */
                return index.ixf().membership_agent();
        }();
        std::string result_string{};
        std::vector<uint64_t> minimiser;

        // TODO: choose between minimizer and syncmers
        auto hash_adaptor = seqan3::views::minimiser_hash(arguments.shape,
                                                          seqan3::window_size{arguments.window_size},
                                                          seqan3::seed{adjust_seed(arguments.shape_weight)});

        for (auto && [id, seq] : records | seqan3::views::slice(start, end))
        {
            result_string.clear();
            result_string += id;
            result_string += '\t';

            auto minimiser_view = seq | hash_adaptor | std::views::common;
            minimiser.assign(minimiser_view.begin(), minimiser_view.end());

            size_t const minimiser_count{minimiser.size()};
            size_t const threshold = thresholder.get(minimiser_count);

            /*if constexpr (is_ibf)
            {
                auto & result = counter.bulk_count(minimiser);
                size_t current_bin{0};
                for (auto && count : result)
                {
                    if (count >= threshold)
                    {
                        result_string += std::to_string(current_bin);
                        result_string += ',';
                    }
                    ++current_bin;
                }
            }
            else
            {*/
                auto & result = counter.bulk_contains(minimiser, threshold); // Results contains user bin IDs
                for (auto && count : result)
                {
                    result_string += std::to_string(count);
                    result_string += ',';
                }
            //}

            if (auto & last_char = result_string.back(); last_char == ',')
                last_char = '\n';
            else
                result_string += '\n';
            synced_out.write(result_string);
        }
    };

    for (auto && chunked_records : fin | seqan3::views::chunk((1ULL << 20) * 10))
    {
        records.clear();
        auto start = std::chrono::high_resolution_clock::now();
        std::ranges::move(chunked_records, std::back_inserter(records));
        auto end = std::chrono::high_resolution_clock::now();
        reads_io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

        cereal_handle.wait();

        do_parallel(worker, records.size(), arguments.threads, compute_time);
    }

    // GCOVR_EXCL_START
    if (arguments.write_time)
    {
        std::filesystem::path file_path{arguments.out_file};
        file_path += ".time";
        std::ofstream file_handle{file_path};
        file_handle << "Index I/O\tReads I/O\tCompute\n";
        file_handle << std::fixed << std::setprecision(2) << index_io_time << '\t' << reads_io_time << '\t'
                    << compute_time;
    }
    // GCOVR_EXCL_STOP
}

} 
