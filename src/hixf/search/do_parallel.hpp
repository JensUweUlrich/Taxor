// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <chrono>
#include <future>
#include <vector>

/*!\file do_parallel.hpp
 * \brief A small helper that splits a range of record indices into contiguous chunks and dispatches one
 *        chunk per thread via std::async, used to parallelise search work across query/reference records.
 */
namespace hixf
{

/*!\brief Splits [0, num_records) into \p threads contiguous chunks and runs \p worker on each chunk in
 *        parallel, accumulating the elapsed wall-clock time.
 * \tparam algorithm_t Callable type invoked as `worker(start, end)` for each chunk.
 * \param worker       The callable to run per chunk; receives the chunk's [start, end) index range.
 * \param num_records  Total number of records to distribute across threads.
 * \param threads      Number of chunks/threads to split the work into.
 * \param compute_time Accumulator that the elapsed wall-clock time of this call is added to.
 *
 * Chunk sizes are `num_records / threads`, except the last chunk which additionally absorbs the remainder
 * so that every record is covered even when \p num_records is not evenly divisible by \p threads.
 */
template <typename algorithm_t>
void do_parallel(algorithm_t && worker, size_t const num_records, size_t const threads, double & compute_time)
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<decltype(std::async(std::launch::async, worker, size_t{}, size_t{}))> tasks;
    size_t const records_per_thread = num_records / threads;

    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = records_per_thread * i;
        size_t const end = i == (threads - 1) ? num_records : records_per_thread * (i + 1);
        tasks.emplace_back(std::async(std::launch::async, worker, start, end));
    }

    // get() (not wait()) so an exception raised in a worker propagates to the caller
    // instead of being silently swallowed
    for (auto && task : tasks)
        task.get();

    auto end = std::chrono::high_resolution_clock::now();
    compute_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

}
