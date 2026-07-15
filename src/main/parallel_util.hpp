#pragma once

#include <algorithm>
#include <future>
#include <vector>

namespace taxor::util
{

// Pointers into a std::map's (or any associative container's) entries, for O(1)
// random-access partitioning across threads. The container must not be structurally
// mutated (insert/erase) while these pointers are in use; in-place mutation of the
// mapped value through a pointer is fine, since that never invalidates iterators/pointers
// to *other* elements of a std::map.
template <typename Map>
std::vector<typename Map::value_type *> to_pointer_vector(Map & map)
{
    std::vector<typename Map::value_type *> out;
    out.reserve(map.size());
    for (auto & kv : map)
        out.push_back(&kv);
    return out;
}

// Launches exactly `threads` async tasks, each given (thread_index, start, end) over a
// contiguous partition of [0, n). thread_index is in [0, threads) and stable even when
// n < threads (some tasks then get an empty range), so a worker can safely index a
// per-thread partial-result vector (e.g. partials[thread_index]) with no risk of two
// threads aliasing the same slot. Uses future::get() (not wait()) so an exception raised
// in a worker propagates to the caller instead of being silently swallowed.
template <typename Worker>
void parallel_for_indexed(Worker && worker, size_t n, size_t threads)
{
    threads = std::max<size_t>(1, threads);
    size_t const chunk = n / threads;

    std::vector<std::future<void>> tasks;
    tasks.reserve(threads);
    for (size_t i = 0; i < threads; ++i)
    {
        size_t const start = chunk * i;
        size_t const end = (i + 1 == threads) ? n : chunk * (i + 1);
        tasks.emplace_back(std::async(std::launch::async, worker, i, start, end));
    }

    for (auto & t : tasks)
        t.get();
}

} // namespace taxor::util
