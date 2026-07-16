// Timers for use in benchmarking.

#pragma once

#include <cstdint>
#include <chrono>

/*!\file timing.hpp
 * \brief Minimal timing helper providing NowNanos(), a monotonic-clock timestamp
 * in nanoseconds; used by xorfilter.hpp for lightweight benchmarking. This is a
 * separate, independent utility from the StopClock stopwatch in StopClock.hpp.
 */

/*!\brief Return the current time point of a steady (monotonic) clock, in
 * nanoseconds since an unspecified epoch.
 * \return Nanosecond timestamp suitable for measuring elapsed durations by
 * subtracting two calls to this function.
 */
::std::uint64_t NowNanos() {
  return ::std::chrono::duration_cast<::std::chrono::nanoseconds>(
             ::std::chrono::steady_clock::now().time_since_epoch())
      .count();
}
