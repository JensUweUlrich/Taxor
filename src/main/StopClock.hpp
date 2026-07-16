#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>

#ifndef STOPCLOCK_HPP_
#define STOPLCOCK_HPP_

/*!\file StopClock.hpp
 * \brief A simple start/stop stopwatch (StopClock) for accumulating elapsed
 * wall-clock time across one or more timing rounds, plus small helper structs
 * (TimeMeasures, Durations) for aggregating per-read timing measurements.
 */

/*!\brief A start/stop stopwatch that accumulates elapsed wall-clock time across
 * one or more start()/stop() rounds.
 *
 * Time is measured using std::chrono::system_clock. Each call to start() begins
 * (or resumes) a timing round; the very first call additionally records the
 * overall begin() timestamp. Each matching stop() call adds the time elapsed
 * since the last start() to the running total, retrievable via elapsed() or
 * runtime().
 */
class StopClock
{
    public:
        using Clock     = std::chrono::system_clock;
        using TimePoint = std::chrono::system_clock::time_point;
        using Seconds   = double;

        /*!\brief Begin (or resume) timing a round.
         *
         * Records the current time as the start of the next stop() interval. On
         * the very first call to start(), also records this time as the overall
         * begin() timestamp.
         */
        void start()
        {
            m_beginRound = Clock::now();

            if ( m_firstStart )
            {
                m_begin      = m_beginRound;
                m_firstStart = false;
            }
        }

        /*!\brief End the current timing round and add its duration to the
         * accumulated runtime.
         *
         * Records the current time as end() and adds the time elapsed since the
         * matching start() call to the running total returned by elapsed()/
         * runtime().
         */
        void stop()
        {
            m_end = Clock::now();
            m_runTime += m_end - m_beginRound;
        }

        /*!\brief Retroactively add an already-elapsed duration to the accumulated
         * runtime, without requiring a matching stop() call.
         *
         * Shifts the current round's start time backwards by \p elapsed
         * (truncated to whole seconds) and adds \p elapsed directly to the
         * running total.
         * \param elapsed Duration, in seconds, to add to the accumulated runtime.
         */
        void decrementStart(std::chrono::duration< Seconds > elapsed)
        {
            m_beginRound -= std::chrono::duration_cast<std::chrono::seconds>(elapsed);
            m_runTime += elapsed;
        }

        /*!\brief Explicitly set the start time of the current round.
         * \param b Time point to use as the start of the current round.
         */
        void setBegin(TimePoint b)
        {
            m_beginRound = b;
        }

        /*!\brief Return the accumulated runtime.
         * \return Total accumulated runtime, in seconds, across all completed
         * start()/stop() rounds (and any decrementStart() calls).
         */
        Seconds elapsed() const noexcept
        {
            return m_runTime.count();
        }

        /*!\brief Return the time point of the very first start() call.
         * \return The overall begin time point.
         */
        TimePoint begin() const noexcept
        {
            return m_begin;
        }

        /*!\brief Return the time point of the most recent stop() call.
         * \return The most recent end time point.
         */
        TimePoint end() const noexcept
        {
            return m_end;
        }

        /*!\brief Return the accumulated runtime as a duration.
         * \return Total accumulated runtime as a std::chrono::duration<Seconds>.
         */
        std::chrono::duration< Seconds > runtime()
        {
            return m_runTime;
        }

    private:
        bool                             m_firstStart{ true };
        TimePoint                        m_begin;
        TimePoint                        m_beginRound;
        TimePoint                        m_end;
        std::chrono::duration< Seconds > m_runTime{ 0.0 };
};

/*!\brief Stream-insertion operator that formats a StopClock::TimePoint as a
 * local date/time string (using the "%c" format).
 * \param stream Output stream to write to.
 * \param timepoint Time point to format and write.
 * \return The same \p stream, for chaining.
 */
template < typename Stream >
inline Stream& operator<<( Stream& stream, const StopClock::TimePoint& timepoint )
{
    const auto time = std::chrono::system_clock::to_time_t( timepoint );
    stream << std::put_time( std::localtime( &time ), "%c" );

    return stream;
}

/*!\brief Groups the StopClock instances used to time the stages of processing
 * a single read (reading, basecalling, classification).
 */
struct TimeMeasures
{
    StopClock timeCompleteRead{};
    StopClock timeBasecallRead{};
    StopClock timeClassifyRead{};
};


/*!\brief Accumulates summed durations, in seconds, split by classified vs.
 * unclassified reads and by processing stage (basecalling vs. classification).
 */
struct Durations {
    double completeClassified = 0;
    double completeUnclassified = 0;
    double basecalling = 0;
    double classification = 0;
};


#endif /* STOPCLOCK_HPP_ */