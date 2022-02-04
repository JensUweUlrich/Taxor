#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>

#ifndef STOPCLOCK_HPP_
#define STOPLCOCK_HPP_

class StopClock
{
    public:
        using Clock     = std::chrono::system_clock;
        using TimePoint = std::chrono::system_clock::time_point;
        using Seconds   = double;

        void start()
        {
            m_beginRound = Clock::now();

            if ( m_firstStart )
            {
                m_begin      = m_beginRound;
                m_firstStart = false;
            }
        }

        void stop()
        {
            m_end = Clock::now();
            m_runTime += m_end - m_beginRound;
        }

        void decrementStart(std::chrono::duration< Seconds > elapsed)
        {
            m_beginRound -= std::chrono::duration_cast<std::chrono::seconds>(elapsed);
            m_runTime += elapsed;
        }

        void setBegin(TimePoint b)
        {
            m_beginRound = b;
        }

        Seconds elapsed() const noexcept
        {
            return m_runTime.count();
        }

        TimePoint begin() const noexcept
        {
            return m_begin;
        }

        TimePoint end() const noexcept
        {
            return m_end;
        }

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

template < typename Stream >
inline Stream& operator<<( Stream& stream, const StopClock::TimePoint& timepoint )
{
    const auto time = std::chrono::system_clock::to_time_t( timepoint );
    stream << std::put_time( std::localtime( &time ), "%c" );

    return stream;
}

struct TimeMeasures
{
    StopClock timeCompleteRead{};
    StopClock timeBasecallRead{};
    StopClock timeClassifyRead{};
};


struct Durations {
    double completeClassified = 0;
    double completeUnclassified = 0;
    double basecalling = 0;
    double classification = 0;
};


#endif /* STOPCLOCK_HPP_ */