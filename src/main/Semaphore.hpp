/*!\file Semaphore.hpp
 * \brief A counting semaphore (usable as a std::scoped_lock/std::unique_lock-compatible BasicLockable via
 *        lock()/unlock()) plus a small helper for polling std::future readiness without blocking.
 *
 * \details Not currently compiled into the taxor binary — this file is only #include'd by build.hpp, which is
 *          itself not #include'd by anything reachable from main.cpp as of this writing, so this class is
 *          effectively dead code too.
 */

#ifndef semaphore_hpp
#define semaphore_hpp

#include <future>

/*!\brief A classic counting semaphore used to cap the number of concurrently running tasks (e.g. limiting how
 *        many std::async jobs run at once).
 *
 * \details Exposes lock()/unlock() so it satisfies BasicLockable and can be used directly with
 *          std::scoped_lock/std::unique_lock: lock() blocks until \c count is non-zero then decrements it
 *          (acquiring one "slot"), unlock() increments \c count back and wakes one waiter. Unlike a mutex,
 *          \c count may start above 1, allowing more than one holder at a time — that's what makes this a
 *          semaphore rather than a plain lock.
 */
class Semaphore {
public:
    explicit Semaphore(size_t count) : count(count) {}
    //!\brief Current number of free "slots" (does not itself acquire the mutex, so may be stale under contention).
    size_t getCount() const { return count; };
    void lock() { // call before critical section
        std::unique_lock<std::mutex> lock(mutex);
        // block until another thread's unlock() makes a slot available, then claim it
        condition_variable.wait(lock, [this] {
          if (count != 0) // written out for clarity, could just be return (count != 0);
              return true;
          else
              return false;
        });
        --count;
    }
    void unlock() {  // call after critical section
        std::unique_lock<std::mutex> lock(mutex);
        ++count;
        // wake exactly one waiting thread since exactly one slot was just freed
        condition_variable.notify_one();
    }

private:
    std::mutex mutex;
    std::condition_variable condition_variable;
    size_t count;
};

/*!\brief Check whether a std::future's result is already available, without blocking.
 * \tparam T The future's result type.
 * \param f The future to poll.
 * \return true if \p f is valid and its result is ready; false if \p f is invalid (e.g. default-constructed,
 *         or already retrieved via get()) or the result is not yet available.
 */
template<typename T>
bool is_future_ready(const std::future<T>& f) {
    if (f.valid()) { // otherwise you might get an exception (std::future_error: No associated state)
        return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    } else {
        return false;
    }
}

#endif