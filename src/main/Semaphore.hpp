
#ifndef semaphore_hpp
#define semaphore_hpp

#include <future>

class Semaphore {
public:
    explicit Semaphore(size_t count) : count(count) {}
    size_t getCount() const { return count; };     
    void lock() { // call before critical section
        std::unique_lock<std::mutex> lock(mutex);
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
        condition_variable.notify_one();
    }

private:
    std::mutex mutex;
    std::condition_variable condition_variable;
    size_t count;
};

template<typename T>
bool is_future_ready(const std::future<T>& f) {
    if (f.valid()) { // otherwise you might get an exception (std::future_error: No associated state)
        return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
    } else {
        return false;
    }
}

#endif