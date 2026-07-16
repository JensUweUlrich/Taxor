
#pragma once

#include <filesystem>
#include <fstream>
#include <mutex>

/*!\file sync_out.hpp
 * \brief A mutex-guarded output-file stream wrapper that allows multiple search worker threads to write
 *        results to the same output file concurrently without interleaving/corrupting each other's output.
 */
namespace hixf
{

/*!\brief Thread-safe wrapper around std::ofstream; serialises concurrent writes with an internal mutex.
 */
class sync_out
{
public:
    sync_out() = default;
    sync_out(sync_out const &) = default;
    sync_out & operator=(sync_out const &) = default;
    sync_out(sync_out &&) = default;
    sync_out & operator=(sync_out &&) = default;
    ~sync_out() = default;

    /*!\brief Opens the underlying output file at the given path.
     * \param path Path of the file to write to.
     */
    sync_out(std::filesystem::path const & path) : file{path}
    {}

    /*!\brief Writes \p data to the file, holding the write mutex for the duration of the write.
     * \param data Value to write; forwarded to the underlying std::ofstream's operator<<.
     */
    template <typename t>
    void write(t && data)
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << std::forward<t>(data);
    }

    /*!\brief Writes \p data to the file, holding the write mutex for the duration of the write.
     * \param data Value to write; forwarded to the underlying std::ofstream's operator<<.
     */
    template <typename t>
    void operator<<(t && data) // Cannot return a reference to itself since multiple threads write in the meantime.
    {
        std::lock_guard<std::mutex> lock(write_mutex);
        file << std::forward<t>(data);
    }

private:
    std::ofstream file;
    std::mutex write_mutex;
};

}
