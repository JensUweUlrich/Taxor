
#include "temp_hash_file.hpp"
#include <filesystem>
#include <fstream>

namespace hixf
{

void create_temp_hash_file(size_t const ixf_pos, std::vector<robin_hood::unordered_flat_set<size_t>> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + ".tmp";
    auto tmp_file = std::filesystem::temp_directory_path() / ixf_tmp_name;
    std::ofstream tmp_stream{tmp_file};
    for (auto bin : node_hashes)
    {   
        for (size_t h : bin)
            tmp_stream << h << " ";
    }
    tmp_stream.close();
}


void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes)
{
    
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_position) + ".tmp";
    auto tmp_file = std::filesystem::temp_directory_path() / ixf_tmp_name;
    if (!std::filesystem::exists(tmp_file))
    {
        std::cerr << ixf_tmp_name << "does not exist!" << std::endl;
        return;
    }
    std::ifstream tmp_stream{tmp_file};
    robin_hood::unordered_flat_set<size_t> hashset{};
    size_t x;
    while (tmp_stream >> x)
    {
        hashset.insert(x);
    }
    tmp_stream.close();
    std::ranges::move(hashset, std::back_inserter(node_hashes));
    //std::filesystem::remove(tmp_file);
    

}

}