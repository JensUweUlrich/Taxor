
#include "temp_hash_file.hpp"
#include <filesystem>
#include <fstream>

namespace hixf
{

void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + ".tmp";
    auto tmp_file = std::filesystem::temp_directory_path() / ixf_tmp_name;
    std::ofstream tmp_stream(tmp_file,std::ios::out | std::ios::trunc | std::ios::binary);
    for (size_t hash : node_hashes)
    {   
        tmp_stream.write((char*)& hash, sizeof(hash));
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
        //std::cerr << ixf_tmp_name << "does not exist!" << std::endl;
        return;
    }
    std::ifstream tmp_stream(tmp_file, std::ios::in | std::ios::binary);
    size_t x;
    while (tmp_stream.read((char * ) & x, sizeof(x)))
    {
        node_hashes.emplace_back(x);
    }
    tmp_stream.close();
    
    
    //std::filesystem::remove(tmp_file);
    

}

}