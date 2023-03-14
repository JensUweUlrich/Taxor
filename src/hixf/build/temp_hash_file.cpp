
#include "temp_hash_file.hpp"
#include <filesystem>
#include <fstream>

namespace hixf
{

void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + ".tmp";
    auto tmp_dir = std::filesystem::current_path() / "hixf_tmp";
    if (!std::filesystem::exists(tmp_dir))
        std::filesystem::create_directory(tmp_dir);
    auto tmp_file = tmp_dir / ixf_tmp_name;
    std::ofstream tmp_stream(tmp_file,std::ios::out | std::ios::trunc | std::ios::binary);
    for (size_t hash : node_hashes)
    {   
        tmp_stream.write((char*)& hash, sizeof(size_t));
    }
   
    tmp_stream.close();
}

void create_temp_hash_file(size_t const ixf_pos, size_t const bin_index, ankerl::unordered_dense::set<size_t> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + "_" + std::to_string(bin_index) + ".tmp";
    auto tmp_dir = std::filesystem::current_path() / "hixf_tmp";
    if (!std::filesystem::exists(tmp_dir))
        std::filesystem::create_directory(tmp_dir);
    auto tmp_file = tmp_dir / ixf_tmp_name;
    std::ofstream tmp_stream(tmp_file,std::ios::out | std::ios::trunc | std::ios::binary);
    for (size_t hash : node_hashes)
    {   
        tmp_stream.write((char*)& hash, sizeof(size_t));
    }
    tmp_stream.close();
}


void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files)
{
    
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_position) + ".tmp";
    auto tmp_dir = std::filesystem::current_path() / "hixf_tmp";
    auto tmp_file = tmp_dir / ixf_tmp_name;
    if (!std::filesystem::exists(tmp_file))
    {
        std::cerr << ixf_tmp_name << "does not exist!" << std::endl << std::flush;
        return;
    }
    std::ifstream tmp_stream(tmp_file, std::ios::in | std::ios::binary);
    size_t x = 0;
    while (tmp_stream.read((char * ) & x, sizeof(size_t)))
    {
        node_hashes.emplace_back(x);
    }
    
    tmp_stream.close();
    tmp_files.insert(tmp_file.string());
    
    //std::filesystem::remove(tmp_file);
    

}

void read_from_temp_hash_file(size_t const ixf_position,
                              uint16_t const bin_index,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files)
{
    
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_position) + "_" + std::to_string(bin_index) + ".tmp";
    auto tmp_dir = std::filesystem::current_path() / "hixf_tmp";
    auto tmp_file = tmp_dir / ixf_tmp_name;
    if (!std::filesystem::exists(tmp_file))
    {
        //std::cerr << ixf_tmp_name << " does not exist!" << std::endl << std::flush;
        return;
    }
    std::ifstream tmp_stream(tmp_file, std::ios::in | std::ios::binary);
    size_t x = 0;
    while (tmp_stream.read((char * ) & x, sizeof(size_t)))
    {
        node_hashes.emplace_back(x);
    }
    tmp_stream.close();
    tmp_files.insert(tmp_file.string());
    

}

}