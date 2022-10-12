
#include "search_hixf.hpp"
#include "search_single.hpp"

namespace hixf
{

//template <bool compressed>
void search_hixf(search_arguments const & arguments)
{
    //using index_structure_t = std::conditional_t<compressed, index_structure::hibf_compressed, index_structure::hixf>;
    auto index = raptor_index<index_structure::hixf>{};
    search_single(arguments, std::move(index));
}

/*
template void search_hibf<false>(search_arguments const & arguments);

template void search_hibf<true>(search_arguments const & arguments);
*/
} 
