#ifndef DYNEARTHSOL3D_SORTINDEX_HPP
#define DYNEARTHSOL3D_SORTINDEX_HPP

#include <algorithm>
#include <numeric>
#include <vector>

/* sort x in ascending order, x is not modified, the sorted order is stored in idx,
 * which the caller sizes to x */
template< typename T, typename index_t >
void sortindex(const std::vector<T>& x, std::vector<index_t>& idx)
{
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&x](index_t left, index_t right) { return x[left] < x[right]; });
}


/* sort x in descending order, as sortindex */
template< typename T, typename index_t >
void sortindex_reversed(const std::vector<T>& x, std::vector<index_t>& idx)
{
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&x](index_t left, index_t right) { return x[left] > x[right]; });
}


#endif
