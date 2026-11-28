#include <cstdio>
#include "parameters.hpp"
#include "geometry.hpp"

int main()
{
    array_t old(NODES_PER_ELEM), current(NODES_PER_ELEM);
    int nodes[NODES_PER_ELEM];
    for (int j = 0; j < NODES_PER_ELEM; ++j) nodes[j] = j;
    double first;
    while (std::scanf("%lf",&first) == 1) {
        old[0][0] = first;
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            for (int d = 0; d < NDIMS; ++d) {
                if (!j && !d) continue;
                if (std::scanf("%lf",&first) != 1) return 1;
                old[j][d] = first;
            }
        for (int j = 0; j < NODES_PER_ELEM; ++j)
            for (int d = 0; d < NDIMS; ++d) {
                if (std::scanf("%lf",&first) != 1) return 1;
                current[j][d] = first;
            }
        const SimplexOverlap result = simplex_overlap(old.view_const(nodes),current.view_const(nodes));
        std::printf("%.17g",result.volume);
        for (int i = 0; i < NODES_PER_ELEM; ++i)
            for (int j = 0; j < NODES_PER_ELEM; ++j) std::printf(" %.17g",result.product[i][j]);
        std::puts("");
    }
}
