// Smoke test for the pruned NVHPC image in Dockerfile.cuda.
//
// The image build compiles and links this with the same nvc++ flags the CI jobs
// use, so that a prune which removed something the real build needs fails while
// the image is being built instead of after it has been pushed to :latest.

#include <cstdio>
#include <vector>

#ifdef NPROF
#include <nvtx3/nvToolsExt.h>
#endif

int main()
{
    std::vector<double> v(8, 1.0);
    double *p = v.data();
    double sum = 0.0;

#ifdef NPROF
    nvtxRangePushA("acc_check");
#endif

    #pragma acc parallel loop reduction(+:sum) copyin(p[0:8])
    for (int i = 0; i < 8; ++i)
        sum += p[i];

#ifdef NPROF
    nvtxRangePop();
#endif

    std::printf("%f\n", sum);
    return 0;
}
