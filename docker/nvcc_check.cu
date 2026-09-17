// Smoke test for the pruned NVHPC image in Dockerfile.cuda -- nvcc half.
//
// The knn-bvh submodule the openacc=1 jobs build compiles its kernels with the
// SDK's bundled nvcc rather than nvc++, and resolves CUDA_HOME by looking up
// `which nvcc`. Both break quietly if cuda/bin is not on PATH, so the image
// build compiles this the same way knn-bvh does.

#include <cstdio>

__global__ void axpy(int n, double a, const double *x, double *y)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n)
        y[i] += a * x[i];
}

extern "C" void launch_axpy(int n, double a, const double *x, double *y)
{
    axpy<<<(n + 255) / 256, 256>>>(n, a, x, y);
}
