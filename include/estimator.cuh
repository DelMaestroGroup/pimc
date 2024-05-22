#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace cuda_wrapper {
    void gpu_isf_wrapper(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent);
    void gpu_isf_wrapper(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent);

    void gpu_ssf_wrapper(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent);
    void gpu_ssf_wrapper(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent);
}
