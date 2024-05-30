#ifndef ESTIMATOR_GPU_CUH 
#define ESTIMATOR_GPU_CUH
#include "common_gpu.h"
#include "device_launch_parameters.h"

void gpu_isf_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent);
void gpu_isf_launcher(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent);

void gpu_ssf_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent, int n_qvecs);
void gpu_ssf_launcher(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent, int n_qvecs);

void gpu_ssf_cyl_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, double maxR, int M, int N, int N_extent, int n_qvecs);
void gpu_ssf_cyl_launcher(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, double maxR, int M, int N, int N_extent, int n_qvecs);

void gpu_es_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent);
void gpu_es_launcher(cudaStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent);
#endif
