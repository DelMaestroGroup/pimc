/**
 * @file estimator.hip.h
 * @author Nathan Nichols
 * @date 04.19.2021
 *
 * @brief Estimator GPU kernels using HIP.
 */
#ifndef ESTIMATOR_GPU_HIP_H 
#define ESTIMATOR_GPU_HIP_H

#include "common_gpu.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNELS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// GPU Kernel for reduction using warp (uses appropriate warp for NVIDIA vs AMD devices i. e. "portable wave aware code")
__device__ void warp_reduce(volatile double *sdata, unsigned int thread_idx) {
    if (warpSize == 64) { if (GPU_BLOCK_SIZE >= 128) sdata[thread_idx] += sdata[thread_idx + 64]; }
    if (GPU_BLOCK_SIZE >= 64) sdata[thread_idx] += sdata[thread_idx + 32];
    if (GPU_BLOCK_SIZE >= 32) sdata[thread_idx] += sdata[thread_idx + 16];
    if (GPU_BLOCK_SIZE >= 16) sdata[thread_idx] += sdata[thread_idx + 8];
    if (GPU_BLOCK_SIZE >= 8) sdata[thread_idx] += sdata[thread_idx + 4];
    if (GPU_BLOCK_SIZE >= 4) sdata[thread_idx] += sdata[thread_idx + 2];
    if (GPU_BLOCK_SIZE >= 2) sdata[thread_idx] += sdata[thread_idx + 1];
}

__device__ void gpu_reduce(volatile double *sdata, unsigned int thread_idx) {
    if (GPU_BLOCK_SIZE >= 1024) {
        if (thread_idx < 512) {
            sdata[thread_idx] += sdata[thread_idx + 512];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (thread_idx < 256) {
            sdata[thread_idx] += sdata[thread_idx + 256];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (thread_idx < 128) {
            sdata[thread_idx] += sdata[thread_idx + 128];
        }
        __syncthreads();
    } 

    if (warpSize == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (thread_idx < 64) {
                sdata[thread_idx] += sdata[thread_idx + 64];
            }
            __syncthreads();
        } 
    }

    if (thread_idx < warpSize) {
        warp_reduce(sdata, thread_idx);
    }
}

// GPU Kernel for ISF calculation
// M timeslices, N particles, N_extent = N + padding
template <bool ZERO_FREQUENCY = false>
__global__ void gpu_isf(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu


    int NNM = N*N*M;
    if (threadIdx.x < NNM) {
        int bead_idx1 = threadIdx.x/N;     // Bead index for first bead
        int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
        int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

        int m_idx2 = (m_idx1 + blockIdx.x) % M;   // Imaginary time index for second bead
        int n_idx2 = (threadIdx.x - bead_idx1*N); // Particle index for second bead
        int bead_idx2 = m_idx2*N + n_idx2;        // Bead index for second bead

        //Get true bead indices in padded beads array
        int true_bead_idx1 = m_idx1*N_extent + n_idx1;
        int true_bead_idx2 = m_idx2*N_extent + n_idx2;

        double q_dot_sep = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            q_dot_sep += qvecs[k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
        }

        s_isf[threadIdx.x] = cos(q_dot_sep);
    } else {
        s_isf[threadIdx.x] = 0.0;
    }

    for (int local_idx = threadIdx.x + GPU_BLOCK_SIZE; local_idx < NNM + GPU_BLOCK_SIZE; local_idx += GPU_BLOCK_SIZE) {
        if (local_idx < NNM) {
            int bead_idx1 = local_idx/N;       // Bead index for first bead
            int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
            int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

            int m_idx2 = (m_idx1 + blockIdx.x) % M; // Imaginary time index for second bead
            int n_idx2 = (local_idx - bead_idx1*N); // Particle index for second bead
            int bead_idx2 = m_idx2*N + n_idx2;      // Bead index for second bead

            //Get true bead indices in padded beads array
            int true_bead_idx1 = m_idx1*N_extent + n_idx1;
            int true_bead_idx2 = m_idx2*N_extent + n_idx2;

            double q_dot_sep = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                q_dot_sep += qvecs[k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
            }

            s_isf[threadIdx.x] += cos(q_dot_sep);
        }
    }
    __syncthreads();
    
    //FIXME This can be abstracted
    // NEED TO REDUCE isf ON SHARED MEMORY AND ADD TO GLOBAL isf
    if (GPU_BLOCK_SIZE >= 1024) {
        if (threadIdx.x < 512) {
            s_isf[threadIdx.x] += s_isf[threadIdx.x + 512];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (threadIdx.x < 256) {
            s_isf[threadIdx.x] += s_isf[threadIdx.x + 256];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (threadIdx.x < 128) {
            s_isf[threadIdx.x] += s_isf[threadIdx.x + 128];
        }
        __syncthreads();
    } 

    if (warpSize == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (threadIdx.x < 64) {
                s_isf[threadIdx.x] += s_isf[threadIdx.x + 64];
            }
            __syncthreads();
        } 
    }

    if (threadIdx.x < warpSize) {
        warp_reduce(s_isf, threadIdx.x);
    }

    if (threadIdx.x == 0) {
        if (ZERO_FREQUENCY) {
            //FIXME would like to get rid of this atomic operation
            atomicAdd(&isf[0], 2.0*s_isf[0]*inorm);
        } else {
            isf[blockIdx.x] = 2.0*s_isf[0]*inorm;
        }
    }
}

// GPU Kernel Launchers
void gpu_isf_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf, dim3(M/2 + 1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}
void gpu_isf_launcher(hipStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf, dim3(M/2 + 1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}

void gpu_ssf_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf, dim3(1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}
void gpu_ssf_launcher(hipStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf, dim3(1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}

void gpu_es_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf<true>, dim3(M/2 + 1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}
void gpu_es_launcher(hipStream_t s, double* __restrict__ isf, double* __restrict__ qvecs, double *beads, double inorm, int M, int N, int N_extent) {
    hipLaunchKernelGGL(gpu_isf<true>, dim3(M/2 + 1), dim3(GPU_BLOCK_SIZE), 0, 0,
            isf, qvecs, beads, inorm, M, N, N_extent);
}

#endif
