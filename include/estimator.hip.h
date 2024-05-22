/**
 * @file estimator.h
 * @author Nathan Nichols
 * @date 04.19.2021
 *
 * @brief Estimator GPU kernels using HIP.
 */

#include "common_gpu.h"

#ifndef ESTIMATOR_HIP_H 
#define ESTIMATOR_HIP_H

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNELS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#ifdef GPU_BLOCK_SIZE
// GPU Kernel for ISF calculation

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
__global__
void gpu_isf(double *isf, double *qvecs, double *beads, int qvec_idx, int tau_idx, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
    int i, bead1_idx, bead2_idx, true_bead1_idx, true_bead2_idx, tau_idx_first_bead,
        particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    double q_dot_sep;
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu
    i = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
    if (i < NNM) {
        bead1_idx = int(i/number_of_particles);
        tau_idx_first_bead = int(bead1_idx/number_of_particles);
        tau_idx_second_bead = (tau_idx_first_bead + tau_idx) % number_of_timeslices;
        particle_idx_second_bead = (i - bead1_idx*number_of_particles);
        bead2_idx = tau_idx_second_bead*number_of_particles + particle_idx_second_bead;
        if (bead1_idx != bead2_idx) {
            particle_idx_first_bead = bead1_idx - tau_idx_first_bead*number_of_particles;

            true_bead1_idx = tau_idx_first_bead*full_number_of_particles + particle_idx_first_bead;
            true_bead2_idx = tau_idx_second_bead*full_number_of_particles + particle_idx_second_bead;

            q_dot_sep = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                q_dot_sep += qvecs[qvec_idx*NDIM + k]*(beads[true_bead2_idx*NDIM + k] - beads[true_bead1_idx*NDIM + k]);
            }

            s_isf[hipThreadIdx_x] = cos(q_dot_sep)*inorm;
        } else {
            s_isf[hipThreadIdx_x] = inorm;
        }
    } else {
        s_isf[hipThreadIdx_x] = 0.0;
    }
    __syncthreads();

    // NEED TO REDUCE isf ON SHARED MEMORY AND ADD TO GLOBAL isf
    if (GPU_BLOCK_SIZE >= 1024) {
        if (hipThreadIdx_x < 512) {
            s_isf[hipThreadIdx_x] += s_isf[hipThreadIdx_x + 512];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (hipThreadIdx_x < 256) {
            s_isf[hipThreadIdx_x] += s_isf[hipThreadIdx_x + 256];
        }
        __syncthreads();
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (hipThreadIdx_x < 128) {
            s_isf[hipThreadIdx_x] += s_isf[hipThreadIdx_x + 128];
        }
        __syncthreads();
    } 

    if (warpSize == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (hipThreadIdx_x < 64) {
                s_isf[hipThreadIdx_x] += s_isf[hipThreadIdx_x + 64];
            }
            __syncthreads();
        } 
    }

    if (hipThreadIdx_x < warpSize) {
        warp_reduce(s_isf, hipThreadIdx_x);
    }

    if (hipThreadIdx_x == 0) {
        //NOTE: May see some performance gain here if temporarily store results
        // to some global device variable and then launch a separate kernel to
        // again reduce those results i.e.
        // tmp_isf[hipBlockIdx_x] = s_isf[0];
        // ^-- reduce on this, but this may get too bloated with multiple qvecs
        // and multiple kernel lanuches per tau_idx (imaginary time separation)
        if ((tau_idx == 0) || (tau_idx == beta_over_two_idx)) {
            atomicAdd(&isf[tau_idx*number_of_qvecs + qvec_idx], 2.0*s_isf[0]);
        } else {
            int _tau_idx = tau_idx > beta_over_two_idx ? number_of_timeslices - tau_idx : tau_idx; 
            atomicAdd(&isf[_tau_idx*number_of_qvecs + qvec_idx], s_isf[0]);
        }
    }
}

#endif
#endif
