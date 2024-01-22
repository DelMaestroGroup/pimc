/**
 * @file estimator.h
 * @author Nathan Nichols
 * @date 04.19.2021
 *
 * @brief Estimator GPU kernels using HIP.
 */

#include "common_gpu.h"
#include "special_functions.cuh"
#include "estimator.cuh"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNELS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// GPU Kernel for ISF calculation
__global__
void gpu_isf(double *isf, double *qvecs, double *beads, double norm, double norm2, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int beta_over_two_idx) {
    int i, tau_idx, tau_idx_first_bead, particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    double q_dot_sep;
    i = blockDim.x * blockIdx.x + threadIdx.x;
    for (int tile = 0; tile < gridDim.x; tile++) {
        __shared__ double s_beads[GPU_BLOCK_SIZE*NDIM]; // beads on shared memory of gpu
        int s_beads_idx = tile * blockDim.x + threadIdx.x;
        if (s_beads_idx < full_number_of_beads) {
            for (int j = 0; j < NDIM; j++) {
                s_beads[threadIdx.x*NDIM + j] = beads[s_beads_idx*NDIM + j];
            }
        }
        __syncthreads();
  
        if (i < full_number_of_beads) {
            tau_idx_first_bead = int(i/full_number_of_particles);
            particle_idx_first_bead = i - tau_idx_first_bead*full_number_of_particles;
            for (int j = 0; j < GPU_BLOCK_SIZE; j++) {
                s_beads_idx = tile * blockDim.x + j;
                if (s_beads_idx < full_number_of_beads) {
                    tau_idx_second_bead = int(s_beads_idx/full_number_of_particles);
                    particle_idx_second_bead = s_beads_idx - tau_idx_second_bead*full_number_of_particles;
                    if ((particle_idx_second_bead < number_of_particles) &&
                            (particle_idx_first_bead < number_of_particles)) {
                        tau_idx = abs(tau_idx_second_bead - tau_idx_first_bead);
                        tau_idx = min(tau_idx, number_of_timeslices - tau_idx);
                        if (i != s_beads_idx) {
                            if ((tau_idx == 0) || (tau_idx == beta_over_two_idx)) {
                                for (int nq = 0; nq < number_of_qvecs; nq++) {
                                    q_dot_sep = 0.0;
                                    for (int k = 0; k < NDIM; k++) {
                                        q_dot_sep += qvecs[nq*NDIM + k]*(s_beads[j*NDIM + k] - beads[i*NDIM + k]);
                                    }
                                    atomicAdd(&isf[tau_idx*number_of_qvecs + nq], cos(q_dot_sep)/norm);
                                }
                            } else {
                                for (int nq = 0; nq < number_of_qvecs; nq++) {
                                    q_dot_sep = 0.0;
                                    for (int k = 0; k < NDIM; k++) {
                                        q_dot_sep += qvecs[nq*NDIM + k]*(s_beads[j*NDIM + k] - beads[i*NDIM + k]);
                                    }
                                    atomicAdd(&isf[tau_idx*number_of_qvecs + nq], cos(q_dot_sep)/norm2);
                                }
                            } 
                        } else {
                            for (int nq = 0; nq < number_of_qvecs; nq++) {
                                atomicAdd(&isf[tau_idx*number_of_qvecs + nq], inorm);
                            }
                        }
                    }
                }
            }
        }
        __syncthreads();
    }
}

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
void gpu_isf2(double *isf, double *qvecs, double *beads, int qvec_idx, int tau_idx, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
    int i, bead1_idx, bead2_idx, true_bead1_idx, true_bead2_idx, tau_idx_first_bead,
        particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    double q_dot_sep;
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu
    i = blockDim.x * blockIdx.x + threadIdx.x;
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

            s_isf[threadIdx.x] = cos(q_dot_sep)*inorm;
        } else {
            s_isf[threadIdx.x] = inorm;
        }
    } else {
        s_isf[threadIdx.x] = 0.0;
    }
    __syncthreads();

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

// GPU Kernel for ISF calculation
__global__
void gpu_isf3(double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
    int i, bead1_idx, bead2_idx, true_bead1_idx, true_bead2_idx, tau_idx_first_bead,
        particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    double q_dot_sep;
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu
    i = blockDim.x * blockIdx.x + threadIdx.x;
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
                q_dot_sep += pow(beads[true_bead2_idx*NDIM + k] - beads[true_bead1_idx*NDIM + k],2);
            }
            q_dot_sep = qmags[qmag_idx]*sqrt(q_dot_sep);

            s_isf[threadIdx.x] = sin(q_dot_sep)*inorm/q_dot_sep;
        } else {
            s_isf[threadIdx.x] = inorm;
        }
    } else {
        s_isf[threadIdx.x] = 0.0;
    }
    __syncthreads();

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
        //NOTE: May see some performance gain here if temporarily store results
        // to some global device variable and then launch a separate kernel to
        // again reduce those results i.e.
        // tmp_isf[hipBlockIdx_x] = s_isf[0];
        // ^-- reduce on this, but this may get too bloated with multiple qmags
        // and multiple kernel lanuches per tau_idx (imaginary time separation)
        if ((tau_idx == 0) || (tau_idx == beta_over_two_idx)) {
            atomicAdd(&isf[tau_idx*number_of_qmags + qmag_idx], 2.0*s_isf[0]);
        } else {
            int _tau_idx = tau_idx > beta_over_two_idx ? number_of_timeslices - tau_idx : tau_idx; 
            atomicAdd(&isf[_tau_idx*number_of_qmags + qmag_idx], s_isf[0]);
        }
    }
}

// GPU Kernel for ISF calculation
__global__
void gpu_isf_dse_mean(double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
    int i, bead1_idx, bead2_idx, true_bead1_idx, true_bead2_idx, tau_idx_first_bead,
        particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    double r, ar, br, _si_ar, _si_br;
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu
    i = blockDim.x * blockIdx.x + threadIdx.x;
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

            r = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                r += pow(beads[true_bead2_idx*NDIM + k] - beads[true_bead1_idx*NDIM + k],2);
            }
            r = sqrt(r);
            ar = qmags[qmag_idx]*r;
            br = qmags[qmag_idx + 1]*r;
            si(ar,&_si_ar);
            si(br,&_si_br);
            s_isf[threadIdx.x] = inorm*(_si_ar - _si_br)/(ar - br);
        } else {
            s_isf[threadIdx.x] = inorm;
        }
    } else {
        s_isf[threadIdx.x] = 0.0;
    }
    __syncthreads();

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
        //NOTE: May see some performance gain here if temporarily store results
        // to some global device variable and then launch a separate kernel to
        // again reduce those results i.e.
        // tmp_isf[hipBlockIdx_x] = s_isf[0];
        // ^-- reduce on this, but this may get too bloated with multiple qmags
        // and multiple kernel lanuches per tau_idx (imaginary time separation)
        if ((tau_idx == 0) || (tau_idx == beta_over_two_idx)) {
            atomicAdd(&isf[tau_idx*number_of_qmags + qmag_idx], 2.0*s_isf[0]);
        } else {
            int _tau_idx = tau_idx > beta_over_two_idx ? number_of_timeslices - tau_idx : tau_idx; 
            atomicAdd(&isf[_tau_idx*number_of_qmags + qmag_idx], s_isf[0]);
        }
    }
}

// GPU Kernel to utilize beta/2 symmetry
__global__
void gpu_isf_symmetry_reduce(double *isf, int number_of_qvecs, int number_of_timeslices,
        int beta_over_two_idx) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int tau_idx = int(i/number_of_qvecs);
    if ((tau_idx > 0) && (tau_idx < beta_over_two_idx)) {
        int q_idx = i - tau_idx*number_of_qvecs;
        int tau_idx2 = number_of_timeslices - tau_idx;
        isf[i] += isf[tau_idx2*number_of_qvecs + q_idx];
    }
}

//GPU Kerenel to set bead separations
__global__
void gpu_isf_set_separations(
        double *separations, double *beads, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads) {
    int i, j, separation_idx, tau_sep, bead1_idx, bead2_idx, true_bead1_idx, true_bead2_idx, tau_idx_first_bead,
        particle_idx_first_bead, tau_idx_second_bead, particle_idx_second_bead;
    int t_idx = blockDim.x * blockIdx.x + threadIdx.x;
    i = t_idx/number_of_beads;
    j = t_idx - bead1_idx*number_of_beads;
    if (j < i) {
        bead1_idx = i;
        bead2_idx = j;
         
        tau_idx_first_bead = int(bead1_idx/number_of_particles);
        tau_idx_second_bead =  int(bead2_idx/number_of_particles);
        tau_sep = abs(tau_idx_second_bead - tau_idx_first_bead);
        
        particle_idx_first_bead = bead1_idx - tau_idx_first_bead*number_of_particles;
        particle_idx_second_bead = bead2_idx - tau_idx_second_bead*number_of_particles;

        true_bead1_idx = tau_idx_first_bead*full_number_of_particles + particle_idx_first_bead;
        true_bead2_idx = tau_idx_second_bead*full_number_of_particles + particle_idx_second_bead;

        if (tau_sep == 0) {
            int ii = particle_idx_first_bead;
            int jj = particle_idx_second_bead;
            if (ii < jj) { ii, jj = jj, ii; }
            separation_idx = (ii*(ii - 1)/2) + jj;
        }
        if ((tau_sep > 0) && (tau_sep < number_of_timeslices/2)) {
            separation_idx = number_of_beads*(number_of_particles - 1)/2;
            separation_idx += number_of_beads*number_of_particles*(tau_sep - 1);
            separation_idx += i*number_of_particles + particle_idx_second_bead;
        }
        if (tau_sep == number_of_timeslices/2) {
            separation_idx = number_of_beads*(number_of_particles - 1)/2;
            separation_idx += number_of_beads*number_of_particles*(number_of_timeslices/2 - 1);
            int ii = i;
            int jj = j;
            int iiN = particle_idx_first_bead;
            int jjN = particle_idx_second_bead;
            if (ii > jj) {
                ii, jj = jj, ii;
                iiN, jjN = jjN, iiN;
            }
            separation_idx += ii*number_of_particles + jjN;
        }

        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            separations[separation_idx*NDIM + k] = beads[true_bead2_idx*NDIM + k] - beads[true_bead1_idx*NDIM + k];
        }
    }
}

__device__
void gpu_isf_set_isf_piece(
        volatile double *isf, double *qvecs, double *separations, int number_of_separations, double tau_sep_factor) {
    for (int i=0; i<number_of_separations; i++) {
        double q_dot_r = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            q_dot_r += separations[i*NDIM + k]*qvecs[k];
        }
        *isf += tau_sep_factor*cos(q_dot_r);
    }
}

__global__
void gpu_isf_set_isf(
        double *isf, double *qvecs, double *separations, int number_of_qvecs,
        int number_of_timeslices, int number_of_particles, int number_of_beads) {
    __shared__ double s_isf[GPU_BLOCK_SIZE]; // temporarily store isf on shared memory of gpu
    int t_idx = blockDim.x * blockIdx.x + threadIdx.x;
    for (int i=0; i<number_of_qvecs; i++) {
        int separation_idx = 0;
        int number_of_separations = number_of_beads*(number_of_particles - 1)/2;
        if (t_idx < number_of_separations) {
            gpu_isf_set_isf_piece(s_isf + threadIdx.x, qvecs + i*NDIM, separations + separation_idx*NDIM, number_of_separations, 2.0/number_of_beads);
            if (t_idx == 0) {
                 s_isf[threadIdx.x] += 1.0;
            }
        } else {
            s_isf[threadIdx.x] = 0.0;
        }
        __syncthreads();
        gpu_reduce(s_isf, threadIdx.x);
        if (threadIdx.x == 0) {
            atomicAdd(&isf[i*(number_of_timeslices/2 + 1)], s_isf[0]);
        }
        for (int j=1; j<number_of_timeslices/2; j++) {
            separation_idx = number_of_beads*(number_of_particles - 1)/2;
            separation_idx += number_of_beads*number_of_particles*(j - 1);
            number_of_separations = number_of_beads*number_of_particles;
            if (t_idx < number_of_separations) {
                gpu_isf_set_isf_piece(s_isf + threadIdx.x, qvecs + i*NDIM, separations + separation_idx*NDIM, number_of_separations, 1.0/number_of_beads);
            } else {
                s_isf[threadIdx.x] = 0.0;
            }
            __syncthreads();
            gpu_reduce(s_isf, threadIdx.x);
            if (threadIdx.x == 0) {
                atomicAdd(&isf[i*(number_of_timeslices/2 + 1) + j], s_isf[0]);
            }
        }
        separation_idx = number_of_beads*(number_of_particles - 1)/2;
        separation_idx += number_of_beads*number_of_particles*(number_of_timeslices/2 - 1);

        number_of_separations = number_of_beads*number_of_particles/2;
        if (t_idx < number_of_separations) {
            gpu_isf_set_isf_piece(s_isf + threadIdx.x, qvecs + i*NDIM, separations + separation_idx*NDIM, number_of_separations, 2.0/number_of_beads);
            if (t_idx == 0) {
                 s_isf[threadIdx.x] += 1.0;
            }
        } else {
            s_isf[threadIdx.x] = 0.0;
        }
        __syncthreads();
        gpu_reduce(s_isf, threadIdx.x);
        if (threadIdx.x == 0) {
            atomicAdd(&isf[i*(number_of_timeslices/2 + 1) + number_of_timeslices/2], s_isf[0]);
        }
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNEL WRAPPER --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
namespace cuda_wrapper {
    void gpu_isf_wrapper(dim3 grid_size, dim3 group_size, double *isf, double *qvecs, double *beads,
        double norm, double norm2, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int beta_over_two_idx) {
        gpu_isf<<<grid_size, group_size, 0, 0>>> ( 
            isf, qvecs, beads, norm, norm2, inorm, number_of_qvecs, number_of_timeslices, number_of_particles,
            number_of_beads, full_number_of_timeslices, full_number_of_particles, full_number_of_beads, beta_over_two_idx);
    }
    void gpu_isf2_wrapper(dim3 grid_size, dim3 group_size,
        double *isf, double *qvecs, double *beads, int qvec_idx, int tau_idx, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf2<<<grid_size, group_size, 0, 0>>> (
            isf, qvecs, beads, qvec_idx, tau_idx, inorm, number_of_qvecs,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf2_wrapper(dim3 grid_size, dim3 group_size, cudaStream_t stream,
        double *isf, double *qvecs, double *beads, int qvec_idx, int tau_idx, double inorm,
        int number_of_qvecs, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf2<<<grid_size, group_size, 0, stream>>> (
            isf, qvecs, beads, qvec_idx, tau_idx, inorm, number_of_qvecs,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf3_wrapper(dim3 grid_size, dim3 group_size,
        double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf3<<<grid_size, group_size, 0, 0>>> (
            isf, qmags, beads, qmag_idx, tau_idx, inorm, number_of_qmags,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf3_wrapper(dim3 grid_size, dim3 group_size, cudaStream_t stream,
        double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf3<<<grid_size, group_size, 0, stream>>> (
            isf, qmags, beads, qmag_idx, tau_idx, inorm, number_of_qmags,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf_dse_mean_wrapper(dim3 grid_size, dim3 group_size,
        double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf_dse_mean<<<grid_size, group_size, 0, 0>>> (
            isf, qmags, beads, qmag_idx, tau_idx, inorm, number_of_qmags,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf_dse_mean_wrapper(dim3 grid_size, dim3 group_size, cudaStream_t stream,
        double *isf, double *qmags, double *beads, int qmag_idx, int tau_idx, double inorm,
        int number_of_qmags, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads, int NNM, int beta_over_two_idx) {
        gpu_isf_dse_mean<<<grid_size, group_size, 0, stream>>> (
            isf, qmags, beads, qmag_idx, tau_idx, inorm, number_of_qmags,
            number_of_timeslices, number_of_particles, number_of_beads,
            full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads, NNM, beta_over_two_idx);
    }
    void gpu_isf_symmetry_reduce_wrapper(dim3 grid_size, dim3 group_size,
        double *isf, int number_of_qvecs, int number_of_timeslices,
        int beta_over_two_idx) {
        gpu_isf_symmetry_reduce<<<grid_size, group_size, 0, 0>>> (
            isf, number_of_qvecs, number_of_timeslices, beta_over_two_idx);
    }
    void gpu_isf_set_separations_wrapper(dim3 grid_size, dim3 group_size,
        double *separations, double *beads, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads) {
        gpu_isf_set_separations<<<grid_size, group_size, 0, 0>>> (
            separations, beads, number_of_timeslices, number_of_particles,
            number_of_beads, full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads);
    }
    void gpu_isf_set_separations_wrapper(dim3 grid_size, dim3 group_size, cudaStream_t stream,
        double *separations, double *beads, int number_of_timeslices, int number_of_particles,
        int number_of_beads, int full_number_of_timeslices, int full_number_of_particles,
        int full_number_of_beads) {
        gpu_isf_set_separations<<<grid_size, group_size, 0, stream>>> (
            separations, beads, number_of_timeslices, number_of_particles,
            number_of_beads, full_number_of_timeslices, full_number_of_particles,
            full_number_of_beads);
    }
    void gpu_isf_set_isf_wrapper(dim3 grid_size, dim3 group_size,
        double *isf, double *qvecs, double *separations, int number_of_qvecs,
        int number_of_timeslices, int number_of_particles, int number_of_beads) {
        gpu_isf_set_isf<<<grid_size, group_size, 0, 0>>> (
            isf, qvecs, separations, number_of_qvecs,
            number_of_timeslices, number_of_particles, number_of_beads);
    }
    void gpu_isf_set_isf_wrapper(dim3 grid_size, dim3 group_size, cudaStream_t stream,
        double *isf, double *qvecs, double *separations, int number_of_qvecs,
        int number_of_timeslices, int number_of_particles, int number_of_beads) {
        gpu_isf_set_isf<<<grid_size, group_size, 0, stream>>> (
            isf, qvecs, separations, number_of_qvecs,
            number_of_timeslices, number_of_particles, number_of_beads);
    }
}
