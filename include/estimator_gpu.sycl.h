/**
 * @file estimator.h
 * @author Nathan Nichols
 * @date 04.19.2021
 *
 * @brief Estimator GPU kernels using SYCL.
 */

#ifndef ESTIMATOR_GPU_SYCL_H 
#define ESTIMATOR_GPU_SYCL_H

#include "common_gpu.h"


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNELS ---------------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// Atomic add function
template <typename T1, typename T2>
static inline T1
atomicAdd(T1* addr, const T2 val) {
    sycl::atomic_ref<T1,
        sycl::memory_order::relaxed,
        sycl::memory_scope::device,
        sycl::access::address_space::global_space> atom(*addr);
    return atom.fetch_add(static_cast<T1>(val));
}

// GPU Kernel for reduction using warp (uses appropriate warp for NVIDIA vs AMD devices i. e. "portable wave aware code")
void warp_reduce(volatile double *sdata, unsigned int thread_idx) {
    if (SUB_GROUP_SIZE == 64) { if (GPU_BLOCK_SIZE >= 128) sdata[thread_idx] += sdata[thread_idx + 64]; }
    if (GPU_BLOCK_SIZE >= 64) sdata[thread_idx] += sdata[thread_idx + 32];
    if (GPU_BLOCK_SIZE >= 32) sdata[thread_idx] += sdata[thread_idx + 16];
    if (GPU_BLOCK_SIZE >= 16) sdata[thread_idx] += sdata[thread_idx + 8];
    if (GPU_BLOCK_SIZE >= 8) sdata[thread_idx] += sdata[thread_idx + 4];
    if (GPU_BLOCK_SIZE >= 4) sdata[thread_idx] += sdata[thread_idx + 2];
    if (GPU_BLOCK_SIZE >= 2) sdata[thread_idx] += sdata[thread_idx + 1];
}

void gpu_reduce(volatile double *sdata, unsigned int thread_idx) {
    auto grp = sycl::ext::oneapi::this_work_item::get_work_group<1>();
    if (GPU_BLOCK_SIZE >= 1024) {
        if (thread_idx < 512) {
            sdata[thread_idx] += sdata[thread_idx + 512];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (thread_idx < 256) {
            sdata[thread_idx] += sdata[thread_idx + 256];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (thread_idx < 128) {
            sdata[thread_idx] += sdata[thread_idx + 128];
        }
        sycl::group_barrier(grp);
    } 

    if (SUB_GROUP_SIZE == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (thread_idx < 64) {
                sdata[thread_idx] += sdata[thread_idx + 64];
            }
            sycl::group_barrier(grp);
        } 
    }

    if (thread_idx < SUB_GROUP_SIZE) {
        warp_reduce(sdata, thread_idx);
    }
}

// GPU Kernel for ISF calculation
// M timeslices, N particles, N_extent = N + padding
template <bool ZERO_FREQUENCY = false>
void gpu_isf(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    //Get item and group
    auto ndItem = sycl::ext::oneapi::this_work_item::get_nd_item<1>();
    auto grp = sycl::ext::oneapi::this_work_item::get_work_group<1>();

    //Set thread index (work item) and block index (group)
    int threadIdx_x = ndItem.get_local_id(0);
    int blockIdx_x = ndItem.get_group(0);

    //Allocate shared memory
    auto s_isf = *sycl::ext::oneapi::group_local_memory_for_overwrite<double[GPU_BLOCK_SIZE]>(grp);


    int NNM = N*N*M;
    if (threadIdx_x < NNM) {
        int bead_idx1 = threadIdx_x/N;     // Bead index for first bead
        int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
        int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

        int m_idx2 = (m_idx1 + blockIdx_x) % M;   // Imaginary time index for second bead
        int n_idx2 = (threadIdx_x - bead_idx1*N); // Particle index for second bead
        //int bead_idx2 = m_idx2*N + n_idx2;        // Bead index for second bead

        //Get true bead indices in padded beads array
        int true_bead_idx1 = m_idx1*N_extent + n_idx1;
        int true_bead_idx2 = m_idx2*N_extent + n_idx2;

        double q_dot_sep = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            q_dot_sep += qvecs[k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
        }

        s_isf[threadIdx_x] = cos(q_dot_sep);
    } else {
        s_isf[threadIdx_x] = 0.0;
    }

    for (int local_idx = threadIdx_x + GPU_BLOCK_SIZE; local_idx < NNM + GPU_BLOCK_SIZE; local_idx += GPU_BLOCK_SIZE) {
        if (local_idx < NNM) {
            int bead_idx1 = local_idx/N;       // Bead index for first bead
            int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
            int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

            int m_idx2 = (m_idx1 + blockIdx_x) % M; // Imaginary time index for second bead
            int n_idx2 = (local_idx - bead_idx1*N); // Particle index for second bead
            //int bead_idx2 = m_idx2*N + n_idx2;      // Bead index for second bead

            //Get true bead indices in padded beads array
            int true_bead_idx1 = m_idx1*N_extent + n_idx1;
            int true_bead_idx2 = m_idx2*N_extent + n_idx2;

            double q_dot_sep = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                q_dot_sep += qvecs[k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
            }

            s_isf[threadIdx_x] += cos(q_dot_sep);
        }
    }
    sycl::group_barrier(grp);
    
    //FIXME This can be abstracted
    // NEED TO REDUCE isf ON SHARED MEMORY AND ADD TO GLOBAL isf
    if (GPU_BLOCK_SIZE >= 1024) {
        if (threadIdx_x < 512) {
            s_isf[threadIdx_x] += s_isf[threadIdx_x + 512];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (threadIdx_x < 256) {
            s_isf[threadIdx_x] += s_isf[threadIdx_x + 256];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (threadIdx_x < 128) {
            s_isf[threadIdx_x] += s_isf[threadIdx_x + 128];
        }
        sycl::group_barrier(grp);
    } 

    if (SUB_GROUP_SIZE == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (threadIdx_x < 64) {
                s_isf[threadIdx_x] += s_isf[threadIdx_x + 64];
            }
            sycl::group_barrier(grp);
        } 
    }

    if (threadIdx_x < SUB_GROUP_SIZE) {
        warp_reduce(s_isf, threadIdx_x);
    }

    if (threadIdx_x == 0) {
        if (ZERO_FREQUENCY) {
            //FIXME would like to get rid of this atomic operation
            atomicAdd(&isf[0], 2.0*s_isf[0]*inorm);
        } else {
            isf[blockIdx_x] = 2.0*s_isf[0]*inorm;
        }
    }
}

// GPU Kernel for SSF calculation for cylinders
// M timeslices, N particles, N_extent = N + padding, ||bead|| < maxR
void gpu_ssf_cyl(double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, double maxR, int M, int N, int N_extent) {
    //Get item and group
    auto ndItem = sycl::ext::oneapi::this_work_item::get_nd_item<1>();
    auto grp = sycl::ext::oneapi::this_work_item::get_work_group<1>();

    //Set thread index (work item) and block index (group)
    int threadIdx_x = ndItem.get_local_id(0);
    int blockIdx_x = ndItem.get_group(0);

    //Allocate shared memory
    auto s_ssf = *sycl::ext::oneapi::group_local_memory_for_overwrite<double[GPU_BLOCK_SIZE]>(grp);

    int NNM = N*N*M;
    if (threadIdx_x < NNM) {
        int bead_idx1 = threadIdx_x/N;     // Bead index for first bead
        int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
        int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

        int m_idx2 = m_idx1;                      // Imaginary time index for second bead
        int n_idx2 = (threadIdx_x - bead_idx1*N); // Particle index for second bead
        //int bead_idx2 = m_idx2*N + n_idx2;        // Bead index for second bead

        //Get true bead indices in padded beads array
        int true_bead_idx1 = m_idx1*N_extent + n_idx1;
        int true_bead_idx2 = m_idx2*N_extent + n_idx2;

        double q_dot_sep = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            double _bead1 = beads[true_bead_idx1*NDIM + k];
            double _bead2 = beads[true_bead_idx2*NDIM + k];
            q_dot_sep += qvecs[NDIM*blockIdx_x + k]*(_bead2 - _bead1);
        }

        double mag_bead1 = 0.0;
        double mag_bead2 = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM - 1; k++) {
            double _bead1 = beads[true_bead_idx1*NDIM + k];
            double _bead2 = beads[true_bead_idx2*NDIM + k];
            mag_bead1 += _bead1*_bead1;
            mag_bead2 += _bead2*_bead2;
        }

        s_ssf[threadIdx_x] = ((mag_bead1 > maxR*maxR) || (mag_bead2 > maxR*maxR)) ? 0.0 : cos(q_dot_sep);
    } else {
        s_ssf[threadIdx_x] = 0.0;
    }

    for (int local_idx = threadIdx_x + GPU_BLOCK_SIZE; local_idx < NNM + GPU_BLOCK_SIZE; local_idx += GPU_BLOCK_SIZE) {
        if (local_idx < NNM) {
            int bead_idx1 = local_idx/N;       // Bead index for first bead
            int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
            int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

            int m_idx2 = m_idx1;                    // Imaginary time index for second bead
            int n_idx2 = (local_idx - bead_idx1*N); // Particle index for second bead
            //int bead_idx2 = m_idx2*N + n_idx2;      // Bead index for second bead

            //Get true bead indices in padded beads array
            int true_bead_idx1 = m_idx1*N_extent + n_idx1;
            int true_bead_idx2 = m_idx2*N_extent + n_idx2;

            double q_dot_sep = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                double _bead1 = beads[true_bead_idx1*NDIM + k];
                double _bead2 = beads[true_bead_idx2*NDIM + k];
                q_dot_sep += qvecs[NDIM*blockIdx_x + k]*(_bead2 - _bead1);
            }

            double mag_bead1 = 0.0;
            double mag_bead2 = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM - 1; k++) {
                double _bead1 = beads[true_bead_idx1*NDIM + k];
                double _bead2 = beads[true_bead_idx2*NDIM + k];
                mag_bead1 += _bead1*_bead1;
                mag_bead2 += _bead2*_bead2;
            }

            s_ssf[threadIdx_x] = ((mag_bead1 > maxR*maxR) || (mag_bead2 > maxR*maxR)) ? 0.0 : cos(q_dot_sep);
        }
    }
    sycl::group_barrier(grp);
    
    //FIXME This can be abstracted
    // NEED TO REDUCE ssf ON SHARED MEMORY AND ADD TO GLOBAL ssf
    if (GPU_BLOCK_SIZE >= 1024) {
        if (threadIdx_x < 512) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 512];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (threadIdx_x < 256) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 256];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (threadIdx_x < 128) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 128];
        }
        sycl::group_barrier(grp);
    } 

    if (SUB_GROUP_SIZE == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (threadIdx_x < 64) {
                s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 64];
            }
            sycl::group_barrier(grp);
        } 
    }

    if (threadIdx_x < SUB_GROUP_SIZE) {
        warp_reduce(s_ssf, threadIdx_x);
    }

    if (threadIdx_x == 0) {
        ssf[blockIdx_x] = 2.0*s_ssf[0]*inorm;
    }
}

// GPU Kernel for SSF calculation
void gpu_ssf(double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    //Get item and group
    auto ndItem = sycl::ext::oneapi::this_work_item::get_nd_item<1>();
    auto grp = sycl::ext::oneapi::this_work_item::get_work_group<1>();

    //Set thread index (work item) and block index (group)
    int threadIdx_x = ndItem.get_local_id(0);
    int blockIdx_x = ndItem.get_group(0);

    //Allocate shared memory
    auto s_ssf = *sycl::ext::oneapi::group_local_memory_for_overwrite<double[GPU_BLOCK_SIZE]>(grp);

    // Access it like a normal array via reference or pointer:
    //auto s_ssf_ptr = sycl::ext::oneapi::group_local_memory_for_overwrite<double[GPU_BLOCK_SIZE]>(grp);
    //auto& s_ssf = *s_ssf_ptr;  // localMem is double[GPU_BLOCK_SIZE]&

    int NNM = N*N*M;
    if (threadIdx_x < NNM) {
        int bead_idx1 = threadIdx_x/N;     // Bead index for first bead
        int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
        int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

        int m_idx2 = m_idx1;                      // Imaginary time index for second bead
        int n_idx2 = (threadIdx_x - bead_idx1*N); // Particle index for second bead
        //int bead_idx2 = m_idx2*N + n_idx2;        // Bead index for second bead

        //Get true bead indices in padded beads array
        int true_bead_idx1 = m_idx1*N_extent + n_idx1;
        int true_bead_idx2 = m_idx2*N_extent + n_idx2;

        double q_dot_sep = 0.0;
        #pragma unroll
        for (int k = 0; k < NDIM; k++) {
            q_dot_sep += qvecs[NDIM*blockIdx_x + k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
        }

        s_ssf[threadIdx_x] = cos(q_dot_sep);
    } else {
        s_ssf[threadIdx_x] = 0.0;
    }

    for (int local_idx = threadIdx_x + GPU_BLOCK_SIZE; local_idx < NNM + GPU_BLOCK_SIZE; local_idx += GPU_BLOCK_SIZE) {
        if (local_idx < NNM) {
            int bead_idx1 = local_idx/N;       // Bead index for first bead
            int m_idx1 = bead_idx1/N;          // Imaginary time index for first bead
            int n_idx1 = bead_idx1 - m_idx1*N; // Particle index for first bead

            int m_idx2 = m_idx1;                    // Imaginary time index for second bead
            int n_idx2 = (local_idx - bead_idx1*N); // Particle index for second bead
            //int bead_idx2 = m_idx2*N + n_idx2;      // Bead index for second bead

            //Get true bead indices in padded beads array
            int true_bead_idx1 = m_idx1*N_extent + n_idx1;
            int true_bead_idx2 = m_idx2*N_extent + n_idx2;

            double q_dot_sep = 0.0;
            #pragma unroll
            for (int k = 0; k < NDIM; k++) {
                q_dot_sep += qvecs[NDIM*blockIdx_x + k]*(beads[true_bead_idx2*NDIM + k] - beads[true_bead_idx1*NDIM + k]);
            }

            s_ssf[threadIdx_x] += cos(q_dot_sep);
        }
    }
    sycl::group_barrier(grp);
    
    //FIXME This can be abstracted
    // NEED TO REDUCE ssf ON SHARED MEMORY AND ADD TO GLOBAL ssf
    if (GPU_BLOCK_SIZE >= 1024) {
        if (threadIdx_x < 512) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 512];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 512) {
        if (threadIdx_x < 256) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 256];
        }
        sycl::group_barrier(grp);
    } 

    if (GPU_BLOCK_SIZE >= 256) {
        if (threadIdx_x < 128) {
            s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 128];
        }
        sycl::group_barrier(grp);
    } 

    if (SUB_GROUP_SIZE == 32) {
        if (GPU_BLOCK_SIZE >= 128) {
            if (threadIdx_x < 64) {
                s_ssf[threadIdx_x] += s_ssf[threadIdx_x + 64];
            }
            sycl::group_barrier(grp);
        } 
    }

    if (threadIdx_x < SUB_GROUP_SIZE) {
        warp_reduce(s_ssf, threadIdx_x);
    }

    if (threadIdx_x == 0) {
        ssf[blockIdx_x] = 2.0*s_ssf[0]*inorm;
    }
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// GPU KERNEL WRAPPER --------------------------------------------------------
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void gpu_isf_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    sycl::queue q = sycl::queue();
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>((M/2 + 1)*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_isf(isf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}
void gpu_isf_launcher(sycl::queue q, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>((M/2 + 1)*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_isf(isf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}

void gpu_ssf_launcher(double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent, int n_qvecs) {
    sycl::queue q = sycl::queue();
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>(n_qvecs*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_ssf(ssf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}
void gpu_ssf_launcher(sycl::queue q, double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent, int n_qvecs) {
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>(n_qvecs*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_ssf(ssf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}

void gpu_ssf_cyl_launcher(double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, double maxR, int M, int N, int N_extent, int n_qvecs) {
    sycl::queue q = sycl::queue();
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>(n_qvecs*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_ssf_cyl(ssf, qvecs, beads, inorm, maxR, M, N, N_extent);
        });
    });
}
void gpu_ssf_cyl_launcher(sycl::queue q, double* __restrict__ ssf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, double maxR, int M, int N, int N_extent, int n_qvecs) {
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>(n_qvecs*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_ssf_cyl(ssf, qvecs, beads, inorm, maxR, M, N, N_extent);
        });
    });
}

void gpu_es_launcher(double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    sycl::queue q = sycl::queue();
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>((M/2 + 1)*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_isf<true>(isf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}
void gpu_es_launcher(sycl::queue q, double* __restrict__ isf, double* __restrict__ qvecs, double* __restrict__ beads, double inorm, int M, int N, int N_extent) {
    q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(sycl::nd_range<1>(sycl::range<1>((M/2 + 1)*GPU_BLOCK_SIZE), sycl::range<1>(GPU_BLOCK_SIZE)),
                [=](sycl::nd_item<1> item) [[sycl::reqd_sub_group_size(SUB_GROUP_SIZE)]] {
            gpu_isf<true>(isf, qvecs, beads, inorm, M, N, N_extent);
        });
    });
}

#endif
