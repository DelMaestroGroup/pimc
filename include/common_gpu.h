/* We either specify the number of dimensions of the simulation at compile time,
 * or it defaults to 1D. */
#ifndef NDIM
#define NDIM 1 ///< Number of spatial dimnsions
#endif

/* Determine if using GPU acceleration */
#ifdef GPU_BLOCK_SIZE
#ifndef USE_CUDA
#include "hip/hip_runtime.h"
#define HIP_ASSERT(x) (assert((x)==hipSuccess))
#endif
#ifdef USE_CUDA
#include <cuda_runtime.h>
#define CUDA_ASSERT(x) (assert((x)==cudaSuccess))
#endif
#ifndef MAX_GPU_STREAMS
#define MAX_GPU_STREAMS 1 ///< Max number of concurrent GPU streams
#endif
#endif
