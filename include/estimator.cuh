#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace cuda_wrapper {
    void gpu_isf_wrapper(dim3,dim3,double*,double*,double*,double,double,double,int,int,int,int,int,int,int,int);
    void gpu_isf2_wrapper(dim3,dim3,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf2_wrapper(dim3,dim3,cudaStream_t,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf3_wrapper(dim3,dim3,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf3_wrapper(dim3,dim3,cudaStream_t,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf_dse_mean_wrapper(dim3,dim3,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf_dse_mean_wrapper(dim3,dim3,cudaStream_t,double*,double*,double*,int,int,double,int,int,int,int,int,int,int,int,int);
    void gpu_isf_symmetry_reduce_wrapper(dim3,dim3,double*,int,int,int);
    void gpu_isf_set_separations_wrapper(dim3,dim3,double*,double*,int,int,int,int,int,int);
    void gpu_isf_set_separations_wrapper(dim3,dim3,cudaStream_t,double*,double*,int,int,int,int,int,int);
    void gpu_isf_set_isf_wrapper(dim3,dim3,double*,double*,double*,int,int,int,int);
    void gpu_isf_set_isf_wrapper(dim3,dim3,cudaStream_t,double*,double*,double*,int,int,int,int);
}
