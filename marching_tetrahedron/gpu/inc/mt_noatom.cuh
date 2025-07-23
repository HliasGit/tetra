#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include "global.h"
#include "struct_gpu.cuh"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

__global__ void set_flags(  dim_t* d_grid,
                            bool *d_flags,
                            size_t size,
                            double threshold,
                            cube_gpu_SoA* d_relevant_cubes,
                            float *time
                        );

__global__ void compact_scan(int *d_idx, int* d_scan, bool *d_flags);