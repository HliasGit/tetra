#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include "struct_gpu.cuh"
#include "global.h"

__global__ void manage_tetra(int n_tetra, float threshold, float *grid, float4 *d_results, int *d_CD1, int *d_CD2, int* counter);

__global__ void make_triangles( int *d_counter, float4 *d_results, Triangle_GPU *d_triangles, coord_t threshold,
                                int *d_pairs, int* d_7val, int* d_0val);

__device__ void make_triangle_from_tetra(float4 *point1, float4 *point2, float4 *point3, float4 *point4,
                                    Triangle_GPU *triangle, int *pairs);