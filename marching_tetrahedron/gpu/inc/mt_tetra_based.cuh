#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include "struct_gpu.cuh"
#include "global.h"

__global__ void manage_tetra(int n_total_tetra, dim_t threshold, dim_t *d_grid, float4 *d_active_tetrahedra, int *d_CD1, int* d_num_active_tetrahedra);


__global__ void make_triangles( int n_active_tetra,
                                int *d_num_generated_triangles,
                                float4 *d_active_tetrahedra,
                                Triangle_GPU *d_triangles,
                                dim_t threshold,
                                int *d_pairs,
                                int *d_7val,
                                int *d_0val,
                                int *d_act_val_vec);
                                
__device__ void make_triangle_from_tetra(float4 *point1, float4 *point2, float4 *point3, float4 *point4,
                                    Triangle_GPU *triangle, int *pairs);