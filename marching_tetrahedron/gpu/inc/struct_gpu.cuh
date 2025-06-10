#pragma once

#include <cuda_runtime.h>

typedef struct cube_gpu_SoA{
    float4 *coord_idx;
    float4 *one_apex;
    float4 *two_apex;
} cube_gpu_SoA;

typedef struct cube_vertices_points_SoA{
    float4 val;
} cube_vertices_points_SoA;