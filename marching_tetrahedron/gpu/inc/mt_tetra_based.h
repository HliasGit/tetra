#pragma once

#include "global.h"
#include "struct_gpu.cuh"

void load_to_const_tetra(Dimensions *dimensions);

void remove_unnecessary_tetrahedra( dim_t *d_grid,
                                    size_t cubes_in_domain,
                                    int *CD1,
                                    dim_t threshold,
                                    float4 **d_active_tetrahedra,
                                    int **d_num_active_tetrahedra,
                                    float *time);

void remove_unnecessary_cubes_mixed(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                int *number_relevant_cubes,
                                cube_gpu_SoA **d_relevant_cubes,
                                double *time);

void tetra_mt(  int *d_num_active_tetrahedra,
                float4 *d_active_tetrahedra,
                dim_t threshold,
                int* pairs,
                int *act_val_vec,
                float *time);

void parallel_march_tetra_mixed(dim_t *d_grid,
                                int *cube_decomposition,
                                dim_t threshold,
                                int number_relevant_cubes,
                                cube_gpu_SoA *d_relevant_cubes,
                                float4 *d_active_tetrahedra,
                                int* act_val_vec,
                                int *pairs,
                                Triangle_GPU **triangles,
                                int *total_triangles,
                                double *time);