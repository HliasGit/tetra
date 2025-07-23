#pragma once

#include <cuda_runtime.h>
#include "global.h"
#include "struct.h"
#include "struct_gpu.cuh"



__global__ void remove_unnecessary_cubes_kernel(dim_t *grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu *d_relevant_cubes);

__global__ void remove_unnecessary_cubes_SoA_kernel(dim_t *grid, int *counter,
                                                size_t size, double threshold,
                                                cube_gpu_SoA *d_relevant_cubes);
                                                
__global__ void skip_preprocessing_k(size_t size,
                                    Dimensions *dim, cube_gpu* d_relevant_cubes);

__global__ void compute_apex(   dim_t *grid, cube_gpu *d_relevant_cubes, int *number_relevant_cubes,
                                cube_vertices_points *d_cube_points_coordinates,
                                Dimensions *dim);

__global__ void compute_apex_float4(   dim_t *grid, cube_gpu_SoA *d_relevant_cubes, int number_relevant_cubes,
                                cube_vertices_points_SoA *d_cube_points_coordinates);

__global__ void compute_march_tetra(dim_t *d_grid, cube_gpu *d_relevant_cubes,
                                    int number_relevant_cubes, int *cube_deco,
                                    cube_vertices_points *d_cube_points_coordinates,
                                    cube_vertices_points *memory_pool, int *pool_index,
                                    dim_t threshold, int* act_val_vec, int* d_pairs,
                                    Triangle_GPU *d_triangles, int *d_counter, Dimensions *dim);

__global__ void compute_march_tetra_SoA(dim_t *d_grid, cube_gpu_SoA *d_relevant_cubes,
                                        int number_relevant_cubes, int *cube_deco,
                                        cube_vertices_points_SoA *d_cube_points_coordinates,
                                        cube_vertices_points_SoA *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter);

__device__ void sort_points(cube_vertices_points **first, cube_vertices_points **second,
                            cube_vertices_points **third, cube_vertices_points **fourth);

__device__ void sort_points_SoA(cube_vertices_points_SoA **first, cube_vertices_points_SoA **second,
                                cube_vertices_points_SoA **third, cube_vertices_points_SoA **fourth);

__device__ void count_elements( int *less, int *eq, int *gre, cube_vertices_points *first,
                                cube_vertices_points *second, cube_vertices_points *third,
                                cube_vertices_points *fourth, dim_t threshold);

__device__ void count_elements_SoA( int *less, int *eq, int *gre, cube_vertices_points_SoA *first,
                                cube_vertices_points_SoA *second, cube_vertices_points_SoA *third, cube_vertices_points_SoA *fourth,
                                dim_t threshold);

__device__ int get_action_value( int less, int eq, int gre);

__device__ void make_triangle(  cube_vertices_points *first, cube_vertices_points *second,
                                cube_vertices_points *third, cube_vertices_points *fourth,
                                Triangle_GPU *triangle, int *pairs);

__device__ void make_triangle_SoA(  cube_vertices_points_SoA *first, cube_vertices_points_SoA *second,
                                cube_vertices_points_SoA *third, cube_vertices_points_SoA *fourth,
                                Triangle_GPU *triangle, int *pairs);

// FLAG VERSION

__global__ void remove_unnecessary_cubes_flag(  dim_t* d_grid,
                                                bool *d_flags,
                                                size_t size,
                                                double threshold,
                                                cube_gpu_SoA* d_relevant_cubes
                                            );