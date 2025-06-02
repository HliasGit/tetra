#pragma once

#include <cuda_runtime.h>
#include "global.h"
#include "struct.h"

__global__ void remove_unnecessary_cubes_kernel(double *grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu *d_relevant_cubes);

__global__ void compute_apex(   double *grid, cube_gpu *d_relevant_cubes, int *number_relevant_cubes,
                                cube_vertices_points *d_cube_points_coordinates,
                                Dimensions *dim);

__global__ void compute_march_tetra(double *d_grid, cube_gpu *d_relevant_cubes,
                                    int number_relevant_cubes, int *cube_deco,
                                    cube_vertices_points *d_cube_points_coordinates,
                                    cube_vertices_points *memory_pool, int *pool_index,
                                    dim_t threshold, int* act_val_vec, int* d_pairs,
                                    Triangle_GPU *d_triangles, int *d_counter, Dimensions *dim);


__device__ void sort_points(cube_vertices_points **first, cube_vertices_points **second,
                            cube_vertices_points **third, cube_vertices_points **fourth);


__device__ void count_elements( int *less, int *eq, int *gre, cube_vertices_points *first,
                                cube_vertices_points *second, cube_vertices_points *third,
                                cube_vertices_points *fourth, dim_t threshold);

__device__ int get_action_value( int less, int eq, int gre);

__device__ void make_triangle(  cube_vertices_points *first, cube_vertices_points *second,
                                cube_vertices_points *third, cube_vertices_points *fourth,
                                Triangle_GPU *triangle, int *pairs, bool debug);