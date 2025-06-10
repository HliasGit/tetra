#pragma once
#include <stdio.h>
#include "global.h"
#include "utils.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "struct_gpu.cuh"


void remove_unnecessary_cubes(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates,
                                double *time);

void remove_unnecessary_cubes(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu_SoA **d_relevant_cubes,
                                double *time);

void skip_preprocessing(size_t cubes_in_domain,
                        Dimensions *dim, cube_gpu **d_relevant_cubes,
                        double *time);

void parallel_march_tetra   (   Dimensions *dim, dim_t *grid, int *cube_decomposition, dim_t threshold,
                                int number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **cube_points_coordinates,
                                int* act_val_vec, int* pairs, Triangle_GPU **triangles, int *total_triangles,
                                double *time);

void parallel_march_tetra   (Dimensions *dim, dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            int number_relevant_cubes,
                            cube_gpu_SoA *d_relevant_cubes, cube_vertices_points_SoA **d_cube_points_coordinates,
                            int* act_val_vec, int *pairs, Triangle_GPU **triangles, int *total_triangles,
                            double *time);

void allocate_d_grid(dim_t **d_grid, dim_t *grid, size_t size);

void print_cuda_error(cudaError_t err, const char* msg);

void print_relevant_points(cube_gpu *d_relevant_cubes, int *number_relevant_cubes);

void print_triangles(   Triangle_GPU *d_relevant_cubes, int *number_relevant_cubes,
                        char *molecule_name, char *molecule_path);