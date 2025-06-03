#pragma once
#include <stdio.h>
#include "global.h"
#include "struct.h"
#include "triangles.h"
#include "utils.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


#ifdef __cplusplus
    extern "C"{
#endif


typedef struct __align__(64) cube_gpu {
    int idx;
    int x;
    int y;
    int z;
    int one_apex[3];
    int two_apex[3];
    int _padding[6]; // ensure total size is 64 bytes
} cube_gpu;


typedef struct Triangle_GPU{
    TriangleVertex v1;
    TriangleVertex v2;
    TriangleVertex v3;
} Triangle_GPU;

typedef struct cube_vertices_points{
    CubeVertex coord;
    dim_t value;
} cube_vertices_points;

void remove_unnecessary_cubes(  double *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates);

void skip_preprocessing(   double *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates);

void parallel_march_tetra   (Dimensions *dim, dim_t *grid, int *cube_decomposition, dim_t threshold,
                            size_t *triangle_counter, size_t *vertex_counter, int number_relevant_cubes,
                            cube_gpu **d_relevant_cubes, cube_vertices_points **cube_points_coordinates,
                            int* act_val_vec, int* pairs, Triangle_GPU **triangles);

void allocate_d_grid(dim_t **d_grid, dim_t *grid, size_t size);

void print_cuda_error(cudaError_t err, const char* msg);

void print_relevant_points(cube_gpu *d_relevant_cubes, int *number_relevant_cubes);

void print_triangles(   Triangle_GPU *d_relevant_cubes, int *number_relevant_cubes,
                        char *molecule_name, char *molecule_path);


#ifdef __cplusplus
    }
#endif