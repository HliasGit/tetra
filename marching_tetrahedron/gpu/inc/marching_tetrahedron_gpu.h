#pragma once
#include <stdio.h>
#include "global.h"
#include "struct.h"
#include "triangles.h"
#include "utils.h"


#ifdef __cplusplus
    extern "C"{
#endif

typedef struct cube_gpu{
    int idx;
    int x;
    int y;
    int z;
    int one_apex[3];
    int two_apex[3];
} cube_gpu;

typedef struct cube_vertices_points{
    CubeVertex coord;
    dim_t value;
} cube_vertices_points;

void remove_unnecessary_cubes(  double *d_grid, size_t total_size, double threshold,
                                Dimensions *dim, int *relevant_cubes_size,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates);

void parallel_march_tetra   (Dimensions *dim, dim_t *grid, int *cube_decomposition, dim_t threshold,
                            void (*func_ptr)(TriangleVertex *, CubeVertex *, CubeVertex *, dim_t *, dim_t *, dim_t),
                            Polyhedra *p, size_t *triangle_counter, size_t *vertex_counter, int relevant_size,
                            cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates);

void allocate_d_grid(dim_t **d_grid, dim_t *grid, size_t size);


#ifdef __cplusplus
    }
#endif