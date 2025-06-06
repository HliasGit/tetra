#ifndef TETRA
#define TETRA

#ifdef VERBOSE
    #define verbose_print(...) printf(__VA_ARGS__)
    #define verbose_call(...) __VA_ARGS__
#else
    #define verbose_print(...)
    #define verbose_call(...)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef __cplusplus
    extern "C"{
#endif

typedef float dim_t;    // type of the scalar field values
typedef float coord_t;  // type for the coordinates (float or double)

typedef struct dimensions{
    size_t x_dim;
    size_t y_dim;
    size_t z_dim;
} Dimensions;

typedef struct cube{
    size_t x;
    size_t y;
    size_t z;
} CubeVertex;

typedef struct vertex{
    coord_t x;
    coord_t y;
    coord_t z;
} TriangleVertex;

typedef struct triangle{
    TriangleVertex *v1;
    TriangleVertex *v2;
    TriangleVertex *v3;
} Triangle;

typedef struct StackNode{
    CubeVertex coordinate;
    dim_t owned_value;
    struct StackNode *next;
    int point;
} StackNode;

typedef struct TriangleCoordNode {
    double coordinate1;
    double coordinate2;
    double coordinate3;
    int level;
    struct TriangleCoordNode* next_list;
    struct TriangleCoordNode* next_level;
    size_t index;
    size_t local_index;
} TriangleCoordNode;

typedef struct TriangleNode{
    TriangleCoordNode *vert1;
    TriangleCoordNode *vert2;
    TriangleCoordNode *vert3;
    int triangle_index;
    struct TriangleNode *next;
} TriangleNode;

typedef struct Polyhedra{
    struct TriangleNode *triangles;
    TriangleCoordNode *root_vertices;
} Polyhedra;

typedef struct nonunique_triangle_node{
    Triangle *tri;
    struct nonunique_triangle_node *next;
} nonunique_triangle_node;

typedef struct cube_gpu {
    int idx;
    int x;
    int y;
    int z;
    int one_apex[3];
    int two_apex[3];
    int _padding[6]; // ensure total size is 64 bytes
} cube_gpu __attribute__((aligned(64)));

typedef struct Triangle_GPU{
    TriangleVertex v1;
    TriangleVertex v2;
    TriangleVertex v3;
} Triangle_GPU;

typedef struct cube_vertices_points{
    CubeVertex coord;
    dim_t value;
} cube_vertices_points;

#ifdef __cplusplus
    }
#endif

#endif //TETRA