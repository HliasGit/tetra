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

typedef double dim_t;    // type of the scalar field values
typedef double coord_t;  // type for the coordinates (float or double)

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
} StackNode;

typedef struct TriangleNode{
    int vert1;
    int vert2;
    int vert3;
    struct TriangleNode *next;
} TriangleNode;

typedef struct VertexNode{
    TriangleVertex *vertex;
    struct VertexNode *next;
    int idx;
} VertexNode;

typedef struct Polyhedra{
    struct TriangleNode *triangles;
    VertexNode *vertices;
} Polyhedra;

#endif //TETRA