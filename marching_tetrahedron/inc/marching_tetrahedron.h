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
} StackNode;

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor);
void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold);
void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count);
void print_grid(const Dimensions *dim, const dim_t *grid);
void find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k, CubeVertex **coordinates);
void push_into_stack(StackNode **start, dim_t value, CubeVertex vert);
void free_stack(StackNode **start);
void print_stack(StackNode *start);
int get_action_value(StackNode *start);
int *get_pairs(int action_val);
Triangle *make_triangle(StackNode *stack, int *pairs);
CubeVertex *get_coordinate_by_idx(StackNode *start, int idx);
void print_to_file(Triangle *triangle,int* count);

#endif //TETRA