#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef float dim_t;    // type of the scalar field values
typedef float coord_t;  // type for the coordinates (float or double)

typedef struct dimensions{
    size_t x_dim;
    size_t y_dim;
    size_t z_dim;
} Dimensions;

typedef struct vertex{
    coord_t x;
    coord_t y;
    coord_t z;
} Vertex;

typedef struct triangle{
    Vertex *v1;
    Vertex *v2;
    Vertex *v3;
} Triangle;

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor);
void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold);
void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count);
void print_grid(const Dimensions *dim, const dim_t *grid);
void find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k, Vertex **coordinates);