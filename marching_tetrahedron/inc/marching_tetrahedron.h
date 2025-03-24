#include "global.h"
#include "struct.h"
#include "triangles.h"
#include "utils.h"

void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold);
void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count, dim_t threshold, double *origin, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t));
void find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k, CubeVertex **coordinates);
int get_action_value(StackNode *start, dim_t threshold);
int *get_pairs(int action_val);