#ifndef MARCHING_TETRAHEDRON_H
#define MARCHING_TETRAHEDRON_H

#include "global.h"
#include "struct.h"
#include "triangles.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold);
void marching_tetrahedra(   Dimensions *dim, dim_t **grid, int *cube_decomposition, dim_t threshold,
                            double *origin, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*,
                            dim_t*, dim_t*, dim_t), Polyhedra *p, size_t *triangles_count,
                            size_t *vertex_counter);
bool find_coordinates(  int idx, const int point, const size_t i, const size_t j, const size_t k,
                        CubeVertex **coordinates);
int get_action_value(StackNode *start, dim_t threshold);
int *get_pairs(int action_val);
bool parity(size_t i, size_t j, size_t k);
double three_det(double mat[3][3]);
bool tetrahedron_determinant(CubeVertex *coordinates);

#ifdef __cplusplus
}
#endif

#endif // MARCHING_TETRAHEDRON_H
