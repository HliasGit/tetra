#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "global.h"

#ifdef __cplusplus
extern "C" {
#endif

Triangle *make_triangle(StackNode *stack, int *pairs, bool two_triangles, dim_t threshold, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t), bool swap);
void midpoint_interpol(TriangleVertex *vtx,CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold);
void linear_interpol(TriangleVertex *vtx,CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold);
bool coordinate_less_than(TriangleVertex *v1, TriangleVertex *v2);
bool coordinate_equals(TriangleVertex *v1, TriangleVertex *v2);
void characterize_triangle(Triangle *tri, StackNode *stack, int *pairs);

#ifdef __cplusplus
}
#endif

#endif // TRIANGLES_H