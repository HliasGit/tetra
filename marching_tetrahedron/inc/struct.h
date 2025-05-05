#ifndef STRUCT_H
#define STRUCT_H

#include "global.h"
#include "triangles.h"

#ifdef __cplusplus
extern "C" {
#endif

void push_into_stack(StackNode **start, dim_t value, CubeVertex vert, int point, int exp_pos, int* permutations);
void free_stack(StackNode **start);
CubeVertex *get_coordinate_by_idx(StackNode *start, int idx);
dim_t *get_value_by_idx(StackNode *start, int idx);
void push_triangle(Polyhedra **p, Triangle *triangle, size_t *vertex_counter);
TriangleCoordNode *add(TriangleCoordNode **root, double *full_coordinate, size_t *idx);
void free_tree(TriangleCoordNode *TriangleCoordNode);
void free_list(TriangleNode *start);
void reverse_list(TriangleNode **head);  

#ifdef __cplusplus
}
#endif

#endif // STRUCT_H