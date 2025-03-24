#include "global.h"

void push_into_stack(StackNode **start, dim_t value, CubeVertex vert);
void free_stack(StackNode **start);
CubeVertex *get_coordinate_by_idx(StackNode *start, int idx);
dim_t *get_value_by_idx(StackNode *start, int idx);