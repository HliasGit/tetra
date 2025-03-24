#include "global.h"

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, double *origin);
void print_grid(const Dimensions *dim, const dim_t *grid);
void print_to_file(Triangle *triangle,int* count, double *origin);
void print_connections(Triangle *triangle, int*count);
void print_stack(StackNode *start);
void merge_files(char *atoms, char* conn);