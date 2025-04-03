#ifndef UTILS_H
#define UTILS_H

#include "global.h"
#include <sys/stat.h>
#include <string.h>

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, double *origin);
void print_grid(const Dimensions *dim, const dim_t *grid);
void print_to_file(Triangle *triangle,int* count, double *origin);
void print_connections(Triangle *triangle, int*count);
void print_stack(StackNode *start);
void merge_files(char *atoms, char* conn);
void print_on_file(Polyhedra *p, char *name);
void print_for_stats(Polyhedra *p);
void print_on_separate_files(Polyhedra *p, char *molecule_name, char* molecule_path, int num_traingles);
void print_atoms_connections_separated(TriangleNode *curr, char *molecule_name, char* result_path, int num_traingles);

#endif // UTILS_H