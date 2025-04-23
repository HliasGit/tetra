#ifndef UTILS_H
#define UTILS_H

#include "global.h"
#include <sys/stat.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, double *origin);
void print_grid(const Dimensions *dim, const dim_t *grid);
void print_stack(StackNode *start);
void print_for_stats(Polyhedra *p);
void print_on_separate_files(Polyhedra *p, char *molecule_name, char* molecule_path, int num_traingles);
void print_atoms_connections_separated(TriangleNode *curr, char *molecule_name, char* result_path, int num_traingles);
void print_to_console_traingles(TriangleNode* start);


#ifdef __cplusplus
}
#endif

#endif // UTILS_H