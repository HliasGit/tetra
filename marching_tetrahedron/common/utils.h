#pragma once

#include "global.h"
#include "struct.h"


#ifdef __cplusplus
extern "C" {
#endif

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, dim_t *origin);
void print_on_separate_files(Polyhedra *p, char *molecule_name, char *molecule_path, int num_triangles);
void print_atoms_connections_separated(TriangleNode *start_triangles, char *molecule_name, char *result_path, int num_traingles);
void print_triangles_cpu(   nonunique_triangle_node *start, char *molecule_name, char *molecule_path);

#ifdef __cplusplus
}
#endif