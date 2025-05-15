#pragma once

#include "global.h"
#include <string.h>
#include "struct.h"


#ifdef __cplusplus
extern "C" {
#endif

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, double *origin);
void print_on_separate_files(Polyhedra *p, char *molecule_name, char *molecule_path, int num_triangles);


#ifdef __cplusplus
}
#endif