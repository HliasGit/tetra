#pragma once

#include "global.h"
#include "struct_gpu.cuh"

void load_to_const_tetra(Dimensions *dimensions);

void remove_unnecessary_tetrahedra( dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                    int *number_relevant_cubes,
                                    cube_gpu_SoA **d_relevant_cubes,
                                    double *time, int *CD1, int *CD2, int *pairs, int* act_val_vec);