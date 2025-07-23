#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
#include "global.h"
#include "struct_gpu.cuh"

void remove_unnecessary_cubes_flag( dim_t* d_grid,
                                    bool *d_flags,
                                    size_t size,
                                    double threshold,
                                    cube_gpu_SoA** d_relevant_cubes,
                                    float *time
                                );
void load_to_const_tetra_flags(Dimensions *dimensions);
void allocate_d_flags(bool **d_flags, size_t size);