#pragma once

#include <cuda_runtime.h>
#include "global.h"


__global__ void hello_kernel();

__global__ void remove_unnecessary_cubes_kernel(double *grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, int* d_list);