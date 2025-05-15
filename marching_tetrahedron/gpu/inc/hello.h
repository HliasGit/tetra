#pragma once
#include <stdio.h>
#include "global.h"

#ifdef __cplusplus
    extern "C"{
#endif

void hello_f();
int remove_unnecessary_cubes(double *grid, size_t size, double threshold, Dimensions *dim, int ** results);

#ifdef __cplusplus
    }
#endif