#include "marching_tetrahedron.h"

int main(){
    Dimensions dim;

    char folder_name[100] = "/home/elia/tesi/code/marching_tetrahedron/test/";
    char *path = strcat(folder_name, "file.bin");

    dim_t threshold = 0.5;
    dim_t *grid;

    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};

    read_file(path, &dim, &grid);
    // verbose_call(print_grid(&dim, grid));
    normalize_grid(&dim, &grid, threshold);
    // verbose_call(print_grid(&dim, grid));

    int count;
    marching_tetrahedra(&dim, &grid, cube_decomposition, &count);

    verbose_print("Count: %d\n", count);

    return 0;
}