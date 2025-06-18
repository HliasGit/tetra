#include <stdio.h>

#include "mt_tetra_based.h"
#include "global.h"

#include "marching_tetrahedron_gpu.h"


int main(int argc, char *argv[]) {

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <threshold> <molecule_name>\n", argv[0]);
        return -1;
    }

    Dimensions dim;

    char molecule_path[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/float/";
    char molecule_name[100];
    char molecule_name_original[100];
    char molecule_path_original[100];
    strcpy(molecule_name, argv[2]);
    strcpy(molecule_path_original, molecule_path);
    strcpy(molecule_name_original, molecule_name);
    char *path = strcat(molecule_path, strcat(molecule_name, ".bin"));
    
    char *endptr;
    dim_t threshold = strtod(argv[1], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid threshold value: %s\n", argv[1]);
        return -1;
    }
    
    printf("Creating surface from file '");
    printf(molecule_name);
    printf("'\nUsing threshold: %f\n", threshold);
    
    printf("Path to molecule file: %s\n", path);

    dim_t *grid;

    dim_t origin[3];
    
    double time = 0;
    
    int CD1[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};
    int CD2[20] = {1,2,5,3,3,5,7,8,2,3,4,8,2,4,5,8,2,3,5,8};
    
    read_file(path, &dim, &grid, origin);

    printf("Grid dimensions: x = %d, y = %d, z = %d\n", dim.x_dim, dim.y_dim, dim.z_dim);


    int pairs[48] = {1,2,1,3,1,4,2,2,1,3,1,4,2,2,3,3,1,4,2,2,3,3,4,4,1,4,2,4,3,3,1,4,2,4,3,4,1,4,2,4,1,3,2,4,2,3,1,3};

    load_to_const_tetra(&dim);
    
    ////////////////////////// MARCH TETRA //////////////////////////
    
    size_t size = dim.x_dim * dim.y_dim * dim.z_dim;
    printf("Total grid size: %zu\n", size);
    
    dim_t *d_grid;

    allocate_d_grid(&d_grid, grid, size);

    int number_relevant_cubes;
    cube_gpu_SoA *d_relevant_cubes;
    cube_vertices_points_SoA *d_cube_points_coordinates;

    int act_val_vec[25] = {0,1,7,6,0,0,2,5,0,0,0,3,0,6,0,0,4,0,6,0,0,0,0,6,0};

    remove_unnecessary_tetrahedra(  d_grid, size, threshold,
                                    &number_relevant_cubes, &d_relevant_cubes, &time,
                                    CD1, CD2, pairs, act_val_vec);

    printf("\n");
    printf("Total GPU time: %f ms\n", time);


    return 0;
}