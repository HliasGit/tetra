#include <stdio.h>

#include "marching_tetrahedron_gpu.h"
#include "mt_noatom.h"
#include "global.h"

#include <vector_types.h>


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
    
    float time = 0;
    
    int CD1[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};
    int act_val_vec[25] = {0,1,7,6,0,0,2,5,0,0,0,3,0,6,0,0,4,0,6,0,0,0,0,6,0};
    int pairs[48] = {1,2,1,3,1,4,2,2,1,3,1,4,2,2,3,3,1,4,2,2,3,3,4,4,1,4,2,4,3,3,1,4,2,4,3,4,1,4,2,4,1,3,2,4,2,3,1,3};
    
    read_file(path, &dim, &grid, origin);

    printf("Grid dimensions: x = %d, y = %d, z = %d\n", dim.x_dim, dim.y_dim, dim.z_dim);

    load_to_const_tetra_flags(&dim);
    
    //////////////////////// MARCH TETRA //////////////////////////
    
    size_t cubes_in_domain = dim.x_dim * dim.y_dim * dim.z_dim;
    printf("Total grid cubes_in_domain: %zu\n", cubes_in_domain);
    
    dim_t *d_grid;
    bool *d_flags;
    cube_gpu_SoA *d_relevant_cubes;

    allocate_d_grid(&d_grid, grid, cubes_in_domain);
    allocate_d_flags(&d_flags, cubes_in_domain);

    remove_unnecessary_cubes_flag(   d_grid,
                                d_flags,
                                cubes_in_domain,
                                threshold,
                                &d_relevant_cubes,
                                &time);

    // // remove_unnecessary_tetrahedra(  d_grid,
    // //                                 cubes_in_domain,
    // //                                 CD1,
    // //                                 threshold,
    // //                                 &d_active_tetrahedra, 
    // //                                 &d_num_active_tetrahedra,
    // //                                 &time);

    // tetra_mt(d_num_active_tetrahedra, d_active_tetrahedra, threshold, pairs, act_val_vec, &time);

    printf("\n");
    printf("Total GPU time: %f ms\n", time);


    return 0;
}