#include <marching_tetrahedron_gpu.h>
#include <stdio.h>

#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(int argc, char *argv[]) {

    Dimensions dim;

    char molecule_path[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/float/";
    char molecule_name[100] = "adenosine";
    char molecule_path_original[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/float/";
    char molecule_name_original[100] = "adenosine";
    char *path = strcat(molecule_path, strcat(molecule_name, ".bin"));
    
    dim_t threshold = 4e-4;
    
    printf("Creating surface from file '");
    printf(molecule_name);
    printf("'\nUsing threshold: %f\n", threshold);
    
    printf("Path to molecule file: %s\n", path);

    dim_t *grid;

    dim_t origin[3];
    
    
    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};
    
    read_file(path, &dim, &grid, origin);

    printf("Grid dimensions: x = %d, y = %d, z = %d\n", dim.x_dim, dim.y_dim, dim.z_dim);
    
    ////////////////////////// MARCH TETRA //////////////////////////
    
    size_t size = dim.x_dim * dim.y_dim * dim.z_dim;
    
    dim_t *d_grid;

    allocate_d_grid(&d_grid, grid, size);

    cube_gpu *d_relevant_cubes;
    cube_vertices_points *d_cube_points_coordinates;

    double time = 0.0;

    skip_preprocessing(   d_grid, size, threshold, &dim, &d_relevant_cubes,
                                &d_cube_points_coordinates, &time);

    printf("\n");
    
    size_t triangles_count;
    size_t vertex_counter;

    printf("Number of processed cubes: %d\n", size);

    Triangle_GPU *triangles = NULL;

    int act_val_vec[25] = {0,1,7,6,0,0,2,5,0,0,0,3,0,6,0,0,4,0,6,0,0,0,0,6,0};

    int pairs[48] = {1,2,1,3,1,4,2,2,1,3,1,4,2,2,3,3,1,4,2,2,3,3,4,4,1,4,2,4,3,3,1,4,2,4,3,4,1,4,2,4,1,3,2,4,2,3,1,3};
    
    int total_triangles = 0;

    parallel_march_tetra(   &dim, d_grid, cube_decomposition, threshold, &triangles_count,
                            &vertex_counter, size, 
                            &d_relevant_cubes, &d_cube_points_coordinates, act_val_vec,
                            pairs, &triangles, &total_triangles, &time);
 
    printf("Total GPU time: %f ms\n", time);
    // print_triangles(triangles, &total_triangles, molecule_name_original, molecule_path_original);


    free(triangles);


    return 0;
}