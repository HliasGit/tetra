#include <marching_tetrahedron_gpu.h>
#include <stdio.h>

#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(int argc, char *argv[]) {

    Dimensions dim;

    char molecule_path[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name[100] = "9mxc_atom";
    char molecule_path_original[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name_original[100] = "9mxc_atom";
    char *path = strcat(molecule_path, strcat(molecule_name, ".bin"));
    
    double threshold = 4e-4;
    
    printf("Creating surface from file '");
    printf(molecule_name);
    printf("'\nUsing threshold: %f\n", threshold);
    
    printf("Path to molecule file: %s\n", path);

    dim_t *grid;

    double origin[3];
    
    
    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};
    
    read_file(path, &dim, &grid, origin);

    printf("Grid dimensions: x = %d, y = %d, z = %d\n", dim.x_dim, dim.y_dim, dim.z_dim);
    
    ////////////////////////// MARCH TETRA //////////////////////////
    
    size_t size = dim.x_dim * dim.y_dim * dim.z_dim;
    
    dim_t *d_grid;

    allocate_d_grid(&d_grid, grid, size);

    int number_relevant_cubes;
    cube_gpu *d_relevant_cubes;
    cube_vertices_points *d_cube_points_coordinates;

    remove_unnecessary_cubes(   d_grid, size, threshold, &dim,
                                &number_relevant_cubes, &d_relevant_cubes,
                                &d_cube_points_coordinates);

    printf("\n");
    
    size_t triangles_count;
    size_t vertex_counter;

    printf("Number of processed cubes: %d\n", size);
    printf("Number of relevant cubes: %d\n", number_relevant_cubes);

    Triangle_GPU *triangles = NULL;

    int act_val_vec[25] = {0,1,7,6,0,0,2,5,0,0,0,3,0,6,0,0,4,0,6,0,0,0,0,6,0};

    int pairs[48] = {1,2,1,3,1,4,2,2,1,3,1,4,2,2,3,3,1,4,2,2,3,3,4,4,1,4,2,4,3,3,1,4,2,4,3,4,1,4,2,4,1,3,2,4,2,3,1,3};
    
    parallel_march_tetra(   &dim, d_grid, cube_decomposition, threshold, &triangles_count,
                            &vertex_counter, number_relevant_cubes, 
                            &d_relevant_cubes, &d_cube_points_coordinates, act_val_vec,
                            pairs, &triangles);



 

    // int empty = 0;
    // printf("number_relevant_cubes: %d\n", number_relevant_cubes);

    // FILE *pdb_file = fopen("output.pdb", "w");
    // if (pdb_file == NULL) {
    //     fprintf(stderr, "Error opening output.pdb for writing\n");
    // } else {
    //     int atom_serial = 1;
    //     for (int i = 0; i < number_relevant_cubes*7; i++) {
    //         if (triangles[i].v1.x == 0.0 && triangles[i].v1.y == 0.0 && triangles[i].v1.z == 0.0 &&
    //             triangles[i].v2.x == 0.0 && triangles[i].v2.y == 0.0 && triangles[i].v2.z == 0.0 &&
    //             triangles[i].v3.x == 0.0 && triangles[i].v3.y == 0.0 && triangles[i].v3.z == 0.0) {
    //             empty++;
    //             continue;
    //         }
    //         if(atom_serial >= 100000)
    //             continue;
    //         fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
    //         // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
    //         atom_serial++;
    //         fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
    //         // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
    //         atom_serial++;
    //         fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
    //         // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
    //         atom_serial++;
    //         fprintf(pdb_file, "CONECT%5d%5d\n", atom_serial - 3, atom_serial - 2);
    //         fprintf(pdb_file, "CONECT%5d%5d\n", atom_serial - 2, atom_serial - 1);
    //         fprintf(pdb_file, "CONECT%5d%5d\n", atom_serial - 3, atom_serial - 1);

    //         // printf("CONECT%5d%5d\n", atom_serial - 3, atom_serial - 2);
    //         // printf("CONECT%5d%5d\n", atom_serial - 2, atom_serial - 1);
    //         // printf("CONECT%5d%5d\n", atom_serial - 3, atom_serial - 1);
    //     }
    //     fclose(pdb_file);
    // }

    print_triangles(triangles, &number_relevant_cubes, molecule_name_original, molecule_path_original);


    free(triangles);


    return 0;
}