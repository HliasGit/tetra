#include <marching_tetrahedron_gpu.h>
#include <stdio.h>

#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(int argc, char *argv[]) {

    Dimensions dim;

    char molecule_path[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name[100] = "squared_torus";
    char molecule_path_original[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name_original[100] = "squared_torus";
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
    
    ////////////////////////// MARCH TETRA //////////////////////////
    
    size_t size = dim.x_dim * dim.y_dim * dim.z_dim;
    
    dim_t *d_grid;

    allocate_d_grid(&d_grid, grid, size);

    int relevant_cubes_size;
    cube_gpu *d_relevant_cubes;
    cube_vertices_points *d_cube_points_coordinates;

    remove_unnecessary_cubes(   d_grid, size, threshold, &dim,
                                &relevant_cubes_size, &d_relevant_cubes,
                                &d_cube_points_coordinates);

    printf("\n");
    
    void (*interpolation_function)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t);
    
    interpolation_function = &midpoint_interpol;
    
    Polyhedra p;
    p.triangles = NULL;
    p.root_vertices = NULL;
    
    
    size_t triangles_count;
    size_t vertex_counter;

    printf("Number of processed cubes: %d\n", relevant_cubes_size);

    Triangle_GPU *triangles;

    int act_val_vec[25] = {0,1,7,6,0,0,2,5,0,0,0,3,0,6,0,0,4,0,6,0,0,0,0,6,0};

    int pairs[48] = {1,2,1,3,1,4,2,2,1,3,1,4,2,2,3,3,1,4,2,2,3,3,4,4,1,4,2,4,3,3,1,4,2,4,3,4,1,4,2,4,1,3,2,4,2,3,1,3};
    
    parallel_march_tetra(   &dim, d_grid, cube_decomposition, threshold, interpolation_function,
                            &p, &triangles_count, &vertex_counter, relevant_cubes_size, 
                            &d_relevant_cubes, &d_cube_points_coordinates, act_val_vec,
                            pairs, &triangles);

 
    // for(int i=0; i<150000; i++){
    //     if (triangles[i].v1.x == 0.0 && triangles[i].v1.y == 0.0 && triangles[i].v1.z == 0.0 &&
    //         triangles[i].v2.x == 0.0 && triangles[i].v2.y == 0.0 && triangles[i].v2.z == 0.0 &&
    //         triangles[i].v3.x == 0.0 && triangles[i].v3.y == 0.0 && triangles[i].v3.z == 0.0) {
    //         continue;
    //     }n
    //     printf("v1: %lf %lf %lf\n", triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
    //     printf("v2: %lf %lf %lf\n", triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
    //     printf("v3: %lf %lf %lf\n", triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
    // }

    int empty = 0;

    FILE *pdb_file = fopen("output.pdb", "w");
    if (pdb_file == NULL) {
        fprintf(stderr, "Error opening output.pdb for writing\n");
    } else {
        int atom_serial = 1;
        for (int i = 0; i < relevant_cubes_size; i++) {
            if (triangles[i].v1.x == 0.0 && triangles[i].v1.y == 0.0 && triangles[i].v1.z == 0.0 &&
                triangles[i].v2.x == 0.0 && triangles[i].v2.y == 0.0 && triangles[i].v2.z == 0.0 &&
                triangles[i].v3.x == 0.0 && triangles[i].v3.y == 0.0 && triangles[i].v3.z == 0.0) {
                empty++;
                continue;
            }
            fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial++, triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
            fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial++, triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
            fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", atom_serial++, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
        }
        fclose(pdb_file);
    }

    printf("EMPTY:  %d\n", empty);


    free(triangles);
    
    // print_on_separate_files(&p, molecule_name_original, molecule_path_original, triangles_count);

    free_tree(p.root_vertices);
    free_list(p.triangles);

    return 0;
}