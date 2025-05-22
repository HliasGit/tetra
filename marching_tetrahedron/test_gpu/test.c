#include <marching_tetrahedron_gpu.h>

#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(int argc, char *argv[]) {

    Dimensions dim;

    char molecule_path[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name[100] = "adenosina";
    char molecule_path_original[100] = "/home/fs72740/evaglietti/tetra/marching_tetrahedron/data/";
    char molecule_name_original[100] = "adenosina";
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
    
    parallel_march_tetra(   &dim, d_grid, cube_decomposition, threshold, interpolation_function,
                            &p, &triangles_count, &vertex_counter, relevant_cubes_size, 
                            &d_relevant_cubes, &d_cube_points_coordinates);

    // marching_tetrahedra_list(&dim, &grid, cube_decomposition, threshold, origin,
    //     interpolation_function, &p, &triangles_count,
    //     &vertex_counter, relevant_cubes_size, results);
    
    // print_on_separate_files(&p, molecule_name_original, molecule_path_original, triangles_count);

    free_tree(p.root_vertices);
    free_list(p.triangles);

    return 0;
}