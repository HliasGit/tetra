#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(){
    Dimensions dim;

    char folder_name[100] = "/home/elia/tesi/code/marching_tetrahedron/test/data/";
    char name[100] = "ala_eddi";
    char name_original[100];
    strcpy(name_original, name);
    char *path = strcat(folder_name, strcat(name, ".bin"));
    dim_t threshold = 0.02;

    printf("Creating surface from file '");
    printf(name);
    printf("'\nUsing threshold: %f\n", threshold);

    dim_t *grid;
    double origin[3];

    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};

    read_file(path, &dim, &grid, origin);

    verbose_call(print_grid(&dim, grid));

    void (*interpolation_function)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t);

    interpolation_function = &midpoint_interpol;    // Choose among midpoint_interpolation and linear interpolation

    Polyhedra p;
    p.triangles = NULL;
    p.vertices = NULL;
    p.root = NULL;

    TriangleNode *root;

    marching_tetrahedra(&dim, &grid, cube_decomposition, threshold, origin, interpolation_function, &p, &root);

    if (p.triangles == NULL) {
        fprintf(stderr, "No triangles have been generated\n");
        exit(-1);
    }

    if (p.root == NULL) {
        fprintf(stderr, "No vertices have been generated\n");
        exit(-1);
    }

    p.triangles = root;
    // print_on_file(&p, name_original);
    // print_for_stats(&p);
    print_on_separate_files(&p, name_original);

    free(grid),
    free_tree(p.root);
    free_list(p.triangles);

    // Print dimensions
    printf("Grid dimensions:\n\tX = %4ld\n\tY = %4ld\n\tZ = %4ld\n", dim.x_dim, dim.y_dim, dim.z_dim);
    return 0;
}