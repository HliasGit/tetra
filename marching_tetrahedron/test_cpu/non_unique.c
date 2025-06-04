#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

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

    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};

    read_file(path, &dim, &grid, origin);

    verbose_call(print_grid(&dim, grid));

    void (*interpolation_function)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t);

    nonunique_triangle_node *start = NULL;

    size_t triangles_count;
    size_t vertex_counter;

    clock_t start_marching, end_marching;
    start_marching = clock();
    marching_tetrahedra_notunique(&dim, &grid, cube_decomposition, threshold, origin, &triangles_count, &start);
    
    end_marching = clock();
    double time_spent = (double)(end_marching - start_marching) / CLOCKS_PER_SEC;
    printf("Took %f seconds for the marhcing tetrahedron computation\n", time_spent);

    printf("Molecule name original: %s\n", molecule_name_original);

    // nonunique_triangle_node *current = start;
    // while (current != NULL) {
    //     printf("Triangle vertices:\n");
    //     for (int i = 0; i < 3; ++i) {
    //         printf("\tVertex %d: (%f, %f, %f)\n", i,
    //            current->tri->v1->x,
    //            current->tri->v2->y,
    //            current->tri->v3->z);
    //     }
    //     printf("\n");
    //     current = current->next;
    // }
    
    // print_triangles_cpu(start, molecule_name_original, molecule_path_original);

    printf("Grid points have been written to grid.pdb\n");

    free(grid);

    // Print dimensions
    printf("Grid dimensions:\n\tX = %4ld\n\tY = %4ld\n\tZ = %4ld\n", dim.x_dim, dim.y_dim, dim.z_dim);
    return 0;
}
