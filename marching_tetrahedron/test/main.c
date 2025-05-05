#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(int argc, char *argv[]) {

    if (argc < 4) {
        fprintf(stderr, "Usage: %s <threshold> <molecule_name> <midpoint|linear> \n", argv[0]);
        return -1;
    }

    if(strcmp(argv[3], "midpoint") != 0 && strcmp(argv[3], "linear") != 0){
        fprintf(stderr, "Choose among <midpoint|linear> as interpolation method \n");
        return -1;
    }

    Dimensions dim;

    char molecule_path[100] = "/home/elia/tesi/tetra/marching_tetrahedron/data/";
    char molecule_name[100];
    char molecule_name_original[100];
    char molecule_path_original[100];
    strcpy(molecule_name, argv[2]);
    strcpy(molecule_path_original, molecule_path);
    strcpy(molecule_name_original, molecule_name);
    char *path = strcat(molecule_path, strcat(molecule_name, ".bin"));
    
    char *endptr;
    double threshold = strtod(argv[1], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Invalid threshold value: %s\n", argv[1]);
        return -1;
    }
    
    printf("Creating surface from file '");
    printf(molecule_name);
    printf("'\nUsing threshold: %f\n", threshold);
    
    printf("Path to molecule file: %s\n", path);

    dim_t *grid;
    double origin[3];

    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};

    read_file(path, &dim, &grid, origin);

    verbose_call(print_grid(&dim, grid));

    void (*interpolation_function)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t);

    if(strcmp(argv[3], "midpoint") == 0){
        interpolation_function = &midpoint_interpol;
        printf("Using mipoint interpolation\n");
    } else {
        interpolation_function = &linear_interpol;    
        printf("Using linear interpolation\n");
    }

    Polyhedra p;
    p.triangles = NULL;
    p.root_vertices = NULL;

    size_t triangles_count;
    size_t vertex_counter;

    clock_t start_marching, end_marching;
    start_marching = clock();
    marching_tetrahedra(&dim, &grid, cube_decomposition, threshold, origin,
                        interpolation_function, &p, &triangles_count,
                        &vertex_counter);
    end_marching = clock();
    double time_spent = (double)(end_marching - start_marching) / CLOCKS_PER_SEC;
    printf("Took %f seconds for the marhcing tetrahedron computation\n", time_spent);

    if (p.triangles == NULL) {
        fprintf(stderr, "No triangles have been generated\n");
        exit(-1);
    }

    if (p.root_vertices == NULL) {
        fprintf(stderr, "No vertices have been generated\n");
        exit(-1);
    }

    // print_for_stats(&p);
    printf("Molecule name original: %s\n", molecule_name_original);
    clock_t start_printing, end_printing;
    start_printing = clock();
    print_on_separate_files(&p, molecule_name_original, molecule_path_original, triangles_count);
    end_printing = clock();
    double time_spent_printing = (double)(end_printing - start_printing) / CLOCKS_PER_SEC;
    printf("Took %f seconds for the printing\n", time_spent_printing);

    // print_to_console_traingles(p.triangles);

    FILE *grid_pdb_file = fopen("grid.pdb", "w");
    if (!grid_pdb_file)
    {
        fprintf(stderr, "Failed to open grid.pdb for writing\n");
        return -1;
    }

    // Iterate through the grid dimensions and write the points
    int atom_index = 1;
    for (size_t x = 0; x < dim.x_dim; ++x)
    {
        for (size_t y = 0; y < dim.y_dim; ++y)
        {
            for (size_t z = 0; z < dim.z_dim; ++z)
            {
                fprintf(grid_pdb_file, "ATOM  %5d  C   GRD A   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n",
                        atom_index, (double)x, (double)y, (double)z);
                atom_index++;
            }
        }
    }

    fclose(grid_pdb_file);
    printf("Grid points have been written to grid.pdb\n");

    free(grid),
    free_tree(p.root_vertices);
    free_list(p.triangles);

    // Print dimensions
    printf("Grid dimensions:\n\tX = %4ld\n\tY = %4ld\n\tZ = %4ld\n", dim.x_dim, dim.y_dim, dim.z_dim);
    return 0;
}