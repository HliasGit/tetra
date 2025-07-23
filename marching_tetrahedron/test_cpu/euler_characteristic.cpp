#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

#include <unordered_map>
#include <list>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{

    if (argc < 4)
    {
        fprintf(stderr, "Usage: %s <threshold> <molecule_name> <midpoint|linear> \n", argv[0]);
        return -1;
    }

    if (strcmp(argv[3], "midpoint") != 0 && strcmp(argv[3], "linear") != 0)
    {
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
    dim_t threshold = strtod(argv[1], &endptr);
    if (*endptr != '\0')
    {
        fprintf(stderr, "Invalid threshold value: %s\n", argv[1]);
        return -1;
    }

    printf("Creating surface from file '");
    printf("%s", molecule_name);
    printf("'\nUsing threshold: %f\n", threshold);

    printf("Path to molecule file: %s\n", path);


    dim_t *grid;
    dim_t origin[3];

    int cube_decomposition[20] = {4, 6, 7, 8, 1, 5, 6, 7, 1, 3, 4, 7, 1, 2, 4, 6, 1, 4, 6, 7};

    read_file(path, &dim, &grid, origin);

    verbose_call(print_grid(&dim, grid));

    void (*interpolation_function)(TriangleVertex *, CubeVertex *, CubeVertex *, dim_t *, dim_t *, dim_t);

    if (strcmp(argv[3], "midpoint") == 0)
    {
        interpolation_function = &midpoint_interpol;
        printf("Using mipoint interpolation\n");
    }
    else
    {
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

    if (p.triangles == NULL)
    {
        fprintf(stderr, "No triangles have been generated\n");
        exit(-1);
    }

    if (p.root_vertices == NULL)
    {
        fprintf(stderr, "No vertices have been generated\n");
        exit(-1);
    }

    printf("Molecule name original: %s\n", molecule_name_original);
    print_on_separate_files(&p, molecule_name_original, molecule_path_original, triangles_count);

    clock_t start_euler, end_euler;
    start_euler = clock();
    std::unordered_map<int, int> edge_count;

    int counter = 0;

    while (p.triangles != NULL)
    {
        int one = p.triangles->vert1->index;
        int two = p.triangles->vert2->index;
        int thr = p.triangles->vert3->index;

        auto add_edge = [&edge_count, &counter](int a, int b) {
            if (a > b) std::swap(a, b);
            int coupled = a * 100000000 + b;
            if (edge_count.find(coupled) == edge_count.end()) {
                edge_count[coupled] = 1;
                counter++;
            }
        };

        add_edge(one, two);
        add_edge(two, thr);
        add_edge(thr, one);

        p.triangles = p.triangles->next;
    }

    printf("# edges: %d\n", counter);

    printf("Euler characteristics: %lu\n", triangles_count+vertex_counter-counter);
    end_euler = clock();

    double time_taken = double(end_euler - start_euler) / CLOCKS_PER_SEC;
    printf("Time taken to compute Euler characteristics: %f seconds\n", time_taken);

    free(grid),
    free_tree(p.root_vertices);
    free_list(p.triangles);

    // Print dimensions
    printf("Grid dimensions:\n\tX = %4ld\n\tY = %4ld\n\tZ = %4ld\n", dim.x_dim, dim.y_dim, dim.z_dim);
    return 0;
}