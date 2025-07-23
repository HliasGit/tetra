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
#include <string>

void make_cpp_vertex_list(TriangleCoordNode *node, std::list<TriangleVertex *> &list)
{
    if (node == NULL)
    {
        return;
    }

    make_cpp_vertex_list(node->next_list, list);
    make_cpp_vertex_list(node->next_level, list);

    TriangleVertex *vtx = (TriangleVertex *)malloc(sizeof(TriangleVertex));

    vtx->x = node->coordinate1;
    vtx->y = node->coordinate2;
    vtx->z = node->coordinate3;

    list.push_back(vtx);
}

TriangleVertex sub_vertices(TriangleVertex *v1, TriangleVertex *v2)
{
    TriangleVertex result;
    result.x = v1->x - v2->x;
    result.y = v1->y - v2->y;
    result.z = v1->z - v2->z;
    return result;
}

TriangleVertex cross_product_vertices(TriangleVertex *v1, TriangleVertex *v2)
{
    TriangleVertex result;
    result.x = v1->y * v2->z - v1->z * v2->y;
    result.y = v1->z * v2->x - v1->x * v2->z;
    result.z = v1->x * v2->y - v1->y * v2->x;
    return result;
}

TriangleVertex scalar_product(TriangleVertex* v2, double scalar){
    TriangleVertex result;
    result.x = scalar * v2->x;
    result.y = scalar * v2->y;
    result.z = scalar * v2->z;
    return result;
}

double length_squared(TriangleVertex *v1, TriangleVertex *v2)
{
    return pow(v2->x - v1->x, 2) + pow(v2->y - v1->y, 2) + pow(v2->z - v1->z, 2);
}

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
    if (*endptr != '\0')
    {
        fprintf(stderr, "Invalid threshold value: %s\n", argv[1]);
        return -1;
    }

    printf("Creating surface from file '");
    printf(molecule_name);
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

    //////////////////////////////////////////// NORMALS ////////////////////////////////////////////

    clock_t start_normals, end_normals;
    start_normals = clock();

    std::unordered_map<int, TriangleVertex> vertices;
    std::unordered_map<int, TriangleVertex> normals;


    TriangleNode *curr = p.triangles;

    int N = 99998;
    int count = 0;
    int file_number = 0;

    struct stat st;

    std::string folder_name =   std::string(molecule_path_original) + "normals/";

    if (stat(folder_name.c_str(), &st) == -1)
    {
        if (mkdir(folder_name.c_str(), 0700) != 0)
        {
            fprintf(stderr, "Failed to create directory: %s\n", folder_name.c_str());
            return -1;
        }
    }

    std::cout <<  "FOLDER NAME" << folder_name << std::endl;

    FILE *fptr;

    while(curr != NULL){
        TriangleVertex vtx1 = {curr->vert1->coordinate1, curr->vert1->coordinate2, curr->vert1->coordinate3};
        TriangleVertex vtx2 = {curr->vert2->coordinate1, curr->vert2->coordinate2, curr->vert2->coordinate3};
        TriangleVertex vtx3 = {curr->vert3->coordinate1, curr->vert3->coordinate2, curr->vert3->coordinate3};

        TriangleVertex v1 = sub_vertices(&vtx2, &vtx1);
        TriangleVertex v2 = sub_vertices(&vtx3, &vtx1);
        TriangleVertex cross = cross_product_vertices(&v1, &v2);

        double l1 = length_squared(&vtx2, &vtx3);
        double l2 = length_squared(&vtx1, &vtx3);
        double l3 = length_squared(&vtx1, &vtx2);

        double weight = l1 + l2 + l3;
        double weight1 = 1/(weight*(l2+l3));
        double weight2 = 1/(weight*(l3+l1));
        double weight3 = 1/(weight*(l1+l2));
        

        vertices[curr->vert1->index].x += scalar_product(&cross, weight1).x;
        vertices[curr->vert1->index].y += scalar_product(&cross, weight1).y;
        vertices[curr->vert1->index].z += scalar_product(&cross, weight1).z;

        vertices[curr->vert2->index].x += scalar_product(&cross, weight2).x;
        vertices[curr->vert2->index].y += scalar_product(&cross, weight2).y;
        vertices[curr->vert2->index].z += scalar_product(&cross, weight2).z;

        vertices[curr->vert3->index].x += scalar_product(&cross, weight3).x;
        vertices[curr->vert3->index].y += scalar_product(&cross, weight3).y;
        vertices[curr->vert3->index].z += scalar_product(&cross, weight3).z;

        normals[curr->vert1->index].x = curr->vert1->coordinate1;
        normals[curr->vert1->index].y = curr->vert1->coordinate2;
        normals[curr->vert1->index].z = curr->vert1->coordinate3;

        normals[curr->vert2->index].x = curr->vert2->coordinate1;
        normals[curr->vert2->index].y = curr->vert2->coordinate2;
        normals[curr->vert2->index].z = curr->vert2->coordinate3;

        normals[curr->vert3->index].x = curr->vert3->coordinate1;
        normals[curr->vert3->index].y = curr->vert3->coordinate2;
        normals[curr->vert3->index].z = curr->vert3->coordinate3;

        curr = curr->next;
    }

    count = 0;
    int idx = 0;
    int local_idx = 0;
    int local_count = 0;
    file_number = 0;

    while (idx != vertex_counter)
    {
        if(count%N == 0){
            std::string file_name = folder_name + std::to_string(file_number) + ".pdb";
            fptr = fopen(file_name.c_str(), "w");
            local_idx = 0;
            local_count = 0;

            std::cout << file_name << std::endl;
        }        

        fprintf(fptr, "ATOM  %5d  N   NOR A   1    %8.3f%8.3f%8.3f  1.00  0.00           N\n",
            local_count,  normals[idx].x,
                    normals[idx].y,
                    normals[idx].z);

        fprintf(fptr, "ATOM  %5d  N   NOR A   1    %8.3f%8.3f%8.3f  1.00  0.00           N\n",
            local_count+1,  vertices[idx].x+normals[idx].x,
                        vertices[idx].y+normals[idx].y,
                        vertices[idx].z+normals[idx].z);

        fprintf(fptr, "CONECT%5d%5d\n", local_count, local_count + 1);
        
        idx++;
        local_idx++;
        local_count += 2;
        count += 2;
        
        if(((count+2)%N == 0)){
            file_number++;
            fclose(fptr);
        }
    }

    end_normals = clock();
    time_spent = (double)(end_normals - start_normals) / CLOCKS_PER_SEC;
    printf("Took %f seconds for the marhcing tetrahedron computation\n", time_spent);

    free(grid),
    free_tree(p.root_vertices);
    free_list(p.triangles);

    // Print dimensions
    printf("Grid dimensions:\n\tX = %4ld\n\tY = %4ld\n\tZ = %4ld\n", dim.x_dim, dim.y_dim, dim.z_dim);
    return 0;
}