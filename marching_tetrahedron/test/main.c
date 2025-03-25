#include "global.h"
#include "triangles.h"
#include "utils.h"
#include "marching_tetrahedron.h"

int main(){
    Dimensions dim;

    char folder_name[100] = "/home/elia/tesi/code/marching_tetrahedron/test/data/";
    char *path = strcat(folder_name, "h_2.bin");

    dim_t threshold = 0.25;
    dim_t *grid;

    double origin[3];

    int cube_decomposition[20] = {4,6,7,8,1,5,6,7,1,3,4,7,1,2,4,6,1,4,6,7};

    read_file(path, &dim, &grid, origin);

    // verbose_call(print_grid(&dim, grid));
    // normalize_grid(&dim, &grid, threshold);
    // verbose_call(print_grid(&dim, grid));

    void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t);

    func_ptr = &linear_interpol;

    Polyhedra p;
    p.triangles = NULL;
    p.vertices = NULL;


    int count = 0;
    marching_tetrahedra(&dim, &grid, cube_decomposition, &count, threshold, origin, func_ptr, &p);

    merge_files("write.pdb", "conn.pdb");

    printf("Count: %d\n", count);

    if(p.triangles == NULL || p.vertices == NULL){
        fprintf(stderr, "either triangles or vertice null\n");
        exit(-1);
    }

    print_with_unique_indices(&p);

    // Print dimensions
    printf("Dimensions: X=%d, Y=%d, Z=%d\n", dim.x_dim, dim.y_dim, dim.z_dim);



    return 0;
}