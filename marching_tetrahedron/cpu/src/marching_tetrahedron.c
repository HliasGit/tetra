#include "marching_tetrahedron.h"
#include <stdbool.h>

/**
 * @brief Function that modifies the grid values subtracting the threshold
 *
 * @param dim Pointer to the Dimenstion data structure
 * @param grid Pointer to pointer of the grid to be modified
 * @param threshold Isosurface value that we want to evaluate
 */
void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold)
{

    for (size_t k = 0; k < dim->z_dim; k++)
    {
        for (size_t j = 0; j < dim->y_dim; j++)
        {
            for (size_t i = 0; i < dim->x_dim; i++)
            {
                // verbose_print("val: %f\n", *grid[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
                (*grid)[i + dim->x_dim * j + k * dim->x_dim * dim->y_dim] -= threshold;
            }
        }
    }
}

/**
 * @brief function that generate the triangles of the final isosurface
 *
 * @param dim Pointer to the Dimenstions structure containing the dimensions
 * @param grid Pointer to pointer to the grid values (cube vertices)
 * @param cube_decomposition Array containing the vertices of the cube decomposition (4*#decompositions elements)
 * @param threshold Value of the threshold for the isosurface
 * @param origin Pointer to the origin coordinates
 * @param func_ptr Function pointer to invoke dynamically the interpolation function
 * @param p Pointer to a Polyhedra data structure containing the beginning of the Triangles and Vertices data
 * @param triangle_counter Pointer to a variable containing the number of triangles created
 * @param vertex_counter Pointer to a variable containing the number of vertices created
 */

void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, dim_t threshold, dim_t *origin,
                         void (*func_ptr)(TriangleVertex *, CubeVertex *, CubeVertex *, dim_t *, dim_t *, dim_t),
                         Polyhedra *p, size_t *triangle_counter, size_t *vertex_counter)
{

    CubeVertex *coordinates;
    StackNode *stack = NULL;
    (*vertex_counter) = 0;
    (*triangle_counter) = 0;

    size_t count_swap = 0;
    size_t re_count_swap = 0;
    bool res;

    size_t progress = 0;
    size_t tot_scan = dim->x_dim * dim->y_dim * dim->z_dim;

    // for every cube in the space
    for (size_t i = 0; i < dim->x_dim - 1; i++)
    {
        for (size_t j = 0; j < dim->y_dim - 1; j++)
        {
            for (size_t k = 0; k < dim->z_dim - 1; k++)
            { // Cube global coordinate

                progress = (i * dim->y_dim * dim->z_dim + j * dim->z_dim + k) * 100 / tot_scan;
                static size_t last_progress = 0;
                if (progress % 10 == 0 && progress != 0 && progress != last_progress)
                {
                    printf("Progress: %zu%%\r", progress); fflush(stdout);
                    last_progress = progress;
                }

                int idx = 0;
                int point = 0;

                // check if every vertex in the cube is F(x,y,x) - C < threshold and in case skip it
                for (int tetra = 0; tetra < 20; tetra += 4)
                { // for every tetrahedron in a cube
                    coordinates = malloc(4 * sizeof(CubeVertex));
                    // printf("Tetra: %d, cube (%d,%d,%d)\n", tetra / 4 + 1, i, j, k);
                    int permutations = 0;
                    
                    for (idx = tetra; idx < tetra + 4; idx++)
                    { // for every point in a tetrahedra

                        point = cube_decomposition[idx];
                        // printf("point coord (%d,%d,%d): %d\n", i, j, k, point);

                        // find the global tetrahedra coordinates.
                        find_coordinates(idx - tetra, point, i, j, k, &coordinates);
                        // printf("find_coordinates result: %d\n", res);
                        verbose_print("    Point: %d\n", point);
                        verbose_print("        coord x: %d\n", coordinates[idx - tetra].x);
                        verbose_print("        coord y: %d\n", coordinates[idx - tetra].y);
                        verbose_print("        coord z: %d\n", coordinates[idx - tetra].z);

                        // get the # neg, # zero, # pos for each tetrahedron in the stack
                        push_into_stack(&stack,
                                        (*grid)[coordinates[idx - tetra].z +
                                                coordinates[idx - tetra].y * dim->z_dim +
                                                coordinates[idx - tetra].x * dim->z_dim * dim->y_dim],
                                        coordinates[idx - tetra]);
                    }

                    // Print the stack for debugging purposes
                    StackNode *current = stack;
                    int stack_index = 0;
                    // printf("Stack contents:\n");
                    while (current != NULL) {
                        // printf("    Value: %f, Coordinates: (%d, %d, %d), point: %d\n",
                        //        current->owned_value,
                        //        current->coordinate.x,
                        //        current->coordinate.y,
                        //        current->coordinate.z,
                        //         current->point);
                        coordinates[stack_index] = current->coordinate;
                        stack_index++;
                        current = current->next;
                    }

                    bool is_positive_orientation = tetrahedron_determinant(coordinates);
                   
                    // get the action value
                    int action_value = get_action_value(stack, threshold);

                    // get the pairs
                    int *pairs = get_pairs(action_value);

                    if (action_value != 0)
                    {
                        // build the triangle
                        Triangle *triangle = make_triangle(stack, pairs, false, threshold, func_ptr, is_positive_orientation);
                        (*triangle_counter)++;
                        
                        push_triangle(&p, triangle, vertex_counter);

                        free(triangle->v1);
                        free(triangle->v2);
                        free(triangle->v3);
                        free(triangle);

                        // build the second triangle in case the tetrahedra has two of them
                        if (action_value == 7 ? true : false)
                        {
                            Triangle *second_triangle = make_triangle(stack, pairs, true, threshold, func_ptr, is_positive_orientation);
                            (*triangle_counter)++;

                            push_triangle(&p, second_triangle, vertex_counter);

                            free(second_triangle->v1);
                            free(second_triangle->v2);
                            free(second_triangle->v3);
                            free(second_triangle);
                        }
                    }

                    free(pairs);
                    free_stack(&stack);
                    free(coordinates);
                }
            }
        }
    }

    // Reverse the triangle list. Since I make it pushing from the head I have to reverse so the head has index 0
    reverse_list(&p->triangles);

    printf("# of triangles: %8zu\n", (*triangle_counter));
    printf("# of vertices:  %8zu\n", (*vertex_counter));
    printf("# to be swapped: %8zu\n", count_swap);
    printf("# to be checked after swap: %8zu\n", re_count_swap);
}

/**
 * @brief Finds the coordinates of a point in the marching tetrahedron algorithm.
 *
 * This function calculates the coordinates of a point within a tetrahedron based on
 * its index and the global cube coordinates. It determin
 * @param idx The index of the point within the tetrahedron (0 to 3).
 * @param point The vertex index of the cube (1 to 8).
 * @param i The x-coordinate of the cube in the grid.
 * @param j The y-coordinate of the cube in the grid.
 * @param k The z-coordinate of the cube in the grid.
 * @param coordinates Pointer to an array of CubeVertex structures where the calculated
 *                    coordinates will be stored.
 */
bool find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k,
                      CubeVertex **coordinates)
{
    if (point < 1 || point > 8)
    {
        fprintf(stderr, "Point can't exceed 1-8");
    }

    coord_t one_apex[3];
    coord_t two_apex[3];

    if (i % 2 == 0)
    {
        one_apex[0] = i;
    }
    else
    {
        one_apex[0] = i + 1;
    }

    if (j % 2 == 0)
    {
        one_apex[1] = j;
    }
    else
    {
        one_apex[1] = j + 1;
    }

    if (k % 2 == 0)
    {
        one_apex[2] = k;
    }
    else
    {
        one_apex[2] = k + 1;
    }

    two_apex[0] = 2 * i + 1 - one_apex[0];
    two_apex[1] = 2 * j + 1 - one_apex[1];
    two_apex[2] = 2 * k + 1 - one_apex[2];

    switch (point)
    {
    case 1:
        (*coordinates)[idx].x = one_apex[0];
        (*coordinates)[idx].y = one_apex[1];
        (*coordinates)[idx].z = one_apex[2];
        break;
    case 2:
        (*coordinates)[idx].x = two_apex[0];
        (*coordinates)[idx].y = one_apex[1];
        (*coordinates)[idx].z = one_apex[2];
        break;
    case 3:
        (*coordinates)[idx].x = one_apex[0];
        (*coordinates)[idx].y = two_apex[1];
        (*coordinates)[idx].z = one_apex[2];
        break;
    case 4:
        (*coordinates)[idx].x = two_apex[0];
        (*coordinates)[idx].y = two_apex[1];
        (*coordinates)[idx].z = one_apex[2];
        break;
    case 5:
        (*coordinates)[idx].x = one_apex[0];
        (*coordinates)[idx].y = one_apex[1];
        (*coordinates)[idx].z = two_apex[2];
        break;
    case 6:
        (*coordinates)[idx].x = two_apex[0];
        (*coordinates)[idx].y = one_apex[1];
        (*coordinates)[idx].z = two_apex[2];
        break;
    case 7:
        (*coordinates)[idx].x = one_apex[0];
        (*coordinates)[idx].y = two_apex[1];
        (*coordinates)[idx].z = two_apex[2];
        break;
    case 8:
        (*coordinates)[idx].x = two_apex[0];
        (*coordinates)[idx].y = two_apex[1];
        (*coordinates)[idx].z = two_apex[2];
        break;
    default:
        fprintf(stderr, "Invalid point value\n");
        exit(-1);
    }

    bool res = true;
    return res;
}

/**
 * @brief Get the action value
 *
 * @param start Pointer to the beginning of the stack
 * @param threshold Pointer to the threshold of the isosurface
 */
int get_action_value(StackNode *start, dim_t threshold)
{
    int val[3] = {0, 0, 0};

    // printf("Threshold %f\n", threshold);
    // exit(1);

    double epsilon = 1e-10;

    while (start != NULL)
    {
        if (start->owned_value - threshold < epsilon)
        {
            val[0]++;
        }
        if (start->owned_value - threshold == epsilon)
        {
            val[1]++;
        }
        if (start->owned_value - threshold > epsilon)
        {
            val[2]++;
        }
        start = start->next;
    }

    if (val[0] == 0 || (val[0] == 2 && val[1] == 2) || (val[0] == 3 && val[1] == 1) || val[0] == 4)
    {
        return 0;
    }

    if (val[0] == 1)
    {
        if (val[1] == 0)
            return 1;
        if (val[1] == 1)
            return 2;
        if (val[1] == 2)
            return 3;
        if (val[1] == 3)
            return 4;
    }

    if (val[0] == 2)
    {
        if (val[1] == 0)
            return 7;
        if (val[1] == 1)
            return 5;
    }

    if (val[0] == 3)
        return 6;
}

/**
 * @brief Get the pairs of how to access the stack given the action value
 *
 * @param action_val Action value to be considered
 */
int *get_pairs(int action_val)
{
    int *res = NULL;
    switch (action_val)
    {

    case 0:
        return res;
    case 1:
        res = malloc(6 * sizeof(int));
        res[0] = 1;
        res[1] = 2;
        res[2] = 1;
        res[3] = 3;
        res[4] = 1;
        res[5] = 4;
        return res;
    case 2:
        res = malloc(6 * sizeof(int));
        res[0] = 2;
        res[1] = 2;
        res[2] = 1;
        res[3] = 3;
        res[4] = 1;
        res[5] = 4;
        return res;
    case 3:
        res = malloc(6 * sizeof(int));
        res[0] = 2;
        res[1] = 2;
        res[2] = 3;
        res[3] = 3;
        res[4] = 1;
        res[5] = 4;
        return res;
    case 4:
        res = malloc(6 * sizeof(int));
        res[0] = 2;
        res[1] = 2;
        res[2] = 3;
        res[3] = 3;
        res[4] = 4;
        res[5] = 4;
        return res;
    case 5:
        res = malloc(6 * sizeof(int));
        res[0] = 1;
        res[1] = 4;
        res[2] = 2;
        res[3] = 4;
        res[4] = 3;
        res[5] = 3;
        return res;
    case 6:
        res = malloc(6 * sizeof(int));
        res[0] = 1;
        res[1] = 4;
        res[2] = 2;
        res[3] = 4;
        res[4] = 3;
        res[5] = 4;
        return res;
    case 7:
        res = malloc(12 * sizeof(int));
        res[0] = 1;
        res[1] = 4;
        res[2] = 2;
        res[3] = 4;
        res[4] = 1;
        res[5] = 3;
        res[6] = 2;
        res[7] = 4;
        res[8] = 2;
        res[9] = 3;
        res[10] = 1;
        res[11] = 3;
        return res;

    default:
        fprintf(stderr, "Error in chooisng the pairs");
        break;
    }
}

bool tetrahedron_determinant(CubeVertex *coords)
{
    double mat1[3][3] = {
        {coords[1].x, coords[2].x, coords[3].x},
        {coords[1].y, coords[2].y, coords[3].y},
        {coords[1].z, coords[2].z, coords[3].z}};

    double mat2[3][3] = {
        {coords[0].x, coords[2].x, coords[3].x},
        {coords[0].y, coords[2].y, coords[3].y},
        {coords[0].z, coords[2].z, coords[3].z}};

    double mat3[3][3] = {
        {coords[0].x, coords[1].x, coords[3].x},
        {coords[0].y, coords[1].y, coords[3].y},
        {coords[0].z, coords[1].z, coords[3].z}};

    double mat4[3][3] = {
        {coords[0].x, coords[1].x, coords[2].x},
        {coords[0].y, coords[1].y, coords[2].y},
        {coords[0].z, coords[1].z, coords[2].z}};

    double det =
        three_det(mat1) - three_det(mat2) + three_det(mat3) - three_det(mat4);

    return det >= 0;
}

double three_det(double mat[3][3])
{
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
           mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
           mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

void marching_tetrahedra_list(Dimensions *dim, dim_t **grid, int *cube_decomposition, dim_t threshold, dim_t *origin,
                            void (*func_ptr)(TriangleVertex *, CubeVertex *, CubeVertex *, dim_t *, dim_t *, dim_t),
                            Polyhedra *p, size_t *triangle_counter, size_t *vertex_counter, int size,
                            int *results)
{

    CubeVertex *coordinates;
    StackNode *stack = NULL;
    (*vertex_counter) = 0;
    (*triangle_counter) = 0;

    size_t count_swap = 0;
    size_t re_count_swap = 0;
    bool res;

    size_t progress = 0;
    size_t tot_scan = dim->x_dim * dim->y_dim * dim->z_dim;

    // for every cube in the space
    for(int cube = 0; cube<size; cube++)
    { // Cube global coordinate

        int idx_cube = results[cube];
        int DZ = dim->z_dim;
        int DY = dim->y_dim;
        int k = idx_cube % DZ;
        int j = (idx_cube / DZ) % DY;
        int i = idx_cube / (DZ * DY);

        progress = (i * dim->y_dim * dim->z_dim + j * dim->z_dim + k) * 100 / tot_scan;
        static size_t last_progress = 0;
        if (progress % 10 == 0 && progress != 0 && progress != last_progress)
        {
            printf("Progress: %zu%%\r", progress); fflush(stdout);
            last_progress = progress;
        }

        int idx = 0;
        int point = 0;

        // check if every vertex in the cube is F(x,y,x) - C < threshold and in case skip it
        for (int tetra = 0; tetra < 20; tetra += 4)
        { // for every tetrahedron in a cube
            coordinates = malloc(4 * sizeof(CubeVertex));
            // printf("Tetra: %d, cube (%d,%d,%d)\n", tetra / 4 + 1, i, j, k);
            int permutations = 0;
            
            for (idx = tetra; idx < tetra + 4; idx++)
            { // for every point in a tetrahedra

                point = cube_decomposition[idx];
                // printf("point coord (%d,%d,%d): %d\n", i, j, k, point);

                // find the global tetrahedra coordinates.
                find_coordinates(idx - tetra, point, i, j, k, &coordinates);
                // printf("find_coordinates result: %d\n", res);
                verbose_print("    Point: %d\n", point);
                verbose_print("        coord x: %d\n", coordinates[idx - tetra].x);
                verbose_print("        coord y: %d\n", coordinates[idx - tetra].y);
                verbose_print("        coord z: %d\n", coordinates[idx - tetra].z);

                // get the # neg, # zero, # pos for each tetrahedron in the stack
                push_into_stack(&stack,
                                (*grid)[coordinates[idx - tetra].z +
                                        coordinates[idx - tetra].y * dim->z_dim +
                                        coordinates[idx - tetra].x * dim->z_dim * dim->y_dim],
                                coordinates[idx - tetra]);
            }

            // Print the stack for debugging purposes
            StackNode *current = stack;
            int stack_index = 0;
            // printf("Stack contents:\n");
            while (current != NULL) {
                // printf("    Value: %f, Coordinates: (%d, %d, %d), point: %d\n",
                //        current->owned_value,
                //        current->coordinate.x,
                //        current->coordinate.y,
                //        current->coordinate.z,
                //         current->point);
                coordinates[stack_index] = current->coordinate;
                stack_index++;
                current = current->next;
            }

            bool is_positive_orientation = tetrahedron_determinant(coordinates);
            
            // get the action value
            int action_value = get_action_value(stack, threshold);

            // get the pairs
            int *pairs = get_pairs(action_value);

            if (action_value != 0)
            {
                // build the triangle
                Triangle *triangle = make_triangle(stack, pairs, false, threshold, func_ptr, is_positive_orientation);
                (*triangle_counter)++;
                
                push_triangle(&p, triangle, vertex_counter);

                free(triangle->v1);
                free(triangle->v2);
                free(triangle->v3);
                free(triangle);

                // build the second triangle in case the tetrahedra has two of them
                if (action_value == 7 ? true : false)
                {
                    Triangle *second_triangle = make_triangle(stack, pairs, true, threshold, func_ptr, is_positive_orientation);
                    (*triangle_counter)++;

                    push_triangle(&p, second_triangle, vertex_counter);

                    free(second_triangle->v1);
                    free(second_triangle->v2);
                    free(second_triangle->v3);
                    free(second_triangle);
                }
            }

            free(pairs);
            free_stack(&stack);
            free(coordinates);
        }
    }
    // Reverse the triangle list. Since I make it pushing from the head I have to reverse so the head has index 0
    reverse_list(&p->triangles);

    printf("# of triangles: %8zu\n", (*triangle_counter));
    printf("# of vertices:  %8zu\n", (*vertex_counter));
    printf("# to be swapped: %8zu\n", count_swap);
    printf("# to be checked after swap: %8zu\n", re_count_swap);
}