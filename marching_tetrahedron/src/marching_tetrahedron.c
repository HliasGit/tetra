#include "marching_tetrahedron.h"

/**
 * @brief Reads a file and initializes the tensor data (linearized) structure along with its dimensions and origin.
 *
 * @param file_name The path to the file to be read.
 * @param dim Pointer to a Dimensions structure where the dimensions of the tensor will be stored.
 * @param tensor Pointer to a pointer of type dim_t where the tensor data will be allocated and stored.
 * @param origin Pointer to a double array where the origin coordinates will be stored.
 */
void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, double *origin){

    FILE *fptr;
    fptr = fopen(file_name, "rb"); 

    if(fptr == NULL) {
        fprintf(stderr, "Not able to open the file\n");
        exit(-1);
    }

    double dx;
    double dy;
    double dz;
    
    fread(&(dx), sizeof(double), 1, fptr);
    fread(&(dy), sizeof(double), 1, fptr);
    fread(&(dz), sizeof(double), 1, fptr);
    fread(&(origin[0]), sizeof(double), 1, fptr);
    fread(&(origin[1]), sizeof(double), 1, fptr);
    fread(&(origin[2]), sizeof(double), 1, fptr);
    fread(&(dim->x_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->y_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->z_dim), sizeof(size_t), 1, fptr);

    verbose_print("x_dim: %d\n", dim->x_dim);
    verbose_print("y_dim: %d\n", dim->y_dim);
    verbose_print("z_dim: %d\n", dim->z_dim);

    printf("Origin:\n");
    printf("    x: %f\n", origin[0]);
    printf("    y: %f\n", origin[1]);
    printf("    z: %f\n", origin[2]);

    *tensor = (dim_t*)malloc(sizeof(dim_t)*dim->x_dim*dim->y_dim*dim->z_dim);
    
    for (int i=0; i<dim->x_dim; i++){
        for (int j=0; j<dim->y_dim; j++){
            for (int k=0; k<dim->z_dim; k++){
                fread(&(*tensor)[k + j*dim->y_dim + i*dim->y_dim*dim->z_dim], sizeof(dim_t), 1, fptr);
                // verbose_print("%f\n", (*tensor)[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
            }
        }
    }
    fclose(fptr);
}

/**
 * @brief Function that modifies the grid values subtracting the threshold
 * 
 * @param dim Pointer to the Dimenstion data structure
 * @param grid Pointer to pointer of the grid to be modified
 * @param threshold Isosurface value that we want to evaluate
 */
void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold){

    for(size_t k=0; k<dim->z_dim; k++){
        for (size_t j=0; j<dim->y_dim; j++){
            for (size_t i=0; i<dim->x_dim; i++){
                // verbose_print("val: %f\n", *grid[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
                (*grid)[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim] -= threshold;
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
 * @param count Pointer for a counter of triangles generated
 * @param threshold Value of the threshold for the isosurface
 * @param origin Pointer to the origin coordinates 
 * @param func_ptr Function pointer to invoke dynamically the interpolation function
 */
void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count, dim_t threshold, double *origin, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t)){
    CubeVertex *coordinates;
    StackNode *stack = NULL;
    // for every cube in the space
    for(size_t k=0; k<dim->z_dim-1; k++){ // Cube global coordinate
        for (size_t j=0; j<dim->y_dim-1; j++){
            for (size_t i=0; i<dim->x_dim-1; i++){

                // check if every vertex in the cube is F(x,y,x) - C < threshold and in case skip it 
                for (int tetra = 0; tetra<20; tetra+=4){ // for every tetrahedron in a cube
                    coordinates = malloc(4*sizeof(CubeVertex));
                    verbose_print("Tetra: %d\n", tetra/4+1);
                    for (int idx=tetra; idx<tetra+4; idx ++){ // for every point in a tetrahedra
                        int point = cube_decomposition[idx];
                        // verbose_print("point coord (%d,%d,%d): %d\n",
                        //             i, j, k, point);

                        // find the global tetrahedra coordinates.
                        find_coordinates(idx-tetra, point, i, j, k, &coordinates);
                        verbose_print("    Point: %d\n", point);
                        verbose_print("        coord x: %d\n", coordinates[idx-tetra].x);
                        verbose_print("        coord y: %d\n", coordinates[idx-tetra].y);
                        verbose_print("        coord z: %d\n", coordinates[idx-tetra].z);

                        // get the # neg, # zero, # pos for each tetrahedron in the stack
                        push_into_stack(&stack,
                                        (*grid)[coordinates[idx-tetra].x + 
                                                coordinates[idx-tetra].y*dim->x_dim + 
                                                coordinates[idx-tetra].z*dim->x_dim*dim->y_dim],
                                        coordinates[idx-tetra]);

                        

                    }
                    
                    // Print for debug
                    verbose_call(print_stack(stack));

                    // get the action value
                    int action_value = get_action_value(stack, threshold);
                    
                    // Print for debug
                    verbose_print("        action val: %d\n", action_value);
                    
                    // get the pairs
                    int *pairs = get_pairs(action_value);

                    if(*count>3300)
                        exit(-1);
                    
                    // Print pairs for debug
                    if(action_value!=0){
                        verbose_print("    Pairs: ");
                        for (int p = 0; p < (action_value == 7 ? 12 : 6); p += 2) {
                            verbose_print("(%d, %d) ", pairs[p], pairs[p + 1]);
                        }
                        verbose_print("\n");
                        // build the triangle 
                        Triangle *triangle = make_triangle(stack, pairs, false, threshold, func_ptr);

                        // exit(-1);
                        
                        (*count)++;
                        printf("Triangle #%d\n", *count);
                        printf("    vertex 1:\n");
                        printf("        x: %f\n", triangle->v1->x);
                        printf("        y: %f\n", triangle->v1->y);
                        printf("        z: %f\n", triangle->v1->z);
                        printf("    vertex 2:\n");
                        printf("        x: %f\n", triangle->v2->x);
                        printf("        y: %f\n", triangle->v2->y);
                        printf("        z: %f\n", triangle->v2->z);
                        printf("    vertex 3:\n");
                        printf("        x: %f\n", triangle->v3->x);
                        printf("        y: %f\n", triangle->v3->y);
                        printf("        z: %f\n", triangle->v3->z);

                        print_to_file(triangle, count, origin);
                        print_connections(triangle, count);

                        if(action_value==7 ? true:false){
                            Triangle *triangle = make_triangle(stack, pairs, action_value==7 ? true:false, threshold, func_ptr);
                        
                            (*count)++;
                            printf("Triangle #%d\n", *count);
                            printf("    vertex 1:\n");
                            printf("        x: %f\n", triangle->v1->x);
                            printf("        y: %f\n", triangle->v1->y);
                            printf("        z: %f\n", triangle->v1->z);
                            printf("    vertex 2:\n");
                            printf("        x: %f\n", triangle->v2->x);
                            printf("        y: %f\n", triangle->v2->y);
                            printf("        z: %f\n", triangle->v2->z);
                            printf("    vertex 3:\n");
                            printf("        x: %f\n", triangle->v3->x);
                            printf("        y: %f\n", triangle->v3->y);
                            printf("        z: %f\n", triangle->v3->z);

                            print_to_file(triangle, count, origin);
                            print_connections(triangle, count);
                        }
                    }

                    // printf("COUNT: %d\n", *count);
                    // exit(1);

                    free(pairs);
                    free_stack(&stack);
                    free(coordinates);
                }
                    
            }
        }
    }
        
        
        

}

/**
 * @brief Print the grid values
 * 
 * @param dim Pointer to the dimensions structure
 * @param grid Pointer to the grid values
 */
void print_grid(const Dimensions *dim, const dim_t *grid){
    for (int k=0; k<dim->z_dim; k++){
        for (int j=0; j<dim->y_dim; j++){
            for (int i=0; i<dim->x_dim; i++){
                verbose_print("%f\n", grid[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim]);
            }
        }
    }
}

/**
/**
 * @brief Finds the coordinates of a point in the marching tetrahedron algorithm.
 * 
 * This function calculates the coordinates of a point within a tetrahedron based on 
 * its index and the global cube coordinates. It determines the position of the point 
 * relative to the cube vertices and assigns the corresponding coordinates.
 * 
 * @param idx The index of the point within the tetrahedron (0 to 3).
 * @param point The vertex index of the cube (1 to 8).
 * @param i The x-coordinate of the cube in the grid.
 * @param j The y-coordinate of the cube in the grid.
 * @param k The z-coordinate of the cube in the grid.
 * @param coordinates Pointer to an array of CubeVertex structures where the calculated 
 *                    coordinates will be stored.
 */
void find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k, CubeVertex **coordinates){
    if (point<1 || point>8){
        fprintf(stderr, "Point can't exceed 1-8");
    }

    coord_t one_apex[3];
    coord_t two_apex[3];

    if (i%2 == 0){
        one_apex[0] = i;
    } else {
        one_apex[0] = i+1;
    }

    if (j%2 == 0){
        one_apex[1] = j;
    } else {
        one_apex[1] = j+1;
    }

    if (k%2 == 0){
        one_apex[2] = k;
    } else {
        one_apex[2] = k+1;
    }

    two_apex[0] = 2*i+1-one_apex[0];
    two_apex[1] = 2*j+1-one_apex[1];
    two_apex[2] = 2*k+1-one_apex[2];

    switch (point-1) {
        case 0:
            (*coordinates)[idx].x = one_apex[0];
            (*coordinates)[idx].y = one_apex[1];
            (*coordinates)[idx].z = one_apex[2];
            break;
        case 1:
            (*coordinates)[idx].x = two_apex[0];
            (*coordinates)[idx].y = one_apex[1];
            (*coordinates)[idx].z = one_apex[2];
            break;
        case 2:
            (*coordinates)[idx].x = one_apex[0];
            (*coordinates)[idx].y = two_apex[1];
            (*coordinates)[idx].z = one_apex[2];
            break;
        case 3:
            (*coordinates)[idx].x = two_apex[0];
            (*coordinates)[idx].y = two_apex[1];
            (*coordinates)[idx].z = one_apex[2];
            break;
        case 4:
            (*coordinates)[idx].x = one_apex[0];
            (*coordinates)[idx].y = one_apex[1];
            (*coordinates)[idx].z = two_apex[2];
            break;
        case 5:
            (*coordinates)[idx].x = two_apex[0];
            (*coordinates)[idx].y = one_apex[1];
            (*coordinates)[idx].z = two_apex[2];
            break;
        case 6:
            (*coordinates)[idx].x = one_apex[0];
            (*coordinates)[idx].y = two_apex[1];
            (*coordinates)[idx].z = two_apex[2];
            break;
        case 7:
            (*coordinates)[idx].x = two_apex[0];
            (*coordinates)[idx].y = two_apex[1];
            (*coordinates)[idx].z = two_apex[2];
            break;
        default:
            fprintf(stderr, "Invalid point value\n");
            exit(-1);
    }
}

/**
 * @brief Push a new node (vertex) into the stack
 * 
 * @param start Pointer to pointer to the beginning of the stack
 * @param value Value to be added to the stack
 * @param vert Vertex (Cube vertex) to be added to the stack
 */
void push_into_stack(StackNode **start, dim_t value, CubeVertex vert){
    StackNode *new = (StackNode*) malloc(sizeof(StackNode));
    new->owned_value = value;
    new->coordinate = vert;
    new->next = NULL;

    if(*start == NULL || (*start)->owned_value >= value) {
        new->next = *start;
        *start = new;
    } else {
        StackNode *current = *start;
        while(current->next != NULL && current->next->owned_value < value) {
            current = current->next;
        }
        new->next = current->next;
        current->next = new;
    }
}

/**
 * @brief Free the stack that stores the tetrahedra
 * 
 * @param start Pointer to pointer to the beginning of the stack
 */
void free_stack(StackNode **start){
    StackNode *current = *start;
    StackNode *next;

    while(current != NULL){
        next = current->next;
        free(current);
        current = next;
    }

    *start = NULL;
}

/**
 * @brief Print the stack, with the coordinates and the value
 * 
 * @param start Pointer to pointer to the beginning of the stack
 */
void print_stack(StackNode *start){
    verbose_print("    Stack content:\n");
    while(start != NULL){
        verbose_print("        Coord: (%d, %d, %d); Val: %f\n",    start->coordinate.x,
                                                            start->coordinate.y,
                                                            start->coordinate.z,
                                                            start->owned_value);
        start = start->next;
    }
}

/**
 * @brief Get the action value
 * 
 * @param start Pointer to the beginning of the stack
 * @param threshold Pointer to the threshold of the isosurface 
 */
int get_action_value(StackNode *start, dim_t threshold){
    int val[3] = {0,0,0};

    // printf("Threshold %d\n", threshold);
    // exit(1);

    while(start != NULL){
        if(start->owned_value-threshold < 0){
            val[0]++;
        }
        if(start->owned_value-threshold == 0){
            val[1]++;
        }
        if(start->owned_value-threshold > 0){
            val[2]++;
        }
        start = start->next;
    }

    if(val[0]==0 || (val[0] == 2 && val[1] == 2) || (val[0]== 3 && val[1] == 1) || val[0] == 4){
        return 0;
    }

    if(val[0] == 1){
        if(val[1] == 0)
            return 1;
        if(val[1] == 1)
            return 2;
        if(val[1] == 2)
            return 3;
        if(val[1] == 3)
            return 4;
    }

    if(val[0] == 2){
        if(val[1] == 0)
            return 7;
        if(val[1] == 1)
            return 5;
    }

    if(val[0] == 3)
        return 6;
}

/**
 * @brief Get the pairs of how to access the stack given the action value
 * 
 * @param action_val Action value to be considered
 */
int *get_pairs(int action_val){    
    int *res = NULL;
    switch (action_val)
    {
    
    case 0:
        return res;
    case 1:
        res = malloc(6 * sizeof(int));
        res[0] = 1; res[1] = 2; res[2] = 1; res[3] = 3; res[4] = 1; res[5] = 4;
        return res;
    case 2:
        res = malloc(6 * sizeof(int));
        res[0] = 2; res[1] = 2; res[2] = 1; res[3] = 3; res[4] = 1; res[5] = 4;
        return res;
    case 3:
        res = malloc(6 * sizeof(int));
        res[0] = 2; res[1] = 2; res[2] = 3; res[3] = 3; res[4] = 1; res[5] = 4;
        return res;
    case 4:
        res = malloc(6 * sizeof(int));
        res[0] = 2; res[1] = 2; res[2] = 3; res[3] = 3; res[4] = 4; res[5] = 4;
        return res;
    case 5:
        res = malloc(6 * sizeof(int));
        res[0] = 1; res[1] = 4; res[2] = 2; res[3] = 4; res[4] = 3; res[5] = 3;
        return res;
    case 6:
        res = malloc(6 * sizeof(int));
        res[0] = 1; res[1] = 4; res[2] = 2; res[3] = 4; res[4] = 3; res[5] = 4;
        return res;
    case 7:
        res = malloc(12 * sizeof(int));
        res[0] = 1; res[1] = 4; res[2] = 2; res[3] = 4; res[4] = 1; res[5] = 3; res[6] = 2; res[7] = 4; res[8] = 2; res[9] = 3; res[10] = 1; res[11] = 3;
        return res;
    
    default:
        fprintf(stderr, "Error in chooisng the pairs");
        break;
    }
}

/**
 * @brief Returns the pointer to the triangle generated from the interpolation
 * 
 * @param stack Pointer to the stack beginning
 * @param pairs Pointer to the pairs in order to access the stack
 * @param two_triangles Boolean to handle the case of building two triangles
 * @param threshold Threshold of the isosurface
 * @param func_ptr Function pointer to invoke dynamically the function
 */
Triangle *make_triangle(StackNode *stack, int *pairs, bool two_triangles, dim_t threshold, void (*func_ptr)(TriangleVertex*, CubeVertex*, CubeVertex*, dim_t*, dim_t*, dim_t)){

    Triangle *triangle = (Triangle *)malloc(sizeof(Triangle));

    triangle->v1 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v2 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v3 = (TriangleVertex *)malloc(sizeof(TriangleVertex));

    for(int idx=0; idx<6; idx+=2){
        CubeVertex *point1 = get_coordinate_by_idx(stack, pairs[idx]-1);
        CubeVertex *point2 = get_coordinate_by_idx(stack, pairs[idx+1]-1);
        dim_t *val1 = get_value_by_idx(stack, pairs[idx]-1);
        dim_t *val2 = get_value_by_idx(stack, pairs[idx+1]-1);

        if(idx==0){
            func_ptr(triangle->v1, point1, point2, val1, val2, threshold);
        }        
        if(idx==2){
            func_ptr(triangle->v2, point1, point2, val1, val2, threshold);
        }
        if(idx==4){
            func_ptr(triangle->v3, point1, point2, val1, val2, threshold);
        }
    }

    if (two_triangles)
    {
        for(int idx=6; idx<12; idx+=2){
            CubeVertex *point1 = get_coordinate_by_idx(stack, pairs[idx]-1);
            CubeVertex *point2 = get_coordinate_by_idx(stack, pairs[idx+1]-1);
            dim_t *val1 = get_value_by_idx(stack, pairs[idx]-1);
            dim_t *val2 = get_value_by_idx(stack, pairs[idx+1]-1);
    
            if(idx==6){
                func_ptr(triangle->v1, point1, point2, val1, val2, threshold);
            }
            if(idx==8){
                func_ptr(triangle->v2, point1, point2, val1, val2, threshold);
            }
            if(idx==10){
                func_ptr(triangle->v3, point1, point2, val1, val2, threshold);
            }
        }
    }
    return triangle;
}

/**
 * @brief Access the stack to get the coordinates
 * 
 * @param start Pointer to the stack beginning
 * @param idx Index needed to be accessed
 */
CubeVertex *get_coordinate_by_idx(StackNode *start, int idx){
    int i=0;
    StackNode *ptr = start;
    while(i<idx){
        i++;
        if(ptr == NULL){
            fprintf(stderr, "Stack is smaller than idx\n");
            exit(-1);
        }
        verbose_print("Iter for the coordinates\n");
        ptr = ptr->next;
    }

    verbose_print("Found\n");
    return &ptr->coordinate;
}

/**
 * @brief Access the stack to get the value
 * 
 * @param start Pointer to the stack beginning
 * @param idx Index needed to be accessed
 */
dim_t *get_value_by_idx(StackNode *start, int idx){
    int i=0;
    StackNode *ptr = start;
    while(i<idx){
        i++;
        if(ptr == NULL){
            fprintf(stderr, "Stack is smaller than idx\n");
            exit(-1);
        }
        verbose_print("Iter for the coordinates\n");
        ptr = ptr->next;
    }

    verbose_print("Found\n");
    return &ptr->owned_value;
}

/**
 * @brief Print to a file the triangle
 * 
 * @param triangle Pointer to the triangle
 * @param count Pointer to the triangle number
 * @param origin Pointer to the origin of the reference system
 */
void print_to_file(Triangle *triangle, int *count, double *origin){
    FILE *fptr;

    // ATOM      1 0    PSE A   0      60.000  58.000  46.500  1.00  1.00           C  
    char str[500];
    snprintf(str, sizeof(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            (*count)*3+0, (triangle->v1->x), (triangle->v1->y), (triangle->v1->z));
    snprintf(str + strlen(str), sizeof(str) - strlen(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            (*count)*3+1, (triangle->v2->x), (triangle->v2->y), (triangle->v2->z));
    snprintf(str + strlen(str), sizeof(str) - strlen(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            (*count)*3+2, (triangle->v3->x), (triangle->v3->y), (triangle->v3->z));

    if ((*count) == 1) {
        fptr = fopen("write.pdb", "w");
        // printf("QUA\n");
        // exit(1);
    } else {
        fptr = fopen("write.pdb", "a");
    }
    fprintf(fptr, "%s", str);

    
    fclose(fptr);
}

/**
 * @brief Print to a file the triangle connections
 * 
 * @param triangle Pointer to the triangle
 * @param count Pointer to the triangle number
 */
void print_connections(Triangle *triangle, int*count){
    FILE *fptr;

    // ATOM      1 0    PSE A   0      60.000  58.000  46.500  1.00  1.00           C  
    char str[500];
    if ((*count)*3+0 > (*count)*3+1) {
        snprintf(str, sizeof(str), "CONECT%5d%5d\n", 
                (*count)*3+1, (*count)*3+0);
    } else {
        snprintf(str, sizeof(str), "CONECT%5d%5d\n", 
                (*count)*3+0, (*count)*3+1);
    }

    if ((*count)*3+1 > (*count)*3+2) {
        snprintf(str + strlen(str), sizeof(str) - strlen(str), "CONECT%5d%5d\n", 
                (*count)*3+2, (*count)*3+1);
    } else {
        snprintf(str + strlen(str), sizeof(str) - strlen(str), "CONECT%5d%5d\n", 
                (*count)*3+1, (*count)*3+2);
    }

    if ((*count)*3+2 > (*count)*3+0) {
        snprintf(str + strlen(str), sizeof(str) - strlen(str), "CONECT%5d%5d\n", 
                (*count)*3+0, (*count)*3+2);
    } else {
        snprintf(str + strlen(str), sizeof(str) - strlen(str), "CONECT%5d%5d\n", 
                (*count)*3+2, (*count)*3+0);
    }

    if ((*count) == 1) {
        fptr = fopen("conn.pdb", "w");
    } else {
        fptr = fopen("conn.pdb", "a");
    }
    fprintf(fptr, "%s", str);
    fclose(fptr);
}

/**
 * @brief Merges two files into one
 * 
 * @param atoms Pointer to the atoms file
 * @param conn Pointer to the connections file 
 */
void merge_files(char *atoms, char* conn){
    FILE *f_atoms = fopen(atoms, "a");
    FILE *f_conn = fopen(conn, "r");

    if (f_atoms == NULL || f_conn == NULL) {
        fprintf(stderr, "Error opening files for merging\n");
        exit(-1);
    }

    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), f_conn) != NULL) {
        fputs(buffer, f_atoms);
    }

    fclose(f_atoms);
    fclose(f_conn);
}

/**
 * @brief Compute the midpoint interpolation
 * 
 * @param vtx Pointer to the vertex to be computed, to store it
 * @param point1 Pointer to the first CubeVertex
 * @param point1 Pointer to the second CubeVertex
 * @param val1 Pointer to the grid value of the first point. Needed just for the function pointer
 * @param val2 Pointer to the grid value of the second point.Needed just for the function pointer
 * @param threshold Value of the threshold value. Needed just for the function pointer
 */
void midpoint_interpol(TriangleVertex *vtx, CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold){
    vtx->x = (coord_t)(point2->x+point1->x)/2;
    vtx->y = (coord_t)(point2->y+point1->y)/2;
    vtx->z = (coord_t)(point2->z+point1->z)/2;
}

/**
 * @brief Compute the linear interpolation
 * 
 * @param vtx Pointer to the vertex to be computed, to store it
 * @param point1 Pointer to the first CubeVertex
 * @param point1 Pointer to the second CubeVertex
 * @param val1 Pointer to the grid value of the first point
 * @param val2 Pointer to the grid value of the second point
 * @param threshold Value of the threshold value 
 */
void linear_interpol(TriangleVertex *vtx,CubeVertex *point1, CubeVertex *point2, dim_t *val1, dim_t *val2, dim_t threshold){
    
    if ((*val2) - (*val1) == 0) {
        fprintf(stderr, "Error: Division by zero in linear interpolation\n");
        exit(-1);
    }

    vtx->x = ((coord_t)(point1->x) + ((((coord_t)point2->x - (coord_t)point1->x) / ((*val2) - (*val1))) * (threshold - (*val1))));
    vtx->y = ((coord_t)(point1->y) + ((((coord_t)point2->y - (coord_t)point1->y) / ((*val2) - (*val1))) * (threshold - (*val1))));
    vtx->z = ((coord_t)(point1->z) + ((((coord_t)point2->z - (coord_t)point1->z) / ((*val2) - (*val1))) * (threshold - (*val1))));

    
}