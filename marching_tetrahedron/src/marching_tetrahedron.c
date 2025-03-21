#include "marching_tetrahedron.h"

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor){

    FILE *fptr;
    fptr = fopen(file_name, "rb"); 

    if(fptr == NULL) {
        fprintf(stderr, "Not able to open the file\n");
        exit(-1);
    }

    double dx;
    double dy;
    double dz;
    double origin_x;
    double origin_y;
    double origin_z;
    
    fread(&(dx), sizeof(double), 1, fptr);
    fread(&(dy), sizeof(double), 1, fptr);
    fread(&(dz), sizeof(double), 1, fptr);
    fread(&(origin_x), sizeof(double), 1, fptr);
    fread(&(origin_y), sizeof(double), 1, fptr);
    fread(&(origin_z), sizeof(double), 1, fptr);
    fread(&(dim->x_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->y_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->z_dim), sizeof(size_t), 1, fptr);

    verbose_print("x_dim: %d\n", dim->x_dim);
    verbose_print("y_dim: %d\n", dim->y_dim);
    verbose_print("z_dim: %d\n", dim->z_dim);

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

void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count){

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
                    int action_value = get_action_value(stack);
                    
                    // Print for debug
                    verbose_print("        action val: %d\n", action_value);
                    
                    // get the pairs
                    int *pairs = get_pairs(action_value);
                    
                    // Print pairs for debug
                    if(action_value!=0){
                        verbose_print("    Pairs: ");
                        for (int p = 0; p < (action_value == 7 ? 12 : 6); p += 2) {
                            verbose_print("(%d, %d) ", pairs[p], pairs[p + 1]);
                        }
                        verbose_print("\n");
                        // build the triangle 
                        Triangle *triangle = make_triangle(stack, pairs);
                        
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

                        print_to_file(triangle, count);
                    }
                    
                    

                    free(pairs);
                    free_stack(&stack);
                    free(coordinates);
                }
                    
            }
        }
    }
        
        
        

}

void print_grid(const Dimensions *dim, const dim_t *grid){
    for (int k=0; k<dim->z_dim; k++){
        for (int j=0; j<dim->y_dim; j++){
            for (int i=0; i<dim->x_dim; i++){
                verbose_print("%f\n", grid[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim]);
            }
        }
    }
}

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

int get_action_value(StackNode *start){
    int val[3] = {0,0,0};

    while(start != NULL){
        if(start->owned_value < 0){
            val[0]++;
        }
        if(start->owned_value == 0){
            val[1]++;
        }
        if(start->owned_value > 0){
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

Triangle *make_triangle(StackNode *stack, int *pairs){

    Triangle *triangle = (Triangle *)malloc(sizeof(Triangle));

    verbose_print("TRIANGLE:\n");

    triangle->v1 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v2 = (TriangleVertex *)malloc(sizeof(TriangleVertex));
    triangle->v3 = (TriangleVertex *)malloc(sizeof(TriangleVertex));

    for(int idx=0; idx<6; idx+=2){ // TODO: MANAGE CASE OF 6 PAIRS (12 NUMBERS)
        CubeVertex *point1 = get_coordinate_by_idx(stack, pairs[idx]-1);
        CubeVertex *point2 = get_coordinate_by_idx(stack, pairs[idx+1]-1);

        if(idx==0){
            triangle->v1->x = (dim_t)(point2->x+point1->x)/2;
            triangle->v1->y = (dim_t)(point2->y+point1->y)/2;
            triangle->v1->z = (dim_t)(point2->z+point1->z)/2;
        }
        if(idx==2){
            triangle->v2->x = (dim_t)(point2->x+point1->x)/2;
            triangle->v2->y = (dim_t)(point2->y+point1->y)/2;
            triangle->v2->z = (dim_t)(point2->z+point1->z)/2;
        }
        if(idx==4){
            triangle->v3->x = (dim_t)(point2->x+point1->x)/2;
            triangle->v3->y = (dim_t)(point2->y+point1->y)/2;
            triangle->v3->z = (dim_t)(point2->z+point1->z)/2;
        }
    }

    return triangle;
}

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

void print_to_file(Triangle *triangle, int *count){
    FILE *fptr;

    // ATOM      1 0    PSE A   0      60.000  58.000  46.500  1.00  1.00           C  
    char str[500];
    snprintf(str, sizeof(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            *count*3+0, triangle->v1->x, triangle->v1->y, triangle->v1->z);
    snprintf(str + strlen(str), sizeof(str) - strlen(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            *count*3+1, triangle->v2->x, triangle->v2->y, triangle->v2->z);
    snprintf(str + strlen(str), sizeof(str) - strlen(str), "ATOM  %5d 0    PSE A   0      %6.3f  %6.3f  %6.3f  1.00  1.00           C\n", 
            *count*3+2, triangle->v3->x, triangle->v3->y, triangle->v3->z);

    if (*count == 1) {
        fptr = fopen("write.pdb", "w");
    } else {
        fptr = fopen("write.pdb", "a");
    }
    fprintf(fptr, "%s", str);
    fclose(fptr);
}