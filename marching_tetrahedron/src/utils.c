#include "utils.h"
#include "struct.h"

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

    verbose_print("Dimensions:\n");
    verbose_print("    x_dim: %zu\n", dim->x_dim);
    verbose_print("    y_dim: %zu\n", dim->y_dim);
    verbose_print("    z_dim: %zu\n", dim->z_dim);

    verbose_print("Origin:\n");
    verbose_print("    x: %f\n", origin[0]);
    verbose_print("    y: %f\n", origin[1]);
    verbose_print("    z: %f\n", origin[2]);

    verbose_print("Cell dimensions:\n");
    verbose_print("    dx: %f\n", dx);
    verbose_print("    dy: %f\n", dy);
    verbose_print("    dz: %f\n", dz);


    *tensor = (dim_t*)malloc(sizeof(dim_t)*dim->x_dim*dim->y_dim*dim->z_dim);
    
    for (int i=0; i<dim->x_dim; i++){
        for (int j=0; j<dim->y_dim; j++){
            for (int k=0; k<dim->z_dim; k++){
                fread(&(*tensor)[k + j*dim->z_dim + i*dim->y_dim*dim->z_dim], sizeof(dim_t), 1, fptr);
                // verbose_print("%f\n", (*tensor)[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
            }
        }
    }

    fclose(fptr);
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
 * @brief Print the stack, with the coordinates and the value
 * 
 * @param start Pointer to pointer to the beginning of the stack
 */
void print_stack(StackNode *start){
    verbose_print("    Stack content:\n");
    while(start != NULL){
        verbose_print("        Coord: (%d, %d, %d); Val: %f\n",     start->coordinate.x,
                                                                    start->coordinate.y,
                                                                    start->coordinate.z,
                                                                    start->owned_value);
        start = start->next;
    }
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
void print_connections(Triangle *triangle, int *count){
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
 * @brief Print on file the triangles and the vertices
 * 
 * This uses the suffix tree data structure to store the vertices and the list for the triangles
 * 
 * @param Polyhedra p is a Pointer to a Polyhedra data structure contatining the first node for triangle and vertices
 * @param char Pointer to the name of the file
 */
void print_on_file(Polyhedra *p, char *name){
    FILE *fptr;
    TriangleCoordNode *curr = p->root;
    VertexNode *del = NULL;
    double first = 0.0;
    double second = 0.0;

    strcat(name, ".pdb");

    fptr = fopen(name, "w");

    print_vertices(p->root, &first, &second, fptr);

    TriangleNode *curr2 = p->triangles;

    // // // Debug print for triangles
    // // printf("Triangles in polyhedron:\n");
    // // TriangleNode *debug_curr = p->triangles;
    // // int triangle_count = 0;

    // // while(debug_curr != NULL){
    // //     printf("Triangle %d: Vertices %d, %d, %d\n", 
    // //            triangle_count++, debug_curr->vert1, debug_curr->vert2, debug_curr->vert3);
    // //     debug_curr = debug_curr->next;
    // // }
    // // printf("Total triangles: %d\n", triangle_count);

    // // if(curr2->vert1 <= 9999 && curr2->vert2 <= 9999 && curr2->vert3 <= 9999){

    while(curr2 != NULL){
        if (curr2->vert1 < curr2->vert2) {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert1->index, curr2->vert2->index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert2->index, curr2->vert1->index);
        }
        
        if (curr2->vert2 < curr2->vert3) {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert2->index, curr2->vert3->index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert3->index, curr2->vert2->index);
        }
        
        if (curr2->vert3 < curr2->vert1) {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert3->index, curr2->vert1->index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr2->vert1->index, curr2->vert3->index);
        }
        curr2 = curr2->next;
    }
    // }
    
    fclose(fptr);
}

void print_for_stats(Polyhedra *p){

    FILE *fptr;
    fptr = fopen("stats.csv", "w");

    TriangleNode *curr2 = p->triangles;

    while(curr2 != NULL){
        if (curr2->vert1 < curr2->vert2) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert1, curr2->vert2);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert2, curr2->vert1);
        }
        
        if (curr2->vert2 < curr2->vert3) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert2, curr2->vert3);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert3, curr2->vert2);
        }
        
        if (curr2->vert3 < curr2->vert1) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert3, curr2->vert1);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert1, curr2->vert3);
        }
        curr2 = curr2->next;
    }

    fclose(fptr);

}

void print_on_separate_files(Polyhedra *p, char *name, int num_triangles){

    int *idxs = print_atoms_separated(p->triangles, name, num_triangles);
    print_connections_separated(p->triangles, name, idxs);
}

int *print_atoms_separated(TriangleNode *curr, char *name, int num_traingles){

    FILE *fptr;
    char file_name[100];
    int N = 50000;
    int min = 0;
    int count = 0; 
    int div = 0;
    int file_number = 0;
    int *offset = (int *)malloc((num_traingles/N + 1) * sizeof(int));

    if (offset == NULL) {
        fprintf(stderr, "Memory allocation failed for offset array\n");
        exit(-1);
    }

    printf("number of offetes: %d\n", num_traingles/N + 1);

    typedef struct print_list{
        struct print_list *next;
        int idx;
    } print_list;

    TriangleNode *counter = curr;

    while(counter != NULL){

        if((count+1)%N == 0 || counter->next == NULL){
            offset[(int)(count/N)] = min;
            printf("count: %d\n", count);
            min = 0;
        }

        if(count / N != 0){
            if (min == 0 || counter->vert1->index < min) {
            min = counter->vert1->index;
            }
            if (counter->vert2->index < min) {
            min = counter->vert2->index;
            }
            if (counter->vert3->index < min) {
            min = counter->vert3->index;
            }
        }

        count++;
        counter = counter->next;
    }

    for (int i = 0; i < 2; i++) {
        printf("Offset[%d]: %d\n", i, offset[i]);
    }

    count = 0;
    print_list *start = NULL;

    while(curr != NULL){

        if(count%N == 0){
            strcpy(file_name, name);
            file_number = count/N;
            div = count;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "w");
            start = NULL;
        }

        div = offset[count/N];
        
        char str[500];
        print_list *temp = start;
        int found1 = 0, found2 = 0, found3 = 0;

        while (temp != NULL) {
            if (temp->idx == curr->vert1->index) found1 = 1;
            if (temp->idx == curr->vert2->index) found2 = 1;
            if (temp->idx == curr->vert3->index) found3 = 1;
            temp = temp->next;
        }

        if (!found1) {
            print_list *new_node = (print_list *)malloc(sizeof(print_list));
            new_node->idx = curr->vert1->index;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert1->index-div, curr->vert1->coordinate1, curr->vert1->coordinate2, curr->vert1->coordinate3);
            fprintf(fptr, "%s", str);
        }

        if (!found2) {
            print_list *new_node = (print_list *)malloc(sizeof(print_list));
            new_node->idx = curr->vert2->index;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert2->index-div, curr->vert2->coordinate1, curr->vert2->coordinate2, curr->vert2->coordinate3);
            fprintf(fptr, "%s", str);
        }

        if (!found3) {
            print_list *new_node = (print_list *)malloc(sizeof(print_list));
            new_node->idx = curr->vert3->index;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert3->index-div, curr->vert3->coordinate1, curr->vert3->coordinate2, curr->vert3->coordinate3);
            fprintf(fptr, "%s", str);
        }

        count++;
        curr = curr->next;
        if(count%N == 0 || curr == NULL){
            // Free the dynamically allocated start list
            while (start != NULL) {
                print_list *temp = start;
                start = start->next;
                free(temp);
            }
            fclose(fptr);
        }
    }

    printf("Counts in the end: %d\n", count);
    return offset;
}

void print_connections_separated(TriangleNode *curr, char *name, int *offsets){
    int N = 50000;
    char file_name[100];

    int count = 0; 

    FILE *fptr;

    int div = 0;

    int file_number = 0;

    while(curr != NULL){

        if(count%N == 0){
            strcpy(file_name, name);
            file_number = count/N;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "a");
        }

        div = offsets[count/N];

        if (curr->vert1->index-div >= 10000 || curr->vert2->index-div >= 10000 || curr->vert3->index-div >= 10000) {
            count++;
            curr = curr->next;
            continue;
        }
        
        if (curr->vert1->index < curr->vert2->index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert1->index-div, curr->vert2->index-div);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert2->index-div, curr->vert1->index-div);
        }

        if (curr->vert2->index < curr->vert3->index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert2->index-div, curr->vert3->index-div);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert3->index-div, curr->vert2->index-div);
        }

        if (curr->vert3->index < curr->vert1->index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert3->index-div, curr->vert1->index-div);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert1->index-div, curr->vert3->index-div);
        }

        count++;
        curr = curr->next;
        if(count%N == 0 || curr == NULL){
            fclose(fptr);
        }
    }
    printf("Counts in the end: %d\n", count);
}