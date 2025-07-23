#include "utils.h"


/**
 * @brief Reads a file and initializes the tensor data (linearized) structure along with its dimensions and origin.
 *
 * @param file_name The path to the file to be read.
 * @param dim Pointer to a Dimensions structure where the dimensions of the tensor will be stored.
 * @param tensor Pointer to a pointer of type dim_t where the tensor data will be allocated and stored.
 * @param origin Pointer to a double array where the origin coordinates will be stored.
 */
void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, dim_t *origin){

    FILE *fptr;
    fptr = fopen(file_name, "rb"); 

    if(fptr == NULL) {
        fprintf(stderr, "Not able to open the file\n");
        exit(-1);
    }

    dim_t dx;
    dim_t dy;
    dim_t dz;
    
    fread(&(dx), sizeof(dim_t), 1, fptr);
    fread(&(dy), sizeof(dim_t), 1, fptr);
    fread(&(dz), sizeof(dim_t), 1, fptr);
    fread(&(origin[0]), sizeof(dim_t), 1, fptr);
    fread(&(origin[1]), sizeof(dim_t), 1, fptr);
    fread(&(origin[2]), sizeof(dim_t), 1, fptr);
    fread(&(dim->x_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->y_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->z_dim), sizeof(size_t), 1, fptr);

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
    printf("    Stack content:\n");
    while(start != NULL){
        printf("        Coord: (%d, %d, %d); Val: %f\n",     start->coordinate.x,
                                                                    start->coordinate.y,
                                                                    start->coordinate.z,
                                                                    start->owned_value);
        start = start->next;
    }
}

/**
 * @brief Print the connections to a csv file to get the idxs connection locality over the connections
 * 
 * @param p Pointer to the polyhedra data structure 
 */
void print_for_stats(Polyhedra *p){

    FILE *fptr;
    fptr = fopen("stats.csv", "w");
    fprintf(fptr, "first,second\n");

    TriangleNode *curr2 = p->triangles;

    while(curr2 != NULL){
        if (curr2->vert1 < curr2->vert2) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert1->index, curr2->vert2->index);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert2->index, curr2->vert1->index);
        }
        
        if (curr2->vert2 < curr2->vert3) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert2->index, curr2->vert3->index);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert3->index, curr2->vert2->index);
        }
        
        if (curr2->vert3 < curr2->vert1) {
            fprintf(fptr, "%8d,%8d\n", curr2->vert3->index, curr2->vert1->index);
        } else {
            fprintf(fptr, "%8d,%8d\n", curr2->vert1->index, curr2->vert3->index);
        }
        curr2 = curr2->next;
    }

    fclose(fptr);
}

/**
 * @brief Print the molecules on separate files so that it's possible to load the whole molecule on VMD through PDB
 * 
 * @param p Pointer to the data structure
 * @param molecule_name Ptr to the name of the molecule
 * @param molecule_path Ptr to the path of the molecule
 * @param num_triangles Number of total triangles generated
 */
void print_on_separate_files(Polyhedra *p, char *molecule_name, char *molecule_path, int num_triangles){
    struct stat st;

    char *folder = "../../results/";
    strcat(molecule_path, folder);
    strcat(molecule_path, molecule_name);

    if (stat(molecule_path, &st) == 0 && S_ISDIR(st.st_mode)) {
        char command[256];
        snprintf(command, sizeof(command), "rm -rf %s", molecule_path);
        if (system(command) != 0) {
            fprintf(stderr, "Failed to remove existing folder: %s\n", molecule_path);
            exit(-1);
        }
    }

    if (stat(molecule_path, &st) == -1) {
        mkdir(molecule_path, 0700);
    }

    strcat(molecule_path, "/");

    print_atoms_connections_separated(p->triangles, molecule_name, molecule_path, num_triangles);
    // print_connections_separated(p->triangles, molecule_name, molecule_path, idxs);
}

/**
 * @brief Print the atoms and the connections on the splitted file
 * 
 * The local counter is introduced to have a atom indexing for every file going 0-9999
 * beign usable on PDB (Max idx for connections is 10000)
 * 
 * @param custart_trianglesrr Pointer to the triangle list
 * @param molecule_name Ptr to the name of the molecule
 * @param result_path Ptr to the result path of the molecule
 * @param num_triangles Number of total triangles generated
 *
 * @return Pointer to the min offsets
 */
void print_atoms_connections_separated(TriangleNode *start_triangles, char *molecule_name, char *result_path, int num_traingles){

    FILE *fptr;
    char file_name[200];
    int N = 33333;
    int count = 0; 
    int file_number = 0;

    printf("The result path is: %s\n", result_path);
    printf("\n");

    printf("Number of produced files:   %d\n", num_traingles/N + 1);

    int local_counter_1 = 0, local_counter_2 = 0, local_counter_3 = 0;
    int temp_counter_1 = 0, temp_counter_2 = 0, temp_counter_3 = 0;

    TriangleNode *curr = start_triangles;

    while(curr != NULL){
        if(count%N == 0){
            strcpy(file_name, result_path);
            strcat(file_name, molecule_name);
            file_number = count/N;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "w");
            local_counter_1 = 0;
            local_counter_2 = 0;
            local_counter_3 = 0;
        }
        
        char str[500];

        snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
        local_counter_1, curr->vert1->coordinate1, curr->vert1->coordinate2, curr->vert1->coordinate3);
        fprintf(fptr, "%s", str);
        temp_counter_1 = local_counter_1;
        local_counter_2 = ++local_counter_1;
        
        snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
        local_counter_2, curr->vert2->coordinate1, curr->vert2->coordinate2, curr->vert2->coordinate3);
        fprintf(fptr, "%s", str);
        temp_counter_2 = local_counter_2;
        local_counter_3 = ++local_counter_2;
        
        snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
        local_counter_3, curr->vert3->coordinate1, curr->vert3->coordinate2, curr->vert3->coordinate3);
        fprintf(fptr, "%s", str);
        temp_counter_3 = local_counter_3;
        local_counter_1 = ++local_counter_3;

        fprintf(fptr, "CONECT%5d%5d\n", temp_counter_1, temp_counter_2);
        fprintf(fptr, "CONECT%5d%5d\n", temp_counter_2, temp_counter_3);
        fprintf(fptr, "CONECT%5d%5d\n", temp_counter_1, temp_counter_3);

        count++;
        curr = curr->next;

        if((count+1%N == 0)){
            fclose(fptr);
        }
    }
}

/**
 * @brief print to the console the triangles in the data structure
 * 
 * @param start Pointer to the beginning of the list
 */
void print_to_console_traingles(TriangleNode* start){
    int count = 0;
    printf("Triangles:\n");
    while (start != NULL) {
        count++;
        printf("    Triangle: %d\n", count);
        printf("        Vertex 1: (%f, %f, %f)\n", start->vert1->coordinate1, start->vert1->coordinate2, start->vert1->coordinate3);
        printf("        Vertex 2: (%f, %f, %f)\n", start->vert2->coordinate1, start->vert2->coordinate2, start->vert2->coordinate3);
        printf("        Vertex 3: (%f, %f, %f)\n", start->vert3->coordinate1, start->vert3->coordinate2, start->vert3->coordinate3);
        start = start->next;
    }
}

void print_triangles_cpu(   nonunique_triangle_node *start,
                        char *molecule_name, char *molecule_path){
    struct stat st;

    const char *folder = "../../results/";
    strcat(molecule_path, folder);
    strcat(molecule_path, molecule_name);

    if (stat(molecule_path, &st) == 0 && S_ISDIR(st.st_mode)) {
        char command[256];
        snprintf(command, sizeof(command), "rm -rf %s", molecule_path);
        if (system(command) != 0) {
            fprintf(stderr, "Failed to remove existing folder: %s\n", molecule_path);
            exit(-1);
        }
    }

    if (stat(molecule_path, &st) == -1) {
        mkdir(molecule_path, 0700);
    }

    strcat(molecule_path, "/");

    int file_number = 0;
    int local_counter = 0;
    int empty = 0;

    FILE *pdb_file;
    while (start != NULL) {
        
        if(local_counter == 0){
            char file_name[256];
            strcpy(file_name, molecule_path);
            strcat(file_name, molecule_name);
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            // printf("Writing triangles to file: %s\n", file_name);
            pdb_file = fopen(file_name, "w");

            // printf("Opening file: %s\n", file_name);
            
            if (!pdb_file) {
                fprintf(stderr, "Failed to open points.pdb for writing.\n");
            } 
            file_number++;
        }

        // printf("Triangle %d:\n", local_counter/3);
        // printf("    Vertex 1: (%f, %f, %f)\n", start->tri->v1->x, start->tri->v1->y, start->tri->v1->z);
        // printf("    Vertex 2: (%f, %f, %f)\n", start->tri->v2->x, start->tri->v2->y, start->tri->v2->z);
        // printf("    Vertex 3: (%f, %f, %f)\n", start->tri->v3->x, start->tri->v3->y, start->tri->v3->z);
        
        

        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, start->tri->v1->x, start->tri->v1->y, start->tri->v1->z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, start->tri->v1->, start->tri->v1->y, start->tri->v1->z);
        local_counter++;
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, start->tri->v2->x, start->tri->v2->y, start->tri->v2->z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, start->tri->v1->, start->tri->v1->y, start->tri->v1->z);
        local_counter++;
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, start->tri->v3->x, start->tri->v3->y, start->tri->v3->z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
        local_counter++;
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 3, local_counter - 2);
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 2, local_counter - 1);
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 3, local_counter - 1);

        if(local_counter == 33333 || start->next == NULL){
            local_counter = 0;
            fclose(pdb_file);
        }

        start = start->next;
    }
    printf("Relevant points written\n");


    
}