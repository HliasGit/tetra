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
 * @param curr Pointer to the triangle list
 * @param molecule_name Ptr to the name of the molecule
 * @param result_path Ptr to the result path of the molecule
 * @param num_triangles Number of total triangles generated
 *
 * @return Pointer to the min offsets
 */
void print_atoms_connections_separated(TriangleNode *start_triangles, char *molecule_name, char *result_path, int num_traingles){

    FILE *fptr;
    char file_name[200];
    int N = 10000;          // N can be increased somehow
    int min = 0;
    int count = 0; 
    int div = 0;
    int file_number = 0;
    int missed = 0;

    printf(result_path);
    printf("\n");

    printf("Number of files: %d\n", num_traingles/N + 1);

    typedef struct print_list{
        struct print_list *next;
        int idx;
        int local_counter;
    } print_list;

    count = 0;
    int local_counter = 0;
    print_list *start = NULL;

    TriangleNode *curr = start_triangles;

    while(curr != NULL){

        if(count%N == 0){
            strcpy(file_name, result_path);
            strcat(file_name, molecule_name);
            file_number = count/N;
            div = count;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "w");
            // printf("File name: %s", file_name);
            start = NULL;
            local_counter = 0;
        }
        
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
            curr->vert1->local_index = local_counter;
            local_counter++;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert1->local_index, curr->vert1->coordinate1, curr->vert1->coordinate2, curr->vert1->coordinate3);
            fprintf(fptr, "%s", str);
        }

        if (!found2) {
            print_list *new_node = (print_list *)malloc(sizeof(print_list));
            new_node->idx = curr->vert2->index;
            curr->vert2->local_index = local_counter;
            local_counter++;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert2->local_index, curr->vert2->coordinate1, curr->vert2->coordinate2, curr->vert2->coordinate3);
            fprintf(fptr, "%s", str);
        }

        if (!found3) {
            print_list *new_node = (print_list *)malloc(sizeof(print_list));
            new_node->idx = curr->vert3->index;
            curr->vert3->local_index = local_counter;
            local_counter++;
            new_node->next = start;
            start = new_node;

            snprintf(str, sizeof(str), "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", 
                curr->vert3->local_index, curr->vert3->coordinate1, curr->vert3->coordinate2, curr->vert3->coordinate3);
            fprintf(fptr, "%s", str);
        }

        if (curr->vert1->index%N >= 10000 || curr->vert2->index%N >= 10000 || curr->vert3->index%N >= 10000) {
            count++;
            curr = curr->next;
            missed++;
            continue;
        }

        if (curr->vert1->local_index < curr->vert2->local_index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert1->local_index, curr->vert2->local_index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert2->local_index, curr->vert1->local_index);
        }

        if (curr->vert2->local_index < curr->vert3->local_index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert2->local_index, curr->vert3->local_index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert3->local_index, curr->vert2->local_index);
        }

        if (curr->vert3->local_index < curr->vert1->local_index) {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert3->local_index, curr->vert1->local_index);
        } else {
            fprintf(fptr, "CONECT%5d%5d\n", curr->vert1->local_index, curr->vert3->local_index);
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

    printf("MISSED: %d\n", missed);
    printf("Counts in the end: %d\n", count);
}

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