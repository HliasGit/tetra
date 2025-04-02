#include "utils.h"
#include "struct.h"
#include <sys/stat.h>
#include <string.h>

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

    int *idxs = print_atoms_separated(p->triangles, molecule_name, molecule_path, num_triangles);
    print_connections_separated(p->triangles, molecule_name, molecule_path, idxs);
}

/**
 * @brief Print the atoms on the splitted file
 * 
 * @param curr Pointer to the triangle list
 * @param molecule_name Ptr to the name of the molecule
 * @param result_path Ptr to the result path of the molecule
 * @param num_triangles Number of total triangles generated
 *
 * @return Pointer to the min offsets
 */
int *print_atoms_separated(TriangleNode *curr, char *molecule_name, char *result_path, int num_traingles){

    FILE *fptr;
    char file_name[200];
    int N = 2500;
    int min = 0;
    int count = 0; 
    int div = 0;
    int file_number = 0;
    int *offset = (int *)malloc((num_traingles/N + 1) * sizeof(int));

    if (offset == NULL) {
        fprintf(stderr, "Memory allocation failed for offset array\n");
        exit(-1);
    }

    printf(result_path);
    printf("\n");

    printf("Number of files afbafub: %d\n", num_traingles/N + 1);

    typedef struct print_list{
        struct print_list *next;
        int idx;
    } print_list;

    TriangleNode *counter = curr;

    while(counter != NULL){

        if((count+1)%N == 0 || counter->next == NULL){
            offset[(int)(count/N)] = min;
            // printf("count: %d\n", count);
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

    // for (int i = 0; i < (num_traingles/N + 1); i++) {
    //     printf("Offset[%d]: %d\n", i, offset[i]);
    // }

    count = 0;
    print_list *start = NULL;

    while(curr != NULL){

        if(count%N == 0){

            // printf("QUI");
            
            strcpy(file_name, result_path);
            strcat(file_name, molecule_name);
            file_number = count/N;
            div = count;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "w");
            // printf("File name: %s", file_name);
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

/**
 * @brief Print the connections on the splitted file
 * 
 * @param curr Pointer to the triangle list
 * @param molecule_name Ptr to the name of the molecule
 * @param result_path Ptr to the result path of the molecule
 * @param offset Pointer to array containing the min offsets
 */
void print_connections_separated(TriangleNode *curr, char *molecule_name, char *result_path, int *offsets){
    int N = 2500;
    char file_name[200];

    int count = 0; 

    int missed = 0;

    FILE *fptr;

    int div = 0;

    int file_number = 0;

    while(curr != NULL){

        if(count%N == 0){
            strcpy(file_name, result_path);
            strcat(file_name, molecule_name);
            file_number = count/N;
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            fptr = fopen(file_name, "a");
        }

        div = offsets[count/N];

        if (curr->vert1->index-div >= 10000 || curr->vert2->index-div >= 10000 || curr->vert3->index-div >= 10000) {
            count++;
            curr = curr->next;
            missed++;
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

    free(offsets);
    printf("MISSED: %d\n", missed);
    printf("Counts in the end: %d\n", count);
}