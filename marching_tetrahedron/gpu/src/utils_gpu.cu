#include "marching_tetrahedron_gpu.h"
#include "marching_tetrahedron_gpu.cuh"

/// @brief Function to allocate on device the ED scalar field 
/// @param d_grid Pointer to the device memory
/// @param grid Pointer to the host memory
/// @param size Number of elements to be allocated
void allocate_d_grid(dim_t **d_grid, dim_t *grid, size_t size){
    cudaError_t err;
    err = cudaMalloc(d_grid, sizeof(dim_t) * size);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(err));
    }
    err = cudaMemcpy(*d_grid, grid, sizeof(dim_t) * size, cudaMemcpyHostToDevice);
    if (err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(err));
    }
}

/// @brief Debug utility to print the pints after the pre-processing pass
/// @param d_relevant_cubes Pointer to the device allocated DS of cube_gpu 
/// @param number_relevant_cubes 
void print_relevant_points(cube_gpu *d_relevant_cubes, int *number_relevant_cubes){
    // Remove "points" directory if it exists, then create it and change into it
    struct stat st = {0};
    if (stat("points", &st) == 0) {
        // Directory exists, remove all .pdb files inside
        system("rm -rf points");
    }
    mkdir("points", 0700);
    chdir("points");

    printf("PRINTING\n");

    int file_number = 0;
    int local_counter = 0;
    FILE *pdb_file;
    for (int i = 0; i < *number_relevant_cubes; ++i) {
        
        if(i%100000 == 0){
            char filename[64];
            snprintf(filename, sizeof(filename), "points_%d.pdb", file_number);
            pdb_file = fopen(filename, "w");

            if (!pdb_file) {
                fprintf(stderr, "Failed to open points.pdb for writing.\n");
            } 

            file_number++;
            local_counter=0;
        }
            
        int x = d_relevant_cubes[i].x;
        int y = d_relevant_cubes[i].y;
        int z = d_relevant_cubes[i].z;

        // printf("REMARK x=%d y=%d z=%d\n", x, y, z);
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, (float)x, (float)y, (float)z);
        local_counter++;

        if((i+1)%100000 == 0){
            fclose(pdb_file);
        }
    }
    printf("Relevant points written to points.pdb\n");
    chdir("..");
}

/// @brief Debug utility to print the pints after the pre-processing pass
/// @param d_relevant_cubes Pointer to the device allocated DS of cube_gpu 
/// @param number_relevant_cubes 
void print_relevant_points_soa(cube_gpu_SoA *d_relevant_cubes, int *number_relevant_cubes){
    // Remove "points" directory if it exists, then create it and change into it

    float4* tmp = d_relevant_cubes->coord_idx;

    struct stat st = {0};
    if (stat("points", &st) == 0) {
        // Directory exists, remove all .pdb files inside
        system("rm -rf points");
    }
    mkdir("points", 0700);
    chdir("points");

    printf("PRINTING\n");

    int file_number = 0;
    int local_counter = 0;
    FILE *pdb_file;
    for (int i = 0; i < *number_relevant_cubes; ++i) {
        
        if(i%100000 == 0){
            char filename[64];
            snprintf(filename, sizeof(filename), "points_%d.pdb", file_number);
            pdb_file = fopen(filename, "w");

            if (!pdb_file) {
                fprintf(stderr, "Failed to open points.pdb for writing.\n");
            } 

            file_number++;
            local_counter=0;
        }

        float x = tmp[i].x;
        float y = tmp[i].y;
        float z = tmp[i].z;

        // printf("REMARK x=%d y=%d z=%d\n", x, y, z);
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, (float)x, (float)y, (float)z);
        local_counter++;

        if((i+1)%100000 == 0){
            fclose(pdb_file);
        }
    }
    printf("Relevant points written to points.pdb\n");
    chdir("..");
}

/// @brief Utility function to print the GPU triangles
/// @param triangles Pointer to the triangles
/// @param triangles_number 
/// @param molecule_name 
/// @param molecule_path 
void print_triangles(   Triangle_GPU *triangles, int *triangles_number,
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

    printf("Number of triangles to print: %d\n", *triangles_number);

    FILE *pdb_file;
    for (int i = 0; i < *triangles_number; i++) {
        
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
            local_counter++;
        }

            
        if (triangles[i].v1.x == 0.0 && triangles[i].v1.y == 0.0 && triangles[i].v1.z == 0.0 &&
            triangles[i].v2.x == 0.0 && triangles[i].v2.y == 0.0 && triangles[i].v2.z == 0.0 &&
            triangles[i].v3.x == 0.0 && triangles[i].v3.y == 0.0 && triangles[i].v3.z == 0.0) {
            empty++;
            continue;
        }

        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v1.x, triangles[i].v1.y, triangles[i].v1.z);
        local_counter++;
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v2.x, triangles[i].v2.y, triangles[i].v2.z);
        local_counter++;
        fprintf(pdb_file, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
        // printf("ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n", local_counter, triangles[i].v3.x, triangles[i].v3.y, triangles[i].v3.z);
        local_counter++;
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 3, local_counter - 2);
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 2, local_counter - 1);
        fprintf(pdb_file, "CONECT%5d%5d\n", local_counter - 3, local_counter - 1);

        if(local_counter == 99994){
            local_counter = 0;
            // printf("empty %d\n", empty);
            // printf("file_number: %d\n", file_number);
            fclose(pdb_file);
        }
    }
    printf("Relevant points written\n");


    
}

/// @brief Utility function to print CUDA errors with a message
/// @param err CUDA error code
/// @param msg Custom message to print with the error
void print_cuda_error(cudaError_t err, const char* msg){
    if (err != cudaSuccess) {
        fprintf(stderr, "error failed: %s\nMsg: %s\n", cudaGetErrorString(err), msg);
    }
}