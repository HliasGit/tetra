#include "marching_tetrahedron.h"

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor){

    FILE *fptr;
    fptr = fopen(file_name, "rb"); 

    if(fptr == NULL) {
        fprintf(stderr, "Not able to open the file\n");
        exit(-1);
    }
    
    fread(&(dim->x_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->y_dim), sizeof(size_t), 1, fptr);
    fread(&(dim->z_dim), sizeof(size_t), 1, fptr);

    printf("x_dim: %d\n", dim->x_dim);
    printf("y_dim: %d\n", dim->y_dim);
    printf("z_dim: %d\n", dim->z_dim);

    *tensor = (dim_t*)malloc(sizeof(dim_t)*dim->x_dim*dim->y_dim*dim->z_dim);
    
    for (int k=0; k<dim->z_dim; k++){
        for (int j=0; j<dim->y_dim; j++){
            for (int i=0; i<dim->x_dim; i++){
                fread(&(*tensor)[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim], sizeof(dim_t), 1, fptr);
                // printf("%f\n", (*tensor)[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
            }
        }
    }



    fclose(fptr);
}

void normalize_grid(Dimensions *dim, dim_t **grid, dim_t threshold){

    for(size_t k=0; k<dim->z_dim; k++){
        for (size_t j=0; j<dim->y_dim; j++){
            for (size_t i=0; i<dim->x_dim; i++){
                // printf("val: %f\n", *grid[i + dim->x_dim*j + k + dim->x_dim*dim->y_dim]);
                (*grid)[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim] -= threshold;
            }
        }
    }
}

void marching_tetrahedra(Dimensions *dim, dim_t **grid, int *cube_decomposition, int *count){

    Vertex *coordinates;
    // for every cube in the space
    for(size_t k=0; k<dim->z_dim-1; k++){ // Cube global coordinate
        for (size_t j=0; j<dim->y_dim-1; j++){
            for (size_t i=0; i<dim->x_dim-1; i++){

                // check if every vertex in the cube is F(x,y,x) - C < threshold and in case skip it 
                for (int tetra = 0; tetra<20; tetra+=4){ // for every tetrahedron in a cube
                    coordinates = malloc(4*sizeof(Vertex));
                    for (int idx=tetra; idx<tetra+4; idx ++){ // for every point in a tetrahedra
                        int point = cube_decomposition[idx];
                        // printf("point coord (%d,%d,%d): %d\n",
                        //             i, j, k, point);

                        // find the global tetrahedra coordinates.
                        find_coordinates(idx-tetra, point, i, j, k, &coordinates);

                    }
                    

                    for (int idx=0; idx<4; idx++){
                        printf("coord x: %f\n", coordinates[idx].x);
                        printf("coord y: %f\n", coordinates[idx].y);
                        printf("coord z: %f\n", coordinates[idx].z);
                    }
                    // get the # neg, # zero, # pos for each tetrahedron
                    // get the action value
                    // get the grid points
                    // build the triangle 
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
                printf("%f\n", grid[i + dim->x_dim*j + k*dim->x_dim*dim->y_dim]);
            }
        }
    }
}

void find_coordinates(int idx, const int point, const size_t i, const size_t j, const size_t k, Vertex **coordinates){
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
        one_apex[1] = i;
    } else {
        one_apex[1] = i+1;
    }

    if (k%2 == 0){
        one_apex[2] = i;
    } else {
        one_apex[2] = i+1;
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

int reverseBits3(int n) {
    int reversed = 0;
    for (int i = 0; i < 3; i++) {
        reversed <<= 1;
        reversed |= (n & 1);
        n >>= 1;
    }
    return reversed;
}