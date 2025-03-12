#include "marching_tetrahedron.h"

int main(){
    FILE *fptr;

    char folder_name[100] = "/home/elia/tesi/code/marching_tetrahedron/test/";
    char *path = strcat(folder_name, "file.bin");

    fptr = fopen(path, "wb");

    if(fptr==NULL){
        printf("null ptr");
    }

    size_t x_d = 2;
    size_t y_d = 2;
    size_t z_d = 2;

    fwrite(&x_d, sizeof(size_t), 1, fptr);
    fwrite(&y_d, sizeof(size_t), 1, fptr);
    fwrite(&z_d, sizeof(size_t), 1, fptr);

    float var = 0;

    for (int z=0; z<z_d; z++){
        for (int y=0; y<y_d; y++){
            for (int x=0; x<x_d; x++){
                // if(x == x_d/2 && y == y_d/2 && z == z_d/2){
                if(x == 1 && y == 1 && z == 1){
                    var = 1;
                    printf("PRINTATO\n");
                } else {
                    var = 0;
                }
                fwrite(&var, sizeof(float), 1, fptr);
            }
        }
    }
    fclose(fptr);
}