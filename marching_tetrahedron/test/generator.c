#include "global.h"

int main(){
    FILE *fptr;

    char folder_name[100] = "/home/elia/tesi/code/marching_tetrahedron/test/";
    char *path = strcat(folder_name, "file.bin");

    fptr = fopen(path, "wb");

    if(fptr==NULL){
        printf("null ptr");
    }

    double dx = 50;
    double dy = 50;
    double dz = 50;
    double origin_x = 50;
    double origin_y = 50;
    double origin_z = 50;
    size_t x_d = 150;
    size_t y_d = 150;
    size_t z_d = 150;

    fwrite(&dx, sizeof(double), 1, fptr);
    fwrite(&dy, sizeof(double), 1, fptr);
    fwrite(&dz, sizeof(double), 1, fptr);
    fwrite(&origin_x, sizeof(double), 1, fptr);
    fwrite(&origin_y, sizeof(double), 1, fptr);
    fwrite(&origin_z, sizeof(double), 1, fptr);
    fwrite(&x_d, sizeof(size_t), 1, fptr);
    fwrite(&y_d, sizeof(size_t), 1, fptr);
    fwrite(&z_d, sizeof(size_t), 1, fptr);

    float var = 0;

    for (int x=0; x<x_d; x++){
        for (int y=0; y<y_d; y++){
            for (int z=0; z<z_d; z++){
                // if(x == 1 && y == 1 && z == 1){
                //     var = 1;
                //     printf("PRINTATO\n");
                // } else {
                //     var = 0;
                // }
                float distance = sqrt((x - x_d / 2) * (x - x_d / 2) + (y - y_d / 2) * (y - y_d / 2) + (z - z_d / 2) * (z - z_d / 2));
                // float dx = fabs(x - x_d / 2);
                // float dy = fabs(y - y_d / 2);
                // float dz = fabs(z - z_d / 2);
                // float distance = fmax(fmax(dx, dy), dz);
                
                if (distance <= 7) {
                    var = 1;
                } else {
                    var = 0;
                }
                fwrite(&var, sizeof(float), 1, fptr);
            }
        }
    }

    fclose(fptr);
}