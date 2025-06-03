#include "global.h"
#include <math.h>

int main(){
    FILE *fptr;

    char folder_name[100] = "/home/elia/tesi/tetra/marching_tetrahedron/data/";
    char *path = strcat(folder_name, "generated.bin");

    fptr = fopen(path, "wb");

    if(fptr==NULL){
        printf("null ptr");
    }

    dim_t dx = 50;
    dim_t dy = 50;
    dim_t dz = 50;
    dim_t origin_x = 50;
    dim_t origin_y = 50;
    dim_t origin_z = 50;
    size_t x_d = 5;
    size_t y_d = 2;
    size_t z_d = 2;

    fwrite(&dx, sizeof(dim_t), 1, fptr);
    fwrite(&dy, sizeof(dim_t), 1, fptr);
    fwrite(&dz, sizeof(dim_t), 1, fptr);
    fwrite(&origin_x, sizeof(dim_t), 1, fptr);
    fwrite(&origin_y, sizeof(dim_t), 1, fptr);
    fwrite(&origin_z, sizeof(dim_t), 1, fptr);
    fwrite(&x_d, sizeof(size_t), 1, fptr);
    fwrite(&y_d, sizeof(size_t), 1, fptr);
    fwrite(&z_d, sizeof(size_t), 1, fptr);

    dim_t var = 0;

    for (size_t x=0; x<x_d; x++){
        for (size_t y=0; y<y_d; y++){
            for (size_t z=0; z<z_d; z++){
                if( (x == 1 || x == 2 /*|| x == 3*/) && 
                    (y == 1 /*|| y == 2 || y == 3*/) && 
                    (z == 1 /*|| z == 2 || z == 3*/)){
                    var = 1;
                    printf("PRINTATO\n");
                } else {
                    var = 0;
                }
                // float distance = sqrt((x - x_d / 2) * (x - x_d / 2) + (y - y_d / 2) * (y - y_d / 2) + (z - z_d / 2) * (z - z_d / 2));
                
                // float dx = fabs((double)x - (double)x_d / 2);
                // float dy = fabs((double)y - (double)y_d / 2);
                // float dz = fabs((double)z - (double)z_d / 2);
                // float distance = fmax(fmax(dx, dy), dz);

                // float center_x = x_d / 2.0;
                // float center_y = y_d / 2.0;
                // float center_z = z_d / 2.0;

                // float width = 8.0;
                // float height = 10.0;
                // float depth = 12.0;

                // float dx = fabs(x - center_x);
                // float dy = fabs(y - center_y);
                // float dz = fabs(z - center_z);

                // float distance = 0;
                // if (dx <= width/2 && dy <= height/2 && dz <= depth/2) {
                //     distance = 0;  // Inside the rectangle
                // } else {
                //     distance = 10;  // Outside the rectangle, any value > 9 will work
                // }
                
                // if (distance <= 1) {
                //     printf("AAAAAAAAAAa\n");
                //     var = 1;
                // } else {
                //     var = 0;
                // }
                fwrite(&var, sizeof(dim_t), 1, fptr);
            }
        }
    }

    fclose(fptr);
}