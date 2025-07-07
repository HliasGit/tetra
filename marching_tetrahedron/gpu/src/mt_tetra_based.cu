#include "mt_tetra_based.cuh"
#include "mt_tetra_based.h"

#include "marching_tetrahedron_gpu.h"
#include "marching_tetrahedron_gpu.cuh"
#include <cstdio>

__constant__ int4 c_dim;

//TODO MAKE THIS PER CUBE
__global__ void manage_tetra(int n_tetra, float threshold, float *d_grid, float4 *d_results, int *d_CD1, int *d_CD2, int* d_counter){
    // check that you're in the grid
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= n_tetra) return;

    // cube coordinates
    int cube_id = tid / 5;
    int cube_i = cube_id / (c_dim.y * c_dim.z);
    int cube_j = (cube_id / c_dim.z) % c_dim.y;
    int cube_k = cube_id % c_dim.z;
    
    if(cube_k == c_dim.z-1 || cube_j == c_dim.y-1 || cube_i == c_dim.x-1) return;
    
    // cube idx
    int cube_flat = cube_i * (c_dim.y * c_dim.z) + cube_j * c_dim.z + cube_k;
    
    int *lut[2] = {d_CD1, d_CD2};
    int *CD = lut[0];
    
    // offset to get the grid value
    int offset_lut[8] = {  0, 1,
                        c_dim.z, c_dim.z + 1,
                        c_dim.z*c_dim.y + 0, c_dim.z*c_dim.y + 1,
                        c_dim.z*c_dim.y + c_dim.z, c_dim.z*c_dim.y + c_dim.z + 1};

    // id of the thread that aligns to one tetrahedron in the cube
    int tetra_local_id = tid % 5;

    
    //offset to know the value of the decomposition
    int offset_CD = tetra_local_id * 4;
    // Print tid % 5 and offset for tid from 0 to 20
    // if (tid < 21) {
    //     printf("%d %% 5 = %d, offset_CD = %d\n", tid, tid % 5, offset_CD);
    // }

    // int skip_up = 0;
    // int skip_down = 0;

    // int vertex_offset_1 = offset_lut[CD[offset_CD + 0] - 1];
    // int vertex_offset_2 = offset_lut[CD[offset_CD + 1] - 1];
    // int vertex_offset_3 = offset_lut[CD[offset_CD + 2] - 1];
    // int vertex_offset_4 = offset_lut[CD[offset_CD + 3] - 1];

    // float val_1 = d_grid[cube_flat + vertex_offset_1];
    // float val_2 = d_grid[cube_flat + vertex_offset_2];
    // float val_3 = d_grid[cube_flat + vertex_offset_3];
    // float val_4 = d_grid[cube_flat + vertex_offset_4];


    // bool skip = (val_1 < threshold && val_2 < threshold && val_3 < threshold && val_4 < threshold) ||
    //         (val_1 > threshold && val_2 > threshold && val_3 > threshold && val_4 > threshold);
    
    // if (skip) return;
    
    const int lut_one_apex[2] = {0,1};
    float4 one_apex;
    float4 two_apex;

    one_apex.x = cube_i + lut_one_apex[((int)cube_i) & 1];
    one_apex.y = cube_j + lut_one_apex[((int)cube_j) & 1];
    one_apex.z = cube_k + lut_one_apex[((int)cube_k) & 1];
    
    two_apex.x = 2 * cube_i + 1 - one_apex.x;
    two_apex.y = 2 * cube_j + 1 - one_apex.y;
    two_apex.z = 2 * cube_k + 1 - one_apex.z;
    
    float val[4];
    float4 coords[8];
    
    for (int i = 0; i < 4; i++) {
        int cd_idx = CD[offset_CD + i] - 1;
        
        coords[0].x = one_apex.x;  coords[0].y = one_apex.y;  coords[0].z = one_apex.z;  coords[0].w = 0.0f;
        coords[1].x = two_apex.x;  coords[1].y = one_apex.y;  coords[1].z = one_apex.z;  coords[1].w = 0.0f;
        coords[2].x = one_apex.x;  coords[2].y = two_apex.y;  coords[2].z = one_apex.z;  coords[2].w = 0.0f;
        coords[3].x = two_apex.x;  coords[3].y = two_apex.y;  coords[3].z = one_apex.z;  coords[3].w = 0.0f;
        coords[4].x = one_apex.x;  coords[4].y = one_apex.y;  coords[4].z = two_apex.z;  coords[4].w = 0.0f;
        coords[5].x = two_apex.x;  coords[5].y = one_apex.y;  coords[5].z = two_apex.z;  coords[5].w = 0.0f;
        coords[6].x = one_apex.x;  coords[6].y = two_apex.y;  coords[6].z = two_apex.z;  coords[6].w = 0.0f;
        coords[7].x = two_apex.x;  coords[7].y = two_apex.y;  coords[7].z = two_apex.z;  coords[7].w = 0.0f;
        
        val[i] = d_grid[(int)(coords[cd_idx].z + coords[cd_idx].y * c_dim.z +
            coords[cd_idx].x * c_dim.z * c_dim.y)];
    }
    
    bool skip = (val[0] < threshold && val[1] < threshold && val[2] < threshold && val[3] < threshold) ||
    (val[0] > threshold && val[1] > threshold && val[2] > threshold && val[3] > threshold);
    
    if (skip) return;

    int ins_pos = atomicAdd(d_counter, 1);
    
    for (int i = 0; i < 4; i++) {
        int cd_idx = CD[offset_CD + i] - 1;
        
        d_results[ins_pos * 4 + i].x = coords[cd_idx].x;
        d_results[ins_pos * 4 + i].y = coords[cd_idx].y;
        d_results[ins_pos * 4 + i].z = coords[cd_idx].z;
        d_results[ins_pos * 4 + i].w = d_grid[(int)(coords[cd_idx].z + coords[cd_idx].y * c_dim.z +
                 coords[cd_idx].x * c_dim.z * c_dim.y)];
    }
}

__global__ void make_triangles( int n_tetra, int *d_counter, float4 *d_results, Triangle_GPU *d_triangles, coord_t threshold,
                                int *d_pairs, int *d_7val, int *d_0val, int *d_act_val_vec){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= n_tetra) return;

    float4 *point1 = &d_results[tid*4 + 0];
    float4 *point2 = &d_results[tid*4 + 1];
    float4 *point3 = &d_results[tid*4 + 2];
    float4 *point4 = &d_results[tid*4 + 3];
    
    // Sort the four points by their .w value (ascending order)
    float4 *points[4] = {point1, point2, point3, point4};
    // Simple bubble sort for 4 elements
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3 - i; ++j) {
            if (points[j]->w > points[j + 1]->w) {
                float4* tmp = points[j];
                points[j] = points[j + 1];
                points[j + 1] = tmp;
            }
        }
    }

    int count_lt = 0, count_eq = 0, count_gt = 0;
    for (int i = 0; i < 4; ++i) {
        if (points[i]->w <  threshold) ++count_lt;
        else if (points[i]->w == threshold) ++count_eq;
        else ++count_gt;
    }
    
    int act_val = d_act_val_vec[count_lt + count_eq*5];
        
    if (act_val!= 0){

        int *use_pairs = &d_pairs[6*(act_val-1)];
        int ins_pos = atomicAdd(d_counter, 1);
        make_triangle_from_tetra(points[0], points[1], points[2], points[3], &d_triangles[ins_pos], use_pairs);
        // printf("d_counter after triangle: %d\n", ins_pos);
        
        if(act_val == 7){
            
            use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 
            int ins_pos = atomicAdd(d_counter, 1);
            make_triangle_from_tetra(points[0], points[1], points[2], points[3], &d_triangles[ins_pos], use_pairs);
            
            // atomicAdd(d_7val, 1);
        }
    }
    // } else {
    //     // printf("SHOULDNT BE HERE\n");
    //     // atomicAdd(d_0val, 1);
    // }
}

__device__ void make_triangle_from_tetra(float4 *point1, float4 *point2, float4 *point3, float4 *point4,
                                    Triangle_GPU *triangle, int *pairs){

    float4 arr[4] = {*point1, *point2, *point3, *point4};
    
    int idx1 = pairs[0] - 1;
    int idx2 = pairs[1] - 1;

    triangle->v1.x = ((coord_t)arr[idx1].x + (coord_t)arr[idx2].x) / 2.0;
    triangle->v1.y = ((coord_t)arr[idx1].y + (coord_t)arr[idx2].y) / 2.0;
    triangle->v1.z = ((coord_t)arr[idx1].z + (coord_t)arr[idx2].z) / 2.0;

    idx1 = pairs[2] - 1;
    idx2 = pairs[3] - 1;


    triangle->v2.x = ((coord_t)arr[idx1].x + (coord_t)arr[idx2].x) / 2.0;
    triangle->v2.y = ((coord_t)arr[idx1].y + (coord_t)arr[idx2].y) / 2.0;
    triangle->v2.z = ((coord_t)arr[idx1].z + (coord_t)arr[idx2].z) / 2.0;

    idx1 = pairs[4] - 1;
    idx2 = pairs[5] - 1;

    triangle->v3.x = ((coord_t)arr[idx1].x + (coord_t)arr[idx2].x) / 2.0;
    triangle->v3.y = ((coord_t)arr[idx1].y + (coord_t)arr[idx2].y) / 2.0;
    triangle->v3.z = ((coord_t)arr[idx1].z + (coord_t)arr[idx2].z) / 2.0;
}

void remove_unnecessary_tetrahedra( dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                    int *number_relevant_cubes,
                                    cube_gpu_SoA **d_relevant_cubes,
                                    double *time, int *CD1, int *CD2, int *pairs, int *act_val_vec){
    
    // //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //


    int n_tetra = cubes_in_domain * 5;
    int n_threads = 1024;
    int n_blocks = (n_tetra + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to manage tetra\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

    int *d_CD1;
    print_cuda_error(cudaMalloc(&d_CD1, sizeof(int)*20), "CD1 malloc");
    print_cuda_error(cudaMemcpy(d_CD1, CD1, sizeof(int)*20, cudaMemcpyHostToDevice), "CD1 memcpy");

    int *d_CD2;
    print_cuda_error(cudaMalloc(&d_CD2, sizeof(int)*20), "CD2 malloc managed");
    print_cuda_error(cudaMemcpy(d_CD1, CD1, sizeof(int)*20, cudaMemcpyHostToDevice), "CD2 memcpy");

    int *d_counter;
    print_cuda_error(cudaMallocManaged(&d_counter, sizeof(int)), "d_counter malloc managed");
    print_cuda_error(cudaMemset(d_counter, 0, sizeof(int)), "d_counter memset managed");

    float4 *d_results;
    print_cuda_error(cudaMallocManaged(&d_results, sizeof(float4) * n_tetra * 4), "results malloc managed");

    int *d_pairs;
    print_cuda_error(cudaMalloc(&d_pairs, sizeof(int) * 48), "pairs malloc"); // 8*6 = 48
    print_cuda_error(cudaMemcpy(d_pairs, pairs, sizeof(int) * 48, cudaMemcpyHostToDevice), "pairs memcpy");

    int *d_7val;
    print_cuda_error(cudaMallocManaged(&d_7val, sizeof(int)), "d_7val malloc managed");
    print_cuda_error(cudaMemset(d_7val, 0, sizeof(int)), "d_7val memset managed");

    int *d_0val;
    print_cuda_error(cudaMallocManaged(&d_0val, sizeof(int)), "d_0val malloc managed");
    print_cuda_error(cudaMemset(d_0val, 0, sizeof(int)), "d_0val memset managed");

    int *d_act_val_vec;
    print_cuda_error(cudaMalloc(&d_act_val_vec, sizeof(int)*25), "cudaMalloc failed for d_act_val_vec alloc: %s\n");
    print_cuda_error(cudaMemcpy(d_act_val_vec, act_val_vec, sizeof(int) * 25, cudaMemcpyHostToDevice), "cudaMemcpy failed for d_act_val_vec copy: %s\n");



    /////////////// KERNEL TO DIVIDE IN TETRA THE SPACE ///////////////
    
    // Take start time
    cudaEvent_t start, stop;
    float elapsedTime = 0.0f;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    manage_tetra<<<n_blocks, n_threads>>>(  n_tetra, threshold, d_grid, d_results, d_CD1, d_CD2, d_counter);
    print_cuda_error(cudaDeviceSynchronize(), "cuda device synch");

    // Take end time
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error after manage_tetra: %s\n", cudaGetErrorString(err));
    }

    // WRITE TO FILE THE DOMAIN DECOMP
    // int local_counter  = 0;
    // FILE *fp = fopen("tetra.pdb", "w");
    // if (!fp) {
    //     printf("Failed to open tetra.pdb for writing\n");
    // } else {
    //     for (int i = 4500000; i < + 4500000 + 100000 -4; i = i+4) {

    //         int id1 = i;
    //         int id2 = i+1;
    //         int id3 = i+2;
    //         int id4 = i+3;

    //         int id1_ = local_counter;
    //         int id2_ = local_counter+1;
    //         int id3_ = local_counter+2;
    //         int id4_ = local_counter+3;

    //         char atom_type1 = (d_results[id1].w < threshold) ? 'C' : 'O';
    //         char atom_type2 = (d_results[id2].w < threshold) ? 'C' : 'O';
    //         char atom_type3 = (d_results[id3].w < threshold) ? 'C' : 'O';
    //         char atom_type4 = (d_results[id4].w < threshold) ? 'C' : 'O';

    //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
    //                 id1_, atom_type1, d_results[id1].x, d_results[id1].y, d_results[id1].z, atom_type1);
    //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
    //                 id2_, atom_type2, d_results[id2].x, d_results[id2].y, d_results[id2].z, atom_type2);
    //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
    //                 id3_, atom_type3, d_results[id3].x, d_results[id3].y, d_results[id3].z, atom_type3);
    //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
    //                 id4_, atom_type4, d_results[id4].x, d_results[id4].y, d_results[id4].z, atom_type4);

    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id1_, id2_);
    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id1_, id3_);
    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id1_, id4_);

    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id2_, id3_);
    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id2_, id4_);

    //         fprintf(fp, "CONECT%5d%5d\n",
    //                 id3_, id4_);

    //         local_counter+=4;
    //     }
    //     fclose(fp);
    //     printf("Wrote %d atoms to tetra.pdb\n", (*d_counter));
    // }
    
    /////////////// KERNEL TO MAKE THE TRIANGLES ///////////////

    n_tetra = *d_counter;
    n_threads = 512;
    n_blocks = (n_tetra + n_threads - 1) / n_threads;

    Triangle_GPU *d_triangles;
    print_cuda_error(cudaMallocManaged(&d_triangles, sizeof(Triangle_GPU) * (*d_counter)*2), "d_triangles malloc managed");

    *d_counter = 0;

    // Take start time for triangle kernel
    cudaEvent_t tri_start, tri_stop;
    float tri_elapsedTime = 0.0f;
    cudaEventCreate(&tri_start);
    cudaEventCreate(&tri_stop);
    cudaEventRecord(tri_start, 0);

    make_triangles<<<n_blocks, n_threads>>>(n_tetra, d_counter, d_results, d_triangles, threshold,
                                            d_pairs, d_7val, d_0val, d_act_val_vec);
    print_cuda_error(cudaDeviceSynchronize(), "cuda device synch");

    // Take end time for triangle kernel
    cudaEventRecord(tri_stop, 0);
    cudaEventSynchronize(tri_stop);
    cudaEventElapsedTime(&tri_elapsedTime, tri_start, tri_stop);

    cudaEventDestroy(tri_start);
    cudaEventDestroy(tri_stop);

    printf("Triangle kernel time: %.3f ms\n", tri_elapsedTime);

    printf("Total # of triangles:                           %d\n", n_tetra);

    cudaError_t last_err = cudaGetLastError();
    if (last_err != cudaSuccess) {
        printf("CUDA error after make_triangles: %s\n", cudaGetErrorString(last_err));
    }

    // WRITE TO FILE THE TRIANGLES
    int final_triangles_number = *d_counter;
    // int local_counter  = 0;
    
    // FILE *fp_tri = fopen("triangle.pdb", "w");
    // if (!fp_tri) {
    //     printf("Failed to open triangle.pdb for writing\n");
    // } else {
    //     for (int i = 0; i < final_triangles_number; ++i) {

    //         int id1 = local_counter * 3 + 1;
    //         int id2 = local_counter * 3 + 2;
    //         int id3 = local_counter * 3 + 3;

    //         fprintf(fp_tri, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n",
    //                 id1, d_triangles[i].v1.x, d_triangles[i].v1.y, d_triangles[i].v1.z);
    //         fprintf(fp_tri, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n",
    //                 id2, d_triangles[i].v2.x, d_triangles[i].v2.y, d_triangles[i].v2.z);
    //         fprintf(fp_tri, "ATOM  %5d C    PSE A   1    %8.2f%8.2f%8.2f 1.00  1.00           C\n",
    //                 id3, d_triangles[i].v3.x, d_triangles[i].v3.y, d_triangles[i].v3.z);

    //         fprintf(fp_tri, "CONECT%5d%5d\n", id1, id2);
    //         fprintf(fp_tri, "CONECT%5d%5d\n", id2, id3);
    //         fprintf(fp_tri, "CONECT%5d%5d\n", id1, id3);

    //         local_counter++;
    //     }
    //     fclose(fp_tri);
    // }
    printf("Wrote %d triangles to triangle.pdb\n", final_triangles_number);

    *time = tri_elapsedTime + elapsedTime;
}

void load_to_const_tetra(Dimensions *dimensions){

    printf("Dimensions: x=%d, y=%d, z=%d\n", dimensions->x_dim, dimensions->y_dim, dimensions->z_dim);
    int4 tmp;
    tmp.x = dimensions->x_dim;
    tmp.y = dimensions->y_dim;
    tmp.z = dimensions->z_dim;

    tmp.w = 0.0f;
    cudaMemcpyToSymbol(c_dim, &tmp, sizeof(int4));

    printf("Dim copied to const mem\n");
}