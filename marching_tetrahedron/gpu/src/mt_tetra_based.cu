#include "mt_tetra_based.cuh"
#include "mt_tetra_based.h"

#include "marching_tetrahedron_gpu.h"
#include "marching_tetrahedron_gpu.cuh"
#include <cstdio>

__constant__ int3 c_dim;


__global__ void remove_unnecessary_cubes_mixed(dim_t* grid, int *counter,
                                                size_t size, double threshold,
                                                cube_gpu_SoA* d_relevant_cubes) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;

    int i = idx / (c_dim.y * c_dim.z);
    int j = (idx / c_dim.z) % c_dim.y;
    int k = idx % c_dim.z;
    
    extern __shared__ float shared_mem[];
    float* s_mem_1 = shared_mem;
    float* s_mem_2 = s_mem_1 + blockDim.x + 1;
    float* s_mem_3 = s_mem_2 + blockDim.x + 1;
    float* s_mem_4 = s_mem_3 + blockDim.x + 1;

    bool all_in = 0;
    bool all_out = 0;

    if(k != c_dim.z || j != c_dim.y || i != c_dim.x){
        s_mem_1[threadIdx.x] = grid[idx]; //caricato la prima fila
        if(threadIdx.x == blockDim.x-1){ // ultimo
            s_mem_1[blockDim.x] = grid[idx + 1]; //caricato ultimo della prima fila
        }

        s_mem_2[threadIdx.x] = grid[idx + c_dim.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_2[blockDim.x] = grid[idx + c_dim.z + 1];
        }

        s_mem_3[threadIdx.x] = grid[idx + c_dim.z * c_dim.y];
        if(threadIdx.x == blockDim.x-1){
            s_mem_3[blockDim.x] = grid[idx + c_dim.z * c_dim.y + 1];
        }

        s_mem_4[threadIdx.x] = grid[idx + c_dim.z * c_dim.y + c_dim.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_4[blockDim.x] = grid[idx + c_dim.z * c_dim.y + c_dim.z + 1];
        }

        __syncthreads();
            
        if( s_mem_1[threadIdx.x] > threshold && s_mem_1[threadIdx.x+1] > threshold &&
            s_mem_2[threadIdx.x] > threshold && s_mem_2[threadIdx.x+1] > threshold &&
            s_mem_3[threadIdx.x] > threshold && s_mem_3[threadIdx.x+1] > threshold &&
            s_mem_4[threadIdx.x] > threshold && s_mem_4[threadIdx.x+1] > threshold)
            {
            all_out = 1;
        }
        if( s_mem_1[threadIdx.x] < threshold && s_mem_1[threadIdx.x+1] < threshold &&
            s_mem_2[threadIdx.x] < threshold && s_mem_2[threadIdx.x+1] < threshold &&
            s_mem_3[threadIdx.x] < threshold && s_mem_3[threadIdx.x+1] < threshold &&
            s_mem_4[threadIdx.x] < threshold && s_mem_4[threadIdx.x+1] < threshold)
            {
            all_in = 1;
        }
    }

    

    if (all_out == 0 && all_in == 0){
        int insert_pos = atomicAdd(counter, 1);
        d_relevant_cubes->coord_idx[insert_pos].w = idx;
        d_relevant_cubes->coord_idx[insert_pos].x = i;
        d_relevant_cubes->coord_idx[insert_pos].y = j;
        d_relevant_cubes->coord_idx[insert_pos].z = k;
    }
}

//TODO MAKE THIS PER CUBE
__global__ void manage_tetra(int n_total_tetra, dim_t threshold, dim_t *d_grid, float4 *d_active_tetrahedra, int *d_CD1, int* d_num_active_tetrahedra){
    // check that you're in the grid
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= n_total_tetra) return;

    // cube coordinates
    int cube_id = tid / 5;
    int cube_i = cube_id / (c_dim.y * c_dim.z);
    int cube_j = (cube_id / c_dim.z) % c_dim.y;
    int cube_k = cube_id % c_dim.z;
    
    if(cube_k == c_dim.z-1 || cube_j == c_dim.y-1 || cube_i == c_dim.x-1) return;
    
    // cube idx
    int cube_flat = cube_i * (c_dim.y * c_dim.z) + cube_j * c_dim.z + cube_k;
    
    // offset to get the grid value
    int offset_lut[8] = {  0, 1,
                        c_dim.z, c_dim.z + 1,
                        c_dim.z*c_dim.y + 0, c_dim.z*c_dim.y + 1,
                        c_dim.z*c_dim.y + c_dim.z, c_dim.z*c_dim.y + c_dim.z + 1};

    // id of the thread that aligns to one tetrahedron in the cube
    int tetra_local_id = tid % 5;
    
    //offset to know the value of the decomposition
    int offset_CD = tetra_local_id * 4;
    
    const int lut_one_apex[2] = {0,1};
    float4 one_apex;
    float4 two_apex;

    one_apex.x = cube_i + lut_one_apex[((int)cube_i) & 1];
    one_apex.y = cube_j + lut_one_apex[((int)cube_j) & 1];
    one_apex.z = cube_k + lut_one_apex[((int)cube_k) & 1];
    
    two_apex.x = 2 * cube_i + 1 - one_apex.x;
    two_apex.y = 2 * cube_j + 1 - one_apex.y;
    two_apex.z = 2 * cube_k + 1 - one_apex.z;
    
    dim_t val[4];
    float4 coords[8];
    
    for (int i = 0; i < 4; i++) {
        int cd_idx = d_CD1[offset_CD + i] - 1;
        
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

    int ins_pos = atomicAdd(d_num_active_tetrahedra, 1);
    
    for (int i = 0; i < 4; i++) {
        int cd_idx = d_CD1[offset_CD + i] - 1;
        
        d_active_tetrahedra[ins_pos * 4 + i].x = coords[cd_idx].x;
        d_active_tetrahedra[ins_pos * 4 + i].y = coords[cd_idx].y;
        d_active_tetrahedra[ins_pos * 4 + i].z = coords[cd_idx].z;
        d_active_tetrahedra[ins_pos * 4 + i].w = d_grid[(int)(coords[cd_idx].z + coords[cd_idx].y * c_dim.z +
                 coords[cd_idx].x * c_dim.z * c_dim.y)];
    }
}

__global__ void compute_apex_tetra_based(   dim_t *grid,
                                            cube_gpu_SoA *d_relevant_cubes,
                                            int number_relevant_cubes,
                                            int *d_cube_decomposition,
                                            float4 *d_active_tetrahedra){
                        
    float3 one_apex;
    float3 two_apex;
    float3 coord_idx;

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int local_tid = tid % 5;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes->coord_idx[tid].x >= c_dim.x-1 ||
        d_relevant_cubes->coord_idx[tid].y >= c_dim.y-1 ||
        d_relevant_cubes->coord_idx[tid].z >= c_dim.z-1 ){
            return;
        }

    coord_idx.x = d_relevant_cubes->coord_idx[tid].x;
    coord_idx.y = d_relevant_cubes->coord_idx[tid].y;
    coord_idx.z = d_relevant_cubes->coord_idx[tid].z;

    const int lut_one_apex[2] = {0,1};

    one_apex.x = coord_idx.x + lut_one_apex[((int)coord_idx.x) & 1];
    one_apex.y = coord_idx.y + lut_one_apex[((int)coord_idx.y) & 1];
    one_apex.z = coord_idx.z + lut_one_apex[((int)coord_idx.z) & 1];
    
    two_apex.x = 2 * coord_idx.x + 1 - one_apex.x;
    two_apex.y = 2 * coord_idx.y + 1 - one_apex.y;
    two_apex.z = 2 * coord_idx.z + 1 - one_apex.z;

    float4 coords[8];
    coords[0] = make_float4(one_apex.x, one_apex.y, one_apex.z, 0.0f);
    coords[1] = make_float4(two_apex.x, one_apex.y, one_apex.z, 0.0f);
    coords[2] = make_float4(one_apex.x, two_apex.y, one_apex.z, 0.0f);
    coords[3] = make_float4(two_apex.x, two_apex.y, one_apex.z, 0.0f);
    coords[4] = make_float4(one_apex.x, one_apex.y, two_apex.z, 0.0f);
    coords[5] = make_float4(two_apex.x, one_apex.y, two_apex.z, 0.0f);
    coords[6] = make_float4(one_apex.x, two_apex.y, two_apex.z, 0.0f);
    coords[7] = make_float4(two_apex.x, two_apex.y, two_apex.z, 0.0f);
    
    for(int v=0; v<8; v++){
        coords[v].w = grid[(int)(   coords[v].z +
                                    coords[v].y * c_dim.z +
                                    coords[v].x * c_dim.z * c_dim.y)];
    }

    for(int tetra=0; tetra<5; tetra++){
        for(int point=0; point<4; point++){
            d_active_tetrahedra[tid] = coords[d_cube_decomposition[tetra*local_tid+point]];
        }
    }
}

__global__ void from_cubes_to_tetra( int number_relevant_cubes,
                                cube_vertices_points_SoA *d_cube_points_coordinates,
                                float4 *d_active_tetrahedra,
                                int *d_cube_decomposition)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= number_relevant_cubes) return;

    for(int tetra=0; tetra<5; tetra++){
        for(int point=0; point<4; point++){
            d_active_tetrahedra[20*tid + tetra*4 + point] = d_cube_points_coordinates[tid*8 + d_cube_decomposition[tetra*4 + point]-1].val;
            // d_active_tetrahedra[20*tid + tetra*4 + point] = make_float4(2.0,2.0,2.0,2.0);
        }
    }
}

__global__ void make_triangles( int n_active_tetra,
                                int *d_num_generated_triangles,
                                float4 *d_active_tetrahedra,
                                Triangle_GPU *d_triangles,
                                dim_t threshold,
                                int *d_pairs,
                                int *d_7val,
                                int *d_0val,
                                int *d_act_val_vec){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= n_active_tetra) return;

    float4 *point1 = &d_active_tetrahedra[tid*4 + 0];
    float4 *point2 = &d_active_tetrahedra[tid*4 + 1];
    float4 *point3 = &d_active_tetrahedra[tid*4 + 2];
    float4 *point4 = &d_active_tetrahedra[tid*4 + 3];
    
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
        int ins_pos = atomicAdd(d_num_generated_triangles, 1);
        make_triangle_from_tetra(points[0], points[1], points[2], points[3], &d_triangles[ins_pos], use_pairs);
        // printf("d_num_generated_triangles after triangle: %d\n", ins_pos);
        
        if(act_val == 7){
            
            use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 
            int ins_pos = atomicAdd(d_num_generated_triangles, 1);
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

void remove_unnecessary_tetrahedra( dim_t *d_grid,
                                    size_t cubes_in_domain,
                                    int *CD1,
                                    dim_t threshold,
                                    float4 **d_active_tetrahedra,
                                    int **d_num_active_tetrahedra,
                                    float *time){
    
    // //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

    int n_total_tetra = cubes_in_domain * 5;
    int n_threads = 1024;
    int n_blocks = (n_total_tetra + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to manage tetra\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

    int *d_CD1;
    print_cuda_error(cudaMalloc(&d_CD1, sizeof(int)*20), "CD1 malloc");
    print_cuda_error(cudaMemcpy(d_CD1, CD1, sizeof(int)*20, cudaMemcpyHostToDevice), "CD1 memcpy");

    print_cuda_error(cudaMallocManaged(d_num_active_tetrahedra, sizeof(int)), "d_num_active_tetrahedra malloc managed");
    print_cuda_error(cudaMemset(*d_num_active_tetrahedra, 0, sizeof(int)), "d_num_active_tetrahedra memset managed");

    print_cuda_error(cudaMallocManaged(d_active_tetrahedra, sizeof(float4) * n_total_tetra * 4), "results malloc managed");

    /////////////// RESOURCES NEEDED FOT THE FOLLOWING KERNEL ///////////////

    /*  int n_total_tetra                     number of total tetrahedra
        float threshold         
        float *d_grid                   grid on device
        float4 *d_active_tetrahedra     active tetrahedra on device
        int *d_CD1                      cube decomp 1 on device
        int *d_num_active_tetrahedra    number of active tetrahedrons on device
        double *time
    */

    /////////////// KERNEL TO DIVIDE IN TETRA THE SPACE ///////////////
    
    // Take start time
    cudaEvent_t start, stop;
    float elapsedTime = 0.0f;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);

    manage_tetra<<<n_blocks, n_threads>>>(  n_total_tetra, threshold, d_grid, *d_active_tetrahedra, d_CD1, *d_num_active_tetrahedra);
    print_cuda_error(cudaDeviceSynchronize(), "cuda device synch");

    // Take end time
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    printf("Preprocessing (tetrahedra selection) kernel time: %.3f ms\n", elapsedTime);

    *time += elapsedTime;

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error after manage_tetra: %s\n", cudaGetErrorString(err));
    }

}

void tetra_mt(  int *d_num_active_tetrahedra,
                float4 *d_active_tetrahedra,
                dim_t threshold,
                int* pairs,
                int *act_val_vec,
                float *time){

    // Here I write the data needed for the computation

    /*  int *d_num_active_tetrahedra                     number of tetra
        float4 *d_active_tetrahedra     store the tetrahedron coordinates and values
        dim_t threshold
        int *pairs                      pairs on GPU
        int *act_val_vec,
        float *time                     overall kernel time
    */

    /*  What needs to be common in the 2 kernels?
        float4 *d_active_tetrahedra
        float *time
    */

    /////////////// KERNEL TO MAKE THE TRIANGLES ///////////////

    int n_active_tetra = *d_num_active_tetrahedra;
    int n_threads = 512;
    int n_blocks = (n_active_tetra + n_threads - 1) / n_threads;

    // Allocate the Triangles GPU
    Triangle_GPU *d_triangles;
    print_cuda_error(cudaMallocManaged(&d_triangles, sizeof(Triangle_GPU) * (n_active_tetra*2)), "d_triangles malloc managed");

    // Allocate the counter for the generated triangles
    int *d_num_generated_triangles;
    print_cuda_error(cudaMallocManaged(&d_num_generated_triangles, sizeof(int)), "d_num_generated_triangles malloc managed");
    print_cuda_error(cudaMemset(d_num_generated_triangles, 0, sizeof(int)), "d_num_generated_triangles memset managed");

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


    // Take start time for triangle kernel
    cudaEvent_t tri_start, tri_stop;
    float tri_elapsedTime = 0.0f;
    cudaEventCreate(&tri_start);
    cudaEventCreate(&tri_stop);
    cudaEventRecord(tri_start, 0);

    make_triangles<<<n_blocks, n_threads>>>(n_active_tetra,
                                            d_num_generated_triangles,
                                            d_active_tetrahedra,
                                            d_triangles,
                                            threshold,
                                            d_pairs,
                                            d_7val,
                                            d_0val,
                                            d_act_val_vec);

    print_cuda_error(cudaDeviceSynchronize(), "cuda device synch");

    // Take end time for triangle kernel
    cudaEventRecord(tri_stop, 0);
    cudaEventSynchronize(tri_stop);
    cudaEventElapsedTime(&tri_elapsedTime, tri_start, tri_stop);

    cudaEventDestroy(tri_start);
    cudaEventDestroy(tri_stop);

    printf("Generate triangle kernel time: %.3f ms\n", tri_elapsedTime);

    printf("Total # of triangles:                           %d\n", n_active_tetra);

    cudaError_t last_err = cudaGetLastError();
    if (last_err != cudaSuccess) {
        printf("CUDA error after make_triangles: %s\n", cudaGetErrorString(last_err));
    }

    // WRITE TO FILE THE TRIANGLES
    int final_triangles_number = *d_num_generated_triangles;

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

    // printf("Wrote %d triangles to triangle.pdb\n", final_triangles_number);

    *time += tri_elapsedTime;
}

void load_to_const_tetra(Dimensions *dimensions){

    printf("Dimensions: x=%d, y=%d, z=%d\n", dimensions->x_dim, dimensions->y_dim, dimensions->z_dim);
    int3 tmp;
    tmp.x = dimensions->x_dim;
    tmp.y = dimensions->y_dim;
    tmp.z = dimensions->z_dim;

    cudaMemcpyToSymbol(c_dim, &tmp, sizeof(int3));

    printf("c_dim copied to const mem\n");
}

// MIXED VERSION

void parallel_march_tetra_mixed(dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                                int number_relevant_cubes,
                                cube_gpu_SoA *d_relevant_cubes, float4 *d_active_tetrahedra,
                                int* act_val_vec, int *pairs, Triangle_GPU **triangles, int *total_triangles,
                                double *time)
    {               

    //      //      //      // GENERAL INFO KERNEL APEX //     //      //      //

    int n_threads = 1024;
    int n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    printf("\nLaunching kernel to write apex\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND MANAGE //     //      //      //

    cube_vertices_points_SoA *d_cube_points_coordinates;
    print_cuda_error(cudaMallocManaged(&d_cube_points_coordinates, sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8), "cudaMalloc failed for d_cube_points_coordinates");

    printf("Allocating %zu bytes for d_cube_points_coordinates\n", sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8);

    printf("Number of cube points coordinates:  %d\n", (number_relevant_cubes) * 8);

    // Take the times for the kernel compute apex
    cudaEvent_t apex_start, apex_stop;
    float apex_elapsedTime;
    cudaEventCreate(&apex_start);
    cudaEventCreate(&apex_stop);
    cudaEventRecord(apex_start, 0);

    compute_apex_float4<<<n_blocks, n_threads>>>(d_grid, d_relevant_cubes, number_relevant_cubes, d_cube_points_coordinates);
    cudaDeviceSynchronize();

    cudaEventRecord(apex_stop, 0);
    cudaEventSynchronize(apex_stop);
    cudaEventElapsedTime(&apex_elapsedTime, apex_start, apex_stop);

    printf("Compute_apex kernel execution time: %f ms\n", apex_elapsedTime);
    *time += apex_elapsedTime;
    
    print_cuda_error(cudaGetLastError(), "CUDA error in compute_apex: %s\n");

    cudaEventDestroy(apex_start);
    cudaEventDestroy(apex_stop);
    cudaDeviceSynchronize();


    printf("Number of relevant cubes: %d\n", number_relevant_cubes);

    // KENREL TO INTERFACE

    int n_tetra = 5;
    n_threads = 1024;
    n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    int n_active_tetra = number_relevant_cubes * n_tetra;
    print_cuda_error(cudaMallocManaged(&d_active_tetrahedra, sizeof(float4) * n_active_tetra * 4), "d_active_tetrahedra malloc managed");

    int *d_cube_decomposition;
    print_cuda_error(cudaMallocManaged(&d_cube_decomposition, sizeof(int) * 20), "d_cube_decomposition malloc managed");
    print_cuda_error(cudaMemcpy(d_cube_decomposition, cube_decomposition, sizeof(int) * 20, cudaMemcpyHostToDevice), "d_cube_decomposition memcpy");
    // Take the times for the from_cubes_to_tetra kernel
    cudaEvent_t tetra_start, tetra_stop;
    float tetra_elapsedTime = 0.0f;
    cudaEventCreate(&tetra_start);
    cudaEventCreate(&tetra_stop);
    cudaEventRecord(tetra_start, 0);

    from_cubes_to_tetra<<<n_blocks, n_threads>>>(number_relevant_cubes, d_cube_points_coordinates, d_active_tetrahedra, d_cube_decomposition);
    cudaDeviceSynchronize();

    cudaEventRecord(tetra_stop, 0);
    cudaEventSynchronize(tetra_stop);
    cudaEventElapsedTime(&tetra_elapsedTime, tetra_start, tetra_stop);

    cudaEventDestroy(tetra_start);
    cudaEventDestroy(tetra_stop);

    printf("from_cubes_to_tetra kernel execution time: %f ms\n", tetra_elapsedTime);
    *time += tetra_elapsedTime;

    
    // printf("Printing d_active_tetrahedra:\n");
    // for (int i = 0; i < n_active_tetra; ++i) {

    //     float x = d_active_tetrahedra[i].x;
    //     float y = d_active_tetrahedra[i].y;
    //     float z = d_active_tetrahedra[i].z;
    //     float w = d_active_tetrahedra[i].w;

    //     printf("(%f, %f, %f), value=%f\n", x, y, z, w);
    // }

    // FILE *fp_dump = fopen("d_active_tetrahedra.bin", "wb");
    // if (!fp_dump) {
    //     printf("Failed to open d_active_tetrahedra.bin for writing\n");
    // } else {
    //     size_t written = fwrite(d_active_tetrahedra, sizeof(float4), n_active_tetra * 4, fp_dump);
    //     if (written != n_active_tetra * 4) {
    //         printf("Warning: wrote only %zu of %d float4s\n", written, n_active_tetra * 4);
    //     }
    //     fclose(fp_dump);
    //     printf("Dumped %d float4s to d_active_tetrahedra.bin\n", n_active_tetra * 4);
    // }
        
    /*  Here I have:
        number of relevant cubes
        cubes points coordinates and values
    */

    /*  To use make_triangles I need
        
        number of active tetra = number of relevant cubes * 5
        pointer to var to store the # of triangles
        pointer to store the active tetra
        pointer to store the triangles
        threshold
        pointer to pairs
        pointer to 7val
        pointer to 0val
        pointer to d_act_val
    
    */

    //      //      //      // GENERAL INFO KERNEL MT //        //      //      //

    n_threads = 512;
    n_blocks = (n_active_tetra + n_threads - 1) / n_threads;

    // Allocate the Triangles GPU
    Triangle_GPU *d_triangles;
    print_cuda_error(cudaMallocManaged(&d_triangles, sizeof(Triangle_GPU) * (n_active_tetra*2)), "d_triangles malloc managed");

    // Allocate the counter for the generated triangles
    int *d_num_generated_triangles;
    print_cuda_error(cudaMallocManaged(&d_num_generated_triangles, sizeof(int)), "d_num_generated_triangles malloc managed");
    print_cuda_error(cudaMemset(d_num_generated_triangles, 0, sizeof(int)), "d_num_generated_triangles memset managed");

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


    // Take start time for triangle kernel
    cudaEvent_t tri_start, tri_stop;
    float tri_elapsedTime = 0.0f;
    cudaEventCreate(&tri_start);
    cudaEventCreate(&tri_stop);
    cudaEventRecord(tri_start, 0);

    make_triangles<<<n_blocks, n_threads>>>(n_active_tetra,
                                            d_num_generated_triangles,
                                            d_active_tetrahedra,
                                            d_triangles,
                                            threshold,
                                            d_pairs,
                                            d_7val,
                                            d_0val,
                                            d_act_val_vec);

    print_cuda_error(cudaDeviceSynchronize(), "cuda device synch");

    // Take end time for triangle kernel
    cudaEventRecord(tri_stop, 0);
    cudaEventSynchronize(tri_stop);
    cudaEventElapsedTime(&tri_elapsedTime, tri_start, tri_stop);

    cudaEventDestroy(tri_start);
    cudaEventDestroy(tri_stop);

    printf("Generate triangle kernel time: %.3f ms\n", tri_elapsedTime);

    printf("Total # of triangles:                           %d\n", *d_num_generated_triangles);

    cudaError_t last_err = cudaGetLastError();
    if (last_err != cudaSuccess) {
        printf("CUDA error after make_triangles: %s\n", cudaGetErrorString(last_err));
    }

    // WRITE TO FILE THE TRIANGLES
    int final_triangles_number = *d_num_generated_triangles;

    int local_counter  = 0;
    
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

    *time += tri_elapsedTime;
}

// void print_after_prprocessing_using_tetra(){
//      // WRITE TO FILE THE DOMAIN DECOMP
//     // int local_counter  = 0;
//     // FILE *fp = fopen("tetra.pdb", "w");
//     // if (!fp) {
//     //     printf("Failed to open tetra.pdb for writing\n");
//     // } else {
//     //     for (int i = 4500000; i < + 4500000 + 100000 -4; i = i+4) {

//     //         int id1 = i;
//     //         int id2 = i+1;
//     //         int id3 = i+2;
//     //         int id4 = i+3;

//     //         int id1_ = local_counter;
//     //         int id2_ = local_counter+1;
//     //         int id3_ = local_counter+2;
//     //         int id4_ = local_counter+3;

//     //         char atom_type1 = (d_results[id1].w < threshold) ? 'C' : 'O';
//     //         char atom_type2 = (d_results[id2].w < threshold) ? 'C' : 'O';
//     //         char atom_type3 = (d_results[id3].w < threshold) ? 'C' : 'O';
//     //         char atom_type4 = (d_results[id4].w < threshold) ? 'C' : 'O';

//     //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
//     //                 id1_, atom_type1, d_results[id1].x, d_results[id1].y, d_results[id1].z, atom_type1);
//     //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
//     //                 id2_, atom_type2, d_results[id2].x, d_results[id2].y, d_results[id2].z, atom_type2);
//     //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
//     //                 id3_, atom_type3, d_results[id3].x, d_results[id3].y, d_results[id3].z, atom_type3);
//     //         fprintf(fp, "ATOM  %5d  %c   TET A   1    %8.3f%8.3f%8.3f  1.00  0.00           %c\n",
//     //                 id4_, atom_type4, d_results[id4].x, d_results[id4].y, d_results[id4].z, atom_type4);

//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id1_, id2_);
//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id1_, id3_);
//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id1_, id4_);

//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id2_, id3_);
//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id2_, id4_);

//     //         fprintf(fp, "CONECT%5d%5d\n",
//     //                 id3_, id4_);

//     //         local_counter+=4;
//     //     }
//     //     fclose(fp);
//     //     printf("Wrote %d atoms to tetra.pdb\n", (*d_num_generated_triangles));
//     // }
// }

void remove_unnecessary_cubes_mixed(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                int *number_relevant_cubes,
                                cube_gpu_SoA **d_relevant_cubes,
                                double *time)
{
    // //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

    int n_threads = 1024;
    int n_blocks = (cubes_in_domain + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to remove unnecessary cubes\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

    // Number of relevant cubes (the ones that are not all in or all out)
    int *d_number_relevant_cubes;
    print_cuda_error(cudaMallocManaged(&d_number_relevant_cubes, sizeof(int)), "cudaMallocManaged failed for d_number_relevant_cubes");
    *d_number_relevant_cubes = 0;    // Initialize to 0

    float4 *d_coord_idx;
    float4 *d_one_apex;
    float4 *d_two_apex;

    print_cuda_error(cudaMallocManaged(d_relevant_cubes, sizeof(cube_gpu_SoA)), "cudaMallocManaged failed for d_relevant_cubes:");

    print_cuda_error(cudaMallocManaged(&d_coord_idx, sizeof(float4)*cubes_in_domain), "cudaMallocManaged failed for d_coord_idx:");
    print_cuda_error(cudaMallocManaged(&d_one_apex, sizeof(float4)*cubes_in_domain), "cudaMallocManaged failed for d_one_apex:");
    print_cuda_error(cudaMallocManaged(&d_two_apex, sizeof(float4)*cubes_in_domain), "cudaMallocManaged failed for d_two_apex:");

    (*d_relevant_cubes)->coord_idx = d_coord_idx;
    (*d_relevant_cubes)->one_apex = d_one_apex;
    (*d_relevant_cubes)->two_apex = d_two_apex;
    
    printf("HERE\n");

    // Setup time reader for the kernel
    cudaEvent_t start, stop;
    float elapsedTime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    remove_unnecessary_cubes_mixed<<<n_blocks, n_threads, (n_threads + 1) * 4 * sizeof(float)>>>(   d_grid, d_number_relevant_cubes,
                                                                cubes_in_domain, threshold,
                                                                *d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);
    *time = elapsedTime;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaDeviceSynchronize();

    *number_relevant_cubes = *d_number_relevant_cubes;

    printf("Number of relevant cubes:           %d\n", (*number_relevant_cubes));
    printf("Total number of cubes:              %zu\n", cubes_in_domain);
    
    // Take the potential error
    print_cuda_error(cudaGetLastError(), "CUDA error");
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");

    // print_relevant_points_soa(*d_relevant_cubes, number_relevant_cubes);
}