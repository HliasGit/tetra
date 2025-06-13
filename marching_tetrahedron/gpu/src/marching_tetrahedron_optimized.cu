#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__global__ void remove_unnecessary_cubes_SoA_kernel(dim_t* grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu_SoA* d_relevant_cubes) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;
    
    int DZ = dim->z_dim;
    int DY = dim->y_dim;
    int DX = dim->x_dim;

    int i = idx / (DZ * DY);

    int rem = idx % (DZ * DY);

    int j = rem / DZ;
    int k = rem % DZ;

    bool all_in;
    bool all_out;

    if(k == DZ-1 || j == DY-1 || i == DX-1){
        all_in = 0;
        all_out = 0;

    } else {
        all_out =
            (grid[idx] < threshold) &&
            (grid[idx+1] < threshold) &&
            (grid[idx + DZ] < threshold) &&
            (grid[idx + DZ+1] < threshold) &&
            (grid[idx + DZ * DY] < threshold) &&
            (grid[idx + DZ * DY + 1] < threshold) &&
            (grid[idx + DZ * DY + DZ] < threshold) &&
            (grid[idx + DZ * DY + DZ+1] < threshold);

        all_in =
            (grid[idx] > threshold) &&
            (grid[idx+1] > threshold) &&
            (grid[idx + DZ] > threshold) &&
            (grid[idx + DZ+1] > threshold) &&
            (grid[idx + DZ * DY] > threshold) &&
            (grid[idx + DZ * DY + 1] > threshold) &&
            (grid[idx + DZ * DY + DZ] > threshold) &&
            (grid[idx + DZ * DY + DZ+1] > threshold);
    }

    if (all_out == 0 && all_in == 0){
        int insert_pos = atomicAdd(counter, 1);
        d_relevant_cubes->coord_idx[insert_pos].w = idx;
        // printf("coord_idx[%d].w = %f\n", insert_pos, d_relevant_cubes->coord_idx[insert_pos].w);
        d_relevant_cubes->coord_idx[insert_pos].x = i;
        d_relevant_cubes->coord_idx[insert_pos].y = j;
        d_relevant_cubes->coord_idx[insert_pos].z = k;
    }

}

__global__ void compute_apex_float4(   dim_t *grid, cube_gpu_SoA *d_relevant_cubes, int number_relevant_cubes,
                                cube_vertices_points_SoA *d_cube_points_coordinates, Dimensions *dim){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes->coord_idx[tid].x >= dim->x_dim-1 ||
        d_relevant_cubes->coord_idx[tid].y >= dim->y_dim-1 ||
        d_relevant_cubes->coord_idx[tid].z >= dim->z_dim-1 ){
            return;
        }

    // printf("coord_idx[%d].w = %f\n", tid, d_relevant_cubes->coord_idx[tid].w);
    
    if (((int)d_relevant_cubes->coord_idx[tid].x & 1) == 0)
    {
        d_relevant_cubes->one_apex[tid].x = d_relevant_cubes->coord_idx[tid].x;
    }
    else
    {
        d_relevant_cubes->one_apex[tid].x = d_relevant_cubes->coord_idx[tid].x + 1;
    }
    
    if (((int)d_relevant_cubes->coord_idx[tid].y & 1) == 0)
    {
        d_relevant_cubes->one_apex[tid].y = d_relevant_cubes->coord_idx[tid].y;
    }
    else
    {
        d_relevant_cubes->one_apex[tid].y = d_relevant_cubes->coord_idx[tid].y + 1;
    }
    
    if (((int)d_relevant_cubes->coord_idx[tid].z & 1) == 0)
    {
        d_relevant_cubes->one_apex[tid].z = d_relevant_cubes->coord_idx[tid].z;
    }
    else
    {
        d_relevant_cubes->one_apex[tid].z = d_relevant_cubes->coord_idx[tid].z + 1;
    }
    
    d_relevant_cubes->two_apex[tid].x = 2 * d_relevant_cubes->coord_idx[tid].x + 1 - d_relevant_cubes->one_apex[tid].x;
    d_relevant_cubes->two_apex[tid].y = 2 * d_relevant_cubes->coord_idx[tid].y + 1 - d_relevant_cubes->one_apex[tid].y;
    d_relevant_cubes->two_apex[tid].z = 2 * d_relevant_cubes->coord_idx[tid].z + 1 - d_relevant_cubes->one_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 0].val.x = d_relevant_cubes->one_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 0].val.y = d_relevant_cubes->one_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 0].val.z = d_relevant_cubes->one_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 1].val.x = d_relevant_cubes->two_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 1].val.y = d_relevant_cubes->one_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 1].val.z = d_relevant_cubes->one_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 2].val.x = d_relevant_cubes->one_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 2].val.y = d_relevant_cubes->two_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 2].val.z = d_relevant_cubes->one_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 3].val.x = d_relevant_cubes->two_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 3].val.y = d_relevant_cubes->two_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 3].val.z = d_relevant_cubes->one_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 4].val.x = d_relevant_cubes->one_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 4].val.y = d_relevant_cubes->one_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 4].val.z = d_relevant_cubes->two_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 5].val.x = d_relevant_cubes->two_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 5].val.y = d_relevant_cubes->one_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 5].val.z = d_relevant_cubes->two_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 6].val.x = d_relevant_cubes->one_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 6].val.y = d_relevant_cubes->two_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 6].val.z = d_relevant_cubes->two_apex[tid].z;

    d_cube_points_coordinates[tid * 8 + 7].val.x = d_relevant_cubes->two_apex[tid].x;
    d_cube_points_coordinates[tid * 8 + 7].val.y = d_relevant_cubes->two_apex[tid].y;
    d_cube_points_coordinates[tid * 8 + 7].val.z = d_relevant_cubes->two_apex[tid].z;

    
    for(int v=0; v<8; v++){
        d_cube_points_coordinates[tid * 8 + v].val.w = grid[(int)(d_cube_points_coordinates[tid * 8 + v].val.z +
                                                            d_cube_points_coordinates[tid * 8 + v].val.y * dim->z_dim +
                                                            d_cube_points_coordinates[tid * 8 + v].val.x * dim->z_dim * dim->y_dim)];
    }
}

__global__ void compute_march_tetra_SoA(dim_t *d_grid, cube_gpu_SoA *d_relevant_cubes,
                                        int number_relevant_cubes, int *cube_deco,
                                        cube_vertices_points_SoA *d_cube_points_coordinates,
                                        cube_vertices_points_SoA *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter, Dimensions *dim){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes->coord_idx[tid].x >= dim->x_dim-1 ||
        d_relevant_cubes->coord_idx[tid].y >= dim->y_dim-1 ||
        d_relevant_cubes->coord_idx[tid].z >= dim->z_dim-1 ){
            return;
        }

    for (int tetra = 0; tetra < 5; tetra++){

        cube_vertices_points_SoA *first     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+0]-1];
        cube_vertices_points_SoA *second    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+1]-1];
        cube_vertices_points_SoA *third     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+2]-1];
        cube_vertices_points_SoA *fourth    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+3]-1];

        sort_points_SoA(&first, &second, &third, &fourth);
        
        int less = 0, eq = 0, gre = 0;
        
        count_elements_SoA(&less, &eq, &gre, first, second, third, fourth, threshold);
        
        int act_val = act_val_vec[less + eq*5];
        // int act_val2 = get_action_value(less, eq, gre);
        
        if (act_val!= 0){
            int *use_pairs = &d_pairs[6*(act_val-1)];

            make_triangle_SoA(first, second, third, fourth, &d_triangles[tid*5+tetra], use_pairs);

            if(act_val == 7){
                use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 
                int insert_pos = atomicAdd(d_counter, 1);

                make_triangle_SoA(first, second, third, fourth, &d_triangles[number_relevant_cubes*5 + insert_pos], use_pairs);
            }
        }
    }
}

__device__ void make_triangle_SoA(  cube_vertices_points_SoA *first, cube_vertices_points_SoA *second,
                                cube_vertices_points_SoA *third, cube_vertices_points_SoA *fourth,
                                Triangle_GPU *triangle, int *pairs){

    cube_vertices_points_SoA* arr[4] = {first, second, third, fourth};
    
    int idx1 = pairs[0] - 1;
    int idx2 = pairs[1] - 1;

    triangle->v1.x = ((coord_t)arr[idx1]->val.x + (coord_t)arr[idx2]->val.x) / 2.0;
    triangle->v1.y = ((coord_t)arr[idx1]->val.y + (coord_t)arr[idx2]->val.y) / 2.0;
    triangle->v1.z = ((coord_t)arr[idx1]->val.z + (coord_t)arr[idx2]->val.z) / 2.0;

    idx1 = pairs[2] - 1;
    idx2 = pairs[3] - 1;

    triangle->v2.x = ((coord_t)arr[idx1]->val.x + (coord_t)arr[idx2]->val.x) / 2.0;
    triangle->v2.y = ((coord_t)arr[idx1]->val.y + (coord_t)arr[idx2]->val.y) / 2.0;
    triangle->v2.z = ((coord_t)arr[idx1]->val.z + (coord_t)arr[idx2]->val.z) / 2.0;

    idx1 = pairs[4] - 1;
    idx2 = pairs[5] - 1;
    
    triangle->v3.x = ((coord_t)arr[idx1]->val.x + (coord_t)arr[idx2]->val.x) / 2.0;
    triangle->v3.y = ((coord_t)arr[idx1]->val.y + (coord_t)arr[idx2]->val.y) / 2.0;
    triangle->v3.z = ((coord_t)arr[idx1]->val.z + (coord_t)arr[idx2]->val.z) / 2.0;
}

__device__ void count_elements_SoA( int *less, int *eq, int *gre, cube_vertices_points_SoA *first,
                                cube_vertices_points_SoA *second, cube_vertices_points_SoA *third, cube_vertices_points_SoA *fourth,
                                dim_t threshold){

    cube_vertices_points_SoA* arr[4] = {first, second, third, fourth};
    for (int i = 0; i < 4; ++i) {
        if (arr[i]->val.w < threshold) {
            (*less)++;
        } else if (arr[i]->val.w == threshold) {
            (*eq)++;
        } else {
            (*gre)++;
        }
    }
}

__device__ void sort_points_SoA(cube_vertices_points_SoA **first, cube_vertices_points_SoA **second, cube_vertices_points_SoA **third, cube_vertices_points_SoA **fourth){

    cube_vertices_points_SoA* arr[4] = {*first, *second, *third, *fourth};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3 - i; ++j) {
            if (arr[j]->val.w > arr[j + 1]->val.w) {
                cube_vertices_points_SoA* tmp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = tmp;
            }
        }
    }
    
    *first = arr[0];
    *second = arr[1];
    *third = arr[2];
    *fourth = arr[3];
}

void remove_unnecessary_cubes(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu_SoA **d_relevant_cubes,
                                double *time)
{
    // //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

    int n_threads = 512;
    int n_blocks = (cubes_in_domain + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to remove unnecessary cubes\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

    // Number of relevant cubes (the ones that are not all in or all out)
    int *d_number_relevant_cubes;
    print_cuda_error(cudaMallocManaged(&d_number_relevant_cubes, sizeof(int)), "cudaMallocManaged failed for d_number_relevant_cubes");
    *d_number_relevant_cubes = 0;    // Initialize to 0

    // Data structure for the dimensions
    Dimensions *d_dim;
    print_cuda_error(cudaMallocManaged(&d_dim, sizeof(Dimensions)), "cudaMallocManaged failed for d_dim: %s");
    *d_dim = *dim;

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
    
    remove_unnecessary_cubes_SoA_kernel<<<n_blocks, n_threads>>>(   d_grid, d_number_relevant_cubes,
                                                                cubes_in_domain, threshold, d_dim,
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

void parallel_march_tetra   (Dimensions *dim, dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            int number_relevant_cubes,
                            cube_gpu_SoA *d_relevant_cubes, cube_vertices_points_SoA **d_cube_points_coordinates,
                            int* act_val_vec, int *pairs, Triangle_GPU **triangles, int *total_triangles,
                            double *time)
    {               

    //      //      //      // GENERAL INFO KERNEL APEX //     //      //      //

    int n_threads = 512;
    int n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    printf("\nLaunching kernel to write apex\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND MANAGE //     //      //      //

    // Data structure for the dimensions
    Dimensions *d_dim;
    print_cuda_error(cudaMallocManaged(&d_dim, sizeof(Dimensions)), "cudaMallocManaged failed for d_dim: %s");
    *d_dim = *dim;      // Copy the DS from the serial code

    print_cuda_error(cudaMallocManaged(d_cube_points_coordinates, sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8), "cudaMalloc failed for d_cube_points_coordinates");

    printf("Allocating %zu bytes for d_cube_points_coordinates\n", sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8);

    printf("Number of cube points coordinates:  %d\n", (number_relevant_cubes) * 8);

    // Take the times for the kernel compute apex
    cudaEvent_t apex_start, apex_stop;
    float apex_elapsedTime;
    cudaEventCreate(&apex_start);
    cudaEventCreate(&apex_stop);
    cudaEventRecord(apex_start, 0);

    compute_apex_float4<<<n_blocks, n_threads>>>(d_grid, d_relevant_cubes, number_relevant_cubes, *d_cube_points_coordinates, d_dim);
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

    //      //      //      // GENERAL INFO KERNEL MT //        //      //      //

    n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    
    printf("\nLaunching kernel to compute MT algo\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);
    
    //      //      //      // MALLOC AND MANAGE //     //      //      //

    int *d_cube_decomposition;
    print_cuda_error(cudaMalloc(&d_cube_decomposition, sizeof(int)*20), "cudaMalloc failed for d_cube_decomposition: %s\n");
    print_cuda_error(cudaMemcpy(d_cube_decomposition, cube_decomposition, sizeof(int) * 20, cudaMemcpyHostToDevice), "cudaMemcpy failed for d_cube_decomposition: %s\n");
    
    int *d_act_val_vec;
    print_cuda_error(cudaMalloc(&d_act_val_vec, sizeof(int)*25), "cudaMalloc failed for d_act_val_vec alloc: %s\n");
    print_cuda_error(cudaMemcpy(d_act_val_vec, act_val_vec, sizeof(int) * 25, cudaMemcpyHostToDevice), "cudaMemcpy failed for d_act_val_vec copy: %s\n");

    int *d_pairs;
    print_cuda_error(cudaMalloc(&d_pairs, sizeof(int)*48), "cudaMalloc failed for d_pairs: %s\n");
    print_cuda_error(cudaMemcpy(d_pairs, pairs, sizeof(int) * 48, cudaMemcpyHostToDevice), "cudaMemcpy failed for d_pairs: %s\n");
    
    Triangle_GPU *d_triangles;
    print_cuda_error(cudaMalloc(&d_triangles, sizeof(Triangle_GPU)*12*number_relevant_cubes), "cudaMalloc failed for d_triangles: %s\n");
    
    int *d_counter;
    print_cuda_error(cudaMallocManaged(&d_counter, sizeof(int)), "cudaMallocManaged failed for d_counter: %s\n");
    *d_counter = 0;

    cube_vertices_points_SoA *d_stack_pool;
    print_cuda_error(cudaMallocManaged(&d_stack_pool, sizeof(cube_vertices_points_SoA)*20*number_relevant_cubes), "cudaMallocManaged failed for d_stack_pool: %s\n");
    
    int *pool_index;
    print_cuda_error(cudaMallocManaged(&pool_index, sizeof(int)), "cudaMallocManaged failed for pool_index: %s\n");

    
    // Take the times
    cudaEvent_t mt_start, mt_stop;
    float mt_elapsedTime;
    cudaEventCreate(&mt_start);
    cudaEventCreate(&mt_stop);
    cudaEventRecord(mt_start, 0);

    compute_march_tetra_SoA<<<n_blocks, n_threads>>>(   d_grid, d_relevant_cubes,
                                                    number_relevant_cubes, d_cube_decomposition,
                                                    *d_cube_points_coordinates, d_stack_pool, pool_index,
                                                    threshold, d_act_val_vec, d_pairs,
                                                    d_triangles, d_counter, d_dim);
    cudaDeviceSynchronize();

    printf("Kernel MT terminated\n");
    printf("7 act-val:                          %d\n", *d_counter);
    
    cudaEventRecord(mt_stop, 0);
    cudaEventSynchronize(mt_stop);
    cudaEventElapsedTime(&mt_elapsedTime, mt_start, mt_stop);
    
    // Ensure host memory is allocated for triangles using the right number of elements
    *total_triangles = number_relevant_cubes * 5 + *d_counter; // 5 tetrahedra per cube plus possible extra triangles
    printf("Total triangles:                    %d\n", (*total_triangles));
    if (*triangles == NULL) {
        *triangles = (Triangle_GPU*)malloc(sizeof(Triangle_GPU) * (*total_triangles));
        if (*triangles == NULL) {
            fprintf(stderr, "Failed to allocate host memory for triangles.\n");
            // Prevent further use of *triangles if allocation failed
            return;
        }
        printf("Host allocation completed\n");
    }

    cudaError_t memcpy_err = cudaMemcpy(*triangles, d_triangles, sizeof(Triangle_GPU) * (*total_triangles), cudaMemcpyDeviceToHost);
    if (memcpy_err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed for triangles: %s\n", cudaGetErrorString(memcpy_err));
    }

    cudaError_t mt_err = cudaGetLastError();
    if (mt_err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error in compute_march_tetra: %s\n", cudaGetErrorString(mt_err));
    }
    printf("compute_march_tetra k exe time:     %f ms\n", mt_elapsedTime);
    *time += mt_elapsedTime;

    cudaEventDestroy(mt_start);
    cudaEventDestroy(mt_stop);
    cudaDeviceSynchronize();

    cudaFree(d_act_val_vec);
    cudaFree(d_pairs);
    cudaFree(d_triangles);
    cudaFree(d_cube_decomposition);
    cudaFree(d_counter);
    cudaFree(d_stack_pool);
    cudaFree(pool_index);
}