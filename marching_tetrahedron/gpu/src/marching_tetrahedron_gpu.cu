#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__global__ void remove_unnecessary_cubes_kernel(double* grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu* d_relevant_cubes) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;
    
    int DZ = dim->z_dim;
    int DY = dim->y_dim;
    int DX = dim->x_dim;

    int i = idx / (DZ * DY);

    int rem = idx % (DZ * DY);

    int j = rem / DZ;
    int k = rem % DZ;

    // if(idx == 39360)
    //     printf("2 idx=%d, i=%d, j=%d, k=%d, rem=%d\n", idx, i, j, k, rem);
    // if(idx == 39360)
    //     printf("3 idx=%d, i=%d, j=%d, k=%d, rem=%d\n", idx, i, j, k, rem);
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
        d_relevant_cubes[insert_pos].idx = idx;
        d_relevant_cubes[insert_pos].x = i;
        d_relevant_cubes[insert_pos].y = j;
        d_relevant_cubes[insert_pos].z = k;
    }

}

__global__ void skip_preprocessing_k(double* grid,size_t size, double threshold,Dimensions *dim, cube_gpu* d_relevant_cubes) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;
    
    int DZ = dim->z_dim;
    int DY = dim->y_dim;

    int i = idx / (DZ * DY);

    int rem = idx % (DZ * DY);

    int j = rem / DZ;
    int k = rem % DZ;

    d_relevant_cubes[idx].idx = idx;
    d_relevant_cubes[idx].x = i;
    d_relevant_cubes[idx].y = j;
    d_relevant_cubes[idx].z = k;
}

__global__ void compute_apex(   double *grid, cube_gpu *d_relevant_cubes, int number_relevant_cubes,
                                cube_vertices_points *d_cube_points_coordinates, Dimensions *dim){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes[tid].x >= dim->x_dim-1 ||
        d_relevant_cubes[tid].y >= dim->y_dim-1 ||
        d_relevant_cubes[tid].z >= dim->z_dim-1 ){
            return;
        }
    
    if ((d_relevant_cubes[tid].x & 1) == 0)
    {
        d_relevant_cubes[tid].one_apex[0] = d_relevant_cubes[tid].x;
    }
    else
    {
        d_relevant_cubes[tid].one_apex[0] = d_relevant_cubes[tid].x + 1;
    }
    
    if ((d_relevant_cubes[tid].y & 1) == 0)
    {
        d_relevant_cubes[tid].one_apex[1] = d_relevant_cubes[tid].y;
    }
    else
    {
        d_relevant_cubes[tid].one_apex[1] = d_relevant_cubes[tid].y + 1;
    }
    
    if ((d_relevant_cubes[tid].z & 1) == 0)
    {
        d_relevant_cubes[tid].one_apex[2] = d_relevant_cubes[tid].z;
    }
    else
    {
        d_relevant_cubes[tid].one_apex[2] = d_relevant_cubes[tid].z + 1;
    }
    
    d_relevant_cubes[tid].two_apex[0] = 2 * d_relevant_cubes[tid].x + 1 - d_relevant_cubes[tid].one_apex[0];
    d_relevant_cubes[tid].two_apex[1] = 2 * d_relevant_cubes[tid].y + 1 - d_relevant_cubes[tid].one_apex[1];
    d_relevant_cubes[tid].two_apex[2] = 2 * d_relevant_cubes[tid].z + 1 - d_relevant_cubes[tid].one_apex[2];

    d_cube_points_coordinates[tid * 8 + 0].coord.x = d_relevant_cubes[tid].one_apex[0];
    d_cube_points_coordinates[tid * 8 + 0].coord.y = d_relevant_cubes[tid].one_apex[1];
    d_cube_points_coordinates[tid * 8 + 0].coord.z = d_relevant_cubes[tid].one_apex[2];

    d_cube_points_coordinates[tid * 8 + 1].coord.x = d_relevant_cubes[tid].two_apex[0];
    d_cube_points_coordinates[tid * 8 + 1].coord.y = d_relevant_cubes[tid].one_apex[1];
    d_cube_points_coordinates[tid * 8 + 1].coord.z = d_relevant_cubes[tid].one_apex[2];

    d_cube_points_coordinates[tid * 8 + 2].coord.x = d_relevant_cubes[tid].one_apex[0];
    d_cube_points_coordinates[tid * 8 + 2].coord.y = d_relevant_cubes[tid].two_apex[1];
    d_cube_points_coordinates[tid * 8 + 2].coord.z = d_relevant_cubes[tid].one_apex[2];

    d_cube_points_coordinates[tid * 8 + 3].coord.x = d_relevant_cubes[tid].two_apex[0];
    d_cube_points_coordinates[tid * 8 + 3].coord.y = d_relevant_cubes[tid].two_apex[1];
    d_cube_points_coordinates[tid * 8 + 3].coord.z = d_relevant_cubes[tid].one_apex[2];

    d_cube_points_coordinates[tid * 8 + 4].coord.x = d_relevant_cubes[tid].one_apex[0];
    d_cube_points_coordinates[tid * 8 + 4].coord.y = d_relevant_cubes[tid].one_apex[1];
    d_cube_points_coordinates[tid * 8 + 4].coord.z = d_relevant_cubes[tid].two_apex[2];

    d_cube_points_coordinates[tid * 8 + 5].coord.x = d_relevant_cubes[tid].two_apex[0];
    d_cube_points_coordinates[tid * 8 + 5].coord.y = d_relevant_cubes[tid].one_apex[1];
    d_cube_points_coordinates[tid * 8 + 5].coord.z = d_relevant_cubes[tid].two_apex[2];

    d_cube_points_coordinates[tid * 8 + 6].coord.x = d_relevant_cubes[tid].one_apex[0];
    d_cube_points_coordinates[tid * 8 + 6].coord.y = d_relevant_cubes[tid].two_apex[1];
    d_cube_points_coordinates[tid * 8 + 6].coord.z = d_relevant_cubes[tid].two_apex[2];

    d_cube_points_coordinates[tid * 8 + 7].coord.x = d_relevant_cubes[tid].two_apex[0];
    d_cube_points_coordinates[tid * 8 + 7].coord.y = d_relevant_cubes[tid].two_apex[1];
    d_cube_points_coordinates[tid * 8 + 7].coord.z = d_relevant_cubes[tid].two_apex[2];

    
    for(int v=0; v<8; v++){
        d_cube_points_coordinates[tid * 8 + v].value = grid[d_cube_points_coordinates[tid * 8 + v].coord.z +
                                                            d_cube_points_coordinates[tid * 8 + v].coord.y * dim->z_dim +
                                                            d_cube_points_coordinates[tid * 8 + v].coord.x * dim->z_dim * dim->y_dim];
    }
}

__global__ void compute_march_tetra(    double *d_grid, cube_gpu *d_relevant_cubes,
                                        int number_relevant_cubes, int *cube_deco,
                                        cube_vertices_points *d_cube_points_coordinates,
                                        cube_vertices_points *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter, Dimensions *dim){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    int  bello = 0;


    if( d_relevant_cubes[tid].x >= dim->x_dim-1 ||
        d_relevant_cubes[tid].y >= dim->y_dim-1 ||
        d_relevant_cubes[tid].z >= dim->z_dim-1 ){
            return;
        }

    for (int tetra = 0; tetra < 5; tetra++){

        cube_vertices_points *first     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+0]-1];
        cube_vertices_points *second    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+1]-1];
        cube_vertices_points *third     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+2]-1];
        cube_vertices_points *fourth    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+3]-1];

        sort_points(&first, &second, &third, &fourth);
        
        int less = 0, eq = 0, gre = 0;
        
        count_elements(&less, &eq, &gre, first, second, third, fourth, threshold);
        
        int act_val = act_val_vec[less + eq*5];
        // int act_val2 = get_action_value(less, eq, gre);
        
        if (act_val!= 0){
            int *use_pairs = &d_pairs[6*(act_val-1)];

            make_triangle(first, second, third, fourth, &d_triangles[tid*5+tetra], use_pairs, d_relevant_cubes[tid].idx == bello);

            if(act_val == 7){
                use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 
                int insert_pos = atomicAdd(d_counter, 1);

                make_triangle(first, second, third, fourth, &d_triangles[number_relevant_cubes*5 + insert_pos], use_pairs, d_relevant_cubes[tid].idx == bello);
            }
        }
    }
}

__device__ void make_triangle(  cube_vertices_points *first, cube_vertices_points *second,
                                cube_vertices_points *third, cube_vertices_points *fourth,
                                Triangle_GPU *triangle, int *pairs, bool debug){

    cube_vertices_points* arr[4] = {first, second, third, fourth};
    
    int idx1 = pairs[0] - 1;
    int idx2 = pairs[1] - 1;

    triangle->v1.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v1.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v1.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;

    if (debug) {
        printf("Triangle v1: (%f, %f, %f)\n", triangle->v1.x, triangle->v1.y, triangle->v1.z);
    }

    idx1 = pairs[2] - 1;
    idx2 = pairs[3] - 1;

    triangle->v2.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v2.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v2.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;

    if (debug) {
        printf("Triangle v2: (%f, %f, %f)\n", triangle->v2.x, triangle->v2.y, triangle->v2.z);
    }

    idx1 = pairs[4] - 1;
    idx2 = pairs[5] - 1;
    
    triangle->v3.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v3.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v3.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;

    if (debug) {
        printf("Triangle v3: (%f, %f, %f)\n\n", triangle->v3.x, triangle->v3.y, triangle->v3.z);
    }
}

__device__ int get_action_value( int less, int eq, int gre){
    
    if (less == 0 || (less == 2 && eq == 2) || (less == 3 && eq == 1) || less == 4)
    {
        return 0;
    }

    if (less == 3)
        return 6;

    if (less == 2)
    {
        if (eq == 0)
            return 7;
        if (eq == 1)
            return 5;
    }

    if (less == 1)
    {
        if (eq == 0)
            return 1;
        if (eq == 1)
            return 2;
        if (eq == 2)
            return 3;
        if (eq == 3)
            return 4;
    }
}

__device__ void count_elements( int *less, int *eq, int *gre, cube_vertices_points *first,
                                cube_vertices_points *second, cube_vertices_points *third, cube_vertices_points *fourth,
                                dim_t threshold){

    cube_vertices_points* arr[4] = {first, second, third, fourth};
    for (int i = 0; i < 4; ++i) {
        if (arr[i]->value < threshold) {
            (*less)++;
        } else if (arr[i]->value == threshold) {
            (*eq)++;
        } else {
            (*gre)++;
        }
    }
}

// TODO SISTEMA IL DOPPIO ASTERISCO
__device__ void sort_points(cube_vertices_points **first, cube_vertices_points **second, cube_vertices_points **third, cube_vertices_points **fourth){

    cube_vertices_points* arr[4] = {*first, *second, *third, *fourth};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3 - i; ++j) {
            if (arr[j]->value > arr[j + 1]->value) {
                cube_vertices_points* tmp = arr[j];
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



void remove_unnecessary_cubes(  double *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates)
{
    //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

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
    *d_dim = *dim;      // Copy the DS from the serial code

    print_cuda_error(cudaMallocManaged(d_relevant_cubes, sizeof(cube_gpu)*cubes_in_domain), "cudaMallocManaged failed for d_relevant_cubes:");
    
    
    // Setup time reader for the kernel
    cudaEvent_t start, stop;
    float elapsedTime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    remove_unnecessary_cubes_kernel<<<n_blocks, n_threads>>>(   d_grid, d_number_relevant_cubes,
                                                                cubes_in_domain, threshold, d_dim,
                                                                *d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaDeviceSynchronize();

    *number_relevant_cubes = *d_number_relevant_cubes;

    printf("Number of relevant cubes:           %d\n", (*number_relevant_cubes));
    printf("Total number of cubes:              %zu\n", cubes_in_domain);
    
    // Take the potential error
    print_cuda_error(cudaGetLastError(), "CUDA error");
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");

    print_relevant_points(*d_relevant_cubes, number_relevant_cubes);
}

void parallel_march_tetra   (Dimensions *dim, dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            size_t *triangle_counter, size_t *vertex_counter, int number_relevant_cubes,
                            cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates,
                            int* act_val_vec, int *pairs, Triangle_GPU **triangles)
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

    print_cuda_error(cudaMallocManaged(d_cube_points_coordinates, sizeof(cube_vertices_points)*(number_relevant_cubes)*8), "cudaMalloc failed for d_cube_points_coordinates");

    printf("Allocating %zu bytes for d_cube_points_coordinates\n", sizeof(cube_vertices_points)*(number_relevant_cubes)*8);

    printf("Number of cube points coordinates:  %d\n", (number_relevant_cubes) * 8);

    // Take the times for the kernel compute apex
    cudaEvent_t apex_start, apex_stop;
    float apex_elapsedTime;
    cudaEventCreate(&apex_start);
    cudaEventCreate(&apex_stop);
    cudaEventRecord(apex_start, 0);

    compute_apex<<<n_blocks, n_threads>>>(d_grid, *d_relevant_cubes, number_relevant_cubes, *d_cube_points_coordinates, d_dim);
    cudaDeviceSynchronize();

    cudaEventRecord(apex_stop, 0);
    cudaEventSynchronize(apex_stop);
    cudaEventElapsedTime(&apex_elapsedTime, apex_start, apex_stop);

    printf("Compute_apex kernel execution time: %f ms\n", apex_elapsedTime);

    print_cuda_error(cudaGetLastError(), "CUDA error in compute_apex: %s\n");

    cudaEventDestroy(apex_start);
    cudaEventDestroy(apex_stop);
    cudaDeviceSynchronize();
    //      //      //      // GENERAL INFO KERNEL MT //        //      //      //

    n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    
    printf("\nLaunching kernel to compute MT algo\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);
    
    //      //      //      // MALLOC AND MANAGE //     //      //      //

    (*vertex_counter) = 0;
    (*triangle_counter) = 0;

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

    cube_vertices_points *d_stack_pool;
    print_cuda_error(cudaMallocManaged(&d_stack_pool, sizeof(cube_vertices_points)*20*number_relevant_cubes), "cudaMallocManaged failed for d_stack_pool: %s\n");
    
    int *pool_index;
    print_cuda_error(cudaMallocManaged(&pool_index, sizeof(int)), "cudaMallocManaged failed for pool_index: %s\n");

    
    // Take the times
    cudaEvent_t mt_start, mt_stop;
    float mt_elapsedTime;
    cudaEventCreate(&mt_start);
    cudaEventCreate(&mt_stop);
    cudaEventRecord(mt_start, 0);

    compute_march_tetra<<<n_blocks, n_threads>>>(   d_grid, *d_relevant_cubes,
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
    int total_triangles = number_relevant_cubes * 5 + *d_counter; // 5 tetrahedra per cube plus possible extra triangles
    printf("Total triangles:                    %d\n", total_triangles);
    if (*triangles == NULL) {
        *triangles = (Triangle_GPU*)malloc(sizeof(Triangle_GPU) * total_triangles);
        if (*triangles == NULL) {
            fprintf(stderr, "Failed to allocate host memory for triangles.\n");
            // Prevent further use of *triangles if allocation failed
            return;
        }
        printf("Host allocation completed\n");
    }

    cudaError_t memcpy_err = cudaMemcpy(*triangles, d_triangles, sizeof(Triangle_GPU) * total_triangles, cudaMemcpyDeviceToHost);
    if (memcpy_err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed for triangles: %s\n", cudaGetErrorString(memcpy_err));
    }

    cudaError_t mt_err = cudaGetLastError();
    if (mt_err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error in compute_march_tetra: %s\n", cudaGetErrorString(mt_err));
    }
    printf("compute_march_tetra k exe time:     %f ms\n", mt_elapsedTime);

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

void skip_preprocessing(  double *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates)
{
    //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

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
    *d_dim = *dim;      // Copy the DS from the serial code

    print_cuda_error(cudaMallocManaged(d_relevant_cubes, sizeof(cube_gpu)*cubes_in_domain), "cudaMallocManaged failed for d_relevant_cubes:");
    
    
    // Setup time reader for the kernel
    cudaEvent_t start, stop;
    float elapsedTime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    skip_preprocessing_k<<<n_blocks, n_threads>>>(   d_grid, cubes_in_domain, threshold, d_dim,*d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaDeviceSynchronize();

    *number_relevant_cubes = *d_number_relevant_cubes;

    printf("Number of relevant cubes:           %d\n", (*number_relevant_cubes));
    printf("Total number of cubes:              %zu\n", cubes_in_domain);
    
    // Take the potential error
    print_cuda_error(cudaGetLastError(), "CUDA error");
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");

    print_relevant_points(*d_relevant_cubes, number_relevant_cubes);
}

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

void print_relevant_points(cube_gpu *d_relevant_cubes, int *number_relevant_cubes){


    // Remove "points" directory if it exists, then create it and change into it
    struct stat st = {0};
    if (stat("points", &st) == 0) {
        // Directory exists, remove all .pdb files inside
        system("rm -rf points");
    }
    mkdir("points", 0700);
    chdir("points");

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

void print_triangles(   Triangle_GPU *triangles, int *number_relevant_cubes,
                        char *molecule_name, char *molecule_path){
    struct stat st;

    const char *folder = "../results/";
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

    FILE *pdb_file;
    for (int i = 0; i < *number_relevant_cubes; ++i) {
        
        if(i%33333 == 0){
            char file_name[256];
            strcpy(file_name, molecule_path);
            strcat(file_name, molecule_name);
            sprintf(file_name + strlen(file_name), "_%d", file_number);
            strcat(file_name, ".pdb");
            // printf("Writing triangles to file: %s\n", file_name);
            pdb_file = fopen(file_name, "w");

            if (!pdb_file) {
                fprintf(stderr, "Failed to open points.pdb for writing.\n");
            } 

            file_number++;
            local_counter=0;
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

        if((i+1)%33333 == 0){
            fclose(pdb_file);
        }
    }
    printf("Relevant points written\n");
    
}


void print_cuda_error(cudaError_t err, const char* msg){
    if (err != cudaSuccess) {
        fprintf(stderr, "error failed: %s\nMsg: %s\n", cudaGetErrorString(err), msg);
    }
}

