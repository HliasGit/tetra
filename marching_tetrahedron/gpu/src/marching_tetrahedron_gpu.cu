#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__global__ void remove_unnecessary_cubes_kernel(dim_t* grid, int *counter,
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

__global__ void skip_preprocessing_k(size_t size, Dimensions *dim, cube_gpu* d_relevant_cubes) {

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

__global__ void compute_apex(   dim_t *grid, cube_gpu *d_relevant_cubes, int number_relevant_cubes,
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

__global__ void compute_march_tetra(    dim_t *d_grid, cube_gpu *d_relevant_cubes,
                                        int number_relevant_cubes, int *cube_deco,
                                        cube_vertices_points *d_cube_points_coordinates,
                                        cube_vertices_points *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter, Dimensions *dim){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

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

            make_triangle(first, second, third, fourth, &d_triangles[tid*5+tetra], use_pairs);

            if(act_val == 7){
                use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 
                int insert_pos = atomicAdd(d_counter, 1);

                make_triangle(first, second, third, fourth, &d_triangles[number_relevant_cubes*5 + insert_pos], use_pairs);
            }
        }
    }
}

__device__ void make_triangle(  cube_vertices_points *first, cube_vertices_points *second,
                                cube_vertices_points *third, cube_vertices_points *fourth,
                                Triangle_GPU *triangle, int *pairs){

    cube_vertices_points* arr[4] = {first, second, third, fourth};
    
    int idx1 = pairs[0] - 1;
    int idx2 = pairs[1] - 1;

    triangle->v1.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v1.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v1.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;

    idx1 = pairs[2] - 1;
    idx2 = pairs[3] - 1;

    triangle->v2.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v2.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v2.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;

    idx1 = pairs[4] - 1;
    idx2 = pairs[5] - 1;
    
    triangle->v3.x = ((coord_t)arr[idx1]->coord.x + (coord_t)arr[idx2]->coord.x) / 2.0;
    triangle->v3.y = ((coord_t)arr[idx1]->coord.y + (coord_t)arr[idx2]->coord.y) / 2.0;
    triangle->v3.z = ((coord_t)arr[idx1]->coord.z + (coord_t)arr[idx2]->coord.z) / 2.0;
}

__device__ int get_action_value(int less, int eq, int gre){
    
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

    return 100;
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

/// @brief Pre-proccessing phase to remove cubes that will not be part of the computation
/// @param d_grid Pointer to device memorized grid
/// @param cubes_in_domain Total size of the domain
/// @param threshold 
/// @param dim 
/// @param number_relevant_cubes 
/// @param d_relevant_cubes Pointer to device that will be allocated in the function
/// @param d_cube_points_coordinates Pointer to device that will be allocated in the function
/// @param time 
void remove_unnecessary_cubes(  dim_t *d_grid, size_t cubes_in_domain, double threshold,
                                Dimensions *dim, int *number_relevant_cubes,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates,
                                double *time)
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
    *d_dim = *dim;

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

    // print_relevant_points(*d_relevant_cubes, number_relevant_cubes);
}

/// @brief Function that given the cubes in the space can compute the tetrahedra and the tre triangles
/// @param dim 
/// @param d_grid Pointer to the device memorized domain
/// @param cube_decomposition Decomposition of a cube in set of tetrahedra
/// @param threshold 
/// @param number_relevant_cubes 
/// @param d_relevant_cubes 
/// @param d_cube_points_coordinates 
/// @param act_val_vec Vector to act_val fast access
/// @param pairs Vectot to access pairs
/// @param triangles Pointer that will keep the generated triangles on host
/// @param total_triangles 
/// @param time 
void parallel_march_tetra   (Dimensions *dim, dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            int number_relevant_cubes,
                            cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates,
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

/// @brief Skip the pre-processing pass and prepare the DS needed to use the MT algorithm
/// @param cubes_in_domain Total size of the domain (x_dim * y_dim * z_dim)
/// @param dim 
/// @param d_relevant_cubes Device pointer to the DS that keeps all the cubes
/// @param time 
void skip_preprocessing(  size_t cubes_in_domain,
                                Dimensions *dim, cube_gpu **d_relevant_cubes,
                                double *time)
{
    //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

    int n_threads = 512;
    int n_blocks = (cubes_in_domain + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to skip preprocessing\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

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
    
    skip_preprocessing_k<<<n_blocks, n_threads>>>(cubes_in_domain, d_dim,*d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);

    *time += elapsedTime;
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaDeviceSynchronize();

    printf("Total number of cubes:              %zu\n", cubes_in_domain);
    
    // Take the potential error
    print_cuda_error(cudaGetLastError(), "CUDA error");
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");
}

void read_file(const char* file_name, Dimensions *dim, dim_t **tensor, dim_t *origin){

    FILE *fptr;
    fptr = fopen(file_name, "rb"); 

    if(fptr == NULL) {
        fprintf(stderr, "Not able to open the file\n");
        exit(-1);
    }

    dim_t dx;
    dim_t dy;
    dim_t dz;
    
    fread(&(dx), sizeof(dim_t), 1, fptr);
    fread(&(dy), sizeof(dim_t), 1, fptr);
    fread(&(dz), sizeof(dim_t), 1, fptr);
    fread(&(origin[0]), sizeof(dim_t), 1, fptr);
    fread(&(origin[1]), sizeof(dim_t), 1, fptr);
    fread(&(origin[2]), sizeof(dim_t), 1, fptr);
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