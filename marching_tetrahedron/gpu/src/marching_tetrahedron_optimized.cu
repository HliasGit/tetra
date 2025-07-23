#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__constant__ int4 dim;
__constant__ int c_pairs[48];

__global__ void remove_unnecessary_cubes_SoA_kernel(dim_t* grid, int *counter,
                                                size_t size, double threshold,
                                                cube_gpu_SoA* d_relevant_cubes) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;

    int i = idx / (dim.y * dim.z);
    int j = (idx / dim.z) % dim.y;
    int k = idx % dim.z;
    
    extern __shared__ float shared_mem[];
    float* s_mem_1 = shared_mem;
    float* s_mem_2 = s_mem_1 + blockDim.x + 1;
    float* s_mem_3 = s_mem_2 + blockDim.x + 1;
    float* s_mem_4 = s_mem_3 + blockDim.x + 1;

    bool all_in = 0;
    bool all_out = 0;

    if(k != dim.z || j != dim.y || i != dim.x){
        s_mem_1[threadIdx.x] = grid[idx]; //caricato la prima fila
        if(threadIdx.x == blockDim.x-1){ // ultimo
            s_mem_1[blockDim.x] = grid[idx + 1]; //caricato ultimo della prima fila
        }

        s_mem_2[threadIdx.x] = grid[idx + dim.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_2[blockDim.x] = grid[idx + dim.z + 1];
        }

        s_mem_3[threadIdx.x] = grid[idx + dim.z * dim.y];
        if(threadIdx.x == blockDim.x-1){
            s_mem_3[blockDim.x] = grid[idx + dim.z * dim.y + 1];
        }

        s_mem_4[threadIdx.x] = grid[idx + dim.z * dim.y + dim.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_4[blockDim.x] = grid[idx + dim.z * dim.y + dim.z + 1];
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

__global__ void compute_apex_float4(   dim_t *grid, cube_gpu_SoA *d_relevant_cubes, int number_relevant_cubes,
                                cube_vertices_points_SoA *d_cube_points_coordinates){
                        
    float3 one_apex;
    float3 two_apex;
    float3 coord_idx;

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes->coord_idx[tid].x >= dim.x-1 ||
        d_relevant_cubes->coord_idx[tid].y >= dim.y-1 ||
        d_relevant_cubes->coord_idx[tid].z >= dim.z-1 ){
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
                                    coords[v].y * dim.z +
                                    coords[v].x * dim.z * dim.y)];
        
        d_cube_points_coordinates[tid * 8 + v].val = coords[v];
    }
}

__global__ void compute_march_tetra_SoA(dim_t *d_grid, cube_gpu_SoA *d_relevant_cubes,
                                        int number_relevant_cubes, int *cube_deco,
                                        cube_vertices_points_SoA *d_cube_points_coordinates,
                                        cube_vertices_points_SoA *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= number_relevant_cubes) return;

    if( d_relevant_cubes->coord_idx[tid].x >= dim.x-1 ||
        d_relevant_cubes->coord_idx[tid].y >= dim.y-1 ||
        d_relevant_cubes->coord_idx[tid].z >= dim.z-1 ){
            return;
        }

    #pragma unroll
    for (int tetra = 0; tetra < 5; tetra++){

        cube_vertices_points_SoA *first     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+0]-1];
        cube_vertices_points_SoA *second    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+1]-1];
        cube_vertices_points_SoA *third     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+2]-1];
        cube_vertices_points_SoA *fourth    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+3]-1];

        sort_points_SoA(&first, &second, &third, &fourth);
        
        int less = 0, eq = 0, gre = 0;
        
        count_elements_SoA(&less, &eq, &gre, first, second, third, fourth, threshold);
        
        // int act_val = act_val_vec[less + eq*5];
        int act_val = get_action_value(less, eq, gre);
        
        if (act_val!= 0){
            int *use_pairs = &c_pairs[6*(act_val-1)];

            make_triangle_SoA(first, second, third, fourth, &d_triangles[tid*5+tetra], use_pairs);

            if(act_val == 7){
                use_pairs = &c_pairs[42]; // 42 beacuse it's 6*7. 
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

inline __device__ void count_elements_SoA( int *less, int *eq, int *gre, cube_vertices_points_SoA *first,
                                cube_vertices_points_SoA *second, cube_vertices_points_SoA *third, cube_vertices_points_SoA *fourth,
                                dim_t threshold){

    cube_vertices_points_SoA* arr[4] = {first, second, third, fourth};
    #pragma unroll
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

inline __device__ void sort_points_SoA(cube_vertices_points_SoA **first, cube_vertices_points_SoA **second, cube_vertices_points_SoA **third, cube_vertices_points_SoA **fourth){

    cube_vertices_points_SoA* arr[4] = {*first, *second, *third, *fourth};
    #pragma unroll
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
    
    remove_unnecessary_cubes_SoA_kernel<<<n_blocks, n_threads, (n_threads + 1) * 4 * sizeof(float)>>>(   d_grid, d_number_relevant_cubes,
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

void parallel_march_tetra   (dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            int number_relevant_cubes,
                            cube_gpu_SoA *d_relevant_cubes, cube_vertices_points_SoA **d_cube_points_coordinates,
                            int* act_val_vec, int *pairs, Triangle_GPU **triangles, int *total_triangles,
                            double *time)
    {               

    //      //      //      // GENERAL INFO KERNEL APEX //     //      //      //

    int n_threads = 1024;
    int n_blocks = (number_relevant_cubes + n_threads - 1) / n_threads;
    printf("\nLaunching kernel to write apex\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND MANAGE //     //      //      //

    print_cuda_error(cudaMallocManaged(d_cube_points_coordinates, sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8), "cudaMalloc failed for d_cube_points_coordinates");

    printf("Allocating %zu bytes for d_cube_points_coordinates\n", sizeof(cube_vertices_points_SoA)*(number_relevant_cubes)*8);

    printf("Number of cube points coordinates:  %d\n", (number_relevant_cubes) * 8);

    // Take the times for the kernel compute apex
    cudaEvent_t apex_start, apex_stop;
    float apex_elapsedTime;
    cudaEventCreate(&apex_start);
    cudaEventCreate(&apex_stop);
    cudaEventRecord(apex_start, 0);

    compute_apex_float4<<<n_blocks, n_threads>>>(d_grid, d_relevant_cubes, number_relevant_cubes, *d_cube_points_coordinates);
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

    n_threads = 1024;
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
                                                    d_triangles, d_counter);
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

void load_to_const(Dimensions *dimensions, int *pairs){
    int4 tmp;
    tmp.x = dimensions->x_dim;
    tmp.y = dimensions->y_dim;
    tmp.z = dimensions->z_dim;

    tmp.w = 0.0f;
    cudaMemcpyToSymbol(dim, &tmp, sizeof(int4));

    cudaMemcpyToSymbol(c_pairs, pairs, sizeof(int)*48);
}

void print_EC(Triangle_GPU *triangles, int total_triangles){ //TODO: FINIRE
    // process triangle vertices in a C++ list
    std::unordered_set<TriangleVertex_GPU> vertices_list;

    // Timing vertex processing using standard C++ chrono (serial code)
    auto vertex_start = std::chrono::high_resolution_clock::now();

    for(int t=0; t<total_triangles; t++){
        if (t % 100000 == 0) {
            printf("\rProcessed %d / %d triangles", t, total_triangles);
            fflush(stdout);
        }
        vertices_list.insert(triangles[t].v1);
        vertices_list.insert(triangles[t].v2);
        vertices_list.insert(triangles[t].v3);
    }

    auto vertex_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> vertex_time = vertex_end - vertex_start;
    printf("\nVertex processing time: %f ms\n", vertex_time.count());
}