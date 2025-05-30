#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__global__ void remove_unnecessary_cubes_kernel(double* grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu* d_relevant_cubes) {

    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;
    
    int DZ = dim->z_dim;
    int DY = dim->y_dim;
    int DX = dim->x_dim;
    int k = idx % DZ;
    int j = (idx / DZ) % DY;
    int i = idx / (DZ * DY);

    bool all_in;
    bool all_out;

    if(k == DZ-1 || j == DY-1 || i == DX-1){
        all_in = 0;
        all_out = 0;
    } else {
        all_out = //Wrong when on the border of the domain
            (grid[idx] < threshold) &&
            (grid[idx+1] < threshold) &&
            (grid[idx + dim->z_dim] < threshold) &&
            (grid[idx + dim->z_dim+1] < threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim] < threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + 1] < threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim] < threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim+1] < threshold);

        all_in = //Wrong when on the border of the domain
            (grid[idx] > threshold) &&
            (grid[idx+1] > threshold) &&
            (grid[idx + dim->z_dim] > threshold) &&
            (grid[idx + dim->z_dim+1] > threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim] > threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + 1] > threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim] > threshold) &&
            (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim+1] > threshold);
    }

    if (all_out == 0 && all_in == 0){
        // To store indices dynamically, you need to allocate a list with a size at least as large as the maximum possible number of insertions.
        // In this case, the maximum is 'size', since at most each thread could insert once.
        int insert_pos = atomicAdd(counter, 1);
        d_relevant_cubes[insert_pos].idx = idx;
        d_relevant_cubes[insert_pos].x = i;
        d_relevant_cubes[insert_pos].y = j;
        d_relevant_cubes[insert_pos].z = k;
    }

}

__global__ void compute_apex(   double *grid, cube_gpu *d_relevant_cubes, int *relevant_size,
                                cube_vertices_points *d_cube_points_coordinates, Dimensions *dim){
    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= *relevant_size) return;

    
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
    
    for (int d = 0; d < 3; ++d) {
        d_relevant_cubes[tid].two_apex[d] = 2 * (&d_relevant_cubes[tid].x)[d] + 1 - d_relevant_cubes[tid].one_apex[d];
    }
    
    // Fill the 8 cube vertex coordinates using a for loop
    for (int v = 0; v < 8; ++v) {
        d_cube_points_coordinates[tid * 8 + v].coord.x = (v & 1) ? d_relevant_cubes[tid].two_apex[0] : d_relevant_cubes[tid].one_apex[0];
        d_cube_points_coordinates[tid * 8 + v].coord.y = (v & 2) ? d_relevant_cubes[tid].two_apex[1] : d_relevant_cubes[tid].one_apex[1];
        d_cube_points_coordinates[tid * 8 + v].coord.z = (v & 4) ? d_relevant_cubes[tid].two_apex[2] : d_relevant_cubes[tid].one_apex[2];
        d_cube_points_coordinates[tid * 8 + v].value = grid[  d_cube_points_coordinates[tid * 8 + v].coord.z +
                                                                    d_cube_points_coordinates[tid * 8 + v].coord.y * dim->z_dim +
                                                                    d_cube_points_coordinates[tid * 8 + v].coord.z * dim->z_dim * dim->y_dim];
    }


}

__global__ void compute_march_tetra(    double *d_grid, cube_gpu *d_relevant_cubes,
                                        int relevant_size, int *cube_deco,
                                        cube_vertices_points *d_cube_points_coordinates,
                                        cube_vertices_points *memory_pool, int *pool_index,
                                        dim_t threshold, int* act_val_vec, int *d_pairs,
                                        Triangle_GPU *d_triangles, int *d_counter){

    size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

    if(tid >= relevant_size) return;

    // printf("tid: %d\n", tid);

    for (int tetra = 0; tetra < 5; tetra++){

        cube_vertices_points *first     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+0]];
        cube_vertices_points *second    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+1]];
        cube_vertices_points *third     = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+2]];
        cube_vertices_points *fourth    = &d_cube_points_coordinates[tid*8+cube_deco[tetra*4+3]];

        sort_points(&first, &second, &third, &fourth);

        int less = 0, eq = 0, gre = 0;

        count_elements(&less, &eq, &gre, first, second, third, fourth, threshold);

        // int act_val2 = get_action_value(less, eq, gre);
        int act_val = act_val_vec[less + eq*5];

        if (act_val!= 0){
            
            int *use_pairs = &d_pairs[6*(act_val-1)]; 

            make_triangle(first, second, third, fourth, &d_triangles[tid*5+tetra], use_pairs);

            if(act_val == 7){
                use_pairs = &d_pairs[42]; // 42 beacuse it's 6*7. 

                int insert_pos = atomicAdd(d_counter, 1);

                make_triangle(first, second, third, fourth, &d_triangles[relevant_size*5 + insert_pos], use_pairs);
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
    // Assign sorted pointers back
    *first = arr[0];
    *second = arr[1];
    *third = arr[2];
    *fourth = arr[3];
}



void remove_unnecessary_cubes(  double *d_grid, size_t total_size, double threshold,
                                Dimensions *dim, int *relevant_cubes_size,
                                cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates)
{
    //      //      //      // GENERAL INFO KERNEL 1 //      //      //      //

    int n_threads = 512;
    int n_blocks = (total_size + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to remove unnecessary cubes\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND COPY //       //      //      //

    cudaError_t alloc_err;

    // Number of relevant cubes (the ones that are not all in or all out)
    int *d_relevant;
    alloc_err = cudaMallocManaged(&d_relevant, sizeof(int));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_relevant: %s\n", cudaGetErrorString(alloc_err));
        return;
    }
    *d_relevant = 0;    // Initialize to 0

    // Data structure for the dimensions
    Dimensions *d_dim;
    alloc_err = cudaMallocManaged(&d_dim, sizeof(Dimensions));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_dim: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_relevant);
        return;
    }
    *d_dim = *dim;      // Copy the DS from the serial code

    // Big array that keeps every point in the grid with some metadata
    alloc_err = cudaMallocManaged(d_relevant_cubes, sizeof(cube_gpu)*total_size);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_relevant_cubes: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_relevant);
        cudaFree(d_dim);
        return;
    }
    


    // Setup time reader for the kernel
    cudaEvent_t start, stop;
    float elapsedTime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    remove_unnecessary_cubes_kernel<<<n_blocks, n_threads>>>(   d_grid, d_relevant,
                                                                total_size, threshold, d_dim,
                                                                *d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);
    
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    // Copy back the # of relevant cubes
    *relevant_cubes_size = *d_relevant;

    printf("Number of unnecessary cubes:        %d\n", (*relevant_cubes_size));
    printf("Total number of cubes:              %zu\n", total_size);
    
    // Take the potential error
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
    }
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");
    
    //      //      //      // GENERAL INFO KERNEL 2 //     //      //      //

    n_blocks = (*d_relevant + n_threads - 1) / n_threads;
    printf("\nLaunching kernel to write apex\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    //      //      //      // MALLOC AND MANAGE //     //      //      //

    cudaMallocManaged(d_cube_points_coordinates, sizeof(cube_vertices_points)*(*d_relevant)*8);

    // Take the times for the kernel compute apex
    cudaEvent_t apex_start, apex_stop;
    float apex_elapsedTime;
    cudaEventCreate(&apex_start);
    cudaEventCreate(&apex_stop);
    cudaEventRecord(apex_start, 0);

    compute_apex<<<n_blocks, n_threads>>>(d_grid, *d_relevant_cubes, d_relevant, *d_cube_points_coordinates, d_dim);
    cudaDeviceSynchronize();

    cudaEventRecord(apex_stop, 0);
    cudaEventSynchronize(apex_stop);
    cudaEventElapsedTime(&apex_elapsedTime, apex_start, apex_stop);

    printf("Compute_apex kernel execution time: %f ms\n", apex_elapsedTime);

    cudaError_t apex_err = cudaGetLastError();
    if (apex_err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error in compute_apex: %s\n", cudaGetErrorString(apex_err));
    }

    cudaEventDestroy(apex_start);
    cudaEventDestroy(apex_stop);
}

void parallel_march_tetra   (Dimensions *dim, dim_t *d_grid, int *cube_decomposition, dim_t threshold,
                            void (*func_ptr)(TriangleVertex *, CubeVertex *, CubeVertex *, dim_t *, dim_t *, dim_t),
                            Polyhedra *p, size_t *triangle_counter, size_t *vertex_counter, int relevant_size,
                            cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates,
                            int* act_val_vec, int *pairs, Triangle_GPU **triangles)
{               
    //      //      //      // GENERAL INFO KERNEL MT //        //      //      //

    int n_threads = 512;
    int n_blocks = (relevant_size + n_threads - 1) / n_threads;
    
    printf("\nLaunching kernel to compute MT algo\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);
    
    //      //      //      // MALLOC AND MANAGE //     //      //      //

    cudaError_t alloc_err;

    (*vertex_counter) = 0;
    (*triangle_counter) = 0;

    int *d_cube_decomposition;
    alloc_err = cudaMalloc(&d_cube_decomposition, sizeof(int)*20);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed for d_cube_decomposition: %s\n", cudaGetErrorString(alloc_err));
        return;
    }
    alloc_err = cudaMemcpy(d_cube_decomposition, cube_decomposition, sizeof(int) * 20, cudaMemcpyHostToDevice);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed for d_cube_decomposition: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_cube_decomposition);
        return;
    }

    int *d_act_val_vec;
    alloc_err = cudaMalloc(&d_act_val_vec, sizeof(int)*25);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed for d_act_val_vec alloc: %s\n", cudaGetErrorString(alloc_err));
        return;
    }
    alloc_err = cudaMemcpy(d_act_val_vec, act_val_vec, sizeof(int) * 25, cudaMemcpyHostToDevice);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed for d_act_val_vec copy: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        return;
    }

    int *d_pairs;
    alloc_err = cudaMalloc(&d_pairs, sizeof(int)*48);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed for d_pairs: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        return;
    }
    alloc_err = cudaMemcpy(d_pairs, pairs, sizeof(int) * 48, cudaMemcpyHostToDevice);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed for d_pairs: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        cudaFree(d_pairs);
        return;
    }

    Triangle_GPU *d_triangles;
    alloc_err = cudaMalloc(&d_triangles, sizeof(Triangle_GPU)*12*relevant_size);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed for d_triangles: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        cudaFree(d_pairs);
        return;
    }

    int *d_counter;
    alloc_err = cudaMallocManaged(&d_counter, sizeof(int));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_counter: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        cudaFree(d_pairs);
        cudaFree(d_triangles);
        cudaFree(d_cube_decomposition);
        return;
    }
    *d_counter = 0;

    cube_vertices_points *d_stack_pool;
    alloc_err = cudaMallocManaged(&d_stack_pool, sizeof(cube_vertices_points)*20*relevant_size);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_stack_pool: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        cudaFree(d_pairs);
        cudaFree(d_triangles);
        cudaFree(d_cube_decomposition);
        cudaFree(d_counter);
        return;
    }
    
    int *pool_index;
    alloc_err = cudaMallocManaged(&pool_index, sizeof(int));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for pool_index: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_act_val_vec);
        cudaFree(d_pairs);
        cudaFree(d_triangles);
        cudaFree(d_cube_decomposition);
        cudaFree(d_counter);
        cudaFree(d_stack_pool);
        return;
    }
    
    // Take the times
    cudaEvent_t mt_start, mt_stop;
    float mt_elapsedTime;
    cudaEventCreate(&mt_start);
    cudaEventCreate(&mt_stop);
    cudaEventRecord(mt_start, 0);

    compute_march_tetra<<<n_blocks, n_threads>>>(   d_grid, *d_relevant_cubes,
                                                    relevant_size, d_cube_decomposition,
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
    int total_triangles = relevant_size * 5 + *d_counter; // 5 tetrahedra per cube plus possible extra triangles
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