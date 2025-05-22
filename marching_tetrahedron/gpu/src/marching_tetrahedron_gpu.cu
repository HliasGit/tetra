#include <marching_tetrahedron_gpu.h>
#include <marching_tetrahedron_gpu.cuh>

__global__ void remove_unnecessary_cubes_kernel(double* grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, cube_gpu* d_relevant_cubes) {

    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    
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
                                        dim_t threshold){

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

        int act_val = get_action_value(less, eq, gre);

        if (act_val!= 0)
            printf("act_val: %d\n", act_val);
        
    }
}

__device__ int get_action_value( int less, int eq, int gre){
    // if(less == 3){
    //     printf("AAAAAAAAAAAAA\n");
    // }

    if (less == 0 || (less == 2 && eq == 2) || (less == 3 && eq == 1) || less == 4)
    {
        return 0;
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

    if (less == 2)
    {
        if (eq == 0)
            return 7;
        if (eq == 1)
            return 5;
    }

    if (less == 3)
        return 6;
}

__device__ void count_elements( int *less, int *eq, int *gre, cube_vertices_points *first,
                                cube_vertices_points *second, cube_vertices_points *third, cube_vertices_points *fourth,
                                dim_t threshold){

    int less_ = 0, eq_ = 0, gre_ = 0;

    // if (first != nullptr && first->value < threshold)
    //     printf("first value: %f\n", first->value);
    // if (second != nullptr && second->value < threshold)
    //     printf("second value: %f\n", second->value);
    // if (third != nullptr && third->value < threshold)
    //     printf("third value: %f\n", third->value);
    // if (fourth != nullptr && fourth->value < threshold)
    //     printf("fourth value: %f\n", fourth->value);

    cube_vertices_points* arr[4] = {first, second, third, fourth};
    for (int i = 0; i < 4; ++i) {
        if (arr[i]->value < threshold) {
            (*less)++;
            less_++;
        } else if (arr[i]->value < threshold) {
            (*eq)++;
            eq_++;
        } else {
            (*gre)++;
            gre_++;
        }
    }

    // Print only if less, eq, and gre are all not 0 or 4
    if ((*less != 0 && *less != 4) && (*eq != 0 && *eq != 4) && (*gre != 0 && *gre != 4)) {
        printf("less: %d, eq: %d, gre: %d\n", *less, *eq, *gre);
    }
}

__device__ void sort_points(cube_vertices_points **first, cube_vertices_points **second, cube_vertices_points **third, cube_vertices_points **fourth){
    // printf("AAAA: %f\n", (*first)->value);

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
    int n_threads = 512;
    int n_blocks = (total_size + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to remove unnecessary cuvbes\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    int *d_relevant;
    cudaError_t alloc_err;
    alloc_err = cudaMallocManaged(&d_relevant, sizeof(int));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_relevant: %s\n", cudaGetErrorString(alloc_err));
        return;
    }
    *d_relevant = 0;

    Dimensions *d_dim;
    alloc_err = cudaMallocManaged(&d_dim, sizeof(Dimensions));
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_dim: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_relevant);
        return;
    }
    *d_dim = *dim;

    alloc_err = cudaMallocManaged(d_relevant_cubes, sizeof(cube_gpu)*total_size);
    if (alloc_err != cudaSuccess) {
        fprintf(stderr, "cudaMallocManaged failed for d_relevant_cubes: %s\n", cudaGetErrorString(alloc_err));
        cudaFree(d_relevant);
        cudaFree(d_dim);
        return;
    }

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

    printf("Kernel rm cubes execution time:     %f ms\n", elapsedTime);

    printf("Number of elements <= threshold:    %d\n", (*d_relevant));
    printf("Total number of cubes:              %zu\n", total_size);

    *relevant_cubes_size = *d_relevant;

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
    }
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    n_blocks = (*d_relevant + n_threads - 1) / n_threads;
    printf("\nLaunching kernel to write apex\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);

    cudaMallocManaged(d_cube_points_coordinates, sizeof(cube_vertices_points)*(*d_relevant)*8);

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

    printf("compute_apex kernel execution time: %f ms\n", apex_elapsedTime);

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
                            cube_gpu **d_relevant_cubes, cube_vertices_points **d_cube_points_coordinates)
{
    (*vertex_counter) = 0;
    (*triangle_counter) = 0;


    int n_threads = 512;
    int n_blocks = (relevant_size + n_threads - 1) / n_threads;

    printf("\nLaunching kernel to compute MT algo\n");
    printf("# blocks                            %d \nand # threads                       %d \n", n_blocks, n_threads);


    cudaMallocManaged(&cube_decomposition, sizeof(int)*20);

    cube_vertices_points *d_stack_pool;
    cudaMallocManaged(&d_stack_pool, sizeof(cube_vertices_points)*20*relevant_size);
    
    cudaEvent_t mt_start, mt_stop;
    float mt_elapsedTime;

    cudaEventCreate(&mt_start);
    cudaEventCreate(&mt_stop);

    cudaEventRecord(mt_start, 0);

    int *pool_index = 0;
    cudaMallocManaged(&pool_index, sizeof(int));

    compute_march_tetra<<<n_blocks, n_threads>>>(   d_grid, *d_relevant_cubes,
                                                    relevant_size, cube_decomposition,
                                                    *d_cube_points_coordinates, d_stack_pool, pool_index,
                                                    threshold);
    cudaDeviceSynchronize();

    cudaError_t mt_err = cudaGetLastError();
    if (mt_err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error in compute_march_tetra: %s\n", cudaGetErrorString(mt_err));
    }

    cudaEventRecord(mt_stop, 0);
    cudaEventSynchronize(mt_stop);

    cudaEventElapsedTime(&mt_elapsedTime, mt_start, mt_stop);

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