#include "mt_noatom.h"
#include "mt_noatom.cuh"
#include "marching_tetrahedron_gpu.h"

__constant__ int3 dim_flags;

__global__ void set_flags(  dim_t* d_grid,
                            bool *d_flags,
                            size_t size,
                            double threshold,
                            cube_gpu_SoA* d_relevant_cubes
                        )
    {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if(idx >= size) return;

    int i = idx / (dim_flags.y * dim_flags.z);
    int j = (idx / dim_flags.z) % dim_flags.y;
    int k = idx % dim_flags.z;
    
    extern __shared__ float shared_mem[];
    float* s_mem_1 = shared_mem;
    float* s_mem_2 = s_mem_1 + blockDim.x + 1;
    float* s_mem_3 = s_mem_2 + blockDim.x + 1;
    float* s_mem_4 = s_mem_3 + blockDim.x + 1;

    bool all_in = 0;
    bool all_out = 0;

    if(k != dim_flags.z-1 || j != dim_flags.y-1 || i != dim_flags.x-1){
        s_mem_1[threadIdx.x] = d_grid[idx]; //caricato la prima fila
        if(threadIdx.x == blockDim.x-1){ // ultimo
            s_mem_1[blockDim.x] = d_grid[idx + 1]; //caricato ultimo della prima fila
        }

        s_mem_2[threadIdx.x] = d_grid[idx + dim_flags.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_2[blockDim.x] = d_grid[idx + dim_flags.z + 1];
        }

        s_mem_3[threadIdx.x] = d_grid[idx + dim_flags.z * dim_flags.y];
        if(threadIdx.x == blockDim.x-1){
            s_mem_3[blockDim.x] = d_grid[idx + dim_flags.z * dim_flags.y + 1];
        }

        s_mem_4[threadIdx.x] = d_grid[idx + dim_flags.z * dim_flags.y + dim_flags.z];
        if(threadIdx.x == blockDim.x-1){
            s_mem_4[blockDim.x] = d_grid[idx + dim_flags.z * dim_flags.y + dim_flags.z + 1];
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
        d_flags[idx] = true;
        d_relevant_cubes->coord_idx[idx].w = idx;
        d_relevant_cubes->coord_idx[idx].x = i;
        d_relevant_cubes->coord_idx[idx].y = j;
        d_relevant_cubes->coord_idx[idx].z = k;
    }
}

__global__ void compact_scan(int *d_idx, int* d_scan, bool *d_flags, int size){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx >= size) return;

    if (d_flags[idx]) {
        int out_idx = d_scan[idx];
        d_idx[out_idx] = idx;
    }
}

void remove_unnecessary_cubes_flag( dim_t* d_grid,
                                    bool *d_flags,
                                    size_t cubes_in_domain,
                                    double threshold,
                                    cube_gpu_SoA** d_relevant_cubes,
                                    float *time
                                    )
{
    int n_threads = 1024;
    int n_blocks = (cubes_in_domain + n_threads - 1) / n_threads;

    printf("Launching set_flags kernel with %d blocks and %d threads per block\n", n_blocks, n_threads);

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


    // Setup time reader for the kernel
    cudaEvent_t start, stop;
    float elapsedTime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    set_flags<<<n_blocks, n_threads, (n_threads + 1) * 4 * sizeof(float)>>>(    d_grid,
                                                                                d_flags,
                                                                                cubes_in_domain,
                                                                                threshold,
                                                                                *d_relevant_cubes);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    cudaDeviceSynchronize();
    
    printf("Kernel remove cubes execution time: %f ms\n", elapsedTime);
    *time = elapsedTime;
    

    // EXCLUSIVE SCAN WITH THRUST
    // Timing device_vector allocation
    cudaEvent_t alloc_start, alloc_stop;
    float allocTime = 0;
    cudaEventCreate(&alloc_start);
    cudaEventCreate(&alloc_stop);
    cudaEventRecord(alloc_start, 0);

    thrust::device_vector<int> thrust_flags(cubes_in_domain);
    thrust::device_vector<int> thrust_scan(cubes_in_domain);

    cudaEventRecord(alloc_stop, 0);
    cudaEventSynchronize(alloc_stop);
    cudaEventElapsedTime(&allocTime, alloc_start, alloc_stop);

    cudaEventDestroy(alloc_start);
    cudaEventDestroy(alloc_stop);

    printf("Thrust device_vector allocation time: %f ms\n", allocTime);

    // Timing exclusive scan
    cudaEvent_t scan_start, scan_stop;
    float scanTime = 0;
    cudaEventCreate(&scan_start);
    cudaEventCreate(&scan_stop);
    cudaEventRecord(scan_start, 0);

    thrust::transform(
        thrust::device_pointer_cast(d_flags),
        thrust::device_pointer_cast(d_flags + cubes_in_domain),
        thrust_flags.begin(),
        thrust::identity<int>()  // or a custom lambda if needed
    );

    // Now perform scan on thrust_flags
    thrust::exclusive_scan(
        thrust_flags.begin(),
        thrust_flags.end(),
        thrust_scan.begin()
    );

    cudaEventRecord(scan_stop, 0);
    cudaEventSynchronize(scan_stop);
    cudaEventElapsedTime(&scanTime, scan_start, scan_stop);

    cudaEventDestroy(scan_start);
    cudaEventDestroy(scan_stop);

    printf("Thrust exclusive scan execution time: %f ms\n", scanTime);
    *time += allocTime + scanTime;

    // KENREL TO COMPACT THE SCAN
    int *d_idx;
    print_cuda_error(cudaMallocManaged(&d_idx, sizeof(int) * cubes_in_domain), "cudaMallocManaged failed for d_idx:");

    int* d_scan;
    print_cuda_error(cudaMallocManaged(&d_scan, sizeof(int) * cubes_in_domain), "cudaMallocManaged failed for d_scan:");
    cudaMemcpy(d_scan, thrust::raw_pointer_cast(thrust_scan.data()), sizeof(int) * cubes_in_domain, cudaMemcpyDeviceToDevice);
    
    // Timing compact_scan kernel
    cudaEvent_t compact_start, compact_stop;
    float compactTime = 0;
    cudaEventCreate(&compact_start);
    cudaEventCreate(&compact_stop);
    cudaEventRecord(compact_start, 0);

    compact_scan<<<n_blocks, n_threads>>>(d_idx, d_scan, d_flags, cubes_in_domain);
    cudaDeviceSynchronize();

    cudaEventRecord(compact_stop, 0);
    cudaEventSynchronize(compact_stop);
    cudaEventElapsedTime(&compactTime, compact_start, compact_stop);

    cudaEventDestroy(compact_start);
    cudaEventDestroy(compact_stop);

    printf("compact_scan kernel execution time: %f ms\n", compactTime);
    *time += compactTime;

    // // Print d_idx values for debugging
    // printf("Printing d_idx values:\n");
    // for (size_t i = 0; i < cubes_in_domain; ++i) {
    //     if(d_idx[i] != 0)
    //     printf("d_idx[%zu] = %d\n", i, d_idx[i]);
    // }

}

void load_to_const_tetra_flags(Dimensions *dimensions){

    printf("Dimensions: x=%d, y=%d, z=%d\n", dimensions->x_dim, dimensions->y_dim, dimensions->z_dim);
    int3 tmp;
    tmp.x = dimensions->x_dim;
    tmp.y = dimensions->y_dim;
    tmp.z = dimensions->z_dim;

    cudaMemcpyToSymbol(dim_flags, &tmp, sizeof(int3));

    printf("dim_flags copied to const mem\n");
}

void allocate_d_flags(bool **d_flags, size_t size){
    cudaMalloc(d_flags, size * sizeof(bool));
    cudaMemset(*d_flags, 0, size * sizeof(bool));
}