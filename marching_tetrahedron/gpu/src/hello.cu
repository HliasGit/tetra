#include <hello.h>
#include <hello.cuh>

__global__ void hello_kernel()
{
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    printf("Hello from thread #%ld!\n", idx);
}

__global__ void remove_unnecessary_cubes_kernel(double* grid, int *counter,
                                                size_t size, double threshold,
                                                Dimensions *dim, int* d_list) {

    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if ((idx + dim->z_dim * dim->y_dim + dim->z_dim+1) >= size) return;

    bool all_out = //Wrong when on the border of the domain
        (grid[idx] < threshold) &&
        (grid[idx+1] < threshold) &&
        (grid[idx + dim->z_dim] < threshold) &&
        (grid[idx + dim->z_dim+1] < threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim] < threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + 1] < threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim] < threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim+1] < threshold);

    bool all_in = //Wrong when on the border of the domain
        (grid[idx] > threshold) &&
        (grid[idx+1] > threshold) &&
        (grid[idx + dim->z_dim] > threshold) &&
        (grid[idx + dim->z_dim+1] > threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim] > threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + 1] > threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim] > threshold) &&
        (grid[idx + dim->z_dim * dim->y_dim + dim->z_dim+1] > threshold);

    if (all_out == 0 && all_in == 0){
        // To store indices dynamically, you need to allocate a list with a size at least as large as the maximum possible number of insertions.
        // In this case, the maximum is 'size', since at most each thread could insert once.
        int insert_pos = atomicAdd(counter, 1);
        d_list[insert_pos] = idx;
    }

}

void hello_f()
{
    hello_kernel<<<1, 256>>>();
    cudaDeviceSynchronize();
}

int remove_unnecessary_cubes(  double *grid, size_t size, double threshold,
                                Dimensions *dim, int **results)
{
    int n_threads = 512;
    int n_blocks = (size + n_threads - 1) / n_threads;
    int counter = 0;
    int *list = (int*)malloc(sizeof(int)*size);

    double *d_grid;
    cudaMalloc(&d_grid, size*sizeof(double));
    cudaMemcpy(d_grid, grid, size*sizeof(double), cudaMemcpyHostToDevice);

    int *d_counter;
    cudaMallocManaged(&d_counter, sizeof(int));
    *d_counter = counter;

    Dimensions *d_dim;
    cudaMallocManaged(&d_dim, sizeof(Dimensions));
    *d_dim = *dim;

    int *d_list;
    cudaMallocManaged(&d_list, sizeof(int)*size);
    *d_list = *list;

    cudaEvent_t start, stop;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    remove_unnecessary_cubes_kernel<<<n_blocks, n_threads>>>(   d_grid, d_counter,
                                                                size, threshold, d_dim,
                                                                d_list);
    cudaDeviceSynchronize();

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);

    cudaEventElapsedTime(&elapsedTime, start, stop);

    printf("Kernel execution time: %f ms\n", elapsedTime);

    printf("Number of elements <= threshold: %d\n", (*d_counter));
    printf("Total points: %zu\n", size);

    // for (int i = 0; i < *d_counter; ++i) {
    //     printf("list[%d] = %d\n", i, d_list[i]);
    // }

    // FILE *fp = fopen("out.pdb", "w");
    // if (fp == NULL) {
    //     fprintf(stderr, "Error opening out.pdb for writing.\n");
    // } else {
    //     for (int i = 0; i < *d_counter; ++i) {
    //         int idx = d_list[i];
    //         int DZ = dim->z_dim;
    //         int DY = dim->y_dim;
    //         int k = idx % DZ;
    //         int j = (idx / DZ) % DY;
    //         int ii = idx / (DZ * DY);
    //         // PDB ATOM line: ATOM  serial  name resName chainID resSeq x y z
    //         // We'll use serial = i+1, name = "C", resName = "RES", chainID = 'A', resSeq = 1, x/y/z = ii/j/k
    //         fprintf(fp, "ATOM  %5d  C   RES A   1    %8.3f%8.3f%8.3f\n", i+1, (float)ii, (float)j, (float)k);
    //     }
    //     fclose(fp);
    // }

    // Assign the device pointer to results for use outside this function
    *results = d_list;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
    }
    printf("remove_unnecessary_cubes_kernel executed successfully.\n");

    return (*d_counter);
}

