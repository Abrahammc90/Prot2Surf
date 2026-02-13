#include <stdio.h>
#include <float.h>
#include <cuda_runtime.h>

/*
 * CUDA kernel for finding minimum distance in distance matrix
 * Uses proper parallel reduction with warp-level and block-level synchronization
 */

// Helper for atomic min on doubles using CAS
__device__ double atomicMin(double* address, double val) {
    unsigned long long *address_as_ull = (unsigned long long *) address;
    unsigned long long old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void find_min_pair_kernel(
    const double *matrix,
    const int *active_points,
    int n,
    double *temp_dist,
    int *temp_i,
    int *temp_j
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (j >= n) return;
    if (!active_points[j]) return;
    
    // Each thread processes one column
    double thread_min = DBL_MAX;
    int thread_i = -1;
    int thread_j = -1;
    
    // Search for minimum in this column
    for (int i = 0; i < j; i++) {
        if (active_points[i]) {
            // Column-major: matrix(i,j) is at matrix[j*n + i]
            double dist = matrix[j * n + i];
            if (dist < thread_min) {
                thread_min = dist;
                thread_i = i;
                thread_j = j;
            }
        }
    }
    
    // Store this thread's result in global memory
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_i >= 0 && idx < n) {
        temp_dist[idx] = thread_min;
        temp_i[idx] = thread_i;
        temp_j[idx] = thread_j;
    }
}

__global__ void reduce_min_kernel(
    double *temp_dist,
    int *temp_i,
    int *temp_j,
    int n,
    int *result_i,
    int *result_j,
    double *result_dist
) {
    // Find global minimum from temporary results
    double global_min = DBL_MAX;
    int global_i = -1;
    int global_j = -1;
    
    for (int idx = 0; idx < n; idx++) {
        if (temp_dist[idx] < global_min) {
            global_min = temp_dist[idx];
            global_i = temp_i[idx];
            global_j = temp_j[idx];
        }
    }
    
    result_dist[0] = global_min;
    result_i[0] = global_i;
    result_j[0] = global_j;
}

extern "C" {

void cuda_find_min_pair_c(
    const double *h_matrix,
    const int *h_active_points,
    int n,
    int *h_result_i,
    int *h_result_j,
    double *h_result_dist
) {
    double *d_matrix = nullptr;
    int *d_active_points = nullptr;
    int *d_result_i = nullptr;
    int *d_result_j = nullptr;
    double *d_result_dist = nullptr;
    double *d_temp_dist = nullptr;
    int *d_temp_i = nullptr;
    int *d_temp_j = nullptr;
    
    // Allocate device memory
    cudaMalloc((void**)&d_matrix, n * n * sizeof(double));
    cudaMalloc((void**)&d_active_points, n * sizeof(int));
    cudaMalloc((void**)&d_result_i, sizeof(int));
    cudaMalloc((void**)&d_result_j, sizeof(int));
    cudaMalloc((void**)&d_result_dist, sizeof(double));
    cudaMalloc((void**)&d_temp_dist, n * sizeof(double));
    cudaMalloc((void**)&d_temp_i, n * sizeof(int));
    cudaMalloc((void**)&d_temp_j, n * sizeof(int));
    
    // Copy input data to device
    cudaMemcpy(d_matrix, h_matrix, n * n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_active_points, h_active_points, n * sizeof(int), cudaMemcpyHostToDevice);
    
    // Initialize temporary and result arrays
    double init_dist = 1e308;
    int init_i = -1, init_j = -1;
    cudaMemset(d_temp_dist, 0xFF, n * sizeof(double));  // Fill with max double
    cudaMemset(d_temp_i, 0xFF, n * sizeof(int));
    cudaMemset(d_temp_j, 0xFF, n * sizeof(int));
    cudaMemcpy(d_result_i, &init_i, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_result_j, &init_j, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_result_dist, &init_dist, sizeof(double), cudaMemcpyHostToDevice);
    
    // Configure kernel
    int threads_per_block = 256;
    int blocks = (n + threads_per_block - 1) / threads_per_block;
    
    // Launch first kernel to find minimums in each column
    find_min_pair_kernel<<<blocks, threads_per_block>>>(
        d_matrix, d_active_points, n,
        d_temp_dist, d_temp_i, d_temp_j
    );
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA first kernel error: %s\n", cudaGetErrorString(err));
        goto cleanup;
    }
    
    cudaDeviceSynchronize();
    
    // Launch reduction kernel to find global minimum
    reduce_min_kernel<<<1, 1>>>(
        d_temp_dist, d_temp_i, d_temp_j, n,
        d_result_i, d_result_j, d_result_dist
    );
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA reduction kernel error: %s\n", cudaGetErrorString(err));
        goto cleanup;
    }
    
    cudaDeviceSynchronize();
    
    // Copy results back
    cudaMemcpy(h_result_i, d_result_i, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_result_j, d_result_j, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_result_dist, d_result_dist, sizeof(double), cudaMemcpyDeviceToHost);

cleanup:
    // Cleanup
    if (d_matrix) cudaFree(d_matrix);
    if (d_active_points) cudaFree(d_active_points);
    if (d_result_i) cudaFree(d_result_i);
    if (d_result_j) cudaFree(d_result_j);
    if (d_result_dist) cudaFree(d_result_dist);
    if (d_temp_dist) cudaFree(d_temp_dist);
    if (d_temp_i) cudaFree(d_temp_i);
    if (d_temp_j) cudaFree(d_temp_j);
}

}
