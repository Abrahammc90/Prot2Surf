#include <stdio.h>
#include <float.h>
#include <cuda_runtime.h>

/*
 * CUDA kernels for GPU-resident hierarchical clustering
 * Matrix stays on GPU throughout all iterations to avoid transfer overhead
 */

// GPU state structure for persistent memory
struct GPUClusterState {
    double *d_matrix;
    int *d_active_points;
    int *d_cluster_size;
    double *d_temp_dist;
    int *d_temp_i;
    int *d_temp_j;
    int *d_result_i;
    int *d_result_j;
    double *d_result_dist;
    int n;
};

static GPUClusterState* gpu_state = nullptr;

// Kernel to find minimum in each column
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
    if (!active_points[j]) {
        temp_dist[j] = DBL_MAX;
        temp_i[j] = -1;
        temp_j[j] = -1;
        return;
    }
    
    // Each thread processes one column
    double thread_min = DBL_MAX;
    int thread_i = -1;
    
    // Search for minimum in this column (lower triangle only)
    for (int i = 0; i < j; i++) {
        if (active_points[i]) {
            // Column-major: matrix(i,j) is at matrix[j*n + i]
            double dist = matrix[j * n + i];
            if (dist < thread_min) {
                thread_min = dist;
                thread_i = i;
            }
        }
    }
    
    // Store this thread's result
    temp_dist[j] = thread_min;
    temp_i[j] = thread_i;
    temp_j[j] = j;
}

// Kernel to reduce to global minimum
__global__ void reduce_min_kernel(
    const double *temp_dist,
    const int *temp_i,
    const int *temp_j,
    int n,
    int *result_i,
    int *result_j,
    double *result_dist
) {
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
    
    *result_dist = global_min;
    *result_i = global_i;
    *result_j = global_j;
}

// Kernel to update distances after merging clusters i and j
// linkage_type: 0=min, 1=max, 2=mean
__global__ void update_distances_kernel(
    double *matrix,
    int n,
    int merge_i,
    int merge_j,
    int linkage_type,
    const int *cluster_size
) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (k >= n || k == merge_i || k == merge_j) return;
    
    // Get distances from k to clusters i and j
    double dist_ki, dist_kj;
    
    // Access matrix element (k, merge_i)
    if (k < merge_i) {
        dist_ki = matrix[merge_i * n + k];  // k is row, merge_i is column
    } else {
        dist_ki = matrix[k * n + merge_i];  // merge_i is row, k is column
    }
    
    // Access matrix element (k, merge_j)
    if (k < merge_j) {
        dist_kj = matrix[merge_j * n + k];
    } else {
        dist_kj = matrix[k * n + merge_j];
    }
    
    // Compute new distance based on linkage
    double new_dist;
    if (linkage_type == 0) {  // min linkage
        new_dist = fmin(dist_ki, dist_kj);
    } else if (linkage_type == 1) {  // max linkage
        new_dist = fmax(dist_ki, dist_kj);
    } else {  // mean linkage - weighted average by cluster size
        new_dist = (dist_ki * cluster_size[merge_i] + dist_kj * cluster_size[merge_j]) / 
                   (cluster_size[merge_i] + cluster_size[merge_j]);
    }
    
    // Update distance in BOTH triangles to maintain symmetry
    if (k < merge_i) {
        matrix[merge_i * n + k] = new_dist;
        matrix[k * n + merge_i] = new_dist;  // Symmetric update
    } else {
        matrix[k * n + merge_i] = new_dist;
        matrix[merge_i * n + k] = new_dist;  // Symmetric update
    }
}

// Kernel to mark cluster as inactive
__global__ void deactivate_cluster_kernel(
    int *active_points,
    int cluster_idx
) {
    active_points[cluster_idx] = 0;
}

// Kernel to update cluster size after merge
__global__ void update_cluster_size_kernel(
    int *cluster_size,
    int merge_i,
    int merge_j
) {
    cluster_size[merge_i] += cluster_size[merge_j];
}

extern "C" {

// Initialize GPU clustering with matrix transfer
int cuda_init_clustering_c(
    const double *h_matrix,
    const int *h_active_points,
    const int *h_cluster_size,
    int n
) {
    if (gpu_state != nullptr) {
        fprintf(stderr, "Error: GPU clustering already initialized\n");
        return -1;
    }
    
    gpu_state = new GPUClusterState();
    gpu_state->n = n;
    
    // Allocate device memory
    cudaMalloc((void**)&gpu_state->d_matrix, n * n * sizeof(double));
    cudaMalloc((void**)&gpu_state->d_active_points, n * sizeof(int));
    cudaMalloc((void**)&gpu_state->d_cluster_size, n * sizeof(int));
    cudaMalloc((void**)&gpu_state->d_temp_dist, n * sizeof(double));
    cudaMalloc((void**)&gpu_state->d_temp_i, n * sizeof(int));
    cudaMalloc((void**)&gpu_state->d_temp_j, n * sizeof(int));
    cudaMalloc((void**)&gpu_state->d_result_i, sizeof(int));
    cudaMalloc((void**)&gpu_state->d_result_j, sizeof(int));
    cudaMalloc((void**)&gpu_state->d_result_dist, sizeof(double));
    
    // Transfer data to GPU once
    cudaMemcpy(gpu_state->d_matrix, h_matrix, n * n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_state->d_active_points, h_active_points, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_state->d_cluster_size, h_cluster_size, n * sizeof(int), cudaMemcpyHostToDevice);
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA init error: %s\n", cudaGetErrorString(err));
        return -1;
    }
    
    return 0;
}

// Find minimum pair and perform merge on GPU
int cuda_find_and_merge_c(
    int *h_result_i,
    int *h_result_j,
    double *h_result_dist,
    int linkage_type
) {
    if (gpu_state == nullptr) {
        fprintf(stderr, "Error: GPU clustering not initialized\n");
        return -1;
    }
    
    int n = gpu_state->n;
    int threads_per_block = 256;
    int blocks = (n + threads_per_block - 1) / threads_per_block;
    
    // Find minimum pair
    find_min_pair_kernel<<<blocks, threads_per_block>>>(
        gpu_state->d_matrix,
        gpu_state->d_active_points,
        n,
        gpu_state->d_temp_dist,
        gpu_state->d_temp_i,
        gpu_state->d_temp_j
    );
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA find kernel error: %s\n", cudaGetErrorString(err));
        return -1;
    }
    
    // Reduce to global minimum
    reduce_min_kernel<<<1, 1>>>(
        gpu_state->d_temp_dist,
        gpu_state->d_temp_i,
        gpu_state->d_temp_j,
        n,
        gpu_state->d_result_i,
        gpu_state->d_result_j,
        gpu_state->d_result_dist
    );
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA reduce kernel error: %s\n", cudaGetErrorString(err));
        return -1;
    }
    
    cudaDeviceSynchronize();
    
    // Copy results back
    cudaMemcpy(h_result_i, gpu_state->d_result_i, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_result_j, gpu_state->d_result_j, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_result_dist, gpu_state->d_result_dist, sizeof(double), cudaMemcpyDeviceToHost);
    
    // If valid pair found, update distances and deactivate merged cluster
    if (*h_result_i >= 0 && *h_result_j >= 0) {
        // Update cluster size before computing new distances
        update_cluster_size_kernel<<<1, 1>>>(
            gpu_state->d_cluster_size,
            *h_result_i,
            *h_result_j
        );
        
        cudaDeviceSynchronize();
        
        update_distances_kernel<<<blocks, threads_per_block>>>(
            gpu_state->d_matrix,
            n,
            *h_result_i,
            *h_result_j,
            linkage_type,
            gpu_state->d_cluster_size
        );
        
        deactivate_cluster_kernel<<<1, 1>>>(
            gpu_state->d_active_points,
            *h_result_j
        );
        
        cudaDeviceSynchronize();
    }
    
    return 0;
}

// Cleanup GPU resources
void cuda_finalize_clustering_c() {
    if (gpu_state == nullptr) return;
    
    cudaFree(gpu_state->d_matrix);
    cudaFree(gpu_state->d_active_points);
    cudaFree(gpu_state->d_cluster_size);
    cudaFree(gpu_state->d_temp_dist);
    cudaFree(gpu_state->d_temp_i);
    cudaFree(gpu_state->d_temp_j);
    cudaFree(gpu_state->d_result_i);
    cudaFree(gpu_state->d_result_j);
    cudaFree(gpu_state->d_result_dist);
    
    delete gpu_state;
    gpu_state = nullptr;
}

}
