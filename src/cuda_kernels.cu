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

// Kernel to merge cluster indexes - copies members from cluster_j to cluster_i
// Fortran uses column-major order: array(row, col) -> offset = (row-1) + (col-1)*n
__global__ void merge_cluster_indexes_kernel(
    int *cluster_indexes,
    int *cluster_count,
    int n,
    int merge_i,
    int merge_j
) {
    int old_count_i = cluster_count[merge_i];
    int count_j = cluster_count[merge_j];
    
    // Copy all members from cluster_j to cluster_i
    // Fortran column-major: cluster_indexes(row, col) = cluster_indexes[(row-1) + (col-1)*n]
    // Note: merge_i and merge_j are 0-based C indices, but stored values are 1-based for Fortran
    for (int k = 0; k < count_j; k++) {
        // Source: cluster_indexes(merge_j+1, k+1) in Fortran
        // Dest:   cluster_indexes(merge_i+1, old_count_i+k+1) in Fortran
        cluster_indexes[merge_i + (old_count_i + k) * n] = cluster_indexes[merge_j + k * n];
    }
    
    // Update cluster count
    cluster_count[merge_i] = old_count_i + count_j;
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

// Complete clustering on GPU - runs all iterations without CPU interaction
int cuda_complete_clustering_c(
    const double *h_matrix,
    const int *h_active_points,
    const int *h_cluster_size,
    int *h_cluster_indexes,        // Output: cluster membership (n*n)
    int *h_cluster_count,          // Output: number of elements per cluster (n)
    int *h_active_clusters,        // Output: which clusters are still active (n)
    int n,
    double dist_threshold,
    int linkage_type
) {
    // Initialize GPU state
    GPUClusterState state;
    state.n = n;
    
    // Allocate device memory
    cudaMalloc((void**)&state.d_matrix, n * n * sizeof(double));
    cudaMalloc((void**)&state.d_active_points, n * sizeof(int));
    cudaMalloc((void**)&state.d_cluster_size, n * sizeof(int));
    cudaMalloc((void**)&state.d_temp_dist, n * sizeof(double));
    cudaMalloc((void**)&state.d_temp_i, n * sizeof(int));
    cudaMalloc((void**)&state.d_temp_j, n * sizeof(int));
    cudaMalloc((void**)&state.d_result_i, sizeof(int));
    cudaMalloc((void**)&state.d_result_j, sizeof(int));
    cudaMalloc((void**)&state.d_result_dist, sizeof(double));
    
    // Additional memory for cluster tracking
    int *d_cluster_indexes, *d_cluster_count;
    cudaMalloc((void**)&d_cluster_indexes, n * n * sizeof(int));
    cudaMalloc((void**)&d_cluster_count, n * sizeof(int));
    
    // Transfer initial data to GPU
    cudaMemcpy(state.d_matrix, h_matrix, n * n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(state.d_active_points, h_active_points, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(state.d_cluster_size, h_cluster_size, n * sizeof(int), cudaMemcpyHostToDevice);
    
    // Initialize cluster_indexes: each cluster starts with just itself
    // Fortran uses column-major and 1-based indexing
    // cluster_indexes(i, 1) = i in Fortran -> h_init_indexes[i-1] = i in C
    int *h_init_indexes = new int[n * n]();
    int *h_init_count = new int[n];
    for (int i = 0; i < n; i++) {
        h_init_indexes[i] = i + 1;  // Store 1-based index for Fortran (column 0, row i)
        h_init_count[i] = 1;
    }
    cudaMemcpy(d_cluster_indexes, h_init_indexes, n * n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cluster_count, h_init_count, n * sizeof(int), cudaMemcpyHostToDevice);
    delete[] h_init_indexes;
    delete[] h_init_count;
    
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA complete clustering init error: %s\n", cudaGetErrorString(err));
        return -1;
    }
    
    // Main clustering loop on CPU (coordinates GPU work)
    int threads_per_block = 256;
    int blocks = (n + threads_per_block - 1) / threads_per_block;
    int remaining_clusters = n;
    int merge_counter = 0;
    int max_iterations = n - 1;  // Maximum possible merges
    
    // Track merge history for CPU-side cluster index reconstruction
    int *merge_history_i = new int[max_iterations];
    int *merge_history_j = new int[max_iterations];
    
    while (remaining_clusters > 1 && merge_counter < max_iterations) {
        // Find minimum pair on GPU
        find_min_pair_kernel<<<blocks, threads_per_block>>>(
            state.d_matrix,
            state.d_active_points,
            n,
            state.d_temp_dist,
            state.d_temp_i,
            state.d_temp_j
        );
        
        reduce_min_kernel<<<1, 1>>>(
            state.d_temp_dist,
            state.d_temp_i,
            state.d_temp_j,
            n,
            state.d_result_i,
            state.d_result_j,
            state.d_result_dist
        );
        
        cudaDeviceSynchronize();
        
        // Copy results to check threshold
        int result_i, result_j;
        double result_dist;
        cudaMemcpy(&result_i, state.d_result_i, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&result_j, state.d_result_j, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&result_dist, state.d_result_dist, sizeof(double), cudaMemcpyDeviceToHost);
        
        // Check termination conditions
        if (result_i < 0 || result_j < 0 || result_dist > dist_threshold) {
            break;
        }
        
        // Store merge history for cluster index reconstruction
        merge_history_i[merge_counter] = result_i;
        merge_history_j[merge_counter] = result_j;
        
        // Merge clusters on GPU
        update_cluster_size_kernel<<<1, 1>>>(
            state.d_cluster_size,
            result_i,
            result_j
        );
        
        // Merge cluster membership on GPU
        merge_cluster_indexes_kernel<<<1, 1>>>(
            d_cluster_indexes,
            d_cluster_count,
            n,
            result_i,
            result_j
        );
        
        update_distances_kernel<<<blocks, threads_per_block>>>(
            state.d_matrix,
            n,
            result_i,
            result_j,
            linkage_type,
            state.d_cluster_size
        );
        
        deactivate_cluster_kernel<<<1, 1>>>(
            state.d_active_points,
            result_j
        );
        
        cudaDeviceSynchronize();
        
        remaining_clusters--;
        merge_counter++;
        
        // Optional progress reporting
        if (merge_counter % 100 == 0) {
            printf("GPU clustering: %d clusters remaining\n", remaining_clusters);
        }
    }
    
    printf("GPU clustering completed: %d merges, %d clusters remaining\n", merge_counter, remaining_clusters);
    
    // Transfer results back to host
    cudaMemcpy(h_active_clusters, state.d_active_points, n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cluster_count, d_cluster_count, n * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cluster_indexes, d_cluster_indexes, n * n * sizeof(int), cudaMemcpyDeviceToHost);
    
    // Copy back cluster sizes for use in post-processing
    int *h_cluster_size_out = new int[n];
    cudaMemcpy(h_cluster_size_out, state.d_cluster_size, n * sizeof(int), cudaMemcpyDeviceToHost);
    
    delete[] h_cluster_size_out;
    delete[] merge_history_i;
    delete[] merge_history_j;
    
    // Free device memory
    cudaFree(state.d_matrix);
    cudaFree(state.d_active_points);
    cudaFree(state.d_cluster_size);
    cudaFree(state.d_temp_dist);
    cudaFree(state.d_temp_i);
    cudaFree(state.d_temp_j);
    cudaFree(state.d_result_i);
    cudaFree(state.d_result_j);
    cudaFree(state.d_result_dist);
    cudaFree(d_cluster_indexes);
    cudaFree(d_cluster_count);
    
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA complete clustering error: %s\n", cudaGetErrorString(err));
        return -1;
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
