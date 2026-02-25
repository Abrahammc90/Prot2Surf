#include <stdio.h>
#include <float.h>
#include <string.h>
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
    double local_min_dist = DBL_MAX;
    int local_min_i = -1;
    
    // Search for minimum in this column (lower triangle only)
    for (int i = 0; i < j; i++) {
        if (active_points[i]) {
            // Column-major: matrix(i,j) is at matrix[j*n + i]
            double dist = matrix[j * n + i];
            if (dist < local_min_dist) {
                local_min_dist = dist;
                local_min_i = i;
            }
        }
    }
    
    // Store this thread's result
    temp_dist[j] = local_min_dist;
    temp_i[j] = local_min_i;
    temp_j[j] = j;
}

// Kernel to find minimum in each column for a batch of columns
__global__ void find_min_pair_batch_kernel(
    const double *matrix,
    const int *active_points,
    int n,
    int col_start,
    int cols,
    double *temp_dist,
    int *temp_i,
    int *temp_j
) {
    int local_col = blockIdx.x * blockDim.x + threadIdx.x;
    if (local_col >= cols) return;

    int j = col_start + local_col;
    if (j >= n) return;

    if (!active_points[j]) {
        temp_dist[local_col] = DBL_MAX;
        temp_i[local_col] = -1;
        temp_j[local_col] = -1;
        return;
    }

    const double *col_ptr = matrix + ((size_t)local_col * (size_t)n);
    double local_min_dist = DBL_MAX;
    int local_min_i = -1;

    for (int i = 0; i < j; i++) {
        if (active_points[i]) {
            double dist = col_ptr[i];
            if (dist < local_min_dist) {
                local_min_dist = dist;
                local_min_i = i;
            }
        }
    }

    temp_dist[local_col] = local_min_dist;
    temp_i[local_col] = local_min_i;
    temp_j[local_col] = j;
}

// Kernel to reduce to global minimum
__global__ void reduce_min_kernel(
    const double *temp_dist,
    const int *temp_i,
    const int *temp_j,
    int n,
    int *min_i_out,
    int *min_j_out,
    double *min_dist_out
) {
    double min_dist = DBL_MAX;
    int min_i = -1;
    int min_j = -1;
    
    for (int idx = 0; idx < n; idx++) {
        if (temp_dist[idx] < min_dist) {
            min_dist = temp_dist[idx];
            min_i = temp_i[idx];
            min_j = temp_j[idx];
        }
    }
    
    *min_dist_out = min_dist;
    *min_i_out = min_i;
    *min_j_out = min_j;
}

// Kernel to update distances after merging clusters min_i and min_j
// linkage_type: 0=min, 1=max, 2=mean
__global__ void update_distances_kernel(
    double *matrix,
    int n,
    int min_i,
    int min_j,
    int linkage_type,
    const int *cluster_size
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n || i == min_i || i == min_j) return;
    
    // Get distances from i to clusters min_i and min_j
    double dist_to_min_i, dist_to_min_j;
    
    // Access matrix element (i, min_i) - column-major
    if (i < min_i) {
        dist_to_min_i = matrix[min_i * n + i];  // i is row, min_i is column
    } else {
        dist_to_min_i = matrix[i * n + min_i];  // min_i is row, i is column
    }
    
    // Access matrix element (i, min_j) - column-major
    if (i < min_j) {
        dist_to_min_j = matrix[min_j * n + i];
    } else {
        dist_to_min_j = matrix[i * n + min_j];
    }
    
    // Compute new distance based on linkage
    double dist;
    if (linkage_type == 0) {  // min linkage
        dist = fmin(dist_to_min_i, dist_to_min_j);
    } else if (linkage_type == 1) {  // max linkage
        dist = fmax(dist_to_min_i, dist_to_min_j);
    } else {  // mean linkage - weighted average by cluster size
        dist = (dist_to_min_i * cluster_size[min_i] + dist_to_min_j * cluster_size[min_j]) / 
               (cluster_size[min_i] + cluster_size[min_j]);
    }
    
    // Update distance in BOTH triangles to maintain symmetry
    if (i < min_i) {
        matrix[min_i * n + i] = dist;
        matrix[i * n + min_i] = dist;  // Symmetric update
    } else {
        matrix[i * n + min_i] = dist;
        matrix[min_i * n + i] = dist;  // Symmetric update
    }
}

// Kernel to mark cluster as inactive
__global__ void deactivate_cluster_kernel(
    int *active_points,
    int min_j
) {
    active_points[min_j] = 0;
}

// Kernel to update cluster size after merge
__global__ void update_cluster_size_kernel(
    int *cluster_size,
    int min_i,
    int min_j
) {
    cluster_size[min_i] += cluster_size[min_j];
}

// Kernel to update cluster parent array after merge - O(n) tracking
__global__ void update_cluster_parent_kernel(
    int *cluster_parent,
    int n,
    int min_i,
    int min_j
) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < n && cluster_parent[k] == min_j) {
        cluster_parent[k] = min_i;
    }
}

extern "C" {

// Complete clustering on GPU - runs all iterations without CPU interaction
int cuda_matrix_clustering_c(
    const double *h_matrix,
    const int *h_active_points,
    const int *h_cluster_size,
    int *h_cluster_parent,         // Output: cluster parent for each point (n) - O(n) tracking
    int *h_cluster_count,          // Output: number of elements per cluster (n)
    int *h_active_clusters,        // Output: which clusters are still active (n)
    int n,
    double dist_threshold,
    int linkage_type
) {
    size_t free_bytes = 0, total_bytes = 0;
    cudaError_t err = cudaMemGetInfo(&free_bytes, &total_bytes);
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA memory query error: %s\n", cudaGetErrorString(err));
        return -1;
    }

    size_t available_bytes = (size_t)(free_bytes * 0.80);

    size_t n_size = (size_t)n;
    size_t n2_size = n_size * n_size;

    size_t full_bytes = n2_size * sizeof(double) +
                        n_size * sizeof(int) * 2 +
                        n_size * sizeof(double) +
                        n_size * sizeof(int) * 2 +
                        sizeof(int) * 2 + sizeof(double) +
                        n_size * sizeof(int);  // cluster_parent array - O(n) tracking

    
    if (full_bytes > available_bytes) {
        fprintf(stderr, "Error: Not enough GPU memory for clustering. Required: %zu bytes, Available: %zu bytes\n", full_bytes, available_bytes);
        return -1;
    }

    // use_batched = true;  // Disable batching for now - focus on GPU-resident path
    // Initialize GPU state
    GPUClusterState state;
    state.n = n;

    // Allocate device memory
    cudaMalloc((void**)&state.d_matrix, n_size * n_size * sizeof(double));
    cudaMalloc((void**)&state.d_active_points, n_size * sizeof(int));
    cudaMalloc((void**)&state.d_cluster_size, n_size * sizeof(int));
    cudaMalloc((void**)&state.d_temp_dist, n_size * sizeof(double));
    cudaMalloc((void**)&state.d_temp_i, n_size * sizeof(int));
    cudaMalloc((void**)&state.d_temp_j, n_size * sizeof(int));
    cudaMalloc((void**)&state.d_result_i, sizeof(int));
    cudaMalloc((void**)&state.d_result_j, sizeof(int));
    cudaMalloc((void**)&state.d_result_dist, sizeof(double));

    // Additional memory for cluster tracking - O(n) parent array
    int *d_cluster_parent;
    cudaMalloc((void**)&d_cluster_parent, n_size * sizeof(int));

    // Transfer initial data to GPU
    cudaMemcpy(state.d_matrix, h_matrix, n_size * n_size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(state.d_active_points, h_active_points, n_size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(state.d_cluster_size, h_cluster_size, n_size * sizeof(int), cudaMemcpyHostToDevice);

    // Initialize cluster_parent: each point starts in its own cluster
    // Fortran uses 1-based indexing
    int *h_init_parent = new int[n_size];
    for (int i = 0; i < n; i++) {
        h_init_parent[i] = i + 1;  // Store 1-based index for Fortran
    }
    cudaMemcpy(d_cluster_parent, h_init_parent, n_size * sizeof(int), cudaMemcpyHostToDevice);
    delete[] h_init_parent;

    err = cudaGetLastError();
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
        int min_i, min_j;
        double min_dist;
        cudaMemcpy(&min_i, state.d_result_i, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&min_j, state.d_result_j, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&min_dist, state.d_result_dist, sizeof(double), cudaMemcpyDeviceToHost);

        // Check termination conditions
        if (min_i < 0 || min_j < 0 || min_dist > dist_threshold) {
            break;
        }

        // Store merge history for cluster index reconstruction
        merge_history_i[merge_counter] = min_i;
        merge_history_j[merge_counter] = min_j;

        // Merge clusters on GPU
        // Update cluster parent on GPU - O(n) operation
        // NOTE: min_i/min_j are 0-based (from CUDA kernels),
        //       but cluster_parent stores 1-based indices (for Fortran).
        //       Convert to 1-based before updating.
        update_cluster_parent_kernel<<<blocks, threads_per_block>>>(
            d_cluster_parent,
            n,
            min_i + 1,  // Convert to 1-based for Fortran
            min_j + 1   // Convert to 1-based for Fortran
        );

        // Update distances BEFORE updating cluster_size so that
        // the mean linkage formula uses the original (pre-merge) sizes.
        update_distances_kernel<<<blocks, threads_per_block>>>(
            state.d_matrix,
            n,
            min_i,
            min_j,
            linkage_type,
            state.d_cluster_size
        );

        // Update cluster size AFTER distance update
        update_cluster_size_kernel<<<1, 1>>>(
            state.d_cluster_size,
            min_i,
            min_j
        );

        deactivate_cluster_kernel<<<1, 1>>>(
            state.d_active_points,
            min_j
        );

        cudaDeviceSynchronize();

        remaining_clusters--;
        merge_counter++;

        // Optional progress reporting
        if (remaining_clusters % 100 == 0) {
            printf("GPU clustering: %d clusters remaining\n", remaining_clusters);
        }
    }

    printf("GPU clustering completed: %d merges, %d clusters remaining\n", merge_counter, remaining_clusters);

    // Transfer results back to host
    cudaMemcpy(h_active_clusters, state.d_active_points, n_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cluster_count, state.d_cluster_size, n_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cluster_parent, d_cluster_parent, n_size * sizeof(int), cudaMemcpyDeviceToHost);

    // Copy back cluster sizes for use in post-processing
    int *h_cluster_size_out = new int[n_size];
    cudaMemcpy(h_cluster_size_out, state.d_cluster_size, n_size * sizeof(int), cudaMemcpyDeviceToHost);

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
    cudaFree(d_cluster_parent);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA complete clustering error: %s\n", cudaGetErrorString(err));
        return -1;
    }

    return 0;
}

}

