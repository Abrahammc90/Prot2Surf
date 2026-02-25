# CUDA Implementation for Clustering

## Overview
CUDA support has been implemented in the linkage subroutine to accelerate the minimum distance finding step during hierarchical clustering. This is the computational bottleneck in the clustering algorithm.

## Building with CUDA Support

### Prerequisites
- NVIDIA CUDA Toolkit (version 11.0 or later recommended)
- NVIDIA GPU with compute capability 3.0 or higher
- gfortran compiler with OpenMP support
- cuBLAS library (optional, for future enhancements)

### Building without CUDA (Default)
```bash
cd clustering_program/src
make clean
make all
```

### Building with CUDA
```bash
cd clustering_program/src
make clean
make all USE_CUDA=1
```

Note: The CUDA implementation will automatically detect if your GPU has the required compute capability.

## Usage

### Running with CUDA Acceleration
```bash
./bin/clust -complexes <file> -matrix <file> -output_name <name> -cuda
```

### Running without CUDA (CPU-based with OpenMP)
```bash
./bin/clust -complexes <file> -matrix <file> -output_name <name>
```

## Command-Line Arguments

- `-cuda`: Enable CUDA acceleration for minimum distance finding
  - Optional parameter
  - If not specified, CPU-based OpenMP acceleration is used instead
  - Requires the program to have been compiled with `USE_CUDA=1`

## Implementation Details

### Fortran Interface (`mod_cuda.f90`)
- Provides Fortran-callable interface to CUDA kernels
- Handles GPU memory allocation and deallocation
- Manages host-device data transfer
- Implements `cuda_find_min_pair()` subroutine

### CUDA Kernels (`cuda_kernels.cu`)
- `find_min_distance_kernel()`: GPU kernel to find minimum distance between active clusters
- `find_min_distance_cuda_()`: Fortran wrapper that:
  - Allocates GPU memory
  - Copies data to GPU
  - Launches kernel with appropriate grid/block configuration
  - Transfers results back to host
  - Cleans up GPU memory

### Modification to Clustering Algorithm (`mod_clust_algorithm.f90`)
- Added optional `use_cuda` parameter to `linkage_clustering()` subroutine
- Conditional execution: uses CUDA kernel if `use_cuda=.true.`, else uses OpenMP
- Maintains numerical equivalence with CPU version

### Main Program (`clust.f90`)
- Added `-cuda` command-line flag
- Passes `use_cuda` flag to clustering subroutine
- Updated help documentation

## Performance Considerations

### When to Use CUDA
- Large matrices (n > 10,000)
- Many iterations of clustering
- GPU is available on the system

### When CUDA May Not Help
- Small matrices (n < 1,000)
- Single use cases (startup/transfer overhead dominates)
- CPU cores are already fully utilized

## Troubleshooting

### Error: "CUDA files not found"
- Ensure CUDA Toolkit is installed
- Check CUDA_INCLUDES path in Makefile (default: `/usr/local/cuda/include`)
- Run: `which nvcc` to verify CUDA installation

### Error: "No CUDA device found" or "Device not compatible"
- Ensure NVIDIA GPU is connected
- Run: `nvidia-smi` to verify GPU is detected
- Check GPU compute capability matches `-arch` setting in Makefile

### Compilation Failures
- Start with `make clean` to remove old object files
- Try building without CUDA first: `make all`
- Check gfortran version: `gfortran --version`
- Check CUDA version: `nvcc --version`

## Future Enhancements

1. **Improved Global Reduction**: Current implementation uses a simple approach. Can be optimized with NVIDIA's reduction templates or Thrust library.
2. **Memory Optimization**: Use GPU memory optimization techniques (coalescing, shared memory patterns)
3. **Asynchronous Transfers**: Use CUDA streams for overlapping computation and communication
4. **Multi-GPU Support**: Extend to multiple GPUs for very large datasets
5. **Batch Processing**: Process multiple clustering iterations concurrently
6. **cuBLAS Integration**: Use cuBLAS for matrix operations if needed

## Files Modified/Created

### New Files
- `mod_cuda.f90` - Fortran module with CUDA interface
- `cuda_kernels.cu` - CUDA kernel implementation

### Modified Files
- `mod_clust_algorithm.f90` - Added CUDA support in linkage_clustering()
- `clust.f90` - Added `-cuda` command-line flag
- `Makefile` - Added CUDA compilation rules and linking

## Author
Abraham Muñiz-Chicharro

## References
- NVIDIA CUDA Programming Guide
- NVIDIA Fortran Interoperability Guide
- Thrust Library documentation
