# Unit Testing Guide for Clustering Program

## Overview

This directory contains comprehensive unit tests for the Fortran subroutines in the clustering program. The tests verify correctness of mathematical operations, data structures, and algorithmic components.

## Test Files

### 1. **test_maths.f90**
Tests mathematical operations in the `maths` module:
- **Cross product calculation** (`cross` subroutine) - 4 tests
  - Standard basis vector products (i×j=k)
  - Anti-commutativity (v×w = -w×v)
  - Self-product is zero (v×v=0)
  - Arbitrary vector products
- **Coordinate transformations** (`update_complex` subroutine) - 7 tests
  - Identity transformations
  - Pure translation
  - XC1-only transformation
  - XC2-only transformation
  - Rotation-only transformation
  - Combined transformations
  - Verification of matrix operations
- **Angle calculations** (`vectors_angle_3D`, `vectors_angle_2D`) - 8 tests
  - 3D angle calculations with various vector combinations
  - 2D angle calculations with perpendicular and parallel vectors
- **RMSD calculations** (`rmsd` subroutine) - 3 tests
  - Identical structures (zero RMSD)
  - Translated structures
  - Rotated structures
- **Center of geometry** (`calculate_cog` subroutine) - 5 tests
  - Simple geometries
  - Symmetric point clouds
  - Offset coordinates
- **Distance calculations** (`calculate_distance` subroutine) - 4 tests
  - Unit vectors
  - Arbitrary vectors
  - Zero distance cases

**Total: 31 tests**

### 2. **test_matrix.f90**
Tests matrix operations in the `mod_matrix` module:
- **Matrix file I/O** (`write_matrix`, `read_matrix`) - 2 tests
  - Small matrices (5×5)
  - Symmetric matrix properties
- **Array file I/O** (`write_array`, `read_array`) - 2 tests
  - Single and double precision arrays
  - Large array serialization
- **Z-coordinate matrix calculation** (`matrix_z_coord`) - 2 tests
  - Identity transformation (coordinates unchanged)
  - Translation transformations
- **RMSD matrix calculation** (`matrix_rmsd_calc`) - 2 tests
  - Symmetry verification
  - Diagonal properties (zeros on diagonal)

**Total: 8 tests**

### 3. **test_threshold.f90**
Tests threshold array operations in the `mod_threshold` module:
- **Z-coordinate array calculation** (`array_z_coord`) - 2 tests
  - Identity transformation
  - Translation effects
- **Minimum distance calculation** (`array_atoms_dist`) - 2 tests
  - Distance computations
  - Minimum value detection
- **Angle calculations** (`array_angle`) - 2 tests
  - 3D angle calculations across arrays
  - Multiple vector pair angles
- **Array sorting** (`sort_array`) - 3 tests
  - Ascending order sorting
  - Index tracking
  - Large array handling

**Total: 9 tests**

### 4. **test_pdb.f90**
Tests PDB file structures and operations in the `mod_pdb` module:
- **PDB type initialization** - 3 tests
  - File structure allocation
  - Atom array creation
  - Residue array creation
- **Atom properties** - 3 tests
  - Coordinate assignment
  - Name and residue information
  - Occupancy and B-factor values
- **Residue management** - 2 tests
  - Residue property assignment
  - Atom array allocation within residues

**Total: 8 tests**

### 5. **test_clustering.f90**
Tests clustering algorithms in the `mod_clust_algorithm` module:
- **Clustering statistics** (`calculate_mean_distance`, `calculate_std_deviation`) - 2 tests
  - Mean distance calculations from clusters
  - Standard deviation computations
- **Simple clustering algorithm** (`simple_clustering`) - 2 tests
  - Tight cluster pair detection
  - Identical point handling
- **Cluster output** (`write_clusters`) - 1 test
  - File output verification
  - Cluster format validation

**Total: 5 tests**

### 6. **test_assoc.f90** (NEW)
Tests association/complexes file reading in the `mod_assoc` module:
- **File reading** (`read_assoc_file`) - 6 tests
  - Successfully read complexes from files
  - Correct structure allocation
  - Data verification (natom, residues, chains)
  - Header line parsing
  - Comment line handling
  - Sequential file reading
- **Object allocation** (`allocate_assoc_object`) - 3 tests
  - Memory allocation
  - Array initialization
  - Size verification
- **File parsing** (`size_assoc`) - 2 tests
  - Line counting with comment handling
  - Empty file handling

**Total: 11 tests**

### 7. **test_read_input.f90** (NEW)
Tests read_input module wrapper functions:
- **Association file reading** (`read_input_assoc`) - 7 tests
  - File reading and structure creation
  - Data verification (natom, residues, chains)
  - Header and comment handling
  - Sequential file operations
  - Memory allocation
  - File existence checking
  - Data integrity validation
- **Error handling** - 3 tests
  - Non-existent file handling
  - Invalid file format detection
  - Sequential read verification

**Total: 10 tests**
  - Occupancy and B-factor values
- **Residue management**
  - Residue property assignment
  - Atom array allocation within residues

## Compilation

### Compile individual tests:

```bash
# Test mathematical operations (requires maths.f90 and mod_pdb.f90)
gfortran -fopenmp -o test_maths test_maths.f90 ../src/maths.o ../src/mod_pdb.o

# Test matrix operations (requires mod_matrix.f90, maths.f90, mod_pdb.f90)
gfortran -fopenmp -o test_matrix test_matrix.f90 ../src/mod_matrix.o ../src/maths.o ../src/mod_pdb.o

# Test threshold operations (requires mod_threshold.f90, maths.f90)
gfortran -fopenmp -o test_threshold test_threshold.f90 ../src/mod_threshold.o ../src/maths.o

# Test PDB operations (requires mod_pdb.f90)
gfortran -o test_pdb test_pdb.f90 ../src/mod_pdb.o

# Test clustering operations (requires mod_clust_algorithm.f90, read_input.f90, etc.)
gfortran -fopenmp -o test_clustering test_clustering.f90 ../src/mod_clust_algorithm.o ../src/read_input.o ../src/mod_pdb.o ../src/mod_assoc.o ../src/maths.o

# Test association file reading (requires mod_assoc.f90)
gfortran -o test_assoc test_assoc.f90 ../src/mod_assoc.o

# Test read_input module (requires read_input.f90, mod_pdb.f90, mod_assoc.f90)
gfortran -o test_read_input test_read_input.f90 ../src/read_input.o ../src/mod_pdb.o ../src/mod_assoc.o
```

### Compile all tests at once:

```bash
# Using make (from src directory)
make test

# Or using the test runner script (from unit_tests directory)
bash run_all_tests.sh
```

## Running Tests

### Run individual tests:

```bash
./test_maths
./test_matrix
./test_threshold
./test_pdb
```

### Run all tests:

```bash
# Create a simple test runner script
bash run_all_tests.sh
```

## Example Output

```
=========================================
UNIT TESTS FOR MATHS MODULE
=========================================

Testing cross product...
  [PASS] Test 1: i x j = k
  [PASS] Test 2: j x i = -k
  [PASS] Test 3: v x v = 0
  [PASS] Test 4: Arbitrary vectors

Testing coordinate transformation (update_complex)...
  [PASS] Test 1: Identity transformation
  [PASS] Test 2: Pure translation

=========================================
TEST SUMMARY
=========================================
Total tests:  6
Passed:       6
Failed:       0
=========================================
```

## Test Coverage

| Module | Subroutine | Tests | Status |
|--------|-----------|-------|--------|
| maths | cross | 4 | Yes |
| maths | update_complex | 7 | Yes |
| maths | vectors_angle_3D | 4 | Yes |
| maths | vectors_angle_2D | 4 | Yes |
| maths | rmsd | 3 | Yes |
| maths | calculate_cog | 5 | Yes |
| maths | calculate_distance | 4 | Yes |
| mod_matrix | write_matrix / read_matrix | 2 | Yes |
| mod_matrix | write_array / read_array | 2 | Yes |
| mod_matrix | matrix_z_coord | 2 | Yes |
| mod_matrix | matrix_rmsd_calc | 2 | Yes |
| mod_threshold | array_z_coord | 2 | Yes |
| mod_threshold | array_atoms_dist | 2 | Yes |
| mod_threshold | array_angle | 2 | Yes |
| mod_threshold | sort_array | 3 | Yes |
| mod_pdb | PDB allocation | 3 | Yes |
| mod_pdb | atom properties | 3 | Yes |
| mod_pdb | residue management | 2 | Yes |
| mod_clust_algorithm | clustering statistics | 2 | Yes |
| mod_clust_algorithm | simple_clustering | 2 | Yes |
| mod_clust_algorithm | write_clusters | 1 | Yes |
| mod_assoc | read_assoc_file | 6 | Yes |
| mod_assoc | allocate_assoc_object | 3 | Yes |
| mod_assoc | size_assoc | 2 | Yes |
| read_input | read_input_assoc | 7 | Yes |
| read_input | error_handling | 3 | Yes |

**Total: 82 tests**

## Adding New Tests

To add tests for additional subroutines:

1. **Create a new test program** following the template:
   ```fortran
   program test_my_module
       implicit none
       integer :: num_tests, num_passed, num_failed
       
       num_tests = 0
       num_passed = 0
       num_failed = 0
       
       ! Call test subroutines
       call test_my_subroutine(num_tests, num_passed, num_failed)
       
       ! Print summary
       if (num_failed > 0) stop 1
   end program
   ```

2. **Implement test subroutines** with clear assertions:
   ```fortran
   subroutine test_my_subroutine(num_tests, num_passed, num_failed)
       ! Test implementation
   end subroutine
   ```

3. **Add compilation rule** to Makefile

4. **Document test** in this README

## Best Practices

1. **Test isolation**: Each test should be independent
2. **Clear naming**: Use descriptive test names
3. **Tolerance handling**: Use appropriate floating-point tolerances (typically 1d-10 or 1d-8)
4. **Edge cases**: Test boundary conditions and special cases
5. **Performance**: Consider testing with realistic data sizes
6. **Documentation**: Include test purpose in comments

## Troubleshooting

### Compilation errors:
- Ensure all required modules are compiled first
- Check for missing `-fopenmp` flag if using OpenMP
- Verify Fortran standard compatibility

### Runtime failures:
- Check for NaN/Inf values in calculations
- Verify array bounds and allocation
- Review floating-point tolerance settings

### Performance issues:
- Monitor memory usage for large matrices
- Use OpenMP pragmas for parallelized loops
- Profile with `gprof` or similar tools

## References

- [Fortran Standard Library](https://stdlib.fortran-lang.org/)
- [OpenMP Fortran API](https://www.openmp.org/)
- [PDB File Format](https://www.wwpdb.org/documentation/file-format)

---

**Author**: Abraham Muñiz-Chicharro  
**Version**: 1.0  
**Date**: February 2026
