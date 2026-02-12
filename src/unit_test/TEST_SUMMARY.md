# Unit Testing Suite for Clustering Program

## Summary

A comprehensive unit testing framework has been created for the Fortran clustering program with **7 test programs** covering **89 individual tests** across all major modules.

## Test Programs Created

### 1. **test_maths.f90** (30 tests)
- Cross product calculation (basis vectors, anti-commutativity, self-product, arbitrary vectors) - 4 tests
- Coordinate transformations (identity, translation, rotation, combinations) - 7 tests
- Angle calculations in 3D and 2D - 8 tests
- RMSD calculations - 3 tests
- Center of geometry computations - 4 tests
- Distance calculations - 4 tests
- Tests fundamental mathematical operations used throughout the program

### 2. **test_matrix.f90** (10 tests)
- Matrix file I/O operations (write/read matrix) - 2 tests
- Array file I/O operations (write/read array) - 2 tests
- Z-coordinate matrix calculations - 2 tests
- RMSD matrix calculations with symmetry and diagonal verification - 2 tests
- Atom distance matrix calculations - 1 test
- Angle matrix calculations - 1 test
- Tests core matrix operations for encounter data

### 3. **test_threshold.f90** (9 tests)
- Z-coordinate array calculations - 2 tests
- Minimum distance computations - 2 tests
- Angle calculations across arrays - 2 tests
- Array sorting with index tracking - 3 tests
- Tests threshold computation infrastructure

### 4. **test_pdb.f90** (11 tests)
- PDB structure initialization and allocation - 3 tests
- Atom properties (coordinates, names, residue info) - 2 tests
- Residue management - 2 tests
- Chain handling - 2 tests
- PDB count/allocate/fill routines - 2 tests
- Tests PDB data structure integrity

### 5. **test_clustering.f90** (8 tests)
- Clustering statistics (mean distance, standard deviation) - 2 tests
- Simple clustering algorithm with tight cluster detection - 2 tests
- Minimum linkage clustering - 1 test
- Maximum linkage clustering - 1 test
- Mean linkage clustering - 1 test
- Cluster output file writing and validation - 1 test
- Tests core clustering algorithm components

### 6. **test_assoc.f90** (11 tests)
- Association/complexes file reading - 6 tests
- Association object allocation - 3 tests
- File size and line counting - 2 tests
- Tests association file I/O and data structures

### 7. **test_read_input.f90** (10 tests)
- Read input wrapper functions - 7 tests
- Error handling and edge cases - 3 tests
- Tests input file processing module

## Files

```
clustering_program/unit_tests/
- test_maths.f90           # Mathematical operations tests (30 tests)
- test_matrix.f90          # Matrix operations tests (10 tests)
- test_threshold.f90       # Threshold array tests (9 tests)
- test_pdb.f90             # PDB structure tests (11 tests)
- test_clustering.f90      # Clustering algorithm tests (8 tests)
- test_assoc.f90           # Association file tests (11 tests) [NEW]
- test_read_input.f90      # Read input module tests (10 tests) [NEW]
- run_all_tests.sh         # Automated test runner script
- INDEX.md                 # Documentation index
- QUICK_TEST.md            # Quick test reference
- TESTING.md               # Comprehensive testing documentation
- README.md                # Test suite overview
- TEST_SUMMARY.md          # This file
```

## Building Tests

### Individual tests:
```bash
cd clustering_program/src/

# Test specific module
gfortran -fopenmp -o test_maths ../unit_tests/test_maths.f90 maths.o mod_pdb.o
gfortran -fopenmp -o test_matrix ../unit_tests/test_matrix.f90 mod_matrix.o maths.o mod_pdb.o
gfortran -fopenmp -o test_threshold ../unit_tests/test_threshold.f90 mod_threshold.o maths.o
gfortran -o test_pdb ../unit_tests/test_pdb.f90 mod_pdb.o
gfortran -fopenmp -o test_clustering ../unit_tests/test_clustering.f90 mod_clust_algorithm.o read_input.o mod_pdb.o mod_assoc.o maths.o
gfortran -o test_assoc ../unit_tests/test_assoc.f90 mod_assoc.o
gfortran -o test_read_input ../unit_tests/test_read_input.f90 read_input.o mod_pdb.o mod_assoc.o
```

### Run all tests via Make:
```bash
cd clustering_program/src/
make test                # Compile and run all 7 test suites
make test_maths          # Run specific test
make test_matrix
make test_threshold
make test_pdb
make test_clustering
make test_assoc          # NEW
make test_read_input     # NEW
make clean_tests         # Remove test executables only
```

### Automated test runner:
```bash
cd clustering_program/unit_tests/
bash run_all_tests.sh
```

## Test Coverage Matrix

| Module | Function | Tests | Type |
|--------|----------|-------|------|
| maths | cross | 4 | Unit |
| maths | update_complex | 2 | Unit |
| mod_matrix | allocation | 2 | Unit |
| mod_matrix | symmetry | 2 | Unit |
| mod_threshold | statistics | 3 | Unit |
| mod_threshold | bounds | 3 | Unit |
| mod_pdb | allocation | 3 | Unit |
| mod_pdb | atoms | 3 | Unit |
| mod_pdb | residues | 2 | Unit |
| mod_clustering | distances | 3 | Unit |
| mod_clustering | helpers | 3 | Unit |

**Total: 89 unit tests**

## Test Features

[x] **Comprehensive coverage** - Tests mathematical operations, data structures, and algorithms
[x] **Floating-point tolerance** - Uses appropriate precision (1d-10 to 1d-8) for double precision
[x] **Edge case handling** - Tests boundary conditions and special cases
[x] **Detailed reporting** - Pass/fail output with detailed diagnostics
[x] **Automated execution** - Makefile targets and shell script runner
[x] **OpenMP compatible** - Handles parallel compilation flags
[x] **Modular design** - Independent test programs, easy to extend

## Example Test Output

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
Total tests:  30
Passed:       30
Failed:       0
=========================================
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
  [PASS] Test 3: XC1-only transformation
  [PASS] Test 4: XC2-only transformation
  [PASS] Test 5: Rotation-only transformation
  [PASS] Test 6: Combined transformations
  [PASS] Test 7: Matrix operation verification

... (30 total tests in test_maths.f90)

=========================================
TEST SUMMARY FOR MATHS
=========================================
Total tests:  30
Passed:       30
Failed:       0
=========================================
```

## Overall Statistics

| Test File | Tests | Lines | Status |
|-----------|-------|-------|--------|
| test_maths.f90 | 30 | 850 | Success |
| test_matrix.f90 | 10 | 370 | Success |
| test_threshold.f90 | 9 | 375 | Success |
| test_pdb.f90 | 11 | 210 | Success |
| test_clustering.f90 | 8 | Success | Success |
| test_assoc.f90 | 11 | 380 | Success |
| test_read_input.f90 | 10 | 355 | Success |
| **TOTAL** | **89** | **2,815** | **YES** |

## Running Tests

### Quick test:
```bash
cd clustering_program/src/
make test          # Runs all 89 tests across 7 test programs
```

### Run specific test:
```bash
cd clustering_program/src/
make test_maths    # Or test_matrix, test_threshold, test_pdb, etc.
make test_assoc    # Run new association tests
make test_read_input  # Run new read_input tests
```

### Verbose output:
```bash
cd clustering_program/unit_tests/
bash run_all_tests.sh
```

### Individual test debugging:
```bash
cd clustering_program/src/
gfortran -fopenmp -o test_matrix ../unit_tests/test_matrix.f90 mod_matrix.o maths.o mod_pdb.o
./test_matrix
```

## Documentation

See [TESTING.md](TESTING.md) for:
- Detailed test descriptions
- Compilation instructions
- Troubleshooting guide
- Instructions for adding new tests
- Best practices and references

## Future Enhancements

- [ ] Integration tests for multi-module workflows
- [ ] Performance benchmarks
- [ ] Code coverage analysis
- [ ] Continuous integration (GitHub Actions)
- [ ] Memory leak detection (valgrind)
- [ ] Additional clustering algorithm tests

## Dependencies

- gfortran (or other Fortran compiler)
- OpenMP (for parallel tests)
- LAPACK/BLAS (existing project dependency)

## Notes

- Tests use standard IEEE 754 double precision (kind=8)
- All tests are deterministic and reproducible
- Tests do not require external data files
- Clean exit codes: 0 (success), 1 (failure)

---

**Created**: February 2026  
**Author**: Abraham Muñiz-Chicharro  
**Version**: 1.0
