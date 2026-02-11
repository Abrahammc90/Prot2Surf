# Quick Reference: Running Unit Tests

## TL;DR

```bash
cd clustering_program/src/
make test          # Compile and run all 82 tests
```

## Test Commands

| Command | Description |
|---------|-------------|
| `make test` | Compile and run all 82 tests (7 test programs) |
| `make test_maths` | Test mathematical operations (31 tests) |
| `make test_matrix` | Test matrix operations (8 tests) |
| `make test_threshold` | Test threshold arrays (9 tests) |
| `make test_pdb` | Test PDB structures (8 tests) |
| `make test_clustering` | Test clustering algorithms (5 tests) |
| `make test_assoc` | Test association file reading (11 tests) ✓ NEW |
| `make test_read_input` | Test read_input module (10 tests) ✓ NEW |
| `make clean_tests` | Remove test executables |
| `bash run_all_tests.sh` | Run all tests with detailed output |

## Test Files

| File | Tests | Lines | Purpose |
|------|-------|-------|---------|
| `test_maths.f90` | 31 | 850 | Cross product, transformations, RMSD, angles, COG, distance |
| `test_matrix.f90` | 8 | 370 | Matrix I/O, array I/O, Z-coord, RMSD |
| `test_threshold.f90` | 9 | 375 | Z-coord, min distance, angles, sorting |
| `test_pdb.f90` | 8 | 210 | PDB types, atoms, residues |
| `test_clustering.f90` | 5 | 275 | Statistics, clustering, output writing |
| `test_assoc.f90` | 11 | 380 | File reading, allocation, parsing ✓ NEW |
| `test_read_input.f90` | 10 | 355 | Input wrapper, error handling ✓ NEW |

## Expected Output

All tests pass with output like:

```
[PASS] Test 1: Test description
[PASS] Test 2: Another test
=========================================
TEST SUMMARY
=========================================
Total tests:  N
Passed:       N
Failed:       0
=========================================
```

## Individual Compilation

```bash
# Manual compilation (if not using Make)
gfortran -o test_maths test_maths.f90 maths.f90 mod_pdb.f90
./test_maths

gfortran -fopenmp -o test_matrix test_matrix.f90 mod_matrix.f90 maths.f90 mod_pdb.f90
./test_matrix

gfortran -fopenmp -o test_threshold test_threshold.f90
./test_threshold

gfortran -o test_pdb test_pdb.f90 mod_pdb.f90
./test_pdb

gfortran -fopenmp -o test_clustering test_clustering.f90
./test_clustering
```

## Documentation

- [TEST_SUMMARY.md](TEST_SUMMARY.md) - Overview of test suite
- [TESTING.md](TESTING.md) - Detailed testing guide

## Key Features

[x] 82 comprehensive tests (updated from 28)
[x] 7 test programs (added test_assoc.f90 and test_read_input.f90)
[x] All major modules covered
[x] OpenMP compatible
[x] No external data files needed
[x] Fast execution (< 2 seconds)
[x] Clear pass/fail reporting
[x] Easy to extend

---

**Pro Tip**: Use `make test 2>&1 | tee test_results.log` to save test output to a file.
