# Unit Testing Framework

This directory contains the comprehensive unit testing framework for the clustering program.

## Contents

### Test Programs
- **test_maths.f90** - Tests for mathematical operations (31 tests)
- **test_matrix.f90** - Tests for matrix operations (8 tests)
- **test_threshold.f90** - Tests for threshold calculations (9 tests)
- **test_pdb.f90** - Tests for PDB file handling (8 tests)
- **test_clustering.f90** - Tests for clustering algorithms (5 tests)
- **test_assoc.f90** - Tests for association file reading (11 tests) ✓ NEW
- **test_read_input.f90** - Tests for read_input module (10 tests) ✓ NEW

### Test Automation
- **run_all_tests.sh** - Automated script to run all unit tests

### Documentation

| Document | Purpose | Best For |
|----------|---------|----------|
| **QUICK_TEST.md** | Quick reference guide | Getting started quickly |
| **TESTING.md** | Comprehensive testing guide | Understanding all details |
| **TEST_SUMMARY.md** | Test suite overview | High-level overview |
| **UNIT_TESTING_COMPLETE.md** | Completion summary | Project status |
| **INDEX.md** | File index and statistics | Navigation |

## Quick Start

### Run all tests
```bash
cd ../src
make test
```

### Run specific test
```bash
cd ../src
make test_maths    # or test_matrix, test_threshold, test_pdb, test_clustering
```

### Clean test executables
```bash
cd ../src
make clean_tests
```

## Where to Start?

1. **New to testing?** → Start with [QUICK_TEST.md](QUICK_TEST.md)
2. **Need details?** → Read [TESTING.md](TESTING.md)
3. **Want overview?** → See [TEST_SUMMARY.md](TEST_SUMMARY.md)
4. **Need file list?** → Check [INDEX.md](INDEX.md)
5. **Project status?** → Read [UNIT_TESTING_COMPLETE.md](UNIT_TESTING_COMPLETE.md)

## Test Coverage

| Module | Tests | Status |
|--------|-------|--------|
| maths.f90 | 31 | Yes |
| mod_matrix.f90 | 8 | Yes |
| mod_threshold.f90 | 9 | Yes |
| mod_pdb.f90 | 8 | Yes |
| mod_clust_algorithm.f90 | 5 | Yes |
| mod_assoc.f90 | 11 | Yes ✓ NEW |
| read_input.f90 | 10 | Yes ✓ NEW |
| **Total** | **82** | **Yes** |

## Build Requirements

- gfortran compiler
- LAPACK and BLAS libraries
- OpenMP support

## Verification

All 82 unit tests are defined and ready to run. The framework provides:
- Comprehensive test coverage across 7 test files
- Assertion helpers for test validation
- Pass/fail reporting with detailed output
- Error message display
- Command-line integration via Makefile
- Standalone bash script runner

## Notable Behavior

**Array Parameter in Clustering:**
- The `-array` option is **mandatory** for all clustering methods **except RMSD**
- When using `-matrix_type rmsd`, the array parameter is **optional**
- When an array is provided, complexes are sorted by the array values
- If array is not provided for non-RMSD methods, an error is raised
- Without an array, cluster averages and standard deviations are computed from the matrix

**Test File Coverage Updates (February 2026):**
- Added comprehensive tests for mod_assoc module (association file reading)
- Added comprehensive tests for read_input module (input file processing)
- Expanded maths module tests from 2 to 31 tests with detailed coordinate transformation coverage
- Updated all test documentation to reflect current test count (82 total tests)

For detailed information, see [TESTING.md](TESTING.md) and [UNIT_TESTING_COMPLETE.md](UNIT_TESTING_COMPLETE.md).
