# Clustering Program - Unit Tests

## Organization

The unit tests have been organized as follows:

### Fortran Unit Tests
**Location:** `src/unit_test/`
- Contains 7 Fortran test programs:
  - `test_maths.f90` - Mathematical operations (cross product, coordinate transformation, etc.)
  - `test_matrix.f90` - Matrix operations
  - `test_threshold.f90` - Threshold calculations
  - `test_pdb.f90` - PDB file handling
  - `test_clustering.f90` - Clustering algorithm
  - `test_assoc.f90` - Association file reading
  - `test_read_input.f90` - Input file parsing
- Documentation files:
  - `README.md` - Test overview
  - `TESTING.md` - Detailed testing guide
  - `INDEX.md` - Test index
  - `TEST_SUMMARY.md` - Summary of test coverage
  - `run_all_tests.sh` - Test runner script

### Python Unit Tests
**Location:** `tools/unit_test/`
- Contains 5 Python test modules:
  - `test_cluster_kon.py` - Tests for NAM algorithm and cluster analysis
  - `test_combine.py` - Tests for file merging
  - `test_plot_angle.py` - Tests for angle plotting
  - `test_plot_dist.py` - Tests for distance plotting
  - `test_utilities.py` - Tests for utility functions
- Configuration files:
  - `conftest.py` - pytest fixtures
  - `pytest.ini` - pytest configuration
  - `README.md` - Python test documentation
- Test runners:
  - `run_tests.sh` - Comprehensive test runner with multiple options
  - `test.sh` - Quick test shortcuts

## Build System Integration

### Makefile (src/)
- `make test` - Runs all Fortran unit tests
- `make fortran_tests` - Runs Fortran tests explicitly
- `make test_maths` - Runs specific Fortran test

Note: Python tests have their own test runners in `tools/unit_test/`

### Doxygen Configuration
- **Root Doxyfile** (for Python and auxiliary tools)
  - FILE_PATTERNS includes: *.py *.sh *.md
  - EXTENSION_MAPPING: py=Python, sh=C
  
- **src/Doxyfile** (for Fortran and unit tests)
  - FILE_PATTERNS includes: *.f90 *.F90 *.f *.F *.py *.sh *.md
  - INPUT = . (includes unit_test directory)
  - RECURSIVE = YES (discovers all test files)

## Summary

✓ **Fortran tests** are located in `src/unit_test/`
✓ **Python tests** are located in `tools/unit_test/`
✓ **Makefile** points to Fortran tests in src
✓ **Doxygen** configured to discover all test files through recursive INPUT

This structure ensures:
- Clear separation of unit tests (Fortran tests in src, Python tests in tools)
- Proper build system coordination
- Complete documentation generation
- Easy test execution from appropriate directories
