# Unit Testing Suite - File Index

## Test Source Files

### Core Test Programs (7 files) - UPDATED

1. **test_maths.f90** (29.4 KB) - EXPANDED
   - Tests: 31 tests (updated from 4)
   - Modules tested: `maths`
   - Features: Cross product, coordinate transformations, angle calculations, RMSD, center of geometry, distance calculations
   - Dependencies: maths.f90, mod_pdb.f90
   - Lines: 850

2. **test_matrix.f90** (13.1 KB)
   - Tests: 8 tests
   - Modules tested: `mod_matrix`
   - Features: Matrix I/O, array I/O, Z-coordinate calculations, RMSD matrix properties
   - Dependencies: mod_matrix.f90, maths.f90, mod_pdb.f90
   - Lines: 370

3. **test_threshold.f90** (13.3 KB) - UPDATED
   - Tests: 9 tests (updated from 6)
   - Modules tested: `mod_threshold`
   - Features: Z-coordinate arrays, minimum distance, angle calculations, array sorting
   - Dependencies: mod_threshold.f90, maths.f90
   - Lines: 375

4. **test_pdb.f90** (7.4 KB)
   - Tests: 8 tests
   - Modules tested: `mod_pdb`
   - Features: PDB structures, atoms, residues, chains
   - Dependencies: mod_pdb.f90
   - Lines: 210

5. **test_clustering.f90** (9.8 KB)
   - Tests: 5 tests
   - Modules tested: `mod_clust_algorithm`
   - Features: Clustering statistics, simple clustering, cluster output writing
   - Dependencies: mod_clust_algorithm.f90, read_input.f90, mod_pdb.f90, mod_assoc.f90, maths.f90
   - Lines: 275

6. **test_assoc.f90** (13.5 KB) ✓ NEW
   - Tests: 11 tests
   - Modules tested: `mod_assoc`
   - Features: Association file reading, object allocation, file parsing
   - Dependencies: mod_assoc.f90
   - Lines: 380

7. **test_read_input.f90** (12.6 KB) ✓ NEW
   - Tests: 10 tests
   - Modules tested: `read_input` (wrapper module)
   - Features: Input file reading, structure creation, error handling
   - Dependencies: read_input.f90, mod_pdb.f90, mod_assoc.f90
   - Lines: 355

## Documentation Files

1. **TESTING.md** (7.2 KB) - UPDATED
   - Comprehensive testing guide with all 7 test files
   - Compilation instructions for each test
   - Updated test coverage table (82 tests total)
   - Best practices and troubleshooting
   - Instructions for adding new tests

2. **TEST_SUMMARY.md** (6.8 KB) - UPDATED
   - Overview of entire test suite with 82 tests
   - Summary of all 7 test programs
   - Test coverage matrix with new modules
   - Example output
   - Overall statistics table

3. **QUICK_TEST.md** (3.2 KB) - UPDATED
   - Quick reference for running all 7 test suites
   - Command summary table with new tests
   - File overview with test counts
   - Pro tips

4. **README.md** (3.5 KB) - UPDATED
   - Test framework overview
   - Directory structure with all 7 test programs
   - Quick start guide
   - Test coverage table (82 tests)
   - Build requirements and verification

5. **INDEX.md** (This file) - UPDATED
   - File listing and descriptions for all 7 test programs
   - Documentation index with 5 files

## Automation Files

1. **run_all_tests.sh** (2.8 KB) - UPDATED
   - Bash script to compile and run all 7 tests
   - Automatic dependency tracking for new tests
   - Summary reporting for all test programs
   - Exit code handling (0 = success, 1 = failure)

## Modified Files

1. **Makefile** - UPDATED
   - Added test targets for test_assoc and test_read_input
   - Updated `test` target to run all 7 tests
   - Updated `clean_tests` target for new executables
   - Integrated with existing build system

## Statistics

| Metric | Value | Updated |
|--------|-------|---------|
| Total test files | 7 | ↑ from 5 |
| Total tests | 82 | ↑ from 28 |
| Documentation files | 5 | ↑ from 4 |
| Total source lines | 2,815 | ↑ from ~1,500 |
| Total size | ~120 KB | ↑ from ~60 KB |
| Modules covered | 7 | ↑ from 5 |
| Functions/subroutines tested | 26 | ↑ from 11 |

## Quick Start

```bash
cd clustering_program/src/

# Run all tests
make test

# Or use the automated script
bash run_all_tests.sh

# Individual test
make test_maths
make test_pdb
```

## Test Distribution

```
Mathematical Operations (4 tests)
├─ cross product
└─ coordinate transformations

Matrix Operations (4 tests)
├─ allocation
└─ symmetry properties

Array/Threshold Operations (6 tests)
├─ statistics
└─ bounds checking

PDB Data Structures (8 tests)
├─ file allocation
├─ atom properties
└─ residue management

Clustering Algorithms (6 tests)
├─ distance matrix properties
└─ clustering helpers
```

## Testing Workflow

```
1. Build Tests (make test)
   ↓
2. Compile each test program
   ├─ Dependency resolution
   └─ Link with required modules
   ↓
3. Execute tests sequentially
   ├─ Each test reports pass/fail
   ├─ Detailed error messages on failure
   └─ Summary statistics
   ↓
4. Report results
   ├─ Total tests run
   ├─ Passed count
   ├─ Failed count
   └─ Exit code (0=success, 1=failure)
```

## Next Steps

1. **Run tests**: `make test`
2. **Read results**: Check PASS/FAIL output
3. **Add tests**: Follow instructions in TESTING.md
4. **Integrate CI**: Set up GitHub Actions (optional)
5. **Performance**: Run with profiling tools if needed

## Files Generated by Tests

When tests run, the following executables are created (removed by `make clean_tests`):

- `test_maths`
- `test_matrix`
- `test_threshold`
- `test_pdb`
- `test_clustering`

## Support & References

- **TESTING.md** - Detailed guide with troubleshooting
- **TEST_SUMMARY.md** - Full overview and future enhancements
- **QUICK_TEST.md** - Common commands and tips

---

**Created**: February 2026
**Total Tests**: 28
**Coverage**: 5 modules, 11 functions
