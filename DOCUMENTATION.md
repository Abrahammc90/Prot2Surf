"""
Clustering Program Documentation

Complete documentation for the cellulose clustering analysis pipeline.

## Overview

The clustering program is part of a comprehensive pipeline for analyzing
protein-cellulose interactions through molecular docking simulations.
It includes Fortran programs for high-performance clustering calculations
and Python tools for data analysis and visualization.

## Directory Structure

```
clustering_program/
├── bin/                            # Compiled executable binaries
│   ├── make_matrix                # Generate encounter complex matrices
│   ├── clust                      # Cluster encounter complexes
│   ├── analyze_residues           # Analyze residue-level interactions
│   └── threshold                  # Apply threshold filtering
│
├── src/                            # Source code and build files
│   ├── Makefile                   # Build configuration and test targets
│   ├── Doxyfile                   # Documentation generation config
│   ├── *.f90                      # Fortran source files
│   ├── *.mod                      # Compiled Fortran modules
│   ├── *.o                        # Object files
│   │
│   └── unit_test/                 # Fortran unit test suite
│       ├── test_*.f90             # 7 Fortran test programs
│       ├── run_all_tests.sh       # Fortran test runner
│       ├── README.md              # Test documentation
│       ├── TESTING.md             # Complete testing guide
│       └── INDEX.md               # Test index
│
├── tools/                          # Python auxiliary analysis tools
│   ├── README.md                  # Tools documentation
│   ├── cluster_kon.py             # NAM algorithm (KON calculations)
│   ├── combine.py                 # Merge complex files
│   ├── regioselectivity.py        # Binding site classification
│   ├── plot_clust_array_*.py      # Distribution visualization (2 scripts)
│   ├── extract_C1_C4_filtered.py  # Data extraction utilities (4 scripts)
│   ├── filter_neg_enc.py          # Complex filtering (1 script)
│   ├── gen_pdbs.py                # PDB generation
│   ├── upper_C1_C4.py             # Extract key binding atoms
│   └── unit_test/                 # Python unit tests
│       ├── test_*.py              # 5 Python test modules (120+ tests)
│       ├── conftest.py            # Pytest configuration
│       ├── pytest.ini             # Pytest settings
│       ├── run_tests.sh           # Comprehensive test runner
│       ├── test.sh                # Quick test shortcuts
│       ├── README.md              # Python test documentation
│       └── TESTING.md             # Testing guide
│
├── docs/                           # Documentation and figures
│   ├── doxygen/                   # Auto-generated API documentation
│   ├── latex/                     # LaTeX documentation sources
│   └── [other documentation]
│
├── Doxyfile                       # Main Doxygen configuration
├── Makefile                       # [Top-level build coordination]
└── README                         # [Quick start guide]
```

## Quick Start

### Building the Fortran Programs

```bash
cd src
make clean        # Remove compiled files
make all          # Build all executables and documentation
```

### Running Tests

```bash
# Run Fortran tests (in src/)
cd src && make test
cd src && make fortran_tests

# Run individual Fortran test
cd src && make test_maths

# Run Python tests (in tools/unit_test/)
cd tools/unit_test && ./run_tests.sh
cd tools/unit_test && ./test.sh quick
cd tools/unit_test && ./test.sh verbose
cd tools/unit_test && ./test.sh coverage
```

## Documentation

### Generating Documentation

```bash
# Generate all documentation
cd src
make docs

# View HTML documentation
open ../docs/doxygen/html/index.html
```

### Documentation Includes

- Fortran source code documentation
- Python module API documentation
- Unit test documentation
- Architecture and design guides
- Tool usage examples

## Fortran Programs

### make_matrix
Generates encounter complex distance matrices from PDB files and association data.

**Input:**
- Receptor PDB file
- Ligand PDB file
- Association file with complex configurations

**Output:**
- Distance/contact matrix file

**Compilation:** `make make_matrix`

### clust
Performs hierarchical clustering on encounter complexes.

**Input:**
- Distance matrix
- Cluster parameters

**Output:**
- Cluster assignments
- Cluster statistics

**Compilation:** `make clust`

### threshold
Filters encounter complexes based on distance and energy thresholds.

**Input:**
- Association file
- Distance/energy thresholds

**Output:**
- Filtered association file

**Compilation:** `make threshold`

### analyze_residues
Analyzes residue-level contributions to binding.

**Input:**
- Encounter complex configurations
- Residue definitions

**Output:**
- Residue interaction statistics

**Compilation:** `make analyze_residues`

## Python Tools

See [tools/README.md](tools/README.md) for detailed tool documentation.

### Key Tools

1. **cluster_kon.py** - Calculate reaction rate constants (NAM algorithm)
2. **combine.py** - Merge multiple encounter complex files
3. **regioselectivity.py** - Classify binding site preferences
4. **plot_clust_array_*.py** - Generate distribution plots (angle/distance)
5. **gen_pdbs.py** - Generate PDB files from encounter complexes

## Unit Testing

### Test Structure

- **Fortran Tests:** 7 compiled test programs covering core algorithms
- **Python Tests:** 5 test modules with 50+ test cases
- **Total Coverage:** ~200 test cases for comprehensive validation

### Running Tests

```bash
# Run Fortran tests
cd src && make test
cd src && make fortran_tests

# Run Python tests
cd tools/unit_test && ./run_tests.sh
cd tools/unit_test && ./test.sh quick
cd tools/unit_test && ./test.sh coverage
```

### Test Documentation

- [Fortran Tests](src/unit_test/) - Comprehensive Fortran unit tests
- [Python Tests](tools/unit_test/) - Complete Python test suite

## Build System

### Makefile Targets (src/)

```bash
# Core targets
make all              # Build everything and generate docs
make clean            # Remove all compiled files
make docs             # Generate Doxygen documentation

# Test targets (Fortran)
make test             # Run all Fortran tests
make fortran_tests    # Run Fortran tests explicitly
make test_maths       # Run specific Fortran test

# Individual program targets
make make_matrix      # Build make_matrix executable
make clust            # Build clust executable
make threshold        # Build threshold executable
make analyze_residues # Build analyze_residues executable
```

### Python Test Runners (tools/unit_test/)

```bash
./run_tests.sh        # Full test suite with all options
./test.sh quick       # Quick test run
./test.sh verbose     # Verbose output
./test.sh coverage    # Coverage report
```

### Compiler Configuration

- **Compiler:** gfortran
- **Flags:** `-O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math`
- **Parallelization:** OpenMP support enabled
- **Math Libraries:** LAPACK and BLAS linked

## Dependencies

### Fortran Programs
- gfortran compiler
- LAPACK library
- BLAS library
- OpenMP support

### Python Tools
- Python 3.7+
- numpy (numerical operations)
- matplotlib (visualization)

### Testing
- pytest >= 6.0
- pytest-cov (coverage reports)
- pytest-xdist (parallel execution)

## Performance Optimization

- OpenMP parallelization for matrix operations
- LAPACK/BLAS for linear algebra
- Memory-efficient algorithms for large datasets
- Native architecture optimizations

## Architecture

### Module Dependencies

```
read_input.f90
├── mod_pdb.f90
└── mod_assoc.f90

mod_matrix.f90
└── maths.f90

mod_threshold.f90
└── maths.f90

make_matrix.f90
├── mod_pdb.f90
├── mod_assoc.f90
├── mod_matrix.f90
└── read_input.f90

clust.f90
├── mod_clust_algorithm.f90
├── mod_matrix.f90
└── mod_pdb.f90

threshold.f90
├── mod_threshold.f90
├── mod_pdb.f90
├── mod_assoc.f90
└── maths.f90
```

## Continuous Integration

Tests can be integrated into CI/CD pipelines:

```bash
cd src
make test
```

Exit code 0 indicates all tests passed.

## Contributing

When adding new features:

1. Update corresponding source files
2. Add unit tests (Fortran or Python as appropriate)
3. Update documentation
4. Run full test suite: `make test`
5. Generate documentation: `make docs`

## Authors

- Abraham Muñiz-Chicharro: Python tools and testing infrastructure
- Original Fortran implementation team

## Version

Clustering Program v1.0 - February 2026

## License

See LICENSE file in project root.

## References

- SDA (Simulation of Diffusional Association) methodology
- NAM algorithm for rate constant calculations
- Hierarchical clustering techniques
"""
