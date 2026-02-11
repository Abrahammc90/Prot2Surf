# Cellulose Clustering Program

## Overview

The Cellulose Clustering Program is a comprehensive Fortran-based tool designed for analyzing and clustering molecular structures, particularly focused on cellulose proteins and their interactions. The program provides powerful algorithms for hierarchical clustering, residue analysis, and distance matrix generation.

## Main Components

### Executables

- **make_matrix**: Generates distance matrices from molecular structures
- **clust**: Performs hierarchical clustering analysis on distance matrices
- **analyze_residues**: Analyzes individual residues and their properties
- **threshold**: Applies thresholding operations for clustering results

### Key Modules

- `mod_pdb`: Handles PDB file reading and structure manipulation
- `mod_matrix`: Matrix operations and management
- `mod_clust_algorithm`: Core clustering algorithms implementation
- `mod_threshold`: Thresholding operations and utilities
- `mod_assoc`: Association and relationship tracking
- `maths`: Mathematical utilities and linear algebra operations

## Features

- High-performance Fortran implementation with OpenMP support
- Linear algebra operations via LAPACK and BLAS
- Efficient clustering algorithms for large molecular datasets
- Residue-level analysis capabilities
- Distance matrix generation and manipulation
- Thresholding and filtering operations

## Compilation

To compile the program:

```bash
cd src
make
```

This will:
1. Generate all executable binaries in the `bin/` directory
2. Create HTML and PDF documentation in the `docs/` directory

## Usage

For detailed usage information, see the individual module and executable documentation.

## Author & License

Part of the Cellulose Project analysis suite.
