# Auxiliary Tools for Clustering Analysis

This directory contains Python auxiliary tools for processing and analyzing encounter complexes from SDA (Simulation of Diffusional Association) simulations.

**Note:** Unit tests for these Python tools are located in `unit_test/`. See [TESTING.md](unit_test/TESTING.md) for test documentation.

## Overview

These tools are part of the cellulose clustering analysis pipeline and provide functionality for:
- Kinetic rate calculations
- Complex file manipulation and filtering
- Visualization of cluster distributions
- Regioselectivity analysis
- PDB file generation and processing

## Tools Description

### Core Analysis Tools

#### cluster_kon.py
**Purpose:** Calculate Smoluchowski rate constants using the NAM (Northrup-Allison-McCammon) algorithm.

**Features:**
- NAM algorithm implementation for diffusion-limited reactions
- SDA input file parsing
- Threshold energy evaluation
- Automatic rate constant calculation

**Usage:**
```bash
python cluster_kon.py <sda_input_file> [options]
```

**Dependencies:** numpy, sys

---

#### regioselectivity.py
**Purpose:** Classify encounter complexes by binding sites (C1/C4 carbon positions).

**Features:**
- Binding site classification
- Statistical analysis of regioselectivity
- Reference complex comparison
- Output files for C1 and C4 classifications

**Usage:**
```bash
python regioselectivity.py <complexes_file> <reference_C1> <reference_C4> [options]
```

**Dependencies:** argparse, sys

---

### File Processing Tools

#### combine.py
**Purpose:** Merge multiple encounter complex association files into a single file.

**Features:**
- Preserves 4-line SDA header format
- Handles multiple input files
- Automatic duplicate handling
- Memory-efficient processing

**Usage:**
```bash
python combine.py <output_file> <input_file1> <input_file2> [...]
```

**Dependencies:** sys

---

#### extract_C1_C4_filtered.py
**Purpose:** Filter encounter complexes based on reference configurations.

**Features:**
- Reference-based filtering
- Multiple reference support
- Maintains association file format
- Selective complex extraction

**Usage:**
```bash
python extract_C1_C4_filtered.py -complexes <file1> <file2> ... -references <ref1> <ref2> ...
```

**Dependencies:** argparse, sys

---

#### filter_neg_enc.py
**Purpose:** Remove encounter complexes with negative z-coordinates.

**Features:**
- PDB file parsing
- Geometric filtering
- Mass-weighted calculations
- Multi-chain support

**Usage:**
```bash
python filter_neg_enc.py <input_complexes> [options]
```

**Dependencies:** numpy, argparse

---

### Visualization Tools

#### plot_clust_array_angle.py
**Purpose:** Generate angle distribution plots for encounter complex clusters.

**Features:**
- Multi-cluster visualization
- Color-coded cluster lines
- High-resolution output (300 DPI)
- Publication-quality figures

**Usage:**
```bash
python plot_clust_array_angle.py <input_data> [options]
```

**Dependencies:** matplotlib, sys

---

#### plot_clust_array_dist.py
**Purpose:** Generate distance distribution plots for encounter complex clusters.

**Features:**
- Multi-cluster visualization
- Color-coded cluster lines
- High-resolution output (300 DPI)
- Publication-quality figures

**Usage:**
```bash
python plot_clust_array_dist.py <input_data> [options]
```

**Dependencies:** matplotlib, sys

---

### PDB Generation and Processing Tools

#### gen_pdbs.py
**Purpose:** Generate individual PDB files from encounter complex configurations.

**Features:**
- Coordinate transformation (rotation + translation)
- Mass-weighted center of mass calculations
- Selective residue filtering
- Multi-chain protein support

**Usage:**
```bash
python gen_pdbs.py <receptor_pdb> <ligand_pdb> <assoc_file> [options]
```

**Options:**
- `-n, --num`: Number of complexes to process (default: all)
- `--resnames`: Filter by residue names
- `-o, --output`: Output file prefix

**Dependencies:** numpy, argparse

---

#### extract_grid.py
**Purpose:** Extract electrostatic potential grid data from UHBD format files.

**Features:**
- Binary UHBD format parsing
- Grid dimension extraction
- 3D numpy array reshaping
- UTF-8 binary data handling

**Usage:**
```bash
python extract_grid.py <input.uhbd> <Nx> <Ny> <Nz> <output.uhbd>
```

**Dependencies:** numpy, sys

---

#### upper_C1_C4.py
**Purpose:** Extract uppermost C1 and C4 atoms from cellulose structures.

**Features:**
- Maximum z-coordinate identification (0.1 Å tolerance)
- Selective C1/C4 atom filtering
- PDB format preservation
- Terminal record handling

**Usage:**
```bash
python upper_C1_C4.py <input_pdb_file>
```

**Output:** upper_C1_C4.pdb

**Dependencies:** sys

---

## Integration with Main Pipeline

These tools are called by the main clustering analysis script (`03_run_all_cluster_analysis.sh`):

1. **combine.py** - Merges aligned encounter complexes
2. **regioselectivity.py** - Analyzes binding site preferences
3. **plot_clust_array_*.py** - Visualizes cluster distributions

## Common Data Formats

### Association File Format
```
# Line 1: Number of complexes
# Line 2: Header information
# Line 3: Column descriptions
# Line 4: Additional metadata
<complex_data_lines>
```

### PDB File Format
Standard Protein Data Bank format with ATOM/HETATM records.

## Author

Abraham Muñiz-Chicharro

## Version

All tools: Version 1.0

## Documentation

Full API documentation is available via Doxygen. Generate documentation by running:
```bash
cd /home/abraham/cellulose_project
doxygen Doxyfile
```

The generated HTML documentation will be in `docs/doxygen/html/`.
