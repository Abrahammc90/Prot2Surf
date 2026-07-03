# Clustering Program (Fortran + Python Tools)

This project provides a complete encounter-analysis and clustering workflow for SDA-style molecular simulations.

It includes:

- High-performance Fortran executables for array clustering and encounter filtering.
- Python utilities for preprocessing, filtering, plotting, and post-analysis.
- Fortran and Python unit test suites.

---

## Repository Layout

- `src/`: Fortran source, Makefile, Fortran test code
- `bin/`: compiled executables
- `tools/`: auxiliary python scripts for data-processing
- `tools/unit_test/`: Python tests (pytest)
- `docs/`: generated Doxygen outputs

---

## Fortran Executables

Built from `src/Makefile` into `bin/`:

- `make_data`: data preparation helpers
- `clust`: hierarchical clustering from array input
- `threshold`: threshold-based encounter selection
- `contact_map`: residue proximity analysis across encounters

Program-specific documentation:

- `src/make_data.md`
- `src/clust.md`
- `src/threshold.md`
- `src/contact_map.md`

---

## Core Fortran Modules

- `mod_clust_algorithm.f90`: array-based mean-linkage hierarchical clustering, OpenMP (CPU)
- `mod_array.f90`: array construction/read/write and clustering helpers
- `mod_threshold.f90`: threshold arrays and sorting utilities
- `mod_pdb.f90`: PDB structures and parsing utilities
- `mod_assoc.f90`: association/complexes file structures and parsing
- `read_input.f90`: high-level loaders for PDB and association inputs
- `maths.f90`: numerical and geometric helper routines

---

## Clustering Execution Model

Clustering is CPU-only in the Fortran implementation:

- array-based clustering runs in `mod_clust_algorithm.f90`
- OpenMP parallel sections are used where implemented

---

## Python Tools (`tools/`)

- `combine.py`: merges association files and removes duplicate complexes
- `extract_C1_C4_filtered.py`: filters complexes by C1/C4 references
- `extract_grid.py`: parses and extracts UHBD grid data
- `gen_pdbs.py`: generates transformed PDB structures from complexes
- `plot_closest_residues.py`: plots closest-residue encounter counts from residue reports
- `plot_clust_array_angle.py`: plots angle distributions per cluster
- `plot_clust_array_dist.py`: plots distance distributions per cluster
- `regioselectivity.py`: classification and percentages by reference matches
- `upper_C1_C4.py`: extracts uppermost C1/C4 atoms from cellulose structures

---

## Testing

### Fortran Tests (`src/unit_test/`)

- `test_maths.f90`
- `test_array.f90`
- `test_threshold.f90`
- `test_pdb.f90`
- `test_clustering.f90`
- `test_assoc.f90`
- `test_read_input.f90`

From `src/`:

```bash
make test
```

For standalone runner:

```bash
cd unit_test
bash run_all_tests.sh
```

### Python Tests (`tools/unit_test/`)

Pytest-based suite covering all tools scripts:

- `test_combine.py`
- `test_extract_C1_C4_filtered.py`
- `test_extract_grid.py`
- `test_gen_pdbs.py`
- `test_plot_closest_residues.py`
- `test_plot_angle.py`
- `test_plot_dist.py`
- `test_regioselectivity.py`
- `test_upper_C1_C4.py`
- `test_utilities.py`

From `tools/unit_test/`:

```bash
./run_tests.sh
```

See `tools/unit_test/README.md` for details.

---

## Build and Documentation

### Documentation Dependencies (Linux/Debian/Ubuntu)

To generate Doxygen HTML and LaTeX/PDF docs (`make docs`), install:

```bash
sudo apt-get update
sudo apt-get install -y doxygen graphviz make \
	texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```

Notes:

- `doxygen`: parses source and generates documentation outputs.
- `graphviz` (`dot`): generates dependency/call graphs used by Doxygen.
- `texlive-latex-extra`: provides `float.sty` (required by Doxygen's `doxygen.sty`).
- `texlive-latex-base` and `texlive-fonts-recommended`: core LaTeX packages for PDF build.

From `src/`:

```bash
make        # build binaries
make docs   # generate Doxygen + PDF docs
make clean  # remove binaries/objects/tests/docs outputs
```

---

## Notes

- CLI usage/help for each Fortran executable is implemented in program `print_help` routines.
- Python tools require Python 3 and scientific dependencies used by each script (`numpy`, `matplotlib`, etc.).

---

## Author

Abraham Muñiz-Chicharro
