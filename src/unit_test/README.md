# Fortran Unit Tests (`src/unit_test`)

This directory contains standalone Fortran test programs for core clustering modules.

## Included Test Programs

- `test_maths.f90`: numerical/geometry operations in `maths`
- `test_array.f90`: array I/O, array-based encounter metrics, and array/metadata reordering in `mod_array`
- `test_threshold.f90`: threshold operations in `mod_threshold`
- `test_pdb.f90`: PDB parsing/structures in `mod_pdb`
- `test_clustering.f90`: clustering logic in `mod_clust_algorithm`
- `test_assoc.f90`: association file handling in `mod_assoc`
- `test_read_input.f90`: integrated input loading in `read_input`

## Current Coverage Notes

- `test_array.f90` currently covers:
	- `write_array()` / `read_array()`
	- `array_z_coord()`
	- `array_atoms_dist()`
	- `array_angle()`
	- `merge_sorted_segments()`
- `test_threshold.f90` covers threshold-oriented array calculations in `mod_threshold`.

## Default Test Execution Paths

### From `src/` (Makefile target)

```bash
make test
```

This is the authoritative default path. It runs the standard Fortran test set defined in `src/Makefile`, including `test_array`.

### From `src/unit_test/` (script)

```bash
bash run_all_tests.sh
```

This script compiles/runs the list currently defined in `run_all_tests.sh`.

## Author

Abraham Muñiz-Chicharro
