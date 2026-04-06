# Python Unit Tests (`tools/unit_test`)

This directory contains pytest-based unit tests for Python scripts in `tools/`.

## Test Files and Scope

- `test_cluster_kon.py` → `cluster_kon.py`
    - NAM rate calculation, SDA input parsing, cluster beta routines
- `test_combine.py` → `combine.py`
    - header parsing, complex parsing, duplicate removal, integration flow
- `test_plot_angle.py` → `plot_clust_array_angle.py`
    - angle data parsing, cluster parsing, encounter indexing
- `test_plot_dist.py` → `plot_clust_array_dist.py`
    - distance data parsing, cluster parsing, encounter indexing
- `test_extract_grid.py` → `extract_grid.py`
    - UHBD header parsing and related grid parsing behavior
- `test_extract_C1_C4_filtered.py` → `extract_C1_C4_filtered.py`
    - complex loading and filtering helpers
- `test_filter_neg_enc.py` → `filter_neg_enc.py`
    - atom object behavior and related logic checks
- `test_gen_pdbs.py` → `gen_pdbs.py`
    - atom/data-structure level checks and transformation-related behavior
- `test_regioselectivity.py` → `regioselectivity.py`
    - CLI/main-path behavior checks
- `test_upper_C1_C4.py` → `upper_C1_C4.py`
    - atom-record parsing and z-coordinate helpers
- `test_utilities.py`
    - shared file-handling and parsing utility-style validations

## Test Infrastructure

- `pytest.ini`: test discovery, markers, output settings
- `conftest.py`: shared fixtures and temporary test data
- `run_tests.sh`: full-featured runner wrapper
- `test.sh`: quick runner wrapper

## Running Tests

From `tools/unit_test/`:

```bash
./run_tests.sh
```

Quick wrapper:

```bash
./test.sh
```

Examples:

```bash
./run_tests.sh --combine
./run_tests.sh --cluster -v
./run_tests.sh -c
```

## Notes

- Tests are designed to be self-contained and use temporary files/fixtures.
- Add new tests as `test_<script_name>.py` when new scripts/functions are introduced.

_Last updated: April 5, 2026_
