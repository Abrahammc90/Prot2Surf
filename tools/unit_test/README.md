"""README for unit tests of clustering analysis tools."""

# Unit Tests for Clustering Analysis Tools

This directory contains comprehensive unit tests for all Python scripts in the `tools/` directory.

## Test Structure

```
unit_test/
- __init__.py                 # Package initialization
- conftest.py                 # Pytest fixtures and configuration
- test_cluster_kon.py        # Tests for cluster_kon.py (NAM algorithm, SDA parsing)
- test_combine.py            # Tests for combine.py (file merging)
- test_plot_angle.py         # Tests for plot_clust_array_angle.py (angle visualization)
- test_plot_dist.py          # Tests for plot_clust_array_dist.py (distance visualization)
- test_utilities.py          # Tests for utility functions (PDB, file handling)
- README.md                  # This file
```

## Test Coverage

### test_cluster_kon.py
Tests the NAM (Northrup-Allison-McCammon) algorithm and SDA input parsing:
- `TestNAMAlgorithm`: Rate constant calculations with various parameters
- `TestSDAInputParser`: SDA configuration file parsing with comments and whitespace
- `TestClusterBeta`: Cluster fraction calculations

**Classes Tested:** NAM_algorithm, sda_input_parser, cluster_beta

### test_combine.py
Tests encounter complex file merging functionality:
- `TestLoadHeaderFromFile`: Header extraction from association files
- `TestLoadComplexesFromFile`: Complex entry loading (skipping headers)
- `TestCombineComplexesRemovingDuplicates`: Duplicate removal and merging
- `TestIntegration`: Full workflow integration

**Functions Tested:** load_header_from_file, load_complexes_from_file, combine_complexes_removing_duplicates

### test_plot_angle.py
Tests angle distribution visualization and data parsing:
- `TestGenerateEncounterXValues`: Sequential x-axis value generation
- `TestParseAngleData`: Angle data file parsing with various formats
- `TestParseClusterFile`: Cluster assignment file parsing
- `TestIntegration`: Complete plotting workflow

**Functions Tested:** _generate_encounter_x_values, _parse_angle_data, _parse_cluster_file

### test_plot_dist.py
Tests distance distribution visualization and data parsing:
- `TestGenerateEncounterXValues`: Sequential x-axis value generation
- `TestParseDistanceData`: Distance data file parsing
- `TestParseClusterFile`: Cluster assignment file parsing with distance mapping
- `TestIntegration`: Complete plotting workflow

**Functions Tested:** _generate_encounter_x_values, _parse_distance_data, _parse_cluster_file

### test_utilities.py
Tests common file handling and data processing:
- `TestPDBFileHandling`: PDB file parsing and coordinate extraction
- `TestAssociationFileHandling`: Association file format handling
- `TestClusterFileHandling`: Cluster file parsing and indexing
- `TestDataTypeConversions`: String to numeric conversions
- `TestFileIOErrors`: File I/O error handling
- `TestFileFormatConsistency`: Format preservation and whitespace handling

## Running Tests

### Install pytest
```bash
pip install pytest
```

### Run all tests
```bash
pytest unit_test/
```

### Run specific test file
```bash
pytest unit_test/test_cluster_kon.py
```

### Run specific test class
```bash
pytest unit_test/test_cluster_kon.py::TestNAMAlgorithm
```

### Run specific test function
```bash
pytest unit_test/test_cluster_kon.py::TestNAMAlgorithm::test_nam_algorithm_basic
```

### Run with verbose output
```bash
pytest unit_test/ -v
```

### Run with coverage report
```bash
pip install pytest-cov
pytest unit_test/ --cov=.. --cov-report=html
```

### Run with output
```bash
pytest unit_test/ -s
```

## Fixtures

Common fixtures provided by `conftest.py`:

- `temp_dir`: Temporary directory for test files
- `sample_sda_input_file`: Sample SDA input file with parameters
- `sample_association_file`: Sample association file with header and entries
- `sample_cluster_file`: Sample cluster assignment file
- `sample_pdb_file`: Sample PDB protein structure file
- `sample_distance_data_file`: Sample distance data file

## Test Design Principles

1. **Isolation**: Each test is independent and uses temporary files
2. **Clarity**: Test names clearly describe what is being tested
3. **Coverage**: Tests cover normal cases, edge cases, and error conditions
4. **Fixtures**: Reusable test data through pytest fixtures
5. **Integration**: Integration tests verify complete workflows

## Expected Behavior

### Passing Tests
All tests should pass by default when tools are implemented correctly.

### Test Organization
- Basic functionality tests
- Edge case tests
- Error handling tests
- Integration tests

## Continuous Integration

These tests can be integrated into CI/CD pipelines:

```bash
# Run tests and exit with failure code if any fail
pytest unit_test/ --tb=short
```

## Contributing

When adding new tools or modifying existing ones:

1. Add corresponding test file (test_toolname.py)
2. Create test classes for major functionality
3. Include tests for:
   - Basic functionality
   - Edge cases
   - Error handling
   - Data validation
4. Use descriptive test names
5. Run full test suite before committing

## Dependencies

- pytest >= 6.0
- numpy (for numerical tests)
- tempfile (standard library)

## Author

Test suite created for cellulose clustering analysis project.
Abraham Muñiz-Chicharro

## Version

Version 1.0 - February 2026
