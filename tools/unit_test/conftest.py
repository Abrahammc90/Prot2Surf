"""Pytest configuration and shared fixtures for unit tests.

@version 1.0
"""

import os
import tempfile
import warnings
import pytest

try:
    from pyparsing import PyparsingDeprecationWarning
except Exception:
    PyparsingDeprecationWarning = Warning

warnings.filterwarnings(
    "ignore",
    category=PyparsingDeprecationWarning,
    module=r"matplotlib\._fontconfig_pattern",
)
warnings.filterwarnings(
    "ignore",
    category=PyparsingDeprecationWarning,
    module=r"matplotlib\._mathtext",
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_sda_input_file(temp_dir):
    """Create a sample SDA input file for testing."""
    filepath = os.path.join(temp_dir, "test_sda.in")
    content = """# SDA Input Configuration
start_pos = 50.0
c = 1.5
nrun = 1000
diffusion_trans = 0.05
diffusion_rot = 0.08
"""
    with open(filepath, 'w') as f:
        f.write(content)
    return filepath


@pytest.fixture
def sample_association_file(temp_dir):
    """Create a sample association/complex file with header and entries."""
    filepath = os.path.join(temp_dir, "test_assoc.txt")
    content = """# Header line 1: Total complexes
# Header line 2: Description
# Header line 3: Parameters
# Header line 4: Column info
  1.0  2.0  3.0  1  0.5
  1.5  2.5  3.5  2  0.6
  2.0  3.0  4.0  3  0.7
"""
    with open(filepath, 'w') as f:
        f.write(content)
    return filepath


@pytest.fixture
def sample_cluster_file(temp_dir):
    """Create a sample cluster assignment file."""
    filepath = os.path.join(temp_dir, "test_clusters.txt")
    content = """Cluster: 0
  1 2 3
Cluster: 1
  4 5 6
Cluster: 2
  7 8
"""
    with open(filepath, 'w') as f:
        f.write(content)
    return filepath


@pytest.fixture
def sample_pdb_file(temp_dir):
    """Create a sample PDB file for testing."""
    filepath = os.path.join(temp_dir, "test.pdb")
    content = """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
TER
END
"""
    with open(filepath, 'w') as f:
        f.write(content)
    return filepath


@pytest.fixture
def sample_distance_data_file(temp_dir):
    """Create a sample distance data file."""
    filepath = os.path.join(temp_dir, "test_distances.txt")
    content = """5.2 6.1 7.3 8.2 5.9
6.5 7.1 8.0 6.8 7.2
"""
    with open(filepath, 'w') as f:
        f.write(content)
    return filepath
