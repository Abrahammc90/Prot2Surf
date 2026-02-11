"""Unit tests for cluster_kon.py - KON rate constant calculation."""

import sys
import os
import pytest
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cluster_kon import NAM_algorithm, sda_input_parser, cluster_beta


class TestNAMAlgorithm:
    """Tests for NAM_algorithm function."""

    def test_nam_algorithm_basic(self):
        """Test NAM algorithm with basic parameters."""
        b_surface = 10.0  # Ų
        q_surface = 15.0  # Ų
        D = 0.05  # Ų/ps
        beta = 0.4  # 1/(kcal/mol)
        
        kon = NAM_algorithm(b_surface, q_surface, D, beta)
        
        assert isinstance(kon, (float, np.floating))
        assert kon > 0  # Rate constant should be positive
        assert np.isfinite(kon)  # Should be a finite number

    def test_nam_algorithm_physical_bounds(self):
        """Test NAM algorithm returns physically reasonable values."""
        # Test with typical cellulose-protein values
        kon1 = NAM_algorithm(8.0, 12.0, 0.05, 0.4)
        kon2 = NAM_algorithm(8.0, 12.0, 0.10, 0.4)
        
        # Higher diffusion should give higher rate constant
        assert kon2 > kon1

    def test_nam_algorithm_edge_cases(self):
        """Test NAM algorithm with edge cases."""
        # Test with very small values
        kon_small = NAM_algorithm(1.0, 1.5, 0.01, 0.1)
        assert kon_small > 0
        
        # Test with larger values
        kon_large = NAM_algorithm(20.0, 30.0, 0.1, 1.0)
        assert kon_large > 0

    def test_nam_algorithm_symmetry(self):
        """Test NAM algorithm with different but proportional surfaces."""
        kon1 = NAM_algorithm(10.0, 15.0, 0.05, 0.4)
        kon2 = NAM_algorithm(20.0, 30.0, 0.05, 0.4)
        
        # Rate should scale with surface area
        ratio = kon2 / kon1
        assert 1.9 < ratio < 2.1  # Should be approximately 2.0


class TestSDAInputParser:
    """Tests for sda_input_parser function."""

    def test_parser_valid_file(self, sample_sda_input_file):
        """Test parser with valid SDA input file."""
        params = sda_input_parser(sample_sda_input_file)
        
        assert 'start_pos' in params
        assert 'c' in params
        assert 'nrun' in params
        assert 'D' in params
        
        assert params['start_pos'] == 50.0
        assert params['c'] == 1.5
        assert params['nrun'] == 1000
        assert params['D'] == 0.05

    def test_parser_skips_comments(self, temp_dir):
        """Test that parser correctly skips comment lines."""
        filepath = os.path.join(temp_dir, "test_comments.in")
        content = """# This is a comment
start_pos = 100.0
# Another comment
c = 2.0
nrun = 2000
diffusion_trans = 0.1
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        params = sda_input_parser(filepath)
        assert params['start_pos'] == 100.0
        assert params['c'] == 2.0

    def test_parser_handles_whitespace(self, temp_dir):
        """Test parser handles various whitespace configurations."""
        filepath = os.path.join(temp_dir, "test_whitespace.in")
        content = """start_pos   =   75.0
c=1.8
  nrun  =  1500  
diffusion_trans    =    0.06
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        params = sda_input_parser(filepath)
        assert params['start_pos'] == 75.0
        assert params['c'] == 1.8
        assert params['nrun'] == 1500

    def test_parser_multiple_diffusion_terms(self, temp_dir):
        """Test parser accumulates multiple diffusion contributions."""
        filepath = os.path.join(temp_dir, "test_diffusion.in")
        content = """start_pos = 50.0
c = 1.5
nrun = 1000
diffusion_trans = 0.03
diffusion_trans = 0.02
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        params = sda_input_parser(filepath)
        assert params['D'] == pytest.approx(0.05)

    def test_parser_file_not_found(self):
        """Test parser raises error for non-existent file."""
        with pytest.raises(FileNotFoundError):
            sda_input_parser("nonexistent_file.in")


class TestClusterBeta:
    """Tests for cluster_beta function."""

    def test_cluster_beta_basic(self, sample_cluster_file):
        """Test cluster_beta calculation with sample file."""
        beta = cluster_beta(sample_cluster_file, 8)
        
        assert 0.0 <= beta <= 1.0
        assert isinstance(beta, (float, np.floating))

    def test_cluster_beta_all_clustered(self, temp_dir):
        """Test when all trajectories are in clusters."""
        filepath = os.path.join(temp_dir, "test_all_clustered.txt")
        content = """Cluster: 0
  1 2 3
Cluster: 1
  4 5 6
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        beta = cluster_beta(filepath, 6)
        # All 6 trajectories should be accounted for
        assert beta == pytest.approx(1.0)

    def test_cluster_beta_partial_clustering(self, temp_dir):
        """Test when only some trajectories are clustered."""
        filepath = os.path.join(temp_dir, "test_partial.txt")
        content = """Cluster: 0
  1 2
Cluster: 1
  3 4
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        beta = cluster_beta(filepath, 10)
        # 4 out of 10 trajectories in clusters
        assert 0.0 < beta < 1.0

    def test_cluster_beta_no_clusters(self, temp_dir):
        """Test with empty cluster file."""
        filepath = os.path.join(temp_dir, "test_no_clusters.txt")
        content = ""
        with open(filepath, 'w') as f:
            f.write(content)
        
        beta = cluster_beta(filepath, 100)
        assert beta == 0.0
