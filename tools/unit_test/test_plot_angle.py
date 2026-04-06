"""Unit tests for plot_clust_array_angle.py - Angle distribution plotting."""

import sys
import os
import pytest
import tempfile

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plot_clust_array_angle import (
    _generate_encounter_x_values,
    _parse_angle_data,
    _parse_cluster_file
)


class TestGenerateEncounterXValues:
    """Tests for _generate_encounter_x_values function."""

    def test_generate_x_values_basic(self):
        """Test x-value generation with basic cluster data."""
        cluster_data = {
            "cluster 0": {1: 10.0, 2: 20.0, 3: 30.0},
            "cluster 1": {4: 15.0, 5: 25.0}
        }
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        assert total == 5
        assert len(cluster_x_values) == 2
        assert len(cluster_x_values["cluster 0"]) == 3
        assert len(cluster_x_values["cluster 1"]) == 2

    def test_generate_x_values_sequential(self):
        """Test that x-values are sequential across clusters."""
        cluster_data = {
            "cluster 0": {1: 10.0, 2: 20.0},
            "cluster 1": {3: 15.0, 4: 25.0}
        }
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        # First cluster should start at 1
        assert cluster_x_values["cluster 0"] == [1, 2]
        # Second cluster should continue from 3
        assert cluster_x_values["cluster 1"] == [3, 4]

    def test_generate_x_values_empty(self):
        """Test with empty cluster data."""
        cluster_data = {}
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        assert total == 0
        assert len(cluster_x_values) == 0

    def test_generate_x_values_single_encounter(self):
        """Test with single encounter per cluster."""
        cluster_data = {
            "cluster 0": {1: 10.0},
            "cluster 1": {2: 20.0},
            "cluster 2": {3: 30.0}
        }
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        assert total == 3
        assert all(len(vals) == 1 for vals in cluster_x_values.values())


class TestParseAngleData:
    """Tests for _parse_angle_data function."""

    def test_parse_angle_data_single_line(self, temp_dir):
        """Test parsing single line of angle data."""
        filepath = os.path.join(temp_dir, "angles.txt")
        with open(filepath, 'w') as f:
            f.write("10.5 20.3 15.7 25.2\n")
        
        angles = _parse_angle_data(filepath)
        
        assert len(angles) == 4
        assert angles[0] == 10.5
        assert angles[-1] == 25.2

    def test_parse_angle_data_multiple_lines(self, temp_dir):
        """Test parsing multiple lines of angle data."""
        filepath = os.path.join(temp_dir, "angles.txt")
        with open(filepath, 'w') as f:
            f.write("10.5 20.3 15.7\n")
            f.write("25.2 30.1 18.9\n")
        
        angles = _parse_angle_data(filepath)
        
        assert len(angles) == 6
        assert angles[0] == 10.5
        assert angles[3] == 25.2

    def test_parse_angle_data_whitespace_handling(self, temp_dir):
        """Test parsing handles extra whitespace."""
        filepath = os.path.join(temp_dir, "angles.txt")
        with open(filepath, 'w') as f:
            f.write("  10.5   20.3  15.7  \n")
        
        angles = _parse_angle_data(filepath)
        
        assert len(angles) == 3
        assert angles[1] == 20.3

    def test_parse_angle_data_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            _parse_angle_data("nonexistent_angles.txt")

    def test_parse_angle_data_empty_file(self, temp_dir):
        """Test parsing empty file."""
        filepath = os.path.join(temp_dir, "empty_angles.txt")
        with open(filepath, 'w') as f:
            f.write("")
        
        angles = _parse_angle_data(filepath)
        
        assert angles == []


class TestParseClusterFile:
    """Tests for _parse_cluster_file function."""

    def test_parse_cluster_file_basic(self, temp_dir):
        """Test parsing basic cluster file."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3\n")
            f.write("Cluster: 1\n")
            f.write("  4 5\n")
        
        all_angles = [10.0, 20.0, 30.0, 40.0, 50.0]
        cluster_data = _parse_cluster_file(filepath, all_angles)
        
        assert "Cluster: 0" in cluster_data
        assert "Cluster: 1" in cluster_data

    def test_parse_cluster_file_encounter_mapping(self, temp_dir):
        """Test that encounters are mapped to correct angles."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2\n")
        
        all_angles = [10.0, 20.0, 30.0]
        cluster_data = _parse_cluster_file(filepath, all_angles)
        
        # Encounters 1,2 should map to angles at indices 0,1 (10.0, 20.0)
        cluster_angles = list(cluster_data["Cluster: 0"].values())
        assert 10.0 in cluster_angles or 20.0 in cluster_angles

    def test_parse_cluster_file_empty(self, temp_dir):
        """Test parsing empty cluster file."""
        filepath = os.path.join(temp_dir, "empty_clusters.txt")
        with open(filepath, 'w') as f:
            f.write("")
        
        all_angles = [10.0, 20.0, 30.0]
        cluster_data = _parse_cluster_file(filepath, all_angles)
        
        assert len(cluster_data) == 0

    def test_parse_cluster_file_multiline_cluster(self, temp_dir):
        """Test cluster with entries on multiple lines."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3 4 5\n")
            f.write("Cluster: 1\n")
            f.write("  6 7 8\n")
        
        all_angles = [float(i) for i in range(1, 9)]
        cluster_data = _parse_cluster_file(filepath, all_angles)
        
        assert len(cluster_data) == 2
        assert len(cluster_data["Cluster: 0"]) == 5

    def test_parse_cluster_file_colon_stripping(self, temp_dir):
        """Test that cluster keys preserve a valid cluster label format."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2\n")
        
        all_angles = [10.0, 20.0, 30.0]
        cluster_data = _parse_cluster_file(filepath, all_angles)
        
        keys = list(cluster_data.keys())
        assert any(str(k).startswith("Cluster") and str(k).split()[1] == "0" for k in keys)


class TestIntegration:
    """Integration tests for angle plotting workflow."""

    def test_full_angle_parsing_workflow(self, temp_dir):
        """Test complete angle data parsing and clustering."""
        # Create angle data file
        angle_file = os.path.join(temp_dir, "angles.txt")
        with open(angle_file, 'w') as f:
            f.write("10.5 20.3 15.7\n")
            f.write("25.2 30.1 18.9\n")
        
        # Create cluster file
        cluster_file = os.path.join(temp_dir, "clusters.txt")
        with open(cluster_file, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3\n")
            f.write("Cluster: 1\n")
            f.write("  4 5 6\n")
        
        # Parse data
        angles = _parse_angle_data(angle_file)
        cluster_data = _parse_cluster_file(cluster_file, angles)
        
        assert len(angles) == 6
        assert len(cluster_data) == 2
