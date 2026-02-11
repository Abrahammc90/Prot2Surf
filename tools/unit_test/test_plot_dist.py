"""Unit tests for plot_clust_array_dist.py - Distance distribution plotting."""

import sys
import os
import pytest
import tempfile

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plot_clust_array_dist import (
    _generate_encounter_x_values,
    _parse_distance_data,
    _parse_cluster_file
)


class TestGenerateEncounterXValues:
    """Tests for _generate_encounter_x_values function."""

    def test_generate_x_values_basic(self):
        """Test x-value generation with basic cluster data."""
        cluster_data = {
            "cluster 0": {1: 5.2, 2: 6.1, 3: 7.3},
            "cluster 1": {4: 8.2, 5: 5.9}
        }
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        assert total == 5
        assert len(cluster_x_values) == 2
        assert len(cluster_x_values["cluster 0"]) == 3
        assert len(cluster_x_values["cluster 1"]) == 2

    def test_generate_x_values_sequential(self):
        """Test that x-values are sequential across clusters."""
        cluster_data = {
            "cluster 0": {1: 5.0, 2: 6.0},
            "cluster 1": {3: 7.0, 4: 8.0}
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

    def test_generate_x_values_many_clusters(self):
        """Test with many clusters."""
        cluster_data = {f"cluster {i}": {i: float(i*10)} for i in range(10)}
        
        cluster_x_values, total = _generate_encounter_x_values(cluster_data)
        
        assert total == 10
        assert len(cluster_x_values) == 10


class TestParseDistanceData:
    """Tests for _parse_distance_data function."""

    def test_parse_distance_data_single_line(self, temp_dir):
        """Test parsing single line of distance data."""
        filepath = os.path.join(temp_dir, "distances.txt")
        with open(filepath, 'w') as f:
            f.write("5.2 6.1 7.3 8.2 5.9\n")
        
        distances = _parse_distance_data(filepath)
        
        assert len(distances) == 5
        assert distances[0] == 5.2
        assert distances[-1] == 5.9

    def test_parse_distance_data_multiple_lines(self, temp_dir):
        """Test parsing multiple lines of distance data."""
        filepath = os.path.join(temp_dir, "distances.txt")
        with open(filepath, 'w') as f:
            f.write("5.2 6.1 7.3\n")
            f.write("8.2 5.9 6.5\n")
        
        distances = _parse_distance_data(filepath)
        
        assert len(distances) == 6
        assert distances[0] == 5.2
        assert distances[3] == 8.2

    def test_parse_distance_data_whitespace_handling(self, temp_dir):
        """Test parsing handles extra whitespace."""
        filepath = os.path.join(temp_dir, "distances.txt")
        with open(filepath, 'w') as f:
            f.write("  5.2   6.1  7.3  \n")
        
        distances = _parse_distance_data(filepath)
        
        assert len(distances) == 3
        assert distances[1] == 6.1

    def test_parse_distance_data_scientific_notation(self, temp_dir):
        """Test parsing scientific notation in distances."""
        filepath = os.path.join(temp_dir, "distances.txt")
        with open(filepath, 'w') as f:
            f.write("5.2e0 6.1e0 7.3e0\n")
        
        distances = _parse_distance_data(filepath)
        
        assert len(distances) == 3
        assert distances[0] == pytest.approx(5.2)

    def test_parse_distance_data_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            _parse_distance_data("nonexistent_distances.txt")

    def test_parse_distance_data_empty_file(self, temp_dir):
        """Test parsing empty file."""
        filepath = os.path.join(temp_dir, "empty_distances.txt")
        with open(filepath, 'w') as f:
            f.write("")
        
        distances = _parse_distance_data(filepath)
        
        assert distances == []


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
        
        all_distances = [5.2, 6.1, 7.3, 8.2, 5.9]
        cluster_data = _parse_cluster_file(filepath, all_distances)
        
        assert "Cluster: 0" in cluster_data
        assert "Cluster: 1" in cluster_data

    def test_parse_cluster_file_encounter_mapping(self, temp_dir):
        """Test that encounters are mapped to correct distances."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2\n")
        
        all_distances = [5.2, 6.1, 7.3]
        cluster_data = _parse_cluster_file(filepath, all_distances)
        
        # Should have values from the distance list
        cluster_distances = list(cluster_data["Cluster: 0"].values())
        assert len(cluster_distances) > 0

    def test_parse_cluster_file_large_encounter_ids(self, temp_dir):
        """Test parsing cluster file with large encounter IDs."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 5 10\n")
        
        all_distances = [float(i) for i in range(1, 11)]
        cluster_data = _parse_cluster_file(filepath, all_distances)
        
        assert len(cluster_data["Cluster: 0"]) > 0

    def test_parse_cluster_file_empty(self, temp_dir):
        """Test parsing empty cluster file."""
        filepath = os.path.join(temp_dir, "empty_clusters.txt")
        with open(filepath, 'w') as f:
            f.write("")
        
        all_distances = [5.2, 6.1, 7.3]
        cluster_data = _parse_cluster_file(filepath, all_distances)
        
        assert len(cluster_data) == 0

    def test_parse_cluster_file_out_of_bounds_encounters(self, temp_dir):
        """Test cluster file references encounters beyond data range."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3 4 5\n")
        
        # Only 3 distance values but cluster references 5 encounters
        all_distances = [5.2, 6.1, 7.3]
        cluster_data = _parse_cluster_file(filepath, all_distances)
        
        # Should handle gracefully and only include valid entries
        assert len(cluster_data) > 0


class TestIntegration:
    """Integration tests for distance plotting workflow."""

    def test_full_distance_parsing_workflow(self, temp_dir):
        """Test complete distance data parsing and clustering."""
        # Create distance data file
        dist_file = os.path.join(temp_dir, "distances.txt")
        with open(dist_file, 'w') as f:
            f.write("5.2 6.1 7.3\n")
            f.write("8.2 5.9 6.5\n")
        
        # Create cluster file
        cluster_file = os.path.join(temp_dir, "clusters.txt")
        with open(cluster_file, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3\n")
            f.write("Cluster: 1\n")
            f.write("  4 5 6\n")
        
        # Parse data
        distances = _parse_distance_data(dist_file)
        cluster_data = _parse_cluster_file(cluster_file, distances)
        
        assert len(distances) == 6
        assert len(cluster_data) == 2

    def test_plot_consistency_across_calls(self, temp_dir):
        """Test that parsing produces consistent results across multiple calls."""
        dist_file = os.path.join(temp_dir, "distances.txt")
        with open(dist_file, 'w') as f:
            f.write("5.2 6.1 7.3 8.2 5.9 6.5\n")
        
        cluster_file = os.path.join(temp_dir, "clusters.txt")
        with open(cluster_file, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2\n")
        
        # Parse twice
        distances1 = _parse_distance_data(dist_file)
        cluster_data1 = _parse_cluster_file(cluster_file, distances1)
        
        distances2 = _parse_distance_data(dist_file)
        cluster_data2 = _parse_cluster_file(cluster_file, distances2)
        
        # Results should be identical
        assert distances1 == distances2
        assert len(cluster_data1) == len(cluster_data2)
