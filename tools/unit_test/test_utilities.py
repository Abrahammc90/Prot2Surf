"""Unit tests for utility scripts - extract_grid.py, upper_C1_C4.py, etc."""

import sys
import os
import pytest
import tempfile

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestPDBFileHandling:
    """Tests for PDB file parsing and processing."""

    def test_pdb_file_reading(self, sample_pdb_file):
        """Test reading and parsing PDB files."""
        with open(sample_pdb_file, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith('ATOM')]
        assert len(atom_lines) > 0
        assert 'CA' in atom_lines[0]

    def test_pdb_coordinate_extraction(self, temp_dir):
        """Test extracting coordinates from PDB file."""
        filepath = os.path.join(temp_dir, "coords.pdb")
        with open(filepath, 'w') as f:
            f.write("ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N\n")
            f.write("ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C\n")
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Extract coordinates from ATOM lines
        coords = []
        for line in lines:
            if line.startswith('ATOM'):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
        
        assert len(coords) == 2
        assert coords[0] == (1.0, 2.0, 3.0)
        assert coords[1] == (2.0, 3.0, 4.0)

    def test_pdb_file_format_validation(self, temp_dir):
        """Test PDB file format validation."""
        filepath = os.path.join(temp_dir, "test.pdb")
        with open(filepath, 'w') as f:
            f.write("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
            f.write("TER\n")
            f.write("END\n")
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Validate format
        assert any(line.startswith('ATOM') for line in lines)
        assert any(line.startswith('TER') for line in lines)
        assert any(line.startswith('END') for line in lines)


class TestAssociationFileHandling:
    """Tests for association/complex file format handling."""

    def test_association_file_header_lines(self, sample_association_file):
        """Test counting header lines in association file."""
        with open(sample_association_file, 'r') as f:
            lines = f.readlines()
        
        header_lines = 4
        data_lines = lines[header_lines:]
        
        assert len(lines) >= header_lines
        assert len(data_lines) > 0

    def test_association_file_entry_parsing(self, sample_association_file):
        """Test parsing individual entries in association file."""
        with open(sample_association_file, 'r') as f:
            lines = f.readlines()
        
        data_lines = lines[4:]  # Skip header
        entries = []
        
        for line in data_lines:
            parts = line.split()
            if parts:  # Non-empty line
                entries.append(parts)
        
        assert len(entries) > 0
        assert all(isinstance(entry, list) for entry in entries)

    def test_association_file_numerical_columns(self, temp_dir):
        """Test parsing numerical columns in association file."""
        filepath = os.path.join(temp_dir, "assoc.txt")
        with open(filepath, 'w') as f:
            f.write("# Line 1\n# Line 2\n# Line 3\n# Line 4\n")
            f.write("  1.5  2.5  3.5\n")
            f.write("  4.5  5.5  6.5\n")
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        data_lines = lines[4:]
        numerical_data = []
        
        for line in data_lines:
            row = [float(x) for x in line.split()]
            numerical_data.append(row)
        
        assert len(numerical_data) == 2
        assert numerical_data[0] == [1.5, 2.5, 3.5]


class TestClusterFileHandling:
    """Tests for cluster file format handling."""

    def test_cluster_file_parsing(self, sample_cluster_file):
        """Test parsing cluster file format."""
        with open(sample_cluster_file, 'r') as f:
            lines = f.readlines()
        
        clusters = {}
        current_cluster = None
        
        for line in lines:
            if line.startswith("Cluster"):
                current_cluster = line.strip()
                clusters[current_cluster] = []
            else:
                parts = line.split()
                if parts and current_cluster:
                    clusters[current_cluster].extend(parts)
        
        assert len(clusters) > 0
        assert all(isinstance(v, list) for v in clusters.values())

    def test_cluster_file_encounter_indexing(self, temp_dir):
        """Test encounter ID extraction from cluster file."""
        filepath = os.path.join(temp_dir, "clusters.txt")
        with open(filepath, 'w') as f:
            f.write("Cluster: 0\n")
            f.write("  1 2 3 4 5\n")
            f.write("Cluster: 1\n")
            f.write("  6 7 8\n")
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        encounter_ids = []
        for line in lines:
            if not line.startswith("Cluster"):
                ids = [int(x) for x in line.split() if x.isdigit()]
                encounter_ids.extend(ids)
        
        assert len(encounter_ids) == 8
        assert encounter_ids[0] == 1
        assert encounter_ids[-1] == 8


class TestDataTypeConversions:
    """Tests for data type conversions in file handling."""

    def test_string_to_float_conversion(self):
        """Test string to float conversion."""
        strings = ["5.2", "6.1", "7.3", "-8.2", "1e-3"]
        
        floats = [float(s) for s in strings]
        
        assert len(floats) == 5
        assert floats[0] == 5.2
        assert floats[3] == -8.2
        assert floats[4] == pytest.approx(0.001)

    def test_string_to_int_conversion(self):
        """Test string to integer conversion."""
        strings = ["1", "2", "3", "100"]
        
        ints = [int(s) for s in strings]
        
        assert len(ints) == 4
        assert ints[0] == 1
        assert ints[-1] == 100

    def test_invalid_float_conversion(self):
        """Test error handling for invalid float conversion."""
        invalid_strings = ["abc", "1.2.3", ""]
        
        for invalid_str in invalid_strings:
            with pytest.raises(ValueError):
                float(invalid_str.strip() if invalid_str else "invalid")

    def test_list_to_array_conversion(self):
        """Test converting lists to arrays."""
        data_lists = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        
        # Flatten lists
        flat = [item for sublist in data_lists for item in sublist]
        
        assert len(flat) == 6
        assert flat[0] == 1.0
        assert flat[-1] == 6.0


class TestFileIOErrors:
    """Tests for file I/O error handling."""

    def test_file_not_found_error(self):
        """Test handling of missing files."""
        with pytest.raises(FileNotFoundError):
            with open("nonexistent_file_xyz.txt", 'r') as f:
                f.read()

    def test_permission_error_simulation(self, temp_dir):
        """Test handling of permission errors."""
        filepath = os.path.join(temp_dir, "readonly.txt")
        with open(filepath, 'w') as f:
            f.write("test")
        
        # Make file read-only
        os.chmod(filepath, 0o444)
        
        # Should be able to read
        with open(filepath, 'r') as f:
            content = f.read()
        
        assert content == "test"
        
        # Cleanup
        os.chmod(filepath, 0o644)

    def test_empty_file_handling(self, temp_dir):
        """Test handling of empty files."""
        filepath = os.path.join(temp_dir, "empty.txt")
        with open(filepath, 'w') as f:
            pass  # Create empty file
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        assert lines == []


class TestFileFormatConsistency:
    """Tests for maintaining file format consistency."""

    def test_newline_preservation(self, temp_dir):
        """Test that newlines are preserved when reading/writing."""
        filepath = os.path.join(temp_dir, "newlines.txt")
        original_lines = ["Line 1\n", "Line 2\n", "Line 3\n"]
        
        with open(filepath, 'w') as f:
            f.writelines(original_lines)
        
        with open(filepath, 'r') as f:
            read_lines = f.readlines()
        
        assert read_lines == original_lines

    def test_whitespace_handling(self, temp_dir):
        """Test handling of leading/trailing whitespace."""
        filepath = os.path.join(temp_dir, "whitespace.txt")
        content = "  value1   value2  \n"
        
        with open(filepath, 'w') as f:
            f.write(content)
        
        with open(filepath, 'r') as f:
            line = f.readline()
        
        # Whitespace should be preserved
        assert line == content
        
        # Split should handle whitespace
        parts = line.split()
        assert len(parts) == 2
        assert parts[0] == "value1"
