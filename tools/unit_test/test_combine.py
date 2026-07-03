"""Unit tests for combine.py - Combining encounter complex files.

@version 1.0
"""

import sys
import os
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from combine import load_header_from_file, load_complexes_from_file, combine_complexes_removing_duplicates


class TestLoadHeaderFromFile:
    """Tests for load_header_from_file function."""

    def test_load_header_basic(self, sample_association_file):
        """Test loading header from association file."""
        header = load_header_from_file(sample_association_file)
        
        assert len(header) == 4
        assert header[0].startswith('# Header line 1')
        assert all(isinstance(line, str) for line in header)

    def test_load_header_preserves_newlines(self, sample_association_file):
        """Test that header lines preserve newline characters."""
        header = load_header_from_file(sample_association_file)
        
        assert all(line.endswith('\n') for line in header)

    def test_load_header_file_not_found(self):
        """Test error handling for non-existent file."""
        with pytest.raises(FileNotFoundError):
            load_header_from_file("nonexistent_file.txt")

    def test_load_header_short_file(self, temp_dir):
        """Test loading header from file with fewer than 4 lines."""
        filepath = os.path.join(temp_dir, "short_file.txt")
        with open(filepath, 'w') as f:
            f.write("Line 1\nLine 2\n")
        
        header = load_header_from_file(filepath)
        assert len(header) == 2  # Should return available lines


class TestLoadComplexesFromFile:
    """Tests for load_complexes_from_file function."""

    def test_load_complexes_basic(self, sample_association_file):
        """Test loading complexes (skip header)."""
        complexes = load_complexes_from_file(sample_association_file)
        
        assert len(complexes) == 3
        assert all(isinstance(line, str) for line in complexes)

    def test_load_complexes_no_header_lines(self, temp_dir):
        """Test loading complexes from file without 4-line header."""
        filepath = os.path.join(temp_dir, "no_header.txt")
        content = """  1.0  2.0  3.0
  1.5  2.5  3.5
"""
        with open(filepath, 'w') as f:
            f.write(content)
        
        complexes = load_complexes_from_file(filepath)
        # Should skip first 4 lines even if fewer lines exist
        assert len(complexes) == 0 or complexes[0].startswith('  ')

    def test_load_complexes_empty_file(self, temp_dir):
        """Test loading complexes from empty file."""
        filepath = os.path.join(temp_dir, "empty.txt")
        with open(filepath, 'w') as f:
            f.write("")
        
        complexes = load_complexes_from_file(filepath)
        assert complexes == []


class TestCombineComplexesRemovingDuplicates:
    """Tests for combine_complexes_removing_duplicates function."""

    def test_combine_no_duplicates(self, sample_association_file):
        """Test combining files with no duplicates."""
        list1 = ["  1.0  2.0  3.0\n"]
        list2 = ["  1.5  2.5  3.5\n"]
        
        combined = combine_complexes_removing_duplicates([list1, list2])
        
        assert len(combined) == 2
        assert "  1.0  2.0  3.0\n" in combined
        assert "  1.5  2.5  3.5\n" in combined

    def test_combine_with_duplicates(self):
        """Test combining files that contain duplicates."""
        list1 = ["  1.0  2.0  3.0\n", "  2.0  3.0  4.0\n"]
        list2 = ["  2.0  3.0  4.0\n", "  3.0  4.0  5.0\n"]
        
        combined = combine_complexes_removing_duplicates([list1, list2])
        
        assert len(combined) == 3
        assert combined.count("  2.0  3.0  4.0\n") == 1  # Should appear only once

    def test_combine_empty_lists(self):
        """Test combining empty lists."""
        combined = combine_complexes_removing_duplicates([[], []])
        
        assert len(combined) == 0

    def test_combine_single_list(self):
        """Test combining single list."""
        list1 = ["  1.0  2.0  3.0\n", "  1.5  2.5  3.5\n"]
        
        combined = combine_complexes_removing_duplicates([list1])
        
        assert len(combined) == 2

    def test_combine_all_duplicates(self):
        """Test combining lists where all entries are duplicates."""
        list1 = ["  1.0  2.0  3.0\n", "  1.0  2.0  3.0\n"]
        list2 = ["  1.0  2.0  3.0\n"]
        
        combined = combine_complexes_removing_duplicates([list1, list2])
        
        assert len(combined) == 1

    def test_combine_preserves_order(self):
        """Test that combining preserves encounter order."""
        list1 = ["  1.0  2.0  3.0\n", "  2.0  3.0  4.0\n"]
        list2 = ["  3.0  4.0  5.0\n", "  4.0  5.0  6.0\n"]
        
        combined = combine_complexes_removing_duplicates([list1, list2])
        
        # First list entries should appear before second list entries
        indices = {line: i for i, line in enumerate(combined)}
        assert indices["  1.0  2.0  3.0\n"] < indices["  3.0  4.0  5.0\n"]


class TestIntegration:
    """Integration tests for combine workflow."""

    def test_full_combine_workflow(self, temp_dir):
        """Test complete workflow: load headers, load complexes, combine."""
        # Create test files
        file1_path = os.path.join(temp_dir, "file1.txt")
        file2_path = os.path.join(temp_dir, "file2.txt")
        
        header = "# Line 1\n# Line 2\n# Line 3\n# Line 4\n"
        
        with open(file1_path, 'w') as f:
            f.write(header)
            f.write("  1.0  2.0\n")
            f.write("  1.5  2.5\n")
        
        with open(file2_path, 'w') as f:
            f.write(header)
            f.write("  1.5  2.5\n")
            f.write("  2.0  3.0\n")
        
        # Load and combine
        header_lines = load_header_from_file(file1_path)
        complexes1 = load_complexes_from_file(file1_path)
        complexes2 = load_complexes_from_file(file2_path)
        
        combined = combine_complexes_removing_duplicates([complexes1, complexes2])
        
        assert len(header_lines) == 4
        assert len(combined) == 3  # 3 unique complexes
