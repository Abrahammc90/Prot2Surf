"""
Unit tests for extract_C1_C4_filtered.py - Filtering encounter complexes by reference.
"""
import os
import tempfile
import pytest

from extract_C1_C4_filtered import load_complexes_from_files

def test_load_complexes_from_files(tmp_path):
    # Create two mock association files with headers and data
    file1 = tmp_path / "assoc1.txt"
    file2 = tmp_path / "assoc2.txt"
    content1 = """# Header1\n# Header2\n# Header3\n# Header4\n  1.0  2.0  3.0\n  1.5  2.5  3.5\n"""
    content2 = """# Header1\n# Header2\n# Header3\n# Header4\n  2.0  3.0  4.0\n  1.5  2.5  3.5\n"""
    file1.write_text(content1)
    file2.write_text(content2)
    result = load_complexes_from_files([str(file1), str(file2)])
    # Should contain all unique complex lines (excluding headers)
    assert "  1.0  2.0  3.0\n" in result
    assert "  1.5  2.5  3.5\n" in result
    assert "  2.0  3.0  4.0\n" in result
    assert len(result) == 3
