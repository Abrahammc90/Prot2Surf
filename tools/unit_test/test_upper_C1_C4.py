"""
Unit tests for upper_C1_C4.py - Extraction of uppermost C1/C4 atoms from PDB files.
"""
import os
import tempfile
import pytest

from upper_C1_C4 import is_atom_record, get_atom_name, get_z_coordinate

def test_is_atom_record():
    assert is_atom_record("ATOM      1  C1  GLU A   1   ...")
    assert is_atom_record("HETATM    2  C4  GLU A   2   ...")
    assert not is_atom_record("REMARK   ...")

def test_get_atom_name():
    line = "ATOM      1  C1  GLU A   1       0.000   0.000   0.000  1.00  0.00           C"
    assert get_atom_name(line) == "C1"
    line2 = "ATOM      2  C4  GLU A   2       1.000   2.000   3.000  1.00  0.00           C"
    assert get_atom_name(line2) == "C4"

def test_get_z_coordinate():
    line = "ATOM      1  C1  GLU A   1       0.000   0.000   5.123  1.00  0.00           C"
    assert get_z_coordinate(line) == pytest.approx(5.123, rel=1e-3)
