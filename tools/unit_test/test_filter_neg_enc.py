"""
Unit tests for filter_neg_enc.py - Atom class and filtering logic.
"""
import pytest
from filter_neg_enc import Atom

def test_atom_initialization():
    atom = Atom()
    assert atom.crds == []
    assert atom.mass == 0
    assert atom.resname == ""
    assert atom.resid == 0
    assert atom.id == 0
    assert atom.name == ""
    assert atom.element == ""
    assert atom.chainID == ""

# Additional tests for Atom methods can be added if methods are imported in __init__
