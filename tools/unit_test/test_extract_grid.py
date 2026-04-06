"""
Unit tests for extract_grid.py - Grid extraction and parsing utilities.
"""
import os
import tempfile
import numpy as np
import pytest

from extract_grid import parse_uhbd_header

def test_parse_uhbd_header_basic():
    header = (
        "   65   65   65   1.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00"
    )
    result = parse_uhbd_header(header)
    assert result['Nx'] == 65
    assert result['Ny'] == 65
    assert result['Nz'] == 65
    assert result['delta'] == pytest.approx(1.0)
    assert result['origin_x'] == pytest.approx(0.0)
    assert result['origin_y'] == pytest.approx(0.0)
    assert result['origin_z'] == pytest.approx(0.0)

# Additional tests for other functions can be added as needed.
