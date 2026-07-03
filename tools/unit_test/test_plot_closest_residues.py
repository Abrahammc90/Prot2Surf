"""Unit tests for plot_closest_residues.py.

@version 1.0
"""

import os
import sys

import matplotlib

# Use a non-interactive backend for headless test environments.
matplotlib.use("Agg")

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from plot_closest_residues import contacts, plot


class TestContacts:
    """Tests for closest-residue contact parsing."""

    def test_read_contacts_and_residue_range(self, temp_dir):
        """Parse contact counts and first/last residues from input files."""
        input_file = os.path.join(temp_dir, "closest.txt")
        pdb_file = os.path.join(temp_dir, "protein.pdb")

        with open(input_file, "w") as f:
            f.write("Residue 27: 101 102 103\n")
            f.write("Residue 30: 200\n")

        with open(pdb_file, "w") as f:
            f.write("ATOM      1  N   ALA A  25       0.000   0.000   0.000  1.00  0.00           N\n")
            f.write("ATOM      2  CA  SER A  27       1.000   0.000   0.000  1.00  0.00           C\n")
            f.write("ATOM      3  C   TYR A  30       1.500   0.500   0.000  1.00  0.00           C\n")
            f.write("ATOM      4  O   ALA A  40       2.000   1.000   0.000  1.00  0.00           O\n")
            f.write("END\n")

        residue_contacts = contacts()
        residue_contacts.read(input_file, pdb_file)

        assert residue_contacts.first == 25
        assert residue_contacts.last == 40
        assert residue_contacts.contacts[27] == 3
        assert residue_contacts.contacts[30] == 1


class TestPlotOutputs:
    """Tests for plotting and output generation."""

    def test_generate_writes_png(self, temp_dir):
        """Generate should save a non-empty PNG file."""
        output_base = os.path.join(temp_dir, "closest_plot.dat")
        expected_png = os.path.join(temp_dir, "closest_plot.png")
        data = {27: 3, 30: 1, 31: 5}
        residue_names = {27: "SER", 30: "TYR", 31: "ASP"}

        plotter = plot()
        plotter.generate(data, residue_names, 25, 40, "Test Plot", output_base)

        assert os.path.exists(expected_png)
        assert os.path.getsize(expected_png) > 0

    def test_write_data_file(self, temp_dir):
        """write_data should export residue/count lines."""
        output_data = os.path.join(temp_dir, "closest_plot.dat")
        data = {27: 3, 30: 1}
        residue_names = {27: "SER", 30: "TYR"}

        plotter = plot()
        plotter.write_data(data, residue_names, output_data)

        with open(output_data, "r") as f:
            content = f.read()

        assert "#ResidueNumber" in content
        assert "27" in content
        assert "3" in content
        assert "30" in content
        assert "1" in content

    def test_generate_with_empty_contacts_does_not_create_png(self, temp_dir):
        """Generate should return early when no valid contacts are provided."""
        output_png = os.path.join(temp_dir, "empty_plot.png")

        plotter = plot()
        plotter.generate({}, {}, 1, 1, "Empty", output_png)

        assert not os.path.exists(output_png)
