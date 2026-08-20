"""
@file plot_closest_residues.py
@brief Plot closest-residue encounter counts and export tabular summaries.
@version 2.0
@date 2026-08-11
@par Usage
python plot_closest_residues.py Cu_z_cluster*_complexes_closest_residues.txt p2_noh.pdb Cu_z_cluster*_complexes.dat <plot_title>

@par Output Files:
    <output_prefix>.png: Publication-quality plot figure
    <output_prefix>.dat: Data file with residue numbers, names, and occurrences

@par Description
Reads closest-residue encounter data and protein residue labels from a PDB file,
then creates a color-coded vertical-line plot and an optional text summary.
"""
import sys
import math
import matplotlib.pyplot as plt
import numpy as np


class contacts:
    """Container for residue-contact counts and residue-name mapping."""
    def __init__(self):
        self.contacts = {}
        self.resid_to_name = {}

    def read(self, input_file, pdb_file):
        """Read residue metadata from PDB and contact counts from input file.

        Args:
            input_file (str): Closest-residue input file.
            pdb_file (str): Protein PDB file with residue IDs and names.
        """
        # --- Read PDB and store residue names ---
        with open(pdb_file, "r") as pdb:
            self.first = None
            self.last = None

            for line in pdb:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    resid = int(line[22:26])
                    resname = line[17:20].strip()

                    if resid not in self.resid_to_name:
                        self.resid_to_name[resid] = resname

                    if self.first is None:
                        self.first = resid

                    self.last = resid

        # --- Read contact file (only lines starting with "Residue") ---
        with open(input_file, "r") as inp:
            for line in inp:
                if not line.startswith("Residue"):
                    continue

                line_splitted = line.split()
                resid_str = line_splitted[1]           # e.g., "152:11367" or "74:"
                resid = int(resid_str.split(":")[0])   # take only the number before colon

                # Count occurrences: if no numbers after colon, count as 1
                encounters_list = line_splitted[2:]
                tot_encounters = len(encounters_list) if encounters_list else 1

                # Only count residues that exist in PDB
                if resid in self.resid_to_name:
                    self.contacts[resid] = tot_encounters


class plot:
    """Plotting and data-export utilities for closest-residue analysis."""

    def generate(self, residue_contacts, residue_names,
                 first_resid, last_resid, title, plot_filename):
        """Generate and save the closest-residue distribution plot.

        Args:
            residue_contacts (dict): Residue ID to encounter count mapping.
            residue_names (dict): Residue ID to residue name mapping.
            first_resid (int): First residue ID observed in the PDB.
            last_resid (int): Last residue ID observed in the PDB.
            title (str): Plot title.
            plot_filename (str): Output base filename ('.png' is enforced).
        """

        if not residue_contacts:
            print("No valid residue contact data found.")
            return

        keys = sorted(residue_contacts.keys())
        values = [residue_contacts[k] for k in keys]
        labels = [residue_names.get(k, "UNK") + str(k) for k in keys]

        max_value = max(values)

        # ---------- Spacing system ----------
        spacing = 1.2
        x_pos = np.arange(len(keys)) * spacing

        distinct_colors = [
            "#1f77b4",  # blue
            "#ff7f0e",  # orange
            "#2ca02c",  # green
            "#d62728",  # red
            "#9467bd",  # purple
            "#8c564b",  # brown
            "#e377c2",  # pink
            "#7f7f7f",  # gray
            "#bcbd22",  # olive
            "#17becf",  # cyan

            "#393b79",  # dark blue
            "#637939",  # dark green
            "#8c6d31",  # dark gold
            "#843c39",  # dark red
            "#7b4173",  # dark purple
            "#3182bd",  # medium blue
            "#e6550d",  # dark orange
            "#31a354",  # medium green
            "#756bb1",  # medium purple
            "#636363",  # dark gray

            "#6baed6",  # light blue
            "#fd8d3c",  # light orange
            "#74c476",  # light green
            "#fb6a6a",  # light red
            "#9e9ac8",  # light purple
            "#a6761d",  # gold/brown
            "#e78ac3",  # light pink
            "#969696",  # light gray
            "#bdbf6f",  # light olive
            "#5bc0de",  # light cyan
        ]

        colors = [distinct_colors[i % len(distinct_colors)] for i in range(len(keys))]

        # Plot vertical lines using spaced positions
        for x, value, color, label in zip(
            x_pos, values, colors, labels
        ):
            plt.vlines(
                x=x,
                ymin=0,
                ymax=value,
                label=label,
                color=color,
                linewidth=3.0
            ) # Adjusted linewidth width for better visibility 

        # Axis limits
        plt.xlim(x_pos[0] - spacing,x_pos[-1] + spacing)

        plt.ylim(0, max_value * 1.05)

        # Labels and formatting
#        ncol = (len(keys) - 1) // 20 + 1

        plt.xlabel("Residue", fontsize=20)

        plt.ylabel("Number of encounters", fontsize=20)

        plt.xticks(x_pos,labels, rotation=90, ha="center", fontsize=14) # if more vertical lines are observed(overlapped residue name), decrease fontsize to 12

        plt.yticks(fontsize=16)

        plt.title(title, fontsize=24)

        plt.grid(True)

        plt.tight_layout()

        # Save figure: remove .dat from filename if exists
        png_filename = plot_filename.replace(".dat", "") + ".png"

        plt.savefig(png_filename, dpi=300)

        plt.close()

    def write_data(self, residue_contacts, residue_labels, filename):
        """Write residue-contact counts to a .dat file.

        Args:
            residue_contacts (dict): Residue ID to encounter count mapping.
            residue_labels (dict): Residue ID to residue name mapping.
            filename (str): Output filename or prefix for '.dat' generation.
        """
        # Fix filename: add .dat only if missing
        filename_out = filename if filename.endswith(".dat") else filename + ".dat"

        with open(filename_out, "w") as f:
            f.write("#ResidueNumber ResidueName Occurrences\n")
            for resid in sorted(residue_contacts.keys()):
                count = residue_contacts[resid]
                resname = residue_labels.get(resid, "UNK")
                f.write(f"{resid:>7d} {resname}{resid:<4d} {count:>7d}\n")


if __name__ == "__main__":

    input_filename = sys.argv[1]
    pdb_filename = sys.argv[2]
    plot_filename = sys.argv[3]
    title = " ".join(sys.argv[4:])

    residues = contacts()

    residues.read(input_filename, pdb_filename)

    p = plot()

    p.generate(
        residues.contacts,
        residues.resid_to_name,
        residues.first,
        residues.last,
        title,
        plot_filename
    )

    p.write_data(
        residues.contacts,
        residues.resid_to_name,
        plot_filename
    )
