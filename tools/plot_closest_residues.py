"""!
@file plot_closest_residues.py
@brief Plot closest-residue encounter counts and export tabular summaries.
@version 1.0
@date 2026-06-07
@par Usage
python plot_closest_residues.py Cu_z_cluster*_complexes_closest_residues.txt p2_noh.pdb Cu_z_cluster*_complexes.dat <plot_title>

@par Output Files
- <output_prefix>.png: publication-quality plot figure.
- <output_prefix>.dat: residue numbers, names, and occurrence counts.

@par Description
Reads closest-residue encounter data and protein residue labels from a PDB file,
then creates a color-coded vertical-line plot and an optional text summary.
"""
#!/usr/bin/env python3
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
                resid_str = line_splitted[1]             # e.g., "152:11367" or "74:"
                resid = int(resid_str.split(":")[0])     # take only the number before colon

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
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
        ]
        colors = [distinct_colors[i % len(distinct_colors)] for i in range(len(keys))]

        # Precompute legend grid and scale figure height with legend rows.
        base_width = 10.0
        base_height = 6.0
        longest_label = max(len(label) for label in labels)
        approx_label_px = max(90, longest_label * 7 + 35)
        usable_width_px = base_width * 100.0 * 0.9
        max_cols = max(1, int(usable_width_px // approx_label_px))
        ncol = min(len(keys), max_cols)
        nrows = int(math.ceil(len(keys) / float(ncol)))

        fig_height = base_height + max(0, nrows - 1) * 1.1
        fig, ax = plt.subplots(figsize=(base_width, fig_height))

        # Plot vertical lines using spaced positions
        bar_width = 1.2 # Adjusted bar width for better visibility
        for x, value, color, label in zip(x_pos, values, colors, labels):
            ax.vlines(x=x, ymin=0, ymax=value, label=label, color=color, linewidth=bar_width)

        # Axis limits
        ax.set_xlim(x_pos[0] - spacing, x_pos[-1] + spacing)
        ax.set_ylim(0, max_value * 1.05)

        # Labels and formatting
        bottom_margin = min(0.72, 0.24 + 0.06 * nrows)
        ax.set_xlabel("Residue", fontsize=20)
        ax.set_ylabel("Number of encounters", fontsize=20)
        ax.set_xticks(x_pos, labels, rotation=90, ha='center', fontsize=14)
        ax.tick_params(axis='y', labelsize=16)
        ax.set_title(title, fontsize=24)

        handles, legend_labels = ax.get_legend_handles_labels()
        fig.legend(
            handles,
            legend_labels,
            ncol=ncol,
            loc="lower left",
            bbox_to_anchor=(0.02, 0.01, 0.96, 0.1),
            mode="expand",
            fontsize=12,
            frameon=True,
            columnspacing=1.0,
            handletextpad=0.4,
            handlelength=2.0,
            borderaxespad=0.0,
        )
        ax.grid(True)
        fig.subplots_adjust(bottom=bottom_margin)

        # Save figure: remove .dat from filename if exists
        png_filename = plot_filename.replace(".dat", "") + ".png"
        fig.savefig(png_filename, dpi=300)
        plt.close(fig)

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


if __name__ == '__main__':
    input_filename = sys.argv[1]
    pdb_filename = sys.argv[2]
    plot_filename = sys.argv[3]
    title = " ".join(sys.argv[4:])

    residues = contacts()
    residues.read(input_filename, pdb_filename)

    p = plot()
    p.generate(residues.contacts,
               residues.resid_to_name,
               residues.first,
               residues.last,
               title,
               plot_filename)

    p.write_data(residues.contacts,
                 residues.resid_to_name,
                 plot_filename)
