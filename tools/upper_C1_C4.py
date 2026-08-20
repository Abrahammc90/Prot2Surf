"""!
@file upper_C1_C4.py
@brief Extract the top-layer C1/C4 atoms bonded to the "down" O4 atoms.
@author Abraham Muñiz-Chicharro
@version 2.0
@date 2026-08-20
@par Usage
python upper_C1_C4.py [input_pdb_file]

@par Usage flags
- Positional @c input_pdb_file : input PDB structure (cellulose or chitin).

Extract Upper C1/C4 Atoms Script

The surface is made of several stacked layers of chains. Within a layer, the
O4 atoms alternate up/down along the chain. This script:
    1. Identifies the topmost layer (highest Z), grouping chains by the
       Z-coordinate of their ROH (reducing-end) atom. Chains are delimited
       by TER records.
    2. Within that top layer, finds every O4 atom that sits "down" - i.e.
       its Z-coordinate is lower than its own residue's C4.
    3. Finds the C1 and C4 atoms bonded (within BOND_CUTOFF) to those
       down-O4 atoms and writes them to the output PDB file.

Input:
    PDB file with atom coordinates

Output:
    upper_C1_C4.pdb - New PDB file containing the C1/C4 atoms of the top
                      layer that are bonded to a "down" O4 atom

Usage:
    python upper_C1_C4.py [input_pdb_file]

Example:
    python upper_C1_C4.py cellulose_noh.pdb

Author: Abraham Muñiz-Chicharro
Version: 2.0
"""

import sys
import math

ATOM_RECORD_TYPES = ("ATOM", "HETATM")
BOND_CUTOFF = 2.5      # Å, max distance to consider a C1/C4 bonded to a "down" O4
LAYER_GAP = 1.0         # Å, gap between ROH Z-values used to split chains into layers
OUTPUT_FILE = "upper_C1_C4.pdb"

def is_atom_record(line):
    """Check if line is an ATOM or HETATM record."""
    return line.startswith(ATOM_RECORD_TYPES)

def get_atom_name(line):
    """Extract atom name from PDB line."""
    return line[12:16].strip()

def get_res_name(line):
    """Extract residue name from PDB line."""
    return line[17:20].strip()

def get_res_seq(line):
    """Extract residue sequence number from PDB line."""
    return line[22:26].strip()

def get_coordinates(line):
    """Extract (x, y, z) coordinates from PDB line."""
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))

def distance(a, b):
    """Euclidean distance between two coordinate tuples."""
    return math.sqrt(sum((i - j) ** 2 for i, j in zip(a, b)))

def split_into_chains(pdb_lines):
    """Split PDB lines into chains, delimited by TER records."""
    chains = []
    current = []
    for line in pdb_lines:
        if line.startswith("END"):
            continue
        current.append(line)
        if line.startswith("TER"):
            chains.append(current)
            current = []
    if current:
        chains.append(current)
    return chains

def chain_reference_z(chain_lines):
    """Z-coordinate of the chain's ROH (reducing-end) atom, used to place it in a layer."""
    for line in chain_lines:
        if is_atom_record(line) and get_res_name(line) == "ROH":
            return get_coordinates(line)[2]
    return None

def top_layer_chains(chains):
    """Group chains into Z layers by their ROH atom and return the chains in the topmost layer."""
    dated = [(chain_reference_z(chain), chain) for chain in chains]
    dated = [item for item in dated if item[0] is not None]
    dated.sort(key=lambda item: item[0])

    layers = [[dated[0]]]
    for z, chain in dated[1:]:
        if z - layers[-1][-1][0] > LAYER_GAP:
            layers.append([])
        layers[-1].append((z, chain))

    top_layer = max(layers, key=lambda layer: max(z for z, _ in layer))
    return [chain for _, chain in top_layer]

def parse_residues(chains):
    """Map (chain_index, resSeq) -> {atom_name: (line, coords)} for a set of chains."""
    residues = {}
    for chain_index, chain_lines in enumerate(chains):
        for line in chain_lines:
            if is_atom_record(line):
                res_id = (chain_index, get_res_seq(line))
                residues.setdefault(res_id, {})[get_atom_name(line)] = (line, get_coordinates(line))
    return residues

def find_down_O4(residues):
    """O4 atoms whose Z is below their own residue's C4 (the 'down' side of the pattern)."""
    down_o4 = []
    for atoms in residues.values():
        if "C4" in atoms and "O4" in atoms:
            _, c4_coord = atoms["C4"]
            o4_line, o4_coord = atoms["O4"]
            if o4_coord[2] < c4_coord[2]:
                down_o4.append((o4_line, o4_coord))
    return down_o4

def find_bonded_C1_C4(chains, down_o4):
    """C1/C4 atoms within BOND_CUTOFF of any 'down' O4, restricted to the given chains."""
    bonded = []
    seen = set()
    for chain_lines in chains:
        for line in chain_lines:
            if not is_atom_record(line) or get_atom_name(line) not in ("C1", "C4"):
                continue
            coord = get_coordinates(line)
            if line in seen:
                continue
            for _, o4_coord in down_o4:
                if distance(coord, o4_coord) <= BOND_CUTOFF:
                    bonded.append(line)
                    seen.add(line)
                    break
    return bonded

def write_output(chains, bonded_atoms, output_file):
    """Write the bonded C1/C4 atoms, preserving TER markers per chain and a final END."""
    bonded_set = set(bonded_atoms)
    output_lines = []

    for chain_lines in chains:
        wrote_atom = False
        for line in chain_lines:
            if line in bonded_set:
                output_lines.append(line)
                wrote_atom = True
            elif line.startswith("TER") and wrote_atom:
                output_lines.append(line)

    output_lines.append("END\n")

    with open(output_file, "w") as f:
        f.writelines(output_lines)

def main():
    """
    Main execution function for extracting the top-layer C1/C4 atoms.

    Command-line arguments:
        argv[1]: Input PDB file path

    Workflow:
        1. Split the structure into chains (delimited by TER)
        2. Group chains into layers by their ROH Z-coordinate; keep the top layer
        3. Within the top layer, find O4 atoms below their own residue's C4
        4. Find C1/C4 atoms bonded (within BOND_CUTOFF) to those down-O4 atoms
        5. Write the result to upper_C1_C4.pdb
    """
    input_pdb_file = sys.argv[1]

    with open(input_pdb_file) as f:
        pdb_lines = f.readlines()

    chains = split_into_chains(pdb_lines)
    top_chains = top_layer_chains(chains)

    residues = parse_residues(top_chains)
    down_o4 = find_down_O4(residues)

    bonded_atoms = find_bonded_C1_C4(top_chains, down_o4)

    write_output(top_chains, bonded_atoms, OUTPUT_FILE)

    print(f"Top layer: {len(top_chains)} chains, {len(down_o4)} down-O4 atoms, "
          f"{len(bonded_atoms)} bonded C1/C4 atoms written to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
