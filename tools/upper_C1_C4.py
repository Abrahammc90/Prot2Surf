"""!
@file upper_C1_C4.py
@brief Extract uppermost C1/C4 atoms from cellulose PDB structures.
@author Abraham Muñiz-Chicharro
@version 1.0
@date 2026-04-05
@par Usage
python upper_C1_C4.py [input_pdb_file]

@par Usage flags
- Positional @c input_pdb_file : input cellulose PDB structure.

Extract Upper C1/C4 Atoms Script

This script extracts the uppermost C1 and C4 atoms from cellulose PDB structures.
These atoms are critical binding sites in cellulose. The script identifies atoms
near the maximum Z-coordinate (highest vertical position) and saves them to a new PDB file.

Features:
    - PDB file parsing
    - Atom coordinate extraction
    - Maximum Z-coordinate identification (0.1 Ų tolerance)
    - Selective atom filtering for C1 and C4
    - Preserves PDB format (ATOM/HETATM records with TER and END markers)
    - Terminal record handling

Input:
    PDB file with atom coordinates

Output:
    upper_C1_C4.pdb - New PDB file containing only C1 and C4 atoms 
                      within 0.1 Ų of the maximum Z-coordinate

Usage:
    python upper_C1_C4.py [input_pdb_file]

Example:
    python upper_C1_C4.py cellulose_structure.pdb

Implementation:
    1. Parse input PDB file
    2. Find maximum Z-coordinate among C1 and C4 atoms
    3. Filter atoms within 0.1 Ų of maximum Z
    4. Preserve PDB record structure (ATOM/HETATM, TER, END)
    5. Write filtered structure to output file

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys

# Constants
ATOM_RECORD_TYPES = ("ATOM", "HETATM")
TARGET_ATOMS = {"C1", "C4"}
Z_THRESHOLD = 0.1
OUTPUT_FILE = "upper_C1_C4.pdb"

def is_atom_record(line):
    """Check if line is an ATOM or HETATM record."""
    return line.startswith(ATOM_RECORD_TYPES)

def get_atom_name(line):
    """Extract atom name from PDB line."""
    return line[12:16].strip()

def get_z_coordinate(line):
    """Extract Z-coordinate from PDB line."""
    return float(line[46:54])

def extract_target_atom_z_coordinates(pdb_lines):
    """Extract Z-coordinates of all C1 and C4 atoms."""
    z_coordinates = []
    for line in pdb_lines:
        if is_atom_record(line):
            atom_name = get_atom_name(line)
            if atom_name in TARGET_ATOMS:
                z_coordinates.append(get_z_coordinate(line))
    return z_coordinates

def filter_upper_atoms(pdb_lines, max_z):
    """Filter atoms to keep only C1/C4 within Z threshold of maximum."""
    filtered_atoms = []
    atoms_found = False
    
    for line in pdb_lines:
        if is_atom_record(line):
            atom_name = get_atom_name(line)
            if atom_name in TARGET_ATOMS:
                z = get_z_coordinate(line)
                if z > max_z - Z_THRESHOLD:
                    filtered_atoms.append(line)
                    atoms_found = True
        elif line.startswith("TER") and atoms_found:
            filtered_atoms.append(line)
            atoms_found = False
        elif line.startswith("END"):
            filtered_atoms.append(line)
    
    return filtered_atoms

def main():
    """
    Main execution function for extracting upper C1/C4 atoms.
    
    Command-line arguments:
        argv[1]: Input PDB file path
        
    Workflow:
        1. Read input PDB file
        2. Extract Z-coordinates of all C1 and C4 atoms
        3. Determine maximum Z-coordinate
        4. Filter atoms: keep C1/C4 atoms with Z >= (max_Z - 0.1)
        5. Preserve terminal records (TER/END)
        6. Write filtered structure to upper_C1_C4.pdb
    """
    input_pdb_file = sys.argv[1]
    
    with open(input_pdb_file) as f:
        pdb_lines = f.readlines()
    
    z_coordinates = extract_target_atom_z_coordinates(pdb_lines)
    max_z = max(z_coordinates)
    
    filtered_atoms = filter_upper_atoms(pdb_lines, max_z)
    
    with open(OUTPUT_FILE, "w") as output:
        output.writelines(filtered_atoms)

if __name__ == "__main__":
    main()