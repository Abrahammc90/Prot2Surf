"""
Filter Negative Encounter Complexes Script

This script processes PDB coordinate files to filter encounter complexes based on
z-coordinate criteria. It identifies and removes complexes with negative
z-coordinates.

Features:
    - PDB file parsing and atom extraction
    - Atomic coordinate manipulation
    - Mass and residue property calculation
    - Geometric filtering of encounter complexes
    - Element identification from atom names
    - Supports multi-chain protein complexes

Classes:
    Atom: Represents a single atom with coordinates, properties, and methods
          for extracting data from PDB format lines

Methods:
    - get_attributes(): Extract all atom properties from PDB line
    - get_crds(): Extract 3D coordinates
    - get_mass(): Calculate atomic mass
    - get_element(): Identify chemical element
    - Distance calculations and filtering

Usage:
    python filter_neg_enc.py <input_complexes> [options]

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import argparse
import numpy as np
import math
import copy

class Atom:
    """
    Represents a single atom from a PDB file with properties and parsing methods.
    
    Attributes:
        crds (tuple): Atomic coordinates (x, y, z) in Ångströms
        mass (float): Atomic mass in atomic mass units
        resname (str): Residue name (e.g., 'ALA', 'GLY')
        resid (int): Residue identification number
        id (int): Atom identification number
        name (str): Atom name (e.g., 'CA', 'CB')
        element (str): Chemical element symbol
        chainID (str): Chain identifier (A, B, etc.)
    """

    def __init__(self):
        """Initialize an empty atom with default values."""
        self.crds = []
        self.mass = 0
        self.resname = ""
        self.resid = 0
        self.id = 0
        self.name = ""

    def get_attributes(self, atom):
        """
        Extract all atom attributes from a PDB format line.
        
        Parameters:
            atom (str): PDB format line containing atom information
            
        Side effects:
            Updates all atom properties (coordinates, name, element, mass, etc.)
        """
        self.get_crds(atom)
        self.get_atomname(atom)
        self.get_element()
        self.get_mass()
        self.get_resname(atom)
        self.get_resid(atom)
        self.get_chainID(atom)
        self.get_atomid(atom)

    def get_resname(self, atom):
        """
        Extract residue name from PDB line (columns 17-20).
        
        Parameters:
            atom (str): PDB format line
        """
        resname = atom[17:20].strip()
        self.resname = resname
    
    def get_chainID(self, atom):
        """
        Extract chain identifier from PDB line (column 21).
        
        Parameters:
            atom (str): PDB format line
        """
        chainID = atom[21]
        self.chainID = chainID

    def get_resid(self, atom):
        """
        Extract residue ID from PDB line (columns 22-26).
        
        Parameters:
            atom (str): PDB format line
        """
        resid = int(atom[22:26])
        self.resid = resid

    def get_atomid(self, atom):
        """
        Extract atom ID from PDB line (columns 6-11).
        
        Parameters:
            atom (str): PDB format line
        """
        atomid = int(atom[6:11])
        self.id = atomid

    def get_atomname(self, atom):
        """
        Extract atom name from PDB line (columns 12-16).
        
        Parameters:
            atom (str): PDB format line
        """
        atomname = atom[12:16].strip(" ")
        self.name = atomname

    def get_crds(self, atom):
        """
        Extract 3D coordinates from PDB line.
        
        Parameters:
            atom (str): PDB format line
            
        Note:
            Coordinates are in columns: X (30-38), Y (38-46), Z (46-54)
        """
        x, y, z = float(atom[30:38]), float(atom[38:46]), float(atom[46:54])
        self.crds = x, y, z

    def get_element(self):
        """
        Identify chemical element from atom name.
        
        Recognizes:
            - Single letter elements: C, N, O, H, F, P, S, I
            - Two letter elements: BR, CL, LI, NA, FE, ZN, CU, MN
            
        Sets self.element to 'X' if element cannot be identified.
        """
        one_element_letters = ["C", "N", "O", "H", "F", "P", "S", "I"]
        two_elements_letters = ["BR", "CL", "LI", "NA", "FE", "ZN", "CU", "MN"]
        self.element = ""
        if self.name[0:2].upper() in two_elements_letters:
            self.element = self.name[0:2].upper()
        elif self.name[0] in one_element_letters:
            self.element = self.name[0]
        else:
            print(self.name+": element not identified")
            self.element = "X"


    def set_crds(self, new_crds):
        """
        Update atom coordinates.
        
        Parameters:
            new_crds (tuple or list): New coordinates (x, y, z)
        """
        self.crds = new_crds

    def get_mass(self):
        """
        Assign atomic mass based on element.
        
        Uses standard atomic masses in atomic mass units.
        Prints warning if element is not recognized.
        """
        atom_masses = {"H": 1, 
                       "C": 12, 
                       "N": 14, 
                       "O": 16, 
                       "F": 19, 
                       "P": 31, 
                       "S": 32, 
                       "CL": 35.5, 
                       "BR": 79.9, 
                       "I": 126.9, 
                       "CU": 63.5, 
                       "ZN": 65.4, 
                       "FE": 55.8, 
                       "MN": 54.9, 
                       "LI": 6.9, 
                       "NA": 23}
        if self.element in atom_masses:
            self.mass = atom_masses[self.element]
        else:
            print(self.name+": element not identified")
            print("Assigning mass of 0")

class Residue:
    """
    Represents a single residue containing multiple atoms.
    
    Attributes:
        resname (str): Residue name (e.g., 'ALA', 'GLY')
        id (int): Residue identification number
        atoms (list): List of Atom objects in this residue
        com (tuple): Center of mass coordinates (calculated via get_com)
    """
    
    def __init__(self, resname, id):
        """Initialize a residue."""
        self.resname = resname
        self.atoms = []
        self.id = id

    def add_atoms(self, atoms):
        """Add atoms to this residue."""
        for atom in atoms:
            self.atoms.append(atom)

    def get_com(self):
        """Calculate and store mass-weighted center of mass for this residue."""
        self.com = center_of_mass(self.atoms)


class Chain:
    """
    Represents a protein chain containing residues.
    
    Attributes:
        id (str): Chain identifier (A, B, etc.)
        atoms (list): All atoms in this chain
        residues (dict): Dictionary of Residue objects indexed by residue ID
        com (tuple): Center of mass coordinates (calculated via get_com)
    """
    
    def __init__(self, id):
        """Initialize a chain."""
        self.id = id
        self.atoms = []
        self.residues = {}

    def add_residues(self, residues):
        """Add residues to this chain."""
        for residue in residues:
            self.residues[residue.id] = residue

    def add_atoms(self, atoms):
        """Add atoms to this chain."""
        for atom in atoms:
            self.atoms.append(atom)

    def get_com(self):
        """Calculate and store mass-weighted center of mass for this chain."""
        self.com = center_of_mass(self.atoms)


class PDB:
    """
    Hierarchical PDB structure builder that organizes atoms into residues and chains.
    
    Attributes:
        pdb (list): List of PDB format lines
        atoms (list): All atoms in the structure
        residues (dict): Dictionary of Residue objects indexed by residue ID
        chains (dict): Dictionary of Chain objects indexed by chain ID
        com (tuple): Center of mass of the entire structure
    """
    
    def __init__(self, pdb):
        """Initialize PDB structure."""
        self.pdb = pdb

    def build_hierarchy(self):
        """
        Parse PDB lines and build hierarchical structure of chains, residues, and atoms.
        
        Creates nested organization:
        - Chains contain residues and atoms
        - Residues contain atoms
        - Atoms contain coordinate and property data
        
        Processes ATOM and HETATM records, with TER marking chain termination.
        """
        self.atoms = []
        self.residues = {}
        self.chains = {}

        tmp_residue_atoms = []
        tmp_chain_atoms = []
        tmp_chain_residues = []

        for idx in range(len(self.pdb)):
            pdb_atom = self.pdb[idx]
            if pdb_atom.startswith("ATOM") or pdb_atom.startswith("HETATM"):
                chain_id = self.pdb[idx][21]
                resnum = int(pdb_atom[22:26])

                #Build atom
                atom = Atom()
                atom.get_attributes(pdb_atom)
                self.atoms.append(atom)

                tmp_residue_atoms.append(atom)
                tmp_chain_atoms.append(atom)

                #Build residue
                if resnum not in self.residues:
                    resname = pdb_atom[17:20]
                    self.residues[resnum] = Residue(resname, resnum)
                    if len(self.residues) > 1:
                        keys = list(self.residues.keys())
                        index = keys.index(resnum)
                        previous_resnum = keys[index - 1]
                        self.residues[previous_resnum].add_atoms(tmp_residue_atoms)
                        tmp_chain_residues.append(self.residues[previous_resnum])
                        tmp_residue_atoms = []
                    

                #Build chain
                if chain_id == " ":
                    chain_id = "A"
                if chain_id not in self.chains:
                    self.chains[chain_id] = Chain(chain_id)
                    if len(self.chains) > 1:
                        keys = list(self.chains.keys())
                        index = keys.index(resnum)
                        previous_chain_id = keys[index - 1]
                        self.chains[previous_chain_id].add_atoms(tmp_chain_atoms)
                        self.chains[previous_chain_id].add_residues(tmp_chain_residues)
                        tmp_chain_atoms = []
                        tmp_chain_residues = []

                    
            elif pdb_atom.startswith("TER") or pdb_atom.startswith("END"):

                #Build residue
                self.residues[resnum].add_atoms(tmp_residue_atoms)
                tmp_chain_residues.append(self.residues[resnum])
                tmp_residue_atoms = []

                #Build chain
                self.chains[chain_id].add_atoms(tmp_chain_atoms)
                self.chains[chain_id].add_residues(tmp_chain_residues)
                tmp_chain_atoms = []
                tmp_chain_residues = []

    def get_com(self):
        """Calculate and store mass-weighted center of mass for entire structure."""
        self.com = center_of_mass(self.atoms)


class transformation:
    """
    Applies geometric transformations (translation) to coordinate vectors.
    
    Attributes:
        crds (tuple): 3D coordinates to transform
        vector (tuple): Translation vector [dx, dy, dz]
    """

    def __init__(self, crds):
        """
        Initialize transformation object.
        
        Parameters:
            crds (array-like): 3D coordinates
        """
        self.crds = crds

    def translate(self, translation_vector):
        """
        Translate coordinates by a given vector.
        
        Parameters:
            translation_vector (array-like): Translation vector [dx, dy, dz] in Ångströms
            
        Side effects:
            Updates self.crds with new translated coordinates
        """
        self.vector = translation_vector
        x = self.crds[0]
        y = self.crds[1]
        z = self.crds[2]
        new_x = "{:.3f}".format(x + self.vector[0])
        new_y = "{:.3f}".format(y + self.vector[1])
        new_z = "{:.3f}".format(z + self.vector[2])
        self.crds = np.array([float(new_x), float(new_y), float(new_z)])


class complexes:
    """
    Processes SDA association files and filters encounter complexes.
    
    Attributes:
        encounters (list): List of encounter complex entries from association file
    """

    def __init__(self, assoc_file):
        """
        Initialize complexes processor.
        
        Parameters:
            assoc_file (str): Path to association file
        """
        self.encounters = open(assoc_file, "r").readlines()

    def filter_negative_com(self, xc1, xc2):
        """
        Filter out encounter complexes with negative Z-coordinates.
        
        Parameters:
            xc1 (array): Receptor center of mass coordinates
            xc2 (array): Ligand center of mass coordinates
            
        Side effects:
            Removes lines from self.encounters where transformed ligand COM has negative Z
            
        Note:
            Iterates backwards through encounters to safely delete during iteration
            Skips first 4 lines (header)
        """
        
        for I in range(len(self.encounters)-1, 3, -1):
            lig_com = xc2.copy()
            line = self.encounters[I]
            line_splitted = line.split()
            transx, transy, transz = float(line_splitted[2]), float(line_splitted[3]), float(line_splitted[4])
            trans = np.array([transx, transy, transz])

            new_com = transformation(lig_com)
            new_com.translate(-xc2+xc1+trans)

            if new_com.crds[2] < 0:
                del self.encounters[I]


def center_of_mass(atoms):
    """
    Calculate mass-weighted center of mass for a list of atoms.
    
    Parameters:
        atoms (list): List of Atom objects
        
    Returns:
        np.array: Center of mass coordinates [x, y, z] in Ångströms
        
    Formula:
        COM = Σ(mass_i * coord_i) / Σ(mass_i)
    """
    sum_mass = 0
    sum_massx, sum_massy, sum_massz = 0, 0, 0
    for atom in atoms:

        sum_mass += atom.mass
        sum_massx += atom.crds[0]*atom.mass
        sum_massy += atom.crds[1]*atom.mass
        sum_massz += atom.crds[2]*atom.mass

    massx = sum_massx/sum_mass
    massy = sum_massy/sum_mass
    massz = sum_massz/sum_mass

    massx = "{:.3f}".format(massx)
    massy = "{:.3f}".format(massy)
    massz = "{:.3f}".format(massz)

    com = np.array([float(massx), float(massy), float(massz)])

    return(com)

def write_encounters(outname, encounters):
    """
    Write encounter complexes to an output file.
    
    Parameters:
        outname (str): Output file path
        encounters (list): List of encounter complex lines to write
        
    Side effects:
        Creates output file with encounter complex data
    """

    with open(outname, "w") as fcomplexes:
        for encounter in encounters:
            fcomplexes.write(encounter)

if __name__ == "__main__":
    """
    Main execution block: Parse command-line arguments and filter negative complexes.
    
    Workflow:
        1. Read receptor PDB and extract center of mass
        2. Read ligand PDB and extract center of mass
        3. Read association file with encounter complexes
        4. Filter: remove complexes where transformed ligand COM has negative Z
        5. Write filtered complexes to output file
    """

    parser = argparse.ArgumentParser(description="Filter negative encounter complexes.")
    parser.add_argument("p1", help="Receptor PDB file")
    parser.add_argument("p2", help="Ligand PDB file")
    parser.add_argument("assoc_file", help="Association file")
    parser.add_argument("-o", "--output", dest="output", type=str, default="Min", help="List of residue names to include in output PDB files. Default: all")

    args = parser.parse_args()
    p1 = args.p1
    p2 = args.p2
    assoc_file = args.assoc_file
    output = args.output

    with open(p1, "r") as pdb:
        p1_lines = pdb.readlines()

    with open(p2, "r")  as pdb:
        p2_lines = pdb.readlines()

    #with open(assoc_file, "r") as assoc:
    #    assoc_lines = assoc.readlines()[4:]

    

    prot = PDB(p1_lines)
    prot.build_hierarchy()
    prot.get_com()

    xc1 = prot.com

    lig = PDB(p2_lines)
    lig_pdb_atoms = p2_lines
    lig.build_hierarchy()
    lig.get_com()

    xc2 = lig.com

    encounter_complexes = complexes(assoc_file)
    encounter_complexes.filter_negative_com(xc1, xc2)

    write_encounters(output, encounter_complexes.encounters)
