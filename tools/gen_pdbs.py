"""!
@file gen_pdbs.py
@brief Generate transformed PDB structures from encounter complexes.
@author Abraham Muñiz-Chicharro
@version 1.0
@date 2026-04-05
@par Usage
python gen_pdbs.py [p1.pdb] [p2.pdb] [assoc_file] [-n [N]] [--resnames [RES ...]] [-o [output_prefix]]

@par Usage flags
- Positional @c p1 : receptor PDB file.
- Positional @c p2 : ligand PDB file.
- Positional @c assoc_file : SDA association file.
- @c -n : number of complexes to process (default: all).
- @c --resnames : optional residue-name filter for output structures.
- @c -o, @c --output : output prefix (default: Min).

Generate PDB Files from Complexes Script

This script transforms a ligand according to encounter complex configurations defined
in an SDA association file. It generates individual PDB files for each encounter complex
by applying rotation and translation transformations to the ligand coordinates based on
the association file data.

Features:
    - Transforms ligand coordinates based on association file data
    - Generates individual PDB files for each encounter complex
    - Supports selective residue output filtering
    - Mass-weighted center of mass calculations
    - Rotation and translation operations
    - Multi-chain protein support

Command Line Arguments:
    p1 (positional):           Receptor PDB file path
    p2 (positional):           Ligand PDB file path
    assoc_file (positional):   Association file with encounter complex data
    -n, --num (optional):      Number of encounter complexes to process (default: -1, all)
    --resnames (optional):     List of residue names to include in output (default: all)
    -o, --output (optional):   Output file prefix (default: "Min")

Usage:
    python gen_pdbs.py [receptor_pdb] [ligand_pdb] [assoc_file] [options]

Examples:
    # Generate PDB files for all encounter complexes
    python gen_pdbs.py receptor.pdb ligand.pdb assoc.txt
    
    # Generate only first 100 complexes with specific residues
    python gen_pdbs.py receptor.pdb ligand.pdb assoc.txt -n 100 --resnames GLU ASP
    
    # Custom output prefix
    python gen_pdbs.py receptor.pdb ligand.pdb assoc.txt -o complex

Classes:
    Atom: Represents a single atom with coordinates, properties, and attributes
    Residue: Container for atoms belonging to a single residue
    Chain: Container for residues and atoms of a protein chain
    PDB: Hierarchical PDB structure builder (atoms -> residues -> chains)
    transformation: Applies rotation and translation to atomic coordinates
    complexes: Processes association file and generates encounter complex PDBs

Functions:
    center_of_mass(): Calculate mass-weighted center of mass for atoms
    write_pdbfile(): Write atoms to PDB format file
    get_indices(): Filter atom indices by residue names

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import argparse
import numpy as np
import math
import copy

class Atom:
    """
    Represents a single atom from a PDB file.
    
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
        self.element = ""
        self.chainID = ""

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
        """
        Initialize a residue.
        
        Parameters:
            resname (str): Residue name from PDB
            id (int): Residue identification number
        """
        self.resname = resname
        self.atoms = []
        self.id = id

    def add_atoms(self, atoms):
        """
        Add atoms to this residue.
        
        Parameters:
            atoms (list): List of Atom objects to add
        """
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
        """
        Initialize a chain.
        
        Parameters:
            id (str): Chain identifier from PDB
        """
        self.id = id
        self.atoms = []
        self.residues = {}

    def add_residues(self, residues):
        """
        Add residues to this chain.
        
        Parameters:
            residues (list): List of Residue objects to add
        """
        for residue in residues:
            self.residues[residue.id] = residue

    def add_atoms(self, atoms):
        """
        Add atoms to this chain.
        
        Parameters:
            atoms (list): List of Atom objects to add
        """
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
        """
        Initialize PDB structure.
        
        Parameters:
            pdb (list): List of PDB format lines (typically read from file)
        """
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
    Applies geometric transformations (translation and rotation) to atomic coordinates.
    
    Attributes:
        atoms (list): List of Atom objects to transform
        vector (tuple): Translation vector [dx, dy, dz]
    """

    def __init__(self, atoms):
        """
        Initialize transformation object.
        
        Parameters:
            atoms (list): List of Atom objects to transform
        """
        self.atoms = atoms

    def translate(self, translation_vector):
        """
        Translate all atoms by a given vector.
        
        Parameters:
            translation_vector (array-like): Translation vector [dx, dy, dz] in Ångströms
            
        Side effects:
            Updates coordinates of all atoms in self.atoms
        """
        self.vector = translation_vector
        for I in range(len(self.atoms)):
            x = self.atoms[I].crds[0]
            y = self.atoms[I].crds[1]
            z = self.atoms[I].crds[2]
            new_x = "{:.3f}".format(x + self.vector[0])
            new_y = "{:.3f}".format(y + self.vector[1])
            new_z = "{:.3f}".format(z + self.vector[2])
            self.atoms[I].crds = np.array([float(new_x), float(new_y), float(new_z)])

    def rotate(self, v1, v2):
        """
        Rotate all atoms based on two vectors defining rotation axes.
        
        Parameters:
            v1 (array-like): First reference vector (3D)
            v2 (array-like): Second reference vector (3D)
            
        Returns:
            list: Rotated atoms
            
        Note:
            Constructs rotation matrix from normalized vectors and their cross product.
        """

        def apply_rotation(vector, rotation_matrix):
            """
            Apply rotation matrix to a single vector.
            
            Parameters:
                vector (array): 3D coordinates
                rotation_matrix (array): 3x3 rotation matrix
                
            Returns:
                array: Rotated coordinates
            """
            rotated_vector = np.zeros((3))
            rotated_vector[0] = np.dot(vector, rotation_matrix[0])
            rotated_vector[1] = np.dot(vector, rotation_matrix[1])
            rotated_vector[2] = np.dot(vector, rotation_matrix[2])
            return rotated_vector
        
        idx, idy, idz = 0, 1, 2

        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)
        v3 = np.cross(v1_norm, v2_norm)
        v3_norm = v3 / np.linalg.norm(v2)
        
        #Rotation matrix
        rot_matrix = np.zeros((3, 3))
        rot_matrix[idx] = np.array([v1_norm[idx], v2_norm[idx], v3_norm[idx]])
        rot_matrix[idy] = np.array([v1_norm[idy], v2_norm[idy], v3_norm[idy]])
        rot_matrix[idz] = np.array([v1_norm[idz], v2_norm[idz], v3_norm[idz]])

        
        for I in range(len(self.atoms)):
            self.atoms[I].crds = apply_rotation(self.atoms[I].crds, rot_matrix)

        return self.atoms

class complexes:
    """
    Processes SDA association files and generates transformed PDB structures.
    
    Attributes:
        encounters (list): List of encounter complex entries from association file
    """

    def __init__(self, assoc_file):
        """
        Initialize complexes processor.
        
        Parameters:
            assoc_file (str): Path to association file (skips 4-line header)
        """
        self.encounters = open(assoc_file, "r").readlines()[4:]

    def generate_pdbs(self, n_encounters, indices):
        """
        Generate PDB files for encounter complexes with transformations.
        
        Parameters:
            n_encounters (int): Number of complexes to process (-1 for all)
            indices (list): Atom indices to include in output PDB
            
        Side effects:
            Creates individual PDB files for each encounter complex
            
        Note:
            Reads transformation data (translation + rotation) from association file
            and applies to ligand coordinates relative to receptor center of mass.
        """

        if n_encounters == -1 or n_encounters > len(self.encounters) -4:
            n_encounters = len(self.encounters)

        for I in range(n_encounters):
            outputf = output+"_"+str(I+1)+".pdb"
            line = self.encounters[I]
            line_splitted = line.split()
            transx, transy, transz = float(line_splitted[2]), float(line_splitted[3]), float(line_splitted[4])
            trans = np.array([transx, transy, transz])

            rot1x, rot1y, rot1z = float(line_splitted[5]), float(line_splitted[6]), float(line_splitted[7])
            rot1 =  np.array([rot1x, rot1y, rot1z])
            rot2x, rot2y, rot2z = float(line_splitted[8]), float(line_splitted[9]), float(line_splitted[10])
            rot2 = np.array([rot2x, rot2y, rot2z])

            atoms_copy = copy.deepcopy(lig.atoms)

            transformable_atoms = transformation(atoms_copy)
            transformable_atoms.translate(-xc2)

            transformable_atoms.rotate(rot1, rot2)
            transformable_atoms.translate(xc1+trans)

            write_pdbfile(outputf, atoms_copy, indices)


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

def write_pdbfile(outname, atoms, indices):
    """
    Write selected atoms to a PDB format file.
    
    Parameters:
        outname (str): Output file path
        atoms (list): List of Atom objects
        indices (list): Indices of atoms to write (filtered atoms)
        
    Side effects:
        Creates PDB file with ATOM records for selected atoms
        
    Note:
        Output format follows standard PDB specification with:
        - Occupancy: 1.00
        - B-factor: 0.00
        - Element symbol in the last column
    """

    with open(outname, "w") as pdb:
        for idx in indices:
            atom = atoms[idx]
            line = "{:<6}{:>5}  {:<4}{:>3} {:1}{:>4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2}".format(
                "ATOM",
                atom.id,
                atom.name,
                atom.resname,
                atom.chainID,
                atom.resid,
                atom.crds[0],
                atom.crds[1],
                atom.crds[2],
                1.00,   # occupancy
                0.00,   # b-factor
                atom.element
            )
            pdb.write(line + "\n")
        pdb.write("END\n")

def get_indices(atoms, list):
    """
    Filter atom indices based on residue name criteria.
    
    Parameters:
        atoms (list): List of Atom objects
        list (list): List of residue names to include (empty list = all atoms)
        
    Returns:
        list: Indices of atoms matching the filter criteria
        
    Note:
        Residue name matching is case-insensitive (converted to uppercase)
    """
    indices = []
    for idx in range(len(atoms)):
        if len(list) == 0:
            indices.append(idx)
        elif atoms[idx].resname.upper() in list:
            indices.append(idx)
    return indices

if __name__ == "__main__":
    """
    Main execution block: Parse command-line arguments and generate encounter complex PDB files.
    
    Workflow:
        1. Read receptor PDB and extract center of mass
        2. Read ligand PDB and extract center of mass
        3. Read association file with transformation data
        4. For each encounter complex:
           - Apply transformation (translate to origin, rotate, translate to receptor)
           - Write transformed coordinates to output PDB file
    """

    parser = argparse.ArgumentParser(description="Transform ligand in encounter complex according to association file.")
    parser.add_argument("p1", help="Receptor PDB file")
    parser.add_argument("p2", help="Ligand PDB file")
    parser.add_argument("assoc_file", help="Association file")
    parser.add_argument("-n", dest="n", type=int, default=-1, help="Number of encounter complexes to process. Default: all (-1)")
    parser.add_argument("--resnames", dest="resnames", nargs="+", default=[], help="List of residue names to include in output PDB files. Default: all")
    parser.add_argument("-o", "--output", dest="output", type=str, default="Min", help="List of residue names to include in output PDB files. Default: all")

    args = parser.parse_args()
    p1 = args.p1
    p2 = args.p2
    assoc_file = args.assoc_file
    n_encounters = args.n
    resname_list = [resname.upper() for resname in args.resnames]
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

    indices = get_indices(lig.atoms, resname_list)

    encounter_complexes = complexes(assoc_file)
    encounter_complexes.generate_pdbs(n_encounters, indices)
