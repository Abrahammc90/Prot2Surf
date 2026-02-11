"""
Cluster KON Analysis Script

This script calculates the Smoluchowski rate constant (kon) for diffusing particles
interacting with spherical surfaces, particularly for cellulose protein kinetics analysis.

Features:
    - NAM (Northrup-Allison-McCammon) algorithm implementation for kon calculation
    - SDA input file parsing for diffusion coefficients
    - Cluster beta distribution analysis
    - Spherical surface interaction models
    - Kon rate constant computation with proper unit conversion

Usage:
    python cluster_kon.py <input_file> [options]

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys
import argparse
import numpy as np

def NAM_algorithm(b_surface, q_surface, D, beta):
    """
    Calculate the Smoluchowski rate constant for a diffusing particle to a sphere.
    
    Parameters:
        b_surface (float): radio to b-surface area (Ų)
        q_surface (float): radio to q-surface area (Ų)
        D (float): Diffusion coefficient (Ų/ps)
        beta (float): Inverse of Boltzmann constant * temperature product
        
    Returns:
        float: Rate constant kon (M⁻¹·s⁻¹)
    """
    # Smoluchowski analyticial rate constant
    kb = 4 * np.pi * b_surface * D
    kq = 4 * np.pi * q_surface * D
    
    beta_inf = beta / (1 - (kb / kq))
    kon = kb * beta_inf  # Å^3 / ps
    # Conversion factor: 1 Å^3/(mol·ps) = 6.022e8 L·(mol·s)
    conv_factor = 6.022e8
    kon = kon * conv_factor
    
    return kon

def sda_input_parser(file_path):
    """
    Parse the SDA input file to extract necessary parameters.
    
    Parameters:
        file_path (str): Path to SDA input configuration file
        
    Returns:
        dict: Dictionary containing extracted parameters:
            - 'start_pos' (float): Starting position for simulation
            - 'c' (float): Characteristic parameter c
            - 'nrun' (int): Number of simulation runs
            - 'D' (float): Total diffusion coefficient (Ų/ps)
            
    Note:
        Parser ignores lines starting with '#' (comments) and empty lines.
        Searches for lines starting with parameter names like 'start_pos=', 'c=', etc.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    params = {}
    D = 0.0
    for line in lines:
        if line.strip().startswith('#') or not line.strip():
            continue
        line_stripped = "".join(line.split())
        if line_stripped.startswith('start_pos='):
            params['start_pos'] = float(line.split('=')[1].strip())
        elif line_stripped.startswith('c='):
            params['c'] = float(line.split('=')[1].strip())
        elif line_stripped.startswith('nrun='):
            params['nrun'] = int(line.split('=')[1].strip())
        elif line_stripped.startswith('diffusion_trans'):
            D += float(line.split('=')[1].strip())
    params['D'] = D

    return params

def cluster_beta(cluster_file, total_trajectories):
    """
    Read cluster file and calculate fraction of trajectories contained in clusters.
    
    Parameters:
        cluster_file (str): Path to cluster information file
        total_trajectories (int): Total number of trajectories in simulation
        
    Returns:
        float: Fraction of trajectories that contributed to clusters (0.0 to 1.0)
        
    Note:
        - Reads cluster entries (skips comments starting with '#')
        - Extracts run numbers from first column
        - Removes duplicate run numbers
        - Calculates: len(unique_runs) / total_trajectories
    """
    clusters = []
    with open(cluster_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                clusters.append(line.strip())
    
    trajectory_numbers = []
    for cluster in clusters:
        run_nb = int(cluster.split()[0])
        trajectory_numbers.append(run_nb)

    trajectory_numbers = list(set(trajectory_numbers))

    return len(trajectory_numbers) / total_trajectories
    
    
def main():
    """
    Main execution function for cluster kon analysis.
    
    Workflow:
        1. Parse command-line arguments for SDA input file and cluster files
        2. Extract simulation parameters from SDA input
        3. For each cluster file:
           - Calculate beta (fraction of trajectories in clusters)
           - Calculate kon using NAM_surface formula
        4. Write results to output file
        
    Command-line arguments:
        --sda_input: Path to SDA configuration file
        --clusters: List of cluster file paths
        -o/--output: Output file path for results
    """
    parser = argparse.ArgumentParser(description="Calculate kon for cluster analysis on spherical surfaces")
    parser.add_argument('--sda_input', type=str, dest='sda_input', required=True, help='Path to SDA input file')
    parser.add_argument('--clusters', type=str, dest='clusters', nargs='+', required=True, help='List of cluster names')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True, help='Output file path')

    args = parser.parse_args()
    sda_params = sda_input_parser(args.sda_input)
    b_surface = sda_params['start_pos']
    q_surface = sda_params['c']
    D = sda_params['D']
    total_trajectories = sda_params['nrun']

    with open(args.output, 'w') as out_f:
        out_f.write("# Cluster kon results\n")
        out_f.write(f"# SDA input file: {args.sda_input}\n")
        out_f.write(f"# b surface: {b_surface} Å, q surface: {q_surface} Å, D: {D} Å^2/ps\n")
        out_f.write(f"# Total trajectories: {total_trajectories}\n")
        out_f.write("\n")
        for cluster_file in args.clusters:
            beta = cluster_beta(cluster_file, total_trajectories)
            kon = NAM_algorithm(b_surface, q_surface, D, beta)
            print(f"{cluster_file}: {kon:.2e} M^-1 s^-1")
            out_f.write(f"{cluster_file}: {kon:.2e} M^-1 s^-1\n")

if __name__ == "__main__":
    main()


