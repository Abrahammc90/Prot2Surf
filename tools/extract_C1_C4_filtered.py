"""
Extract and Filter Encounter Complexes Script

This script filters and classifies encounter complexes from SDA simulations based on
reference complex files. It extracts complexes that match reference configurations
for C1 and C4 carbon atoms (critical binding sites in cellulose).

Features:
    - Filters encounter complexes against reference files
    - Generates filtered output files
    - Maintains association file format (4-line header + complexes)
    - Supports multiple reference classifications

Usage:
    python extract_C1_C4_filtered.py -complexes <file1> <file2> ... -references <ref1> <ref2> ...

Arguments:
    -complexes: List of encounter complex association files to filter
    -references: List of reference association files for filtering

Output:
    Creates filtered_<filename> for each reference file containing matching complexes

Implementation:
    1. Parse command-line arguments for complex and reference files
    2. Load all complex entries from complex files (skip 4-line headers)
    3. For each reference file:
        - Extract header information
        - Filter reference complexes that exist in complex files
        - Write filtered results to output file
    4. Output format: filtered_<original_filename>

Example:
    python extract_C1_C4_filtered.py -complexes aligned_complexes -references ref_C1 ref_C4

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys
import argparse
from pathlib import Path

# Constants
HEADER_LINES = 4
OUTPUT_PREFIX = "filtered_"

def load_complexes_from_files(file_list):
    """
    Load all encounter complexes from multiple files, skipping headers.
    
    Parameters:
        file_list (list): List of file paths to read from
        
    Returns:
        set: Set of unique complex entries (excludes header lines)
    """
    complexes = []
    for filepath in file_list:
        with open(filepath, "r") as f:
            lines = f.readlines()
            complexes.extend(lines[HEADER_LINES:])
    return set(complexes)

def extract_filename(filepath):
    """Extract basename from filepath."""
    return Path(filepath).name

def filter_reference_against_complexes(reference_file, encounter_complexes):
    """
    Filter reference complexes against encounter complexes set.
    
    Parameters:
        reference_file (str): Path to reference association file
        encounter_complexes (set): Set of valid encounter complexes
        
    Returns:
        tuple: (filename, header_lines, filtered_complexes)
    """
    with open(reference_file, "r") as f:
        lines = f.readlines()
    
    header = lines[:HEADER_LINES]
    reference_complexes = lines[HEADER_LINES:]
    
    # Filter: keep only complexes that exist in encounter set
    filtered = [complex_line for complex_line in reference_complexes 
                if complex_line in encounter_complexes]
    
    filename = extract_filename(reference_file)
    return filename, header, filtered

def write_filtered_file(output_filename, header_lines, complex_lines):
    """Write filtered complexes to output file."""
    with open(output_filename, "w") as output_file:
        output_file.writelines(header_lines)
        output_file.writelines(complex_lines)

def main():
    """
    Main execution function for filtering encounter complexes.
    
    Workflow:
        1. Parse command-line arguments
        2. Load all encounter complexes from input files (skip headers)
        3. For each reference file:
           - Extract header (first 4 lines)
           - Load all reference complexes
           - Filter: keep only complexes that exist in encounter set
           - Write filtered header + complexes to output file
        4. Output files named: filtered_<basename>
    """
    parser = argparse.ArgumentParser(description="Filter encounter complexes against references.")
    
    parser.add_argument(
        "-complexes", 
        nargs="+", 
        required=True,
        help="List of complex association files"
    )
    parser.add_argument(
        "-references", 
        nargs="+", 
        required=True,
        help="List of reference association files"
    )
    
    args = parser.parse_args()
    
    # Load all encounter complexes into a set for efficient lookup
    encounter_complexes = load_complexes_from_files(args.complexes)
    
    # Process each reference file
    for reference_file in args.references:
        filename, header, filtered_complexes = filter_reference_against_complexes(
            reference_file, 
            encounter_complexes
        )
        
        output_filename = f"{OUTPUT_PREFIX}{filename}"
        write_filtered_file(output_filename, header, filtered_complexes)

if __name__ == "__main__":
    main()