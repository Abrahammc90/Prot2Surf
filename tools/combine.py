"""!
@file combine.py
@brief Combine association files and remove duplicate encounter complexes.
@author Abraham Muñiz-Chicharro
@version 1.0
@date 2026-04-05
@par Usage
python combine.py [assoc_file1] [assoc_file2] ... [output_file]

@par Usage flags
- Positional inputs: one or more association files.
- Final positional argument: output filename.

Combine Encounter Complexes Files Script

This script combines multiple SDA encounter complexes files into a single output file,
removing duplicate entries while preserving the header information.

Features:
    - Merges multiple encounter complexes files
    - Removes duplicate complex entries
    - Preserves header information from first file
    - Maintains file format consistency

Usage:
    python combine.py [assoc_file1] [assoc_file2] ... [output_file]
    
    The last argument is treated as the output filename.
    All preceding arguments are treated as input association files.

Example:
    python combine.py assoc_1 assoc_2 assoc_3 combined_output

Implementation:
    1. Reads first file to extract 4-line header
    2. Iterates through all input files, collecting complex entries
    3. Skips duplicate lines using membership testing
    4. Writes combined header and complexes to output file

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys

# Constants
HEADER_LINES = 4

def load_header_from_file(filepath):
    """Load header lines from association file."""
    with open(filepath) as f:
        return f.readlines()[:HEADER_LINES]

def load_complexes_from_file(filepath):
    """Load complex entries from association file (skip header)."""
    with open(filepath) as f:
        return f.readlines()[HEADER_LINES:]

def combine_complexes_removing_duplicates(complex_lists):
    """
    Combine multiple lists of complexes, removing duplicates.
    
    Parameters:
        complex_lists (list): List of lists containing complex lines
        
    Returns:
        list: Combined list with duplicates removed (order-preserving)
    """
    combined = []
    seen = set()
    
    for complex_list in complex_lists:
        for line in complex_list:
            if line not in seen:
                combined.append(line)
                seen.add(line)
    
    return combined

def main():
    """
    Main execution function for combining association files.
    
    Command-line arguments:
        argv[1:-1]: Input association files
        argv[-1]: Output file path
        
    Workflow:
        1. Extract header from first input file (first 4 lines)
        2. Collect all complex entries from all files
        3. Remove duplicates while preserving order
        4. Write combined results to output file
    """
    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    
    # Load header from first file
    header = load_header_from_file(input_files[0])
    
    # Load and combine complexes from all files
    all_complexes = [load_complexes_from_file(filepath) for filepath in input_files]
    combined_complexes = combine_complexes_removing_duplicates(all_complexes)
    
    # Write combined output
    with open(output_file, "w") as output:
        output.writelines(header)
        output.writelines(combined_complexes)

if __name__ == "__main__":
    main()

