"""
Regioselectivity Analysis Script

This script classifies encounter complexes by their regioselectivity, determining
which binding sites on cellulose are preferentially targeted. It analyzes the
distribution of complexes across different reference binding modes.

Features:
    - Encounter complex classification
    - Reference-based binding site identification
    - Regioselectivity quantification
    - Occurrence-weighted statistics
    - Percentage distribution calculation
    - Multi-reference support

Usage:
    python regioselectivity.py -complexes <file1> <file2> ... -references <ref1> <ref2> ...

Arguments:
    -complexes: List of encounter complex association files to classify
    -references: List of reference association files for classification

Output:
    Console output showing:
    - Classification category: number of occurrences / total (percentage)
    - Supports multiple binding site classifications

Example:
    python regioselectivity.py -complexes aligned_complexes -references filtered_C1_complexes filtered_C4_complexes

Implementation:
    1. Parse command-line arguments for complex and reference files
    2. Load all complex entries from input files
    3. For each complex, cross-reference against all reference files
    4. Count occurrences (weighted by occurrence field)
    5. Print classification statistics with percentages

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys
import argparse

def main():
    """
    Main execution function for regioselectivity classification.
    
    Workflow:
        1. Parse command-line arguments
        2. Read encounter complexes from all input files (skip 4-line headers)
        3. Load reference files for cross-reference matching
        4. For each encounter complex:
           - Extract occurrence count (columns 183-193)
           - Match against all reference files
           - Create classification key from matching references
        5. Aggregate and print results with percentages
    """
    parser = argparse.ArgumentParser(description="Classify encounter complexes.")

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

    encounters_to_classify = []

    for assoc_filename in args.complexes:
        encounters_to_classify += open(assoc_filename, "r").readlines()[4:]

    references = {}
    for assoc_filename in args.references:
        reference_lines = open(assoc_filename, "r").readlines()
        references[assoc_filename] = reference_lines

    classification = {}
    total_encounters = 0

    for encounter in encounters_to_classify:
        occurrency = int(float(encounter[183:193]))
        key = ""
        for reference_key in references:
            if encounter in references[reference_key]:
                key += (reference_key + ", ")
        key = key[:-2]
        if key not in classification:
            classification[key] = 0
        classification[key] += occurrency
        total_encounters += occurrency

    for key in classification:
        print(key+":", str(classification[key])+"/"+str(total_encounters), "({:.2f} %)".format(classification[key]/total_encounters*100))

if __name__ == "__main__":
    main()
