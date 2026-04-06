# clust.f90

## Description

`clust` performs hierarchical clustering on encounter data using `min`, `max`, or `mean` linkage.

The input passed with `-input` is auto-detected:

- single line → treated as a 1D array
- multiple lines → treated as a 2D distance matrix

The program reads encounter complexes, runs clustering, and writes output files using the provided output name.

## Parser description
# Required Arguments

- `-complexes <file>`
- `-input <file>`
- `-output_name <name>`

# Optional arguments

- `-nb_encounters <N>` (default: all)
- `-linkage <min|max|mean>` (default: `mean`)
- `-help`

## Modules Used

- `mod_matrix`
- `mod_array`
- `mod_clust_algorithm`
- `read_input`

## Main Routines Called

- `read_assoc` (from `read_input`)
- `read_array` (from `mod_array`)
- `read_matrix` (from `mod_matrix`)
- `linkage_clustering_from_array` (from `mod_clust_algorithm`)
- `linkage_clustering_from_matrix` (from `mod_clust_algorithm`)
- `print_help` (internal)

## Notes

- Clustering runs on CPU for both array and matrix inputs.
- If `-nb_encounters` is omitted, all encounters in the complexes file are used.

## Usage Example

```bash
./clust -complexes assoc_complexes -input matrix_z.txt -nb_encounters 5000 -output_name Cu_z -linkage mean
```

## Author

Abraham Muñiz-Chicharro
