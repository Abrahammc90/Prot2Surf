# clust_all.f90

## Description

`clust_all` builds array/matrix data in memory from structural and encounter inputs, then performs hierarchical clustering in a single workflow.

It is designed to avoid intermediate disk I/O for data construction steps while still writing clustering outputs.

## Parser description
# Required Arguments

- `-pdb2 <file>`
- `-complexes <file>`
- `-atoms2 <atom ...>` (or `-atoms2a`)
- `-data_type <type>`
- `-datadist <file>`

# Optional arguments

- `-pdb1 <file>` (required for `atoms_dist` and angle types)
- `-atoms1`, `-atoms1a`, `-atoms1b`, `-atoms2b` (as needed by data type)
- `-nb_encounters <N>`
- `-linkage <min|max|mean>` (default: `mean`)
- `-output_name <name>`
- `-help`

## Modules Used

- `read_input`
- `mod_matrix`
- `mod_array`
- `mod_clust_algorithm`
- `mod_assoc`

## Workflow Summary

1. Parse CLI options and validate required inputs.
2. Read PDB and complexes files.
3. Build requested data representation according to `-data_type`.
4. Run clustering with selected linkage method.
5. Write clustering outputs using `-output_name`.

## Usage

```bash
./clust_all -help
```

```bash
./clust_all -data_type rmsd -help
```

## Author

Abraham Muñiz-Chicharro
