# make_data.f90

## Description

`make_data` generates encounter-derived numeric data and writes it to disk.

The selected mode is controlled by `-data_type`:

- `rmsd`
- `z_coord`
- `atoms_dist`
- `2D_angle`
- `3D_angle`

## Parser description
# Required Arguments

- `-pdb2 <file>`
- `-atoms2 <atom ...>` (or `-atoms2a`)
- `-complexes <file>`
- `-data_type <type>`
- `-input <file>`

# Optional arguments

- `-pdb1 <file>` (required for `atoms_dist` and angle types)
- `-atoms1`, `-atoms1a`, `-atoms1b`, `-atoms2b` (as needed by data type)
- `-nb_encounters <N>`
- `-help`

## Modules Used

- `read_input`
- `mod_matrix`
- `mod_array`
- `mod_assoc`

## Usage

```bash
./make_data -help
```

```bash
./make_data -data_type rmsd -help
```

## Author

Abraham Muñiz-Chicharro
