# make_data.f90

## Description

`make_data` generates encounter-derived numeric arrays and writes them to disk.

The selected mode is controlled by `-data_type`:

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
- `-output <file>`

# Optional arguments

- `-pdb1 <file>` (required for `atoms_dist` and angle types)
- `-atoms1`, `-atoms1a`, `-atoms1b`, `-atoms2b` (as needed by data type)
- `-nb_encounters <N>`
- `-help`

## Modules Used

- `read_input`
- `mod_array`
- `mod_assoc`

## Usage

```bash
./make_data -help
```

```bash
./make_data -data_type z_coord -help
```

```bash
./make_data -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes -data_type z_coord -output array_z.txt
```

## Author

Abraham Muñiz-Chicharro
