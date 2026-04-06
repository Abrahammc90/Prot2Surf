# threshold.f90

## Description

`threshold` selects encounter complexes based on a threshold criterion and writes both the computed value array and the filtered complexes.

Supported threshold modes include:

- `z_coord`
- `atoms_dist`
- `2D_angle`
- `3D_angle`

Type-specific help can be shown with:

```bash
./threshold -threshold_type <type> -help
```

## Parser description
# Required Arguments

- `-pdb2 <file>`
- `-atoms2 <atom ...>` (or `-atoms2a`)
- `-array <file>`
- `-complexes <file>`
- `-complexes_output <file>`
- `-threshold_type <type>`
- `-cutoff <value>`

# Optional arguments

- `-pdb1 <file>`
- `-atoms1`, `-atoms1a`, `-atoms1b`, `-atoms2b`
- `-nb_encounters <N>`
- `-help`

## Modules Used

- `read_input`
- `mod_threshold`

## Main Routines Called

- `read_pdb` (from `read_input`)
- `read_assoc` (from `read_input`)
- `array_z_coord`, `array_angle`, `array_atoms_dist` (from `mod_threshold`)
- `sort_array` (from `mod_threshold`)
- `write_array`, `write_cutoff_complexes` (internal)
- `print_help` (internal)

## Usage Example

```bash
./threshold -pdb2 p2_noh.pdb -atoms2 Cu -complexes assoc_complexes -threshold_type z_coord -complexes_output threshold_z.txt -array z_values.txt -cutoff 5.0
```

```bash
./threshold -threshold_type atoms_dist -help
```

## Author

Abraham Muñiz-Chicharro
