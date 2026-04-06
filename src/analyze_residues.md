# analyze_residues.f90

## Description

`analyze_residues` evaluates residue proximity to a surface across encounter complexes.

For each residue in `pdb2`, it computes a center-of-geometry, applies encounter transformations, and records encounters where the transformed residue position is below a distance threshold.

## Parser description
# Required Arguments

- `-pdb2 <file>`
- `-complexes <file>`

# Optional arguments

- `-nb_encounters <N>` (default: all)
- `-threshold <value>` (default: `6.0`)
- `-help`

## Modules Used

- `read_input`
- `maths`

## Main Routines Called

- `read_pdb` (from `read_input`)
- `read_assoc` (from `read_input`)
- `calculate_cog` (from `maths`)
- `update_complex` (from `maths`)
- `print_help` (internal)

## Output

- `closest_residues.txt`: each residue line includes encounter indexes that satisfy the threshold condition.

## Usage Example

```bash
./analyze_residues -pdb2 p2_noh.pdb -complexes assoc_complexes -nb_encounters 5000 -threshold 6.0
```

## Author

Abraham Muñiz-Chicharro
