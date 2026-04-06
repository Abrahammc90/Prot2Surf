"""!
@file extract_grid.py
@brief Parse and extract UHBD electrostatic grid data.
@author Abraham Muñiz-Chicharro
@version 1.0
@date 2026-04-05
@par Usage
python extract_grid.py [input.uhbd] [Nx] [Ny] [Nz] [output.uhbd]

@par Usage flags
- Positional @c input.uhbd : source UHBD grid file.
- Positional @c Nx @c Ny @c Nz : target output dimensions.
- Positional @c output.uhbd : destination UHBD file.

Extract Grid Data from UHBD Files Script

This script loads and extracts electrostatic potential grid data from UHBD format files. 
It's used for processing electrostatic grid data generated with APBS.

Features:
    - Parses binary UHBD format files
    - Extracts grid dimensions, spacing, and origin
    - Reads electrostatic potential values
    - Returns 3D numpy array with proper reshaping
    - Handles UTF-8 binary data decoding
    - Refactored for improved readability and maintainability

Usage:
    python extract_grid.py [input.uhbd] [Nx] [Ny] [Nz] [output.uhbd]

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import numpy as np
import sys

# UHBD Format Constants
HEADER_LINE_COUNT = 2
GRID_VALUE_PER_LINE = 6
FLOAT_FIELD_WIDTH = 13


def parse_uhbd_header(header_line):
    """
    Parse UHBD header line containing grid parameters.
    
    UHBD format uses fixed-width fields:
    - Columns 0-6: Nx
    - Columns 7-13: Ny
    - Columns 14-20: Nz
    - Columns 21-32: delta
    - Columns 33-44: origin_x
    - Columns 45-56: origin_y
    - Columns 57-68: origin_z
    
    Returns:
        dict: Contains keys 'Nx', 'Ny', 'Nz', 'delta', 'origin_x', 'origin_y', 'origin_z'
    """
    try:
        return {
            'Nx': int(float(header_line[0:7])),
            'Ny': int(float(header_line[7:14])),
            'Nz': int(float(header_line[14:21])),
            'delta': float(header_line[21:33]),
            'origin_x': float(header_line[33:45]),
            'origin_y': float(header_line[45:57]),
            'origin_z': float(header_line[57:69])
        }
    except ValueError:
        parts = header_line.split()
        if len(parts) < 7:
            raise
        return {
            'Nx': int(float(parts[0])),
            'Ny': int(float(parts[1])),
            'Nz': int(float(parts[2])),
            'delta': float(parts[3]),
            'origin_x': float(parts[4]),
            'origin_y': float(parts[5]),
            'origin_z': float(parts[6])
        }


def is_coordinate_header(line_parts):
    """Check if line is a coordinate header (exactly three integers)."""
    if len(line_parts) != 3:
        return False
    return all(part.lstrip('-').isdigit() for part in line_parts)


def read_uhbd_file(filename):
    """
    Load UHBD file and extract header and grid data.
    
    Parameters:
        filename (str): Path to UHBD format file
    
    Returns:
        tuple: (header_lines, grid_data, grid_params) where:
            - header_lines (list): File header information
            - grid_data (np.ndarray): 3D array of grid values reshaped as (Nz, Nx, Ny)
            - grid_params (dict): Grid metadata (Nx, Ny, Nz, delta, origin coordinates)
    """
    with open(filename, "rb") as f:
        header = []
        
        # Read initial header lines
        for _ in range(HEADER_LINE_COUNT):
            line = f.readline()
            header.append(line.decode('utf-8'))
        
        # Read grid parameters from third header line
        header_line = f.readline().decode('utf-8')
        header.append(header_line)
        grid_params = parse_uhbd_header(header_line)
        
        # Read remaining header lines
        for _ in range(2):
            line = f.readline()
            header.append(line.decode('utf-8'))
        
        # Read grid data, filtering out coordinate headers
        data = []
        for line in f.readlines():
            decoded_line = line.decode('utf-8').strip()
            parts = decoded_line.split()
            
            # Skip coordinate header lines (z-slice indicators)
            if not is_coordinate_header(parts):
                data.extend(parts)
        
        # Reshape data into 3D array
        data_array = np.array(data, dtype=np.float32)
        Nz = grid_params['Nz']
        Nx = grid_params['Nx']
        Ny = grid_params['Ny']
        data_array = data_array.reshape((Nz, Nx, Ny))
    
    return header, data_array, grid_params


def extract_center_region(grid_data, current_dims, target_dims):
    """
    Extract central region from grid with new dimensions.
    
    Parameters:
        grid_data (np.ndarray): Original 3D grid
        current_dims (tuple): Original dimensions (Nx, Ny, Nz)
        target_dims (tuple): Target dimensions (Nx, Ny, Nz)
        
    Returns:
        tuple: (extracted_grid, offsets) where offsets are [dx, dy, dz]
    """
    Nx, Ny, Nz = current_dims
    new_Nx, new_Ny, new_Nz = target_dims
    
    # Calculate symmetric offsets from grid edges
    dx = (Nx - new_Nx) // 2
    dy = (Ny - new_Ny) // 2
    dz = (Nz - new_Nz) // 2
    
    # Extract center region
    extracted = grid_data[
        dz:dz + new_Nz,
        dx:dx + new_Nx,
        dy:dy + new_Ny
    ]
    
    return extracted, [dx, dy, dz]


def format_grid_value(value):
    """Format single grid value in UHBD scientific notation."""
    return "{:13.5e}".format(value)


def write_coordinate_header(grd_file, z_index, Nx, Ny):
    """Write z-slice coordinate header line."""
    header_str = "{:7s}{:7s}{:7s}\n".format(str(z_index), str(Nx), str(Ny))
    grd_file.write(header_str)


def write_grid_slice(grd_file, grid_slice):
    """Write grid data for one z-slice (6 values per line in UHBD format)."""
    flat_data = grid_slice.flatten()
    
    for i in range(0, len(flat_data), GRID_VALUE_PER_LINE):
        row = flat_data[i:i + GRID_VALUE_PER_LINE]
        formatted_row = "".join(format_grid_value(val) for val in row)
        grd_file.write(formatted_row + "\n")


def update_header_for_new_grid(header, original_dims, new_dims, offsets, delta):
    """
    Update header lines with new grid dimensions and adjusted origin.
    
    Parameters:
        header (list): Original header lines (modified in-place)
        original_dims (tuple): Original (Nx, Ny, Nz)
        new_dims (tuple): New (Nx, Ny, Nz)
        offsets (list): [dx, dy, dz] pixel offsets from original origin
        delta (float): Grid spacing (unchanged)
    """
    new_Nx, new_Ny, new_Nz = new_dims
    dx, dy, dz = offsets
    
    # Update line 1 with new Z dimension
    h1_parts = header[1].split()
    header[1] = header[1][:38] + "{:7s}{:7s}{:7s}\n".format(
        str(new_Nz), h1_parts[1], str(new_Nz))
    
    # Update line 2: new dimensions + adjusted origin
    h2 = header[2]
    origin_x = float(h2[33:45])
    origin_y = float(h2[45:57])
    origin_z = float(h2[57:69])
    
    # Adjust origin based on pixel offset
    new_origin_x = origin_x + dx * delta
    new_origin_y = origin_y + dy * delta
    new_origin_z = origin_z + dz * delta
    
    header[2] = "{:7s}{:7s}{:7s}{:12.5e}{:12.5e}{:12.5e}{:12.5e}\n".format(
        str(new_Nx), str(new_Ny), str(new_Nz), delta, 
        new_origin_x, new_origin_y, new_origin_z)


def write_uhbd_file(filename, header, grid_data, new_dims):
    """
    Write grid data to UHBD format file.
    
    Parameters:
        filename (str): Output file path
        header (list): Updated header lines
        grid_data (np.ndarray): 3D array of grid data
        new_dims (tuple): Grid dimensions (Nx, Ny, Nz)
    """
    new_Nx, new_Ny, new_Nz = new_dims
    
    with open(filename, "w") as grd:
        # Write header
        for line in header:
            grd.write(line)
        
        # Write each z-slice with coordinate header and data
        for z_index in range(new_Nz):
            write_coordinate_header(grd, z_index + 1, new_Nx, new_Ny)
            write_grid_slice(grd, grid_data[z_index])


def main():
    """
    Extract and resize grid from UHBD file.
    
    Command-line arguments:
        argv[1]: Input UHBD file path
        argv[2]: New grid dimension Nx
        argv[3]: New grid dimension Ny
        argv[4]: New grid dimension Nz
        argv[5]: Output file path
        
    Workflow:
        1. Load UHBD grid file
        2. Extract central region with new dimensions
        3. Update header with new dimensions and origin
        4. Write to output file in UHBD format
    """
    input_filename = sys.argv[1]
    target_dims = [int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])]
    output_filename = sys.argv[5]
    
    # Load original grid
    header, grid_data, grid_params = read_uhbd_file(input_filename)
    current_dims = [grid_params['Nx'], grid_params['Ny'], grid_params['Nz']]
    
    # Extract center region
    extracted_grid, offsets = extract_center_region(
        grid_data, current_dims, target_dims
    )
    
    # Update header with new parameters
    update_header_for_new_grid(
        header, current_dims, target_dims, offsets, grid_params['delta']
    )
    
    # Write output file
    write_uhbd_file(output_filename, header, extracted_grid, target_dims)


if __name__ == "__main__":
    main()
