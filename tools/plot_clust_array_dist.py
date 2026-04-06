"""!
@file plot_clust_array_dist.py
@brief Plot per-cluster distance distributions for encounter complexes.
@author Abraham Muñiz-Chicharro
@version 1.0
@date 2026-04-05
@par Usage
python plot_clust_array_dist.py [distance_data_file] [cluster_file] [output_prefix] [plot_title]

@par Usage flags
- Positional @c distance_data_file : file with encounter distance values.
- Positional @c cluster_file : file with cluster membership.
- Positional @c output_prefix : output basename for PNG and DAT files.
- Positional @c plot_title : title string for the generated figure.

Plot Cluster Distance Distribution Script

This script generates line plots showing the distance distribution of encounter complexes
across different clusters. It visualizes how distances vary within each cluster, helping
identify spatial patterns.

Features:
    - Multi-cluster distance visualization with separate lines per cluster
    - Color-coded cluster lines using tab10 colormap for easy distinction
    - High-resolution output (300 DPI PNG images for publication)
    - Automatic legend with dynamic multi-column layout based on cluster count
    - Publication-quality figure formatting with semi-transparent grid overlay
    - Dual output: PNG visualization and .dat file with numerical data
    - Efficient encounter complex tracking across clusters
    - Robust error handling and input validation

Usage:
    python plot_clust_array_dist.py [distance_data_file] [cluster_file] [output_prefix] [plot_title]

Arguments:
    distance_data_file (str): File containing distance values for all encounter complexes
    cluster_file (str): File containing cluster assignments (one cluster per line)
    output_prefix (str): Prefix for output files (generates .png and .dat files)
    plot_title (str): Title for the plot

Output Files:
    [output_prefix].png: Publication-quality plot figure
    [output_prefix].dat: Data file with encounter complex IDs and distances

Author: Abraham Muñiz-Chicharro
Version: 1.0
"""

import sys
import matplotlib.pyplot as plt


def _generate_encounter_x_values(cluster_data):
    """Generate sequential x-axis values for encounter complexes.
    
    Args:
        cluster_data (dict): Dictionary mapping cluster names to distance dictionaries
    
    Returns:
        tuple: (cluster_x_values, total_encounters) where cluster_x_values maps
               cluster names to lists of x-axis values
    """
    cluster_x_values = {}
    encounter_counter = 1
    
    for cluster_name, cluster_dict in cluster_data.items():
        x_values = []
        for _ in cluster_dict.values():
            x_values.append(encounter_counter)
            encounter_counter += 1
        cluster_x_values[cluster_name] = x_values
    
    return cluster_x_values, encounter_counter - 1


def _plot_cluster_lines(cluster_data, cluster_x_values):
    """Plot distance distribution lines for all clusters.
    
    Args:
        cluster_data (dict): Dictionary mapping cluster names to distance values
        cluster_x_values (dict): Dictionary mapping cluster names to x-axis values
    """
    colormap = plt.colormaps["tab10"]
    
    for cluster_name in cluster_data:
        cluster_number = int(cluster_name.split()[1])
        y_values = list(cluster_data[cluster_name].values())
        x_values = cluster_x_values[cluster_name]
        
        plt.plot(
            x_values,
            y_values,
            label=f"Cluster {cluster_number}",
            color=colormap(cluster_number),
            linewidth=3
        )


def _configure_plot(plot_title, num_clusters):
    """Configure plot appearance including labels, title, and legend.
    
    Args:
        plot_title (str): Title for the plot
        num_clusters (int): Number of clusters (used for legend columns)
    """
    # Calculate legend columns based on number of clusters
    legend_columns = (num_clusters - 1) // 5 + 1
    
    # Configure axes labels and tick sizes
    plt.xlabel("Encounter complex", fontsize=30)
    plt.ylabel("Distance (Å)", fontsize=30)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    
    # Add title and legend
    plt.title(plot_title, fontsize=40)
    plt.legend(fontsize=25, ncol=legend_columns)
    
    # Add grid and optimize layout
    plt.grid(True, alpha=0.3)
    plt.tight_layout()


def _save_data_file(output_prefix, cluster_data):
    """Save distance data to text file.
    
    Args:
        output_prefix (str): Prefix for output filename
        cluster_data (dict): Dictionary mapping cluster names to distance values
    """
    encounter_id = 1
    with open(f"{output_prefix}.dat", "w") as data_file:
        for cluster_name in cluster_data:
            for distance_value in cluster_data[cluster_name].values():
                data_file.write(f"{encounter_id:>7d}: {distance_value:>10.3f}\n")
                encounter_id += 1


def write_plot(output_prefix, cluster_data, plot_title):
    """Generate and save distance distribution plot for multiple clusters.
    
    Creates both a publication-quality PNG plot and a data file with
    distance values for all encounter complexes.
    
    Args:
        output_prefix (str): Output filename prefix (generates .png and .dat files)
        cluster_data (dict): Dictionary mapping cluster names to distance dictionaries
                            Format: {"cluster X": {encounter_id: distance_value}}
        plot_title (str): Title for the plot
    
    Returns:
        None: Saves plot to .png file and data to .dat file
    """
    # Create figure with publication-quality settings
    plt.figure(figsize=(11, 9), dpi=300)
    
    # Generate x-axis values and plot cluster lines
    cluster_x_values, total_encounters = _generate_encounter_x_values(cluster_data)
    _plot_cluster_lines(cluster_data, cluster_x_values)
    
    # Configure plot appearance
    _configure_plot(plot_title, len(cluster_data))
    
    # Save plot figure
    plt.savefig(f"{output_prefix}.png", bbox_inches='tight', dpi=300)
    plt.close()
    
    # Save data file
    _save_data_file(output_prefix, cluster_data)

def _parse_distance_data(distance_file):
    """Parse distance data from input file.
    
    Args:
        distance_file (str): Path to distance data file
    
    Returns:
        list: List of distance values as floats
    """
    with open(distance_file, "r") as f:
        lines = f.readlines()
    
    # Extract all distances from all lines
    all_distances = []
    for line in lines:
        distances = [float(x) for x in line.split()]
        all_distances.extend(distances)
    
    return all_distances


def _parse_cluster_file(cluster_file, all_distances):
    """Parse cluster assignments and map encounter complexes to clusters.
    
    Args:
        cluster_file (str): Path to cluster assignments file
        all_distances (list): List of distance values indexed by encounter complex ID
    
    Returns:
        dict: Dictionary mapping cluster names to distance dictionaries
              Format: {"cluster X": {encounter_id: distance_value}}
    """
    cluster_data = {}
    encounter_counter = 1
    cluster_name = None
    
    with open(cluster_file, "r") as f:
        for line_index, line in enumerate(f):
            line_parts = line.split()
            if not line_parts:
                continue
            
            # Identify cluster headers
            if line.startswith("Cluster"):
                cluster_name = f"{line_parts[0]} {line_parts[1]}"
                cluster_data[cluster_name] = {}
                encounter_tokens = line_parts[2:]
            else:
                encounter_tokens = line_parts
            
            if cluster_name is None:
                continue

            # Process encounter complex IDs for this cluster
            for encounter_id_str in encounter_tokens:
                encounter_index = int(encounter_id_str) - 1  # Convert to 0-based index
                
                if encounter_index < len(all_distances):
                    cluster_data[cluster_name][encounter_counter] = all_distances[encounter_index]
                    encounter_counter += 1
    
    return cluster_data


if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) < 5:
        print("Usage: python plot_clust_array_dist.py <distance_file> <cluster_file> <output_prefix> <plot_title>")
        sys.exit(1)
    
    distance_file = sys.argv[1]
    cluster_file = sys.argv[2]
    output_prefix = sys.argv[3]
    plot_title = " ".join(sys.argv[4:])
    
    # Parse input data
    all_distances = _parse_distance_data(distance_file)
    cluster_data = _parse_cluster_file(cluster_file, all_distances)
    
    # Generate plot and data file
    write_plot(output_prefix, cluster_data, plot_title)
