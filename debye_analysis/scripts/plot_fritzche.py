#!/usr/bin/env python3
"""
Plot the excitation energies for Be-like ions to reproduce 
Figure 1 from Saha & Fritzsche (2006).
"""

import csv
import pathlib

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: The script cannot import matplotlib. Ensure it is installed.")
    exit(1)

def main():
    # Define file paths
    output_directory = pathlib.Path("fritzsche_2006_results")
    csv_file_path = output_directory / "fritzsche_2006_data.csv"
    
    if not csv_file_path.exists():
        print(f"Error: The file {csv_file_path} does not exist.")
        exit(1)
        
    result_list = []
    
    # Read and parse the CSV data
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            # Skip empty entries (unbound states due to pressure ionization)
            if row["scaled_energy_cm"].strip():
                result_list.append({
                    "atomic_symbol": row["atomic_symbol"],
                    "screening_parameter": float(row["screening_parameter"]),
                    "ee_screen": int(row["ee_screen"]),
                    "scaled_energy_cm": float(row["scaled_energy_cm"])
                })

    # Create the two-panel figure
    # We use a slightly taller figure to match the paper's aspect ratio
    figure, (axis_left, axis_right) = plt.subplots(1, 2, figsize=(10, 8))
    
    # Define marker styles mapping to closely match the paper
    # Format: (marker for ee=0, marker for ee=1)
    style_map = {
        "C":  ('s', 'v'), # Square / Down-triangle
        "N":  ('o', 'D'), # Circle / Diamond
        "O":  ('^', '>'), # Up-triangle / Right-triangle
        "Si": ('s', 'v'), # Square / Down-triangle
        "Fe": ('o', 'D'), # Circle / Diamond
        "Mo": ('^', '>')  # Up-triangle / Right-triangle
    }
    
    # Plot the light ions on the left axis
    plot_ion(result_list, "C", axis_left, "C III", style_map["C"])
    plot_ion(result_list, "N", axis_left, "N IV", style_map["N"])
    plot_ion(result_list, "O", axis_left, "O V", style_map["O"])
    
    # Plot the heavy ions on the right axis
    plot_ion(result_list, "Si", axis_right, "Si XI", style_map["Si"])
    plot_ion(result_list, "Fe", axis_right, "Fe XXIII", style_map["Fe"])
    plot_ion(result_list, "Mo", axis_right, "Mo XXXIX", style_map["Mo"])
    
    # Format the left axis
    axis_left.set_xlabel(r"$\lambda$ (a.u.)", fontsize=12)
    axis_left.set_ylabel(r"$\Delta E / Z^2$ (cm$^{-1}$)", fontsize=12)
    axis_left.set_xlim(0.0, 0.45)
    # The absolute limits are slightly adjusted to frame FAC's local exchange values
    axis_left.set_ylim(1275, 1550) 
    axis_left.tick_params(direction="in", top=True, right=True)
    
    # Format the right axis
    axis_right.set_xlabel(r"$\lambda$ (a.u.)", fontsize=12)
    axis_right.set_ylabel(r"$\Delta E / Z^2$ (cm$^{-1}$)", fontsize=12)
    axis_right.set_xlim(0.0, 0.55)
    axis_right.set_ylim(400, 1000)
    axis_right.tick_params(direction="in", top=True, right=True)
    
    plt.tight_layout()
    
    # Save the output
    plot_file_path = output_directory / "fritzsche_2006_fig1_reproduced.png"
    plt.savefig(plot_file_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot successfully generated and saved to {plot_file_path}")

def plot_ion(result_list, atomic_symbol, axis, label, markers):
    """
    Filter the data for a specific ion and plot the two screening models.
    """
    marker_solid, marker_dashed = markers

    # Extract electron-nucleus only data (solid line)
    en_only_data = [row for row in result_list if row["atomic_symbol"] == atomic_symbol and row["ee_screen"] == 0]
    en_only_data.sort(key=lambda row: row["screening_parameter"])
    
    # Extract full Debye data (dashed line)
    full_debye_data = [row for row in result_list if row["atomic_symbol"] == atomic_symbol and row["ee_screen"] == 1]
    full_debye_data.sort(key=lambda row: row["screening_parameter"])
    
    x_en_only = [row["screening_parameter"] for row in en_only_data]
    y_en_only = [row["scaled_energy_cm"] for row in en_only_data]
    
    x_full_debye = [row["screening_parameter"] for row in full_debye_data]
    y_full_debye = [row["scaled_energy_cm"] for row in full_debye_data]
    
    # Plot the lines
    axis.plot(x_en_only, y_en_only, marker=marker_solid, linestyle='-', color='black', markersize=5, linewidth=1.2)
    axis.plot(x_full_debye, y_full_debye, marker=marker_dashed, linestyle='--', color='black', markersize=5, linewidth=1.2)
    
    # Add the text label near the end of the solid line
    if x_en_only and y_en_only:
        final_x = x_en_only[-1]
        final_y = y_en_only[-1]
        axis.text(final_x + 0.02, final_y, label, verticalalignment='center', fontsize=10)

if __name__ == "__main__":
    main()
