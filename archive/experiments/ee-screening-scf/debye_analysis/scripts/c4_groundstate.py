#!/usr/bin/env python3
"""
Calculate the ground-state energy of He-like C4+ in a Debye-Huckel plasma.
This script uses the native pfac API. It includes full configuration interaction (CI) up to n=6.
"""

import os
import math
import csv
import pathlib

try:
    from pfac import fac
except ImportError:
    print("Error: The script cannot import pfac.fac. Ensure that the 'fac' conda environment is active.")
    exit(1)

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# --- CONSTANTS ---
NPS = 1.0e-2
RBOHR = 0.52917721067
HARTREE_EV = 27.211386018

def calculate_temperature_ev(lambda_d_a0):
    """
    Calculate the electron temperature in eV.
    
    Parameters:
    lambda_d_a0 : float
        The Debye length in a0.
        
    Returns:
    float
        The electron temperature in eV.
    """
    return lambda_d_a0 * lambda_d_a0 * 4.0 * math.pi * (NPS * RBOHR**3) * HARTREE_EV

# --- CHEN ET AL. (2018) TABLE I FAC REFERENCE DATA ---
CHEN_TABLE_I_FAC = {
    "inf": -32.4771,
    100.0: -32.3768,
    20.0:  -31.9081,
    12.5:  -31.5601,
    10.0:  -31.3298,
    5.0:   -30.1973,
    3.333: -29.0950,
    2.5:   -28.0242,
    2.0:   -26.9818
}

def get_ground_state_energy_au(level_file_path):
    """
    Read the ground-state energy (E0) from the FAC level file.
    
    Parameters:
    level_file_path : str
        The file path to the ASCII level file.
        
    Returns:
    float or None
        The ground-state energy in a.u., or None if the file does not exist.
    """
    if not os.path.exists(level_file_path):
        return None
    with open(level_file_path, 'r') as level_file:
        for line in level_file:
            if line.startswith("E0"):
                energy_ev = float(line.split(',')[1])
                return energy_ev / HARTREE_EV
    return None

def run_fac_groundstate(lambda_d_label, lambda_d_a0, max_principal_quantum_number, output_directory):
    """
    Run the FAC calculation with the pfac API.
    
    Parameters:
    lambda_d_label : str or float
        The text label for the Debye length.
    lambda_d_a0 : float
        The Debye length in a0.
    max_principal_quantum_number : int
        The maximum principal quantum number (n) for the configuration-interaction basis.
    output_directory : pathlib.Path
        The directory to save the FAC output files.
        
    Returns:
    float or None
        The ground-state energy in a.u., or None if the calculation fails.
    """
    run_name = f"C4_ground_L{lambda_d_label}_n{max_principal_quantum_number}".replace(".", "p")
    
    # Clear the FAC memory to prevent configuration overlap.
    fac.ReinitRadial(2)
    fac.ClearOrbitalTable(0)
    
    fac.SetAtom('C')
    
    # A large finite Debye length (1.0e6 a0) represents the free-ion limit (infinity).
    if lambda_d_a0 == float('inf'):
        lambda_d_a0 = 1.0e6
        
    temperature_ev = calculate_temperature_ev(lambda_d_a0)
    
    # Enable Debye screening.
    fac.SetOption('orbital:debye_mode', 1)
    fac.SetOption('radial:ee_screen', 1)
    
    group_list = ["g1"]
    fac.Config('g1', '1s2')
    
    # Build the configuration-interaction basis up to the maximum principal quantum number.
    if max_principal_quantum_number >= 2:
        for n1 in range(1, max_principal_quantum_number + 1):
            configuration_list = []
            for n2 in range(n1, max_principal_quantum_number + 1):
                if n1 == 1 and n2 == 1:
                    continue  # The g1 group already includes this configuration.
                
                if n1 != n2:
                    configuration_list.append(f"{n1}*1 {n2}*1")
                else:
                    configuration_list.append(f"{n1}*2")
                    
            if configuration_list:
                group_name = f"g_n{n1}"
                # Unpack the list into arguments for the C API.
                fac.Config(group_name, *configuration_list)
                group_list.append(group_name)

    fac.PlasmaScreen(4.0, NPS, temperature_ev, 1, 0.0, 1)
    
    fac.ConfigEnergy(0)
    fac.OptimizeRadial(['g1'])
    fac.ConfigEnergy(1)
    
    binary_file_path = str(output_directory / f"{run_name}.lev.b")
    ascii_file_path = str(output_directory / f"{run_name}.lev")
    
    fac.Structure(binary_file_path, group_list)
    fac.MemENTable(binary_file_path)
    fac.PrintTable(binary_file_path, ascii_file_path, 1)
    
    return get_ground_state_energy_au(ascii_file_path)

def main():
    output_directory = pathlib.Path("c4_ground_pfac_results")
    output_directory.mkdir(exist_ok=True)
    
    result_list = []
    
    print("Run the FAC calculation with the native pfac API.")
    print(f"Output directory: {output_directory}/\n")
    print("-" * 85)
    print(f"{'lambda_D (a0)':<13} | {'Calculated E0':<15} | {'Chen E0 (Lit)':<15} | {'Abs Diff':<10} | {'Rel Err %':<10} | {'Max n'}")
    print("-" * 85)

    for lambda_d_label, chen_energy_au in CHEN_TABLE_I_FAC.items():
        if lambda_d_label == "inf":
            lambda_d_a0 = float('inf')
        else:
            lambda_d_a0 = float(lambda_d_label)
        
        calculated_energy_au = None
        best_principal_quantum_number = None
        
        # Start at n=6. If FAC does not converge, decrease n by 1 and start the calculation again.
        for current_n in range(6, 0, -1):
            try:
                energy_au = run_fac_groundstate(lambda_d_label, lambda_d_a0, current_n, output_directory)
                if energy_au is not None:
                    calculated_energy_au = energy_au
                    best_principal_quantum_number = current_n
                    break
            except Exception as exception:
                # The calculation failed. Try the next lower basis size.
                pass
        
        if calculated_energy_au is None:
            print(f"{str(lambda_d_label):<13} | {'FAILED':<15} | {chen_energy_au:<15.4f} | {'-':<10} | {'-':<10} | -")
            continue
            
        absolute_error_au = calculated_energy_au - chen_energy_au
        relative_error_percent = 100.0 * absolute_error_au / abs(chen_energy_au)
        
        print(f"{str(lambda_d_label):<13} | {calculated_energy_au:<15.4f} | {chen_energy_au:<15.4f} | {absolute_error_au:<10.4f} | {relative_error_percent:<10.4f} | n={best_principal_quantum_number}")
            
        result_list.append({
            "lambda_D_label": lambda_d_label,
            "lambda_D_a0": lambda_d_a0,
            "E0_fac_au": calculated_energy_au,
            "chen_fac_au": chen_energy_au,
            "absolute_error_au": absolute_error_au,
            "relative_error_percent": relative_error_percent,
            "max_n_used": best_principal_quantum_number
        })

    print("-" * 85)

    if result_list:
        csv_file_path = output_directory / "c4_ground_comparison.csv"
        with open(csv_file_path, 'w', newline='') as csv_file:
            csv_writer = csv.DictWriter(csv_file, fieldnames=result_list[0].keys())
            csv_writer.writeheader()
            csv_writer.writerows(result_list)
        
        max_absolute_error = max(abs(result_row["absolute_error_au"]) for result_row in result_list)
        max_relative_error = max(abs(result_row["relative_error_percent"]) for result_row in result_list)
        print(f"\nMaximum absolute error: {max_absolute_error:.6f} a.u.")
        print(f"Maximum relative error: {max_relative_error:.6f} %")

        if HAS_MATPLOTLIB:
            plot_results(result_list, output_directory)

def plot_results(result_list, output_directory):
    """
    Plot the calculated ground-state energy against the Debye length.
    """
    finite_result_list = [result_row for result_row in result_list if result_row["lambda_D_label"] != "inf"]
    infinity_result = next((result_row for result_row in result_list if result_row["lambda_D_label"] == "inf"), None)
    
    debye_lengths_a0 = [result_row["lambda_D_a0"] for result_row in finite_result_list]
    calculated_energies_au = [result_row["E0_fac_au"] for result_row in finite_result_list]
    chen_energies_au = [result_row["chen_fac_au"] for result_row in finite_result_list]
    
    plt.figure(figsize=(8, 6))
    plt.plot(debye_lengths_a0, calculated_energies_au, 'b-', label="Calculated energy (pfac API)")
    plt.plot(debye_lengths_a0, chen_energies_au, 'ro', label="Chen et al. (Table I)")
    
    if infinity_result:
        plt.axhline(y=infinity_result["E0_fac_au"], color='b', linestyle='--', label=r"Calculated $\lambda_D = \infty$ limit")
        plt.axhline(y=infinity_result["chen_fac_au"], color='r', linestyle=':', label=r"Chen $\lambda_D = \infty$ limit")

    plt.xscale('log')
    plt.gca().invert_xaxis()
    plt.xlabel(r"Debye Length $\lambda_D$ ($a_0$) [Log Scale]")
    plt.ylabel(r"Ground-State Energy $E_0$ (a.u.)")
    plt.title(r"He-like C$^{4+}$ Ground State vs Debye Length")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    
    plot_file_path = output_directory / "c4_ground_energy_ci.png"
    plt.savefig(plot_file_path, dpi=300)
    plt.close()

if __name__ == "__main__":
    main()
