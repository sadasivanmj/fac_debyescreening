#!/usr/bin/env python3
"""
Calculate the excitation energy of the 3P1 intercombination line for Be-like ions.
This script generates .sf files and uses the sfac binary to isolate radial solver failures.
"""

import os
import math
import csv
import pathlib
import subprocess

# --- CONSTANTS ---
SFAC_BINARY = os.environ.get("SFAC", "sfac")
NPS = 1.0e-2
RBOHR = 0.52917721067
HARTREE_EV = 27.211386018
CM_INVERSE_PER_HARTREE = 219474.63

# The target ions: (Element Symbol, Atomic Number)
TARGET_IONS = [("C", 6), ("N", 7), ("O", 8), ("Si", 14), ("Fe", 26), ("Mo", 42)]

# The screening parameters (lambda) from the paper x-axis.
SCREENING_PARAMETERS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]

def calculate_temperature_ev(lambda_d_a0):
    """
    Calculate the electron temperature in eV.
    """
    return lambda_d_a0 * lambda_d_a0 * 4.0 * math.pi * (NPS * RBOHR**3) * HARTREE_EV

def get_excitation_energy_au(level_file_path):
    """
    Read the ground-state energy and the 3P1 excited-state energy from the FAC level file.
    """
    if not level_file_path.exists():
        return None
        
    ground_energy_ev = None
    level_list = []
    
    with open(level_file_path, 'r') as level_file:
        for line in level_file:
            if line.startswith("E0"):
                ground_energy_ev = float(line.split(',')[1])
                continue
                
            parts = line.split()
            if len(parts) > 6 and parts[0].isdigit():
                energy_ev = float(parts[2])
                vnl = parts[4]
                two_j = int(parts[5])
                level_list.append({
                    "energy_ev": energy_ev,
                    "vnl": vnl,
                    "two_j": two_j
                })
                
    if ground_energy_ev is None:
        return None
        
    # Find the 1s2 2s1 2p1 states. The outer shell is 2p (vnl == "201").
    # The J = 1 state has 2J = 2.
    target_levels = [level for level in level_list if level["vnl"] == "201" and level["two_j"] == 2]
    
    if not target_levels:
        return None
        
    # The 3P1 state is the lower energy state of the two J=1 states.
    target_levels.sort(key=lambda level: level["energy_ev"])
    excited_energy_ev = target_levels[0]["energy_ev"]
    
    excitation_energy_au = excited_energy_ev / HARTREE_EV
    return excitation_energy_au

def create_sfac_script(atomic_symbol, atomic_number, screening_parameter, ee_screen, run_name):
    """
    Generate the text for the FAC .sf input file.
    """
    # A screening parameter of 0.0 represents the free-ion limit (infinity).
    if screening_parameter == 0.0:
        lambda_d_a0 = 1.0e6
    else:
        lambda_d_a0 = 1.0 / screening_parameter
        
    temperature_ev = calculate_temperature_ev(lambda_d_a0)
    net_charge = float(atomic_number) - 4.0
    
    lines = [
        f"SetAtom('{atomic_symbol}')",
        "SetOption('orbital:debye_mode', 1)",
        f"SetOption('radial:ee_screen', {ee_screen})",
        "Config('g1', '1s2 2s2')",
        "Config('g2', '1s2 2s1 2p1')",
        "Config('g3', '1s2 2p2')",
        f"PlasmaScreen({net_charge}, {NPS}, {temperature_ev}, 1, 0.0, 1)",
        "ConfigEnergy(0)",
        "OptimizeRadial(['g1', 'g2'])",
        "ConfigEnergy(1)",
        f"Structure('{run_name}.lev.b', ['g1', 'g2', 'g3'])",
        f"MemENTable('{run_name}.lev.b')",
        f"PrintTable('{run_name}.lev.b', '{run_name}.lev', 1)"
    ]
    return "\n".join(lines) + "\n"

def run_fac_calculation(atomic_symbol, atomic_number, screening_parameter, ee_screen, output_directory):
    """
    Run the FAC calculation for a Be-like ion using the sfac subprocess.
    """
    run_name = f"{atomic_symbol}_L{screening_parameter}_ee{ee_screen}".replace(".", "p")
    script_content = create_sfac_script(atomic_symbol, atomic_number, screening_parameter, ee_screen, run_name)
    
    script_file_path = output_directory / f"{run_name}.sf"
    ascii_file_path = output_directory / f"{run_name}.lev"
    
    with open(script_file_path, "w") as script_file:
        script_file.write(script_content)
        
    # Execute the sfac binary.
    process_result = subprocess.run(
        [SFAC_BINARY, script_file_path.name], 
        cwd=output_directory, 
        capture_output=True, 
        text=True
    )
    
    if process_result.returncode != 0:
        return None
        
    return get_excitation_energy_au(ascii_file_path)

def main():
    output_directory = pathlib.Path("fritzsche_2006_results")
    output_directory.mkdir(exist_ok=True)
    
    result_list = []
    
    print("Start the FAC calculations for the Be-like ions.")
    print(f"Output directory: {output_directory}/\n")

    for atomic_symbol, atomic_number in TARGET_IONS:
        for screening_parameter in SCREENING_PARAMETERS:
            for ee_screen in [0, 1]:
                
                print(f"Calculate {atomic_symbol} | lambda = {screening_parameter} | ee_screen = {ee_screen} ...")
                
                excitation_energy_au = run_fac_calculation(
                    atomic_symbol, 
                    atomic_number, 
                    screening_parameter, 
                    ee_screen, 
                    output_directory
                )

                if excitation_energy_au is not None:
                    print("Success.")
                    scaled_energy_cm = (excitation_energy_au * CM_INVERSE_PER_HARTREE) / (atomic_number * atomic_number)
                else:
                    print("The orbital is unbound. The calculation failed.")
                    scaled_energy_cm = None
                    
                result_list.append({
                    "atomic_symbol": atomic_symbol,
                    "atomic_number": atomic_number,
                    "screening_parameter": screening_parameter,
                    "ee_screen": ee_screen,
                    "excitation_energy_au": excitation_energy_au,
                    "scaled_energy_cm": scaled_energy_cm
                })

    csv_file_path = output_directory / "fritzsche_2006_data.csv"
    with open(csv_file_path, 'w', newline='') as csv_file:
        csv_writer = csv.DictWriter(csv_file, fieldnames=result_list[0].keys())
        csv_writer.writeheader()
        csv_writer.writerows(result_list)
        
    print(f"\nCalculations are complete. Data saved to {csv_file_path}")

if __name__ == "__main__":
    main()
