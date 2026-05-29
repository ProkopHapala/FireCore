#!/usr/bin/env python3

import os
import glob
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt

def detect_fe_column(filename: str) -> int:
    """
    Detects the column index for free energy ('TI_F' or 'cumulative_FE').
    Returns the column index or -1 if not found.
    """
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("#"):
                # Stop if we are past the header comments
                break
            # Look for a header line with column names
            if "lambda" in line and ("TI_F" in line or "cumulative_FE" in line):
                parts = line.replace("#", "").split()
                for i, p in enumerate(parts):
                    if p in ("TI_F", "cumulative_FE"):
                        return i
    return -1

def main():
    parser = argparse.ArgumentParser(description="Plot difference of max(F) - min(F) vs Temperature")
    parser.add_argument("--results_dir", default=".", help="Directory containing the temperature subdirectories (e.g., 'results_sweep')")
    args = parser.parse_args()

    results_dir = args.results_dir
    
    # Recursively find all .dat files in the directory
    files = glob.glob(os.path.join(results_dir, "**", "*.dat"), recursive=True)
    
    if not files:
        print(f"No .dat files found recursively in '{results_dir}'.")
        return

    temperatures = []
    delta_Fs = []
    
    print(f"Found {len(files)} '.dat' files. Processing...")

    for f in files:
        # Extract temperature from the directory name (e.g., matching "temp_300K")
        m = re.search(r'temp_(\d+(?:\.\d+)?)K', f)
        if not m:
            print(f"Skipping {f}: could not determine temperature from path.")
            continue
        
        T = float(m.group(1))

        try:
            fe_col_idx = detect_fe_column(f)
            if fe_col_idx == -1:
                print(f"Warning: Could not detect a Free Energy column ('TI_F' or 'cumulative_FE') in the header of '{f}'. Skipping.")
                continue

            data = np.loadtxt(f, comments="#")
            
            if data.size == 0:
                continue
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            if fe_col_idx >= data.shape[1]:
                print(f"Skipping {f}: detected FE column {fe_col_idx} is out of bounds for data shape {data.shape}.")
                continue
                
            # Use the detected column for Free Energy
            F = data[:, fe_col_idx]
            delta_F = np.max(F) - np.min(F)
            
            temperatures.append(T)
            delta_Fs.append(delta_F)
            
        except Exception as e:
            print(f"Skipping {f} due to error: {e}")
            
    if not temperatures:
        print(f"No valid data was processed from the found files in '{results_dir}'.")
        return
        
    # Sort the values by temperature for sequential plotting
    sort_idx = np.argsort(temperatures)
    temperatures = np.array(temperatures)[sort_idx]
    delta_Fs = np.array(delta_Fs)[sort_idx]
    
    print(f"Plotting data for {len(temperatures)} temperatures.")

    # Plotting
    plt.figure(figsize=(8, 5))
    plt.plot(temperatures, delta_Fs, marker='o', linestyle='-', linewidth=2, color='#2196F3')
    # plt.xscale('log')
    plt.xlabel('Temperature [K]', fontsize=12)
    plt.ylabel(r'$\Delta F = \max(F) - \min(F)$ [eV]', fontsize=12)
    plt.title('Maximum Free Energy Difference vs Temperature', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Save the plot
    out_path = os.path.join(results_dir, "delta_F_vs_T.png")
    plt.savefig(out_path, dpi=150)
    print(f"Plot saved to {out_path}")
    
    # Display the plot if not in a headless environment
    if 'DISPLAY' in os.environ or os.name == 'nt':
        try:
            plt.show()
        except Exception as e:
            print(f"Could not display plot interactively: {e}")
    else:
        print("No display found. Plot was saved to file.")


if __name__ == "__main__":
    main()
