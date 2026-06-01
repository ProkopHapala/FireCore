# Free Energy Calculation Ecosystem

This directory contains a high-performance framework for calculating free energy profiles using **Thermodynamic Integration (TI)** and (experimentally) the **Jarzynski Equality (JE)**. The architecture supports everything from quick single-system tests to massive grid-search optimizations and temperature sweeps, all accelerated by an OpenCL Molecular Dynamics engine.

---

## 1. Advanced Sweep Control: `free_energy_calculation.sh`

The primary tool for production runs is the sweep runner. It allows you to define a batch of simulations (e.g., a temperature sweep) in a single configuration file.

### Preparing the Configuration
Configs are located in `configs/` (e.g., `default_free_energy_config.json`). The JSON is structured into two main sections:

1.  **`common`**: Global settings inherited by all runs in the sweep.
    *   `xyz_name`: Path to the molecular structure.
    *   `constraints`: Path to the constraint definition file.
    *   `surf_name`: Path to the surface (e.g., `data/xyz/NaCl_3x3_L3`) or `"none"` for vacuum.
    *   `nSys`, `nLambda`, `nMDsteps`, `nEQsteps`: Control the statistical quality and length of the simulation.
    *   `temperature`, `dt`, `mode` (TI, JE, or BOTH).
2.  **`parameter_sets`**: A list of dictionaries. Each entry represents one independent simulation and overrides the `common` block.
    *   *Example:* To run a temperature sweep, keep the `common` block constant and add entries like `{"name": "300K", "temperature": 300.0}` to this list.

### Usage
```bash
./free_energy_calculation.sh configs/your_config.json
```
*If no argument is provided, it defaults to `configs/default_free_energy_config.json`.*

### Analysis & Outputs
*   Each run generates its own subfolder in `results/<output_root>/<run_name>`.
*   An **interactive HTML plot** is generated immediately for each run.
*   Upon completion of the entire sweep, `plot_temperature_sweep.py` is triggered to create:
    *   `comparison_temperatures.html`: An overlay of all free energy curves.
    *   `delta_F_vs_T.csv`: Extracted $\Delta F$ values for further thermodynamic analysis.

---

## 2. Standalone Scripts: `run_DA.sh`, `run_nHex.sh`, etc.

For rapid testing of specific molecules without a complex JSON config, use the pre-configured root scripts.

*   **Molecules:** Donor-Acceptor (`run_DA.sh`), n-Hexane (`run_nHex.sh`), Combined Systems (`run_combined_systems.sh`), and Entropic Spring (`run_ES.sh`).
*   **CLI Customization:** These scripts support powerful command-line overrides:
    *   `--ff`: Switch between `MMFF` and `UFF`.
    *   `--mode`: `TI`, `JE`, or `BOTH`.
    *   `--surf_name`: Provide a path to a surface or `none`.
    *   `-T <temp>`: Set temperature in Kelvin.
    *   `--dt <val>`: Adjust the MD time step.
    *   `--hard_atoms`, `--soft_atoms`, `--hard_dist`, `--soft_dist`: Choose the constraint variant.

**Example:**
```bash
./run_DA.sh --mode TI --ff MMFF --surf_name none -T 300 --soft_atoms
```
The script will build the necessary libraries, run the MD, and automatically open an interactive analysis report.

---

## 3. Parameter Optimization: `run_all_optimizations.sh`

To find the optimal balance between **Accuracy** (low RMSD) and **Performance** (low Wall Time), use the grid-search benchmarking tool.

### Configuration
Adjust parameters in **`configs/optimize_parameters/opt_all.json`**. Unlike standard sweeps, this config accepts **lists** for parameters (e.g., `"dt_list": [0.01, 0.05, 0.1]`). The runner performs a Cartesian product of all provided lists, testing every possible combination.

### Usage & Pareto Analysis
```bash
./run_all_optimizations.sh
```
The script executes hundreds of runs and then triggers `collect_ES_summary.py`. The resulting analysis is found in `results/bench_ES_all/summary/`:
*   `summary_dashboard.html`: **The main visualization.** It plots **Wall Time vs. RMSD**, allowing you to visually identify the "Pareto Front"—the most efficient parameter sets.
*   `best_pareto.csv`: A filtered list of the most efficient runs.

---

## 4. Critical Technical Notes

### Entropic Spring (ES) Benchmarking
The Entropic Spring model (`run_ES.sh`) is used to validate the engine against analytical solutions.
*   **Requirement:** Currently, you **MUST manually disable** non-bonding and angle interactions inside the C++ engine code to match the analytical assumptions of the ideal chain.
*   **TODO:** Automate interaction toggling via Python flags/JSON configuration.

### Known Limitations
*   **Jarzynski Equality (JE):** The JE method is currently experimental and often fails to produce physically meaningful results due to high work variance. **Thermodynamic Integration (TI) is the recommended production method.**
*   **Potential Energy Scans:** Infrastructure for "Rigid" and "Relaxed" scans exists but is currently broken due to atom indexing and constraint issues in the OpenCL backend.
    *   **TODO:** Fix potential scan logic to properly integrate with the multi-system topology. Feel free to test, but expect non-physical results.

### Features
*   **Surfaces:** Can be added to any run via the `--surf_name` flag.
*   **Forcefields:** Easily switchable (`MMFF` vs `UFF`) to compare accuracy/speed.
*   **Temperature:** Fully parametric, enabling easy calculation of entropy and enthalpy from $\Delta F(T)$ curves.
