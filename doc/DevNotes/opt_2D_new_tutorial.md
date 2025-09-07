# opt_2D_new.py - Comprehensive Tutorial

This is a detailed guide for using `opt_2D_new.py`, the command-line utility for visualizing and comparing 2D energy surfaces from XYZ files. This script supports multiple modes and flexible configuration through command-line arguments.

## 1. Overview

`opt_2D_new.py` processes single XYZ files containing 2D energy scan data and performs:
- **plot**: Visualize reference energies from XYZ comments
- **model**: Compute and compare model energies without fitting
- **fit**: Perform parameter fitting and compare results
- **scan**: Scan individual DOF parameters to analyze sensitivity

## 2. Command-Line Arguments

All boolean flags use integer values (1=enable, 0=disable) for flexible default control.

### Core Arguments
- `--mode`: Operating mode (`plot`, `model`, `fit`, `scan`) - default: `fit`
- `--input` / `-i`: Input XYZ file path - default: `/home/prokophapala/Desktop/CARBSIS/wb97m-split/H2O-A1_H2O-D1-y.xyz`
- `--dof-selection`: DOF selection file - default: `dofSelection_MorseSR_H2O.dat`

### Model Parameters
- `--verbosity`: Verbosity level (0=quiet, 1=minimal, 2=normal, 3=verbose) - default: 2
- `--nstep`: Number of fitting steps - default: 1000
- `--fmax`: Target force maximum for convergence - default: 1e-8
- `--dt`: Integrator timestep - default: 0.01
- `--max-step`: Maximum step size - default: 0.05
- `--damping`: Damping coefficient - default: 0.0

### Model Selection
- `--lj`: Use Lennard-Jones instead of Morse potentials (0=Morse, 1=LJ) - default: 0
- `--kMorse`: Global Morse curvature parameter - default: 1.8
- `--Lepairs`: Global electron-pair scaling - default: 1.0

### Weighting Controls
- `--weight-a`: Exponential weight amplitude - default: 1.0
- `--weight-alpha`: Exponential weight sharpness - default: 4.0
- `--emin-min`: Energy threshold for weighting segments - default: -0.02
- `--n-before-min-morse`: Points before minimum for Morse weighting - default: 100
- `--n-before-min-lj`: Points before minimum for LJ weighting - default: 2

### Output Controls
- `--epairs`: Enable/disable electron-pair terms (1=enable, 0=disable) - default: 1
- `--show`: Show figure window (1=show, 0=hide) - default: 1
- `--line`: Plot r_min(angle) and E_min(angle) lines (1=plot, 0=hide) - default: 1
- `--out-xyz`: Output fitted XYZ file (1=output, 0=no output) - default: 0
- `--save`: Path to save plot (PNG format)

### Scan Mode Arguments
- `--scan_dofs`: List of DOF indices to scan (space-separated integers)
- `--scan_range`: Scan range [min, max, n_steps] - default: [-1.0, 1.0, 100]
- `--soft_clamp`: Enable soft clamp during scan (1=enable, 0=disable) - default: 1
- `--user_weights`: Enable user weights during scan (1=enable, 0=disable) - default: 1
- `--regularization`: Enable regularization during scan (1=enable, 0=disable) - default: 1
- `--clamp_thresholds`: Soft clamp thresholds [start, max] - default: [4.0, 6.0]

## 3. Operating Modes

### Plot Mode (`--mode plot`)
Plots reference energies extracted from XYZ comment lines without any model computation.

```bash
python opt_2D_new.py --mode plot --input H2O-A1_H2O-D1-y.xyz --save reference_plot.png
```

### Model Mode (`--mode model`)
Computes model energies using initial parameters without fitting. Useful for baseline comparison.

```bash
python opt_2D_new.py --mode model --input H2O-A1_H2O-D1-y.xyz --lj 0 --show 1
```

### Fit Mode (`--mode fit`) - Default
Performs full parameter fitting using the selected DOF file and optimization parameters.

```bash
python opt_2D_new.py --mode fit --input H2O-A1_H2O-D1-y.xyz --nstep 500 --fmax 1e-6 --save fitted_plot.png
```

### Scan Mode (`--mode scan`)
Scans specified DOF parameters to analyze parameter sensitivity and optimization landscape.

```bash
# Scan first DOF with custom range
python opt_2D_new.py --mode scan --scan_dofs 0 --scan_range -2.0 2.0 50 --save scan_plot.png

# Scan multiple DOFs
python opt_2D_new.py --mode scan --scan_dofs 0 1 2 --scan_range -1.0 1.0 100
```

## 4. Weighting System

The script uses sophisticated weighting to emphasize important regions:

### Exponential Weighting
- `--weight-a`: Controls the amplitude of exponential weighting
- `--weight-alpha`: Controls the sharpness/falloff rate
- `--emin-min`: Energy threshold below which weighting is applied

### Morse vs LJ Weighting
- Morse potentials: `--n-before-min-morse` points before energy minimum get increased weight
- LJ potentials: `--n-before-min-lj` points before minimum get increased weight

Example configurations:
```bash
# Strong weighting near minimum for Morse
python opt_2D_new.py --weight-a 2.0 --weight-alpha 6.0 --n-before-min-morse 150

# Conservative weighting for LJ
python opt_2D_new.py --lj 1 --weight-a 1.0 --weight-alpha 3.0 --n-before-min-lj 3
```

## 5. Soft Clamp and Regularization

### Soft Clamp
Prevents parameters from exploring unphysical regions:
- `--clamp_thresholds`: [start_threshold, max_threshold]
- `--soft_clamp`: Enable soft clamping (1=enable, 0=disable)

### Regularization
Applies penalties to keep parameters near target values:
- `--regularization`: Enable regularization (1=enable, 0=disable)
- Controlled via `dofSelection_*.dat` file with regularization parameters

## 6. Output Options

### Plotting
- `--show`: Control figure window display (1=show, 0=hide)
- `--line`: Control r_min and E_min line plotting (1=plot, 0=hide)
- `--save`: Save plot to PNG file

### XYZ Export
- `--out-xyz`: Export fitted geometries to XYZ file (1=output, 0=no output)

### Verbosity Control
- `--verbosity`: Control amount of console output
  - 0: Silent operation
  - 1: Basic progress information
  - 2: Normal output with DOF values
  - 3: Verbose debugging information

## 7. Advanced Usage Examples

### Complete Fitting Workflow
```bash
python opt_2D_new.py \
  --mode fit \
  --input H2O-A1_H2O-D1-y.xyz \
  --dof-selection dofSelection_MorseSR_H2O.dat \
  --nstep 2000 \
  --fmax 1e-8 \
  --weight-a 1.5 \
  --weight-alpha 4.0 \
  --save final_fit.png \
  --out-xyz 1 \
  --verbosity 2
```

### Parameter Sensitivity Analysis
```bash
# Scan key parameters
python opt_2D_new.py --mode scan --scan_dofs 0 1 2 3 --scan_range -0.5 0.5 50 --save sensitivity_scan.png

# Scan with custom weighting
python opt_2D_new.py --mode scan --scan_dofs 0 --user_weights 0 --save unweighted_scan.png
```

### LJ Model with Custom Parameters
```bash
python opt_2D_new.py \
  --mode fit \
  --lj 1 \
  --n-before-min-lj 5 \
  --kMorse 2.0 \
  --Lepairs 0.8 \
  --clamp_thresholds 3.0 5.0 \
  --save lj_fit.png
```

## 8. Troubleshooting

### Common Issues
- **No convergence**: Increase `--nstep`, decrease `--fmax`, adjust `--dt`
- **Unphysical parameters**: Enable soft clamp, tighten bounds in DOF file
- **Poor weighting**: Adjust `--weight-a`, `--weight-alpha`, `--n-before-min-*`
- **Memory issues**: Reduce `--nstep`, use smaller scan ranges
- **Plot not showing**: Check `--show 1` and matplotlib backend

### Debug Mode
```bash
python opt_2D_new.py --verbosity 3 --mode model  # Verbose output for debugging
```

This tutorial covers all features of `opt_2D_new.py`. For questions about the underlying FitREQ model, see `FitREQ.h` and the main tutorial.
