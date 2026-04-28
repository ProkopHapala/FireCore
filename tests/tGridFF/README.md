# Conformational Sampling

This program is designed for **sampling molecular conformations on ionic crystal surfaces** with performance analysis of various computational parameters.



## Opening the GUI

```bash
python3 throughput_gui.py
```

## Purpose

This tool enables:
- **Sampling of molecular conformations** on ionic crystal surfaces (e.g., NaCl)
- **Optimization of computational parameters** for molecular dynamics
- **Performance analysis** of different simulation settings
- **Data generation** for studying molecule-surface interactions
## Main GUI Functions

### 1. "Optimization" Tab - Main Simulations

#### Basic Parameters
- **dovdW**: Van der Waals interactions (0/1)
- **doSurfAtoms**: Use surface atoms (0/1)
- **bGridFF**: Grid Force Field (0=disabled, 1=linear, 6=bSpline)
- **bTex**: Texture (0/1)
- **bSaveToDatabase**: Save to database (0/1) (calculation is much slower, mostly without surface)

#### File Paths
- **XYZ File**: Molecule file (e.g., `data/xyz/xylitol_WO_gridFF`)
- **Surface File**: Template for surface files with `${N}` placeholder
  - Example: `data/xyz/surfaces_for_throughput/NaCl_${N}x${N}_L3`
  - `${N}` is automatically replaced with surface size (1x1, 4x4, 8x8, etc.)

#### Convergence Criteria
- **Fconv**: Force convergence criterion (e.g., 1e-4)

#### Advanced Parameters
- **Replicas**: Number of replicas (e.g., "1000,5000")
- **Perframes**: Per frame parameter (e.g., "20,500")
- **PerVF**: Per VF parameter (e.g., "20,50")
- **nPBC**: Periodic boundary conditions (e.g., "(1,1,0)")
- **Ns**: Surface sizes (e.g., "1-16" or "1,4,8,16")

#### Local Memory Parameters
- **nlocMMFFs**: Local memory for MMFF
- **nlocmoves**: Local memory for moves
- **nlocNBFFs**: Local memory for NBFF
- **nlocSurfs**: Local memory for surface
- **nlocGridFFs**: Local memory for GridFF
- **nlocGridFFbSplines**: Local memory for GridFF bSplines

### 2. "Visualization" Tab - Results Analysis

- **Generate interactive performance plots**
- **View results table**
- **Open results folder**
- **Auto-open in browser**

## Usage

### 1. Loading Presets
```
File → Load Preset → Select file (e.g., params_surface_with_gridff.txt)
```

### 2. Manual Setup
1. Set basic parameters according to simulation type
2. Enter paths to XYZ and surface files
3. Set parameter ranges for testing
4. Run optimization

### 3. Results Analysis
1. Switch to "Visualization" tab
2. Select results file
3. Generate interactive plots
4. Analyze performance data

## Output Data

The program generates:
- **minima*.dat**: Individual simulation results
- **results.dat**: Summary results of all runs
- **results_table.csv**: Table for analysis
- **Interactive HTML plots**: For performance visualization


Individual simulation results could be studied if copied to `number_evaluations_vs_time/` folder and run evaluation_vs_time.py script `python3 evaluation_vs_time.py`
