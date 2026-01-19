# Grid Generation Script

This script generates grid for desired substrate using GPU OpenCL acceleration


## Usage
```bash
python3 generate_grid.py --fname <substrate_name> [--desired_voxel <size>]
```

### Arguments
- `--fname`: Substrate XYZ filename (with or without .xyz extension)
- `--desired_voxel`: Voxel size in Angstroms (default: 0.1)

### Examples
1. Basic usage:
   ```bash
   python3 generate_grid.py --fname NaCl
   ```

2. Custom voxel size:
   ```bash
   python3 generate_grid.py --fname NaCl --desired_voxel 0.1
   ```

## Input Files
Place substrate XYZ files in either:
- Current directory
- `FireCore/cpp/common_resources/xyz/` directory

## Output
Grid files will be saved in:
```
FireCore/cpp/common_resources/<substrate_name>/Bspline_PLQd.npy
```

## Notes
- Output includes Coulomb, Pauli, and London components

## Used Substrates for this Manuscripts:
**For Scans (Rigid and Relax) and for defect surface calculations**
- NaCl_perfect_8x8.xyz
- NaCl_perfect_20x20.xyz
- NaCl_45_defect_align_20x20.xyz

**For CPU performance comparison with LAMMPS**
- NaCl_perfect_8x8.xyz
- NaCl_perfect_12x12.xyz
- NaCl_perfect_16x16.xyz
- NaCl_perfect_20x20.xyz

**For GPU performance**
- NaCl_1x1_L3.xyz
- NaCl_2x2_L3.xyz
- NaCl_3x3_L3.xyz
- NaCl_4x4_L3.xyz
- NaCl_5x5_L3.xyz
- NaCl_6x6_L3.xyz
- NaCl_7x7_L3.xyz
- NaCl_8x8_L3.xyz
- NaCl_9x9_L3.xyz
- NaCl_10x10_L3.xyz
- NaCl_11x11_L3.xyz
- NaCl_12x12_L3.xyz
- NaCl_13x13_L3.xyz
- NaCl_14x14_L3.xyz
- NaCl_15x15_L3.xyz
- NaCl_16x16_L3.xyz

# Scan Calculations Script

This script performs molecular scans over substrates to analyze force field interactions.

## Key Features
- 1D and 2D potential energy surface scans
- Support for rigid and relaxed scans
- Comparison with LAMMPS data
- Multiple force field components: total, morse, coulomb

## Basic Usage
```bash
python3 generate_scans.py \
  --molecule <molecule.xyz> \
  --substrate <substrate.xyz> \
  --output-dir <results> \
  --scan-mode <1d|2d> \
  [additional options]
```

## Common Use Cases

### 1. Simple 1D Scan
```bash
python3 generate_scans.py \
  --molecule H2O.xyz \
  --substrate NaCl_8x8.xyz \
  --output-dir scan_results \
  --scan-mode 1d \
  --scan-dir "0,0,1" \
  --nscan 50 \
  --span-min 2.0 \
  --span-max 6.0 \
  --scan-origin "0,0,3" \
  --scan-types total

```

### 2. Relaxed 1D Scan
```bash
python3 generate_scans.py \
  --molecule CO.xyz \
  --substrate Au_111.xyz \
  --output-dir relaxed_scan \
  --scan-mode 1d \
  --scan-dir "1,0,0" \
  --relaxed \
  --niter-max 1000 \
  --fconv 0.01 \
  --dt 0.02 \
  --cons-atom 26 \
  --scan-origin "0,0,3" \
  --scan-types total
```

### 3. 2D Surface Scan
```bash
python3 generate_scans.py \
  --molecule NH3.xyz \
  --substrate graphene.xyz \
  --output-dir 2d_scan \
  --scan-mode 2d \
  --scan-dir1 "1,0,0" \
  --scan-dir2 "0,1,0" \
  --scan-origin "0,0,3" \
  --nscan-1 20 \
  --span-min-1 -5.0 \
  --span-max-1 5.0 \
  --nscan-2 20 \
  --span-min-2 -5.0 \
  --span-max-2 5.0 \
  --scan-types total
```

### 4. Compare with LAMMPS
```bash
python3 generate_scans.py \
  --molecule CH4.xyz \
  --substrate MgO.xyz \
  --output-dir comparison \
  --scan-mode 1d \
  --scan-dir "0,0,1" \
  --compare \
  --lammps-1d-dir lammps_data
```

## Key Options 

| Option | Description |
|--------|-------------|
| `--scan-mode` | `1d` for linear scans, `2d` for surface scans |
| `--scan-dir` | Direction vector for 1D scans (e.g., "0,0,1" for Z-axis) |
| `--scan-dir1/2` | Direction vectors for 2D scans |
| `--scan-origin` | Starting point for any scans |
| `--relaxed` | Enable geometry relaxation at each scan point |
| `--niter-max` | Max relaxation iterations (default: 1000) |
| `--fconv` | Force convergence threshold (default: 0.01 eV/Å) |
| `--dt` | Time step for molecular dynamics relaxation (default: 0.01 fs) |
| `--cons-atom` | Atom index to constrain during relaxation (default: None) |
| `--scan-types` | Force components: `total`, `morse`, `coulomb` |
| `--compare` | Compare with existing LAMMPS data |

## Relaxation Parameters
For relaxed scans, these additional parameters control the optimization process:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--dt` | Time step for molecular dynamics relaxation | 0.01 fs |
| `--fconv` | Force convergence threshold | 0.05 eV/Å |
| `--niter-max` | Maximum relaxation iterations | 1000 |

## Complete Examples

### 1. Rigid Scan (Section 6.1)
```bash
python3 generate_scans.py \
  --scan-types total \
  --scan-mode 1d \
  --scan-origin "0,0,0" \
  --nscan 125 \
  --span-min 2.6 \
  --span-max 15.1 \
  --molecule data/xyz/PTCDA \
  --substrate data/xyz/NaCl_perfect_8x8.xyz \
  --scan-dir "0,0,1"
```

### 2. Relaxed Scan (Section 6.2)
```bash
python3 generate_scans.py \
  --scan-types total \
  --scan-mode 1d \
  --scan-dir "0.0,0.0,1.0" \
  --scan-origin "0.0,0.0,0.0" \
  --nscan 201 \
  --span-min 1.3 \
  --span-max 21.4 \
  --molecule data/xyz/PTCDA \
  --substrate data/xyz/NaCl_perfect_8x8.xyz \
  --relaxed \
  --niter-max 100000 \
  --dt 0.02 \
  --fconv 0.001 \
  --cons-atom 26
```

## Key Notes
1. **Scan Direction**: For 1D scans, use `--scan-dir1` instead of `--scan-dir`
2. **Origin Point**: Crucial for positioning the molecule relative to the substrate
3. **Constrained Atom**: Atom 26 (typically a central atom) is fixed during relaxation
4. **Time Step**: Smaller dt (0.02 fs) provides more stable relaxation
5. **Force Convergence**: Tighter fconv (0.001 eV/Å) yields more accurate results
6. **Iteration Limit**: High niter-max (100000) ensures convergence for complex systems

## Output Files
- `scan_1d.dat`: 1D scan results
- `scan_2d.dat`: 2D scan matrix
- `relaxed_geometries/`: Optimized structures (for relaxed scans)
