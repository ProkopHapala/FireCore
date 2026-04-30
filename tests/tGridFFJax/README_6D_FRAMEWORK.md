# GridFF-JAX 6D Framework

Generic 6D (x, y, z, α, β, γ) sampling, fitting, and MLIP diagnostic framework for molecule-on-surface systems, built on JAX for GPU/TPU/multi-CPU acceleration.

## Overview

This framework addresses the fundamental limitation of the previous z-scan-only fitting pipeline: **parameters fitted using only vertical (z) variation at a few high-symmetry sites cannot reproduce lateral corrugation accurately.** The 6D framework systematically samples the full rigid-body configuration space and fits transferable per-element parameters from this comprehensive dataset.

### Architecture

```
┌────────────────────────────────────────────────────────┐
│  1. 6D Pose Sampler                                     │
│     (u, v, z, α, β, γ) systematic grid + random poses  │
│     Input: molecule list file + substrate density        │
├────────────────────────────────────────────────────────┤
│  2. Teacher Evaluation                                   │
│     MAD-SURF / MACE / any MLIP / DFT (future)          │
│     Output: energies + forces for all poses              │
├────────────────────────────────────────────────────────┤
│  3. 6D Fitting Pipeline                                  │
│     Stratified train/val/test split                      │
│     Multi-molecule joint fitting → transferable params   │
├────────────────────────────────────────────────────────┤
│  4. MLIP Diagnostics                                     │
│     Smoothness, force consistency, corrugation, compare  │
├────────────────────────────────────────────────────────┤
│  5. Export                                               │
│     → Bspline_PLQd.npy (existing C++/OpenCL runtime)    │
│     → element_pairs_Ag.json (transferable parameters)    │
└────────────────────────────────────────────────────────┘
```

## Quick Start

### Step 1: Generate 6D Reference Data

```bash
# Single molecule (CHONH2 = formamide)
python run_6d_sampling.py --molecules CHONH2 \
    --n-u 10 --n-v 10 --n-z 20 --n-orient 8

# Multiple molecules from a file
python run_6d_sampling.py --molecule-file molecules_tier1.txt

# Custom substrate (Cu instead of Ag)
python run_6d_sampling.py --molecules H,CO,H2O \
    --chgcar /path/to/Cu/CHGCAR --locpot /path/to/Cu/LOCPOT \
    --substrate Cu --model-path /path/to/model.model
```

**Output:** `datasets_6d/<molecule>_6d.npz` + `.json` for each molecule.

### Step 2: Fit Transferable Parameters

```bash
# Single molecule fit
python run_6d_fit.py --dataset-dir datasets_6d/ --molecules CHONH2

# Multi-molecule joint fit (recommended for transferability)
python run_6d_fit.py --dataset-dir datasets_6d/ \
    --molecules H,CO,H2O,CHONH2 --two-stage-c6 --raw-r6

# Custom training
python run_6d_fit.py --dataset-dir datasets_6d/ \
    --molecules CHONH2 --max-steps 800 --force-weight 15.0
```

**Output:** `datasets_6d/fit_results/fit_params.json`, `element_pairs_Ag.json`

### Step 3: Validate

```bash
python run_6d_validate.py --fit-dir datasets_6d/fit_results/ \
    --molecules CHONH2 --z-heights 2.0,2.5,3.0,3.5
```

**Output:** Per-molecule z-scan and 2D lateral scan comparisons with error metrics.

### Step 4: Run MLIP Diagnostics

```bash
# Basic diagnostics (force consistency + smoothness + corrugation)
python run_mlip_diagnostics.py --molecules CHONH2

# With GridFF comparison
python run_mlip_diagnostics.py --molecules CHONH2 \
    --fit-dir datasets_6d/fit_results/ --compare-gridff
```

**Output:** `mlip_diagnostics/<molecule>/diagnostic_report.json` + raw data.

## Molecule List File Format

The sampler reads molecules from a plain text file (`molecules_tier1.txt`):

```text
# Comments start with '#'
H                                    # Built-in name
CHONH2                               # Built-in name
CO                                   # Built-in name
/path/to/benzene.xyz  name=C6H6     # XYZ file with custom name
/path/to/mol.xyz  anchor=2 charges=true  # With options
```

Built-in molecules: `H`, `C`, `O`, `N`, `CO`, `H2O`, `CHONH2`.
Custom molecules: provide an XYZ file path.

## 6D Sampling Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_u` | 10 | Lateral grid points along u (fractional a) |
| `n_v` | 10 | Lateral grid points along v (fractional b) |
| `n_z` | 20 | Height points per (u,v) position |
| `z_range` | (1.5, 5.5) | Height range above surface (Å) |
| `z_bias_power` | 1.5 | >1 concentrates samples near surface |
| `n_orient` | 8 | Molecular orientations (Fibonacci sphere + yaw) |
| `tilt_max_deg` | 60.0 | Maximum tilt from surface normal |
| `n_yaw` | 4 | Yaw angles per tilt direction |
| `random_fraction` | 0.1 | Extra random poses (fraction of systematic) |
| `include_high_symmetry_sites` | true | Insert top/bridge/hollow explicitly |

**Total poses per molecule** ≈ `n_u × n_v × n_z × n_orient × (1 + random_fraction)`

Example: 10 × 10 × 20 × 8 × 1.1 ≈ 17,600 poses.

## MLIP Diagnostic Tools

### Force-Energy Consistency
Compares analytical MLIP forces with numerical `dE/dr` via central finite differences.
**Key metric:** RMS |F_analytical − F_numerical| — if this is large, the MLIP energy surface is not smooth.

### Smoothness Analysis
Evaluates energy over a (u,v) grid at multiple z-heights. Computes:
- Gradient magnitude maps
- Laplacian maps (2nd derivative)
- Roughness metric: std(Laplacian) — high values indicate non-physical oscillations

### Corrugation Analysis
Compares lateral energy variation (E_max − E_min over unit cell) between methods.
Reports per-site energies at top/bridge/hollow.

### Multi-Method Comparison
Generic parity statistics between any pair of methods:
- Energy RMSE/MAE/bias/R² (meV)
- Force RMSE/MAE (eV/Å)
- Per-site and z-binned breakdown

## Transferable Parameters

The framework produces **element-substrate pair parameters** (e.g., H-Ag, C-Ag):

```json
{
  "substrate": "Ag",
  "pairs": {
    "H-Ag": {"radius_offset": 0.12, "energy_scale": 1.5, "c6_coeff": 23.0, ...},
    "C-Ag": {"radius_offset": -0.05, "energy_scale": 2.1, "c6_coeff": 45.0, ...},
    ...
  }
}
```

These are analogous to UFF parameters but specific to the substrate and fitted from first-principles (MLIP or DFT) data. Multi-molecule joint fitting ensures they transfer across different adsorbates.

**Future (Phase 2):** Atom-type-aware parameters (C_R-Ag vs C_3-Ag) using connectivity-based typing. The data structures are already in place in `atom_types.py`.

## Code Layout

```
pyBall/gridff_jax/
├── pose_sampling/
│   ├── sampler_6d.py        # 6D systematic pose generator
│   ├── molecules.py         # Built-in molecule definitions
│   ├── rigid.py             # Pose transformation utilities
│   └── sites.py             # Surface site detection
├── fit/
│   ├── fit_6d.py            # 6D dataset, splitting, merging
│   ├── optimize.py          # Parameter optimization (existing)
│   └── dataset.py           # Legacy dataset utilities
├── diagnostics/
│   ├── force_consistency.py # F vs dE/dr check
│   ├── smoothness.py        # Energy surface smoothness
│   ├── corrugation.py       # Lateral corrugation analysis
│   ├── comparison.py        # Multi-method parity statistics
│   └── report.py            # Automated report generation
├── atom_types.py            # Element-pair parameter registry
├── config.py                # Configuration (incl. Sampler6DConfig)
└── interfaces.py            # Shared data structures

tests/tGridFFJax/
├── run_6d_sampling.py       # Generate 6D reference datasets
├── run_6d_fit.py            # Fit parameters from 6D data
├── run_6d_validate.py       # Validate with lateral + z scans
├── run_mlip_diagnostics.py  # MLIP diagnostic suite
├── molecules_tier1.txt      # Example molecule list
└── README_6D_FRAMEWORK.md   # This file
```

## Extensibility

### Adding a new MLIP teacher
Implement `TeacherBackend.evaluate_batch(density, poses) -> TeacherResult` in
`pyBall/gridff_jax/teacher_backends/`. The diagnostics and fitting pipelines
work with any backend that provides this interface.

### Adding a new substrate
1. Provide CHGCAR + LOCPOT (or equivalent density data)
2. Set `--substrate <element>` and optionally tune `metal_bulk_electron_density`
3. Run the same pipeline — all code is substrate-agnostic

### Adding molecular deformations (7D+)
The `Sampler6DConfig.extra_dimensions` field accepts `ExtraDimension` entries:
```python
from pyBall.gridff_jax.pose_sampling.sampler_6d import ExtraDimension
config.extra_dimensions = [
    ExtraDimension(name="bond_stretch", range=(-0.1, 0.1), n_points=5),
]
```
The data structure is ready; the sampling logic for deformations is a future extension.

### Adding DFT as teacher
Implement a `VaspSinglePointBackend` that runs VASP for each pose and returns
energies + forces in the `TeacherResult` format. The fitting and diagnostics
pipelines will work unchanged.

## Performance Notes

- **JAX on GPU/TPU**: Fitting uses JAX autograd + optax for GPU-accelerated optimization. `jax.vmap` parallelizes batch evaluation.
- **Multi-CPU**: JAX XLA backend automatically parallelizes across CPU cores when GPU is unavailable.
- **Memory**: ~17K poses × 6 atoms × 3 coords ≈ 300 KB. PLQ grids (~500 MB) dominate memory.
- **Runtime GridFF**: Fitted parameters export to `Bspline_PLQd.npy` for the existing C++/OpenCL runtime. No changes to the fast evaluation path.
