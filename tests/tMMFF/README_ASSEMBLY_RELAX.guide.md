# Assembly and Relaxation Workflow Pipeline

This directory contains scripts to generate molecular assemblies (e.g., placing multiple molecules in a periodic lattice) and subsequently relax these assemblies on a surface using the PyOpenCL MMFF implementation.

## 1. Generating Assemblies (`test_assembly.py`)

The `test_assembly.py` script explores the configuration space of rigid molecules in a periodic cell by applying rotations and translations. It evaluates steric clashes and intermolecular distances to find viable packed structures.

### Example Usage:
```bash
python test_assembly.py --preset triptycene --rot_mode tilt --nrot 32 --n_tilt 9 --tilt_range 0.35 --nshift 10 --clash_max 5.0 --plot_trans --plot_rot --plot_best_k 0 --z_highlight 0.4 --dump --nPBC_test 2 --nPBC_xyz 0
```

### Key Outputs:
- A directory named `assembly_out_<preset>/` containing:
  - `mol2/`: A subdirectory with individual `.mol2` files for the best configurations (e.g., `rank_0000_idx_73359.mol2`). These files retain the original Tripos atom types (like `C.ar`) which are critical for the downstream relaxation.
  - Visualization plots (`trans_sampling_*.png`, `rot_sampling_*.png`, `pareto_*.png`).
  - A concatenated `.xyz` movie of the best configurations.

## 2. Single Configuration Relaxation (`test_ditetraceno_surface.py`)

To closely examine and relax a single configuration from the generated assemblies, use `test_ditetraceno_surface.py`. This script applies the MMFF force field along with a surface potential (Hamaker or Morse) and performs molecular dynamics relaxation.

### Example Usage:
```bash
python test_ditetraceno_surface.py --mol2 assembly_out_triptycene/mol2/rank_0000_idx_73359.mol2 --nsteps 10000 --stride 100 --dt 0.02 --damp 0.01
```

### Key Outputs:
- `trj_<name>_<mode>.xyz`: The relaxation trajectory (useful for debugging and visualizing the relaxation process).
- `trj_<name>_<mode>_final.mol2`: The final relaxed geometry in MOL2 format.

## 3. Batch Parallel Relaxation (`relax_batch_surface.py`)

When you have generated many candidate configurations and want to relax all of them efficiently, use `relax_batch_surface.py`. This script takes advantage of the GPU by packing multiple systems (`nSystems > 1`) into the OpenCL buffers and relaxing them completely in parallel.

### Example Usage:
```bash
python relax_batch_surface.py --mol2-dir assembly_out_triptycene/mol2 --out-dir assembly_out_triptycene/relaxed_batch --max-systems 50 --nsteps 5000 --stride 1000 --nPBC 1 1 0
```

### Features:
- **Parallel Execution:** Evaluates forces and integrates MD for multiple configurations simultaneously on the GPU.
- **No Trajectories:** To maximize performance, it does not write intermediate `.xyz` trajectories.
- **Final Geometries:** Outputs the final relaxed geometry for each configuration to the specified `--out-dir` as `*_relaxed.mol2`.

## Important Notes on Atom Typing:
The MOL2 format uses specific atom types in column 6 (e.g., `C.ar` for aromatic carbon, `N.2` for sp2 nitrogen). 
- `test_assembly.py` preserves these types during export.
- `MMFF.py` maps these Tripos types to internal parameters (e.g., `C_R`) ensuring that planar constraints (`Ksp`, `Kpp`) are correctly applied during relaxation. **Never manually strip or alter these types in the intermediate MOL2 files**, or the molecules will deform into incorrect sp3 geometries.
