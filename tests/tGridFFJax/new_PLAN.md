# Plan: Plug-and-Play Substrate-Aware DFT/MLIP-to-GridFF Workflow

## Summary
Build a v1 “GridFF auto” workflow that ingests DFT/MLIP outputs, validates them, selects the correct substrate physics recipe, fits GridFF, and exports reproducible artifacts. V1 targets **Ag/metal + NaCl/ionic** first, while structuring the backend so perovskite/ferroelectric recipes can be added later without changing the user interface.

Core principle:

```text
field source = what builds the grid
teacher source = what fits the grid
substrate recipe = how physics is assembled
```

The workflow will **ingest outputs only** in v1. It will not run VASP, Psi4, QE, PySCF, or MLIP jobs.

## Key Changes

- Add a user-facing bundle workflow:
  - Primary entry: `--bundle /path/to/gridff_bundle.json`
  - Keep CLI overrides for expert users: `--substrate-class`, `--density-kind`, `--chgcar`, `--locpot`, `--density-xyz`, `--teacher-npz`, pose-grid args.
  - Bundle records substrate class, field source, teacher labels, adsorbate, pose grid, units, and method metadata.

- Split substrate physics recipes cleanly:
  - `metal`: VASP/cube density + potential, `metal_density_plq`, density-derived Pauli/London, LOCPOT/potential Coulomb, metal screening/image-compatible settings.
  - `ionic_pointcharge`: charged `surface_xyz`, current `parity_core`, pairwise Pauli/London, Ewald/point-charge Coulomb, no image plane, no metal C6 stage.
  - `ionic_dft_volumetric`: density/potential cube or CHGCAR/LOCPOT input, pairwise ion-resolved Pauli/London, DFT potential as Coulomb channel, no metallic image/screening assumptions.
  - Perovskite/ferroelectric recipes remain explicit stubs until physics is defined.

- Standardize teacher labels:
  - One canonical label format: `pose_id`, `pose_params`, `positions`, `energies_eV`, optional `forces_eV_per_A`, `adsorbate_symbols`, `adsorbate_charges`, metadata.
  - Sources may be Psi4, VASP, QE, PySCF, MACE, NequIP, ML-DFT, or manual data, but must be imported into the canonical format.
  - Existing `precomputed` teacher becomes the stable interchange layer.

- Add hard validation before fitting:
  - Field source and teacher labels must agree on substrate geometry, cell, adsorbate, pose grid, units, and coordinate frame.
  - Ionic point-charge mode must fail if charges are missing or all zero.
  - Ionic modes must fail if metal image/metal-density settings are accidentally enabled.
  - Cache signatures must include full config, field source, teacher identity, pose grid, adsorbate charges, and substrate geometry.
  - Store adsorbate charges in all cached datasets and reload paths.
  - Export failure must fail the run, not only warn.

## Implementation Shape

- Add a small “bundle loader + validator” layer in `pyBall/gridff_jax` that converts either a manifest or CLI flags into a locked `RunConfig`.
- Refactor `benchmark_substrate_6d.py` so it can consume the locked config instead of manually mixing substrate class, density kind, and teacher kind.
- Add a new builder mode name for ionic DFT volumetrics, e.g. `ionic_dft_plq`, using pairwise ion PL terms plus DFT potential for Coulomb, not metal `rho^power` logic.
- Write output artifacts consistently:
  - `input_manifest.lock.json`
  - `run_config.lock.json`
  - `validation_report.json`
  - `apples_dataset.npz`
  - `fit_params_full.json`
  - exported GridFF files
  - `SUMMARY.md`

## Test Plan

- Unit tests:
  - recipe selection for `metal`, `ionic_pointcharge`, `ionic_dft_volumetric`
  - manifest parsing and CLI override precedence
  - required-file validation for VASP, cube, xyz, and precomputed teacher labels
  - adsorbate charges preserved through save/load
  - cache invalidates when teacher, density, pose grid, or charges change
  - force shape and unit validation

- Smoke tests:
  - NaCl point-charge + synthetic/precomputed labels
  - NaCl Psi4-style precomputed labels, energy-only and energy+forces
  - Ag volumetric path with existing CHGCAR/LOCPOT-style config
  - cube density + potential path with generated tiny cube fixtures

- Failure tests:
  - ionic xyz without charges fails
  - ionic substrate accidentally using metal image settings fails
  - mismatched pose grid between labels and requested fit fails
  - missing potential in `ionic_dft_volumetric` fails unless explicit fallback is requested

## Assumptions

- V1 supports Ag-like metals and NaCl-like ionic substrates robustly first.
- External DFT/MLIP execution is outside v1; users provide outputs.
- Psi4/PySCF/QE data enter through cube files and/or canonical precomputed labels.
- The framework guarantees consistency, validation, and correct recipe selection; physical accuracy still depends on the supplied DFT/MLIP data quality.
- Existing expert CLI workflows remain available, but the bundle path becomes the recommended user-facing workflow.
