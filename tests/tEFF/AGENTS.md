# EFF (eFF / Reactive Force Field) Tests

## Purpose

CPU vs GPU parity validation for the eFF (electron force field) and RARFF reactive force fields. These are semi-quantum force fields using localized orbital ansätze.

## Ownership

- eFF/RARFF CPU/GPU parity scripts
- Relaxed scan protocols for energy and force comparison
- Component-wise energy decomposition tests

## Local Contracts

- **Run from this directory** — scripts may use relative paths to data and kernels.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **Static-force parity:** Check at `dt=0` first (should be ~1e-5). Relaxation parity is harder due to integrator/path differences.
- **Relaxed scan protocol:** `run_relax_parity_protocol.py` uses fixed ions, `dt=0.01`, damping=0.1, 2000 steps.

## Work Guidance

### Test Scripts
- `test_ocl_vs_cpu.py` — basic OpenCL vs CPU parity
- `run_process_xyz_e.py` / `run_process_xyz.py` — energy evaluation drivers
- `run_relax_parity_protocol.py` — 3-stage parity (short, long, scan)
- `run_scan_batch.sh` — batch relaxed scan parity with plotting
- `debug_singlepoint_cpu_gpu.py` — single-point comparison
- `run_dynamics.py` — dynamics test harness

### Plotting
- `plot_scan_parity.py` — plot CPU/GPU parity curves (CPU dotted lw=1.5, GPU solid lw=0.5)
- `plot_EA.py` / `plot_EE.py` — energy component plots

## Verification

- Run `run_scan_batch.sh` or `python3 test_ocl_vs_cpu.py` from this directory
- Static parity: max|dE| ~1e-5, force norms ~1e-5
- If relaxation diverges, check integrator step update and damping parameters, not just force formulas

## Child DOX Index

- No child subtrees
