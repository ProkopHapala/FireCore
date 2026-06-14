# UFF Force Field Tests

## Purpose

CPU vs GPU parity validation for the Universal Force Field (UFF): bonds, angles, dihedrals, inversions, and non-bonded terms.

## Ownership

- Python test scripts and parity harnesses
- `run.sh` / `make.sh` automation scripts (if present)

## Local Contracts

- **Run from this directory** — scripts use relative paths to `common_resources/cl/*.cl`.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **Multi-system GPU:** `bMMFF=True` is required for UFF initialization (`makeFFs()` path).
- **2D NDRange:** `global.y = nSystems` with per-system buffer offsets (`i0a`, `i0b`, etc.).

## Work Guidance

### Test Scripts
- `test_UFF_multi.py` — CPU vs GPU parity for multi-system evaluation
- `test_UFF_ocl.py` — OpenCL-specific tests
- `test_MMFF_ocl_parity.py` / `test_MMFF_multi_parity.py` — MMFF parity
- `test_parity_suite.py` — consolidated parity suite
- `run.py` — general UFF/MFF test driver

### Key Conventions
- Force assembly uses auxiliary buffers → global `aforce` to avoid atomics
- Exclusion lists (`excl[]`) or subtraction method for 1-2 / 1-3 bonded neighbors
- PBC shifts precomputed in `shifts[npbc]`

## Verification

- Run `run.sh` or `python3 test_UFF_multi.py` from this directory
- Parity tolerance: energy and forces should match CPU reference to ~1e-5
- Full stdout must be visible (no `| tail`, `| head`)

## Child DOX Index

- No child subtrees
