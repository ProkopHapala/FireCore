# Kekule Structure Explorer Tests

## Purpose

Tests for the Kekule Structure Explorer GUI backend (`pyBall.KekuleBackend`): topology building, hydrogen passivation, H-bond detection, DFTB relaxation, and ribbon parity.

## Ownership

- Topology tests (benzene, naphthalene, anthracene, etc.)
- H-bond detection validation
- DFTB relaxation result checks
- Ribbon parity tests (zigzag ribbons with N/NH passivation)

## Local Contracts

- **Run from this directory** — scripts use relative paths to data.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **sp3 H orientation:** Recent fix for tetrahedral hydrogen placement must be preserved.
- **H-bond detection:** PBC-aware X-H...Y detection for molecules.

## Work Guidance

- `test_kekule_topology.py` — topology building for standard aromatic systems
- `test_kekule_hbonds.py` — H-bond detection validation
- `test_kekule_relax.py` — DFTB relaxation tests
- `test_ribbon_parity.py` — ribbon structure parity (zigzag, N/NH passivation)

### Results Directories
- `out_topology/` — 18 geometry snapshots (benzene, naphthalene, anthracene, etc.)
- `out_relax/` — DFTB relaxation results

## Verification

- Run `python3 test_kekule_topology.py` from this directory
- Topology should match expected ring counts and atom types
- H-bond detection should find correct donor/acceptor pairs

## Child DOX Index

- No child subtrees
