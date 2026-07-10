---
type: TopicalAudit
title: Dihedral/Torsion Terms in test_FFfit.py
description: Implementation, status, and open issues for adding UFF/Prokop torsion terms to the Python FFfit pipeline.
tags: [fffit, forcefield, dihedral, torsion, vibration, uff]
timestamp: 2026-07-09
---

# Dihedral/Torsion Terms in `test_FFfit.py`

This file documents the recent attempt to add UFF/Prokop-style torsion (dihedral) terms to the Python force-field fitting pipeline (`tests/tSiNCs/test_FFfit.py`). The goal is to improve the low-frequency and torsion-dominated part of vibrational spectra for Si and C nanocrystals by fitting an additional set of transferable torsion parameters.

## Summary

The implementation is functionally complete: Python mirrors `evalDihedral_Prokop()` from `cpp/common/molecular/UFF.h`, dihedral topology is built from the bond graph, a global dihedral type map is created, and the dihedral sensitivity is added to the mode-basis fit and to `H_model` for frequency comparison. However, the fit becomes numerically unstable and the resulting vibrational spectra are worse than without dihedrals. The root cause is that the dihedral sensitivity `A_p = dH/dV` is not positive semidefinite for the real reference geometries, because the reference `phi` values are not exactly at the minima of the `1 + d cos(n phi)` potential. The resulting `f' * C` term in the torsion Hessian adds negative curvature and is strongly correlated with the angle sensitivity, which causes the angle stiffnesses to collapse to near zero.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| Python | [`tests/tSiNCs/test_FFfit.py`](../../../tests/tSiNCs/test_FFfit.py) | functional but unstable | `dihedral_energy_gradient`, `dihedral_hessian`, `build_dihedrals`, `compute_dihedral_sensitivity`, `add_dihedral_hessian` and CLI flags `--dihedrals`, `--dihedral-n`, `--dihedral-d`, `--dihedral-h` |
| C++ reference | [`cpp/common/molecular/UFF.h`](../../../cpp/common/molecular/UFF.h) | reference | `evalDihedral_Prokop()` (lines ~918-1000) computes energy and forces for a single torsion |
| C++ engine | [`cpp/common/molecular/FFfit.h`](../../../cpp/common/molecular/FFfit.h) | new / untracked | theory and `FFfit` class; no dihedral code yet but documents the Hessian/fitting theory |
| C++ wrapper | [`cpp/libs/Molecular/FFfit_lib.cpp`](../../../cpp/libs/Molecular/FFfit_lib.cpp) | new / untracked | ctypes C-wrapper exposing `FFfit` to Python; `add_dihedral` not yet added |
| Build config | [`cpp/libs/Molecular/CMakeLists.txt`](../../../cpp/libs/Molecular/CMakeLists.txt) | modified | adds `FFfit_lib` shared library |
| Python wrapper | [`pyBall/FFfit.py`](../../../pyBall/FFfit.py) | new / untracked | ctypes Python wrapper for `FFfit_lib` |
| Python plotting | [`tests/tSiNCs/plot_pyscf_vib_results.py`](../../../tests/tSiNCs/plot_pyscf_vib_results.py) | new | plots PySCF vibration results and spectra |
| PySCF driver | [`tests/tSiNCs/run_small_np_pyscf_vib.py`](../../../tests/tSiNCs/run_small_np_pyscf_vib.py) | new | full PySCF relax + vibration pipeline for small nanoparticles |
| Viewer plugin | [`pyBall/GUI/mol_browser_plugins/vibration_spectrum.py`](../../../pyBall/GUI/mol_browser_plugins/vibration_spectrum.py) | modified | supports PySCF `.npy` format in addition to nanocrystal pipeline NPZ |

## Algorithm and Energy Form

The implemented torsion term follows `evalDihedral_Prokop()` in `UFF.h` without the optional non-bonded subtraction (`bSubNonBond = false`):

```
E = V * (1 + d * cos(n * phi))
```

- `V` is the fitted dihedral barrier (in eV).
- `d` is the phase sign (`+1` or `-1`).
- `n` is the periodicity (integer, default `3` for sp3).
- `phi` is the torsion angle defined by the four atoms `i-j-k-l`.

The gradient is computed from the cross-product vectors and the scalar `f = -V * d * n * sin(n * phi)`, exactly as in the C++ code. The Hessian is obtained by central finite differences of the analytical gradient. The per-type sensitivity `A_p` is the sum of all per-dihedral Hessians `d²E/dr²` with `V = 1`.

Key functions in `test_FFfit.py`:

- `dihedral_energy_gradient(pos, d, n)` — analytical energy and `(4, 3)` gradient.
- `dihedral_hessian(pos, d, n, h)` — `(12, 12)` symmetric Hessian via FD.
- `build_dihedrals(symbols, positions, bonds, d, n, dihedral)` — enumerates all `i-j-k-l` sequences around each central bond.
- `dihedral_type_key(symbols, i, j, k, l)` — groups dihedrals by central-bond element pair and sorted outer-atom pair.
- `compute_dihedral_sensitivity(...)` — builds the full `(3N, 3N)` `A_p` matrix for each dihedral type.
- `add_dihedral_hessian(H, k, dihedral_A)` — adds `Σ_p k[p] * A_p` to the model Hessian.

## Integration into the Workflow

1. **Topology** — `build_dihedrals()` is called after `build_topology()` if `--dihedrals` is set. The `d` and `n` values are taken from the CLI or default to `d=1`, `n=3`.
2. **Global map** — `build_global_param_map()` now returns a `global_dihedral_map` in addition to the bond/angle maps. Dihedral parameter indices are offset after bond/angle/3rd-bond indices.
3. **Precomputation** — `compute_dihedral_sensitivity()` is called once per system and stored in `sys['dihedral_A']`.
4. **C++ fitters** — `FFfit` instances are still created with `n_total` parameters (including dihedral slots), but `set_n_free(n_cpp)` is used for the C++ LSQ/GD methods, which ignore dihedrals. `set_n_free(n_total)` is used for the Python mode-basis fit.
5. **Mode-basis fit** — `fit_modes_multi()` accepts `dihedral_A_per_system` and passes it to `get_sensitivity_action()`, which adds `A_p` to the C++ `H` for each one-hot parameter vector.
6. **Evaluation** — `add_dihedral_hessian()` is called before `compare_frequencies()` and in the mode-basis diagnostics so the model Hessian includes the dihedral contribution.

## CLI Flags

```bash
python3 tests/tSiNCs/test_FFfit.py --cases all_Si --dihedrals --dihedral-n 3 --dihedral-d 1.0 --dihedral-h 1e-5
```

- `--dihedrals` — enable dihedral term detection and fitting.
- `--dihedral-n` — periodicity (default `3`).
- `--dihedral-d` — phase sign (default `1.0`).
- `--dihedral-h` — finite-difference step for the Hessian (default `1e-5` Å).

## Test Results

### Without dihedrals (baseline)

- Si systems: `relFrob` ~4-6%.
- C diamond systems: `relFrob` ~10-12%.

### With `--dihedrals` (default `d=1`, `n=3`)

- Si systems:
  - `k_angle` collapses to ~0.001-0.07 eV/rad².
  - `k_dihedral` ~0.004-0.035 eV.
  - `relFrob` degrades to ~10-12%.
  - Low model frequencies appear (e.g. 21.9 cm⁻¹ vs reference 36.9 cm⁻¹).
- C diamond systems:
  - `relFrob` degrades to ~50-60%.
  - Many dihedrals in C diamond have `phi` ≈ 120°, which is a **maximum** of `1 + cos(3 phi)` with `d=1`, so the sensitivity is negative in the torsion direction.

## Root Cause: `A_d` is Not Positive Semidefinite

For a torsion term `E = V * (1 + d cos(n phi))` the Hessian is

```
H = V * [ f'' * b ⊗ b + f' * C ]
```

where `b = ∂phi/∂r` and `C = ∂²phi/∂r²`. The term `f' * C` is non-zero whenever `phi` is not exactly at a minimum of `1 + d cos(n phi)`. It is not positive semidefinite and, for the short Si-H/C-H bonds in the 4-atom torsion, `C` can have very large eigenvalues. Consequences:

1. `A_d = dH/dV` has negative eigenvalues.
2. Fitting `V > 0` adds negative curvature to the model Hessian.
3. `A_d` is highly correlated with the angle sensitivity `A_angle`, so the mode-basis fit transfers stiffness from angles to dihedrals (`k_angle` collapses).
4. For C diamond many dihedrals have `phi` near 120°, i.e. at a maximum when `d=1`, so the torsion direction itself is destabilized.

## Open Issues

1. **Torsion sensitivity sign/stability** — Decide whether to use the exact `A_d` (with `f' * C`), the pure curvature term `f'' * b ⊗ b`, or a phase-selected `d` per dihedral to keep the nearest minimum close to the actual `phi`.
2. **Dihedral/angle correlation** — The fit is ill-conditioned because dihedral and angle sensitivities overlap. A regularization or NNLS prior that keeps `k_angle` near the C++ LSQ value may be needed.
3. **Per-dihedral phase selection** — For C diamond, `d=1` is wrong for `phi ≈ 120°`; `d` should be `-1` for those dihedrals. This may require splitting a single "type" into subtypes by `d` or by `phi` sector.
4. **C++ `FFfit` integration** — `FFfit_lib.cpp` and `FFfit.h` currently have no `add_dihedral` / `set_dihedral_param`. Everything is done in Python. Moving the torsion evaluation to C++ would let us use the exact UFF machinery and the same parameter map for the C++ LSQ/GD methods.
5. **Test harness** — No dedicated unit test for `dihedral_energy_gradient` against finite differences exists in the repository. A quick `pytest`/`test_` script should be added to guard regressions.

## Next Steps

1. Choose the treatment of `A_d` (exact vs. pure-curvature vs. phase-selected).
2. If pure-curvature is chosen, implement `A_d = f'' * b ⊗ b` and update `add_dihedral_hessian()` accordingly.
3. Add per-dihedral phase selection (`d = -sign(cos(n phi))`) or per-type phase selection to avoid maxima for C diamond.
4. Re-run `--dihedrals` on `all_Si` and `all` and compare `relFrob` and frequency residuals.
5. Once stable, consider adding `fffit_add_dihedral` to `FFfit_lib.cpp` and `FFfit.h` so the C++ methods can also fit dihedrals.

## 2026-07 correction: signed angle and inverse bond lengths

The original Python mirror contained two errors that invalidated its torsion Hessian:

- The C++ complex pair is `(cos(phi), sin(phi))`, but Python interpreted the second component as `-sin(phi)`. Reported signed torsions therefore had the opposite sign.
- `UFF::hneigh[].w` stores inverse bond length. Python used the bond length itself in the endpoint forces and used the inverse ratios in the wrong direction.

`dihedral_angle()` now returns the signed UFF/Prokop angle from `atan2(sin(phi), cos(phi))`. `dihedral_energy_gradient()` uses inverse lengths exactly as `evalDihedral_Prokop()`. A dedicated test compares all 12 analytical gradient components with central finite differences of the torsion energy.

With these corrections, exact UFF torsion sensitivities can be compared fairly. For `Si_R3p8`, however, bounded hybrid fitting still drives their amplitudes to zero at the current `d=1, n=3` phase. Thus they do not improve this case; the remaining issue is the physical form/phase and indefinite prestress curvature, not the angle or gradient implementation.
