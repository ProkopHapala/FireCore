#!/usr/bin/env python3
"""Proof: MACE z-scan → PLQ fit → GridFF scan → comparison.

SAME Ag surface used for everything:
  • PLQ fields  built from CHGCAR/LOCPOT (3×3 Ag slab, 36 atoms, 4 layers)
  • MACE scan   uses the SAME 36 Ag atoms from the CHGCAR as the substrate
  • PLQ scan    evaluates E_PLQ at the same positions

Pipeline
--------
1.  Load Ag slab positions + PLQ B-spline fields from CHGCAR/LOCPOT.
2.  Read molecule (C2H7N) from extxyz – first MD frame, orient flat.
3.  MACE z-scan  z_rel = 2 → 12 Å above Ag top layer, step 0.25 Å:
      E_sub   = MACE(Ag only, pbc=T)          ← computed once
      E_mol   = MACE(mol only, same cell/pbc, at vacuum midpoint) ← once
      E_tot   = MACE(Ag + mol @ z, pbc=T)
      E_int   = E_tot − E_sub − E_mol
      Fz_int  = Σ_j Fz_j(Ag+mol, mol-atoms)  ← total z-force on molecule
4.  Build design matrices:
      X_E[k, j*Nel+el]  = Σ_{i∈el} V_j( r_i^k )       (field values)
      X_F[k, j*Nel+el]  = Σ_{i∈el} −∂V_j/∂z( r_i^k )  (−z-gradient)
5.  Linear regression (numpy.linalg.lstsq):
      theta_E   – fit to energies only
      theta_EF  – joint fit to [E_int ; w·Fz_int]
6.  PLQ prediction:
      E_PLQ(z)  = X_E @ theta_EF
      Fz_PLQ(z) = X_F @ theta_EF
7.  Four comparison plots.

Cell constraint: Ag top at z=22.056 Å, cell_z=37.052 Å → max z_rel=12 Å
  (sufficient: MACE cutoff ≈6 Å, E_int → 0 well before 12 Å)

Data files
----------
  CHGCAR : /home/niel/git/ORR_HER_Ag_Colab/.../final_scf_12x12x1/CHGCAR
  LOCPOT : /home/niel/git/ORR_HER_Ag_Colab/.../workfunc_12x12x1/LOCPOT
  MACE   : tests/tGridFFJax/mad-surf_data/models/…/MACE_model.model
  EXTXYZ : tests/tGridFFJax/mad-surf_data/full_train_test_std_config_types.extxyz

Usage
-----
  python3 tests/tGridFFJax/mace_vs_plq_zscan.py
  python3 tests/tGridFFJax/mace_vs_plq_zscan.py --device cuda
  python3 tests/tGridFFJax/mace_vs_plq_zscan.py --force-recompute
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from collections import Counter

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_proof")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ────────────────────────────────────────────────────────────────────
HERE     = os.path.dirname(os.path.abspath(__file__))
FIRECORE = os.path.abspath(os.path.join(HERE, "..", ".."))
EXTXYZ   = os.path.join(HERE, "mad-surf_data",
                         "full_train_test_std_config_types.extxyz")
MACE_MDL = os.path.join(HERE, "mad-surf_data", "models",
                         "full_dataset_config_weights", "MACE_model.model")
CHGCAR   = ("/home/niel/git/ORR_HER_Ag_Colab/results/"
            "Ag_ORR_HER/slab_clean/final_scf_12x12x1/CHGCAR")
LOCPOT   = ("/home/niel/git/ORR_HER_Ag_Colab/results/"
            "Ag_ORR_HER/slab_clean/workfunc_12x12x1/LOCPOT")
MADSURF  = os.path.join(HERE, "MAD-SURF-main")
OUT_DIR  = os.path.join(HERE, "dft_plq_analysis")

for _p in [FIRECORE, MADSURF]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

ELEMENTS = ["C", "H", "N", "O"]
N_EL     = len(ELEMENTS)
EL_MAP   = {e: i for i, e in enumerate(ELEMENTS)}
_RE_CFG  = re.compile(r'config_type=(\S+)')


# ════════════════════════════════════════════════════════════════════════════
# 1.  Ag slab from CHGCAR  +  PLQ B-spline fields
# ════════════════════════════════════════════════════════════════════════════

def load_substrate(chgcar_path, locpot_path, out_dir):
    """Load Ag slab geometry and PLQ B-spline coefficients from VASP files.

    Returns
    -------
    ag_syms  : list[str]
    ag_pos   : (36, 3)  Cartesian Å
    cell     : (3, 3)
    coeffs   : (nz, ny, nx, 3)  B-spline coeffs for [P, L, Q]
    meta     : dict  cell, origin, voxel, grid_xyz, z_ref
    """
    from pyBall.gridff_jax.config import (DensityBackendConfig, GridConfig,
                                           FeatureToggles)
    from pyBall.gridff_jax.density_backends.vasp_volumetric import (
        VaspVolumetricBackend)
    from pyBall.gridff_jax.substrate_fields import build_substrate_grids

    cache_c = os.path.join(out_dir, "_Bspline_PLQ_coeffs.npy")
    cache_m = os.path.join(out_dir, "_Bspline_PLQ_meta.npy")

    cfg     = DensityBackendConfig(kind="vasp_volumetric",
                                   chgcar_path=chgcar_path,
                                   locpot_path=locpot_path)
    density = VaspVolumetricBackend(cfg).load()
    ag_syms = list(density.symbols)
    ag_pos  = np.asarray(density.positions, dtype=float)
    cell    = np.asarray(density.cell,      dtype=float)

    print(f"  Ag slab : {len(ag_syms)} atoms  "
          f"z=[{ag_pos[:,2].min():.3f}, {ag_pos[:,2].max():.3f}] Å")
    print(f"  Cell    : {cell[0]}\n            {cell[1]}\n            {cell[2]}")

    if os.path.exists(cache_c) and os.path.exists(cache_m):
        print(f"  PLQ coeffs: loading cache → {cache_c}")
        coeffs = np.load(cache_c)
        meta   = np.load(cache_m, allow_pickle=True).item()
    else:
        gc    = GridConfig(builder_mode="metal_density_plq")
        mol_el = tuple(sorted(set(ELEMENTS)))
        print(f"  Building PLQ B-spline grids (grid {density.grid_shape_xyz})…")
        t0    = time.perf_counter()
        grids = build_substrate_grids(density, reactive_elements=mol_el,
                                      grid_config=gc, toggles=FeatureToggles(),
                                      prefer_jax=False)
        print(f"  PLQ built in {time.perf_counter()-t0:.1f}s")
        coeffs = np.stack([grids["pauli_coeff_zyx"],
                           grids["london_coeff_zyx"],
                           grids["coulomb_coeff_zyx"]], axis=-1)
        meta   = {"cell":     density.cell,
                  "origin":   density.origin,
                  "voxel":    density.voxel,
                  "grid_xyz": density.grid_shape_xyz,
                  "z_ref":    grids["metadata"]["z_ref"]}
        np.save(cache_c, coeffs)
        np.save(cache_m, meta)  # type: ignore[arg-type]
        print(f"  Cached → {cache_c}")

    print(f"  PLQ shape : {coeffs.shape}   z_ref = {meta['z_ref']:.4f} Å")
    return ag_syms, ag_pos, cell, coeffs, meta


# ════════════════════════════════════════════════════════════════════════════
# 2.  Molecule from extxyz  (one frame, flat orientation)
# ════════════════════════════════════════════════════════════════════════════

def read_molecule(extxyz_path):
    """Extract non-Ag atoms from first MD_interface_Ag frame.  Geometry only."""
    with open(extxyz_path, "r", errors="ignore") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            try:
                n = int(line.strip())
            except ValueError:
                continue
            info = fh.readline()
            atoms = [fh.readline() for _ in range(n)]
            m = _RE_CFG.search(info)
            if not m or m.group(1) != "MD_interface_Ag":
                continue
            syms, pos = [], []
            for a in atoms:
                p = a.split()
                if len(p) >= 4 and p[0] != "Ag":
                    syms.append(p[0])
                    pos.append([float(p[1]), float(p[2]), float(p[3])])
            print(f"  Molecule : {len(syms)} atoms  {dict(Counter(syms))}")
            return syms, np.array(pos, dtype=float)
    raise RuntimeError("no MD_interface_Ag frame found")


def orient_flat(mol_pos):
    """PCA rotation: molecular normal → z.  Returns centered flat positions."""
    com = mol_pos.mean(axis=0)
    c   = mol_pos - com
    _, V = np.linalg.eigh(c.T @ c)   # ascending eigenvalues
    R   = np.column_stack([V[:, 2], V[:, 1], V[:, 0]])
    fp  = c @ R
    if fp[:, 2].mean() < 0:
        fp[:, 2] *= -1
    return fp


# ════════════════════════════════════════════════════════════════════════════
# 3.  MACE z-scan
# ════════════════════════════════════════════════════════════════════════════

def mace_zscan(ag_syms, ag_pos, mol_syms, mol_flat,
               cell, z_scan, mace_path, device):
    """Run MACE for each z in z_scan.  Returns dict of arrays."""
    from mace.calculators import MACECalculator  # type: ignore
    from ase import Atoms  # type: ignore

    print(f"\n  MACE model : {mace_path}")
    calc = MACECalculator(model_paths=mace_path,
                          device=device, default_dtype="float64")

    # Ag top layer: highest atom not near cell top (PBC bottom images live at z≈0)
    cell_z   = float(cell[2, 2])
    real_m   = ag_pos[:, 2] < (cell_z - 1.0)
    z_top_ag = float(ag_pos[real_m, 2].max())
    atop_xy  = ag_pos[real_m][np.argmax(ag_pos[real_m][:, 2]), :2]
    print(f"  z_top_ag = {z_top_ag:.4f} Å   atop site xy = {atop_xy}")

    # E_sub (substrate, once)
    sub  = Atoms(symbols=ag_syms, positions=ag_pos,
                 cell=cell, pbc=[True, True, True])
    sub.calc = calc
    e_sub = float(sub.get_potential_energy())
    print(f"  E_sub    = {e_sub:.4f} eV  ({len(ag_syms)} Ag atoms)")

    # E_mol (molecule at vacuum midpoint, same cell/pbc, no Ag)
    z_bot_img = cell_z - float(ag_pos[real_m, 2].min()) + float(ag_pos[real_m, 2].min())
    # safe reference: halfway between Ag top and (cell_top - small margin)
    z_mol_ref = z_top_ag + (cell_z - 2.0 - z_top_ag) / 2.0
    mol_ref   = mol_flat.copy()
    mol_ref[:, 0] += atop_xy[0]
    mol_ref[:, 1] += atop_xy[1]
    mol_ref[:, 2] += z_mol_ref
    mol_at = Atoms(symbols=mol_syms, positions=mol_ref,
                   cell=cell, pbc=[True, True, True])
    mol_at.calc = calc
    e_mol = float(mol_at.get_potential_energy())
    print(f"  E_mol    = {e_mol:.4f} eV  (z_ref={z_mol_ref:.2f} Å, same cell/pbc)")

    # Scan
    n_ag = len(ag_syms)
    z_list, e_int_list, fz_list = [], [], []
    print(f"\n  z-scan: {len(z_scan)} points  "
          f"[{z_scan[0]:.2f} → {z_scan[-1]:.2f} Å]")
    t0 = time.perf_counter()

    for i, z_rel in enumerate(z_scan):
        mol_z = mol_flat.copy()
        mol_z[:, 0] += atop_xy[0]
        mol_z[:, 1] += atop_xy[1]
        mol_z[:, 2] += z_top_ag + z_rel

        comb = Atoms(symbols=ag_syms + mol_syms,
                     positions=np.vstack([ag_pos, mol_z]),
                     cell=cell, pbc=[True, True, True])
        comb.calc = calc
        e_tot  = float(comb.get_potential_energy())
        f_all  = np.asarray(comb.get_forces(), dtype=float)
        f_mol  = f_all[n_ag:]

        e_int  = e_tot - e_sub - e_mol
        fz_int = float(f_mol[:, 2].sum())   # total z-force on molecule

        z_list.append(z_rel);  e_int_list.append(e_int);  fz_list.append(fz_int)

        if (i + 1) % 10 == 0 or (i + 1) == len(z_scan):
            print(f"    [{i+1:3d}/{len(z_scan)}] z={z_rel:.2f}Å  "
                  f"E_int={e_int*1000:+8.1f} meV  "
                  f"Fz={fz_int:+.4f} eV/Å  "
                  f"({time.perf_counter()-t0:.0f}s)")

    print(f"  Scan done in {time.perf_counter()-t0:.1f}s")
    return {
        "z_rel":    np.array(z_list),
        "e_int":    np.array(e_int_list),    # eV
        "fz_int":   np.array(fz_list),       # eV/Å  (sum over mol atoms)
        "e_sub":    e_sub,
        "e_mol":    e_mol,
        "z_top_ag": z_top_ag,
        "atop_xy":  atop_xy,
        "mol_syms": mol_syms,
    }


# ════════════════════════════════════════════════════════════════════════════
# 4.  PLQ field sampling  +  z-gradient
# ════════════════════════════════════════════════════════════════════════════

def sample_plq(coeffs, meta, positions):
    """Cubic B-spline: (N,3) positions → (N,3) [V_P, V_L, V_Q]."""
    from pyBall.gridff_jax.hybrid_energy.model import _sample_bspline_numpy
    from pyBall.gridff_jax.interfaces import DensityData
    stub = DensityData(
        cell=np.asarray(meta["cell"], dtype=float),
        origin=np.asarray(meta["origin"], dtype=float),
        voxel=meta["voxel"],
        symbols=[], positions=np.zeros((0, 3), dtype=float),
        charges=None, rho_zyx=None, v_loc_zyx=None,
        grid_shape_xyz_hint=meta["grid_xyz"],
    )
    r = _sample_bspline_numpy(coeffs, stub, positions)
    return r if r.ndim == 2 else r[:, None]   # (N, 3)


def sample_plq_grad_z(coeffs, meta, positions, dz=0.01):
    """Numerical z-gradient of PLQ fields: −∂V/∂z via central diff."""
    dv = np.zeros_like(positions); dv[:, 2] = dz
    vp = sample_plq(coeffs, meta, positions + dv)
    vm = sample_plq(coeffs, meta, positions - dv)
    return (vp - vm) / (2 * dz)    # (N, 3)  positive = field increases with z


def positions_in_plq_frame(mol_pos_global, z_top_md, meta):
    """Shift z so z_top_md aligns with PLQ z_ref; fold x,y into PLQ cell."""
    p   = mol_pos_global.copy().astype(float)
    p[:, 2] += float(meta["z_ref"]) - z_top_md   # z alignment

    # Fold x,y into PLQ cell (non-orthogonal PBC)
    c   = np.asarray(meta["cell"], dtype=float)
    a1  = c[0, :2];  a2 = c[1, :2]
    M   = np.column_stack([a1, a2])
    frc = (np.linalg.inv(M) @ p[:, :2].T).T % 1.0
    p[:, :2] = (M @ frc.T).T
    return p


# ════════════════════════════════════════════════════════════════════════════
# 5.  Design matrices  +  linear regression
# ════════════════════════════════════════════════════════════════════════════

def build_matrices(scan, coeffs, meta):
    """Build X_E (energy) and X_F (force) design matrices for all scan points.

    X_E[k, j*N_EL + el] = Σ_{i∈el}   V_j(r_i^k)          eV / coeff
    X_F[k, j*N_EL + el] = Σ_{i∈el}  −∂V_j/∂z(r_i^k)     eV/Å / coeff
    """
    z_arr    = scan["z_rel"]
    atop_xy  = scan["atop_xy"]
    z_top    = scan["z_top_ag"]
    mol_syms = scan["mol_syms"]
    mol_flat = scan["mol_flat"]
    n        = len(z_arr)

    X_E = np.zeros((n, 3 * N_EL), dtype=float)
    X_F = np.zeros((n, 3 * N_EL), dtype=float)

    for k, z_rel in enumerate(z_arr):
        mol_z = mol_flat.copy()
        mol_z[:, 0] += atop_xy[0]
        mol_z[:, 1] += atop_xy[1]
        mol_z[:, 2] += z_top + z_rel

        pos_plq   = positions_in_plq_frame(mol_z, z_top, meta)
        plq_vals  = sample_plq(       coeffs, meta, pos_plq)   # (N_mol, 3)
        plq_grad  = sample_plq_grad_z(coeffs, meta, pos_plq)   # (N_mol, 3)

        for i, sym in enumerate(mol_syms):
            el = EL_MAP.get(sym, -1)
            if el < 0:
                continue
            for j in range(3):
                X_E[k, j * N_EL + el] +=  plq_vals[i, j]
                X_F[k, j * N_EL + el] += -plq_grad[i, j]   # Fz = −∂E/∂z

    return X_E, X_F


def fit_plq(X_E, X_F, e_int, fz_int, force_weight=0.1):
    """Return theta_E (energy only) and theta_EF (joint energy+force)."""
    # energy-only
    theta_E, _, rank_E, _ = np.linalg.lstsq(X_E, e_int, rcond=None)

    # joint: stack energy rows + weighted force rows
    w    = force_weight
    X_jt = np.vstack([X_E,     w * X_F])
    y_jt = np.hstack([e_int,   w * fz_int])
    theta_EF, _, rank_EF, _ = np.linalg.lstsq(X_jt, y_jt, rcond=None)

    e_pred_E   = X_E @ theta_E
    e_pred_EF  = X_E @ theta_EF
    fz_pred_E  = X_F @ theta_E
    fz_pred_EF = X_F @ theta_EF

    rmse_E_e   = float(np.sqrt(np.mean((e_pred_E  - e_int )**2))) * 1000
    rmse_EF_e  = float(np.sqrt(np.mean((e_pred_EF - e_int )**2))) * 1000
    rmse_E_f   = float(np.sqrt(np.mean((fz_pred_E  - fz_int)**2)))
    rmse_EF_f  = float(np.sqrt(np.mean((fz_pred_EF - fz_int)**2)))

    print(f"\n  Fit results (rank E={rank_E}, EF={rank_EF}):")
    print(f"    theta_E  : RMSE_E = {rmse_E_e:.1f} meV   RMSE_F = {rmse_E_f:.4f} eV/Å")
    print(f"    theta_EF : RMSE_E = {rmse_EF_e:.1f} meV   RMSE_F = {rmse_EF_f:.4f} eV/Å")
    return {
        "theta_E":    theta_E,
        "theta_EF":   theta_EF,
        "e_E":        e_pred_E,
        "e_EF":       e_pred_EF,
        "fz_E":       fz_pred_E,
        "fz_EF":      fz_pred_EF,
        "rmse_E_e":   rmse_E_e,
        "rmse_EF_e":  rmse_EF_e,
        "rmse_E_f":   rmse_E_f,
        "rmse_EF_f":  rmse_EF_f,
    }


def print_coeffs(theta, mol_syms, label):
    els = sorted(set(mol_syms) & set(ELEMENTS))
    print(f"\n  {label}:")
    print(f"  {'El':>4}  {'A_Pauli':>14}  {'B_London':>14}  {'C_Coulomb':>14}")
    for el in els:
        i = EL_MAP[el]
        print(f"  {el:>4}  {theta[0*N_EL+i]:>14.5f}  "
              f"{theta[1*N_EL+i]:>14.5f}  {theta[2*N_EL+i]:>14.5f}")


# ════════════════════════════════════════════════════════════════════════════
# 6.  Plots
# ════════════════════════════════════════════════════════════════════════════

def _save(fig, path):
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {path}")


def plot_energy(z, e_mace, e_E, e_EF, fit, out_dir):
    """Plot 1: energy comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(z, e_mace * 1000, "k-o",  ms=5, lw=2,   label="MACE  E_int")
    ax.plot(z, e_E    * 1000, "b--",        lw=1.8,
            label=f"PLQ (E-fit)   RMSE={fit['rmse_E_e']:.0f} meV")
    ax.plot(z, e_EF   * 1000, "r-",         lw=1.8,
            label=f"PLQ (E+F-fit) RMSE={fit['rmse_EF_e']:.0f} meV")
    ax.axhline(0, color="k", lw=0.5, ls=":")
    ax.set_xlabel("z above Ag top layer (Å)", fontsize=12)
    ax.set_ylabel("E_int  (meV)",              fontsize=12)
    ax.set_title(
        "Interaction energy: MACE z-scan  vs  optimal PLQ fit\n"
        f"Ag slab: {CHGCAR}\n"
        f"MACE: {MACE_MDL}",
        fontsize=9)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.25)
    _save(fig, os.path.join(out_dir, "proof_plot1_energy.png"))


def plot_force(z, fz_mace, fz_E, fz_EF, fit, out_dir):
    """Plot 2: force comparison (Fz = sum of z-forces on all mol atoms)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(z, fz_mace, "k-o",  ms=5, lw=2,   label="MACE  Fz_int  (sum over mol)")
    ax.plot(z, fz_E,    "b--",        lw=1.8,
            label=f"PLQ (E-fit)   RMSE={fit['rmse_E_f']:.4f} eV/Å")
    ax.plot(z, fz_EF,   "r-",         lw=1.8,
            label=f"PLQ (E+F-fit) RMSE={fit['rmse_EF_f']:.4f} eV/Å")
    ax.axhline(0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("z above Ag top layer (Å)", fontsize=12)
    ax.set_ylabel("Fz  (eV/Å)",               fontsize=12)
    ax.set_title(
        "Interaction force z-component: MACE  vs  PLQ\n"
        "(Fz > 0 → pushed away from surface; Fz < 0 → attracted)",
        fontsize=10)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.25)
    _save(fig, os.path.join(out_dir, "proof_plot2_force.png"))


def plot_parity(e_mace, e_E, e_EF, fit, out_dir):
    """Plot 3: parity scatter."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, e_plq, label, rmse in [
        (axes[0], e_E,  "E-only fit",  fit["rmse_E_e"]),
        (axes[1], e_EF, "E+F joint",   fit["rmse_EF_e"]),
    ]:
        ym = e_mace * 1000;  yp = e_plq * 1000
        ax.scatter(ym, yp, s=30, color="steelblue", alpha=0.8)
        lo = min(ym.min(), yp.min()) - 20
        hi = max(ym.max(), yp.max()) + 20
        ax.plot([lo, hi], [lo, hi], "k--", lw=1)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("MACE E_int (meV)"); ax.set_ylabel("PLQ E_int (meV)")
        r2 = 1 - np.var(yp - ym) / max(np.var(ym), 1e-30)
        ax.set_title(f"{label}\nRMSE={rmse:.1f} meV   R²={r2:.4f}")
        ax.grid(True, alpha=0.25)
    fig.suptitle("PLQ linear fit parity  —  MACE as target", fontsize=12)
    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "proof_plot3_parity.png"))


def plot_residuals(z, e_mace, e_E, e_EF, fz_mace, fz_E, fz_EF, out_dir):
    """Plot 4: residuals vs z (energy and force)."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(z, (e_E  - e_mace) * 1000, "b--", lw=1.5, label="E-only fit")
    ax1.plot(z, (e_EF - e_mace) * 1000, "r-",  lw=1.5, label="E+F joint")
    ax1.axhline(0, color="k", lw=0.5, ls=":")
    ax1.set_ylabel("ΔE_int  (meV)")
    ax1.set_title("Residuals (PLQ − MACE)\n"
                  "Systematic pattern → PLQ functional form limitation")
    ax1.legend(); ax1.grid(True, alpha=0.25)

    ax2.plot(z, fz_E  - fz_mace, "b--", lw=1.5, label="E-only fit")
    ax2.plot(z, fz_EF - fz_mace, "r-",  lw=1.5, label="E+F joint")
    ax2.axhline(0, color="k", lw=0.5, ls=":")
    ax2.set_xlabel("z above Ag top layer (Å)")
    ax2.set_ylabel("ΔFz  (eV/Å)")
    ax2.legend(); ax2.grid(True, alpha=0.25)

    plt.tight_layout()
    _save(fig, os.path.join(out_dir, "proof_plot4_residuals.png"))


# ════════════════════════════════════════════════════════════════════════════
# 7.  Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device",          default="cpu")
    ap.add_argument("--z-min",   type=float, default=2.0)
    ap.add_argument("--z-max",   type=float, default=12.0)
    ap.add_argument("--z-step",  type=float, default=0.25)
    ap.add_argument("--force-recompute", action="store_true")
    ap.add_argument("--force-weight", type=float, default=0.1,
                    help="Weight for force rows in joint E+F regression (default 0.1)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    scan_cache = os.path.join(OUT_DIR, "_proof_mace_scan.npy")

    print("=" * 68)
    print("PROOF: MACE z-scan → PLQ fit → GridFF scan → comparison")
    print(f"  CHGCAR  : {CHGCAR}")
    print(f"  LOCPOT  : {LOCPOT}")
    print(f"  MACE    : {MACE_MDL}")
    print(f"  EXTXYZ  : {EXTXYZ}")
    print(f"  Out dir : {OUT_DIR}")
    print(f"  Device  : {args.device}")
    print("=" * 68)

    # ── 1. Substrate + PLQ ────────────────────────────────────────────────
    print("\n[1] Loading Ag slab + PLQ fields from CHGCAR/LOCPOT …")
    ag_syms, ag_pos, cell, coeffs, meta = load_substrate(CHGCAR, LOCPOT, OUT_DIR)

    # ── 2. Molecule ───────────────────────────────────────────────────────
    print("\n[2] Reading molecule from extxyz …")
    mol_syms, mol_pos = read_molecule(EXTXYZ)
    mol_flat = orient_flat(mol_pos)
    z_spread = float(np.ptp(mol_flat[:, 2]))
    print(f"  Flat orientation: z-spread = {z_spread:.3f} Å")

    # ── 3. MACE z-scan ────────────────────────────────────────────────────
    z_scan = np.arange(args.z_min,
                       args.z_max + args.z_step / 2,
                       args.z_step)

    if os.path.exists(scan_cache) and not args.force_recompute:
        print(f"\n[3] Loading cached MACE scan: {scan_cache}")
        scan = np.load(scan_cache, allow_pickle=True).item()
        print(f"  {len(scan['z_rel'])} points  "
              f"z=[{scan['z_rel'][0]:.2f}, {scan['z_rel'][-1]:.2f}] Å")
    else:
        print(f"\n[3] Running MACE z-scan  ({len(z_scan)} points) …")
        scan = mace_zscan(ag_syms, ag_pos, mol_syms, mol_flat,
                          cell, z_scan, MACE_MDL, args.device)
        scan["mol_flat"] = mol_flat          # store for design matrix
        np.save(scan_cache, scan, allow_pickle=True)  # type: ignore[arg-type]
        print(f"  Cached → {scan_cache}")

    if "mol_flat" not in scan:
        scan["mol_flat"] = mol_flat

    z      = scan["z_rel"]
    e_mace = scan["e_int"]
    fz_mace = scan["fz_int"]
    i_min   = int(np.argmin(e_mace))
    print(f"\n  MACE summary:")
    print(f"    E_int range : {e_mace.min()*1000:.1f} → {e_mace.max()*1000:.1f} meV")
    print(f"    E_int min   : {e_mace[i_min]*1000:.1f} meV  @ z = {z[i_min]:.2f} Å")
    print(f"    Fz range    : {fz_mace.min():.4f} → {fz_mace.max():.4f} eV/Å")

    # ── 4 & 5. Design matrices + regression ──────────────────────────────
    print("\n[4] Building PLQ design matrices X_E and X_F …")
    X_E, X_F = build_matrices(scan, coeffs, meta)
    print(f"  X_E shape : {X_E.shape}   rank = {np.linalg.matrix_rank(X_E)}")
    print(f"  X_F shape : {X_F.shape}   rank = {np.linalg.matrix_rank(X_F)}")

    print(f"\n[5] Linear regression  (force_weight={args.force_weight}) …")
    fit = fit_plq(X_E, X_F, e_mace, fz_mace, force_weight=args.force_weight)
    print_coeffs(fit["theta_E"],  mol_syms, "theta (E-only fit)")
    print_coeffs(fit["theta_EF"], mol_syms, "theta (E+F joint fit)")

    # ── 6. Plots ──────────────────────────────────────────────────────────
    print("\n[6] Generating comparison plots …")
    plot_energy(z, e_mace, fit["e_E"],  fit["e_EF"],  fit, OUT_DIR)
    plot_force( z, fz_mace, fit["fz_E"], fit["fz_EF"], fit, OUT_DIR)
    plot_parity(e_mace, fit["e_E"], fit["e_EF"], fit, OUT_DIR)
    plot_residuals(z, e_mace, fit["e_E"], fit["e_EF"],
                   fz_mace, fit["fz_E"], fit["fz_EF"], OUT_DIR)

    print(f"\nDone.  All outputs in {OUT_DIR}/")
    print("\n  proof_plot1_energy.png   – MACE vs PLQ_E vs PLQ_EF  energy vs z")
    print("  proof_plot2_force.png    – MACE vs PLQ_E vs PLQ_EF  Fz vs z")
    print("  proof_plot3_parity.png   – E_PLQ vs E_MACE scatter")
    print("  proof_plot4_residuals.png– (PLQ−MACE) energy + force vs z")


if __name__ == "__main__":
    main()
