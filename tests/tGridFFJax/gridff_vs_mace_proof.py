#!/usr/bin/env python3
"""GridFF vs MACE multi-site proof and functional-form diagnosis.

Uses the SAME Ag slab (CHGCAR) for both the PLQ field construction and the
MACE substrate.  Molecule from first MD_interface_Ag extxyz frame, flat orientation.

Pipeline
--------
1.  Load Ag slab + PLQ B-spline coefficients (from mace_vs_plq_zscan.py cache).
2.  Load molecule + atop-site MACE z-scan (from previous cache).
3.  Fit PLQ linearly to the atop MACE scan  → theta_atop.
4.  Find atop / bridge / FCC-hollow surface sites from Ag top layer.
5.  Run MACE z-scans at bridge and hollow (cached separately).
6.  Evaluate GridFF at ALL three sites using theta_atop (same parameters).
7.  Global re-fit PLQ to all three sites simultaneously → theta_global.
8.  Evaluate GridFF with theta_global.
9.  Diagnostic plots:
      proof2_plot1_rmse_table.png     – RMSE summary table, both fits
      proof2_plot2_multisite_E.png    – energy: MACE vs GridFF at 3 sites
      proof2_plot3_multisite_F.png    – Fz:     MACE vs GridFF at 3 sites
      proof2_plot4_residuals.png      – (GridFF - MACE) at all 3 sites
      proof2_plot5_decay.png          – log|E_int| vs z: PLQ vs MACE decay exponent
      proof2_plot6_global_fit.png     – improvement with global multi-site fit

Usage
-----
  python3 tests/tGridFFJax/gridff_vs_mace_proof.py
  python3 tests/tGridFFJax/gridff_vs_mace_proof.py --device cpu --force-recompute
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from collections import Counter

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_proof2")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

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
CACHE_COEFFS = os.path.join(OUT_DIR, "_Bspline_PLQ_coeffs.npy")
CACHE_META   = os.path.join(OUT_DIR, "_Bspline_PLQ_meta.npy")
CACHE_ATOP   = os.path.join(OUT_DIR, "_proof_mace_scan.npy")   # from prev script

for _p in [FIRECORE, MADSURF]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

ELEMENTS = ["C", "H", "N", "O"]
N_EL     = len(ELEMENTS)
EL_MAP   = {e: i for i, e in enumerate(ELEMENTS)}
_RE_CFG  = re.compile(r'config_type=(\S+)')


# ════════════════════════════════════════════════════════════════════════════
# I.  Load cached substrate, PLQ coefficients, molecule
# ════════════════════════════════════════════════════════════════════════════

def load_substrate_cached():
    """Load Ag slab geometry and PLQ B-spline coefficients from cache."""
    from pyBall.gridff_jax.config import DensityBackendConfig
    from pyBall.gridff_jax.density_backends.vasp_volumetric import VaspVolumetricBackend

    cfg     = DensityBackendConfig(kind="vasp_volumetric",
                                   chgcar_path=CHGCAR, locpot_path=LOCPOT)
    density = VaspVolumetricBackend(cfg).load()
    ag_syms = list(density.symbols)
    ag_pos  = np.asarray(density.positions, dtype=float)
    cell    = np.asarray(density.cell,      dtype=float)

    assert os.path.exists(CACHE_COEFFS), \
        f"PLQ cache not found: {CACHE_COEFFS}\nRun mace_vs_plq_zscan.py first."
    coeffs = np.load(CACHE_COEFFS)
    meta   = np.load(CACHE_META, allow_pickle=True).item()
    print(f"  Ag slab : {len(ag_syms)} atoms  "
          f"z=[{ag_pos[:,2].min():.3f}, {ag_pos[:,2].max():.3f}] Å")
    print(f"  PLQ     : shape {coeffs.shape}  z_ref={meta['z_ref']:.4f} Å")
    return ag_syms, ag_pos, cell, coeffs, meta


def read_molecule():
    """Extract non-Ag atoms from first MD_interface_Ag frame."""
    with open(EXTXYZ, "r", errors="ignore") as fh:
        while True:
            line = fh.readline()
            if not line: break
            try: n = int(line.strip())
            except ValueError: continue
            info  = fh.readline()
            atoms = [fh.readline() for _ in range(n)]
            m = _RE_CFG.search(info)
            if not m or m.group(1) != "MD_interface_Ag": continue
            syms, pos = [], []
            for a in atoms:
                p = a.split()
                if len(p) >= 4 and p[0] != "Ag":
                    syms.append(p[0]); pos.append([float(p[1]),float(p[2]),float(p[3])])
            print(f"  Molecule : {len(syms)} atoms  {dict(Counter(syms))}")
            return syms, np.array(pos, dtype=float)
    raise RuntimeError("no MD_interface_Ag frame")


def orient_flat(mol_pos):
    com = mol_pos.mean(axis=0)
    c   = mol_pos - com
    _, V = np.linalg.eigh(c.T @ c)
    R   = np.column_stack([V[:, 2], V[:, 1], V[:, 0]])
    fp  = c @ R
    if fp[:, 2].mean() < 0: fp[:, 2] *= -1
    return fp


# ════════════════════════════════════════════════════════════════════════════
# II.  Surface site geometry
# ════════════════════════════════════════════════════════════════════════════

def get_surface_sites(ag_pos, cell):
    """Compute atop / bridge / FCC-hollow xy coordinates from top Ag layer.

    Returns dict: {site_name: (x, y)} in Cartesian Å.
    """
    # Top layer: atoms within 0.3 Å of the maximum z
    top_mask = ag_pos[:, 2] > ag_pos[:, 2].max() - 0.3
    top_xy   = ag_pos[top_mask, :2].copy()
    cell_2d  = cell[:2, :2]           # 2×2 sub-matrix for 2D PBC

    def min_img_vec(dr):
        frac = np.linalg.solve(cell_2d.T, dr)
        frac -= np.round(frac)
        return cell_2d.T @ frac

    # atop: highest atom in the top layer (same as mace_vs_plq_zscan)
    i_top  = np.argmax(ag_pos[top_mask, 2])
    ref    = top_xy[i_top].copy()

    # nearest neighbours of ref with PBC
    nn_list = []
    for pos in top_xy:
        dr   = min_img_vec(pos - ref)
        dist = float(np.linalg.norm(dr))
        if dist > 0.1:                    # skip self
            nn_list.append((dist, ref + dr))
    nn_list.sort(key=lambda x: x[0])

    # bridge = midpoint between ref and nearest neighbor
    nn_dist = nn_list[0][0]
    nn1 = np.asarray(nn_list[0][1])
    bridge = (ref + nn1) / 2.0

    # FCC-hollow = centroid of ref and two ADJACENT (not opposite) NNs.
    # Adjacent means they also form a NN pair (|nn1 - nn2| ≈ nn_dist).
    best_nn2  = None
    best_dev  = np.inf
    for d2, nn2_raw in nn_list[1:]:
        if abs(d2 - nn_dist) > 0.15:
            continue                      # not a NN of ref
        nn2 = np.asarray(nn2_raw)
        d12 = float(np.linalg.norm(nn2 - nn1))   # direct distance (already PBC-wrapped)
        dev = abs(d12 - nn_dist)
        if dev < best_dev:
            best_dev  = dev
            best_nn2  = nn2
    hollow = (ref + nn1 + best_nn2) / 3.0

    sites = {"atop": ref, "bridge": bridge, "fcc_hollow": hollow}
    print(f"  Surface sites (xy):")
    for name, xy in sites.items():
        print(f"    {name:12s}  x={xy[0]:.3f} y={xy[1]:.3f} Å")
    return sites


# ════════════════════════════════════════════════════════════════════════════
# III.  MACE z-scan at a given lateral site
# ════════════════════════════════════════════════════════════════════════════

def mace_zscan_site(ag_syms, ag_pos, mol_syms, mol_flat,
                    cell, z_scan, site_xy, site_name,
                    z_top_ag, e_sub, e_mol, calc, cache_dir):
    """Run MACE z-scan at given xy site.  Uses pre-computed E_sub and E_mol."""
    from ase import Atoms

    cache = os.path.join(cache_dir, f"_proof2_mace_{site_name}.npy")
    if os.path.exists(cache):
        d = np.load(cache, allow_pickle=True).item()
        if len(d["z_rel"]) == len(z_scan) and np.allclose(d["z_rel"], z_scan, atol=1e-6):
            print(f"    [{site_name}] loaded from cache: {cache}")
            return d

    n_ag   = len(ag_syms)
    z_list, e_list, fz_list = [], [], []
    t0 = time.perf_counter()
    for i, z_rel in enumerate(z_scan):
        mol_z       = mol_flat.copy()
        mol_z[:, 0] += site_xy[0]
        mol_z[:, 1] += site_xy[1]
        mol_z[:, 2] += z_top_ag + z_rel
        comb = Atoms(symbols=ag_syms + mol_syms,
                     positions=np.vstack([ag_pos, mol_z]),
                     cell=cell, pbc=[True, True, True])
        comb.calc = calc
        e_tot  = float(comb.get_potential_energy())
        f_all  = np.asarray(comb.get_forces(), dtype=float)
        z_list.append(z_rel)
        e_list.append(e_tot - e_sub - e_mol)
        fz_list.append(float(f_all[n_ag:, 2].sum()))
        if (i + 1) % 10 == 0 or (i + 1) == len(z_scan):
            print(f"    [{site_name}] [{i+1:3d}/{len(z_scan)}] "
                  f"z={z_rel:.2f}Å  E={e_list[-1]*1000:+8.1f} meV  "
                  f"({time.perf_counter()-t0:.0f}s)")

    result = {"z_rel":  np.array(z_list),
              "e_int":  np.array(e_list),
              "fz_int": np.array(fz_list),
              "site":   site_name}
    np.save(cache, result, allow_pickle=True)
    return result


# ════════════════════════════════════════════════════════════════════════════
# IV.  PLQ design matrix and regression
# ════════════════════════════════════════════════════════════════════════════

def sample_plq(coeffs, meta, positions):
    """B-spline sample: (N,3) pos → (N,3) [V_P, V_L, V_Q]."""
    from pyBall.gridff_jax.hybrid_energy.model import _sample_bspline_numpy
    from pyBall.gridff_jax.interfaces import DensityData
    stub = DensityData(
        cell=np.asarray(meta["cell"], dtype=float),
        origin=np.asarray(meta["origin"], dtype=float),
        voxel=meta["voxel"], symbols=[], positions=np.zeros((0,3),dtype=float),
        charges=None, rho_zyx=None, v_loc_zyx=None,
        grid_shape_xyz_hint=meta["grid_xyz"])
    r = _sample_bspline_numpy(coeffs, stub, positions)
    return r if r.ndim == 2 else r[:, None]


def sample_plq_grad_z(coeffs, meta, positions, dz=0.01):
    dv = np.zeros_like(positions); dv[:, 2] = dz
    return (sample_plq(coeffs, meta, positions + dv) -
            sample_plq(coeffs, meta, positions - dv)) / (2 * dz)


def positions_in_plq_frame(mol_pos_global, z_top_md, meta):
    p   = mol_pos_global.copy().astype(float)
    p[:, 2] += float(meta["z_ref"]) - z_top_md
    c   = np.asarray(meta["cell"], dtype=float)
    a1  = c[0, :2]; a2 = c[1, :2]
    M   = np.column_stack([a1, a2])
    frc = (np.linalg.inv(M) @ p[:, :2].T).T % 1.0
    p[:, :2] = (M @ frc.T).T
    return p


def build_design_matrices(mol_syms, mol_flat, site_xy, z_scan,
                          z_top_ag, coeffs, meta):
    """Build X_E and X_F for a single scan site."""
    n = len(z_scan)
    X_E = np.zeros((n, 3 * N_EL), dtype=float)
    X_F = np.zeros((n, 3 * N_EL), dtype=float)
    for k, z_rel in enumerate(z_scan):
        mol_z       = mol_flat.copy()
        mol_z[:, 0] += site_xy[0]
        mol_z[:, 1] += site_xy[1]
        mol_z[:, 2] += z_top_ag + z_rel
        pp = positions_in_plq_frame(mol_z, z_top_ag, meta)
        vv = sample_plq(      coeffs, meta, pp)   # (N_mol, 3)
        gv = sample_plq_grad_z(coeffs, meta, pp)  # (N_mol, 3)
        for i, sym in enumerate(mol_syms):
            el = EL_MAP.get(sym, -1)
            if el < 0: continue
            for j in range(3):
                X_E[k, j*N_EL + el] +=  vv[i, j]
                X_F[k, j*N_EL + el] += -gv[i, j]
    return X_E, X_F


def fit_plq(X_E_list, X_F_list, e_int_list, fz_list,
            force_weight=0.1, label=""):
    """Fit PLQ parameters to multiple sites simultaneously."""
    X_E_all  = np.vstack(X_E_list)
    X_F_all  = np.vstack(X_F_list)
    e_all    = np.hstack(e_int_list)
    fz_all   = np.hstack(fz_list)

    # Energy-only fit
    theta_E, _, rank_E, _ = np.linalg.lstsq(X_E_all, e_all, rcond=None)

    # Joint E+F fit
    w   = force_weight
    Xjt = np.vstack([X_E_all, w * X_F_all])
    yjt = np.hstack([e_all,   w * fz_all])
    theta_EF, _, rank_EF, _ = np.linalg.lstsq(Xjt, yjt, rcond=None)

    def rmse_e(th): return float(np.sqrt(np.mean((X_E_all@th - e_all)**2)))*1e3
    def rmse_f(th): return float(np.sqrt(np.mean((X_F_all@th - fz_all)**2)))
    def r2(th):
        res = X_E_all@th - e_all
        return 1 - float(np.var(res)) / max(float(np.var(e_all)), 1e-30)

    stats = {
        "E":  {"rmse_e": rmse_e(theta_E),  "rmse_f": rmse_f(theta_E),
               "r2": r2(theta_E),  "rank": rank_E,  "theta": theta_E},
        "EF": {"rmse_e": rmse_e(theta_EF), "rmse_f": rmse_f(theta_EF),
               "r2": r2(theta_EF), "rank": rank_EF, "theta": theta_EF},
    }
    print(f"\n  [{label or 'fit'}]  rank E={rank_E} EF={rank_EF}")
    for tag, s in stats.items():
        print(f"    theta_{tag}  RMSE_E={s['rmse_e']:.1f} meV  "
              f"RMSE_F={s['rmse_f']:.4f} eV/Å  R²={s['r2']:.4f}")
    return stats


# ════════════════════════════════════════════════════════════════════════════
# V.  GridFF z-scan (using fitted theta + B-spline coefficients)
# ════════════════════════════════════════════════════════════════════════════

def gridff_zscan(mol_syms, mol_flat, site_xy, z_scan,
                 z_top_ag, theta, coeffs, meta):
    """Evaluate PLQ GridFF energy and Fz at each z in z_scan."""
    X_E, X_F = build_design_matrices(mol_syms, mol_flat, site_xy, z_scan,
                                     z_top_ag, coeffs, meta)
    e_plq  = X_E @ theta
    fz_plq = X_F @ theta
    return e_plq, fz_plq


# ════════════════════════════════════════════════════════════════════════════
# VI.  Exponential decay analysis  (functional-form diagnosis)
# ════════════════════════════════════════════════════════════════════════════

def fit_exponential_tail(z, e_int, z_min_fit=3.5):
    """Fit |E_int| = A * exp(-alpha * z) to the tail.
    Returns (alpha, A, r2) or None if too few points."""
    mask = (z >= z_min_fit) & (np.abs(e_int) > 1e-6)
    if mask.sum() < 4:
        return None
    zf = z[mask]; ef = np.abs(e_int[mask])
    X  = np.column_stack([np.ones_like(zf), zf])
    c, _, _, _ = np.linalg.lstsq(X, np.log(ef), rcond=None)
    log_A, neg_alpha = c
    alpha = -neg_alpha
    pred  = np.exp(log_A + neg_alpha * zf)
    r2    = 1 - float(np.var(pred - ef)) / max(float(np.var(ef)), 1e-30)
    return {"alpha": float(alpha), "A": float(np.exp(log_A)), "r2": float(r2)}


def print_coeff_table(theta, mol_syms, label):
    els = sorted(set(mol_syms) & set(ELEMENTS))
    print(f"\n  {label}:")
    print(f"  {'El':>4}  {'A_Pauli':>14}  {'B_London':>14}  {'C_Coulomb':>14}")
    for el in els:
        i = EL_MAP[el]
        print(f"  {el:>4}  {theta[0*N_EL+i]:>14.5f}  "
              f"{theta[1*N_EL+i]:>14.5f}  {theta[2*N_EL+i]:>14.5f}")


# ════════════════════════════════════════════════════════════════════════════
# VII.  Plots
# ════════════════════════════════════════════════════════════════════════════

COLORS = {"atop": "#1f77b4", "bridge": "#ff7f0e", "fcc_hollow": "#2ca02c"}
SITE_LABELS = {"atop": "atop", "bridge": "bridge", "fcc_hollow": "FCC hollow"}


def _save(fig, name, out_dir):
    p = os.path.join(out_dir, name)
    fig.savefig(p, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {p}")


def plot_rmse_table(stats_atop, stats_global, mol_syms, theta_atop_EF, theta_glob_EF, out_dir):
    """Plot 1: RMSE summary table + coefficient bar chart."""
    fig = plt.figure(figsize=(14, 7))
    gs  = GridSpec(1, 2, figure=fig, wspace=0.4)

    # ── left: RMSE table ────────────────────────────────────────────────
    ax0 = fig.add_subplot(gs[0])
    ax0.axis("off")
    rows = [
        ["Fit",         "Sites", "RMSE_E (meV)", "RMSE_F (eV/Å)", "R²"],
        ["E-only atop", "atop",
         f"{stats_atop['E']['rmse_e']:.1f}",
         f"{stats_atop['E']['rmse_f']:.4f}",
         f"{stats_atop['E']['r2']:.4f}"],
        ["E+F atop",    "atop",
         f"{stats_atop['EF']['rmse_e']:.1f}",
         f"{stats_atop['EF']['rmse_f']:.4f}",
         f"{stats_atop['EF']['r2']:.4f}"],
        ["E-only global", "3 sites",
         f"{stats_global['E']['rmse_e']:.1f}",
         f"{stats_global['E']['rmse_f']:.4f}",
         f"{stats_global['E']['r2']:.4f}"],
        ["E+F global",  "3 sites",
         f"{stats_global['EF']['rmse_e']:.1f}",
         f"{stats_global['EF']['rmse_f']:.4f}",
         f"{stats_global['EF']['r2']:.4f}"],
    ]
    tbl = ax0.table(cellText=rows[1:], colLabels=rows[0],
                    loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(10)
    tbl.scale(1.2, 2.0)
    ax0.set_title("PLQ Fit Quality Summary", fontsize=12, pad=20)

    # ── right: coefficient bar chart ───────────────────────────────────
    ax1 = fig.add_subplot(gs[1])
    els  = sorted(set(mol_syms) & set(ELEMENTS))
    n_el = len(els)
    x    = np.arange(n_el)
    w    = 0.18
    field_labels = ["Pauli A", "London B", "Coulomb C"]
    offsets = [-w, 0, w]
    for fi, (off, flabel) in enumerate(zip(offsets, field_labels)):
        vals_atop = [theta_atop_EF[fi*N_EL + EL_MAP[e]] for e in els]
        vals_glob = [theta_glob_EF[fi*N_EL + EL_MAP[e]] for e in els]
        ax1.bar(x + off - 0.5*w, vals_atop, w*0.45, label=f"{flabel} (atop)", alpha=0.8)
        ax1.bar(x + off + 0.5*w, vals_glob, w*0.45, label=f"{flabel} (global)", alpha=0.5,
                hatch="/")
    ax1.set_xticks(x); ax1.set_xticklabels(els)
    ax1.set_xlabel("Element"); ax1.set_ylabel("Coefficient")
    ax1.set_title("PLQ Coefficients: atop-fit vs global-fit\n(E+F joint)")
    ax1.legend(fontsize=7, ncol=2); ax1.grid(True, alpha=0.2, axis="y")
    ax1.set_yscale("symlog", linthresh=1.0)

    _save(fig, "proof2_plot1_rmse_table.png", out_dir)


def plot_multisite_energy(scans, gridff_atop, gridff_glob, out_dir):
    """Plot 2: MACE vs GridFF energy at 3 sites."""
    sites = list(scans.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for ax, site in zip(axes, sites):
        zm      = scans[site]["z_rel"]
        em      = scans[site]["e_int"] * 1000
        eg_atop = gridff_atop[site][0] * 1000   # energy
        eg_glob = gridff_glob[site][0] * 1000
        ax.plot(zm, em,      "k-o", ms=5, lw=2, label="MACE (target)")
        ax.plot(zm, eg_atop, "b--", lw=1.8,
                label=f"GridFF atop-fit\nRMSE={float(np.sqrt(np.mean((eg_atop-em)**2))):.0f} meV")
        ax.plot(zm, eg_glob, "r-",  lw=1.8,
                label=f"GridFF global-fit\nRMSE={float(np.sqrt(np.mean((eg_glob-em)**2))):.0f} meV")
        ax.axhline(0, color="k", lw=0.4, ls=":")
        ax.set_xlabel("z above Ag top layer (Å)")
        ax.set_ylabel("E_int (meV)")
        ax.set_title(f"Site: {SITE_LABELS[site]}")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    fig.suptitle("Interaction energy: MACE vs PLQ/GridFF   (same θ for all sites)",
                 fontsize=12)
    plt.tight_layout()
    _save(fig, "proof2_plot2_multisite_E.png", out_dir)


def plot_multisite_force(scans, gridff_atop, gridff_glob, out_dir):
    """Plot 3: MACE vs GridFF Fz at 3 sites."""
    sites = list(scans.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    for ax, site in zip(axes, sites):
        zm      = scans[site]["z_rel"]
        fm      = scans[site]["fz_int"]
        fg_atop = gridff_atop[site][1]
        fg_glob = gridff_glob[site][1]
        ax.plot(zm, fm,      "k-o", ms=5, lw=2, label="MACE Fz")
        ax.plot(zm, fg_atop, "b--", lw=1.8,
                label=f"GridFF atop-fit\nRMSE={float(np.sqrt(np.mean((fg_atop-fm)**2))):.4f} eV/Å")
        ax.plot(zm, fg_glob, "r-",  lw=1.8,
                label=f"GridFF global-fit\nRMSE={float(np.sqrt(np.mean((fg_glob-fm)**2))):.4f} eV/Å")
        ax.axhline(0, color="k", lw=0.7, ls="--")
        ax.set_xlabel("z above Ag top layer (Å)")
        ax.set_ylabel("Fz (eV/Å)")
        ax.set_title(f"Site: {SITE_LABELS[site]}")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.2)
    fig.suptitle("z-force on molecule: MACE vs PLQ/GridFF", fontsize=12)
    plt.tight_layout()
    _save(fig, "proof2_plot3_multisite_F.png", out_dir)


def plot_residuals(scans, gridff_glob, out_dir):
    """Plot 4: (GridFF_global - MACE) energy and force residuals."""
    sites = list(scans.keys())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for site in sites:
        zm = scans[site]["z_rel"]
        de = (gridff_glob[site][0] - scans[site]["e_int"]) * 1000
        df = gridff_glob[site][1] - scans[site]["fz_int"]
        c  = COLORS[site]
        ax1.plot(zm, de, "-o", ms=4, lw=1.5, color=c, label=SITE_LABELS[site])
        ax2.plot(zm, df, "-o", ms=4, lw=1.5, color=c, label=SITE_LABELS[site])

    ax1.axhline(0, color="k", lw=0.5, ls=":")
    ax1.set_ylabel("ΔE_int  (GridFF − MACE)  [meV]")
    ax1.set_title("Residuals: global PLQ fit − MACE\n"
                  "Systematic pattern ⇒ PLQ functional form insufficient")
    ax1.legend(); ax1.grid(True, alpha=0.2)

    ax2.axhline(0, color="k", lw=0.5, ls=":")
    ax2.set_xlabel("z above Ag top layer (Å)")
    ax2.set_ylabel("ΔFz  (GridFF − MACE)  [eV/Å]")
    ax2.legend(); ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    _save(fig, "proof2_plot4_residuals.png", out_dir)


def plot_decay_analysis(scans, out_dir):
    """Plot 5: log|E_int| vs z — compare exponential decay constants."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    decay_info = {}
    for site in scans:
        zm = scans[site]["z_rel"]
        em = scans[site]["e_int"]
        c  = COLORS[site]

        # log-linear plot
        valid = np.abs(em) > 1e-6
        ax1.semilogy(zm[valid], np.abs(em[valid])*1000, "-o", ms=4, lw=1.5,
                     color=c, label=SITE_LABELS[site])

        # exponential fit to tail
        fit = fit_exponential_tail(zm, em, z_min_fit=3.5)
        if fit:
            z_tail = np.linspace(3.5, zm.max(), 60)
            ax1.semilogy(z_tail, fit["A"] * np.exp(-fit["alpha"] * z_tail) * 1000,
                         "--", lw=1.2, color=c, alpha=0.7,
                         label=f"  α={fit['alpha']:.2f} Å⁻¹  R²={fit['r2']:.3f}")
            decay_info[site] = fit

    ax1.set_xlabel("z above Ag top layer (Å)")
    ax1.set_ylabel("|E_int| (meV, log scale)")
    ax1.set_title("Exponential decay analysis\n"
                  "Dashed: fit to z > 3.5 Å tail")
    ax1.legend(fontsize=9); ax1.grid(True, alpha=0.2, which="both")

    # ── right: residuals from exponential fit ──────────────────────────
    for site, fit in decay_info.items():
        zm = scans[site]["z_rel"]
        em = np.abs(scans[site]["e_int"])
        e_exp = fit["A"] * np.exp(-fit["alpha"] * zm)
        rel_resid = (em - e_exp) / (np.abs(em) + 1e-8)
        ax2.plot(zm, rel_resid * 100, "-o", ms=4, lw=1.5,
                 color=COLORS[site], label=SITE_LABELS[site])

    ax2.axhline(0, color="k", lw=0.5, ls=":")
    ax2.set_xlabel("z above Ag top layer (Å)")
    ax2.set_ylabel("Relative deviation from pure exponential (%)")
    ax2.set_title("Non-exponential character of MACE E_int\n"
                  "> 20% deviation ⇒ PLQ single-exponential insufficient")
    ax2.legend(); ax2.grid(True, alpha=0.2)
    ax2.set_ylim(-150, 150)

    # Print decay constants
    print("\n  Exponential decay constants (E_int ~ A·exp(−α·z), fit z > 3.5 Å):")
    for site, fit in decay_info.items():
        print(f"    {SITE_LABELS[site]:12s}  α={fit['alpha']:.3f} Å⁻¹  "
              f"A={fit['A']*1000:.1f} meV  R²={fit['r2']:.4f}")
    print("  (PLQ Pauli field α_Morse ≈ 1.5 Å⁻¹, London ≈ 0.75 Å⁻¹)")

    plt.tight_layout()
    _save(fig, "proof2_plot5_decay.png", out_dir)


def plot_global_fit_comparison(scans, gridff_atop, gridff_glob, out_dir):
    """Plot 6: atop-fit vs global-fit — shows benefit of multi-site training."""
    sites  = list(scans.keys())
    n_sites = len(sites)
    fig, axes = plt.subplots(2, n_sites, figsize=(5*n_sites, 9), sharex=True)

    for col, site in enumerate(sites):
        zm = scans[site]["z_rel"]
        em = scans[site]["e_int"] * 1000
        e_atop = gridff_atop[site][0] * 1000
        e_glob = gridff_glob[site][0] * 1000
        res_atop = e_atop - em
        res_glob = e_glob - em

        # top row: energy
        axes[0, col].plot(zm, em,     "k-o", ms=4, lw=2,   label="MACE")
        axes[0, col].plot(zm, e_atop, "b--", lw=1.5, label="atop-fit PLQ")
        axes[0, col].plot(zm, e_glob, "r-",  lw=1.5, label="global PLQ")
        axes[0, col].axhline(0, color="k", lw=0.4, ls=":")
        axes[0, col].set_title(f"Site: {SITE_LABELS[site]}")
        axes[0, col].set_ylabel("E_int (meV)")
        axes[0, col].legend(fontsize=8); axes[0, col].grid(True, alpha=0.2)

        # bottom row: residuals
        rmse_a = float(np.sqrt(np.mean(res_atop**2)))
        rmse_g = float(np.sqrt(np.mean(res_glob**2)))
        axes[1, col].plot(zm, res_atop, "b--", lw=1.5,
                          label=f"atop-fit  RMSE={rmse_a:.1f} meV")
        axes[1, col].plot(zm, res_glob, "r-",  lw=1.5,
                          label=f"global    RMSE={rmse_g:.1f} meV")
        axes[1, col].axhline(0, color="k", lw=0.5, ls=":")
        axes[1, col].set_xlabel("z above Ag top layer (Å)")
        axes[1, col].set_ylabel("ΔE (meV)")
        axes[1, col].legend(fontsize=8); axes[1, col].grid(True, alpha=0.2)

    fig.suptitle("atop-fit PLQ  vs  global multi-site PLQ\n"
                 "If global ≈ atop → PLQ form is transferable.\n"
                 "If global ≫ atop at other sites → need more functional richness.",
                 fontsize=11)
    plt.tight_layout()
    _save(fig, "proof2_plot6_global_fit.png", out_dir)


# ════════════════════════════════════════════════════════════════════════════
# VIII.  Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device",          default="cpu")
    ap.add_argument("--z-min",   type=float, default=2.0)
    ap.add_argument("--z-max",   type=float, default=8.0,
                    help="Max z for bridge/hollow scans (shorter than atop to save time)")
    ap.add_argument("--z-step",  type=float, default=0.25)
    ap.add_argument("--force-recompute", action="store_true")
    ap.add_argument("--force-weight", type=float, default=0.1)
    ap.add_argument("--skip-new-mace", action="store_true",
                    help="Skip MACE at bridge/hollow (plot atop only)")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    print("=" * 68)
    print("GridFF vs MACE multi-site proof  (PLQ functional form diagnosis)")
    print(f"  CHGCAR  : {CHGCAR}")
    print(f"  LOCPOT  : {LOCPOT}")
    print(f"  MACE    : {MACE_MDL}")
    print(f"  Device  : {args.device}")
    print("=" * 68)

    # ── 1. Substrate ──────────────────────────────────────────────────────
    print("\n[1] Loading Ag slab + PLQ fields …")
    ag_syms, ag_pos, cell, coeffs, meta = load_substrate_cached()

    # ── 2. Molecule ───────────────────────────────────────────────────────
    print("\n[2] Molecule …")
    mol_syms, mol_pos = read_molecule()
    mol_flat = orient_flat(mol_pos)
    print(f"  z-spread flat: {float(np.ptp(mol_flat[:,2])):.3f} Å")

    # ── 3. Load atop MACE scan ────────────────────────────────────────────
    print(f"\n[3] Loading atop MACE scan from cache: {CACHE_ATOP}")
    assert os.path.exists(CACHE_ATOP), \
        f"Atop cache missing: run mace_vs_plq_zscan.py first."
    scan_atop = np.load(CACHE_ATOP, allow_pickle=True).item()
    if "mol_flat" not in scan_atop:
        scan_atop["mol_flat"] = mol_flat
    z_top_ag = float(scan_atop["z_top_ag"])
    atop_xy  = scan_atop["atop_xy"]
    e_sub    = float(scan_atop["e_sub"])
    e_mol    = float(scan_atop["e_mol"])
    z_atop   = scan_atop["z_rel"]
    print(f"  {len(z_atop)} points  z=[{z_atop[0]:.2f},{z_atop[-1]:.2f}] Å  "
          f"E_min={scan_atop['e_int'].min()*1000:.1f} meV")

    # ── 4. Surface sites ──────────────────────────────────────────────────
    print("\n[4] Computing surface sites …")
    sites = get_surface_sites(ag_pos, cell)
    # Override atop_xy with the exact value from MACE scan
    sites["atop"] = atop_xy

    # ── 5. MACE at bridge / hollow  ───────────────────────────────────────
    z_scan_new = np.arange(args.z_min, args.z_max + args.z_step/2, args.z_step)

    scans = {"atop": {"z_rel": z_atop,
                      "e_int": scan_atop["e_int"],
                      "fz_int": scan_atop["fz_int"]}}

    if not args.skip_new_mace:
        print(f"\n[5] MACE bridge / hollow scans  ({len(z_scan_new)} pts each) …")
        from mace.calculators import MACECalculator
        print(f"  Loading MACE: {MACE_MDL}")
        calc = MACECalculator(model_paths=MACE_MDL,
                              device=args.device, default_dtype="float64")
        for site_name in ("bridge", "fcc_hollow"):
            if args.force_recompute:
                cache_path = os.path.join(OUT_DIR, f"_proof2_mace_{site_name}.npy")
                if os.path.exists(cache_path): os.remove(cache_path)
            scans[site_name] = mace_zscan_site(
                ag_syms, ag_pos, mol_syms, mol_flat,
                cell, z_scan_new, sites[site_name], site_name,
                z_top_ag, e_sub, e_mol, calc, OUT_DIR)
    else:
        print("\n[5] Skipping bridge/hollow MACE (--skip-new-mace)")
        # Load from cache if available
        for site_name in ("bridge", "fcc_hollow"):
            c = os.path.join(OUT_DIR, f"_proof2_mace_{site_name}.npy")
            if os.path.exists(c):
                scans[site_name] = np.load(c, allow_pickle=True).item()
                print(f"  Loaded cached {site_name}")

    # Align all scans to the same z range for global fit
    # Use z_scan_new for bridge/hollow, z_atop for atop
    # For design matrices, we use each site's own z array

    # ── 6. Fit PLQ to atop site ───────────────────────────────────────────
    print("\n[6] PLQ fit — atop site only …")
    X_E_atop, X_F_atop = build_design_matrices(
        mol_syms, mol_flat, sites["atop"],
        z_atop, z_top_ag, coeffs, meta)
    stats_atop = fit_plq([X_E_atop], [X_F_atop],
                         [scan_atop["e_int"]], [scan_atop["fz_int"]],
                         force_weight=args.force_weight, label="atop")
    theta_atop_EF = stats_atop["EF"]["theta"]
    print_coeff_table(theta_atop_EF, mol_syms, "theta_EF (atop-fit)")

    # ── 7. Design matrices for other sites ────────────────────────────────
    X_E_list = [X_E_atop]
    X_F_list = [X_F_atop]
    e_list   = [scan_atop["e_int"]]
    fz_list  = [scan_atop["fz_int"]]

    for site_name in ("bridge", "fcc_hollow"):
        if site_name not in scans: continue
        zs = scans[site_name]["z_rel"]
        XE, XF = build_design_matrices(
            mol_syms, mol_flat, sites[site_name],
            zs, z_top_ag, coeffs, meta)
        X_E_list.append(XE); X_F_list.append(XF)
        e_list.append(scans[site_name]["e_int"])
        fz_list.append(scans[site_name]["fz_int"])

    # ── 8. Global multi-site fit ──────────────────────────────────────────
    print("\n[7] PLQ fit — global multi-site …")
    stats_global = fit_plq(X_E_list, X_F_list, e_list, fz_list,
                           force_weight=args.force_weight, label="global")
    theta_glob_EF = stats_global["EF"]["theta"]
    print_coeff_table(theta_glob_EF, mol_syms, "theta_EF (global fit)")

    # ── 9. GridFF scans at all sites with both theta sets ─────────────────
    print("\n[8] GridFF scans at all sites …")
    gridff_atop_results = {}
    gridff_glob_results = {}
    for site_name, scan in scans.items():
        zs = scan["z_rel"]
        ea, fa = gridff_zscan(mol_syms, mol_flat, sites[site_name],
                              zs, z_top_ag, theta_atop_EF, coeffs, meta)
        eg, fg = gridff_zscan(mol_syms, mol_flat, sites[site_name],
                              zs, z_top_ag, theta_glob_EF, coeffs, meta)
        gridff_atop_results[site_name] = (ea, fa)
        gridff_glob_results[site_name] = (eg, fg)

        rmse_e_a = float(np.sqrt(np.mean((ea - scan["e_int"])**2))) * 1000
        rmse_e_g = float(np.sqrt(np.mean((eg - scan["e_int"])**2))) * 1000
        print(f"  {SITE_LABELS[site_name]:12s}  "
              f"RMSE atop-θ = {rmse_e_a:.1f} meV   "
              f"RMSE global-θ = {rmse_e_g:.1f} meV")

    # ── 10. Plots ─────────────────────────────────────────────────────────
    print("\n[9] Generating diagnostic plots …")
    plot_rmse_table(stats_atop, stats_global, mol_syms,
                    theta_atop_EF, theta_glob_EF, OUT_DIR)
    plot_multisite_energy(scans, gridff_atop_results, gridff_glob_results, OUT_DIR)
    plot_multisite_force(scans, gridff_atop_results, gridff_glob_results, OUT_DIR)
    plot_residuals(scans, gridff_glob_results, OUT_DIR)
    plot_decay_analysis(scans, OUT_DIR)
    plot_global_fit_comparison(scans, gridff_atop_results,
                               gridff_glob_results, OUT_DIR)

    # ── Summary diagnosis ─────────────────────────────────────────────────
    print("\n" + "=" * 68)
    print("DIAGNOSIS SUMMARY")
    print("=" * 68)
    atop_rmse  = stats_atop["EF"]["rmse_e"]
    global_rmse = stats_global["EF"]["rmse_e"]
    print(f"\n  PLQ fit quality (E+F joint):")
    print(f"    atop only   : {atop_rmse:.1f} meV  "
          f"(force RMSE = {stats_atop['EF']['rmse_f']:.4f} eV/Å)")
    print(f"    global      : {global_rmse:.1f} meV  "
          f"(force RMSE = {stats_global['EF']['rmse_f']:.4f} eV/Å)")

    if len(scans) >= 3:
        # Per-site RMSE with global theta
        site_rmses = {}
        for i, site_name in enumerate(scans.keys()):
            pred = X_E_list[i] @ theta_glob_EF
            site_rmses[site_name] = float(np.sqrt(np.mean(
                (pred - e_list[i])**2))) * 1000
        print(f"\n  Per-site RMSE with global theta:")
        for s, rmse in site_rmses.items():
            print(f"    {SITE_LABELS[s]:12s}  {rmse:.1f} meV")

        worst = max(site_rmses.values())
        print(f"\n  Worst-site RMSE = {worst:.1f} meV")
        if worst < 50:
            print("  ✓  PLQ functional form is SUFFICIENT for this molecule/surface.")
            print("     Recommendation: use global multi-site fit as production theta.")
        elif worst < 150:
            print("  ⚠  PLQ is adequate but limited.  Consider adding:")
            print("     • More lateral sites in fit (bridge, hollow, bridge2)")
            print("     • Per-element alpha_morse optimization")
            print("     • Reactive channel for N-Ag bonding")
        else:
            print("  ✗  PLQ functional form is INSUFFICIENT.")
            print("     Root cause (check plot5 decay analysis):")
            print("     • If α_MACE ≫ 1.5 Å⁻¹: Pauli repulsion too soft → need larger alpha")
            print("     • If residuals are site-dependent: corrugation mismatch")
            print("     • If residuals grow at z < 2.5 Å: hard-core needs Morse 2-body term")
    else:
        print("\n  (Only atop site computed.  Run without --skip-new-mace for full diagnosis.)")

    print(f"\n  All plots in: {OUT_DIR}/")
    print("  proof2_plot1_rmse_table.png    – RMSE summary + coefficient bars")
    print("  proof2_plot2_multisite_E.png   – energy: MACE vs GridFF at 3 sites")
    print("  proof2_plot3_multisite_F.png   – Fz: MACE vs GridFF at 3 sites")
    print("  proof2_plot4_residuals.png     – (GridFF − MACE) residuals")
    print("  proof2_plot5_decay.png         – exponential decay diagnosis")
    print("  proof2_plot6_global_fit.png    – atop-fit vs global-fit comparison")


if __name__ == "__main__":
    main()
