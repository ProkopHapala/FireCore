#!/usr/bin/env python3
"""MACE z-scan for CHONH2/Ag(111): E_int and F_z vs height.

Generates interaction energy and force data by running the MAD-SURF MACE
model — no DFT energies extracted. Only atom positions (geometry) are read
from the extxyz training file; the DFT REF_energy field is never used.

Method:
  1. Read ONE MD_interface_Ag frame → fixed Ag slab geometry
  2. Read ONE NMS frame (≤12 atoms) → equilibrium molecule geometry
  3. Orient molecule flat via PCA rotation
  4. Fix Ag slab; scan molecule z from z_min to z_max above surface
  5. E_int(z) = E_MACE(slab+mol@z) − E_MACE(slab) − E_MACE(mol_vacuum)

Saves scan cache: dft_plq_analysis/_mace_zscan_cache.npy
Produces 4 plots in dft_plq_analysis/:
  plot1_eint_vs_z.png
  plot2_fz_vs_z.png
  plot3_eint_fz_combined.png
  plot4_sanity_negdEdz_vs_Fz.png

Usage:
  python3 tests/tGridFFJax/analyze_madsurf_dft_zscan.py
  python3 tests/tGridFFJax/analyze_madsurf_dft_zscan.py --device cuda
  python3 tests/tGridFFJax/analyze_madsurf_dft_zscan.py --z-min 1.5 --z-max 8.0 --z-step 0.1
  python3 tests/tGridFFJax/analyze_madsurf_dft_zscan.py --force-recompute
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from collections import Counter

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_zscan")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ── canonical paths ──────────────────────────────────────────────────────────
HERE       = os.path.dirname(os.path.abspath(__file__))
FIRECORE   = os.path.abspath(os.path.join(HERE, "..", ".."))
EXTXYZ     = os.path.join(HERE, "mad-surf_data",
                          "full_train_test_std_config_types.extxyz")
MACE_MODEL = os.path.join(HERE, "mad-surf_data", "models",
                          "full_dataset_config_weights", "MACE_model.model")
MADSURF_MAIN = os.path.join(HERE, "MAD-SURF-main")
OUT_DIR    = os.path.join(HERE, "dft_plq_analysis")

for _p in [FIRECORE, MADSURF_MAIN]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

_RE_CONFIG  = re.compile(r'config_type=(\S+)')
_RE_LATTICE = re.compile(r'Lattice="([^"]+)"')


# ════════════════════════════════════════════════════════════════════════════
# 1. Geometry from extxyz (positions only — DFT energy never read)
# ════════════════════════════════════════════════════════════════════════════

def read_geometry_frames(extxyz_path):
    """Read first MD_interface_Ag frame; split into Ag slab + molecule.

    The molecule is extracted from the MD frame itself (non-Ag atoms),
    guaranteeing consistency with the MD training data.
    REF_energy is never read.

    Returns
    -------
    ag_syms  : list[str]    Ag symbols
    ag_pos   : (N_ag, 3)
    mol_syms : list[str]    molecule symbols
    mol_pos  : (N_mol, 3)
    cell     : (3, 3)
    """
    cell   = None
    n_total = 0

    print(f"Reading geometry from: {extxyz_path}")
    with open(extxyz_path, "r", errors="ignore") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            try:
                n_atoms = int(line)
            except ValueError:
                continue

            info      = fh.readline().rstrip("\n")
            atom_lines = [fh.readline() for _ in range(n_atoms)]
            n_total  += 1

            m = _RE_CONFIG.search(info)
            ctype = m.group(1) if m else "unknown"

            m = _RE_LATTICE.search(info)
            if m:
                vals = list(map(float, m.group(1).split()))
                cell = np.array(vals, dtype=float).reshape(3, 3)

            if ctype != "MD_interface_Ag":
                continue

            syms, pos = [], []
            ok = True
            for al in atom_lines:
                parts = al.split()
                if len(parts) < 4:
                    ok = False; break
                syms.append(parts[0])
                pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
            if not ok:
                continue
            pos = np.array(pos, dtype=float)

            ag_m     = np.array([s == "Ag" for s in syms])
            ag_syms  = [s for s, m_ in zip(syms, ag_m) if m_]
            mol_syms = [s for s, m_ in zip(syms, ag_m) if not m_]
            ag_pos   = pos[ag_m]
            mol_pos  = pos[~ag_m]

            comp = dict(Counter(mol_syms))
            print(f"  Ag slab  : {ag_m.sum()} atoms (MD_interface_Ag frame {n_total})")
            print(f"  Molecule : {len(mol_syms)} atoms {comp}")
            if cell is not None:
                print(f"  Cell     : {cell[0]}, {cell[1]}, {cell[2]}")
            return ag_syms, ag_pos, mol_syms, mol_pos, cell

    raise RuntimeError("No MD_interface_Ag frame found in extxyz")


# ════════════════════════════════════════════════════════════════════════════
# 2. Molecule orientation — PCA flat rotation
# ════════════════════════════════════════════════════════════════════════════

def orient_molecule_flat(mol_pos):
    """Rotate molecule so its molecular plane is in xy (normal axis → z).

    Uses PCA: the eigenvector with smallest eigenvalue is the normal.

    Returns
    -------
    flat_pos : (N, 3)  centered flat positions (z-spread minimised)
    R        : (3, 3)  rotation matrix applied
    """
    com = mol_pos.mean(axis=0)
    centered = mol_pos - com
    cov = centered.T @ centered
    _eigenvalues, eigenvectors = np.linalg.eigh(cov)   # ascending eigenvalue order
    # eigenvectors[:,0] = normal to molecule plane (smallest eigenvalue)
    # eigenvectors[:,2] = dominant spread axis (largest eigenvalue)
    R = np.column_stack([
        eigenvectors[:, 2],   # new x = dominant spread
        eigenvectors[:, 1],   # new y = middle spread
        eigenvectors[:, 0],   # new z = normal (→ becomes "up")
    ])
    flat_pos = centered @ R
    # ensure "up" direction is consistent: most atoms above COM if needed
    if flat_pos[:, 2].mean() < 0:
        flat_pos[:, 2] *= -1
    z_spread = float(np.ptp(flat_pos[:, 2]))
    print(f"  Flat molecule: z-spread = {z_spread:.3f} Å, "
          f"z ∈ [{flat_pos[:,2].min():.3f}, {flat_pos[:,2].max():.3f}] Å")
    return flat_pos, R


def find_z_top_ag(ag_pos, cell):
    """Return z of the Ag top layer (excluding PBC images near cell top)."""
    cell_z = float(cell[2, 2])
    real_m = ag_pos[:, 2] < (cell_z - 1.0)  # exclude PBC images at ~cell_z
    return float(ag_pos[real_m, 2].max())


# ════════════════════════════════════════════════════════════════════════════
# 3. MACE z-scan
# ════════════════════════════════════════════════════════════════════════════

def run_mace_zscan(ag_syms, ag_pos, mol_syms, mol_pos_flat,
                   cell, z_scan, mace_model_path, device, z_top_ag):
    """Scan molecule vertically above fixed Ag slab using MACE.

    Molecule COM is placed above the top Ag atom (atop site) at each z.
    E_int(z) = E_MACE(slab+mol@z) − E_MACE(slab) − E_MACE(mol_vac)

    Returns dict with arrays: z_rel, e_int [eV], fz_mol, fx_mol, fy_mol [eV/Å]
    """
    from mace.calculators import MACECalculator  # type: ignore
    from ase import Atoms  # type: ignore

    print(f"\nLoading MACE: {mace_model_path}")
    calc = MACECalculator(model_paths=mace_model_path,
                          device=device, default_dtype="float64")

    # Fixed lateral position: above the top Ag atom (atop site)
    cell_z   = float(cell[2, 2])
    real_m   = ag_pos[:, 2] < (cell_z - 1.0)
    top_idx  = np.argmax(ag_pos[real_m][:, 2])
    lat_xy   = ag_pos[real_m][top_idx, :2]
    print(f"  Lateral (atop site): x={lat_xy[0]:.3f}, y={lat_xy[1]:.3f} Å")

    # E_slab (computed once, fixed Ag)
    slab_at = Atoms(symbols=ag_syms, positions=ag_pos,
                    cell=cell, pbc=[True, True, True])
    slab_at.calc = calc
    e_slab = float(slab_at.get_potential_energy())
    print(f"  E_slab = {e_slab:.4f} eV  ({len(ag_syms)} Ag atoms, fixed)")

    # E_mol reference: molecule placed at midpoint of vacuum gap, SAME cell and
    # pbc as E_slab/E_comb. This ensures exact cancellation of any cell-dependent
    # MACE energy terms (Ewald, long-range, etc.) that differ between pbc=True
    # and pbc=False calculations.
    #
    # Safe z: equidistant from real Ag top and Ag bottom PBC image.
    #   z_top_real  ≈ z_top_ag (~4.7 Å)
    #   z_bot_image = cell_z - |z_bot_real| ≈ 32.7 - 0.5 ≈ 32.2 Å
    #   midpoint z  ≈ (z_top_ag + z_bot_image) / 2 ≈ 18.5 Å → 14 Å above z_top_ag
    cell_z       = float(cell[2, 2])
    z_bot_image  = cell_z - 0.5                        # PBC image of bottom Ag
    z_mol_ref    = (z_top_ag + z_bot_image) / 2.0      # midpoint of vacuum gap
    z_rel_ref    = z_mol_ref - z_top_ag
    mol_ref_pos  = mol_pos_flat.copy()
    mol_ref_pos[:, 0] += lat_xy[0]
    mol_ref_pos[:, 1] += lat_xy[1]
    mol_ref_pos[:, 2] += z_mol_ref
    mol_vac = Atoms(symbols=mol_syms, positions=mol_ref_pos,
                    cell=cell, pbc=[True, True, True])
    mol_vac.calc = calc
    e_mol = float(mol_vac.get_potential_energy())
    print(f"  E_mol  = {e_mol:.4f} eV  ({len(mol_syms)} atoms, "
          f"z_ref={z_mol_ref:.2f} Å = z_top+{z_rel_ref:.1f} Å, same cell/pbc)")

    n_ag = len(ag_syms)
    z_rel_out, e_int_out = [], []
    fz_out, fx_out, fy_out = [], [], []

    print(f"\nZ-scan: {len(z_scan)} points  "
          f"z ∈ [{z_scan[0]:.2f}, {z_scan[-1]:.2f}] Å above surface")
    t0 = time.perf_counter()

    for i, z_rel in enumerate(z_scan):
        mol_z = mol_pos_flat.copy()
        mol_z[:, 0] += lat_xy[0]
        mol_z[:, 1] += lat_xy[1]
        mol_z[:, 2] += z_top_ag + z_rel

        comb = Atoms(symbols=ag_syms + mol_syms,
                     positions=np.vstack([ag_pos, mol_z]),
                     cell=cell, pbc=[True, True, True])
        comb.calc = calc
        e_comb = float(comb.get_potential_energy())
        f_all  = np.asarray(comb.get_forces(), dtype=float)
        f_mol  = f_all[n_ag:]

        e_int = e_comb - e_slab - e_mol
        e_int_out.append(e_int)
        fz_out.append(float(f_mol[:, 2].mean()))
        fx_out.append(float(f_mol[:, 0].mean()))
        fy_out.append(float(f_mol[:, 1].mean()))
        z_rel_out.append(z_rel)

        if (i + 1) % 10 == 0 or (i + 1) == len(z_scan):
            print(f"  z={z_rel:.2f}Å  E_int={e_int*1000:+8.1f} meV  "
                  f"Fz={fz_out[-1]:+.4f} eV/Å  "
                  f"[{i+1}/{len(z_scan)}, {time.perf_counter()-t0:.1f}s]")

    elapsed = time.perf_counter() - t0
    print(f"Scan done in {elapsed:.1f}s  ({elapsed/len(z_scan)*1000:.0f} ms/point)")

    return {
        "z_rel":      np.array(z_rel_out),
        "e_int":      np.array(e_int_out),
        "fz_mol":     np.array(fz_out),
        "fx_mol":     np.array(fx_out),
        "fy_mol":     np.array(fy_out),
        "e_slab":     e_slab,
        "e_mol":      e_mol,
        "lat_xy":     lat_xy,
        "z_top_ag":   z_top_ag,
        "extxyz":     EXTXYZ,
        "mace_model": mace_model_path,
    }


# ════════════════════════════════════════════════════════════════════════════
# 4. Plots
# ════════════════════════════════════════════════════════════════════════════

def _save(fig, path):
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved : {path}")


def plot1_eint(scan, out_dir):
    z  = scan["z_rel"]
    ei = scan["e_int"] * 1000  # → meV
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z, ei, "b-o", ms=4, lw=1.5, label="MACE E_int")
    ax.axhline(0, color="k", lw=0.5, ls="--")
    i_min = int(np.argmin(ei))
    ax.annotate(f"min = {ei[i_min]:.1f} meV\n@ z = {z[i_min]:.2f} Å",
                xy=(z[i_min], ei[i_min]),
                xytext=(z[i_min]+0.4, ei[i_min]+50),
                arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax.set_xlabel("z above Ag top layer (Å)")
    ax.set_ylabel("E_int (meV)")
    ax.set_title(f"MACE Interaction Energy  E_int = E(slab+mol) − E(slab) − E(mol)\n"
                 f"MACE: {scan['mace_model']}")
    ax.legend(); ax.grid(True, alpha=0.25)
    _save(fig, os.path.join(out_dir, "plot1_eint_vs_z.png"))


def plot2_fz(scan, out_dir):
    z  = scan["z_rel"]
    fz = scan["fz_mol"]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z, fz, "r-o", ms=4, lw=1.5, label="MACE ⟨Fz⟩ (mean over mol atoms)")
    ax.axhline(0, color="k", lw=0.7, ls="--")
    ax.set_xlabel("z above Ag top layer (Å)")
    ax.set_ylabel("Mean Fz on molecule (eV/Å)")
    ax.set_title(f"MACE Force z-component  (positive = pushed away from surface)\n"
                 f"MACE: {scan['mace_model']}")
    ax.legend(); ax.grid(True, alpha=0.25)
    _save(fig, os.path.join(out_dir, "plot2_fz_vs_z.png"))


def plot3_combined(scan, out_dir):
    z  = scan["z_rel"]
    ei = scan["e_int"] * 1000
    fz = scan["fz_mol"]
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    l1, = ax1.plot(z, ei, "b-o", ms=4, lw=1.5, label="E_int (meV)")
    l2, = ax2.plot(z, fz, "r-s", ms=4, lw=1.5, label="⟨Fz⟩ (eV/Å)")
    ax1.axhline(0, color="b", lw=0.4, ls="--")
    ax2.axhline(0, color="r", lw=0.4, ls="--")
    ax1.set_xlabel("z above Ag top layer (Å)")
    ax1.set_ylabel("E_int (meV)", color="b")
    ax2.set_ylabel("Mean Fz (eV/Å)", color="r")
    ax1.tick_params(axis="y", labelcolor="b")
    ax2.tick_params(axis="y", labelcolor="r")
    ax1.set_title("MACE z-scan: Interaction Energy and Force\n"
                  f"Molecule: {Counter(scan.get('mol_syms', []))}  "
                  f"Lateral: atop site")
    lines = [l1, l2]
    ax1.legend(lines, [str(l.get_label()) for l in lines], loc="upper right")
    ax1.grid(True, alpha=0.25)
    _save(fig, os.path.join(out_dir, "plot3_eint_fz_combined.png"))


def plot4_sanity(scan, out_dir):
    """Sanity check: −dE/dz should equal Fz (force = −grad E)."""
    z  = scan["z_rel"]
    ei = scan["e_int"]
    fz = scan["fz_mol"]

    # Numerical derivative of energy
    neg_dEdz = -np.gradient(ei, z)  # eV/Å

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(z, neg_dEdz, "b-",  lw=2,   label="-dE_int/dz  (numerical)")
    ax.plot(z, fz,       "r--", lw=1.5, label="⟨Fz⟩ MACE  (direct from model)")
    ax.axhline(0, color="k", lw=0.5, ls=":")
    ax.set_xlabel("z above Ag top layer (Å)")
    ax.set_ylabel("Force z-component (eV/Å)")
    ax.set_title("Sanity check: −dE/dz vs Fz  (should overlap if MACE is energy-consistent)\n"
                 "Divergence here means MACE forces ≠ −∇E (gradient of MACE energy)")
    ax.legend(); ax.grid(True, alpha=0.25)
    rmse = float(np.sqrt(np.mean((neg_dEdz - fz)**2)))
    ax.text(0.98, 0.05, f"RMSE = {rmse*1000:.1f} meV/Å",
            transform=ax.transAxes, ha="right", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7))
    _save(fig, os.path.join(out_dir, "plot4_sanity_negdEdz_vs_Fz.png"))


# ════════════════════════════════════════════════════════════════════════════
# 5. Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--z-min",   type=float, default=1.5)
    p.add_argument("--z-max",   type=float, default=8.0)
    p.add_argument("--z-step",  type=float, default=0.1)
    p.add_argument("--device",  default="cpu")
    p.add_argument("--force-recompute", action="store_true",
                   help="Ignore cached scan data and recompute")
    args = p.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    cache_path = os.path.join(OUT_DIR, "_mace_zscan_cache.npy")

    print("=" * 64)
    print("MACE z-scan  —  CHONH2/Ag(111)  interaction energy + forces")
    print(f"  EXTXYZ     : {EXTXYZ}")
    print(f"  MACE model : {MACE_MODEL}")
    print(f"  Output dir : {OUT_DIR}")
    print(f"  Cache      : {cache_path}")
    print(f"  Device     : {args.device}")
    print(f"  z range    : {args.z_min:.2f} → {args.z_max:.2f} Å  step {args.z_step:.2f} Å")
    print("=" * 64)

    # Load or compute scan data
    if os.path.exists(cache_path) and not args.force_recompute:
        print(f"\nLoading cached scan from {cache_path}")
        scan = np.load(cache_path, allow_pickle=True).item()
        print(f"  z range: {scan['z_rel'].min():.2f} → {scan['z_rel'].max():.2f} Å"
              f"  ({len(scan['z_rel'])} points)")
    else:
        ag_syms, ag_pos, mol_syms, mol_pos, cell = read_geometry_frames(EXTXYZ)
        mol_pos_flat, _ = orient_molecule_flat(mol_pos)
        z_top_ag = find_z_top_ag(ag_pos, cell)
        print(f"  z_top_ag = {z_top_ag:.4f} Å")

        z_scan = np.arange(args.z_min, args.z_max + args.z_step / 2, args.z_step)
        scan = run_mace_zscan(ag_syms, ag_pos, mol_syms, mol_pos_flat,
                              cell, z_scan, MACE_MODEL, args.device, z_top_ag)
        scan["mol_syms"] = mol_syms
        np.save(cache_path, scan, allow_pickle=True)  # type: ignore[arg-type]
        print(f"\nCached → {cache_path}")

    # Print summary statistics
    ei = scan["e_int"]
    fz = scan["fz_mol"]
    z  = scan["z_rel"]
    i_min = int(np.argmin(ei))
    print(f"\nScan summary:")
    print(f"  E_int range : {ei.min()*1000:.1f} to {ei.max()*1000:.1f} meV")
    print(f"  E_int min   : {ei[i_min]*1000:.1f} meV  @ z = {z[i_min]:.2f} Å")
    print(f"  Fz range    : {fz.min():.4f} to {fz.max():.4f} eV/Å")
    iz0 = int(np.argmin(np.abs(fz)))  # Fz ≈ 0 → equilibrium
    print(f"  Fz ≈ 0 at z ≈ {z[iz0]:.2f} Å  (equilibrium height)")

    # Generate plots
    print("\nGenerating plots...")
    plot1_eint(scan, OUT_DIR)
    plot2_fz(scan, OUT_DIR)
    plot3_combined(scan, OUT_DIR)
    plot4_sanity(scan, OUT_DIR)
    print(f"\nAll plots in {OUT_DIR}/")


if __name__ == "__main__":
    main()
