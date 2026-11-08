#!/usr/bin/env python3
"""plot_hydride_motif_spectra.py — {100} vs {111}, XH vs XH2, 5-ring vs 7-ring.

Thin orchestration. PDOS partition + stacked plot live in pyBall.FFfit_utils / FFfit_plots.
PySCF L1 stacked figures (`--pyscf`) are the Si-NC plot template: chemistry colors, total DOS, rug, xy inset.
Report: doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md

Usage (from repo root):
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --skip-hessian
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --pyscf /home/prokop/SIMULATIONS/SiNCs/pyscf_vib_results
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --pyscf ... --dftb-l1 /path/to/DFTB/L1
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet ridge_110  # ⟨110⟩ extrema = edge (default)
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet miller_111
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet xh_align
  python3 tests/tSiNCs/plot_hydride_motif_spectra.py --l2-only --facet wulff   # COM support
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from pyBall.FFfit_utils import (
    gaussian_spectrum, stretch_mode_weights, weighted_mean_cm1, xh_bonds_from_topology,
    NBHD_TAG_RGB, NBHD_PLOT_STYLE, nbhd_legend_label, primary_face_families_from_name,
    cartesian_modes_from_hessian, modes_as_3N, xh_bonds_from_coords, neighborhood_xh_groups,
    assert_harmonic_spectrum_at_minimum, write_hydride_halogen_xyz, as_signed_wavenumbers_cm1,
    wulff_families_from_name, atom_stretch_in_window, atom_group_mode_weights, atom_tags_from_owner,
    resolve_pyscf_vib_case_dir,
    apply_ring_tags, load_dftb_vibrations_tag, masses_amu_from_Z,
)
from pyBall.FFfit_plots import plot_stacked_method_pdos, plot_ownmin_method_dos
from pyBall.nanocrystal_pipeline import solve_normal_modes_from_hessian_npz, solve_spectrum

ATLAS = REPO / "tests/tSiNCs/OUT_chem_atlas/atlas"
OUT = REPO / "tests/tSiNCs/OUT_motif_spectra"
PLOT = OUT / "plots"
TEST = REPO / "tests/tSiNCs"

# Experimental windows (Kusová). MMFF default peaks are offset; still shade them as the target.
KUSOVA = {"SiH": (2040.0, 2100.0), "SiH2": (2100.0, 2140.0), "SiH3": (2140.0, 2170.0), "CH": (2815.0, 2900.0), "CH2": (2900.0, 2975.0)}
EXPECTED_GROUPS = {  # n(XH), n(XH2), n(XH3) groups — same graph for Si and C
    "L0_adamantane": (4, 6, 0), "L0_Si10H16": (4, 6, 0),
    "cube": (12, 12, 0), "cube_5ring": (12, 11, 0), "cube_7ring": (12, 13, 0),
    "octahedron": (24, 6, 0), "octahedron_5ring": (24, 5, 0), "octahedron_7ring": (24, 7, 0),
}
CRYSTALS = [
    ("L0_adamantane", ATLAS / "L0_ref/C_adamantane", "C", "cage", "Td cage, mixed CH/CH2"),
    ("L0_Si10H16", ATLAS / "L0_ref/Si_Si10H16", "Si", "cage", "Td cage, mixed SiH/SiH2"),
    ("cube_C", ATLAS / "L1_dft/cube/C", "C", "cube", "{100} CH2-rich"),
    ("cube_Si", ATLAS / "L1_dft/cube/Si", "Si", "cube", "{100} SiH2-rich"),
    ("octa_C", ATLAS / "L1_dft/octahedron/C", "C", "octahedron", "{111} CH-rich"),
    ("octa_Si", ATLAS / "L1_dft/octahedron/Si", "Si", "octahedron", "{111} SiH-rich"),
    ("cube5_C", ATLAS / "L1_dft/cube_5ring/C", "C", "cube_5ring", "collapse one CH2"),
    ("cube5_Si", ATLAS / "L1_dft/cube_5ring/Si", "Si", "cube_5ring", "collapse one SiH2"),
    ("cube7_C", ATLAS / "L1_dft/cube_7ring/C", "C", "cube_7ring", "insert one CH2"),
    ("cube7_Si", ATLAS / "L1_dft/cube_7ring/Si", "Si", "cube_7ring", "insert one SiH2"),
    ("octa5_C", ATLAS / "L1_dft/octahedron_5ring/C", "C", "octahedron_5ring", "collapse one CH2"),
    ("octa5_Si", ATLAS / "L1_dft/octahedron_5ring/Si", "Si", "octahedron_5ring", "collapse one SiH2"),
    ("octa7_C", ATLAS / "L1_dft/octahedron_7ring/C", "C", "octahedron_7ring", "insert one CH2"),
    ("octa7_Si", ATLAS / "L1_dft/octahedron_7ring/Si", "Si", "octahedron_7ring", "insert one SiH2"),
]
SIGMA = 8.0
GRID = np.linspace(0.0, 3600.0, 1801)
STRETCH_WIN = {"C": (1800.0, 2400.0), "Si": (1800.0, 2400.0)}  # MMFF default hydride stretch (not Kusová scale)


def _assert_groups(motif, groups):
    exp = EXPECTED_GROUPS[motif]
    got = (int(groups["XH"].shape[0]), int(groups["XH2"].shape[0] // 2), int(groups["XH3"].shape[0] // 3) if groups["XH3"].shape[0] else 0)
    if groups["XH2"].shape[0] % 2:
        raise ValueError(f"{motif}: XH2 bond count {groups['XH2'].shape[0]} not even")
    if got != exp:
        raise ValueError(f"{motif}: hydride groups {got} != expected {exp}")


def ensure_hessian(crystal_dir: Path, skip: bool) -> Path:
    mol2 = crystal_dir / "relaxed.mol2"
    hess = crystal_dir / "04_hessian.npz"
    spec = crystal_dir / "05_spectrum.npz"
    if not mol2.is_file():
        raise FileNotFoundError(mol2)
    def _bonds_ok(path):
        d = np.load(path)
        if "bonds_ij" not in d.files or "Z" not in d.files:
            return False
        Z = np.asarray(d["Z"]).reshape(-1)
        return not any(int(Z[i]) == 1 and int(Z[j]) == 1 for i, j in np.asarray(d["bonds_ij"]))
    if hess.is_file() and _bonds_ok(hess):
        if not spec.is_file():
            solve_spectrum(hess, spec, out_plot=crystal_dir / "spectrum.png")
        return hess
    print(f"  Hessian FD  {mol2.relative_to(REPO)}  (fresh process — MMFF is not re-entrant)")
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO) + os.pathsep + env.get("PYTHONPATH", "")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pyBall.nanocrystal_pipeline", "hessian", "--relaxed-xyz", str(mol2), "--out-npz", str(hess)], cwd=str(REPO), env=env)
    subprocess.check_call([sys.executable, "-m", "pyBall.nanocrystal_pipeline", "spectrum", "--hessian-npz", str(hess), "--out-npz", str(spec), "--out-plot", str(crystal_dir / "spectrum.png")], cwd=str(REPO), env=env)
    if not _bonds_ok(hess):
        raise ValueError(f"{hess}: H–H bonds after mol2 remap")
    return hess


def analyze_crystal(cid, crystal_dir, elem, motif, note, skip_hess):
    hess_path = ensure_hessian(crystal_dir, skip_hess)
    nm = solve_normal_modes_from_hessian_npz(hess_path)
    pos, Z = nm["pos"], nm["Z"]
    if Z is None:
        raise ValueError(f"{cid}: Hessian missing Z")
    d = np.load(hess_path)
    if "bonds_ij" not in d.files:
        raise KeyError(f"{cid}: {hess_path} missing bonds_ij (MMFF order)")
    bonds = np.asarray(d["bonds_ij"], dtype=np.int32)
    groups, _nH = xh_bonds_from_topology(Z, bonds)
    _assert_groups(motif if motif != "cage" else ("L0_adamantane" if elem == "C" else "L0_Si10H16"), groups)
    W = stretch_mode_weights(pos, nm["modes_cart"], groups)
    mask = nm["vib_mask"]
    om = nm["omegas_cm"][mask]
    om = om.copy()
    w_xh, w_xh2 = W["XH"][mask], W["XH2"][mask]
    w_tot = w_xh + w_xh2 + W["XH3"][mask]
    phys = (om > 10.0) & (om < 4000.0)
    om, w_xh, w_xh2, w_tot = om[phys], w_xh[phys], w_xh2[phys], w_tot[phys]
    lo, hi = STRETCH_WIN[elem]
    mxh, sxh = weighted_mean_cm1(om, w_xh, lo, hi)
    mxh2, sxh2 = weighted_mean_cm1(om, w_xh2, lo, hi)
    if elem == "Si":
        mxh_k, _ = weighted_mean_cm1(om, w_xh, *KUSOVA["SiH"])
        mxh2_k, _ = weighted_mean_cm1(om, w_xh2, *KUSOVA["SiH2"])
    else:
        mxh_k = mxh2_k = float("nan")
    om_phys = om[(om > 10.0) & (om < 4000.0)]
    rec = {
        "id": cid, "elem": elem, "motif": motif, "note": note, "dir": str(crystal_dir),
        "natoms": int(nm["natoms"]), "n_XH": int(groups["XH"].shape[0]), "n_XH2": int(groups["XH2"].shape[0] // 2),
        "mean_XH_cm1": mxh, "mean_XH2_cm1": mxh2, "wsum_XH": sxh, "wsum_XH2": sxh2,
        "mean_XH_kusova": mxh_k, "mean_XH2_kusova": mxh2_k,
        "omega_max_cm1": float(om_phys.max()) if om_phys.size else float("nan"), "hessian": str(hess_path),
    }
    cache = OUT / "modes" / f"{cid}.npz"
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, omega_cm=om, w_XH=w_xh, w_XH2=w_xh2, w_stretch=w_tot)
    rec["cache"] = str(cache)
    print(f"  {cid:16s} N={rec['natoms']:3d}  XH={rec['n_XH']:2d} XH2={rec['n_XH2']:2d}  "
          f"<XH>={mxh:7.1f}  <XH2>={mxh2:7.1f}  ωmax={rec['omega_max_cm1']:.0f}")
    return rec


def load_cache(cid):
    d = np.load(OUT / "modes" / f"{cid}.npz")
    return {k: np.asarray(d[k]) for k in d.files}


def dos(cid, key, sigma=SIGMA):
    d = load_cache(cid)
    w = None if key == "count" else d[key]
    return gaussian_spectrum(d["omega_cm"], GRID, sigma, weights=w)


def shade_kusova(ax, elem):
    if elem == "Si":
        ax.axvspan(*KUSOVA["SiH"], color="#9ecae1", alpha=0.25, label="Kusová SiH 2040–2100")
        ax.axvspan(*KUSOVA["SiH2"], color="#fc9272", alpha=0.25, label="Kusová SiH₂ 2100–2140")
        return
    ax.axvspan(*KUSOVA["CH"], color="#9ecae1", alpha=0.25, label="Kusová CH 2815–2900")
    ax.axvspan(*KUSOVA["CH2"], color="#fc9272", alpha=0.25, label="Kusová CH₂ 2900–2975")


def savefig(fig, name):
    PLOT.mkdir(parents=True, exist_ok=True)
    path = PLOT / name
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path.relative_to(REPO)}")
    return path


def plot_overlay(cids, labels, title, fname, elem, keys=("w_stretch",), xlim=None, ylabel="stretch PDOS (a.u.)"):
    fig, ax = plt.subplots(figsize=(9.2, 3.6))
    colors = plt.cm.tab10(np.linspace(0, 0.7, max(len(cids), 1)))
    ls = ["-", "--", "-.", ":"]
    for i, (cid, lab) in enumerate(zip(cids, labels)):
        for k, key in enumerate(keys):
            y = dos(cid, key)
            ax.plot(GRID, y, color=colors[i], ls=ls[k % len(ls)], lw=1.2, label=f"{lab}" + ("" if len(keys) == 1 else f" {key}"))
    shade_kusova(ax, elem)
    lo, hi = xlim if xlim else STRETCH_WIN[elem]
    ax.set_xlim(lo, hi)
    ax.set_xlabel("ω (cm⁻¹)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    return savefig(fig, fname)


def plot_xh_split(cid, title, fname, elem):
    fig, ax = plt.subplots(figsize=(9.2, 3.6))
    ax.plot(GRID, dos(cid, "w_XH"), color="#2171b5", lw=1.4, label="XH stretch")
    ax.plot(GRID, dos(cid, "w_XH2"), color="#cb181d", lw=1.4, label="XH₂ stretch")
    shade_kusova(ax, elem)
    ax.set_xlim(*STRETCH_WIN[elem])
    ax.set_xlabel("ω (cm⁻¹)")
    ax.set_ylabel("Σ q²  (stretch projection)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return savefig(fig, fname)


def plot_delta(parent, child, title, fname, elem):
    fig, ax = plt.subplots(figsize=(9.2, 3.8))
    yp, yc = dos(parent, "w_stretch"), dos(child, "w_stretch")
    ax.plot(GRID, yp, color="#636363", lw=1.0, label="parent")
    ax.plot(GRID, yc, color="#e6550d", lw=1.2, label="defect")
    ax.plot(GRID, yc - yp, color="#31a354", lw=1.0, ls="--", label="ΔS = defect − parent")
    shade_kusova(ax, elem)
    ax.axhline(0.0, color="k", lw=0.4)
    ax.set_xlim(*STRETCH_WIN[elem])
    ax.set_xlabel("ω (cm⁻¹)")
    ax.set_ylabel("stretch PDOS / ΔS")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return savefig(fig, fname)


def plot_hydride_bars(recs, fname):
    rows = [r for r in recs if r["motif"] not in ("cage",)]
    labels = [r["id"] for r in rows]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    ax.bar(x - 0.18, [r["n_XH"] for r in rows], 0.36, label="XH ({111}-like)", color="#2171b5")
    ax.bar(x + 0.18, [r["n_XH2"] for r in rows], 0.36, label="XH₂ ({100}-like)", color="#cb181d")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("group count")
    ax.set_title("L1 Wulff hydrides (topology bonds, not distance)")
    ax.legend()
    fig.tight_layout()
    return savefig(fig, fname)


def plot_older_faces():
    """Reuse OUT_dftb_vs_mmff frequency lists: C octa {111} vs trunc {100+111} vs rhombic {110}."""
    root = TEST / "OUT_dftb_vs_mmff"
    series = [
        ("octahedron_C", "{111} octa"),
        ("truncated_octahedron_C", "{100}+{111} trunc"),
        ("rhombic_dodecahedron_C", "{110} rhombic"),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=False)
    for ax, kind, xlim in ((axes[0], "mmff_default_frequencies_cm1.npy", (1800, 2300)), (axes[1], "mmff_scaled_frequencies_cm1.npy", (2700, 3300))):
        for name, lab in series:
            om = np.load(root / name / kind)
            assert_harmonic_spectrum_at_minimum(om, ctx=f"plot_older_faces {name} {kind}: ")
            om = om[np.isfinite(om) & (om > 10)]
            y = gaussian_spectrum(om, GRID, 10.0)
            ax.plot(GRID, y, lw=1.2, label=lab)
        ax.set_xlim(*xlim)
        ax.set_ylabel("mode DOS")
        ax.legend(fontsize=8)
        ax.set_title(kind.replace(".npy", "") + "  (L2-size C Wulff, total DOS, not typed XH/XH2)")
    axes[1].set_xlabel("ω (cm⁻¹)")
    axes[0].axvspan(*KUSOVA["SiH"], color="#9ecae1", alpha=0.2)
    axes[0].axvspan(*KUSOVA["SiH2"], color="#fc9272", alpha=0.2)
    axes[1].axvspan(*KUSOVA["CH"], color="#9ecae1", alpha=0.2)
    axes[1].axvspan(*KUSOVA["CH2"], color="#fc9272", alpha=0.2)
    fig.tight_layout()
    return savefig(fig, "older_L2_C_faces_total_DOS.png")


def plot_older_ensemble_hydrides():
    """Sphere ensemble mixes hydrides — cannot isolate 100/111."""
    from pyBall.FFfit_utils import hydride_counts_from_coords
    root = TEST / "OUT_nc_ensemble_v2/data/crystals"
    dirs = sorted(p for p in root.iterdir() if (p / "04_hessian.npz").is_file())
    if not dirs:
        return None
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    xs, ch, ch2, nats = [], [], [], []
    for i, d in enumerate(dirs):
        z = np.load(d / "04_hessian.npz")
        pos, Z = np.asarray(z["pos"]), np.asarray(z["Z"]).reshape(-1)
        c = hydride_counts_from_coords(pos, Z)
        xs.append(i)
        ch.append(c["CH"] + c["SiH"])
        ch2.append(c["CH2"] + c["SiH2"])
        nats.append(int(Z.size))
        ax.text(i, max(ch[-1], ch2[-1]) + 0.5, f"N={nats[-1]}", ha="center", fontsize=7)
    ax.bar(np.array(xs) - 0.18, ch, 0.36, label="XH", color="#2171b5")
    ax.bar(np.array(xs) + 0.18, ch2, 0.36, label="XH₂", color="#cb181d")
    ax.set_xlabel("ensemble sphere index")
    ax.set_ylabel("group count (distance)")
    ax.set_title("OUT_nc_ensemble_v2: random diamond spheres — mixed XH/XH2, no facet control")
    ax.legend()
    fig.tight_layout()
    return savefig(fig, "older_ensemble_spheres_hydrides.png")


TAG_COLOR = NBHD_TAG_RGB
NBHD_KEYS = ["XH@111", "XH2@100", "XH@100", "XH2@111", "XH@110", "XH2@110", "XH@edge", "XH2@edge", "ring5", "ring7"]
NBHD_LS = {k: v["ls"] for k, v in NBHD_PLOT_STYLE.items()}
NBHD_PLOT_C = {k: v["color"] for k, v in NBHD_PLOT_STYLE.items()}
L2_C_CASES = (
    ("octahedron_C", "{111} octa, CH faces + ridge edges"),
    ("truncated_octahedron_C", "{100} CH₂ + {111} CH faces, ridge edges"),
    ("rhombic_dodecahedron_C", "{110} faces (rhombic dodecahedron; no {111} facets)"),
)
DFTB_L2 = Path("/home/prokop/SIMULATIONS/SiNCs/DFTB")
MMFF_L2 = TEST / "OUT_dftb_vs_mmff"
KUSOVA_SHADE_C = ((2815.0, 2900.0, "#9ecae1", "Kusová CH"), (2900.0, 2975.0, "#fc9272", "Kusová CH₂"))


def _mmff_crystal_dir(name):
    """MMFF L2 Hessians at that MMFF's own minimum. Do not fall back to WRONG_at_DFTB_geometry (invalid spectra)."""
    p = MMFF_L2 / name
    proto = p / "mmff_ownmin_protocol.json"
    if (p / "mmff_default_hessian.npy").is_file() and proto.is_file():
        return p
    raise FileNotFoundError(
        f"{p}/mmff_ownmin_protocol.json missing. MMFF-at-DFTB-geometry files in WRONG_at_DFTB_geometry/ are not spectra. "
        "See doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md"
    )


def _load_xyz_Zpos(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    smap = {"H": 1, "C": 6, "Si": 14}
    pos = np.zeros((n, 3)); Z = np.zeros(n, dtype=np.int32)
    for i in range(n):
        w = lines[2 + i].split()
        Z[i] = smap[w[0]]
        pos[i] = (float(w[1]), float(w[2]), float(w[3]))
    return pos, Z


def _atom_chem_colors(Z, nH, pos, owner_tag, ring5, ring7):
    n = len(Z)
    col = np.tile(NBHD_TAG_RGB["bulk"], (n, 1))
    sizes = np.full(n, 8.0)
    com = pos[Z > 1].mean(axis=0)
    for i in range(n):
        if int(Z[i]) == 1:
            continue
        sizes[i] = 38.0
        tag = owner_tag.get(i)
        if tag in TAG_COLOR:
            col[i] = TAG_COLOR[tag]
        if tag in ('ring5', 'ring7') or i in ring5 or i in ring7:
            sizes[i] = 70.0
    for i in range(n):
        if int(Z[i]) != 1:
            continue
        sizes[i] = 20.0
        tag = owner_tag.get(i)
        if tag in TAG_COLOR:
            col[i] = TAG_COLOR[tag]
    return col, sizes


def _owner_from_nbhd(nbhd):
    """Atom → neighborhood tag from bond groups (skips ring5/ring7 extras)."""
    owner = {}
    for tag, bij in nbhd.items():
        if tag in ("ring5", "ring7"):
            continue
        for i, j in np.asarray(bij).reshape(-1, 2):
            owner[int(i)] = tag
            owner[int(j)] = tag
    return owner


def _owner_tags(pos, Z, nH, xh_groups, families=('100', '111', '110')):
    """Wulff COM support tags (legacy). Prefer _owner_from_nbhd after neighborhood_xh_groups."""
    com = pos[Z > 1].mean(axis=0)
    owner = {}
    from pyBall.FFfit_utils import facet_kind_from_vec
    for cls, bij in xh_groups.items():
        for i, j in np.asarray(bij).reshape(-1, 2):
            fac = facet_kind_from_vec(pos[int(i)] - com, families=families)
            tag = f"{cls}@{fac}"
            owner[int(i)] = tag
            owner[int(j)] = tag
    return owner


def plot_nbhd_pdos(om, W, title, fname, xlim, shade_elem, face_families=('100', '111')):
    fig, ax = plt.subplots(figsize=(9.4, 3.8))
    any_line = False
    for k in NBHD_KEYS:
        if k not in W or float(np.sum(W[k])) < 1e-18:
            continue
        ax.plot(GRID, gaussian_spectrum(om, GRID, SIGMA, W[k]), color=NBHD_PLOT_C[k], ls=NBHD_LS[k], lw=1.35, label=nbhd_legend_label(k, shade_elem, face_families))
        any_line = True
    if not any_line:
        raise ValueError(f"plot_nbhd_pdos: no neighborhood weight in {title}")
    shade_kusova(ax, shade_elem)
    ax.set_xlim(*xlim)
    ax.set_xlabel("ω (cm⁻¹)")
    ax.set_ylabel("stretch PDOS")
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    return savefig(fig, fname)


def plot_participation(pos, Z, bonds, part, fname, title, highlight=None):
    from pyBall.plotUtils import plot_nc_views
    p = np.asarray(part, dtype=np.float64)
    vmax = float(np.percentile(p[p > 0], 98)) if np.any(p > 0) else 1.0
    if vmax <= 0:
        vmax = 1.0
    t = np.clip(p / vmax, 0.0, 1.0)
    cmap = plt.cm.YlOrRd
    colors = cmap(t)[:, :3]
    colors[t < 1e-6] = np.array([0.85, 0.85, 0.85])
    sizes = 8.0 + 70.0 * t
    sizes[Z > 1] = np.maximum(sizes[Z > 1], 22.0)
    PLOT.mkdir(parents=True, exist_ok=True)
    path = PLOT / fname
    plot_nc_views(path, pos, colors, sizes, bonds_ij=bonds, highlight=highlight, title=title)
    print(f"  wrote {path.relative_to(REPO)}")
    return path


def analyze_spatial_l1(cid, crystal_dir, elem, xlim, motif, facet_mode="ridge_110", ridge_below_A=0.5):
    from pyBall.plotUtils import plot_nc_views
    PLOT.mkdir(parents=True, exist_ok=True)
    hess = crystal_dir / "04_hessian.npz"
    d = np.load(hess)
    pos, Z, bonds = np.asarray(d["pos"]), np.asarray(d["Z"]), np.asarray(d["bonds_ij"], dtype=np.int32)
    xh, nH = xh_bonds_from_topology(Z, bonds)
    ring_lengths = (5,) if "5ring" in motif else ((7,) if "7ring" in motif else ())
    families = wulff_families_from_name(motif)
    faces = ('111',) if facet_mode == 'miller_111' else primary_face_families_from_name(motif)
    nbhd, rings = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=bonds, ring_lengths=ring_lengths, families=families, facet_mode=facet_mode, face_families=faces, ridge_below_A=ridge_below_A)
    if "5ring" in motif and not rings["ring5"]:
        raise ValueError(f"{cid}: 5ring motif but heavy_cycles found no 5-rings")
    if "7ring" in motif and not rings["ring7"]:
        raise ValueError(f"{cid}: 7ring motif but heavy_cycles found no 7-rings")
    owner_h = _owner_from_nbhd(nbhd)
    tags_h = atom_tags_from_owner(len(Z), owner_h)
    nH_bulk = int(np.sum((Z == 1) & (tags_h == 'bulk')))
    if nH_bulk:
        raise ValueError(f"{cid}: {nH_bulk} hydrogens with no hydride tag")
    owner = apply_ring_tags(owner_h, nbhd, rings)
    tags = atom_tags_from_owner(len(Z), owner)
    r5, r7 = rings["ring5"], rings["ring7"]
    col, sizes = _atom_chem_colors(Z, nH, pos, owner, set(r5), set(r7))
    xyzp = OUT / "xyz" / f"l1_halogen_{cid}.xyz"
    write_hydride_halogen_xyz(xyzp, pos, Z, owner_h, faces, comment=f"{cid} facet={facet_mode} ridge_below_A={ridge_below_A}")
    print(f"  wrote {xyzp.relative_to(REPO)}")
    plot_nc_views(PLOT / f"map_chem_{cid}.png", pos, col, sizes, bonds_ij=bonds, highlight=r5 + r7, title=f"{cid}  facet={facet_mode}  blue={{111}}  red={{100}}  orange=edge")
    bits = ",".join(f"{k}={len(v)}" for k, v in nbhd.items() if len(v))
    print(f"  wrote plots/map_chem_{cid}.png  n5={len(r5)} n7={len(r7)}  {bits}")
    out = {"id": cid, "n_ring5": len(r5), "n_ring7": len(r7), "nbhd_n": {k: int(len(nbhd[k])) for k in nbhd}, "spectrum_ok": False}
    try:
        nm = solve_normal_modes_from_hessian_npz(hess)
    except RuntimeError as e:
        print(f"  SKIP PDOS {cid}: Hessian not at a minimum — maps/xyz only\n    {e}")
        return out
    masses = np.diag(np.asarray(d["M"], dtype=np.float64)).reshape(len(Z), 3)[:, 0]
    mask = (nm["omegas_cm"] > 10.0) & (nm["omegas_cm"] < 4000.0)
    om = nm["omegas_cm"][mask]
    U = nm["modes_cart"][:, mask]
    Wst = stretch_mode_weights(pos, U, nbhd)
    plot_nbhd_pdos(om, Wst, f"{cid} neighborhood stretch PDOS  facet={facet_mode}", f"nbhd_{cid}.png", xlim, elem, faces)
    Wdos = atom_group_mode_weights(U, masses, tags)
    plot_stacked_method_pdos(
        GRID, [dict(method="MMFF default (own min)", omega_cm=om, w_groups=Wdos, note="")],
        PLOT / f"l1_stack_{cid}.png", elem=elem, face_families=faces, xlim=(0.0, 3600.0),
        title=f"{cid}  facet={facet_mode}", pos=pos, atom_colors=col, atom_sizes=sizes,
    )
    print(f"  wrote plots/l1_stack_{cid}.png")
    all_xh = np.vstack([xh[k] for k in xh if len(xh[k])])
    wins = [("SiHwin", 2040, 2100), ("SiH2win", 2100, 2140), ("hi2245", 2220, 2270)] if elem == "Si" else [("defCH", 1800, 1950), ("hiCH", 2200, 2350)]
    for wname, lo, hi in wins:
        part = atom_stretch_in_window(pos, U, om, all_xh, lo, hi)
        plot_participation(pos, Z, bonds, part, f"loc_{cid}_{wname}.png", f"{cid}: stretch localization {lo:.0f}–{hi:.0f} cm⁻¹ (yellow→red)", highlight=r5 + r7)
    out["spectrum_ok"] = True
    return out


PYSCF_CASES = [  # job folder name, atlas mol2 dir, elem, motif — Z order matches mol2 (checked)
    ("cube_C", ATLAS / "L1_dft/cube/C", "C", "cube"),
    ("cube_5ring_C", ATLAS / "L1_dft/cube_5ring/C", "C", "cube_5ring"),
    ("cube_7ring_C", ATLAS / "L1_dft/cube_7ring/C", "C", "cube_7ring"),
    ("octahedron_C", ATLAS / "L1_dft/octahedron/C", "C", "octahedron"),
    ("octahedron_5ring_C", ATLAS / "L1_dft/octahedron_5ring/C", "C", "octahedron_5ring"),
    ("octahedron_7ring_C", ATLAS / "L1_dft/octahedron_7ring/C", "C", "octahedron_7ring"),
    ("cube_Si", ATLAS / "L1_dft/cube/Si", "Si", "cube"),
    ("cube_5ring_Si", ATLAS / "L1_dft/cube_5ring/Si", "Si", "cube_5ring"),
    ("cube_7ring_Si", ATLAS / "L1_dft/cube_7ring/Si", "Si", "cube_7ring"),
    ("octahedron_Si", ATLAS / "L1_dft/octahedron/Si", "Si", "octahedron"),
    ("octahedron_5ring_Si", ATLAS / "L1_dft/octahedron_5ring/Si", "Si", "octahedron_5ring"),
    ("octahedron_7ring_Si", ATLAS / "L1_dft/octahedron_7ring/Si", "Si", "octahedron_7ring"),
]


def _nbhd_pdos_row(pos, Z, bonds, name, elem, motif, omega_cm, modes, masses, facet_mode="ridge_110", ridge_below_A=0.5):
    """Neighborhood atom-group PDOS at this geometry. Bonds from atlas mol2 (Z order must match xyz)."""
    xh, nH = xh_bonds_from_topology(Z, bonds)
    _assert_groups(motif, xh)
    ring_lengths = (5,) if "5ring" in motif else ((7,) if "7ring" in motif else ())
    families = wulff_families_from_name(motif)
    faces = primary_face_families_from_name(motif)
    nbhd, rings = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=bonds, ring_lengths=ring_lengths, families=families, facet_mode=facet_mode, face_families=faces, ridge_below_A=ridge_below_A)
    if "5ring" in motif and not rings["ring5"]:
        raise ValueError(f"{name}: 5ring motif but heavy_cycles found no 5-rings")
    if "7ring" in motif and not rings["ring7"]:
        raise ValueError(f"{name}: 7ring motif but heavy_cycles found no 7-rings")
    owner_h = _owner_from_nbhd(nbhd)
    tags_h = atom_tags_from_owner(len(Z), owner_h)
    nH_bulk = int(np.sum((Z == 1) & (tags_h == "bulk")))
    if nH_bulk:
        raise ValueError(f"{name}: {nH_bulk} hydrogens with no hydride tag")
    owner = apply_ring_tags(owner_h, nbhd, rings)
    tags = atom_tags_from_owner(len(Z), owner)
    om = as_signed_wavenumbers_cm1(omega_cm)
    U = modes_as_3N(modes, len(Z))
    masses = np.asarray(masses, dtype=np.float64)
    if U.shape[1] != om.shape[0]:
        raise ValueError(f"{name}: modes {U.shape} vs freq {om.shape}")
    if masses.shape != (len(Z),):
        raise ValueError(f"{name}: masses {masses.shape} vs N={len(Z)}")
    n_imag = int(np.sum(om < -10.0))
    phys = (om > 10.0) & (om < 4000.0)
    om_p, U_p = om[phys], U[:, phys]
    W = atom_group_mode_weights(U_p, masses, tags)
    r5, r7 = rings["ring5"], rings["ring7"]
    col, sizes = _atom_chem_colors(Z, nH, pos, owner, set(r5), set(r7))
    note = f"{n_imag} imag — not a minimum" if n_imag else ""
    return dict(method=name, omega_cm=om_p, omega_all=om, w_groups=W, note=note, pos=pos, atom_colors=col, atom_sizes=sizes,
                bonds_ij=bonds, elem=elem, faces=faces, n_imag=n_imag, nbhd=nbhd, Z=Z)


def _atlas_bonds_and_Z(atlas_dir, xyz_path, name):
    """Bonds from atlas mol2; fail loud if this xyz's Z order differs from atlas relaxed.xyz."""
    from pyBall.io.crystal_npz import load_bonds_ij_from_mol2
    mol2 = atlas_dir / "relaxed.mol2"
    atlas_xyz = atlas_dir / "relaxed.xyz"
    if not mol2.is_file():
        raise FileNotFoundError(mol2)
    if not atlas_xyz.is_file():
        raise FileNotFoundError(atlas_xyz)
    pos, Z = _load_xyz_Zpos(xyz_path)
    _, Z_a = _load_xyz_Zpos(atlas_xyz)
    if Z.shape != Z_a.shape or not np.array_equal(Z, Z_a):
        raise ValueError(f"{name}: Z order of {xyz_path} != {atlas_xyz} (mol2 bonds would be wrong)")
    return pos, Z, load_bonds_ij_from_mol2(mol2)


def _pyscf_nbhd_row(job_dir, atlas_dir, name, elem, motif, facet_mode="ridge_110", ridge_below_A=0.5):
    """Neighborhood atom-group PDOS from PySCF npy. Same tagging as analyze_spatial_l1."""
    pos, Z, bonds = _atlas_bonds_and_Z(atlas_dir, job_dir / "relaxed.xyz", name)
    om = as_signed_wavenumbers_cm1(np.load(job_dir / "frequencies_cm1.npy"))
    U = np.load(job_dir / "modes.npy")
    masses = np.asarray(np.load(job_dir / "masses.npy"), dtype=np.float64)
    return _nbhd_pdos_row(pos, Z, bonds, name, elem, motif, om, U, masses, facet_mode=facet_mode, ridge_below_A=ridge_below_A)


DFTB_L1_SK = {  # elem → (folder suffix, plot label, overlay color)
    "C": (("3ob-3-1", "DFTB+ 3ob-3-1", "#2171b5"), ("mio-1-1", "DFTB+ mio-1-1", "#d94801")),
    "Si": (("pbc-0-3", "DFTB+ pbc-0-3", "#2171b5"), ("matsci-0-3", "DFTB+ matsci-0-3", "#d94801")),
}


def _dftb_l1_nbhd_row(job_dir, atlas_dir, name, elem, motif, facet_mode="wulff", ridge_below_A=0.5):
    """Neighborhood PDOS from DFTB+ L1 vibrations.tag (no npy in that bake). Own-minimum xyz."""
    xyz = job_dir / "relaxed.xyz"
    if not xyz.is_file():
        xyz = job_dir / "geo_end.xyz"
    if not xyz.is_file():
        raise FileNotFoundError(f"{job_dir}: neither relaxed.xyz nor geo_end.xyz")
    tag = job_dir / "vibrations.tag"
    if not tag.is_file():
        raise FileNotFoundError(tag)
    pos, Z, bonds = _atlas_bonds_and_Z(atlas_dir, xyz, name)
    om, U = load_dftb_vibrations_tag(tag, len(Z))
    st_path = job_dir / "status.json"
    if st_path.is_file():
        st = json.loads(st_path.read_text())
        fjs = np.asarray(st["frequencies"], dtype=np.float64)
        if fjs.shape != om.shape:
            raise ValueError(f"{name}: status.json nfreq={fjs.shape} vs tag {om.shape}")
        dmax = float(np.max(np.abs(om - fjs)))
        if dmax > 0.05:
            raise ValueError(f"{name}: tag ω vs status.json max|Δ|={dmax:.3f} cm⁻¹ (Hartree→cm⁻¹ mismatch)")
        if int(st.get("n_imag", -1)) != int(np.sum(om < -10.0)):
            raise ValueError(f"{name}: status n_imag={st.get('n_imag')} vs tag {int(np.sum(om < -10.0))}")
    masses = masses_amu_from_Z(Z)
    rec = _nbhd_pdos_row(pos, Z, bonds, name, elem, motif, om, U, masses, facet_mode=facet_mode, ridge_below_A=ridge_below_A)
    rec["omega_all"] = om
    return rec


def analyze_pyscf_jobs(jobs_dir, out_dir, facet_mode="ridge_110"):
    """Replace ad-hoc per-crystal envelopes with plot_stacked_method_pdos (neighborhood PDOS + eigenstate rug)."""
    jobs_dir = Path(jobs_dir); out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows_c, rows_si = [], []
    for name, atlas_dir, elem, motif in PYSCF_CASES:
        rec = _pyscf_nbhd_row(resolve_pyscf_vib_case_dir(jobs_dir, name), atlas_dir, name, elem, motif, facet_mode=facet_mode)
        plot_stacked_method_pdos(
            GRID, [dict(method=f"PySCF PBE  {name}", omega_cm=rec["omega_cm"], w_groups=rec["w_groups"], note=rec["note"])],
            out_dir / f"pdos_{name}.png", elem=elem, face_families=rec["faces"],
            xlim=(0.0, 3300.0) if elem == "C" else (0.0, 2450.0),
            title=f"{name}  PySCF PBE  facet={facet_mode}", pos=rec["pos"], atom_colors=rec["atom_colors"], atom_sizes=rec["atom_sizes"],
        )
        bits = " ".join(f"{k}={len(rec['nbhd'][k])}" for k in rec["nbhd"] if len(rec["nbhd"][k]))
        print(f"  {name:24s} n_imag={rec['n_imag']:2d}  {bits}  -> pdos_{name}.png")
        (rows_c if elem == "C" else rows_si).append(rec)
    plot_stacked_method_pdos(
        GRID, [dict(method=r["method"], omega_cm=r["omega_cm"], w_groups=r["w_groups"], note=r["note"],
                    pos=r["pos"], atom_colors=r["atom_colors"], atom_sizes=r["atom_sizes"], bonds_ij=r["bonds_ij"]) for r in rows_c],
        out_dir / "stacked_C.png", elem="C", face_families=("100", "111"), xlim=(0.0, 3300.0),
        title="PySCF PBE diamond L1 — neighborhood PDOS (same colors as L1/L2 motif plots)",
    )
    plot_stacked_method_pdos(
        GRID, [dict(method=r["method"], omega_cm=r["omega_cm"], w_groups=r["w_groups"], note=r["note"],
                    pos=r["pos"], atom_colors=r["atom_colors"], atom_sizes=r["atom_sizes"], bonds_ij=r["bonds_ij"]) for r in rows_si],
        out_dir / "stacked_Si.png", elem="Si", face_families=("100", "111"), xlim=(0.0, 2450.0),
        title="PySCF PBE Si L1 — neighborhood PDOS (red notes = leftover imaginaries)",
    )
    print(f"  wrote {out_dir / 'stacked_C.png'}")
    print(f"  wrote {out_dir / 'stacked_Si.png'}")
    return rows_c, rows_si


def _nbhd_count_str(nbhd):
    return " ".join(f"{k}={len(nbhd[k])}" for k in nbhd if len(nbhd[k]))


def _stretch_peaks(om, W, elem):
    lo, hi = (2700.0, 3300.0) if elem == "C" else (1800.0, 2300.0)
    out = {}
    for k, w in W.items():
        m, s = weighted_mean_cm1(om, w, lo, hi)
        if s > 1e-12:
            out[k] = {"mean_cm1": m, "wsum": s}
    return out


def analyze_l1_pyscf_vs_dftb(pyscf_dir, dftb_dir, out_dir, facet_mode="wulff"):
    """Same 12 L1 crystals: PySCF PBE vs DFTB+ SK sets. Stacked neighborhood PDOS + untyped DOS overlay."""
    pyscf_dir, dftb_dir, out_dir = Path(pyscf_dir), Path(dftb_dir), Path(out_dir)
    if not pyscf_dir.is_dir():
        raise FileNotFoundError(pyscf_dir)
    if not dftb_dir.is_dir():
        raise FileNotFoundError(dftb_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for name, atlas_dir, elem, motif in PYSCF_CASES:
        pyscf_job = resolve_pyscf_vib_case_dir(pyscf_dir, name)
        rec_p = _pyscf_nbhd_row(pyscf_job, atlas_dir, name, elem, motif, facet_mode=facet_mode)
        rec_p["method"] = "PySCF PBE"
        rec_p["color"] = "#111111"
        dftb_rows = []
        for sk, lab, col in DFTB_L1_SK[elem]:
            rec_d = _dftb_l1_nbhd_row(dftb_dir / f"{name}_{sk}", atlas_dir, f"{name}_{sk}", elem, motif, facet_mode=facet_mode)
            rec_d["method"] = lab
            rec_d["color"] = col
            dftb_rows.append(rec_d)
        rows = [rec_p] + dftb_rows
        xlim = (0.0, 3300.0) if elem == "C" else (0.0, 2450.0)
        stacked_path = out_dir / f"compare_{name}.png"
        plot_stacked_method_pdos(
            GRID,
            [dict(method=r["method"], omega_cm=r["omega_cm"], w_groups=r["w_groups"], note=r["note"],
                  pos=r["pos"], atom_colors=r["atom_colors"], atom_sizes=r["atom_sizes"], bonds_ij=r["bonds_ij"]) for r in rows],
            stacked_path, elem=elem, face_families=rows[0]["faces"], xlim=xlim,
            title=f"{name}  PySCF PBE vs DFTB+  facet={facet_mode}",
        )
        spec_dir = out_dir / name
        overlay_rows = []
        for r in rows:
            overlay_rows.append(dict(
                label=f"{r['method']}" + (f"  ({r['note']})" if r["note"] else " (own min)"),
                omega_cm=r["omega_all"], color=r["color"],
                require_minimum=(r["n_imag"] == 0 and r["omega_all"].size == 3 * len(r["Z"])),
            ))
        ov = plot_ownmin_method_dos(overlay_rows, spec_dir, f"{name}  — each method at its own minimum", grid=GRID, sigma=SIGMA)
        bits = " | ".join(f"{r['method']}: n_imag={r['n_imag']} {_nbhd_count_str(r['nbhd'])}" for r in rows)
        print(f"  {name:24s} {bits}")
        print(f"    stacked {stacked_path.name}  overlay {ov.relative_to(out_dir)}")
        summary.append({
            "name": name, "elem": elem, "motif": motif,
            "methods": [{
                "method": r["method"], "n_imag": r["n_imag"],
                "nbhd_n": {k: int(len(r["nbhd"][k])) for k in r["nbhd"] if len(r["nbhd"][k])},
                "stretch": _stretch_peaks(r["omega_cm"], r["w_groups"], elem),
                "omega_max": float(np.max(r["omega_cm"])) if len(r["omega_cm"]) else float("nan"),
            } for r in rows],
        })
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    imgs = "\n".join(
        f'<h2 id="{s["name"]}">{s["name"]}</h2>\n'
        f'<p>n_imag: ' + ", ".join(f'{m["method"]}={m["n_imag"]}' for m in s["methods"]) + "</p>\n"
        f'<p><img src="compare_{s["name"]}.png" width="1100" alt="stacked PDOS {s["name"]}"></p>\n'
        f'<p><img src="{s["name"]}/spectra_overlay.png" width="900" alt="DOS overlay {s["name"]}">\n'
        f'<img src="{s["name"]}/spectra_stretch.png" width="900" alt="stretch {s["name"]}"></p>'
        for s in summary
    )
    (out_dir / "index.html").write_text(
        "<!DOCTYPE html><html><head><meta charset='utf-8'><title>L1 PySCF vs DFTB+</title></head><body>\n"
        "<h1>L1 PySCF PBE vs DFTB+ (neighborhood PDOS + DOS)</h1>\n"
        "<p>Same 12 chem-atlas crystals. Each method at its own minimum. Stacked colors = neighborhood groups "
        "(Σ w = 1). Red notes = Hessian not a stationary point. "
        "SSOT: <code>doc/Topics/FTIR_Nanocrystals/PySCF_L1_neighborhood_PDOS.md</code>.</p>\n"
        f"{imgs}\n</body></html>\n"
    )
    print(f"  wrote {out_dir / 'summary.json'}  {out_dir / 'index.html'}")
    return summary


def _l2_nbhd(pos, Z, name, facet_mode="xh_align", align_cos=0.9, ridge_below_A=0.5):
    xh, nH = xh_bonds_from_coords(pos, Z)
    families = wulff_families_from_name(name)
    face_fam = primary_face_families_from_name(name)
    nbhd, rings = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=None, families=families, facet_mode=facet_mode, face_families=face_fam, align_cos=align_cos, ridge_below_A=ridge_below_A)
    return xh, nH, nbhd, rings, families


def _l2_cache_load(path):
    d = np.load(path)
    keys = [str(x) for x in d["keys"]]
    Wst = np.asarray(d["W"])
    if Wst.shape[0] != len(keys):
        raise ValueError(f"{path}: W {Wst.shape} vs keys {len(keys)}")
    return np.asarray(d["omega_cm"]), {k: Wst[i] for i, k in enumerate(keys)}


def _l2_cache_save(path, om, W):
    keys = [k for k in NBHD_KEYS if k in W]
    extra = [k for k in W if k not in keys]
    keys = keys + extra
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, omega_cm=np.asarray(om), keys=np.array(keys), W=np.stack([np.asarray(W[k]) for k in keys]))


def _mmff_modes_project(name, pos, Z, atom_tags, kind, skip_cache, masses, facet_mode="xh_align", align_cos=0.9, ridge_below_A=0.5):
    """kind = 'default' | 'scaled'. Hessian at that MMFF's own minimum, mol2/DFTB atom order (eV/Å²)."""
    cache = OUT / "modes" / f"l2_{name}_mmff_{kind}_atompdos_{facet_mode}_c{int(round(100 * float(align_cos)))}_b{int(round(100 * float(ridge_below_A)))}.npz"
    cdir = _mmff_crystal_dir(name)
    hess_p = cdir / f"mmff_{kind}_hessian.npy"
    freq_p = cdir / f"mmff_{kind}_frequencies_cm1.npy"
    if not hess_p.is_file():
        raise FileNotFoundError(hess_p)
    saved_all = np.asarray(np.load(freq_p), dtype=np.float64)
    assert_harmonic_spectrum_at_minimum(saved_all, ctx=f"{freq_p}: ")
    if skip_cache and cache.is_file():
        return _l2_cache_load(cache)
    n3 = 3 * len(Z)
    H = np.asarray(np.load(hess_p), dtype=np.float64)
    if H.shape != (n3, n3):
        raise ValueError(f"{hess_p}: shape {H.shape} want ({n3},{n3})")
    print(f"  eigh MMFF {kind}  {name}  {n3}×{n3}")
    om, U = cartesian_modes_from_hessian(H, masses)
    saved = saved_all[np.isfinite(saved_all) & (saved_all > 10.0)]
    if saved.size == 0:
        raise ValueError(f"{freq_p}: no frequencies > 10")
    dmax = abs(float(om.max()) - float(saved.max()))
    if dmax > 20.0:
        raise ValueError(f"{name} MMFF {kind}: eigh νmax={om.max():.1f} vs saved {saved.max():.1f} (Δ={dmax:.1f})")
    W = atom_group_mode_weights(U, masses, atom_tags)
    _l2_cache_save(cache, om, W)
    return om, W


def _group_peaks(om, W, lo, hi):
    out = {}
    for k, w in W.items():
        if float(np.sum(w)) < 1e-18:
            continue
        m, s = weighted_mean_cm1(om, w, lo, hi)
        out[k] = {"mean_cm1": m, "wsum": s}
    return out


def analyze_l2_dftb_vs_mmff(skip_cache=False, facet_mode="ridge_110", align_cos=0.9, ridge_below_A=0.5):
    """Same 3 L2 C crystals: DFTB+ 3ob-3-1 vs MMFF. Atom-group PDOS partitions the DOS."""
    PLOT.mkdir(parents=True, exist_ok=True)
    recs = []
    dftb_omegas = []
    for name, lab in L2_C_CASES:
        dftb = DFTB_L2 / name
        if not (dftb / "modes.npy").is_file():
            raise FileNotFoundError(dftb / "modes.npy")
        pos, Z = _load_xyz_Zpos(dftb / "geo_end.xyz")
        masses = np.asarray(np.load(dftb / "masses.npy"), dtype=np.float64)
        if masses.shape != (len(Z),):
            raise ValueError(f"{name}: masses {masses.shape} vs N={len(Z)}")
        xh, nH, nbhd, rings, families = _l2_nbhd(pos, Z, name, facet_mode=facet_mode, align_cos=align_cos, ridge_below_A=ridge_below_A)
        faces = ('111',) if facet_mode == 'miller_111' else primary_face_families_from_name(name)
        owner = _owner_from_nbhd(nbhd)
        tags = atom_tags_from_owner(len(Z), owner)
        nH_bulk = int(np.sum((Z == 1) & (tags == 'bulk')))
        if nH_bulk:
            raise ValueError(f"{name}: {nH_bulk} hydrogens with no hydride tag")
        xyzp = OUT / "xyz" / f"l2_halogen_{name}.xyz"
        write_hydride_halogen_xyz(xyzp, pos, Z, owner, faces, comment=f"{name} facet={facet_mode} align_cos={align_cos} ridge_below_A={ridge_below_A}")
        xyzp_mode = OUT / "xyz" / f"l2_halogen_{name}_{facet_mode}.xyz"
        write_hydride_halogen_xyz(xyzp_mode, pos, Z, owner, faces, comment=f"{name} facet={facet_mode} align_cos={align_cos} ridge_below_A={ridge_below_A}")
        print(f"  wrote {xyzp.relative_to(REPO)}  and {xyzp_mode.name}  Jmol F=CH-face Cl=CH-edge Br=CH2-face I=CH2-edge")
        col, sizes = _atom_chem_colors(Z, nH, pos, owner, set(), set())
        legend = "  ".join(nbhd_legend_label(k, "C", faces) for k in (['bulk'] + NBHD_KEYS) if (k == 'bulk' or (k in nbhd and len(nbhd[k]))))
        def _nbhd_bits(n):
            return " ".join(f"{k}={len(n[k])}" for k in NBHD_KEYS if k in n and len(n[k]))
        print(f"  map {name}  facet={facet_mode}  N={len(Z)}  bulk={int(np.sum(tags=='bulk'))}  {_nbhd_bits(nbhd)}")
        for other in ("wulff", "xh_align", "miller_111", "ridge_110"):
            if other == facet_mode:
                continue
            _, _, nbo, _, _ = _l2_nbhd(pos, Z, name, facet_mode=other, align_cos=align_cos, ridge_below_A=ridge_below_A)
            print(f"    counts {other:12s} {_nbhd_bits(nbo)}")
        om_d = np.asarray(np.load(dftb / "frequencies_cm1.npy"), dtype=np.float64)
        assert_harmonic_spectrum_at_minimum(om_d, ctx=f"DFTB+ {name}: ")
        U_d = modes_as_3N(np.load(dftb / "modes.npy"), len(Z))
        if U_d.shape[1] != om_d.shape[0]:
            raise ValueError(f"{name}: DFTB modes {U_d.shape} vs freq {om_d.shape}")
        W_d = atom_group_mode_weights(U_d, masses, tags)
        om_m, W_m = _mmff_modes_project(name, pos, Z, tags, "default", skip_cache, masses, facet_mode=facet_mode, align_cos=align_cos, ridge_below_A=ridge_below_A)
        om_s, W_s = _mmff_modes_project(name, pos, Z, tags, "scaled", skip_cache, masses, facet_mode=facet_mode, align_cos=align_cos, ridge_below_A=ridge_below_A)
        stacked = PLOT / f"l2_compare_{name}.png"
        plot_stacked_method_pdos(
            GRID,
            [
                dict(method="DFTB+ 3ob-3-1", omega_cm=om_d, w_groups=W_d, note=""),
                dict(method="MMFF default", omega_cm=om_m, w_groups=W_m, note=""),
                dict(method="MMFF CH×1.995", omega_cm=om_s, w_groups=W_s, note=""),
            ],
            stacked, elem="C", face_families=faces, xlim=(0.0, 3600.0),
            title=f"{name}", pos=pos, atom_colors=col, atom_sizes=sizes,
        )
        print(f"  wrote {stacked.relative_to(REPO)}")
        om_m_all = np.asarray(np.load(_mmff_crystal_dir(name) / "mmff_default_frequencies_cm1.npy"), dtype=np.float64)
        om_s_all = np.asarray(np.load(_mmff_crystal_dir(name) / "mmff_scaled_frequencies_cm1.npy"), dtype=np.float64)
        ov = plot_ownmin_method_dos(
            [
                dict(label="DFTB+ 3ob-3-1 (own min)", omega_cm=om_d, color="#111111"),
                dict(label="MMFF default (own min)", omega_cm=om_m_all, color="#1f77b4"),
                dict(label="MMFF CH×1.995 ang×0.30 (own min)", omega_cm=om_s_all, color="#ff7f0e"),
            ],
            MMFF_L2 / name,
            f"{name}  — each method at its own minimum",
            grid=GRID, sigma=8.0,
        )
        print(f"  wrote {ov.relative_to(REPO)}")
        rec = {
            "name": name, "lab": lab, "natoms": int(len(Z)), "nH": int((Z == 1).sum()),
            "nbhd_n": {k: int(len(nbhd[k])) for k in nbhd},
            "n_bulk": int(np.sum(tags == 'bulk')),
            "facet_mode": facet_mode,
            "faces": list(faces),
            "legend": legend,
            "peaks": {
                "DFTB+ 3ob-3-1": _group_peaks(om_d, W_d, 2700, 3300),
                "MMFF default": _group_peaks(om_m, W_m, 1800, 2400),
                "MMFF CH×1.995": _group_peaks(om_s, W_s, 2700, 3300),
            },
        }
        recs.append(rec)
        for method, peaks in rec["peaks"].items():
            bits = "  ".join(f"{nbhd_legend_label(k,'C',faces)}={v['mean_cm1']:.0f}" for k, v in peaks.items() if k != 'bulk' and np.isfinite(v["mean_cm1"]))
            print(f"    {method:18s}  {bits}")
        dftb_omegas.append((name, om_d))
    cols = ("#2171b5", "#cb181d", "#238b45")
    plot_ownmin_method_dos(
        [dict(label=n, omega_cm=om, color=cols[i]) for i, (n, om) in enumerate(dftb_omegas)],
        MMFF_L2 / "dftb_only",
        "DFTB+ 3ob-3-1  — three L2 C Wulff cuts (own minima)",
        grid=GRID, sigma=8.0,
    )
    ov_dftb = MMFF_L2 / "dftb_only" / "spectra_overlay.png"
    (MMFF_L2 / "dftb_only" / "overlay_spectra.png").write_bytes(ov_dftb.read_bytes())
    print(f"  wrote {(MMFF_L2 / 'dftb_only' / 'overlay_spectra.png').relative_to(REPO)}")
    plot_older_faces()
    return recs


def analyze_spatial_l2_dftb():
    """Kept as alias: DFTB maps are written by analyze_l2_dftb_vs_mmff."""
    return analyze_l2_dftb_vs_mmff(skip_cache=True)


def _l2_section_html(l2_recs):
    if not l2_recs:
        return ""
    def img(name, cap):
        return f'<figure><img src="plots/{name}" alt=""/><figcaption>{cap}</figcaption></figure>'
    blocks = []
    for rec in l2_recs:
        name, lab = rec["name"], rec["lab"]
        faces = tuple(rec["faces"]) if rec.get("faces") else primary_face_families_from_name(name)
        tags = [k for k in NBHD_KEYS if rec["nbhd_n"].get(k)]
        head = "<tr><th>method</th>" + "".join(f"<th>{nbhd_legend_label(k,'C',faces)}</th>" for k in tags) + "</tr>"
        rows = []
        for method, peaks in rec["peaks"].items():
            cells = []
            for k in tags:
                v = peaks.get(k)
                cells.append("—" if v is None or not np.isfinite(v["mean_cm1"]) else f"{v['mean_cm1']:.0f}")
            rows.append("<tr><td>" + method + "</td>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
        counts = ", ".join(f"{nbhd_legend_label(k,'C',faces)}={rec['nbhd_n'][k]}" for k in tags)
        fm = rec.get("facet_mode", "xh_align")
        blocks.append(f"""
<h3>{name}</h3>
<p>{lab}. N={rec['natoms']}, nH={rec['nH']}, bulk C={rec.get('n_bulk','?')}. facet={fm}. Bond counts: {counts}.</p>
{img(f'l2_compare_{name}.png', 'Top: atom colors = PDOS colors. Then one panel per method: stacked group PDOS; black line = total DOS = sum of the stack.')}
<table>{head}{''.join(rows)}</table>
<p>Table: stretch-weighted ⟨ω⟩ of each group in that method’s hydride window. The plots themselves are the full DOS partition, not this table.</p>
""")
    return f"""
<h2 id="l2-compare">L2 C crystals — DFTB+ 3ob-3-1 vs MMFF (same groups, methods stacked)</h2>
<p>These are the three diamonds that already have DFTB+ Hessians
(<code>octahedron_C</code>, <code>truncated_octahedron_C</code>, <code>rhombic_dodecahedron_C</code>).
Each method’s Hessian is at <b>that method’s own minimum</b> (DFTB+ SCC relax; MMFF FIRE then Hessian). Neighborhood tags are taken from the DFTB geometry (same atom order as atlas mol2). Do not use <code>WRONG_at_DFTB_geometry/</code>.
Face vs edge (switch <code>--facet</code>):
<code>ridge_110</code> (default) = C/Si within 0.5 Å of a ⟨110⟩ extremity are edge (12 signed / 6 unsigned axes); leftover hydrides get the crystal’s primary Miller face.
<code>miller_111</code> = sitting {{111}} octant from heavy−COM, then X–H within align_cos=0.90 of <b>that</b> ⟨111⟩.
<code>xh_align</code> = hydride-neighbor terrace vs ridge.
<code>wulff</code> = COM support.
Colors: blue = CH@{{111}} face, orange/brown = edge. Jmol: <code>OUT_motif_spectra/xyz/l2_halogen_*.xyz</code> — F=CH face, Cl=CH edge, Br=CH₂ face, I=CH₂ edge.</p>
<div class="ok"><b>How to read:</b> each mode’s weight is partitioned over atoms (p<sub>i</sub> = m<sub>i</sub>|u<sub>i</sub>|²). Groups are those atoms. Stacked PDOS sums to the black total-DOS line. Gray = diamond core; colors = hydride neighborhoods (same as the 3D map).</div>
{''.join(blocks)}
"""


def write_html(recs, plots, l2_recs=None):
    recs = list(recs or [])
    l2_recs = list(l2_recs or [])
    if not recs:
        sp = OUT / "stats.json"
        if sp.is_file():
            raw = json.loads(sp.read_text())
            recs = raw["ok"] if isinstance(raw, dict) and "ok" in raw else (raw if isinstance(raw, list) else [])

    def img(name, cap):
        return f'<figure><img src="plots/{name}" alt=""/><figcaption>{cap}</figcaption></figure>'

    css = """body { font: 16px/1.45 system-ui,sans-serif; max-width: 980px; margin: 24px auto; padding: 0 16px; color:#222; }
h1 { font-size: 22px; } h2 { font-size: 17px; margin-top: 1.8em; border-bottom:1px solid #ddd; }
h3 { font-size: 15px; margin-top: 1.4em; }
.gap { background:#fff8e6; padding:10px 12px; border:1px solid #e6d9a8; }
.ok { background:#eef8ee; padding:10px 12px; border:1px solid #b7ddb7; }
figure { margin: 1em 0 1.6em; } img { max-width:100%; border:1px solid #ddd; }
figcaption { font-size:13px; color:#444; margin-top:6px; }
table { border-collapse:collapse; font-size:13px; } td,th { border:1px solid #ccc; padding:4px 8px; text-align:left; }
code { font-size: 90%; } ul.check li { margin: 0.35em 0; }"""
    head = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<title>Hydride motif spectra — DFTB+ vs MMFF</title>
<style>{css}</style></head><body>
<h1>Hydride motif spectra — DFTB+ 3ob-3-1 vs MMFF</h1>
<p>Write-up: <a href="../../../doc/Topics/FTIR_Nanocrystals/HydrideMotif_spectra.md">HydrideMotif_spectra.md</a>.
L2 Hessians: <code>/home/prokop/SIMULATIONS/SiNCs/DFTB/</code> (DFTB+ own min) and <code>tests/tSiNCs/OUT_dftb_vs_mmff/</code> (MMFF FIRE then Hessian at that MMFF’s own minimum).</p>
{_l2_section_html(l2_recs)}
"""
    if not recs:
        html = head + "</body></html>"
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "index.html").write_text(html)
        print(f"  wrote {(OUT/'index.html').relative_to(REPO)}")
        return
    by_id = {r["id"]: r for r in recs}

    def num(cid, k):
        v = by_id[cid][k]
        return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else (f"{v:.1f}" if isinstance(v, float) else str(v))

    def _fk(v):
        return "" if v is None or not np.isfinite(v) else f"{float(v):.1f}"
    q2_rows = "".join(
        f"<tr><td>{r['id']}</td><td>{r['n_XH']}</td><td>{r['n_XH2']}</td><td>{r['mean_XH_cm1']:.1f}</td><td>{r['mean_XH2_cm1']:.1f}</td>"
        f"<td>{r['mean_XH2_cm1']-r['mean_XH_cm1']:+.1f}</td><td>{_fk(r.get('mean_XH_kusova'))}</td><td>{_fk(r.get('mean_XH2_kusova'))}</td></tr>"
        for r in recs
    )
    gallery = []
    for cid, _d, elem, motif, note in CRYSTALS:
        gallery.append(f"<h3>{cid}</h3><p>{elem}, {motif}: {note}.</p>")
        gallery.append(img(f'map_chem_{cid}.png', f'{cid}: ridge_110 atom colors (legend = PDOS stack).'))
        if (PLOT / f"l1_stack_{cid}.png").is_file():
            gallery.append(img(f'l1_stack_{cid}.png', f'{cid}: stacked atom-group PDOS; black line = total DOS.'))
        if (PLOT / f"nbhd_{cid}.png").is_file():
            gallery.append(img(f'nbhd_{cid}.png', f'{cid}: neighborhood stretch PDOS (same tags).'))
    miss = [c for c, *_ in CRYSTALS if c not in by_id]
    fail_box = ""
    if miss:
        fail_box = "<div class=\"gap\"><b>L1 MMFF Hessians not at a minimum — not spectra:</b> " + ", ".join(miss) + ". Atom maps and halogen XYZ still use <code>ridge_110</code> on the stored geometry. Do not use leftover <code>q1_</code>/<code>q2_</code>/<code>q3_</code> PNGs for those crystals.</div>"
    q1_phys = ""
    if "cube_C" in by_id and "octa_C" in by_id:
        q1_phys = f"<p><b>Physics:</b> cube {by_id['cube_C']['n_XH']} XH + {by_id['cube_C']['n_XH2']} XH₂; octa {by_id['octa_C']['n_XH']} XH + {by_id['octa_C']['n_XH2']} XH₂ (identical graph for Si).</p>"
    html = head + f"""<div class="gap"><b>L1 MMFF atlas below</b> (small Wulff cubes/octa, not the L2 DFTB crystals). Hessian = MMFF FD from <code>relaxed.mol2</code>, NonBonded off. Unscaled frequencies. 5/7-ring bonds from mol2. Face/edge = <code>ridge_110</code>.</div>
{fail_box}

<h2>Checklist</h2>
<ul class="check">
<li><b>Q1. {{100}} vs {{111}}</b> — cube vs octahedron, C and Si. Geometry and hydride counts: yes (atlas). Spectra: L1 MMFF stretch PDOS below; older L2 C total-DOS (octa/trunc/rhombic, no cube).</li>
<li><b>Q2. XH vs XH₂</b> — project each mode onto X–H stretch, split by nH on the heavy from topology. Cube is mixed 12+12, not pure XH₂; octa is 24 XH + 6 XH₂.</li>
<li><b>Q3. 5-ring vs 7-ring</b> — ΔS vs the same parent (collapse or insert one XH₂), not a random ensemble.</li>
</ul>

<h2>0. Where the groups sit (all L0/L1 — <code>ridge_110</code>)</h2>
<p>Atom color = neighborhood of the X–H oscillator, same as the stacked PDOS.
Blue = XH@{{111}} face, red = XH₂@{{100}} face, orange = ⟨110⟩-extremum edge.
Jmol: <code>OUT_motif_spectra/xyz/l1_halogen_*.xyz</code> — F=XH face, Cl=XH edge, Br=XH₂ face, I=XH₂ edge.</p>
{"".join(gallery)}

<h2>2. Where those frequencies live on the particle</h2>
<p>Yellow→red = stretch participation in a named window. Black halo = 5/7-ring atoms. If a spectral feature belongs to {{100}} SiH₂, the red blob must sit on the cube faces, not the core.</p>
{img('loc_cube_Si_SiHwin.png', 'Cube Si, 2040–2100 cm⁻¹ (Kusová SiH): edge/corner XH, not the {100} terraces.')}
{img('loc_cube_Si_SiH2win.png', 'Cube Si, 2100–2140 cm⁻¹ (Kusová SiH₂): {100} face hydrides.')}
{img('loc_cube_Si_hi2245.png', 'Cube Si, 2220–2270 cm⁻¹ (sharp MMFF band): collective {100} terrace — the mode 5/7-rings destroy.')}
{img('loc_octa_Si_SiHwin.png', 'Octa Si, 2040–2100: {111} face XH.')}
{img('loc_cube5_Si_hi2245.png', '5-ring cube, 2220–2270: the coupled terrace band is gone/shifted; remaining weight is not on the 5-ring halo.')}
{img('loc_cube7_Si_hi2245.png', '7-ring cube, same window.')}

<h2>Q1 — {{100}} cube vs {{111}} octahedron (totals)</h2>
{q1_phys}
<p><b>What the plots show:</b> Si cube has a double peak across the Kusová SiH / SiH₂ windows plus a sharp coupled band ~2245 cm⁻¹. Si octa is XH-dominated, strongest weight ~2010 cm⁻¹ (slightly below the Kusová SiH window) with almost no SiH₂-window intensity. C default-MMFF CH stretches sit in this same 1800–2400 window, not at 2815–2975; cube C is much stronger here than octa C (more CH₂ oscillators + a sharp ~2320 cm⁻¹ band).</p>
{img('hydride_counts_L1.png', 'Topology hydride counts on L1 Wulff crystals. Cube ≈ equal XH and XH₂; octa is XH-dominated. 5-ring removes one XH₂; 7-ring adds one.')}
{img('q1_cube_vs_octa_C.png', 'C L1: total X–H stretch PDOS. Cube vs octa. Shaded bands = Kusová C–H windows (MMFF default peaks sit lower).')}
{img('q1_cube_vs_octa_Si.png', 'Si L1: same contrast. Shaded = Kusová SiH / SiH₂ windows.')}
<p>Older L2 C Wulff (OUT_dftb_vs_mmff) had octa / truncated / rhombic Hessians already — total mode DOS, not typed:</p>
{img('older_L2_C_faces_total_DOS.png', 'L2-size C: default MMFF (top) and scaled MMFF (bottom). Truncated = mixed {{100}}+{{111}}; rhombic = {{110}}. No cube in that folder.')}

<h2>Q2 — XH vs XH₂ on the same crystal</h2>
<p><b>Physics:</b> q<sub>b</sub> = r̂ · (u<sub>H</sub> − u<sub>X</sub>), W = Σ q² over bonds of that class.
Experiment: SiH₂ sits above SiH (~2100–2140 vs ~2040–2100 cm⁻¹). MMFF default may or may not reproduce the ordering.</p>
<table><tr><th>crystal</th><th>n XH</th><th>n XH2</th><th>&lt;XH&gt; 1800–2400</th><th>&lt;XH2&gt; 1800–2400</th><th>Δ</th><th>&lt;XH&gt; in 2040–2100</th><th>&lt;XH2&gt; in 2100–2140</th></tr>
{q2_rows}</table>
<p>Full-window means are pulled by the sharp ~2245 cm⁻¹ cube band (mostly XH projection). Inside the Kusová windows on <b>cube Si</b>, XH peaks ~2080 and XH₂ ~2110 — the experimental order. Octa Si has almost no XH₂-window weight (6 XH₂ vs 24 XH). Default MMFF does <b>not</b> globally place XH₂ above XH; FFfit on spheres was the calibration that reduced subtype RMSE 41→32 cm⁻¹.</p>
{img('q2_cube_C_XH_XH2.png', 'Cube C: XH vs XH₂ stretch PDOS on the same Hessian.')}
{img('q2_octa_C_XH_XH2.png', 'Octa C: XH-rich; XH₂ is a minority shoulder.')}
{img('q2_cube_Si_XH_XH2.png', 'Cube Si.')}
{img('q2_octa_Si_XH_XH2.png', 'Octa Si.')}
{img('q2_cages.png', 'L0 cages (adamantane / Si10H16): mixed XH and XH₂ on one Td molecule — the smallest controlled split.')}

<h2>Q3 — 5-ring (collapse XH₂) vs 7-ring (insert XH₂)</h2>
<p><b>Physics:</b> one extra (or missing) surface XH₂ plus a 5- or 7-membered ring. Characteristic features are <i>differences vs parent</i>, not new isolated peaks in a large mix.</p>
<p><b>What the plots show:</b> both collapse (5-ring) and insert (7-ring) on the Si cube <i>break the sharp ~2245 cm⁻¹ coupled {{100}} band</i> (large negative ΔS). That is a ring-defect fingerprint of a collective facet mode, not a new isolated XH₂ oscillator. 5-ring also loses weight near 2085 cm⁻¹ (parent SiH-window peak). Several remote peaks (1965, 2175, …) are unchanged (ΔS ≈ 0).</p>
{img('q3_cube5_C.png', 'Cube C: 5-ring − parent. Expect a loss of XH₂ stretch weight.')}
{img('q3_cube7_C.png', 'Cube C: 7-ring − parent. Expect a gain of XH₂ stretch weight.')}
{img('q3_cube5_Si.png', 'Cube Si 5-ring ΔS: loss at ~2085 and especially ~2245 cm⁻¹.')}
{img('q3_cube7_Si.png', 'Cube Si 7-ring ΔS: the same 2245 cm⁻¹ coupled band collapses; not a simple +XH₂ peak.')}
{img('q3_octa5_C.png', 'Octa C 5-ring ΔS (parent already XH-rich; one XH₂ is a small fraction).')}
{img('q3_octa7_C.png', 'Octa C 7-ring ΔS.')}
{img('q3_octa5_Si.png', 'Octa Si 5-ring ΔS.')}
{img('q3_octa7_Si.png', 'Octa Si 7-ring ΔS.')}

<h2>Older calculations in tests/tSiNCs (what they can and cannot answer)</h2>
<table>
<tr><th>Tree</th><th>What it is</th><th>Q1 100/111</th><th>Q2 XH/XH2</th><th>Q3 5/7-ring</th></tr>
<tr><td>OUT_chem_atlas</td><td>30 Wulff geometries, 14 MMFF-relaxed</td><td>shapes + counts yes; spectra = this page (L1)</td><td>counts yes</td><td>geometries yes; ΔS = this page</td></tr>
<tr><td>OUT_dftb_vs_mmff</td><td>L2-size C octa/trunc/rhombic MMFF+DFTB freqs</td><td>111 vs 110 vs mixed; <b>no cube</b></td><td>total DOS only</td><td>no</td></tr>
<tr><td>OUT_FFfit_plots</td><td>PySCF Hessian fit on SiH4 + Si spheres</td><td>no (spheres)</td><td>subtype RMSE 41→32 cm⁻¹ — closest old typed split</td><td>no</td></tr>
<tr><td>OUT_nc_ensemble_v2</td><td>8 random diamond spheres, full 01→05</td><td>no</td><td>mixed hydrides</td><td>no</td></tr>
<tr><td>OUT_nc_consolidated</td><td>1 sphere + accumulated spectrum</td><td>no</td><td>no</td><td>no</td></tr>
<tr><td>OUT_nc_atlas</td><td>pre-Wulff Miller/sphere C SVG</td><td>shapes only</td><td>no Hessian</td><td>no</td></tr>
<tr><td>OUT_small_nc</td><td>Si/C spheres &lt;100 atoms</td><td>no</td><td>FFfit geometries</td><td>no</td></tr>
<tr><td>OUT_XRD</td><td>Debye powder on diamond spheres</td><td colspan="3">not FTIR</td></tr>
<tr><td>SiH_CH_FTIRs_Kusova</td><td>experiment PNGs</td><td colspan="3">overlay target, not simulation</td></tr>
<tr><td>fixtures/si_1nm_passivation</td><td>README for 01→05</td><td colspan="3">folders empty on this disk</td></tr>
</table>
{img('older_ensemble_spheres_hydrides.png', 'Random spheres mix XH and XH₂; size series is not a facet experiment.')}
<p>FFfit typed Si overlays (spheres, not Wulff):
<a href="../OUT_FFfit_plots/typing_comparison_all_Si_hierarchical_joint_v2/">typing_comparison_all_Si_hierarchical_joint_v2</a>.
Kusová:
<a href="../SiH_CH_FTIRs_Kusova/SiH-standard-v01.png">Si–H standards</a> ·
<a href="../SiH_CH_FTIRs_Kusova/ND-fits-v01.png">ND C–H</a>.</p>

<h2>What is still missing</h2>
<ul>
<li>MMFF default stretch frequencies are offset from Kusová; FFfit / QM Hessians on these same Wulff mol2 files are the next calibration.</li>
<li>No SiH₃ / CH₃ (prune drops corners). No L2/L3 Hessian (too big for this FD pass).</li>
<li>IR intensities are stretch-projection DOS, not bond-dipole FTIR.</li>
</ul>
</body></html>
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(html)
    print(f"  wrote {(OUT/'index.html').relative_to(REPO)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-hessian", action="store_true", help="reuse existing 04_hessian.npz and cached L2 MMFF eigensolves")
    ap.add_argument("--l2-only", action="store_true", help="only DFTB+ vs MMFF stacked PDOS on the 3 L2 C crystals")
    ap.add_argument("--pyscf", default=None, help="PySCF vib results dir (canonical: .../pyscf_vib_results; luna_/tight_ prefixes resolved)")
    ap.add_argument("--pyscf-out", default=None, help="output dir for --pyscf (default tests/tSiNCs/OUT_pyscf_jobs)")
    ap.add_argument("--dftb-l1", default=None, help="DFTB+ L1 job root ({crystal}_{sk}/vibrations.tag); with --pyscf, per-crystal comparison")
    ap.add_argument("--compare-out", default=None, help="output dir for --pyscf --dftb-l1 (default tests/tSiNCs/OUT_dftb_vs_pyscf_l1)")
    ap.add_argument("--facet", choices=("ridge_110", "miller_111", "xh_align", "wulff"), default="ridge_110", help="face/edge: ⟨110⟩ extrema (default), X–H vs sitting ⟨111⟩, hydride-neighbor terrace, or Wulff COM")
    ap.add_argument("--align-cos", type=float, default=0.90, help="miller_111 / xh_align cosine gate (0.90 = 10% from 1, 0.95 = 5%)")
    ap.add_argument("--ridge-below", type=float, default=0.5, help="ridge_110: Å below a ⟨110⟩ extremum counted as edge")
    args = ap.parse_args()
    os.chdir(REPO)
    OUT.mkdir(parents=True, exist_ok=True)
    PLOT.mkdir(parents=True, exist_ok=True)
    if args.pyscf:
        facet = args.facet
        if facet == "ridge_110":
            print("  L1 PySCF: ridge_110 tags every hydride as edge on these sizes; using facet=wulff (COM support)")
            facet = "wulff"
        if args.dftb_l1:
            outc = Path(args.compare_out) if args.compare_out else TEST / "OUT_dftb_vs_pyscf_l1"
            print(f"== L1 PySCF vs DFTB+  {args.pyscf}  vs  {args.dftb_l1}  -> {outc}")
            analyze_l1_pyscf_vs_dftb(args.pyscf, args.dftb_l1, outc, facet_mode=facet)
            print("done", outc)
            return
        outp = Path(args.pyscf_out) if args.pyscf_out else TEST / "OUT_pyscf_jobs"
        print(f"== PySCF neighborhood PDOS  {args.pyscf}  -> {outp}")
        analyze_pyscf_jobs(args.pyscf, outp, facet_mode=facet)
        print("done", outp)
        return
    if args.dftb_l1:
        raise SystemExit("--dftb-l1 requires --pyscf (same 12 crystals, stacked PDOS comparison)")
    if args.l2_only:
        print("== L2 DFTB+ 3ob-3-1 vs MMFF (3 C crystals)")
        l2 = analyze_l2_dftb_vs_mmff(skip_cache=args.skip_hessian, facet_mode=args.facet, align_cos=args.align_cos, ridge_below_A=args.ridge_below)
        write_html([], None, l2_recs=l2)
        print("done", OUT)
        return
    recs = []
    hess_fail = []
    for row in CRYSTALS:
        print(f"== {row[0]}")
        try:
            recs.append(analyze_crystal(*row, args.skip_hessian))
        except RuntimeError as e:
            print(f"  SKIP spectrum {row[0]}: {e}")
            hess_fail.append((row[0], str(e)))
    (OUT / "stats.json").write_text(json.dumps({"ok": recs, "hessian_fail": [{"id": a, "error": b} for a, b in hess_fail]}, indent=2))
    ok_ids = {r["id"] for r in recs}
    if recs:
        plot_hydride_bars(recs, "hydride_counts_L1.png")
    def _have(*cids):
        return all(c in ok_ids for c in cids)
    if _have("cube_C", "octa_C"):
        plot_overlay(["cube_C", "octa_C"], ["cube {100}", "octa {111}"], "Q1  C  cube vs octa — X–H stretch PDOS", "q1_cube_vs_octa_C.png", "C")
    if _have("cube_Si", "octa_Si"):
        plot_overlay(["cube_Si", "octa_Si"], ["cube {100}", "octa {111}"], "Q1  Si  cube vs octa — X–H stretch PDOS", "q1_cube_vs_octa_Si.png", "Si")
    if _have("cube_C"):
        plot_xh_split("cube_C", "Q2  cube C — XH vs XH₂", "q2_cube_C_XH_XH2.png", "C")
    if _have("octa_C"):
        plot_xh_split("octa_C", "Q2  octa C — XH vs XH₂", "q2_octa_C_XH_XH2.png", "C")
    if _have("cube_Si"):
        plot_xh_split("cube_Si", "Q2  cube Si — XH vs XH₂", "q2_cube_Si_XH_XH2.png", "Si")
    if _have("octa_Si"):
        plot_xh_split("octa_Si", "Q2  octa Si — XH vs XH₂", "q2_octa_Si_XH_XH2.png", "Si")
    if _have("L0_adamantane", "L0_Si10H16"):
        plot_overlay(["L0_adamantane", "L0_Si10H16"], ["adamantane C", "Si10H16"], "Q2  L0 cages — total stretch PDOS", "q2_cages.png", "C", keys=("w_XH", "w_XH2"))
    for parent, child, title, fname, elem in (
        ("cube_C", "cube5_C", "Q3  cube C  5-ring − parent", "q3_cube5_C.png", "C"),
        ("cube_C", "cube7_C", "Q3  cube C  7-ring − parent", "q3_cube7_C.png", "C"),
        ("cube_Si", "cube5_Si", "Q3  cube Si  5-ring − parent", "q3_cube5_Si.png", "Si"),
        ("cube_Si", "cube7_Si", "Q3  cube Si  7-ring − parent", "q3_cube7_Si.png", "Si"),
        ("octa_C", "octa5_C", "Q3  octa C  5-ring − parent", "q3_octa5_C.png", "C"),
        ("octa_C", "octa7_C", "Q3  octa C  7-ring − parent", "q3_octa7_C.png", "C"),
        ("octa_Si", "octa5_Si", "Q3  octa Si  5-ring − parent", "q3_octa5_Si.png", "Si"),
        ("octa_Si", "octa7_Si", "Q3  octa Si  7-ring − parent", "q3_octa7_Si.png", "Si"),
    ):
        if _have(parent, child):
            plot_delta(parent, child, title, fname, elem)
        else:
            print(f"  SKIP {fname}: Hessian not at minimum")
    plot_older_faces()
    plot_older_ensemble_hydrides()
    print("== spatial + neighborhood (L0/L1 MMFF, all atlas crystals)")
    for cid, d, elem, motif, _note in CRYSTALS:
        print(f"  spatial {cid}")
        analyze_spatial_l1(cid, d, elem, STRETCH_WIN[elem] if elem == "Si" else (1800, 2400), motif, facet_mode=args.facet, ridge_below_A=args.ridge_below)
    print("== L2 DFTB+ 3ob-3-1 vs MMFF (3 C crystals)")
    l2 = analyze_l2_dftb_vs_mmff(skip_cache=args.skip_hessian, facet_mode=args.facet, align_cos=args.align_cos, ridge_below_A=args.ridge_below)
    write_html(recs, None, l2_recs=l2)
    print("done", OUT)


if __name__ == "__main__":
    main()
