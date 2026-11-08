#!/usr/bin/env python3
# === AUTO-DOC BEGIN ===
"""
plot_pyscf_vib_results.py - Plot and analyze PySCF vibration results (Hessian + spectra).

## Input format contract

Each case directory under ``results/`` must contain:

  frequencies_cm1.npy   (3N,)     float64 or complex128   [cm⁻¹; complex = imaginary freqs]
  hessian.npy           (N,N,3,3) float64                 [Hartree/Bohr²]
  modes.npy             (n_vib,N,3) float64               [n_vib = 3N-6; vibrational only]
  masses.npy            (N,)      int64                   [amu]
  relaxed.xyz           standard XYZ geometry
  status.json           optional metadata (E_hf, status)

The Hessian is plotted in atomic-block layout: [atom0_xyz, atom1_xyz, ...]
with 3×3 blocks on the diagonal, not [all_x, all_y, all_z].

## Usage
  python3 plot_pyscf_vib_results.py /path/to/results --noshow
  python3 plot_pyscf_vib_results.py /path/to/results --overlay --atlas-overlays --no-hessian --noshow
  python3 plot_pyscf_vib_results.py /path/to/results --cases SiH4 Si_R3p8
  python3 plot_pyscf_vib_results.py /path/to/results --size-series Si
"""
# === AUTO-DOC END ===
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from pyBall.FFfit_utils import as_signed_wavenumbers_cm1, gaussian_spectrum, pyscf_vib_job_case_name, list_pyscf_vib_case_dirs  # noqa: E402

# Experimental hydride windows — same numbers as plot_hydride_motif_spectra.py (Kusová).
KUSOVA = {"SiH": (2040.0, 2100.0), "SiH2": (2100.0, 2140.0), "SiH3": (2140.0, 2170.0), "CH": (2815.0, 2900.0), "CH2": (2900.0, 2975.0)}
ATLAS_OVERLAYS = [
    ("overlay_C", "PySCF PBE/ccECP-cc-pVDZ — diamond L1", ["cube_C", "cube_5ring_C", "cube_7ring_C", "octahedron_C", "octahedron_5ring_C", "octahedron_7ring_C"]),
    ("overlay_Si", "PySCF PBE/ccECP-cc-pVDZ — Si L1", ["cube_Si", "cube_5ring_Si", "cube_7ring_Si", "octahedron_Si", "octahedron_5ring_Si", "octahedron_7ring_Si"]),
    ("q1_cube_vs_octa_C", "Q1 {100} cube vs {111} octa — C", ["cube_C", "octahedron_C"]),
    ("q1_cube_vs_octa_Si", "Q1 {100} cube vs {111} octa — Si", ["cube_Si", "octahedron_Si"]),
    ("q3_cube_rings_C", "Q3 cube C: parent / 5-ring / 7-ring", ["cube_C", "cube_5ring_C", "cube_7ring_C"]),
    ("q3_cube_rings_Si", "Q3 cube Si: parent / 5-ring / 7-ring", ["cube_Si", "cube_5ring_Si", "cube_7ring_Si"]),
    ("q3_octa_rings_C", "Q3 octa C: parent / 5-ring / 7-ring", ["octahedron_C", "octahedron_5ring_C", "octahedron_7ring_C"]),
    ("q3_octa_rings_Si", "Q3 octa Si: parent / 5-ring / 7-ring", ["octahedron_Si", "octahedron_5ring_Si", "octahedron_7ring_Si"]),
]


import matplotlib
if "--noshow" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


# ============================================================
# Data loading
# ============================================================

def load_case(case_dir):
    """Load all vibration data for one case directory."""
    d = Path(case_dir)
    data = {}
    for key, fname in [('freqs', 'frequencies_cm1.npy'), ('hessian', 'hessian.npy'),
                       ('modes', 'modes.npy'), ('masses', 'masses.npy')]:
        p = d / fname
        if not p.exists():
            raise FileNotFoundError(f"Missing {fname} in {case_dir}")
        data[key] = np.load(str(p))
    xyz_p = d / 'relaxed.xyz'
    if xyz_p.exists():
        data['symbols'], data['positions'] = read_xyz(str(xyz_p))
    json_p = d / 'status.json'
    if json_p.exists():
        with open(json_p) as f:
            data['status'] = json.load(f)
    data['name'] = pyscf_vib_job_case_name(d)[0]
    return data


def read_xyz(path):
    """Read XYZ file, return (symbols list, positions array)."""
    with open(path) as f:
        n = int(f.readline().strip())
        f.readline()  # comment
        syms, pos = [], []
        for _ in range(n):
            parts = f.readline().split()
            syms.append(parts[0])
            pos.append([float(x) for x in parts[1:4]])
    return syms, np.array(pos, dtype=np.float64)


def discover_cases(results_dir):
    """One directory per canonical crystal (luna_/tight_/modefollow_ collapsed)."""
    return [str(p) for p in list_pyscf_vib_case_dirs(results_dir)]


# ============================================================
# Frequency utilities
# ============================================================

def filter_real_freqs(freqs, threshold=10.0):
    """Return real positive frequencies above threshold. Imaginary freqs (complex) are excluded."""
    w = as_signed_wavenumbers_cm1(freqs)
    return w[w > threshold]


def get_imaginary_freqs(freqs, threshold=10.0):
    """Return |imaginary frequencies| (PySCF stores unstable modes as +i|ν|; signed helper → negative)."""
    w = as_signed_wavenumbers_cm1(freqs)
    return -w[w < -threshold]


def freq_diagnostics(freqs):
    """Signed cm⁻¹ diagnostics. PySCF jobs store 3N-6 (rigid already projected) so n_rigid≈0 is expected."""
    w = as_signed_wavenumbers_cm1(freqs)
    n_imag = int(np.sum(w < -10.0))
    n_soft = int(np.sum((w >= -10.0) & (w <= 10.0)))
    n_vib = int(np.sum(w > 10.0))
    worst = float(w.min()) if w.size else float("nan")
    nu1 = float(w[w > 10.0].min()) if n_vib else float("nan")
    nu_max = float(w.max()) if w.size else float("nan")
    return {"w": w, "n_imag": n_imag, "n_soft": n_soft, "n_vib": n_vib, "worst": worst, "nu1": nu1, "nu_max": nu_max}


def shade_kusova(ax, name):
    """Shade Kusová hydride windows. Cite plot_hydride_motif_spectra.py."""
    keys = ("CH", "CH2") if name.endswith("_C") else ("SiH", "SiH2", "SiH3")
    colors = {"CH": "#93c5fd", "CH2": "#fde68a", "SiH": "#93c5fd", "SiH2": "#fde68a", "SiH3": "#fbcfe8"}
    for k in keys:
        lo, hi = KUSOVA[k]
        ax.axvspan(lo, hi, color=colors[k], alpha=0.22, lw=0, zorder=0)
        ax.text(0.5 * (lo + hi), ax.get_ylim()[1] * 0.92, k, ha="center", va="top", fontsize=7, color="#444")


# ============================================================
# Plotting: Hessian matrix
# ============================================================

def plot_hessian(data, outdir=None, noshow=False):
    """Plot Hessian matrix as 2D imshow with symmetric color scale."""
    hess = data['hessian']  # (N, N, 3, 3)
    N = hess.shape[0]
    hess_2d = hess.reshape(N * 3, N * 3)
    vmax = np.abs(hess_2d).max()
    if vmax < 1e-30:
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(hess_2d, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
    ax.set_title(f'Hessian — {data["name"]}  (N={N}, 3N={3*N})\n'
                 f'Layout: [atom0 xyz, atom1 xyz, ...]  3×3 atomic blocks on diagonal\n'
                 f'Range: [{hess_2d.min():.2e}, {hess_2d.max():.2e}]')
    ax.set_xlabel('Atom j (x,y,z blocks)')
    ax.set_ylabel('Atom i (x,y,z blocks)')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Hartree/Bohr²')

    # Grid lines at atom boundaries (3x3 blocks)
    for i in range(1, N):
        ax.axhline(i * 3 - 0.5, color='gray', lw=0.4, alpha=0.4)
        ax.axvline(i * 3 - 0.5, color='gray', lw=0.4, alpha=0.4)
    # Atom number ticks at block centers
    if N <= 30:
        ax.set_xticks(np.arange(N) * 3 + 1)
        ax.set_xticklabels(np.arange(N), fontsize=7)
        ax.set_yticks(np.arange(N) * 3 + 1)
        ax.set_yticklabels(np.arange(N), fontsize=7)

    plt.tight_layout()
    if outdir:
        out = Path(outdir) / f'{data["name"]}_hessian.png'
        fig.savefig(str(out), dpi=150)
        print(f"  Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


def plot_hessian_offdiag(data, outdir=None, noshow=False):
    """Plot off-diagonal blocks (inter-atom coupling) as a heatmap of block norms."""
    hess = data['hessian']  # (N, N, 3, 3)
    N = hess.shape[0]
    block_norms = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            block_norms[i, j] = np.linalg.norm(hess[i, j])

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(block_norms, cmap='hot_r', aspect='equal')
    ax.set_title(f'Inter-atom Hessian block norms — {data["name"]}  (N={N})\n'
                 f'Each cell = ‖H[i,j]‖ for the 3×3 block')
    ax.set_xlabel('Atom j')
    ax.set_ylabel('Atom i')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='‖H_ij‖')

    plt.tight_layout()
    if outdir:
        out = Path(outdir) / f'{data["name"]}_hessian_blocks.png'
        fig.savefig(str(out), dpi=150)
        print(f"  Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


# ============================================================
# Plotting: Vibration spectrum
# ============================================================

def plot_spectrum(data, xmin=0.0, xmax=3500, width=20.0, outdir=None, noshow=False, suffix="spectrum"):
    """Two-panel figure: stick spectrum + Gaussian-broadened."""
    freqs = data['freqs']
    diag = freq_diagnostics(freqs)
    real_f = filter_real_freqs(freqs, threshold=10.0)
    imag_f = get_imaginary_freqs(freqs, threshold=10.0)
    warn = f"  — NOT a minimum: {diag['n_imag']} imag, worst {diag['worst']:.1f} cm⁻¹" if diag["n_imag"] else ""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(f'Vibrational spectrum — {data["name"]}  '
                 f'(N={len(data.get("symbols", []))}, {len(real_f)} real modes'
                 + (f', {len(imag_f)} imaginary' if len(imag_f) else '') + f'){warn}',
                 fontsize=12, color=("#b91c1c" if diag["n_imag"] else "#111"))
    x_lo = float(xmin)
    if len(imag_f) and xmin < 50.0:
        x_lo = min(x_lo, float(-np.max(imag_f)) - 40.0)
    x = np.arange(x_lo, xmax + 1, 0.5)
    ax1.vlines(real_f, 0, 1, colors='#2563eb', linewidth=0.6, alpha=0.75, label='real')
    if len(imag_f):
        ax1.vlines(-imag_f, 0, 0.5, colors='#dc2626', linewidth=0.8, alpha=0.85, label='imaginary')
        ax1.axvline(0.0, color='#111', lw=0.6, ls='--', alpha=0.5)
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)
    shade_kusova(ax1, data["name"])
    spec = gaussian_spectrum(real_f, x, width)
    if spec.max() > 0:
        spec /= spec.max()
    ax2.plot(x, spec, color='#1d4ed8', lw=1.4, label=f'Gaussian (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec, color='#3b82f6', alpha=0.15)
    shade_kusova(ax2, data["name"])
    ax2.set_xlim(x_lo, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.25, ls=':')
    plt.tight_layout()
    out = None
    if outdir:
        out = Path(outdir) / f'{data["name"]}_{suffix}.png'
        fig.savefig(str(out), dpi=150)
        print(f"  Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)
    return out


# ============================================================
# Plotting: Overlay spectra across cases
# ============================================================

def plot_overlay(all_data, xmin=0.0, xmax=3500, width=20.0, outdir=None, noshow=False, name="overlay_spectra", title=None):
    """Overlay broadened spectra from all cases."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    fig.suptitle(title or "Vibrational spectra overlay — all cases", fontsize=14)
    imag_lo = 0.0
    if xmin < 50.0:
        for data in all_data:
            imag_f = get_imaginary_freqs(data["freqs"], threshold=10.0)
            if len(imag_f):
                imag_lo = min(imag_lo, float(-np.max(imag_f)))
    x_lo = min(float(xmin), imag_lo - 40.0) if imag_lo < 0 else float(xmin)
    x = np.arange(x_lo, xmax + 1, 0.5)
    colors = plt.cm.tab10(np.linspace(0, 0.9, max(len(all_data), 1)))
    elem_name = all_data[0]["name"] if all_data else ""
    for i, data in enumerate(all_data):
        real_f = filter_real_freqs(data['freqs'], threshold=10.0)
        imag_f = get_imaginary_freqs(data['freqs'], threshold=10.0)
        c = colors[i]
        label = data['name']
        n_imag = freq_diagnostics(data["freqs"])["n_imag"]
        if n_imag:
            label = f"{label} ({n_imag} imag)"
        ax1.vlines(real_f, 0, 1, colors=c, linewidth=0.5, alpha=0.6, label=label)
        if len(imag_f):
            ax1.vlines(-imag_f, 0, 0.45, colors=c, linewidth=0.7, alpha=0.9)
        spec = gaussian_spectrum(real_f, x, width)
        if spec.max() > 0:
            spec /= spec.max()
        ax2.plot(x, spec, color=c, lw=1.2, alpha=0.8, label=label)
    if imag_lo < 0:
        ax1.axvline(0.0, color='#111', lw=0.6, ls='--', alpha=0.5)
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=8, loc='upper right', ncol=2)
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)
    shade_kusova(ax1, elem_name)
    shade_kusova(ax2, elem_name)
    ax2.set_xlim(x_lo, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=8, loc='upper right', ncol=2)
    ax2.grid(True, alpha=0.25, ls=':')
    plt.tight_layout()
    out = None
    if outdir:
        out = Path(outdir) / f"{name}.png"
        fig.savefig(str(out), dpi=150)
        print(f"Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)
    return out


# ============================================================
# Plotting: Compare Hessians side-by-side
# ============================================================

def plot_compare_hessians(all_data, outdir=None, noshow=False):
    """Side-by-side Hessian imshow for all cases."""
    n = len(all_data)
    if n == 0:
        return
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, data in zip(axes, all_data):
        hess = data['hessian']
        N = hess.shape[0]
        hess_2d = hess.reshape(N * 3, N * 3)
        vmax = np.abs(hess_2d).max()
        if vmax < 1e-30:
            vmax = 1.0
        im = ax.imshow(hess_2d, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title(f'{data["name"]}\nN={N}, 3N={3*N}', fontsize=10)
        ax.set_xlabel('3N')
        ax.set_ylabel('3N')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle('Hessian matrix comparison', fontsize=14)
    plt.tight_layout()
    if outdir:
        out = Path(outdir) / 'compare_hessians.png'
        fig.savefig(str(out), dpi=150)
        print(f"Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


# ============================================================
# Plotting: Size-series analysis
# ============================================================

def plot_size_series(all_data, prefix, outdir=None, noshow=False):
    """Plot frequency trends vs particle size for cases starting with prefix."""
    series = [d for d in all_data if d['name'].startswith(prefix)]
    if len(series) < 2:
        print(f"Not enough cases with prefix '{prefix}' for size series (found {len(series)})")
        return

    # Extract radius from name (e.g., Si_R3p8 -> 3.8)
    import re
    sizes = []
    for d in series:
        m = re.search(r'R(\d+)p(\d+)', d['name'])
        if m:
            r = float(f"{m.group(1)}.{m.group(2)}")
        else:
            r = float(d.get('status', {}).get('natoms', 0))
        sizes.append((r, d))
    sizes.sort(key=lambda x: x[0])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f'Size-dependent vibration trends — {prefix} series', fontsize=14)

    colors = plt.cm.viridis(np.linspace(0, 0.9, len(sizes)))
    x = np.arange(0, 2501, 0.5)

    for i, (r, data) in enumerate(sizes):
        real_f = filter_real_freqs(data['freqs'], threshold=10.0)
        spec = gaussian_spectrum(real_f, x, 25.0)
        if spec.max() > 0:
            spec /= spec.max()
        ax1.plot(x, spec, color=colors[i], lw=1.0, alpha=0.8, label=f'R={r:.1f}Å (N={len(data.get("symbols", []))})')

    ax1.set_ylabel('Intensity (normalized)')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(True, alpha=0.25, ls=':')

    # Peak frequency vs size
    peak_freqs = []
    rs = []
    for r, data in sizes:
        real_f = filter_real_freqs(data['freqs'], threshold=50.0)
        if len(real_f):
            peak_freqs.append(np.median(real_f))
            rs.append(r)
    ax2.plot(rs, peak_freqs, 'o-', color='#dc2626', lw=1.5, markersize=6)
    ax2.set_xlabel('Radius (Å)')
    ax2.set_ylabel('Median frequency (cm⁻¹)')
    ax2.grid(True, alpha=0.25, ls=':')

    plt.tight_layout()
    if outdir:
        out = Path(outdir) / f'size_series_{prefix}.png'
        fig.savefig(str(out), dpi=150)
        print(f"Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


# ============================================================
# Summary table
# ============================================================

def print_summary(all_data):
    """Print summary table of all cases. Imaginary modes mean the Hessian is not at a stationary point."""
    print(f"\n{'='*110}")
    print(f"  SUMMARY: {len(all_data)} cases   (freqs are 3N-6; n_soft≈0 expected; n_imag>0 = NOT a minimum)")
    print(f"{'='*110}")
    print(f"  {'Case':<24s} {'N':>4s} {'n_vib':>6s} {'n_imag':>7s} {'n_soft':>6s} {'worst':>8s} {'ν1':>8s} {'νmax':>8s} {'E(Ha)':>14s} {'ok?':>4s}")
    print(f"  {'-'*24} {'-'*4} {'-'*6} {'-'*7} {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*14} {'-'*4}")
    n_bad = 0
    for d in all_data:
        N = d['hessian'].shape[0] if d.get('hessian') is not None else len(d.get('symbols', []))
        diag = freq_diagnostics(d['freqs'])
        E = d.get('status', {}).get('energy_Ha', float('nan'))
        ok = "NO" if diag["n_imag"] else "yes"
        if diag["n_imag"]:
            n_bad += 1
        print(f"  {d['name']:<24s} {N:4d} {diag['n_vib']:6d} {diag['n_imag']:7d} {diag['n_soft']:6d} {diag['worst']:8.1f} {diag['nu1']:8.1f} {diag['nu_max']:8.1f} {E:14.8f} {ok:>4s}")
    print()
    if n_bad:
        print(f"  WARNING: {n_bad}/{len(all_data)} cases have |ν|>10 imaginary modes.")
        print("  Those Hessians are not at a stationary point of this energy. See Hessian_at_own_minimum.md")
        print("  Plots still written (red sticks = imaginaries). Do not treat flagged Si/C curves as FTIR peak positions.\n")


def write_index_html(outdir, all_data, overlay_files, source_dir):
    """L2 review page: provenance + per-case spectra + motif overlays."""
    rows = []
    for d in all_data:
        diag = freq_diagnostics(d["freqs"])
        E = d.get("status", {}).get("energy_Ha", float("nan"))
        flag = "NOT min" if diag["n_imag"] else "ok"
        cls = "bad" if diag["n_imag"] else "ok"
        rows.append(f"<tr class='{cls}'><td>{d['name']}</td><td>{len(d.get('symbols', []))}</td><td>{diag['n_vib']}</td>"
                    f"<td>{diag['n_imag']}</td><td>{diag['worst']:.1f}</td><td>{diag['nu1']:.1f}</td><td>{diag['nu_max']:.1f}</td>"
                    f"<td>{E:.6f}</td><td>{flag}</td></tr>")
    figs = []
    for name, cap in (("stacked_C", "Diamond L1 — full spectrum, one row per crystal (0–3500 cm⁻¹)"),
                      ("stacked_Si", "Si L1 — full spectrum (red labels = not a minimum)")):
        p = Path(outdir) / f"{name}.png"
        if p.is_file():
            figs.append(f"<h3>{cap}</h3><figure><img src='{p.name}' alt=''/><figcaption>{p.name}</figcaption></figure>")
    for name, title, _ in ATLAS_OVERLAYS:
        p = Path(outdir) / f"{name}.png"
        if p.is_file():
            figs.append(f"<h3>{title}</h3><figure><img src='{p.name}' alt=''/><figcaption>{p.name}</figcaption></figure>")
        p2 = Path(outdir) / f"{name}_stretch.png"
        if p2.is_file():
            figs.append(f"<figure><img src='{p2.name}' alt=''/><figcaption>{p2.name} — hydride stretch window</figcaption></figure>")
    per = []
    for d in all_data:
        diag = freq_diagnostics(d["freqs"])
        tag = " class='bad'" if diag["n_imag"] else ""
        per.append(f"<h3{tag}>{d['name']}</h3>")
        for suf in ("spectrum", "stretch"):
            p = Path(outdir) / f"{d['name']}_{suf}.png"
            if p.is_file():
                per.append(f"<figure><img src='{p.name}' alt=''/><figcaption>{p.name}</figcaption></figure>")
    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<title>PySCF L1 nanocrystal spectra</title>
<style>body {{ font: 16px/1.45 system-ui,sans-serif; max-width: 980px; margin: 24px auto; padding: 0 16px; color:#222; }}
h1 {{ font-size: 22px; }} h2 {{ font-size: 17px; margin-top: 1.8em; border-bottom:1px solid #ddd; }}
h3 {{ font-size: 15px; margin-top: 1.4em; }} h3.bad {{ color:#b91c1c; }}
.gap {{ background:#fff8e6; padding:10px 12px; border:1px solid #e6d9a8; }}
.bad-box {{ background:#fde8e8; padding:10px 12px; border:1px solid #e8b4b4; }}
figure {{ margin: 1em 0 1.6em; }} img {{ max-width:100%; border:1px solid #ddd; }}
figcaption {{ font-size:13px; color:#444; margin-top:6px; }}
table {{ border-collapse:collapse; font-size:13px; }} td,th {{ border:1px solid #ccc; padding:4px 8px; }}
tr.bad td {{ background:#fde8e8; }} tr.ok td {{ background:#eef8ee; }}
code {{ font-size: 90%; }}</style></head><body>
<h1>PySCF PBE nanocrystal vibrational spectra (L1 chem-atlas)</h1>
<p>Source: <code>{source_dir}</code>. Method: PBE / ccECP-cc-pVDZ (from chembook.json). Frequencies are 3N−6 (rigid modes already projected). Shaded bands = Kusová hydride windows (same as <code>plot_hydride_motif_spectra.py</code>).</p>
<div class="bad-box"><b>Stationary-point check:</b> several Si jobs (and <code>cube_C</code>) have imaginary modes with |ν|&gt;10 cm⁻¹. Those Hessians are <b>not</b> at a PySCF minimum — leftover forces. Do not cite flagged stretch peaks as FTIR positions. Protocol: <code>doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md</code>. Red sticks = imaginary frequencies (not dropped).</div>
<h2>Diagnostics</h2>
<table><tr><th>case</th><th>N</th><th>n_vib</th><th>n_imag</th><th>worst</th><th>ν₁</th><th>νmax</th><th>E (Ha)</th><th>ok?</th></tr>
{''.join(rows)}
</table>
<h2>Motif overlays</h2>
{''.join(figs)}
<h2>Per-case spectra</h2>
{''.join(per)}
</body></html>
"""
    out = Path(outdir) / "index.html"
    out.write_text(html)
    print(f"REVIEW: {out}")


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description='Plot and analyze PySCF vibration results')
    ap.add_argument('results_dir', help='Directory with case subdirectories (each containing frequencies_cm1.npy, hessian.npy, etc.)')
    ap.add_argument('--cases', nargs='+', default=None, help='Specific case names to process (default: all)')
    ap.add_argument('--outdir', default=None, help='Output directory for PNG files (default: <results_dir>/plots)')
    ap.add_argument('--xmax', type=float, default=3500, help='Max frequency axis (cm^-1); 3500 covers PBE C–H stretches')
    ap.add_argument('--width', type=float, default=15.0, help='Gaussian broadening width (cm^-1)')
    ap.add_argument('--overlay', action='store_true', help='Overlay spectra from all cases')
    ap.add_argument('--atlas-overlays', action='store_true', help='Chem-atlas motif overlays (cube vs octa, 5/7-ring) + hydride-stretch zooms + index.html')
    ap.add_argument('--compare-hessian', action='store_true', help='Compare Hessians side-by-side')
    ap.add_argument('--size-series', nargs='+', default=None, help='Plot size-dependent trends for cases with given prefix (e.g., Si C_diamond)')
    ap.add_argument('--noshow', action='store_true', help='Save figures only, do not open GUI')
    ap.add_argument('--no-hessian', action='store_true', help='Skip per-case Hessian plots')
    ap.add_argument('--no-spectrum', action='store_true', help='Skip per-case spectrum plots')
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    if not os.path.isdir(results_dir):
        print(f"ERROR: not a directory: {results_dir}")
        sys.exit(1)

    outdir = args.outdir or os.path.join(results_dir, 'plots')
    os.makedirs(outdir, exist_ok=True)

    # Discover and load cases
    case_dirs = discover_cases(results_dir)
    if args.cases:
        wanted = set(args.cases)
        by_case = {pyscf_vib_job_case_name(c)[0]: c for c in case_dirs}
        case_dirs = [by_case[c] for c in wanted if c in by_case]
        missing = wanted - set(by_case)
        if missing:
            raise KeyError(f"cases not found: {sorted(missing)}")

    if not case_dirs:
        print(f"No cases found in {results_dir}")
        sys.exit(1)

    all_data = []
    for cd in case_dirs:
        data = load_case(cd)
        all_data.append(data)
        print(f"  Loaded: {data['name']}  (N={data['hessian'].shape[0]})")

    print_summary(all_data)

    by_name = {d["name"]: d for d in all_data}

    # Per-case plots
    for data in all_data:
        if not args.no_hessian:
            print(f"  Plotting Hessian: {data['name']}")
            plot_hessian(data, outdir=outdir, noshow=args.noshow)
            plot_hessian_offdiag(data, outdir=outdir, noshow=args.noshow)
        if not args.no_spectrum:
            print(f"  Plotting spectrum: {data['name']}")
            plot_spectrum(data, xmax=args.xmax, width=args.width, outdir=outdir, noshow=args.noshow)
            lo, hi = (2600.0, 3300.0) if data["name"].endswith("_C") else (1800.0, 2400.0)
            plot_spectrum(data, xmin=lo, xmax=hi, width=args.width, outdir=outdir, noshow=args.noshow, suffix="stretch")

    # Overlay
    if args.overlay:
        print("  Plotting overlay spectra...")
        plot_overlay(all_data, xmax=args.xmax, width=args.width, outdir=outdir, noshow=args.noshow)

    overlay_files = []
    if args.atlas_overlays:
        for name, title, cases in ATLAS_OVERLAYS:
            subset = [by_name[c] for c in cases if c in by_name]
            missing = [c for c in cases if c not in by_name]
            if missing:
                raise KeyError(f"atlas overlay {name}: missing cases {missing}")
            print(f"  Overlay {name} ({len(subset)} cases)")
            plot_overlay(subset, xmax=args.xmax, width=args.width, outdir=outdir, noshow=args.noshow, name=name, title=title)
            stretch_lo, stretch_hi = (2600.0, 3300.0) if cases[0].endswith("_C") else (1800.0, 2400.0)
            plot_overlay(subset, xmin=stretch_lo, xmax=stretch_hi, width=args.width, outdir=outdir, noshow=args.noshow,
                         name=f"{name}_stretch", title=title + " — hydride stretch")
            overlay_files.append(name)
        write_index_html(outdir, all_data, overlay_files, results_dir)

    # Compare Hessians
    if args.compare_hessian:
        print("  Plotting Hessian comparison...")
        plot_compare_hessians(all_data, outdir=outdir, noshow=args.noshow)

    # Size series
    if args.size_series:
        for prefix in args.size_series:
            print(f"  Plotting size series: {prefix}")
            plot_size_series(all_data, prefix, outdir=outdir, noshow=args.noshow)

    print(f"\nAll plots saved to: {outdir}")


if __name__ == '__main__':
    main()
