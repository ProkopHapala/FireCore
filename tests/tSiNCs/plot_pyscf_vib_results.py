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
  python3 plot_pyscf_vib_results.py /path/to/results
  python3 plot_pyscf_vib_results.py /path/to/results --cases SiH4 Si_R3p8
  python3 plot_pyscf_vib_results.py /path/to/results --overlay
  python3 plot_pyscf_vib_results.py /path/to/results --size-series Si
"""
# === AUTO-DOC END ===
import argparse
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


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
    data['name'] = d.name
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
    """Find all case subdirectories with frequencies_cm1.npy."""
    cases = []
    for d in sorted(os.listdir(results_dir)):
        p = os.path.join(results_dir, d)
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, 'frequencies_cm1.npy')):
            cases.append(p)
    return cases


# ============================================================
# Frequency utilities
# ============================================================

def filter_real_freqs(freqs, threshold=10.0):
    """Return real positive frequencies above threshold. Imaginary freqs (complex) are excluded."""
    if np.iscomplexobj(freqs):
        mask = (freqs.imag == 0) & (freqs.real > threshold)
        return freqs[mask].real
    return freqs[freqs > threshold].real


def get_imaginary_freqs(freqs, threshold=10.0):
    """Return |imaginary frequencies| (stored as pure imaginary in PySCF)."""
    if np.iscomplexobj(freqs):
        mask = freqs.imag > threshold
        return freqs[mask].imag
    return np.array([])


def broaden(freqs, x, width=20.0):
    """Sum of Gaussians centered at freqs, evaluated on x grid."""
    spec = np.zeros_like(x)
    for f in freqs:
        spec += np.exp(-((x - f)**2) / (2 * width**2))
    return spec


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

def plot_spectrum(data, xmax=2500, width=20.0, outdir=None, noshow=False):
    """Two-panel figure: stick spectrum + Gaussian-broadened."""
    freqs = data['freqs']
    real_f = filter_real_freqs(freqs, threshold=10.0)
    imag_f = get_imaginary_freqs(freqs, threshold=10.0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(f'Vibrational spectrum — {data["name"]}  '
                 f'(N={len(data.get("symbols", []))}, {len(real_f)} real modes'
                 + (f', {len(imag_f)} imaginary' if len(imag_f) else '') + ')',
                 fontsize=13)

    x = np.arange(0, xmax + 1, 0.5)

    # Stick panel
    ax1.vlines(real_f, 0, 1, colors='#2563eb', linewidth=0.6, alpha=0.75, label='real')
    if len(imag_f):
        ax1.vlines(-imag_f, 0, 0.5, colors='#dc2626', linewidth=0.6, alpha=0.6, label='imaginary')
    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)

    # Broadened panel
    spec = broaden(real_f, x, width)
    if spec.max() > 0:
        spec /= spec.max()
    ax2.plot(x, spec, color='#1d4ed8', lw=1.4, label=f'Gaussian (σ={width:.0f} cm⁻¹)')
    ax2.fill_between(x, 0, spec, color='#3b82f6', alpha=0.15)
    ax2.set_xlim(0, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, alpha=0.25, ls=':')

    plt.tight_layout()
    if outdir:
        out = Path(outdir) / f'{data["name"]}_spectrum.png'
        fig.savefig(str(out), dpi=150)
        print(f"  Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


# ============================================================
# Plotting: Overlay spectra across cases
# ============================================================

def plot_overlay(all_data, xmax=2500, width=20.0, outdir=None, noshow=False):
    """Overlay broadened spectra from all cases."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    fig.suptitle('Vibrational spectra overlay — all cases', fontsize=14)
    x = np.arange(0, xmax + 1, 0.5)
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(all_data)))

    for i, data in enumerate(all_data):
        freqs = data['freqs']
        real_f = filter_real_freqs(freqs, threshold=10.0)
        c = colors[i]
        label = data['name']
        ax1.vlines(real_f, 0, 1, colors=c, linewidth=0.5, alpha=0.6, label=label)
        spec = broaden(real_f, x, width)
        if spec.max() > 0:
            spec /= spec.max()
        ax2.plot(x, spec, color=c, lw=1.2, alpha=0.8, label=label)

    ax1.set_ylim(0, 1.1)
    ax1.set_ylabel('Stick (arb.)')
    ax1.legend(fontsize=8, loc='upper right', ncol=2)
    ax1.grid(True, alpha=0.25, ls=':')
    ax1.tick_params(labelbottom=False)

    ax2.set_xlim(0, xmax)
    ax2.set_xlabel('Frequency (cm⁻¹)')
    ax2.set_ylabel('Intensity (normalized)')
    ax2.legend(fontsize=8, loc='upper right', ncol=2)
    ax2.grid(True, alpha=0.25, ls=':')

    plt.tight_layout()
    if outdir:
        out = Path(outdir) / 'overlay_spectra.png'
        fig.savefig(str(out), dpi=150)
        print(f"Saved: {out}")
    if not noshow:
        plt.show()
    plt.close(fig)


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
        spec = broaden(real_f, x, width=25.0)
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
    """Print summary table of all cases."""
    print(f"\n{'='*80}")
    print(f"  SUMMARY: {len(all_data)} cases")
    print(f"{'='*80}")
    print(f"  {'Case':<20s} {'N':>4s} {'3N':>5s} {'n_real':>7s} {'n_imag':>7s} {'E(Ha)':>14s} {'status':>8s}")
    print(f"  {'-'*20} {'-'*4} {'-'*5} {'-'*7} {'-'*7} {'-'*14} {'-'*8}")
    for d in all_data:
        N = d['hessian'].shape[0]
        real_f = filter_real_freqs(d['freqs'], threshold=10.0)
        imag_f = get_imaginary_freqs(d['freqs'], threshold=10.0)
        E = d.get('status', {}).get('energy_Ha', float('nan'))
        status = d.get('status', {}).get('status', '?')
        print(f"  {d['name']:<20s} {N:4d} {3*N:5d} {len(real_f):7d} {len(imag_f):7d} {E:14.8f} {status:>8s}")
    print()


# ============================================================
# Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(description='Plot and analyze PySCF vibration results')
    ap.add_argument('results_dir', help='Directory with case subdirectories (each containing frequencies_cm1.npy, hessian.npy, etc.)')
    ap.add_argument('--cases', nargs='+', default=None, help='Specific case names to process (default: all)')
    ap.add_argument('--outdir', default=None, help='Output directory for PNG files (default: <results_dir>/plots)')
    ap.add_argument('--xmax', type=float, default=2500, help='Max frequency axis (cm^-1)')
    ap.add_argument('--width', type=float, default=20.0, help='Gaussian broadening width (cm^-1)')
    ap.add_argument('--overlay', action='store_true', help='Overlay spectra from all cases')
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
        case_dirs = [c for c in case_dirs if os.path.basename(c) in wanted]
        missing = wanted - {os.path.basename(c) for c in case_dirs}
        if missing:
            print(f"WARNING: cases not found: {missing}")

    if not case_dirs:
        print(f"No cases found in {results_dir}")
        sys.exit(1)

    all_data = []
    for cd in case_dirs:
        try:
            data = load_case(cd)
            all_data.append(data)
            print(f"  Loaded: {data['name']}  (N={data['hessian'].shape[0]})")
        except Exception as e:
            print(f"  FAILED to load {cd}: {e}")

    if not all_data:
        print("No data loaded successfully")
        sys.exit(1)

    print_summary(all_data)

    # Per-case plots
    for data in all_data:
        if not args.no_hessian:
            print(f"  Plotting Hessian: {data['name']}")
            plot_hessian(data, outdir=outdir, noshow=args.noshow)
            plot_hessian_offdiag(data, outdir=outdir, noshow=args.noshow)
        if not args.no_spectrum:
            print(f"  Plotting spectrum: {data['name']}")
            plot_spectrum(data, xmax=args.xmax, width=args.width, outdir=outdir, noshow=args.noshow)

    # Overlay
    if args.overlay:
        print("  Plotting overlay spectra...")
        plot_overlay(all_data, xmax=args.xmax, width=args.width, outdir=outdir, noshow=args.noshow)

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
