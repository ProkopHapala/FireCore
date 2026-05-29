#!/usr/bin/env python3
"""
plot_vib_spectra.py - Overlay vibrational spectra from multiple methods for one molecule.

Usage:
    python plot_vib_spectra.py adamantane
    python plot_vib_spectra.py sila_adamantane
    python plot_vib_spectra.py adamantane --xmax 3500 --width 20
    python plot_vib_spectra.py adamantane --noshow  # save only, no GUI
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from vib_utils import freq_npy_path, filter_real_freqs

COLORS = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown']

# Method tags that may exist per molecule (in preferred display order)
METHOD_TAGS = {
    'adamantane':      ['dftb_mio-1-1', 'dftb_3ob-3-1', 'cp2k_pbe_szv_molopt_sr_gth', 'gpaw_lcao_dzp_pbe', 'psi4_b3lyp_cc-pvdz', 'pyscf_hf_sto-3g', 'pyscf_b3lyp_sto-3g'],
    'sila_adamantane': ['dftb_matsci-0-3', 'dftb_pbc-0-3', 'cp2k_pbe_szv_molopt_sr_gth', 'gpaw_lcao_dzp_pbe', 'psi4_b3lyp_cc-pvdz', 'pyscf_hf_sto-3g', 'pyscf_b3lyp_sto-3g'],
}

METHOD_LABELS = {
    'dftb_mio-1-1':              'DFTB+ mio-1-1',
    'dftb_3ob-3-1':              'DFTB+ 3ob-3-1',
    'dftb_matsci-0-3':           'DFTB+ matsci-0-3',
    'dftb_pbc-0-3':              'DFTB+ pbc-0-3',
    'cp2k_pbe_szv_molopt_sr_gth': 'CP2K PBE/SZV-MOLOPT',
    'gpaw_lcao_dzp_pbe':         'GPAW PBE/dzp',
    'psi4_b3lyp_cc-pvdz':        'Psi4 B3LYP/cc-pVDZ',
    'pyscf_hf_sto-3g':           'HF / STO-3G',
    'pyscf_b3lyp_sto-3g':        'B3LYP / STO-3G',
}


def broaden(freqs, x, width):
    """Sum of Gaussians centered at freqs, evaluated on x grid."""
    spec = np.zeros_like(x)
    for f in freqs:
        spec += np.exp(-((x - f)**2) / (2 * width**2))
    return spec


def load_all(mol_name, workdir='.'):
    """Load all available frequency arrays for this molecule. Returns dict {tag: freqs}."""
    data = {}
    for tag in METHOD_TAGS.get(mol_name, []):
        p = freq_npy_path(mol_name, tag, workdir)
        if p.exists():
            data[tag] = np.load(str(p))
            print(f"  Loaded {tag}: {len(data[tag])} modes from {p}")
        else:
            print(f"  Missing: {p}")
    return data


def plot_overlay(mol_name, data, xmax=3500, width=20, workdir='.', noshow=False):
    """Two-panel figure: top=stick, bottom=broadened, one curve per method."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    fig.suptitle(f'Vibrational spectra — {mol_name}', fontsize=14)
    x = np.arange(0, xmax + 1, 0.5)

    for i, (tag, freqs) in enumerate(data.items()):
        color = COLORS[i % len(COLORS)]
        label = METHOD_LABELS.get(tag, tag)
        rf = filter_real_freqs(freqs, threshold=10.0)

        # Stick panel
        ax1.vlines(rf, 0, 1, color=color, linewidth=0.8, alpha=0.7, label=label)

        # Broadened panel
        spec = broaden(rf, x, width)
        if spec.max() > 0:
            spec /= spec.max()
        ax2.plot(x, spec, color=color, linewidth=1.5, label=label, alpha=0.85)

    ax1.set_ylabel('Stick intensity (arb.)', fontsize=11)
    ax1.set_ylim(0, 1.3)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, alpha=0.25)
    ax1.tick_params(labelbottom=False)

    ax2.set_xlabel('Frequency (cm$^{-1}$)', fontsize=12)
    ax2.set_ylabel('Intensity (normalized)', fontsize=11)
    ax2.set_xlim(0, xmax)
    ax2.legend(fontsize=9, loc='upper left')
    ax2.grid(True, alpha=0.25)

    plt.tight_layout()
    out = Path(workdir) / f'{mol_name}_vib_overlay.png'
    fig.savefig(str(out), dpi=150)
    print(f"Saved: {out}")
    if not noshow:
        plt.show()
    return fig


def main():
    parser = argparse.ArgumentParser(description='Plot vibrational spectra overlay (multiple methods, one molecule)')
    parser.add_argument('molecule', choices=list(METHOD_TAGS.keys()))
    parser.add_argument('--workdir', default='.', help='Directory with cached frequency files')
    parser.add_argument('--xmax',   type=float, default=3500, help='Max frequency axis (cm^-1)')
    parser.add_argument('--width',  type=float, default=20,   help='Gaussian broadening width (cm^-1)')
    parser.add_argument('--noshow', action='store_true', help='Save figure only, do not open GUI')
    args = parser.parse_args()

    print(f"\nLoading frequencies for: {args.molecule}")
    data = load_all(args.molecule, workdir=args.workdir)
    assert data, f"No frequency data found for {args.molecule}. Run run_vib_spectra.py first."

    plot_overlay(args.molecule, data, xmax=args.xmax, width=args.width,
                 workdir=args.workdir, noshow=args.noshow)


if __name__ == '__main__':
    main()
