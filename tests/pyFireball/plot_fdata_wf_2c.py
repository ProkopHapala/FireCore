import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Utility to visualize Fdata contents:
# - Wavefunctions: one-body radial (basis/*.wf*)
# - 2-center integrals: overlap/kinetic/vna/vnl/etc. (columns are mu-nu matrix elements)
# - vnl: non-local pseudopotential projectors; each column is one projector-channel pair
#   (channels argument controls how many columns to plot for readability)
# - PP cutoff: rc_PP read from info.dat (no channels; just a radius marker)

# Add pyBall to path
sys.path.append(os.path.abspath('../../'))
from pyBall.FireballOCL.FdataParser import FdataParser


def plot_wavefunctions(parser: FdataParser, nz: int, shell: int, out: str):
    """
    Plot radial wavefunction for one species/shell.
    parser : FdataParser with fdata_dir set
    nz     : nuclear charge to select basis files (e.g., 1 for H, 6 for C)
    shell  : 1-based shell index within available wf files for that nz
    out    : output image filename
    """
    # loud debug per guidelines
    wf_files = parser.find_wf(nz)
    print(f"[DEBUG] wavefunction files for nz={nz}: {wf_files}")
    if not wf_files:
        raise RuntimeError(f"No wavefunction files found for nz={nz}")
    if shell < 1 or shell > len(wf_files):
        print(f"[DEBUG] shell index requested={shell}, available={len(wf_files)}")
    fname = wf_files[shell - 1]
    rec = parser.read_wf(fname)
    mesh = rec['mesh']
    rmax = rec['rcmax']
    r = np.linspace(0.0, rmax, mesh)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(r, rec['data'], label=f"nz={nz} shell={shell} l={rec['l']}")
    ax.set_xlabel("r (Angstrom)")
    ax.set_ylabel("psi(r)")
    ax.legend()
    ax.set_title(os.path.basename(fname))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved wavefunction plot to {out}")


def plot_wavefunctions_all(parser: FdataParser, nz: int, out_prefix: str):
    """
    Plot all wavefunction shells for a species.
    nz         : nuclear charge
    out_prefix : filename prefix; suffix _shell{idx}.png is added
    """
    wf_files = parser.find_wf(nz)
    print(f"[DEBUG] wavefunction files for nz={nz}: {wf_files}")
    if not wf_files:
        raise RuntimeError(f"No wavefunction files found for nz={nz}")
    for i, fname in enumerate(wf_files, start=1):
        rec = parser.read_wf(fname)
        mesh = rec['mesh']
        rmax = rec['rcmax']
        r = np.linspace(0.0, rmax, mesh)
        fig, ax = plt.subplots(figsize=(6,4))
        ax.plot(r, rec['data'], label=f"nz={nz} shell={i} l={rec['l']}")
        ax.set_xlabel("r (Angstrom)")
        ax.set_ylabel("psi(r)")
        ax.legend()
        ax.set_title(os.path.basename(fname))
        fig.tight_layout()
        out = f"{out_prefix}_shell{i}.png"
        fig.savefig(out, dpi=200)
        print(f"Saved wavefunction plot to {out}")


def plot_2c(parser: FdataParser, root: str, nz1: int, nz2: int, channels: int, out: str):
    """
    Plot 2-center integral columns for a pair of species.
    parser   : FdataParser
    root     : 2c basename (overlap, kinetic, vna, vnl, vxc, etc.)
    nz1,nz2  : nuclear charges of species 1 and 2
    channels : how many mu-nu columns to draw (limits clutter)
    out      : output image filename
    """
    path = parser.find_2c(root, nz1, nz2)
    if not os.path.exists(path):
        raise RuntimeError(f"2c file not found: {path}")
    rec = parser.read_2c(path)
    z = np.linspace(0.0, rec['zmax'], rec['numz'])
    data = rec['data']
    nchan = min(channels, data.shape[1])
    print(f"[DEBUG] plotting first {nchan} / {data.shape[1]} channels from {path}")
    fig, ax = plt.subplots(figsize=(7,4))
    for i in range(nchan):
        ax.plot(z, data[:, i], label=f"mu-nu {i}")
    ax.set_xlabel("distance z (Angstrom)")
    ax.set_ylabel(f"{root} integral")
    ax.legend(ncol=2, fontsize=8)
    ax.set_title(f"{root} nz1={rec['nucz1']} nz2={rec['nucz2']}")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved 2c plot to {out}")


def plot_2c_list(parser: FdataParser, root: str, specs, channels: int, out_prefix: str):
    """
    Plot 2c integrals for a list of species pairs.
    specs     : list of (nz1, nz2) tuples
    channels  : how many mu-nu columns to draw for each
    out_prefix: filename prefix; suffix _{nz1}_{nz2}.png is added
    """
    for nz1, nz2 in specs:
        out = f"{out_prefix}_{nz1}_{nz2}.png"
        plot_2c(parser, root, nz1, nz2, channels, out)


def plot_2c_compare(parser: FdataParser, root: str, specs, channels: int, out: str):
    """
    Plot multiple species pairs into a single figure for easy comparison.
    specs    : list of (nz1, nz2)
    channels : how many mu-nu columns per pair to overlay
    out      : output image filename
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    for nz1, nz2 in specs:
        path = parser.find_2c(root, nz1, nz2)
        if not os.path.exists(path):
            print(f"[DEBUG] skipping missing {root} {nz1}-{nz2} at {path}")
            continue
        rec = parser.read_2c(path)
        z = np.linspace(0.0, rec['zmax'], rec['numz'])
        data = rec['data']
        nchan = min(channels, data.shape[1])
        for i in range(nchan):
            ax.plot(z, data[:, i], label=f"{nz1}-{nz2} ch{i}")
    ax.set_xlabel("distance z (Angstrom)")
    ax.set_ylabel(f"{root} integral")
    ax.legend(ncol=2, fontsize=8)
    ax.set_title(f"{root} comparison")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved 2c comparison plot to {out}")


def plot_vnl_pp(parser: FdataParser, nz1: int, nz2: int, channels: int, out: str):
    """
    Plot non-local PP integrals (vnl).
    vnl files are projector-channel tables: each column is a mu-nu projector pair.
    channels controls how many columns to draw (helps readability when dozens exist).
    """
    root = "vnl"
    path = parser.find_2c(root, nz1, nz2)
    if not os.path.exists(path):
        print(f"[DEBUG] vnl file not found for nz1={nz1} nz2={nz2}: {path}")
        return
    rec = parser.read_2c(path)
    z = np.linspace(0.0, rec['zmax'], rec['numz'])
    data = rec['data']
    nchan = min(channels, data.shape[1])
    print(f"[DEBUG] plotting first {nchan} / {data.shape[1]} vnl channels from {path}")
    fig, ax = plt.subplots(figsize=(7,4))
    for i in range(nchan):
        ax.plot(z, data[:, i], label=f"PP mu-nu {i}")
    ax.set_xlabel("distance z (Angstrom)")
    ax.set_ylabel("vnl integral")
    ax.legend(ncol=2, fontsize=8)
    ax.set_title(f"vnl nz1={rec['nucz1']} nz2={rec['nucz2']}")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved vnl plot to {out}")


def plot_pp_radii(parser: FdataParser, nz: int, out: str):
    """
    Plot PP cutoff radius rc_PP if available.
    This comes from info.dat and is the radial cutoff used by pseudopotential projectors.
    """
    if not hasattr(parser, "species_info"):
        parser.parse_info()
    info = parser.species_info.get(nz, {})
    rc = info.get("rc_PP", None)
    if rc is None:
        print(f"[DEBUG] rc_PP missing for nz={nz}")
        return
    fig, ax = plt.subplots(figsize=(4,2))
    ax.axvline(rc, color="C1", label=f"rc_PP={rc}")
    ax.set_xlabel("r (Angstrom)")
    ax.set_yticks([])
    ax.legend()
    ax.set_title(f"PP cutoff nz={nz}")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved PP cutoff plot to {out}")


def read_pp(path: str):
    """
    Minimal parser for radial pseudopotential basis/*.pp.
    Returns r (Angstrom) and V(r) arrays using lines that contain two floats.
    """
    if not os.path.exists(path):
        raise RuntimeError(f"PP file not found: {path}")
    r_list, v_list = [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                r_val = float(parts[0].replace('D','E'))
                v_val = float(parts[1].replace('D','E'))
            except ValueError:
                continue
            r_list.append(r_val)
            v_list.append(v_val)
    if not r_list:
        raise RuntimeError(f"Failed to parse any (r,V) pairs from {path}")
    return np.array(r_list), np.array(v_list)


def plot_pp(parser: FdataParser, nz: int, out: str):
    """
    Plot one-body radial pseudopotential from basis/{nz}.pp (no channels).
    """
    path = os.path.join(parser.fdata_dir, "basis", f"{nz:03d}.pp")
    r, v = read_pp(path)
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(r, v, label=f"{nz:03d}.pp")
    ax.set_xlabel("r (Angstrom)")
    ax.set_ylabel("V_pp(r)")
    ax.legend()
    ax.set_title(os.path.basename(path))
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved PP plot to {out}")


def main():
    ap = argparse.ArgumentParser(description="Plot wavefunctions and 2c integrals from Fdata")
    ap.add_argument("--fdata", default="./Fdata", help="Fdata directory")
    ap.add_argument("--nz", type=int, default=1, help="nuclear charge for wavefunction")
    ap.add_argument("--shell", type=int, default=1, help="which wavefunction shell (1-based)")
    ap.add_argument("--wf_out", default="wf.png", help="output image for wavefunction")
    ap.add_argument("--root", default="overlap", help="2c root (overlap, kinetic, vna, vnl, ...)")
    ap.add_argument("--nz1", type=int, default=6, help="nuclear charge 1 for 2c")
    ap.add_argument("--nz2", type=int, default=6, help="nuclear charge 2 for 2c")
    ap.add_argument("--channels", type=int, default=6, help="how many 2c channels to plot")
    ap.add_argument("--c2_out", default="c2.png", help="output image for 2c integrals")
    ap.add_argument("--c2_list", default='1:1,6:6,1:6', help="comma-separated list of nz1:nz2 pairs (e.g., '1:1,6:6') to batch plot 2c")
    ap.add_argument("--vnl_channels", type=int, default=6, help="how many vnl channels to plot")
    ap.add_argument("--vnl_out", default="vnl.png", help="output image for vnl integrals")
    ap.add_argument("--pp_cutoff_out", default="pp_cutoff.png", help="output image for PP cutoff")
    ap.add_argument("--pp_out", default="pp.png", help="output image for one-body PP")
    ap.add_argument("--plot_all_wf", action="store_true", help="plot all wavefunction shells for nz")
    args = ap.parse_args()

    parser = FdataParser(args.fdata)
    # this populates species_info and keeps loud if missing
    parser.parse_info()

    if args.plot_all_wf:
        plot_wavefunctions_all(parser, args.nz, os.path.splitext(args.wf_out)[0])
    else:
        plot_wavefunctions(parser, args.nz, args.shell, args.wf_out)
    if args.c2_list:
        specs = []
        for pair in args.c2_list.split(","):
            if ":" in pair:
                a, b = pair.split(":")
                try:
                    specs.append((int(a), int(b)))
                except ValueError:
                    print(f"[DEBUG] skipping malformed pair '{pair}'")
        if specs:
            plot_2c_list(parser, args.root, specs, args.channels, os.path.splitext(args.c2_out)[0])
    else:
        plot_2c(parser, args.root, args.nz1, args.nz2, args.channels, args.c2_out)
    # Pseudopotential helpers: vnl integrals and rc_PP marker
    plot_vnl_pp(parser, args.nz1, args.nz2, args.vnl_channels, args.vnl_out)
    plot_pp_radii(parser, args.nz1, args.pp_cutoff_out)
    plot_pp(parser, args.nz1, args.pp_out)
    plt.show()


if __name__ == "__main__":
    main()
