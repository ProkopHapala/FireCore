import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from pyBall.FireballOCL.FdataParser import FdataParser, read_pp

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


def plot_na(parser: FdataParser, nz: int, out: str):
    """
    Plot neutral-atom potentials basis/{nz}_*.na* (may have multiple shells).
    """
    pattern = os.path.join(parser.fdata_dir, "basis", f"{nz:03d}_*.na*")
    paths = sorted(glob.glob(pattern))
    if not paths:
        print(f"[DEBUG] NA files not found for nz={nz} pattern={pattern}")
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    for path in paths:
        r, v = read_pp(path)
        ax.plot(r, v, label=os.path.basename(path))
    ax.set_xlabel("r (Angstrom)")
    ax.set_ylabel("V_na(r)")
    ax.legend(fontsize=8)
    ax.set_title(f"Neutral atom potentials nz={nz}")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"Saved NA plot to {out}")


# --- Wavefunction comparison and polynomial fits ---

def _load_wf(parser: FdataParser, nz: int, shell: int):
    wf_files = parser.find_wf(nz)
    print(f"[DEBUG] wavefunction files for nz={nz}: {wf_files}")
    if not wf_files or shell < 1 or shell > len(wf_files):
        raise RuntimeError(f"Missing wf shell={shell} for nz={nz}")
    rec = parser.read_wf(wf_files[shell - 1])
    r = np.linspace(0.0, rec["rcmax"], rec["mesh"])
    return r, rec["data"], rec["l"], rec["rcmax"], os.path.basename(wf_files[shell - 1])


def _fit_poly_best(r, y, degs):
    best = None
    best_sse = None
    for d in degs:
        coefs = np.polynomial.polynomial.polyfit(r, y, d)
        y_fit = np.polynomial.polynomial.Polynomial(coefs)(r)
        sse = float(np.sum((y - y_fit) ** 2))
        if (best_sse is None) or (sse < best_sse):
            best_sse = sse
            best = (d, coefs, y_fit)
    return best_sse, best


def _wf_fit_summary(parser: FdataParser, elements, degs, shells=(1, 2)):
    summary = []
    for nz in elements:
        for shell in shells:
            try:
                r, psi, l, rc, name = _load_wf(parser, nz, shell=shell)
                sse, best = _fit_poly_best(r, psi, degs)
                d, coefs, fit = best
                summary.append((nz, l, shell, d, sse, coefs, r, psi, fit))
            except Exception as e:
                print(f"[DEBUG] fit skip shell={shell} Z{nz}: {e}")
    return summary


def _print_fit_summary(summary):
    print("Fit summary (nz, l, shell, degree, SSE, coefs):")
    for nz, l, shell, d, sse, coefs, _, _, _ in summary:
        print(f"Z{nz} l={l} shell={shell} deg{d} SSE={sse:.3e} coefs={coefs}")


def plot_wf_compare_by_L(parser: FdataParser, elements=(1, 6, 7, 8), degs=range(4, 9)):
    """
    Overlay s-shell (panel A) and p-shell (panel B) across elements (data only).
    """
    colors = plt.cm.tab10.colors
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    ax_s, ax_p = axs
    for i, nz in enumerate(elements):
        col = colors[i % len(colors)]
        try:
            r_s, psi_s, l_s, rc_s, name_s = _load_wf(parser, nz, shell=1)
            ax_s.plot(r_s, psi_s, color=col, label=f"Z{nz} s")
        except Exception as e:
            print(f"[DEBUG] skip s for Z{nz}: {e}")
        try:
            r_p, psi_p, l_p, rc_p, name_p = _load_wf(parser, nz, shell=2)
            ax_p.plot(r_p, psi_p, color=col, label=f"Z{nz} p")
        except Exception as e:
            print(f"[DEBUG] skip p for Z{nz}: {e}")

    ax_s.set_title("s (shell1)")
    ax_p.set_title("p (shell2)")
    for ax in axs:
        ax.set_xlabel("r (Å)")
        ax.set_ylabel("ψ(r)")
        ax.grid(True)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("wf_compare_sp.png", dpi=200)
    print("Saved wf_compare_sp.png")
    return None


def plot_wf_fits(parser: FdataParser, elements=(1, 6, 7, 8), degs=range(4, 9)):
    """
    Fit and plot wavefunctions (s and p) across elements, showing data + fit; prints summary.
    """
    summary = []
    colors = plt.cm.tab10.colors
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    ax_s, ax_p = axs
    for i, nz in enumerate(elements):
        col = colors[i % len(colors)]
        try:
            r_s, psi_s, l_s, rc_s, name_s = _load_wf(parser, nz, shell=1)
            sse_s, best_s = _fit_poly_best(r_s, psi_s, degs)
            d_s, coefs_s, fit_s = best_s
            ax_s.plot(r_s, psi_s, color=col, label=f"Z{nz} s data")
            ax_s.plot(r_s, fit_s, "--", color=col, lw=0.8, label=f"Z{nz} s fit d{d_s}")
            summary.append((nz, "s", 1, d_s, sse_s, coefs_s, r_s, psi_s, fit_s))
        except Exception as e:
            print(f"[DEBUG] fit skip s for Z{nz}: {e}")
        try:
            r_p, psi_p, l_p, rc_p, name_p = _load_wf(parser, nz, shell=2)
            sse_p, best_p = _fit_poly_best(r_p, psi_p, degs)
            d_p, coefs_p, fit_p = best_p
            ax_p.plot(r_p, psi_p, color=col, label=f"Z{nz} p data")
            ax_p.plot(r_p, fit_p, "--", color=col, lw=0.8, label=f"Z{nz} p fit d{d_p}")
            summary.append((nz, "p", 2, d_p, sse_p, coefs_p, r_p, psi_p, fit_p))
        except Exception as e:
            print(f"[DEBUG] fit skip p for Z{nz}: {e}")
    ax_s.set_title("s fits")
    ax_p.set_title("p fits")
    for ax in axs:
        ax.set_xlabel("r (Å)")
        ax.set_ylabel("ψ(r)")
        ax.grid(True)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("wf_fits.png", dpi=200)
    print("Saved wf_fits.png")
    _print_fit_summary(summary)
    return summary


def plot_wf_compare_by_elem(parser: FdataParser, elements=(1, 6, 7, 8), degs=range(4, 9)):
    """
    Stacked per-element plots (s and p where available), data only.
    """
    colors = plt.cm.tab10.colors
    fig, axes = plt.subplots(len(elements), 1, figsize=(8, 2.5 * len(elements)), sharex=False)
    if len(elements) == 1:
        axes = [axes]
    for ax, nz in zip(axes, elements):
        col = colors[elements.index(nz) % len(colors)]
        has_any = False
        try:
            r_s, psi_s, l_s, rc_s, name_s = _load_wf(parser, nz, shell=1)
            ax.plot(r_s, psi_s, color=col, label="s")
            has_any = True
        except Exception as e:
            print(f"[DEBUG] stack skip s for Z{nz}: {e}")
        try:
            r_p, psi_p, l_p, rc_p, name_p = _load_wf(parser, nz, shell=2)
            ax.plot(r_p, psi_p, color="k", label="p", alpha=0.7)
            has_any = True
        except Exception as e:
            print(f"[DEBUG] stack skip p for Z{nz}: {e}")
        ax.set_title(f"Z{nz}" if has_any else f"Z{nz} (no wf)")
        ax.set_xlabel("r (Å)")
        ax.set_ylabel("ψ(r)")
        ax.grid(True)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("wf_per_element.png", dpi=200)
    print("Saved wf_per_element.png")
    return None


def wf_fit_summary(parser: FdataParser, elements=(1, 6, 7, 8), degs=range(4, 9)):
    """
    Compute fits and print summary (no plots).
    """
    summary = _wf_fit_summary(parser, elements, degs, shells=(1, 2))
    _print_fit_summary(summary)
    return summary