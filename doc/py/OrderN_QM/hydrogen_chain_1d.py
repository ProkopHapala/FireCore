import argparse
import numpy as np
import matplotlib.pyplot as plt

from OrderN import (
    find_neighbor_pairs,
    build_matrix_from_pairs,
    plot_matrix,
    canonical_reference,
    plot_spectrum_and_map,
)


def build_dimer_chain(n_cells=10, cell_a=2.5, dimer_distance=0.74):
    """Return positions (N,3) and cell indices for a 1D chain of dimers."""
    atoms = []
    cell_ids = []
    for ic in range(n_cells):
        x0 = ic * cell_a
        atoms.append([x0, 0.0, 0.0])
        atoms.append([x0 + dimer_distance, 0.0, 0.0])
        cell_ids.extend([ic, ic])
    pos = np.array(atoms)
    cell_ids = np.array(cell_ids, dtype=int)
    print(f"#DEBUG build_dimer_chain n_atoms={len(atoms)} cell_a={cell_a} d={dimer_distance}")
    return pos, cell_ids


def build_hamiltonian(pos, cell_ids, onsite=0.0, decay=1.0, t0=-1.0, wrap_ends=False, conf_amp=0.0):
    # include both intra-dimer and inter-cell neighbors
    rmax = 1.05 * max(np.diff(np.sort(np.unique(pos[:,0]))).max(), 1e-6)
    pairs = find_neighbor_pairs(pos, cell_ids=cell_ids, tol=1.05, max_distance=rmax)
    if wrap_ends and len(pos) > 2:
        inter_r = pos[2,0] - pos[1,0] if pos.shape[0] >= 3 else rmax
        pairs.append((len(pos)-1, 0, inter_r))
    H = build_matrix_from_pairs(len(pos), pairs, diag_value=onsite, prefactor=t0, decay=decay)
    if conf_amp != 0.0:
        xs = pos[:,0]
        x0, x1 = xs.min(), xs.max()
        L = max(x1 - x0, 1e-9)
        xc = 0.5*(x0 + x1)
        shape = 0.5 * (1.0 - np.cos(np.pi * (xs - xc) / (0.5*L)))  # 0 at center, 1 at ends
        np.fill_diagonal(H, onsite + conf_amp * shape)
    return H


def build_overlap(pos, cell_ids, diag=1.0, decay=1.0, s0=1.0, wrap_ends=False):
    rmax = 1.05 * max(np.diff(np.sort(np.unique(pos[:,0]))).max(), 1e-6)
    pairs = find_neighbor_pairs(pos, cell_ids=cell_ids, tol=1.05, max_distance=rmax)
    if wrap_ends and len(pos) > 2:
        inter_r = pos[2,0] - pos[1,0] if pos.shape[0] >= 3 else rmax
        pairs.append((len(pos)-1, 0, inter_r))
    S = build_matrix_from_pairs(len(pos), pairs, diag_value=diag, prefactor=s0, decay=decay)
    return S


def plot_hamiltonian(H, title="1D dimer chain Hamiltonian"):
    plot_matrix(H, title=title, symmetric=True)


def plot_overlap(S, title="1D overlap matrix"):
    plot_matrix(S, title=title, symmetric=True)



def plot_geometry(pos, cell_a, title="1D geometry"):
    fig, ax = plt.subplots(figsize=(7, 2.5))
    ax.scatter(pos[:, 0], np.zeros_like(pos[:, 0]), c="k", s=30, label="H")
    # Draw bonds within cell and to next cell (visual only)
    for i in range(0, len(pos), 2):
        ax.plot(pos[i : i + 2, 0], [0, 0], "-", c="tab:blue", lw=2, label="dimer" if i == 0 else None)
    # guide lines for cells
    xmax = pos[:, 0].max()
    for x in np.arange(0, xmax + 1e-6, cell_a):
        ax.axvline(x, color="gray", ls="--", lw=0.5)
    ax.set_xlabel("x [Å]")
    ax.set_yticks([])
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="1D hydrogen dimer chain tight-binding reference")
    ap.add_argument("--n_cells",  type=int, default=6)
    ap.add_argument("--cell_a",   type=float, default=2.0)
    ap.add_argument("--dimer_distance", type=float, default=0.74)
    # Defaults aligned to Harrison-like scale in eV/Å
    ap.add_argument("--onsite",   type=float, default=-13.6)   # H 1s energy in eV
    ap.add_argument("--decay",    type=float, default=0.60)   # decay length in Å
    ap.add_argument("--t0",       type=float, default=-24.5)  # ~-18 eV at 0.74 Å with decay 0.60
    ap.add_argument("--s0",              type=float, default=1.0)    # overlap scale
    ap.add_argument("--TS_decay_factor", type=float, default=1.2)
    ap.add_argument("--conf_amp",        type=float, default=1.0, help="cosine confinement amplitude added to onsite")
    ap.add_argument("--wrap_ends",    type=int, default=1, help="periodic end-to-end coupling")
    ap.add_argument("--plot_geom",    type=int, default=1, help="show geometry scatter")
    ap.add_argument("--plot_mats",    type=int, default=1, help="show H and S matrices")
    ap.add_argument("--plot_map",     type=int, default=1, help="show MO coefficient map")
    ap.add_argument("--plot_density", type=int, default=1, help="show Mulliken density")
    ap.add_argument("--no_show",      type=int, default=0, help="skip plt.show (for batch)")
    ap.add_argument("--verbosity",    type=int, default=1, help="printing verbosity; states printed if >2")
    args = ap.parse_args()

    pos, cell_ids = build_dimer_chain(n_cells=args.n_cells, cell_a=args.cell_a, dimer_distance=args.dimer_distance)
    H = build_hamiltonian(pos, cell_ids, onsite=args.onsite, decay=args.decay, t0=args.t0, wrap_ends=bool(args.wrap_ends), conf_amp=args.conf_amp)
    S = build_overlap(pos, cell_ids, diag=1.0, decay=args.decay*args.TS_decay_factor, s0=args.s0, wrap_ends=bool(args.wrap_ends))
    print(f"#DEBUG |H| max={np.abs(H).max():.3f} min={np.abs(H).min():.3f}")
    print(f"#DEBUG |S| max={np.abs(S).max():.3f} min={np.abs(S).min():.3f}")

    if args.plot_geom:
        plot_geometry(pos, cell_a=args.cell_a, title="1D hydrogen dimer chain geometry")

    if args.plot_mats:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        vmax_h = np.abs(H).max()
        im0 = axes[0].imshow(H, origin="lower", cmap="bwr", vmin=-vmax_h, vmax=vmax_h)
        axes[0].set_title("Hamiltonian")
        axes[0].set_xlabel("j"); axes[0].set_ylabel("i")
        vmax_s = np.abs(S).max()
        im1 = axes[1].imshow(S, origin="lower", cmap="bwr", vmin=-vmax_s, vmax=vmax_s)
        axes[1].set_title("Overlap")
        axes[1].set_xlabel("j"); axes[1].set_ylabel("i")
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
        plt.tight_layout()

    # canonical O(N^3) reference
    e, C, rho = canonical_reference(H, S, n_occ=len(pos)//2)
    gap = e[len(pos)//2] - e[len(pos)//2 - 1]
    homo = e[len(pos)//2 - 1]
    lumo = e[len(pos)//2]
    valence_width = homo - e[0]
    print(f"#DEBUG eigen min/max {e.min():.3f}/{e.max():.3f} gap={gap:.4f} homo={homo:.4f} lumo={lumo:.4f} valence_width={valence_width:.4f}")
    print(f"#INFO  BandGap={gap:.4f}  HOMO={homo:.4f}  LUMO={lumo:.4f}  ValenceWidth={valence_width:.4f}")
    # report representative hoppings/overlaps
    if len(pos) >= 4:
        t_intra = H[0, 1]
        s_intra = S[0, 1]
        t_inter = H[1, 2]
        s_inter = S[1, 2]
        print(f"#INFO  T on/off ( {t_intra:.4f} , {t_inter:.4f} )")
        print(f"#INFO  S on/off ( {s_intra:.4f} , {s_inter:.4f} )")
    # dump eigenstates (controlled by verbosity)
    if args.verbosity > 2:
        for idx, ei in enumerate(e):
            coeffs = " ".join(f"{c:.3f}" for c in C[:, idx])
            print(f"#STATE {idx:02d} E={ei:.6f}  coeffs: {coeffs}")

    if args.plot_map:
        info = f"gap={gap:.3f}\nHOMO={homo:.3f}\nLUMO={lumo:.3f}\nvalence width={valence_width:.3f}"
        plot_spectrum_and_map(e, C, n_occ=len(pos)//2, positions=pos[:,0], plot_spectrum=False, info_text=info)

    if args.plot_density:
        fig, ax = plt.subplots(figsize=(7, 2.5))
        xs = pos[:,0]
        order = np.argsort(xs)
        ax.plot(xs[order], rho[order], "-k", lw=3, label="density (Mulliken)")
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("ρ")
        ax.set_title("Total density along chain")
        ax.legend()
        plt.tight_layout()

    if not args.no_show:
        plt.show()
