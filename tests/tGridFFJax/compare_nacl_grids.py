#!/usr/bin/env python3
"""Compare two NaCl GridFF Bspline_PLQd.npy grids per channel.

For each of the three channels (P=Pauli, L=London, Q=Coulomb), plot:
  - z-profile (mean over the xy plane vs z) for both grids on one axis
  - 2D slice (xy heatmap at the equilibrium-z plane) for each grid
  - shape / range / |max| summary

This is the meaningful apples-to-apples comparison between two GridFF
grids that may have different voxel spacings — it bypasses voxel-by-voxel
arithmetic (which is meaningless across different grids) and looks at
physical observables: where each channel has structure along z, and what
the lateral pattern looks like at the equilibrium plane.

Usage
-----
  python compare_nacl_grids.py \\
      --ref  cpp/common_resources/NaCl_1x1_L1/Bspline_PLQd.npy \\
      --new  tests/tGridFFJax/runs/nacl_ionic_smoke/Bspline_PLQd.npy \\
      --ref-cell 4.0 4.0 20.0 --new-cell 4.0 4.0 20.0 \\
      --out tests/tGridFFJax/runs/nacl_ionic_smoke/compare_vs_pointcharge.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CHANNEL_NAMES = ["Pauli", "London", "Coulomb"]


def load_grid(path: Path):
    arr = np.load(path)
    # File layout (from export_firecore_artifacts):
    #   arr shape (Nx, Ny, Nz, 3)  with channels [Pauli, London, Coulomb]
    # Verify and return.
    if arr.ndim != 4 or arr.shape[-1] != 3:
        raise ValueError(f"Unexpected shape {arr.shape} for {path}; "
                         f"expected (Nx, Ny, Nz, 3)")
    return arr  # (Nx, Ny, Nz, 3)


def z_profile(grid: np.ndarray) -> np.ndarray:
    """Mean over xy for each z slice and each channel. Returns (Nz, 3)."""
    return grid.mean(axis=(0, 1))


def xy_slice(grid: np.ndarray, z_index: int) -> np.ndarray:
    """Returns (Nx, Ny, 3) slice at z_index."""
    return grid[:, :, z_index, :]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True, help="Reference Bspline_PLQd.npy")
    p.add_argument("--new", required=True, help="Candidate Bspline_PLQd.npy")
    p.add_argument("--ref-cell", nargs=3, type=float, default=(4.0, 4.0, 20.0))
    p.add_argument("--new-cell", nargs=3, type=float, default=(4.0, 4.0, 20.0))
    p.add_argument("--z-min", type=float, default=2.0,
                   help="z height (Å) to use for the xy heatmap slices.")
    p.add_argument("--out", required=True, help="Output PNG path")
    return p.parse_args()


def main():
    args = parse_args()
    ref = load_grid(Path(args.ref))
    new = load_grid(Path(args.new))

    print(f"reference: shape={ref.shape}  cell={args.ref_cell}  "
          f"voxel_z={args.ref_cell[2]/ref.shape[2]:.3f} Å")
    print(f"new      : shape={new.shape}  cell={args.new_cell}  "
          f"voxel_z={args.new_cell[2]/new.shape[2]:.3f} Å")
    print()
    print("Per-channel range and |max|:")
    print(f"{'channel':10s}  {'ref_range':>20s}  {'new_range':>20s}  {'ref_|max|':>12s}  {'new_|max|':>12s}")
    for i, name in enumerate(CHANNEL_NAMES):
        rrng = f"[{ref[..., i].min():+.3g}, {ref[..., i].max():+.3g}]"
        nrng = f"[{new[..., i].min():+.3g}, {new[..., i].max():+.3g}]"
        print(f"{name:10s}  {rrng:>20s}  {nrng:>20s}  "
              f"{np.abs(ref[..., i]).max():>12.4g}  {np.abs(new[..., i]).max():>12.4g}")

    z_ref = np.linspace(0.0, args.ref_cell[2], ref.shape[2], endpoint=False)
    z_new = np.linspace(0.0, args.new_cell[2], new.shape[2], endpoint=False)
    iz_ref = int(np.argmin(np.abs(z_ref - args.z_min)))
    iz_new = int(np.argmin(np.abs(z_new - args.z_min)))

    prof_ref = z_profile(ref)
    prof_new = z_profile(new)

    fig, axes = plt.subplots(3, 3, figsize=(13.5, 11.0), constrained_layout=True)
    for i, name in enumerate(CHANNEL_NAMES):
        # Column 0: z-profile (mean over xy)
        ax = axes[i, 0]
        ax.plot(z_ref, prof_ref[:, i], "k-", lw=1.5, label=f"reference (n_z={ref.shape[2]})")
        ax.plot(z_new, prof_new[:, i], "C0--", lw=1.5,
                label=f"new (n_z={new.shape[2]})")
        ax.set_xlabel("z [Å]")
        ax.set_ylabel(f"<{name}>_xy")
        ax.set_title(f"{name} z-profile (mean over xy)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
        ax.axvline(args.z_min, color="r", lw=0.5, ls=":")

        # Columns 1+2: xy heatmap at z=z_min for each grid
        sl_ref = xy_slice(ref, iz_ref)[..., i].T
        sl_new = xy_slice(new, iz_new)[..., i].T
        vmin = min(sl_ref.min(), sl_new.min())
        vmax = max(sl_ref.max(), sl_new.max())
        for ax2, sl, label in [(axes[i, 1], sl_ref, f"ref @ z={z_ref[iz_ref]:.2f} Å"),
                               (axes[i, 2], sl_new, f"new @ z={z_new[iz_new]:.2f} Å")]:
            im = ax2.imshow(sl, origin="lower", aspect="equal",
                            extent=[0, args.new_cell[0], 0, args.new_cell[1]],
                            vmin=vmin, vmax=vmax, cmap="RdBu_r" if name == "Coulomb" else "viridis")
            ax2.set_title(f"{name}  {label}", fontsize=9)
            ax2.set_xlabel("x [Å]"); ax2.set_ylabel("y [Å]")
            plt.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)

    fig.suptitle(f"NaCl GridFF channels — reference (point-charge) vs new (JAX fit)",
                 fontsize=12)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
