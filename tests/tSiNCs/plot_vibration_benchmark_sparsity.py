#!/usr/bin/env python3
"""Plot log10|H_ij| imshow for vibration benchmark .npz files (no Hessian recomputation)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DIR = REPO / "tests/tSiNCs/fixtures/vibration_benchmarks"


def reconstruct_from_npz_blocks(d: np.lib.npyio.NpzFile) -> np.ndarray | None:
    if "blocks" not in d.files:
        return None
    neigh_idx = d["neigh_idx"]
    neigh_count = d["neigh_count"]
    blocks = d["blocks"]
    natoms, max_neigh = neigh_idx.shape
    dim = natoms * 3
    H = np.zeros((dim, dim), dtype=np.float64)
    for p in range(natoms):
        for j in range(int(neigh_count[p])):
            o = int(neigh_idx[p, j])
            if o < 0:
                continue
            H[o * 3:(o + 1) * 3, p * 3:(p + 1) * 3] = blocks[p, j]
    return 0.5 * (H + H.T)


def load_K_projected(d: np.lib.npyio.NpzFile) -> np.ndarray:
    if "H_dense_projected" in d.files:
        return np.asarray(d["H_dense_projected"], dtype=np.float64)
    K = sp.csr_matrix((d["K_csr_data"], d["K_csr_indices"], d["K_csr_indptr"]), shape=tuple(d["K_csr_shape"]))
    return K.toarray()


def plot_logabs_H(H: np.ndarray, out_png: Path, *, title: str, meta: dict, label: str):
    ndof = H.shape[0]
    absH = np.abs(H)
    nnz = int(np.count_nonzero(absH))
    density = nnz / (ndof * ndof)
    floor = max(float(absH[absH > 0].min()) if nnz else 1.0, 1e-30)
    logH = np.full_like(absH, np.nan, dtype=np.float64)
    mask = absH > 0.0
    logH[mask] = np.log10(absH[mask])
    vmin = np.log10(floor)
    vmax = float(np.nanmax(logH))

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    im = ax.imshow(logH, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")
    ax.set_xlabel("j (DOF)")
    ax.set_ylabel("i (DOF)")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\log_{10}|H_{ij}|$  (NaN = exact zero)")
    txt = (f"{label}\nDOF={ndof}  nnz={nnz}  density={100*density:.2f}%\n"
           f"n_shells={meta.get('n_shells', '?')}  eigh={meta.get('eigh_seconds', '?')}s")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    return nnz, density


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=str, default=str(DEFAULT_DIR))
    args = ap.parse_args()
    bench_dir = Path(args.dir)
    if not bench_dir.is_dir():
        raise FileNotFoundError(bench_dir)
    for npz_path in sorted(bench_dir.glob("*.npz")):
        d = np.load(npz_path, allow_pickle=True)
        meta = json.loads(str(d["meta_json"]))
        stem = npz_path.stem
        Hproj = load_K_projected(d)
        out_proj = bench_dir / f"{stem}_Kproj_logabs.png"
        nnz, density = plot_logabs_H(Hproj, out_proj, title=f"{stem}: projected K", meta=meta, label="projected K (solver matrix)")
        print(f"{stem}: projected DOF={Hproj.shape[0]} nnz={nnz} density={100*density:.3f}% -> {out_proj.name}")
        Hblk = reconstruct_from_npz_blocks(d)
        if Hblk is not None:
            out_blk = bench_dir / f"{stem}_Hblocks_logabs.png"
            nnz_b, density_b = plot_logabs_H(Hblk, out_blk, title=f"{stem}: shell blocks (pre-projection)", meta=meta, label="shell-block H (n_shells from .npz)")
            print(f"{stem}: blocks   DOF={Hblk.shape[0]} nnz={nnz_b} density={100*density_b:.3f}% -> {out_blk.name}")


if __name__ == "__main__":
    main()
