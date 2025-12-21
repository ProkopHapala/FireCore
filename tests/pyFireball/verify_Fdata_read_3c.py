import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add pyBall to path
sys.path.append(os.path.abspath('../../'))
import pyBall.FireCore as fc
from pyBall.FireballOCL.FdataParser import FdataParser


def main():
    ap = argparse.ArgumentParser(description="Compare 3c bcna table as read by Python vs exported from Fortran.")
    ap.add_argument("--root", default="bcna", help="table root (bcna/den3)")
    ap.add_argument("--interaction", type=int, default=1, help="interaction id (1 for bcna)")
    ap.add_argument("--isorp", type=int, default=0, help="spin channel")
    ap.add_argument("--nz1", type=int, default=1)
    ap.add_argument("--nz2", type=int, default=1)
    ap.add_argument("--nz3", type=int, default=1)
    ap.add_argument("--itheta", type=int, default=1, help="theta component 1..5")
    ap.add_argument("--iME", type=int, default=1, help="matrix element index (1-based)")
    ap.add_argument("--i_nz", type=int, default=0, help="nonzero slot in Python-loaded table (0-based)")
    ap.add_argument("--fdata", default="./Fdata", help="Fdata directory")
    ap.add_argument("--out", default="compare_bcna_tables.png", help="output image filename")
    args = ap.parse_args()

    # Ensure Fortran tables are initialized (sets icon3c, bcna_* arrays)
    # Minimal dummy system: 1 atom type 1 at origin
    fc.initialize(atomType=np.array([1], dtype=np.int32),
                  atomPos=np.zeros((3,1), dtype=np.float64),
                  verbosity=0)

    parser = FdataParser(args.fdata)
    _, d3c = parser.load_species_data([args.nz1, args.nz2, args.nz3])

    fig, axes = plt.subplots(5, 3, figsize=(12, 16))
    fig.suptitle(f"bcna tables nz=({args.nz1},{args.nz2},{args.nz3}) isorp={args.isorp} iME={args.iME}")

    for itheta in range(1, 6):
        key = (args.root, itheta, args.isorp, args.nz1, args.nz2, args.nz3)
        rec = d3c.get(key)
        if rec is None:
            raise RuntimeError(f"Missing record {key} in Fdata")
        data_py = rec["data"]
        if args.i_nz < 0 or args.i_nz >= data_py.shape[2]:
            raise RuntimeError(f"i_nz out of range for data shape {data_py.shape}")
        slice_py = data_py[:, :, args.i_nz]
        maxsize = slice_py.size

        arr_f, hx_f, hy_f, status = fc.export_bcna_table(args.interaction, args.isorp, args.nz1, args.nz2, args.nz3, itheta, args.iME, maxsize)
        if arr_f.shape != slice_py.shape:
            raise RuntimeError(f"Shape mismatch Fortran {arr_f.shape} vs Python {slice_py.shape} for itheta={itheta}")

        vlim = max(np.max(np.abs(slice_py)), np.max(np.abs(arr_f)))
        diff = arr_f - slice_py
        dlim = np.max(np.abs(diff))

        r = itheta - 1
        ax_f = axes[r, 0]; ax_p = axes[r, 1]; ax_d = axes[r, 2]
        imf = ax_f.imshow(arr_f, origin="lower", cmap="coolwarm", vmin=-vlim, vmax=vlim)
        imp = ax_p.imshow(slice_py, origin="lower", cmap="coolwarm", vmin=-vlim, vmax=vlim)
        imd = ax_d.imshow(diff, origin="lower", cmap="bwr", vmin=-dlim, vmax=dlim)
        ax_f.set_title(f"Fortran bcna_0{itheta}")
        ax_p.set_title(f"Python bcna_0{itheta}")
        ax_d.set_title(f"Diff (F-P) bcna_0{itheta}")
        ax_f.set_xlabel("x"); ax_p.set_xlabel("x"); ax_d.set_xlabel("x")
        ax_f.set_ylabel("y")
        fig.colorbar(imf, ax=ax_f, fraction=0.046, pad=0.04)
        fig.colorbar(imp, ax=ax_p, fraction=0.046, pad=0.04)
        fig.colorbar(imd, ax=ax_d, fraction=0.046, pad=0.04)

        # quick stats: Fortran min/max, Python min/max, max abs diff
        print(f"[itheta={itheta}] F[min,max]=({arr_f.min(): .9e},{arr_f.max(): .9e}) "
              f"P[min,max]=({slice_py.min(): .9e},{slice_py.max(): .9e}) "
              f"max|F-P|={dlim: .9e}")

    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.savefig(args.out, dpi=200)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
