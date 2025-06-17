#!/usr/bin/env python3
"""faf_func.py – ultra-minimal functional re-implementation of
FoldedAtomicFunctions.main_example with no classes.

The script:
1. builds a 2-D grid
2. adds two atoms with periodic replicas
3. sums Coulomb + Morse potentials
4. fits the raw potential to a plane-wave×exp(z) basis
5. shows three concise plots (raw slice, fit vs ref, singular values)

All helpers are one-liners or short blocks, re-using *utils*, *basis* and
*optimize* modules where possible.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from plot_utils import (
    plot1D,
    plot_SV,
    imshow_grid,
    plot2Dapprox,
    plot2Dbasis,
)
from basis_utils import (
    coulomb2D,
    morse2D,
    cos_exp_basis,
)

# -----------------------------------------------------------------------------
# 0. Tiny helper functions (mostly one-liners)                                  
# -----------------------------------------------------------------------------
mk_grid  = lambda Lx, Lz, dx, z0: (np.arange(0, Lx, dx), z0 + np.arange(0, Lz, dx))
fit_ls   = lambda V, phi: np.linalg.lstsq(phi.T, V.ravel(), rcond=None)[0]
reco     = lambda c, phi, shp: (c @ phi).reshape(shp)
rmse     = lambda A, B: float(np.sqrt(((A - B) ** 2).mean()))

# -----------------------------------------------------------------------------
# 1. Main driver                                                                
# -----------------------------------------------------------------------------

def main():
    # parameters (kept very short names)
    Lx, Lz, dx, z0 = 10.0, 10.0, 0.2, 2.0
    n_img = 3  # ± images in x

    # atoms (two examples, charges zero – coulomb part will be zero)
    ats = [
        dict(x=0.0, y=0.0, z=0.0, q=0.0, r0=1.6, D=0.01, a=1.6),
        dict(x=Lx / 2, y=0.0, z=0.0, q=0.0, r0=2.3, D=0.01, a=1.6),
    ]

    xs, zs = mk_grid(Lx, Lz, dx, z0)
    X, Z = np.meshgrid(xs, zs, indexing="ij")

    Vc = coulomb2D(X, Z, ats, n_img)
    Vm = morse2D(X, Z, ats, n_img)
    V = Vc + Vm  # reference potential within one unit cell

    # extent for imshow helpers
    extent = [0, Lx, z0, z0 + Lz]

    # plot reference potential with atoms
    colors = ["blue", "red"]
    atoms_vis = [ {"x": at["x"], "z": at["z"], "r0": at["r0"], "color": colors[i % len(colors)]} for i, at in enumerate(ats) ]
    imshow_grid(V, extent, "Reference Morse potential", atoms=atoms_vis, fname="ref_potential.png")

    # basis + fit (plane waves × exp)
    nx_h, nz_f = 3, 4
    phi = cos_exp_basis(X, Z, nx=nx_h, nz=nz_f)
    coeffs = fit_ls(V, phi)
    V_fit = reco(coeffs, phi, X.shape)

    # compare fit vs ref
    plot2Dapprox(V, V_fit, extent, fname="fit_vs_ref.png")

    # singular values of phi (library only)
    plot_SV(np.linalg.svd(phi, compute_uv=False), K_opt=0)

    # visualise basis rows (with coefficients)
    labels = [f"cos({k})*e^-{j}z" for k in range(nx_h + 1) for j in range(1, nz_f + 1)]
    plot2Dbasis(phi, X.shape, extent, coeffs=coeffs, labels=labels, fname="basis_set.png")

    # 1-D central slice quick check
    ix = X.shape[0] // 2
    plot1D(zs, np.vstack((V[ix], V_fit[ix])), "Central slice ref vs fit")

    print("RMSE (eV):", rmse(V, V_fit))
    return dict(grid=(xs, zs), V=V, V_fit=V_fit, coeffs=coeffs)

if __name__ == "__main__":
    main()
