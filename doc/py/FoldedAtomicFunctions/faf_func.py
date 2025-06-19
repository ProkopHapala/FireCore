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

from plot_utils import plot1D, plot_SV, imshow_grid, plot_atoms, plot2Dapprox, plot2Dbasis

from basis_utils import cos_exp_basis, mult_atoms, getMorseQ

# -----------------------------------------------------------------------------
# 0. Tiny helper functions (mostly one-liners)                                  
# -----------------------------------------------------------------------------
#mk_grid  = lambda Lx, Lz, dx, z0 : (np.arange(0, Lx, dx), z0 + np.arange(0, Lz, dx))
#fit_ls   = lambda V, phi         :  np.linalg.lstsq(phi.T, V.ravel(), rcond=None)[0]
#reco     = lambda c, phi, shp    : (c @ phi).reshape(shp)
#rmse     = lambda A, B           : float(np.sqrt(((A - B) ** 2).mean()))


def makePotentialXZ( apos, apars, X, Z, apar=(1.4,0.01,0.1,1.6) ):
    apars = apars.copy()
    apars[:,0] =          apars[:,0] + apar[0]
    apars[:,1] = np.sqrt( apars[:,1] * apar[1] )
    apars[:,2] =          apars[:,2] * apar[2]
    apars[:,3] = np.sqrt( apars[:,3] * apar[3] )
    #print( "xs ", xs)
    # print("apos:\n",   apos)
    # print("apars:\n",  apars)
    # print("colors:\n", colors)
    # print("X.shape:", X.shape)
    V = getMorseQ(apos, apars, X, Z)
    V0 = np.average(V[:,-1]); V-=V0; print("V0:", V0)
    return V, extent, X, Z

def scan_Qs( apos, apars, X, Z, Lx=5.0, nbx=8, nbz=5, a0=0.5 ):
    
    phi, labels = cos_exp_basis(X, Z, nbx, nbz, a0=a0, Lx=Lx)

    for i,l in enumerate(labels): print(i,l)

    Qs = [0.0,-0.1,0.1] 
    #Qs = [-0.0,-0.1] 
    #Qs = [-0.1] 
    #Qs = [+0.1] 
    for Qsond in Qs:
        V, extent, X, Z = makePotentialXZ( apos, apars, X, Z, apar=(1.4,0.01,Qsond,1.6))
        
        #ax = imshow_grid(V, extent=extent, title="Reference Morse potential", figsize=(20,6))
        #plot_atoms(apos, colors, sz=100, ax=ax, bEqual=True)

        #plt.show()
       
        coeffs = np.linalg.lstsq(phi.T, V.ravel(), rcond=None)[0]     ;print("coeffs.shape:", coeffs.shape)
        V_fit = (coeffs @ phi).reshape(X.shape)                       ;print("V_fit.shape:", V_fit.shape)

        err = V - V_fit; err[Z<2.0] = np.nan
        plot2Dapprox(V, V_fit, err=err, extent=extent, scErr=0.001, title=f"Fit vs Ref Q={Qsond}")
        plt.savefig(f"fit_vs_ref_Q{Qsond:.2f}.png")

        #plot2Dbasis(phi, X.shape, extent, coeffs=coeffs, labels=labels, nrow=nbz)


def scan_basis_a0( apos, apars, X, Z, Q_list=[-0.1, 0.0, 0.1], a0_list = [0.20,0.40,0.60], nbx=8, nbz=5, Lx=5.0 ):
    # Containers: results[Q]['max'|'avg'] -> list indexed by a0
    results = {Q: {"max": [], "avg": []} for Q in Qs}

    for a0 in a0_list:
        phi, labels = cos_exp_basis(X, Z, nbx, nbz, a0=a0, Lx=Lx)
        for Qval in Q_list:
            V, _extent, _X, _Z = makePotentialXZ(apos, apars, X, Z, apar=(1.4, 0.01, Qval, 1.6))
            coeffs = np.linalg.lstsq(phi.T, V.ravel(), rcond=None)[0]
            V_fit = (coeffs @ phi).reshape(X.shape)

            # error metrics in region z > 2.0
            mask = Z > 2.0
            err_region = np.abs(V - V_fit)[mask]
            results[Qval]["max"].append(np.nanmax(err_region))
            results[Qval]["avg"].append(np.nanmean(err_region))

    # Plot the scan – six curves (max & avg for each charge)
    fig, ax = plt.subplots(figsize=(6, 4))
    markers = {"max": "o", "avg": "s"}
    colors  = {0.0: "k", -0.1: "b", 0.1: "r"}

    for Qval in Q_list:
        ax.plot(a0_list, results[Qval]["max"], '.-', color=colors[Qval], label=f"Q={Qval:+.2f} max")
        ax.plot(a0_list, results[Qval]["avg"], '.:', color=colors[Qval], label=f"Q={Qval:+.2f} avg")

    ax.set_xlabel("a0 (Å⁻¹)")
    ax.set_ylabel("Error (eV)")
    ax.set_title("Fit error vs a0 (region z > 2 Å)")
    ax.grid(True, ls=":", lw=0.5)
    ax.legend()
    plt.tight_layout()



if __name__ == "__main__":
     # parameters (kept very short names)
    Lx, Lz, dx, z0 = 5.0, 10.0, 0.1, 1.5
    npbc = 2  # ± images in x

    # atoms (two examples, charges zero – coulomb part will be zero)
    Q = 1.0
    atoms = [
    #   (x,y,z)      (R0, E0, Q, aMorse)    color    
        ((0.0 , 0.0, 0.0), (1.4, 0.01, +Q, 1.6), "m" ),  # Na
        ((Lx/2, 0.0, 0.0), (2.2, 0.01, -Q, 1.6), "g" ),  # Cl
        ((0.0 , 0.0,-2.5), (1.4, 0.01, -Q, 1.6), "g" ),
        ((Lx/2, 0.0,-2.5), (2.2, 0.01, +Q, 1.6), "m" ),
    ]
    apos   = np.array([at[0] for at in atoms])
    apars  = np.array([at[1] for at in atoms])
    colors = [at[2] for at in atoms]

    apos, apars, colors = mult_atoms(apos, apars, colors, nPBC=(npbc,0,0), Ls=(Lx,1.,1.), corners={1,3} )
    extent = [-Lx/2, Lx/2, z0, z0 + Lz]   ;print("extent:", extent)
    xs, zs = np.arange(extent[0], extent[1], dx), np.arange(extent[2], extent[3], dx)
    X, Z   = np.meshgrid(xs, zs, indexing="ij")

    nbx, nbz, a0 = 8, 3, 0.4

    a0_list = [0.20, 0.25, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]
    Qs    = [0.0, -0.1, 0.1]

    scan_Qs(apos, apars, X, Z, Lx=Lx, nbx=nbx, nbz=nbz, a0=a0)
    #scan_basis_a0(apos, apars, X, Z, Q_list=Qs, a0_list=a0_list, nbx=nbx, nbz=nbz, Lx=Lx)


    plt.show()
