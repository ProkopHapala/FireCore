#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

# Allow running from tests/pyFireball
sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import FireCore as fc
from pyBall import plotUtils as pu

def fit_molecular_plane(positions):
    """Return (origin, normal) for best-fit plane through atomic positions.

    origin  – centroid of all positions
    normal  – unit vector normal to the plane (from PCA, smallest-variance direction)
    """
    pos = np.asarray(positions, dtype=float)
    origin = pos.mean(axis=0)
    X = pos - origin
    # PCA via SVD: columns of Vt are principal directions
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    normal = Vt[-1]  # direction with smallest variance
    # Normalize and fix orientation to have positive z if possible (for consistency)
    normal /= np.linalg.norm(normal)
    if normal[2] < 0.0:
        normal = -normal
    return origin, normal

def build_plane_grid(origin, normal, size_x, size_y, nx, ny, z_height):
    """Build a 2D grid of points in a plane parallel to the molecular plane.

    The plane passes through `origin + z_height * normal` and is spanned by
    two in-plane orthonormal vectors (x_axis, y_axis).

    Returns
    -------
    points : (nx*ny, 3) array of 3D coordinates
    X_local, Y_local : (ny, nx) arrays of in-plane coordinates (for reference)
    """
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)

    # Choose an arbitrary vector not parallel to n to build an in-plane basis
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial, n)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])

    x_axis = trial - np.dot(trial, n) * n
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(n, x_axis)
    y_axis /= np.linalg.norm(y_axis)

    xs = np.linspace(-0.5 * size_x, 0.5 * size_x, nx)
    ys = np.linspace(-0.5 * size_y, 0.5 * size_y, ny)
    X_local, Y_local = np.meshgrid(xs, ys)  # shape (ny, nx)

    plane_origin = origin + z_height * n
    points = plane_origin[None, :] + X_local[..., None] * x_axis[None, None, :] + Y_local[..., None] * y_axis[None, None, :]
    points = points.reshape(-1, 3)
    return points, X_local, Y_local

def project_atoms_to_plane(apos, origin, normal, z_height):
    """Project atomic positions onto the same plane coordinates used for the sampling grid."""
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial, n)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    x_axis = trial - np.dot(trial, n) * n
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(n, x_axis)
    y_axis /= np.linalg.norm(y_axis)
    plane_origin = origin + z_height * n
    apos = np.asarray(apos, float)
    vecs = apos - plane_origin[None, :]
    xs = np.dot(vecs, x_axis)
    ys = np.dot(vecs, y_axis)
    apos_plane = np.stack([xs, ys, np.zeros_like(xs)], axis=1)
    return apos_plane

def evaluate_density_on_plane(points, f_den=1.0, f_den0=0.0):
    """Evaluate electronic density at given 3D points using FireBall (dens2points)."""
    return fc.dens2points(points, f_den=f_den, f_den0=f_den0)

def evaluate_orbital_on_plane(points, iMO=1, ikpoint=1):
    """Evaluate molecular orbital iMO at given 3D points using FireBall (orb2points)."""
    return fc.orb2points(points, iMO=iMO, ikpoint=ikpoint)

def plot_field_2d(values, X_local, Y_local, title, outfile, cmap="viridis", overlay=None):
    """Plot a scalar field defined on a regular (X_local, Y_local) grid and save to PNG."""
    ny, nx = X_local.shape
    field2d = values.reshape(ny, nx)

    plt.figure(figsize=(6, 5))
    extent = [X_local.min(), X_local.max(), Y_local.min(), Y_local.max()]
    # Use symmetric color range around zero so that mid-color (e.g. white) is at 0
    vmax = np.max(np.abs(field2d))
    vmin = -vmax
    plt.imshow(field2d, origin="lower", extent=extent, cmap=cmap, aspect="equal", vmin=vmin, vmax=vmax)
    plt.colorbar()
    if overlay is not None:
        overlay()
    plt.xlabel("x in plane [Å]")
    plt.ylabel("y in plane [Å]")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate FireBall density and MOs on a plane above a planar molecule.")
    parser.add_argument("--xyz",         type=str,   default=os.path.join("..", "..", "cpp", "common_resources", "xyz", "PTCDA.xyz"), help="Path to input XYZ file (planar molecule). Default: PTCDA.xyz in cpp/common_resources/xyz.")
    parser.add_argument("--z",           type=float, default=1.0,   help="Height above the molecular plane (in Å).")
    parser.add_argument("--size",        type=float, default=20.0,  help="Size of the sampling square in the plane (in Å). Applied to both x and y.")
    parser.add_argument("--nx",          type=int,   default=200,   help="Number of grid points along x in the plane.")
    parser.add_argument("--ny",          type=int,   default=200,   help="Number of grid points along y in the plane.")
    parser.add_argument("--orbital-min", type=int,   default=-10,   help="Index of first MO to sample (1-based), or offset from HOMO if --homo-range is set.")
    parser.add_argument("--orbital-max", type=int,   default=4,     help="Index of last MO to sample (1-based), or offset from HOMO if --homo-range is set. If 0 in absolute mode, a small range near the top of occupied states is chosen heuristically.")
    parser.add_argument("--homo-range",  type=int,   default=1,     help="Interpret orbital-min/max as offsets relative to HOMO (0/1).")
    parser.add_argument("--overlay-atoms",  type=int,default=0,     help="Plot atoms and bonds overlay on top of MO images (0/1).")
    parser.add_argument("--overlay-bonds",  type=int,default=1,     help="Plot atoms and bonds overlay on top of MO images (0/1).")
    parser.add_argument("--overlay-labels", type=int,default=0,     help="Show atom labels in overlay (0/1).")
    parser.add_argument("--nmax-scf",    type=int,   default=200,   help="Maximum number of SCF iterations for evalForce.")
    parser.add_argument("--verbosity",   type=int,   default=0,     help="FireBall verbosity level (0 = quiet).")
    parser.add_argument("--outdir",      type=str,   default="mo_plane_outputs", help="Directory to store PNG images.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("# Loading molecule from", args.xyz)
    mol = AtomicSystem(fname=args.xyz)
    print("Loaded molecule with", len(mol.apos), "atoms")

    # Fit plane to atomic coordinates
    origin, normal = fit_molecular_plane(mol.apos)
    print("Plane origin:", origin)
    print("Plane normal:", normal)

    # Initialize FireBall and run SCF via evalForce to converge density
    print("# Initializing FireBall (FireCore) ...")
    fc.setVerbosity(args.verbosity)
    fc.initialize(atomType=mol.atypes, atomPos=mol.apos, verbosity=args.verbosity)
    # Query number of orbitals and other dimensions from HS export helper
    dims = fc.get_HS_dims()
    norb = dims.norbitals
    print("================================== ")
    print("Number of orbitals (from get_HS_dims):", norb)
    print("# HS dimensions from get_HS_dims():")
    print("  natoms     =", dims.natoms)
    print("  norbitals  =", dims.norbitals)
    print("  nspecies   =", dims.nspecies)
    print("  neigh_max  =", dims.neigh_max)
    print("  numorb_max =", dims.numorb_max)
    print("  nsh_max    =", dims.nsh_max)
    print("  ME2c_max   =", dims.ME2c_max)
    print("  max_mu_dim =", (dims.max_mu_dim1, dims.max_mu_dim2, dims.max_mu_dim3))
    print("  mbeta_max  =", dims.mbeta_max)
    print("  nspecies_fdata =", dims.nspecies_fdata)
    print(" ================================== ")

    # Electron-count / HOMO estimate from FireBall (ztot) exposed via get_HS_dims
    nelec = dims.nelec
    homo_index_est = nelec // 2
    print("# Total electrons from FireBall (nelec):", nelec)
    print("# Estimated HOMO index (closed-shell, 1-based):", homo_index_est)

    print("# Running SCF via evalForce (nmax_scf = %d) ..." % args.nmax_scf)
    forces, energies = fc.evalForce(mol.apos, nmax_scf=args.nmax_scf)
    print("SCF done. Total energy Etot = %.6f eV" % energies[0])

    # Get orbital energies for the chosen k-point (assume ikp=1 / Gamma)
    eigen = fc.get_eigen(ikp=1, norb=norb)

    # Build sampling grid in plane at height z
    print("# Building sampling grid in plane: size = %.2f Å, z = %.2f Å, nx = %d, ny = %d" % (args.size, args.z, args.nx, args.ny))
    points, X_local, Y_local = build_plane_grid(origin, normal, args.size, args.size,args.nx, args.ny, args.z)

    # Evaluate density on the grid
    # print("# Evaluating density on plane using dens2points ...")
    # density_values = evaluate_density_on_plane(points, f_den=1.0, f_den0=0.0)
    # dens_title = "Density at z = %.2f Å" % args.z
    # dens_png = os.path.join(args.outdir, "density_z%.2f.png" % args.z)
    # plot_field_2d(density_values, X_local, Y_local, dens_title, dens_png, cmap="viridis")
    # print("Saved density map to", dens_png)

    # Decide which MOs to sample
    if norb is None or norb <= 0:
        print("WARNING: norbitals not available or zero; skipping MO projections.")
        i_min, i_max = 1, 0
    elif args.homo_range:
        # Interpret orbital-min/max as offsets with respect to HOMO index estimate
        i_min = max(1, homo_index_est + args.orbital_min)
        i_max = min(norb, homo_index_est + args.orbital_max)
    elif args.orbital_max <= 0:
        # Heuristic: sample a small set of highest-energy occupied / nearby orbitals
        i_max = norb
        i_min = max(1, norb - 3)
    else:
        i_min = max(1, args.orbital_min)
        i_max = min(norb, args.orbital_max)

    print("# Sampling molecular orbitals in range [%d, %d]" % (i_min, i_max))

    apos_plane = None
    bOverlay = args.overlay_atoms or args.overlay_bonds or args.overlay_labels
    if bOverlay:
        if getattr(mol, "bonds", None) is None:
            mol.findBonds(Rcut=3.0, RvdwCut=0.5)
        apos_plane = project_atoms_to_plane(mol.apos, origin, normal, args.z)

    for iMO in range(i_min, i_max + 1):
        print("  Evaluating MO", iMO, "on plane using orb2points ...")
        mo_values = evaluate_orbital_on_plane(points, iMO=iMO, ikpoint=1)
        rel_idx = iMO - homo_index_est
        E_MO = eigen[iMO-1] if (iMO-1) < len(eigen) else np.nan
        mo_title = "MO %d (HOMO%+d)  E = %.3f eV  z = %.2f Å" % (iMO, rel_idx, E_MO, args.z)
        mo_png = os.path.join(args.outdir, "MO%04d_z%.2f.png" % (iMO, args.z))
        def overlay_cb():
            if not bOverlay or apos_plane is None:
                return
            class _Sys:
                pass
            sys2 = _Sys()
            sys2.apos = apos_plane
            sys2.enames = mol.enames
            sys2.bonds = mol.bonds
            pu.plotSystem(sys2, bBonds=args.overlay_bonds, bAtoms=args.overlay_atoms, axes=(0,1), bLabels=bool(args.overlay_labels), labels=mol.enames)
        # Use a diverging colormap for signed orbital amplitude
        plot_field_2d(mo_values, X_local, Y_local, mo_title, mo_png, cmap="seismic", overlay=overlay_cb)
        print("    Saved", mo_png)

    print("# Done.")
