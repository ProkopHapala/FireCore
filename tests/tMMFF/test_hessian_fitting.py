#!/usr/bin/env python3
"""
Test Hessian parameter fitting for UFF force field.
Workflow:
  1. Init UFF (bMMFF=True, bUFF=True)
  2. getBuffs_UFF() - print all topology/params
  3. setSwitchesUFF to enable bonds+angles, disable non-bonded
  4. Relax geometry with run()
  5. getHessianContext() - print all basis matrices, atom indices, params
  6. assembleHessianFromParams() - print assembled Hessian
  7. Fit parameters (linear LS) and compare

Usage:
    cd tests/tMMFF && ./run.sh
    or: cd tests/tMMFF && python3 test_hessian_fitting.py
"""

import os, sys
import numpy as np
sys.path.append("../../")
from pyBall import MMFF
from pyBall import FTIR

data_dir = "../../cpp/common_resources"

def test_hessian_fitting(xyz_path="../../web/common_resources/mol/H2O.mol2", lam=1e-6):
    print("=" * 70)
    print("Hessian Parameter Fitting Test")
    print("=" * 70)

    # 1. Init UFF --- bMMFF=True is REQUIRED to use UFF (see test_UFF_multi.py)
    print(f"\n--- 1. Init UFF from {xyz_path} ---")
    MMFF.setVerbosity(verbosity=1, idebug=0)  # idebug=1 triggers buggy bonParams[atom_idx] in evalBonds debug print
    MMFF.init(
        xyz_name       = xyz_path,
        surf_name      = None,
        smile_name     = None,
        sElementTypes  = os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes     = os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes     = os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes    = os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes = os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF  = True,   # MUST be True to use UFF
        bUFF   = True,
        bEpairs= False,
        nPBC   = (1,1,1),
    )

    # 2. Get UFF buffers (topology + params) and print everything
    print("\n--- 2. getBuffs_UFF() ---")
    MMFF.getBuffs_UFF()
    print(f"natoms={MMFF.natoms}  nbonds={MMFF.nbonds}  nangles={MMFF.nangles}  ndihedrals={MMFF.ndihedrals}  ninversions={MMFF.ninversions}")
    print(f"apos:\n{MMFF.apos}")
    print(f"bonAtoms:\n{MMFF.bonAtoms}")
    print(f"bonParams (k, l0):\n{MMFF.bonParams}")
    print(f"angAtoms:\n{MMFF.angAtoms}")
    print(f"angParams (k, theta0, ...):\n{MMFF.angParams}")
    print(f"neighs:\n{MMFF.neighs}")
    print(f"neighBs:\n{MMFF.neighBs}")

    # 3. Set switches: bonds+angles ON, non-bonded OFF
    print("\n--- 3. setSwitchesUFF ---")
    MMFF.setSwitches(NonBonded=-1, SurfAtoms=-1, GridFF=-1)
    MMFF.setSwitchesUFF(DoBond=1, DoAngle=1, DoDihedral=-1, DoInversion=-1, DoAssemble=1, SubtractBondNonBond=-1, ClampNonBonded=-1)

    # 4. Relax geometry
    print("\n--- 4. Relax geometry ---")
    MMFF.setTrjName("trj_hessian_test.xyz", savePerNsteps=10)
    nconv = MMFF.run(nstepMax=1000, dt=0.05, Fconv=1e-6, ialg=2, damping=0.1)
    print(f"Converged after {nconv} steps")
    print(f"apos after relax:\n{MMFF.apos}")
    print(f"fapos after relax (should be ~0):\n{MMFF.fapos}")
    fmax = np.max(np.abs(MMFF.fapos))
    print(f"max|F| = {fmax:.3e}")
    assert fmax < 0.1, f"Relaxation did not converge! max|F|={fmax:.3e}"

    # 5. Get Hessian context - print all intermediate arrays
    print("\n--- 5. getHessianContext ---")
    ctx = MMFF.getHessianContext(get_bases=True, get_atoms=True, get_params=True)
    nb, na, nat = ctx['n_bonds'], ctx['n_angles'], ctx['n_atoms']
    print(f"n_atoms={nat}  n_bonds={nb}  n_angles={na}")
    print(f"bond_atoms:\n{ctx['bond_atoms']}")
    print(f"bond_params (k):\n{ctx['bond_params']}")
    print(f"angle_atoms:\n{ctx['angle_atoms']}")
    print(f"angle_params (k):\n{ctx['angle_params']}")
    print(f"bond_bases shape: {ctx['bond_bases'].shape}")
    for ib in range(nb):
        print(f"  bond[{ib}] atoms={ctx['bond_atoms'][ib]}  k={ctx['bond_params'][ib]:.6f}")
        B = ctx['bond_bases'][ib]   # (4,3,3): [Bii, Bij, Bji, Bjj]
        for label, block in zip(['Bii','Bij','Bji','Bjj'], B):
            print(f"    {label}:\n{block}")
    print(f"angle_bases shape: {ctx['angle_bases'].shape}")
    for ia in range(na):
        print(f"  angle[{ia}] atoms={ctx['angle_atoms'][ia]}  k={ctx['angle_params'][ia]:.6f}")
        print(f"    B[{ia}] (9x9):\n{ctx['angle_bases'][ia]}")

    # 6. Assemble Hessian from current params and print
    print("\n--- 6. assembleHessianFromParams ---")
    k_current = np.concatenate([ctx['bond_params'], ctx['angle_params']])
    print(f"k vector: {k_current}")
    H = MMFF.assembleHessianFromParams(k_current, nat)
    H = 0.5 * (H + H.T)
    print(f"Assembled Hessian ({H.shape}):\n{H}")
    print(f"Hessian norm: {np.linalg.norm(H,'fro'):.6e}")
    assert np.all(np.isfinite(H)), "Hessian contains NaN/Inf!"

    # 7. Fit parameters
    print(f"\n--- 7. Fit parameters (λ={lam}) ---")
    result = FTIR.fit_hessian_parameters_ctx(H, ctx, lam=lam, method='ridge')
    k_bonds_fit  = result['k_bonds']
    k_angles_fit = result['k_angles']
    print(f"Fitted bond k:   {k_bonds_fit}")
    print(f"Original bond k: {ctx['bond_params']}")
    print(f"Bond k error:    {np.linalg.norm(k_bonds_fit - ctx['bond_params']):.3e}")
    print(f"Fitted angle k:  {k_angles_fit}")
    print(f"Original angle k:{ctx['angle_params']}")
    print(f"Angle k error:   {np.linalg.norm(k_angles_fit - ctx['angle_params']):.3e}")
    print(f"Residual: {result['residual']:.3e}  R²: {result['r2']:.6f}  cond#: {result['condition_number']:.2e}")

    MMFF.clear()
    print("\n" + "=" * 70)
    print("Test DONE")
    print("=" * 70)
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Test Hessian parameter fitting")
    parser.add_argument("-i", "--input", default="../../web/common_resources/mol/H2O.mol2")
    parser.add_argument("--lam", type=float, default=1e-6)
    args = parser.parse_args()
    test_hessian_fitting(args.input, args.lam)
