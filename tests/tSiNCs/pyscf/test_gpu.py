#!/usr/bin/env python3
"""Minimal GPU4PySCF test on H2O and CH4."""

import time
import numpy as np
from pyscf import gto, dft

def run_cpu(mol, xc='b3lyp'):
    """Run CPU DFT single-point."""
    mf = dft.RKS(mol)
    mf.xc = xc
    t0 = time.time()
    e = mf.kernel()
    t = time.time() - t0
    return e, t

def run_gpu(mol, xc='b3lyp'):
    """Run GPU DFT single-point via gpu4pyscf."""
    from gpu4pyscf.dft import rks
    mf = rks.RKS(mol)
    mf.xc = xc
    t0 = time.time()
    e = mf.kernel()
    t = time.time() - t0
    return e, t

def make_h2o():
    mol = gto.M(
        atom='''
        O  0.000000  0.000000  0.117300
        H  0.756950  0.000000 -0.469200
        H -0.756950  0.000000 -0.469200
        ''',
        basis='6-31g(d)',
        verbose=0,
    )
    return mol

def make_ch4():
    mol = gto.M(
        atom='''
        C  0.000000  0.000000  0.000000
        H  0.626700  0.626700  0.626700
        H -0.626700 -0.626700  0.626700
        H -0.626700  0.626700 -0.626700
        H  0.626700 -0.626700 -0.626700
        ''',
        basis='6-31g(d)',
        verbose=0,
    )
    return mol

if __name__ == '__main__':
    print("=" * 60)
    print("GPU4PySCF smoke test")
    print("=" * 60)

    # Check GPU availability
    try:
        import gpu4pyscf
        from gpu4pyscf.dft import rks
        print("gpu4pyscf imported OK")
    except Exception as exc:
        print(f"gpu4pyscf import FAILED: {exc}")
        raise

    for name, make_mol in [("H2O", make_h2o), ("CH4", make_ch4)]:
        print(f"\n--- {name} ---")
        mol = make_mol()

        e_cpu, t_cpu = run_cpu(mol)
        print(f"CPU  E = {e_cpu:14.8f} Ha   time = {t_cpu:.3f} s")

        e_gpu, t_gpu = run_gpu(mol)
        print(f"GPU  E = {e_gpu:14.8f} Ha   time = {t_gpu:.3f} s")

        diff = abs(e_cpu - e_gpu)
        print(f"|E_cpu - E_gpu| = {diff:.2e} Ha")
        assert diff < 1e-6, "Energy mismatch between CPU and GPU!"

    print("\n" + "=" * 60)
    print("All tests passed.")
    print("=" * 60)
