#!/usr/bin/env python3
"""GPU4PySCF single-point test on adamantane (C10H16) and Si10H16."""

import time
from pyscf import gto, dft

def run_cpu(mol, xc='b3lyp'):
    mf = dft.RKS(mol); mf.xc = xc
    t0 = time.time()
    e = mf.kernel()
    return e, time.time() - t0

def run_gpu(mol, xc='b3lyp'):
    from gpu4pyscf.dft import rks
    mf = rks.RKS(mol); mf.xc = xc
    t0 = time.time()
    e = mf.kernel()
    return e, time.time() - t0

def make_adamantane():
    """C10H16 from mol2 file."""
    atom = """
    C  -1.5847   1.1071  -0.2162
    C  -0.6083   0.8543   0.9552
    H  -3.6133  -0.5923   0.6075
    H  -2.9194  -0.7530   3.0217
    C  -2.9750   1.4862   0.3421
    H  -3.4095   1.6322   3.6622
    H  -1.7269   1.1832   4.1010
    C  -1.1444  -0.3003   1.8315
    C  -2.8525   2.7679   1.1972
    H  -3.6769   1.6676  -0.5007
    C  -3.5094   0.3309   1.2192
    H  -2.4807   3.6074   0.5694
    H  -3.8519   3.0576   1.5905
    C  -1.8782   2.5162   2.3700
    C  -2.5346   0.0779   2.3912
    C  -2.4146   1.3606   3.2452
    C  -0.4883   2.1354   1.8114
    H  -1.2200  -1.2307   1.2266
    H  -0.4397  -0.4983   2.6689
    H  -4.5166   0.5894   1.6139
    H  -1.1992   1.9278  -0.8604
    H  -1.6660   0.1937  -0.8456
    H   0.2229   1.9650   2.6495
    H  -0.0876   2.9682   1.1928
    H   0.3914   0.5817   0.5527
    H  -1.7912   3.4383   2.9848
    """
    return gto.M(atom=atom, basis='6-31g(d)', verbose=0)

def make_si10h():
    """Si10H16 from XYZ file."""
    atom = """
    Si  -2.21321   1.84157  -0.20380
    Si  -2.96215   2.08274   1.88455
    Si  -2.70958  -0.20602  -0.93874
    Si  -0.00022   2.12451  -0.23747
    Si   0.96718   0.60524   1.08035
    Si   0.21709   0.85048   3.16725
    Si   0.46687  -1.44249   0.34833
    Si  -1.74615  -1.72677   0.38008
    Si  -2.49429  -1.48146   2.46786
    Si  -1.99564   0.56510   3.20434
    H   -2.84952   2.84252  -1.07178
    H   -2.62994   3.43025   2.36845
    H   -4.41887   1.88895   1.90561
    H   -4.16725  -0.39232  -0.91295
    H   -2.20957  -0.36451  -2.31176
    H    0.49078   1.95976  -1.61291
    H    0.32607   3.47137   0.25176
    H    2.42476   0.79224   1.05745
    H    0.54234   2.20186   3.64563
    H    0.85245  -0.15302   4.03369
    H    0.95852  -1.60257  -1.02748
    H    1.09761  -2.44040   1.22391
    H   -2.07477  -3.07522  -0.10390
    H   -1.85358  -2.48062   3.33477
    H   -3.95265  -1.66405   2.48740
    H   -2.48844   0.72553   4.57958
    """
    return gto.M(atom=atom, basis='6-31g(d)', verbose=0)

if __name__ == '__main__':
    print("=" * 60)
    print("GPU4PySCF single-point: larger systems")
    print("=" * 60)

    for name, make_mol in [("Adamantane C10H16", make_adamantane),
                           ("Si10H16", make_si10h)]:
        print(f"\n--- {name} ---")
        mol = make_mol()
        nao = mol.nao
        print(f"Atoms = {mol.natm}, basis functions (nao) = {nao}")

        e_cpu, t_cpu = run_cpu(mol)
        print(f"CPU  E = {e_cpu:14.8f} Ha   time = {t_cpu:.3f} s")

        e_gpu, t_gpu = run_gpu(mol)
        print(f"GPU  E = {e_gpu:14.8f} Ha   time = {t_gpu:.3f} s")

        diff = abs(e_cpu - e_gpu)
        print(f"|ΔE| = {diff:.2e} Ha   speedup = {t_cpu/t_gpu:.2f}x")
        assert diff < 1e-6, "Energy mismatch!"

    print("\n" + "=" * 60)
    print("All tests passed.")
    print("=" * 60)
