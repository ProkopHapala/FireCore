import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import FireCore as fc


def build_nccn_cell():
    a_CC = 1.42
    xa = a_CC * np.cos(np.pi/6)
    ya = a_CC * np.sin(np.pi/6)
    Lx = 2.0 * xa
    Ly = 20.0
    Lz = 20.0

    apos = np.array([
        [0.0,      0.0, 0.0],
        [xa,       ya,  0.0],
        [xa,       ya + a_CC, 0.0],
        [0.0, 2.0*ya + a_CC, 0.0],
    ], dtype=np.float64)
    apos[:,1] -= np.mean(apos[:,1])
    apos[:,1] += 0.5 * Ly
    apos[:,2] += 0.5 * Lz

    atypes = np.array([7, 6, 6, 7], dtype=np.int32)
    lvs = np.array([
        [Lx, 0.0, 0.0],
        [0.0, Ly,  0.0],
        [0.0, 0.0, Lz ],
    ], dtype=np.float64)
    return atypes, apos, lvs


def main():
    atypes, apos, lvs = build_nccn_cell()
    Lx = lvs[0,0]
    nk_scf = 1
    k_scf = np.zeros((nk_scf, 3), dtype=np.float64)
    w_scf = np.ones(nk_scf, dtype=np.float64) / nk_scf

    print("NCCN minimal ribbon cell")
    print("atypes=", atypes)
    print("apos=\n", apos)
    print("lvs=\n", lvs)
    print("Lx=", Lx, "BZ edge pi/Lx=", np.pi/Lx)

    fc.initialize(atomType=atypes, atomPos=apos, verbosity=1, lvs=lvs, kpoints=k_scf, kweights=w_scf)
    dims = fc.get_HS_dims(force_refresh=True)
    norb = int(dims.norbitals)
    print("norb=", norb)

    forces, energies = fc.evalForce(apos, nmax_scf=100)
    print("SCF energies=", energies)

    fc.assembleH(apos, iforce=0, Kscf=0)

    nk = 101
    kxs = np.linspace(-np.pi/Lx, np.pi/Lx, nk)
    bands = np.zeros((nk, norb), dtype=np.float64)
    for ik, kx in enumerate(kxs):
        fc.solveH(np.array([kx, 0.0, 0.0], dtype=np.float64), ikpoint=1)
        bands[ik, :] = fc.get_eigen(ikp=1, norb=norb)
        if (ik % 20) == 0:
            print("ik", ik, "kx", kx, "E[:4]", bands[ik,:min(4,norb)])

    out_npz = "NCCN_ribbon_bands.npz"
    out_png = "NCCN_ribbon_bands.png"
    np.savez(out_npz, kxs=kxs, bands=bands, atypes=atypes, apos=apos, lvs=lvs)

    plt.figure(figsize=(7,5))
    for ib in range(norb):
        plt.plot(kxs * Lx / np.pi, bands[:,ib], 'k-', lw=1.0)
    plt.xlabel(r"$k_x L_x / \pi$")
    plt.ylabel("E (eV)")
    plt.title("NCCN minimal graphene ribbon band structure")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    print("saved", out_npz, out_png)


if __name__ == "__main__":
    main()
