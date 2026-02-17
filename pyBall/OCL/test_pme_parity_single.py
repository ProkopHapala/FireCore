import os
import sys
import numpy as np

_repo_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "../.."))
if _repo_root not in sys.path: sys.path.insert(0, _repo_root)

from pyBall.OCL.pauli_ocl import PauliSolverCL
from pyBall.OCL.HubbardSolver import HubbardSolver


def eval_manybody_energies_dense(Esite, Wij=None, W_scalar=0.0):
    """Match PME.cl energy convention: sum_i n_i*Esite[i] + sum_{i<j} n_i n_j * (Wij[i,j] or W_scalar)."""
    n_sites = len(Esite)
    n_states = 1 << n_sites
    Es = np.zeros(n_states, dtype=np.float64)
    for s in range(n_states):
        e = 0.0
        for i in range(n_sites):
            if (s >> i) & 1:
                e += float(Esite[i])
        for i in range(n_sites):
            if (s >> i) & 1:
                for j in range(i + 1, n_sites):
                    if (s >> j) & 1:
                        e += float(Wij[i, j] if Wij is not None else W_scalar)
        Es[s] = e
    return Es


def main():
    np.set_printoptions(linewidth=np.inf, precision=8, suppress=True)

    # ------------------------------------------------------------------
    # Shared 4-site embedding (use only first 2 sites, other 2 are "spectators")
    # ------------------------------------------------------------------
    n_sites = 4
    pSites = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0, 0.0],
        [1e3, 0.0, 0.0, 1e3],
        [0.0, 1e3, 0.0, 1e3],
    ], dtype=np.float32)

    # tip position and bias
    tip_pos = np.array([[2.5, 0.0, 3.0]], dtype=np.float32)
    Vbias = 0.50
    Vtips = np.array([Vbias], dtype=np.float32)

    # rotations identity for PME.cl
    rots = np.tile(np.eye(3, dtype=np.float32), (n_sites, 1, 1))

    # multipole coefficients: PME.cl expects 10 coefficients, use monopole only
    cs = np.zeros(10, dtype=np.float32)
    cs[0] = 1.0
    order = 0

    # parameters in pauli_ocl.py: [Rtip,zV0,zVd,Esite_ref,beta,Gamma,W,bMirror,bRamp]
    Rtip = 3.0
    zV0 = -0.5
    zVd = 0.0
    Esite_ref = 0.0
    beta = 0.5
    Gamma = 1.0
    W_scalar = 0.0
    bMirror = 1.0
    bRamp = 1.0
    params = np.array([Rtip, zV0, zVd, Esite_ref, beta, Gamma, W_scalar, bMirror, bRamp], dtype=np.float32)

    # ------------------------------------------------------------------
    # Full PME path (PME.cl)
    # ------------------------------------------------------------------
    pauli = PauliSolverCL(nSingle=n_sites)
    pauli.set_lead(0, mu=0.0, temp=4.0 * 8.61733326214511e-5)  # eV
    pauli.set_lead(1, mu=Vbias, temp=4.0 * 8.61733326214511e-5)  # eV

    curr_full, Es_full, Ts_full, P_full, StateEs_full, K_full, CurMat_full = pauli.scan_current_tip(
        pTips=tip_pos,
        Vtips=Vtips,
        pSites=pSites,
        params=params,
        order=order,
        cs=cs,
        rots=rots,
        return_probs=True,
        return_state_energies=True,
        return_curmat=True,
    )

    Es_full = Es_full.reshape(1, n_sites)[0]
    Ts_full = Ts_full.reshape(1, n_sites)[0]
    P_full = P_full.reshape(1, 16)[0]
    StateEs_full = StateEs_full.reshape(1, 16)[0]

    # ------------------------------------------------------------------
    # Hubbard precalc path (hubbard.cl)
    # ------------------------------------------------------------------
    # Hubbard precalc expects posE float4 with E0 in w, and pTips float4 with Vbias in w
    posE = pSites.copy()
    pTips4 = np.zeros((1, 4), dtype=np.float32)
    pTips4[0, :3] = tip_pos[0]
    pTips4[0, 3] = Vbias

    build_opts = ["-DDBG_PRECALC=1"]
    hub = HubbardSolver(build_options=build_opts)

    # Hubbard multipoleCoefs are per-site, nMulti is number of coeffs per site
    multipoleCoefs = np.zeros((n_sites, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0

    hub_params = {
        "Rtip": float(Rtip),
        "zMirror": float(zV0),
        "zRampOffset": float(zVd),
        "Thop_decay": float(beta),
        "Thop_amp": 1.0,
    }

    Esite_hub, Tsite_hub = hub.precalc_esite_thop(
        posE=posE,
        pTips=pTips4,
        rots=rots,
        multipoleCoefs=multipoleCoefs,
        bMirror=True,
        bRamp=True,
        params=hub_params,
    )

    Esite_hub = Esite_hub.reshape(1, n_sites)[0]
    Tsite_hub = Tsite_hub.reshape(1, n_sites)[0]

    # ------------------------------------------------------------------
    # Comparisons: Esite/Thop
    # ------------------------------------------------------------------
    print("\n=== Tip interaction parity (PME.cl vs hubbard.cl precalc) ===")
    print("Es_full:", Es_full)
    print("Esite_hub:", Esite_hub)
    print("dE:", Esite_hub - Es_full)
    print("max|dE|:", np.max(np.abs(Esite_hub - Es_full)))

    print("\nTs_full:", Ts_full)
    print("Tsite_hub:", Tsite_hub)
    print("dT:", Tsite_hub - Ts_full)
    print("max|dT|:", np.max(np.abs(Tsite_hub - Ts_full)))

    # ------------------------------------------------------------------
    # Many-body energies parity (derived from Esite only, Wij=0)
    # ------------------------------------------------------------------
    E_many = eval_manybody_energies_dense(Es_full, Wij=None, W_scalar=0.0)
    print("\n=== Many-body energy parity (derived vs PME.cl StateEs) ===")
    print("StateEs_full:", StateEs_full)
    print("E_many:", E_many)
    print("dE_state:", E_many - StateEs_full)
    print("max|dE_state|:", np.max(np.abs(E_many - StateEs_full)))

    # ground state index from both
    i_gs_full = int(np.argmin(StateEs_full))
    i_gs_derived = int(np.argmin(E_many))
    print("\n=== Ground state (argmin energy) ===")
    print("i_gs_full   :", i_gs_full, "E=", StateEs_full[i_gs_full])
    print("i_gs_derived:", i_gs_derived, "E=", E_many[i_gs_derived])

    # show probabilities/current from full solver (hubbard PME parity not implemented yet)
    print("\n=== Full PME outputs (reference) ===")
    print("current_full:", curr_full)
    print("probs:", P_full)


if __name__ == "__main__":
    main()
