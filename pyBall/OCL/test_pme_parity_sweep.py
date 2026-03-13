import os
import sys
import numpy as np
import matplotlib.pyplot as plt

_repo_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "../.."))
if _repo_root not in sys.path: sys.path.insert(0, _repo_root)

from pyBall.OCL.pauli_ocl import PauliSolverCL
from pyBall.OCL.HubbardSolver import HubbardSolver


def main():
    np.set_printoptions(linewidth=np.inf, precision=6, suppress=True)

    # Effective 2-site embedded into 4-site (PME.cl hardcoded N_SITES=4)
    n_sites = 4
    pSites = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0, 0.0],
        [1e3, 0.0, 0.0, 1e3],
        [0.0, 1e3, 0.0, 1e3],
    ], dtype=np.float32)

    # tip line scan (x sweep) with fixed bias
    zTip = 3.0
    Vbias = 0.50
    xs = np.linspace(-2.0, 7.0, 91).astype(np.float32)
    tip_pos = np.zeros((len(xs), 3), dtype=np.float32)
    tip_pos[:, 0] = xs
    tip_pos[:, 1] = 0.0
    tip_pos[:, 2] = zTip
    Vtips = np.full(len(xs), Vbias, dtype=np.float32)

    rots = np.tile(np.eye(3, dtype=np.float32), (n_sites, 1, 1))

    cs = np.zeros(10, dtype=np.float32)
    cs[0] = 1.0
    order = 0

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

    # Full solver
    pauli = PauliSolverCL(nSingle=n_sites)
    pauli.set_lead(0, mu=0.0, temp=4.0 * 8.61733326214511e-5)
    pauli.set_lead(1, mu=Vbias, temp=4.0 * 8.61733326214511e-5)

    curr_full, Es_full, Ts_full, P_full, StateEs_full, K_full, CurMat_full = pauli.scan_current_tip(
        pTips=tip_pos,
        Vtips=Vtips,
        pSites=pSites,
        params=params,
        order=order,
        cs=cs,
        rots=rots,
        return_probs=False,
        return_state_energies=True,
    )

    # Hubbard precalc
    pTips4 = np.zeros((len(xs), 4), dtype=np.float32)
    pTips4[:, :3] = tip_pos
    pTips4[:, 3] = Vbias

    hub = HubbardSolver(build_options=["-DDBG_PRECALC=0"])
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
        posE=pSites,
        pTips=pTips4,
        rots=rots,
        multipoleCoefs=multipoleCoefs,
        bMirror=True,
        bRamp=True,
        params=hub_params,
    )

    # Compare only first 2 sites (the relevant ones)
    Es_full = Es_full.reshape(len(xs), n_sites)
    Ts_full = Ts_full.reshape(len(xs), n_sites)

    fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)

    axs[0, 0].plot(xs, Es_full[:, 0], label="PME E0")
    axs[0, 0].plot(xs, Esite_hub[:, 0], '--', label="Hub E0")
    axs[0, 0].set_ylabel("Esite[0] (eV)")
    axs[0, 0].legend()

    axs[0, 1].plot(xs, Es_full[:, 1], label="PME E1")
    axs[0, 1].plot(xs, Esite_hub[:, 1], '--', label="Hub E1")
    axs[0, 1].set_ylabel("Esite[1] (eV)")
    axs[0, 1].legend()

    axs[1, 0].plot(xs, Ts_full[:, 0], label="PME T0")
    axs[1, 0].plot(xs, Tsite_hub[:, 0], '--', label="Hub T0")
    axs[1, 0].set_ylabel("T[0]")
    axs[1, 0].set_xlabel("tip x (A)")
    axs[1, 0].legend()

    axs[1, 1].plot(xs, Ts_full[:, 1], label="PME T1")
    axs[1, 1].plot(xs, Tsite_hub[:, 1], '--', label="Hub T1")
    axs[1, 1].set_ylabel("T[1]")
    axs[1, 1].set_xlabel("tip x (A)")
    axs[1, 1].legend()

    fig.suptitle(f"Tip interaction parity (x sweep, V={Vbias:.2f})")
    fig.tight_layout()

    out = "/home/prokop/git/FireCore/figs/parity_tip_xsweep.png"
    fig.savefig(out, dpi=150)
    print("saved", out)

    # quick numeric summary
    dE = Esite_hub[:, :2] - Es_full[:, :2]
    dT = Tsite_hub[:, :2] - Ts_full[:, :2]
    print("max|dE|=", np.max(np.abs(dE)))
    print("max|dT|=", np.max(np.abs(dT)))


if __name__ == "__main__":
    main()
