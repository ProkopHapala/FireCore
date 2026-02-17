import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

_repo_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "../.."))
if _repo_root not in sys.path: sys.path.insert(0, _repo_root)

from pyBall.OCL.pauli_ocl import PauliSolverCL
from pyBall.OCL.HubbardSolver import HubbardSolver

kB_eV = 8.61733326214511e-5


def make_basis_u4(nTips, nSite, states):
    """states: iterable of int bitmasks (up to 128 bits). Returns int32 (nTips,N,4)."""
    states = np.asarray(states, dtype=np.uint32)
    N = states.size
    if nSite > 128:
        raise ValueError("nSite>128 not supported")
    out = np.zeros((nTips, N, 4), dtype=np.int32)
    # only first word used for nSite<=32
    out[:, :, 0] = states.reshape(1, N)
    return out


def make_zero_sparse_W(nSite, nMaxNeigh=1):
    W_val = np.zeros((nSite, nMaxNeigh), dtype=np.float32)
    W_idx = np.zeros((nSite, nMaxNeigh), dtype=np.int32)
    nNeigh = np.zeros(nSite, dtype=np.int32)
    return (W_val, W_idx, nNeigh)


def eval_state_energies(Esite, W_scalar=0.0):
    nSite = Esite.size
    nStates = 1 << nSite
    Es = np.zeros(nStates, dtype=np.float64)
    for s in range(nStates):
        e = 0.0
        for i in range(nSite):
            if (s >> i) & 1:
                e += float(Esite[i])
        # no Wij for these parity tests by default
        Es[s] = e
    return Es


def curmat_to_isite(CurMat, nSite):
    """Sum PME.cl curmat contributions (only b->c added edges stored) into per-site array."""
    nStates = CurMat.shape[0]
    isite = np.zeros(nSite, dtype=np.float64)
    for b in range(nStates):
        mb = b
        for c in range(nStates):
            v = CurMat[b, c]
            if v == 0.0:
                continue
            diff = b ^ c
            if diff == 0:
                continue
            site = int(np.log2(diff & -diff))
            if site < nSite:
                isite[site] += v
    return isite


def write_summary(path, meta, spreads):
    with open(path, "w") as f:
        f.write("# Parity test summary\n")
        f.write("\n## Parameters\n")
        f.write(json.dumps(meta, indent=2))
        f.write("\n\n## Spreads (min/max)\n")
        for k, (mn, mx) in spreads.items():
            f.write(f"{k}: {mn:.8g} {mx:.8g}\n")


def save_csv(path, header, arr):
    np.savetxt(path, arr, delimiter=",", header=header, comments="")


def run_scan(out_dir, tag, pSites, rots, multipole_cs, order, params_pme, hub_params, scan_var, scan_values, tip_base_pos, zTip, fixed_bias=0.5, T_K=4.0, Gamma=1.0, do_curmat=True):
    os.makedirs(out_dir, exist_ok=True)

    nSite = pSites.shape[0]
    nStates = 1 << nSite

    # prepare tips
    nTips = len(scan_values)
    tip_pos = np.repeat(np.array(tip_base_pos, dtype=np.float32).reshape(1, 3), nTips, axis=0)
    tip_pos[:, 2] = zTip

    if scan_var == "x":
        tip_pos[:, 0] = np.asarray(scan_values, dtype=np.float32)
        Vtips = np.full(nTips, float(fixed_bias), dtype=np.float32)
    elif scan_var == "V":
        Vtips = np.asarray(scan_values, dtype=np.float32)
    else:
        raise ValueError("scan_var must be 'x' or 'V'")

    # Full solver setup
    pauli = PauliSolverCL(nSingle=nSite)
    T_eV = T_K * kB_eV
    pauli.set_lead(0, mu=0.0, temp=T_eV)

    # set lead1 mu per point inside kernel via v_tips (PME.cl uses v_tips as mu1)

    # PME params: [Rtip, zV0, zVd, Esite_ref, beta, Gamma, W, bMirror, bRamp]
    # params_pme passed as np.float32 array

    # Hubbard setup
    hub = HubbardSolver(build_options=["-DDBG_PRECALC=0"])
    multipoleCoefs = np.zeros((nSite, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0

    # Run full solver: need Es/Ts always and probs/current
    curr_full, Es_full, Ts_full, P_full, StateEs_full, K_full, CurMat_full = pauli.scan_current_tip(
        pTips=tip_pos,
        Vtips=Vtips,
        pSites=pSites,
        params=params_pme,
        order=order,
        cs=multipole_cs,
        rots=rots,
        return_probs=True,
        return_state_energies=True,
        return_curmat=do_curmat,
    )

    Es_full = Es_full.reshape(nTips, nSite)
    Ts_full = Ts_full.reshape(nTips, nSite)
    P_full = P_full.reshape(nTips, nStates)
    StateEs_full = StateEs_full.reshape(nTips, nStates)
    if CurMat_full is not None:
        CurMat_full = CurMat_full.reshape(nTips, nStates, nStates)

    # Run hubbard precalc
    pTips4 = np.zeros((nTips, 4), dtype=np.float32)
    pTips4[:, :3] = tip_pos
    pTips4[:, 3] = Vtips

    Esite_hub, Tsite_hub = hub.precalc_esite_thop(
        posE=pSites,
        pTips=pTips4,
        rots=rots,
        multipoleCoefs=multipoleCoefs,
        bMirror=True,
        bRamp=True,
        params=hub_params,
    )

    # Hubbard PME dense with full basis
    states = np.arange(nStates, dtype=np.uint32)
    basis_u4 = make_basis_u4(nTips, nSite, states)
    W_sparse = make_zero_sparse_W(nSite, nMaxNeigh=1)

    curr_hub, probs_hub, isite_hub = hub.solve_pme_dense_fullbasis(
        W_sparse=W_sparse,
        Esite=Esite_hub.astype(np.float32),
        Tsite=Tsite_hub.astype(np.float32),
        basis_u4=basis_u4,
        nTips=nTips,
        nSite=nSite,
        mu_tip=0.0,          # NOTE: in hubbard kernels mu_tip is passed explicitly, but PME.cl uses Vtips as mu1
        T_env=T_eV,
        Gamma_sub=(Gamma / np.pi) ** 2,
        Gamma_tip=(Gamma / np.pi) ** 2,
        bAlloc=True,
    )

    probs_hub = probs_hub[:, :nStates]

    # Full solver per-site current
    isite_full = np.zeros((nTips, nSite), dtype=np.float64)
    if CurMat_full is not None:
        for it in range(nTips):
            isite_full[it] = curmat_to_isite(CurMat_full[it], nSite)

    # Ground states
    gsE_full = np.min(StateEs_full, axis=1)
    gsE_hub = np.min(np.array([eval_state_energies(Esite_hub[it]) for it in range(nTips)]), axis=1)

    # spreads
    spreads = {
        "Esite_full": (float(np.min(Es_full)), float(np.max(Es_full))),
        "Esite_hub": (float(np.min(Esite_hub)), float(np.max(Esite_hub))),
        "T_full": (float(np.min(Ts_full)), float(np.max(Ts_full))),
        "T_hub": (float(np.min(Tsite_hub)), float(np.max(Tsite_hub))),
        "P_full": (float(np.min(P_full)), float(np.max(P_full))),
        "P_hub": (float(np.min(probs_hub)), float(np.max(probs_hub))),
        "I_full": (float(np.min(curr_full)), float(np.max(curr_full))),
        "I_hub": (float(np.min(curr_hub)), float(np.max(curr_hub))),
    }

    meta = {
        "tag": tag,
        "scan_var": scan_var,
        "scan_min": float(np.min(scan_values)),
        "scan_max": float(np.max(scan_values)),
        "nTips": int(nTips),
        "nSite": int(nSite),
        "T_K": float(T_K),
        "T_eV": float(T_eV),
        "params_pme": [float(x) for x in params_pme],
        "hub_params": hub_params,
        "tip_base_pos": [float(x) for x in tip_base_pos],
        "zTip": float(zTip),
    }

    write_summary(os.path.join(out_dir, f"summary_{tag}.txt"), meta, spreads)

    # CSV export: scan + Esite/T + gsE + currents
    cols = [scan_var]
    data = [np.asarray(scan_values, dtype=np.float64)]
    for i in range(nSite):
        cols += [f"Efull_{i}", f"Ehub_{i}", f"Tfull_{i}", f"Thub_{i}", f"Isite_full_{i}", f"Isite_hub_{i}"]
        data += [Es_full[:, i], Esite_hub[:, i], Ts_full[:, i], Tsite_hub[:, i], isite_full[:, i], isite_hub[:, i]]
    cols += ["gsE_full", "gsE_hub", "I_full", "I_hub"]
    data += [gsE_full, gsE_hub, curr_full, curr_hub]

    out_mat = np.vstack(data).T
    save_csv(os.path.join(out_dir, f"scan_{tag}.csv"), ",".join(cols), out_mat)

    # Probabilities export (top few states by avg prob full)
    pavg = np.mean(P_full, axis=0)
    top = np.argsort(-pavg)[:min(8, nStates)]
    pdat = [np.asarray(scan_values, dtype=np.float64)]
    pcols = [scan_var]
    for s in top:
        pcols += [f"Pfull_{s}", f"Phub_{s}"]
        pdat += [P_full[:, s], probs_hub[:, s]]
    save_csv(os.path.join(out_dir, f"probs_{tag}.csv"), ",".join(pcols), np.vstack(pdat).T)

    # Plot Esite and T
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    for i in range(min(4, nSite)):
        axs[0, 0].plot(scan_values, Es_full[:, i], label=f"full E{i}")
        axs[0, 0].plot(scan_values, Esite_hub[:, i], '--', label=f"hub E{i}")
        axs[1, 0].plot(scan_values, Ts_full[:, i], label=f"full T{i}")
        axs[1, 0].plot(scan_values, Tsite_hub[:, i], '--', label=f"hub T{i}")
    axs[0, 0].set_ylabel("Esite")
    axs[1, 0].set_ylabel("T")
    axs[1, 0].set_xlabel(scan_var)
    axs[0, 0].legend(fontsize=8)
    axs[1, 0].legend(fontsize=8)

    # Plot probabilities (top)
    for s in top:
        axs[0, 1].plot(scan_values, P_full[:, s], label=f"full P{s}")
        axs[0, 1].plot(scan_values, probs_hub[:, s], '--', label=f"hub P{s}")
    axs[0, 1].set_ylabel("Prob")
    axs[0, 1].legend(fontsize=8)

    # Plot currents
    axs[1, 1].plot(scan_values, curr_full, label="full Itot")
    axs[1, 1].plot(scan_values, curr_hub, '--', label="hub Itot")
    for i in range(min(4, nSite)):
        axs[1, 1].plot(scan_values, isite_full[:, i], ':', label=f"full Isite{i}")
        axs[1, 1].plot(scan_values, isite_hub[:, i], '-.', label=f"hub Isite{i}")
    axs[1, 1].set_ylabel("Current")
    axs[1, 1].set_xlabel(scan_var)
    axs[1, 1].legend(fontsize=7)

    fig.suptitle(f"{tag} parity scan ({scan_var})")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"plot_{tag}.png"), dpi=150)
    plt.close(fig)

    # report numeric parity
    err = {
        "max_dE": float(np.max(np.abs(Esite_hub - Es_full))),
        "max_dT": float(np.max(np.abs(Tsite_hub - Ts_full))),
        "max_dP": float(np.max(np.abs(probs_hub - P_full))),
        "max_dI": float(np.max(np.abs(curr_hub - curr_full))),
        "max_dIsite": float(np.max(np.abs(isite_hub - isite_full))),
    }
    with open(os.path.join(out_dir, f"errors_{tag}.json"), "w") as f:
        json.dump(err, f, indent=2)

    return err, spreads


def main():
    np.set_printoptions(linewidth=np.inf, precision=6, suppress=True)

    # Common physical parameters
    Rtip = 3.0
    zV0 = -0.5
    zVd = 0.0
    beta = 0.5
    Gamma = 1.0
    bMirror = 1.0
    bRamp = 1.0
    Esite_ref = 0.0
    W_scalar = 0.0
    T_K = 4.0
    zTip = 3.0

    # PME.cl multipole expects 10 coeffs
    multipole_cs = np.zeros(10, dtype=np.float32)
    multipole_cs[0] = 1.0
    order = 0

    # hubbard precalc uses per-site coefs, pass nMulti=1 (monopole)
    hub_params = {
        "Rtip": float(Rtip),
        "zMirror": float(zV0),
        "zRampOffset": float(zVd),
        "Thop_decay": float(beta),
        "Thop_amp": 1.0,
    }

    base_tip = (2.5, 0.0, zTip)

    # Output roots
    out2 = os.path.join(_repo_root, "figs", "parity_2site")
    out4 = os.path.join(_repo_root, "figs", "parity_4site")

    # ------------------------------------------------------------------
    # 2-site embedded in 4-site
    # ------------------------------------------------------------------
    pSites2 = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0, 0.0],
        [1e3, 0.0, 0.0, 1e3],
        [0.0, 1e3, 0.0, 1e3],
    ], dtype=np.float32)
    rots2 = np.tile(np.eye(3, dtype=np.float32), (4, 1, 1))

    # PME params for full solver
    params_pme = np.array([Rtip, zV0, zVd, Esite_ref, beta, Gamma, W_scalar, bMirror, bRamp], dtype=np.float32)

    # x scan at fixed V
    xs = np.linspace(-2.0, 7.0, 91)
    Vfix = 0.50
    err_x2, spr_x2 = run_scan(os.path.join(out2, "xscan"), "2site_xscan", pSites2, rots2, multipole_cs, order, params_pme, hub_params, "x", xs, base_tip, zTip, fixed_bias=Vfix, T_K=T_K, Gamma=Gamma)

    # V scan at fixed x
    Vs = np.linspace(0.0, 1.0, 101)
    err_v2, spr_v2 = run_scan(os.path.join(out2, "Vscan"), "2site_Vscan", pSites2, rots2, multipole_cs, order, params_pme, hub_params, "V", Vs, base_tip, zTip, T_K=T_K, Gamma=Gamma)

    # ------------------------------------------------------------------
    # true 4-site (a simple line)
    # ------------------------------------------------------------------
    pSites4 = np.array([
        [0.0, 0.0, 0.0, 0.0],
        [5.0, 0.0, 0.0, 0.0],
        [10.0, 0.0, 0.0, 0.0],
        [15.0, 0.0, 0.0, 0.0],
    ], dtype=np.float32)
    rots4 = np.tile(np.eye(3, dtype=np.float32), (4, 1, 1))

    err_x4, spr_x4 = run_scan(os.path.join(out4, "xscan"), "4site_xscan", pSites4, rots4, multipole_cs, order, params_pme, hub_params, "x", xs, base_tip, zTip, fixed_bias=Vfix, T_K=T_K, Gamma=Gamma)
    err_v4, spr_v4 = run_scan(os.path.join(out4, "Vscan"), "4site_Vscan", pSites4, rots4, multipole_cs, order, params_pme, hub_params, "V", Vs, base_tip, zTip, T_K=T_K, Gamma=Gamma)

    print("\n=== DONE ===")
    print("2-site outputs:", out2)
    print("4-site outputs:", out4)


if __name__ == "__main__":
    main()
