#!/usr/bin/env python3
"""
Parity tests: PauliSolverCL (PME.cl full solver) vs HubbardSolver (hubbard.cl dense PME).
Uses realistic parameters from Ruslan's JSON configs and geometry files.
Outputs CSV, TXT summaries, PNG plots into dedicated directories.
Supports 1D scans (x-sweep and V-sweep) and 2D scans (x-V maps).
"""
import os, sys, json, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_repo_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "../.."))
if _repo_root not in sys.path: sys.path.insert(0, _repo_root)

from pyBall.OCL.pauli_ocl import PauliSolverCL
from pyBall.OCL.HubbardSolver import HubbardSolver, make_sparse_W

kB_eV = 8.61733326214511e-5

# =====================================================================
#  Helpers
# =====================================================================

def load_geometry(path):
    """Load site geometry from tab-separated file: x y z E0 per line."""
    data = np.loadtxt(path, dtype=np.float32)
    if data.ndim == 1: data = data.reshape(1, -1)
    assert data.shape[1] >= 4, f"Expected >=4 columns in {path}, got {data.shape[1]}"
    return data[:, :4]  # (nSite, 4) = x,y,z,E0


def embed_2site_in_4(sites2):
    """Embed 2 real sites into 4-site array (PME.cl requires N_SITES=4). Spectators far away, E0=+100."""
    sites4 = np.zeros((4, 4), dtype=np.float32)
    sites4[:2] = sites2
    sites4[2] = [9999.0, 9999.0, 0.0, 100.0]
    sites4[3] = [9999.0, -9999.0, 0.0, 100.0]
    return sites4


def make_basis_u4(nTips, nSite, active_sites=None):
    """Basis states as int32 (nTips, nStates, 4). Only active sites vary; others stay 0."""
    if active_sites is None:
        active_sites = list(range(nSite))
    nActive = len(active_sites)
    nStates = 1 << nActive
    states = np.zeros(nStates, dtype=np.int32)
    for k in range(nStates):
        s = 0
        for b, site_idx in enumerate(active_sites):
            if k & (1 << b): s |= (1 << site_idx)
        states[k] = s
    out = np.zeros((nTips, nStates, 4), dtype=np.int32)
    out[:, :, 0] = states.reshape(1, nStates)
    return out, nStates


def make_uniform_Wij(nSite, W):
    """Dense uniform Wij: W for all i!=j pairs."""
    Wij = np.full((nSite, nSite), W, dtype=np.float32)
    np.fill_diagonal(Wij, 0.0)
    return Wij


def curmat_to_isite(CurMat, nSite):
    nStates = CurMat.shape[0]
    isite = np.zeros(nSite, dtype=np.float64)
    for b in range(nStates):
        for c in range(nStates):
            v = CurMat[b, c]
            if v == 0.0: continue
            diff = b ^ c
            if diff == 0 or bin(diff).count('1') != 1: continue
            site = (diff & -diff).bit_length() - 1
            if site < nSite: isite[site] += v
    return isite


def save_csv(path, header, arr):
    np.savetxt(path, arr, delimiter=",", header=header, comments="")


def spread(arr):
    return (float(np.min(arr)), float(np.max(arr)))


# =====================================================================
#  1D scan runner
# =====================================================================

def run_1d_scan(out_dir, tag, pSites, rots, cs, order, params_pme, hub_params, multipoleCoefs,
                scan_var, scan_values, tip_base_pos, zTip, fixed_bias, T_K, GammaS, GammaT,
                Wij_dense, real_sites, all_params_dict, active_sites=None):
    """Run a 1D parity scan (x or V) and produce CSV/PNG/TXT artifacts."""
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    nSite = pSites.shape[0]
    nStates_full = 1 << nSite  # full solver always uses 2^nSite
    nTips = len(scan_values)
    T_eV = T_K * kB_eV
    if active_sites is None: active_sites = list(range(nSite))

    # Prepare tip positions and biases
    tip_pos = np.tile(np.array(tip_base_pos, dtype=np.float32), (nTips, 1))
    if scan_var == "x":
        tip_pos[:, 0] = np.asarray(scan_values, dtype=np.float32)
        Vtips = np.full(nTips, float(fixed_bias), dtype=np.float32)
    elif scan_var == "V":
        Vtips = np.asarray(scan_values, dtype=np.float32)
    else:
        raise ValueError(f"Unknown scan_var={scan_var}")

    # --- Full solver (PME.cl) ---
    pauli = PauliSolverCL(nSingle=nSite)
    pauli.set_lead(0, mu=0.0, temp=T_eV)
    pauli.set_lead(1, mu=0.0, temp=T_eV)  # mu1 overridden by v_tips in kernel

    curr_full, Es_full, Ts_full, P_full, StateEs_full, K_full, CurMat_full = pauli.scan_current_tip(
        pTips=tip_pos[:, :3], Vtips=Vtips, pSites=pSites, params=params_pme,
        order=order, cs=cs, rots=rots,
        return_probs=True, return_state_energies=True, return_curmat=True, Wij=Wij_dense,
    )
    Es_full    = Es_full.reshape(nTips, nSite)
    Ts_full    = Ts_full.reshape(nTips, nSite)
    P_full     = P_full.reshape(nTips, nStates_full)
    StateEs_full = StateEs_full.reshape(nTips, nStates_full)
    CurMat_full  = CurMat_full.reshape(nTips, nStates_full, nStates_full) if CurMat_full is not None else None

    # --- Hubbard precalc ---
    hub = HubbardSolver(build_options=["-DDBG_PRECALC=0"])
    pTips4 = np.zeros((nTips, 4), dtype=np.float32)
    pTips4[:, :3] = tip_pos[:, :3]
    pTips4[:, 3]  = Vtips

    Esite_hub, Tsite_hub = hub.precalc_esite_thop(
        posE=pSites, pTips=pTips4, rots=rots, multipoleCoefs=multipoleCoefs,
        bMirror=True, bRamp=True, params=hub_params,
    )

    # --- Hubbard PME dense (active-sites basis, per-V loop for correct mu_tip) ---
    W_sparse = make_sparse_W(Wij_dense, nMaxNeigh=max(1, nSite - 1))
    basis_u4, nStates_hub = make_basis_u4(1, nSite, active_sites=active_sites)
    gamma_sub = (GammaS / np.pi) ** 2
    gamma_tip = (GammaT / np.pi) ** 2

    curr_hub   = np.zeros(nTips, dtype=np.float32)
    probs_hub_raw = np.zeros((nTips, nStates_hub), dtype=np.float32)
    isite_hub  = np.zeros((nTips, nSite), dtype=np.float32)

    for it in range(nTips):
        c, p, si = hub.solve_pme_dense_fullbasis(
            W_sparse=W_sparse,
            Esite=Esite_hub[it:it+1].astype(np.float32),
            Tsite=Tsite_hub[it:it+1].astype(np.float32),
            basis_u4=basis_u4,
            nTips=1, nSite=nSite,
            mu_tip=float(Vtips[it]),
            T_env=T_eV,
            Gamma_sub=gamma_sub, Gamma_tip=gamma_tip,
            bAlloc=(it == 0),
        )
        curr_hub[it]   = c[0]
        probs_hub_raw[it] = p[0, :nStates_hub]
        isite_hub[it]  = si[0]

    # Map hub probs to full-solver state indexing for comparison
    hub_state_indices = basis_u4[0, :nStates_hub, 0]  # actual state masks
    probs_hub = np.zeros((nTips, nStates_full), dtype=np.float32)
    for k in range(nStates_hub):
        si_full = int(hub_state_indices[k])
        if si_full < nStates_full: probs_hub[:, si_full] = probs_hub_raw[:, k]

    nStates = nStates_full  # for CSV/plot compatibility

    # Per-site current from full solver curmat
    isite_full = np.zeros((nTips, nSite), dtype=np.float64)
    if CurMat_full is not None:
        for it in range(nTips): isite_full[it] = curmat_to_isite(CurMat_full[it], nSite)

    # Ground state energies
    gsE_full = np.min(StateEs_full, axis=1)
    gs_idx_full = np.argmin(StateEs_full, axis=1)

    elapsed = time.time() - t0
    ns = real_sites  # indices of real (non-spectator) sites

    # --- Spreads ---
    spreads = {}
    for i in ns:
        spreads[f"Esite_full[{i}]"] = spread(Es_full[:, i])
        spreads[f"Esite_hub[{i}]"]  = spread(Esite_hub[:, i])
        spreads[f"T_full[{i}]"]     = spread(Ts_full[:, i])
        spreads[f"T_hub[{i}]"]      = spread(Tsite_hub[:, i])
    spreads["I_full"]  = spread(curr_full)
    spreads["I_hub"]   = spread(curr_hub)
    spreads["P_full"]  = spread(P_full)
    spreads["P_hub"]   = spread(probs_hub)
    spreads["gsE_full"] = spread(gsE_full)

    # --- Summary ---
    meta = dict(all_params_dict)
    meta.update({"tag": tag, "scan_var": scan_var, "scan_min": float(np.min(scan_values)),
                 "scan_max": float(np.max(scan_values)), "nPts": int(nTips), "nSite": int(nSite),
                 "nStates": int(nStates), "real_sites": [int(x) for x in ns],
                 "elapsed_sec": round(elapsed, 2)})
    with open(os.path.join(out_dir, f"summary_{tag}.txt"), "w") as f:
        f.write(f"# Parity test: {tag}\n\n## ALL Parameters\n")
        f.write(json.dumps(meta, indent=2))
        f.write("\n\n## Spreads (min, max)\n")
        for k, (mn, mx) in spreads.items(): f.write(f"  {k}: [{mn:.8g}, {mx:.8g}]\n")

    # --- CSV ---
    cols = [scan_var]; data = [np.asarray(scan_values, dtype=np.float64)]
    for i in ns:
        cols += [f"Efull_{i}", f"Ehub_{i}", f"Tfull_{i}", f"Thub_{i}", f"Isite_full_{i}", f"Isite_hub_{i}"]
        data += [Es_full[:, i], Esite_hub[:, i], Ts_full[:, i], Tsite_hub[:, i], isite_full[:, i], isite_hub[:, i]]
    cols += ["gsE_full", "I_full", "I_hub"]
    data += [gsE_full, curr_full.astype(np.float64), curr_hub.astype(np.float64)]
    save_csv(os.path.join(out_dir, f"scan_{tag}.csv"), ",".join(cols), np.vstack(data).T)

    # Probs CSV (top states)
    pavg = np.mean(np.abs(P_full), axis=0)
    top = np.argsort(-pavg)[:min(8, nStates)]
    pcols = [scan_var]; pdat = [np.asarray(scan_values, dtype=np.float64)]
    for s in top:
        pcols += [f"Pfull_{s}", f"Phub_{s}"]
        pdat += [P_full[:, s], probs_hub[:, s]]
    save_csv(os.path.join(out_dir, f"probs_{tag}.csv"), ",".join(pcols), np.vstack(pdat).T)

    # --- PLOTS (only real sites) ---
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    sv = np.asarray(scan_values)
    for i in ns:
        axs[0, 0].plot(sv, Es_full[:, i], label=f"full E{i}")
        axs[0, 0].plot(sv, Esite_hub[:, i], '--', label=f"hub E{i}")
        axs[1, 0].plot(sv, Ts_full[:, i], label=f"full T{i}")
        axs[1, 0].plot(sv, Tsite_hub[:, i], '--', label=f"hub T{i}")
    axs[0, 0].set_ylabel("Esite (eV)"); axs[0, 0].legend(fontsize=7)
    axs[1, 0].set_ylabel("T_hop");      axs[1, 0].legend(fontsize=7); axs[1, 0].set_xlabel(scan_var)

    for s in top[:6]:
        axs[0, 1].plot(sv, P_full[:, s], label=f"full P{s:04b}")
        axs[0, 1].plot(sv, probs_hub[:, s], '--', label=f"hub P{s:04b}")
    axs[0, 1].set_ylabel("Probability"); axs[0, 1].legend(fontsize=6, ncol=2)

    axs[1, 1].plot(sv, curr_full, 'b-', lw=1.5, label="full I_tot")
    axs[1, 1].plot(sv, curr_hub, 'r--', lw=1.5, label="hub I_tot")
    for i in ns:
        axs[1, 1].plot(sv, isite_full[:, i], ':', alpha=0.7, label=f"full Isite{i}")
        axs[1, 1].plot(sv, isite_hub[:, i], '-.', alpha=0.7, label=f"hub Isite{i}")
    axs[1, 1].set_ylabel("Current"); axs[1, 1].set_xlabel(scan_var); axs[1, 1].legend(fontsize=6, ncol=2)

    fig.suptitle(f"{tag}")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"plot_{tag}.png"), dpi=150)
    plt.close(fig)

    # --- Error report ---
    err = {}
    for i in ns:
        err[f"max_dE[{i}]"] = float(np.max(np.abs(Esite_hub[:, i] - Es_full[:, i])))
        err[f"max_dT[{i}]"] = float(np.max(np.abs(Tsite_hub[:, i] - Ts_full[:, i])))
    err["max_dP"] = float(np.max(np.abs(probs_hub - P_full)))
    err["max_dI"] = float(np.max(np.abs(curr_hub - curr_full)))
    with open(os.path.join(out_dir, f"errors_{tag}.json"), "w") as f:
        json.dump(err, f, indent=2)

    print(f"  [{tag}] done in {elapsed:.1f}s  max_dE={max(err[k] for k in err if 'dE' in k):.2e}  max_dI={err['max_dI']:.2e}  I_spread={spreads['I_full']}")
    return err, spreads


# =====================================================================
#  2D scan runner (x-V map)
# =====================================================================

def run_2d_scan(out_dir, tag, pSites, rots, cs, order, params_pme, hub_params, multipoleCoefs,
                xs, Vs, tip_base_pos, zTip, T_K, GammaS, GammaT,
                Wij_dense, real_sites, all_params_dict, active_sites=None):
    """Run a 2D x-V parity scan and produce imshow PNG + CSV."""
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()

    nSite = pSites.shape[0]
    nStates_full = 1 << nSite
    nx, nv = len(xs), len(Vs)
    T_eV = T_K * kB_eV
    ns = real_sites
    if active_sites is None: active_sites = list(range(nSite))

    # Build flat tip arrays for full solver (row-major: V changes fastest)
    nTips = nx * nv
    tip_pos_flat = np.zeros((nTips, 3), dtype=np.float32)
    Vtips_flat   = np.zeros(nTips, dtype=np.float32)
    for ix, x in enumerate(xs):
        for iv, v in enumerate(Vs):
            idx = ix * nv + iv
            tip_pos_flat[idx] = [x, tip_base_pos[1], zTip]
            Vtips_flat[idx] = v

    # --- Full solver ---
    pauli = PauliSolverCL(nSingle=nSite)
    pauli.set_lead(0, mu=0.0, temp=T_eV)
    pauli.set_lead(1, mu=0.0, temp=T_eV)

    curr_full, Es_full, Ts_full, P_full, StateEs_full, _, _ = pauli.scan_current_tip(
        pTips=tip_pos_flat, Vtips=Vtips_flat, pSites=pSites, params=params_pme,
        order=order, cs=cs, rots=rots,
        return_probs=True, return_state_energies=True, return_curmat=False, Wij=Wij_dense,
    )
    Es_full = Es_full.reshape(nTips, nSite)
    Ts_full = Ts_full.reshape(nTips, nSite)
    P_full  = P_full.reshape(nTips, nStates_full)
    StateEs_full = StateEs_full.reshape(nTips, nStates_full)

    # --- Hubbard precalc + PME (loop over tips for varying mu_tip) ---
    hub = HubbardSolver(build_options=["-DDBG_PRECALC=0"])
    pTips4 = np.zeros((nTips, 4), dtype=np.float32)
    pTips4[:, :3] = tip_pos_flat
    pTips4[:, 3]  = Vtips_flat

    Esite_hub, Tsite_hub = hub.precalc_esite_thop(
        posE=pSites, pTips=pTips4, rots=rots, multipoleCoefs=multipoleCoefs,
        bMirror=True, bRamp=True, params=hub_params,
    )

    W_sparse = make_sparse_W(Wij_dense, nMaxNeigh=max(1, nSite - 1))
    basis_u4, nStates_hub = make_basis_u4(1, nSite, active_sites=active_sites)
    gamma_sub = (GammaS / np.pi) ** 2
    gamma_tip = (GammaT / np.pi) ** 2

    curr_hub  = np.zeros(nTips, dtype=np.float32)

    for it in range(nTips):
        c, p, _ = hub.solve_pme_dense_fullbasis(
            W_sparse=W_sparse,
            Esite=Esite_hub[it:it+1].astype(np.float32),
            Tsite=Tsite_hub[it:it+1].astype(np.float32),
            basis_u4=basis_u4, nTips=1, nSite=nSite,
            mu_tip=float(Vtips_flat[it]), T_env=T_eV,
            Gamma_sub=gamma_sub, Gamma_tip=gamma_tip,
            bAlloc=(it == 0),
        )
        curr_hub[it] = c[0]

    # Reshape to 2D maps (nx, nv)
    I_full_2d = curr_full.reshape(nx, nv)
    I_hub_2d  = curr_hub.reshape(nx, nv)
    dI_2d     = (I_hub_2d - I_full_2d)

    elapsed = time.time() - t0
    print(f"  [{tag}] 2D scan {nx}x{nv} done in {elapsed:.1f}s  max|dI|={np.max(np.abs(dI_2d)):.2e}")

    # --- Save CSV of current maps ---
    np.savetxt(os.path.join(out_dir, f"I_full_{tag}.csv"), I_full_2d, delimiter=",")
    np.savetxt(os.path.join(out_dir, f"I_hub_{tag}.csv"), I_hub_2d, delimiter=",")

    # --- Summary ---
    meta = dict(all_params_dict)
    meta.update({"tag": tag, "scan": "2D x-V", "x_min": float(xs[0]), "x_max": float(xs[-1]),
                 "V_min": float(Vs[0]), "V_max": float(Vs[-1]), "nx": nx, "nV": nv,
                 "nSite": nSite, "real_sites": [int(x) for x in ns], "elapsed_sec": round(elapsed, 2)})
    with open(os.path.join(out_dir, f"summary_{tag}.txt"), "w") as f:
        f.write(f"# 2D Parity test: {tag}\n\n## ALL Parameters\n")
        f.write(json.dumps(meta, indent=2))
        f.write(f"\n\n## Spreads\n  I_full: {spread(curr_full)}\n  I_hub: {spread(curr_hub)}\n  max|dI|: {np.max(np.abs(dI_2d)):.8g}\n")

    # --- Plot ---
    ext = [float(Vs[0]), float(Vs[-1]), float(xs[0]), float(xs[-1])]
    vmax = max(np.max(np.abs(I_full_2d)), np.max(np.abs(I_hub_2d)), 1e-15)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    im0 = axs[0].imshow(I_full_2d, extent=ext, origin='lower', aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axs[0].set_title("Full PME current"); axs[0].set_xlabel("V"); axs[0].set_ylabel("x"); fig.colorbar(im0, ax=axs[0])
    im1 = axs[1].imshow(I_hub_2d, extent=ext, origin='lower', aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axs[1].set_title("Hubbard PME current"); axs[1].set_xlabel("V"); fig.colorbar(im1, ax=axs[1])
    dvmax = max(np.max(np.abs(dI_2d)), 1e-15)
    im2 = axs[2].imshow(dI_2d, extent=ext, origin='lower', aspect='auto', cmap='RdBu_r', vmin=-dvmax, vmax=dvmax)
    axs[2].set_title("Difference (hub-full)"); axs[2].set_xlabel("V"); fig.colorbar(im2, ax=axs[2])
    fig.suptitle(f"{tag}: 2D current map (x vs V)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"plot2d_{tag}.png"), dpi=150)
    plt.close(fig)


# =====================================================================
#  Main
# =====================================================================

def main():
    np.set_printoptions(linewidth=np.inf, precision=6, suppress=True)

    # === Parameters from Ruslan's JSON configs ===
    Rtip    = 3.0
    zV0     = -1.0
    zVd     = 15.0
    decay   = 0.3     # beta / Thop_decay
    GammaS  = 0.01
    GammaT  = 0.01
    W       = 0.05    # inter-site Coulomb
    T_K     = 3.0
    zTip    = 5.0
    Q0      = 1.0     # monopole
    Qzz     = 10.0    # zz quadrupole
    bMirror = 1.0
    bRamp   = 1.0
    Esite_ref = -0.09  # reference on-site (overridden by geometry file per-site)

    # Multipole coefficients for PME.cl: [Q0, px, py, pz, Qxx, Qyy, Qzz, Qyz, Qzx, Qxy]
    cs = np.zeros(10, dtype=np.float32)
    cs[0] = Q0; cs[6] = Qzz
    order = 2  # include quadrupole

    # PME params: [Rtip, zV0, zVd, Esite_ref, beta, Gamma, W, bMirror, bRamp]
    # Gamma here feeds (Gamma/pi)^2 in pauli_ocl.py
    params_pme = np.array([Rtip, zV0, zVd, Esite_ref, decay, GammaS, W, bMirror, bRamp], dtype=np.float32)

    # Hubbard precalc params
    hub_params = {"Rtip": Rtip, "zMirror": zV0, "zRampOffset": zVd, "Thop_decay": decay, "Thop_amp": 1.0}

    # Dict of ALL params for reporting
    all_params = {
        "Rtip": Rtip, "zV0": zV0, "zVd": zVd, "decay": decay, "GammaS": GammaS, "GammaT": GammaT,
        "W": W, "T_K": T_K, "T_eV": T_K * kB_eV, "zTip": zTip, "Q0": Q0, "Qzz": Qzz,
        "bMirror": bool(bMirror), "bRamp": bool(bRamp), "Esite_ref": Esite_ref,
        "multipole_order": order, "multipole_cs": [float(x) for x in cs],
        "params_pme_vec": [float(x) for x in params_pme],
        "hub_params": hub_params,
    }

    _ocl_dir = os.path.dirname(os.path.abspath(__file__))

    # ==================================================================
    # 2-SITE system (from Ruslan_long.txt, embedded in 4-site for PME.cl)
    # ==================================================================
    geom2_path = os.path.join(_ocl_dir, "Ruslan_long.txt")
    sites2_raw = load_geometry(geom2_path)  # (2, 4)
    pSites2 = embed_2site_in_4(sites2_raw)
    nSite2 = 4
    real2 = [0, 1]
    rots2 = np.tile(np.eye(3, dtype=np.float32), (nSite2, 1, 1))
    Wij2 = make_uniform_Wij(nSite2, W)
    Wij2[2, :] = 0; Wij2[:, 2] = 0; Wij2[3, :] = 0; Wij2[:, 3] = 0  # no W for spectators

    # Per-site multipole coefficients (nSite, nMulti)
    multipoleCoefs2 = np.zeros((nSite2, 10), dtype=np.float32)
    multipoleCoefs2[:, 0] = Q0
    multipoleCoefs2[:, 6] = Qzz

    VBias2 = 1.1  # from JSON
    tip_base2 = [0.0, 0.0, zTip]
    out2 = os.path.join(_repo_root, "figs", "parity_2site_ruslan")
    all_params_2 = dict(all_params, geometry=geom2_path, VBias_fixed=VBias2, nSite_real=2, nSite_embed=4)

    print("=== 2-SITE (Ruslan_long.txt, embedded in 4) ===")
    xs2 = np.linspace(-10.0, 10.0, 201)
    Vs2 = np.linspace(0.0, 1.5, 151)

    run_1d_scan(os.path.join(out2, "xscan"), "2site_xscan", pSites2, rots2, cs, order, params_pme, hub_params, multipoleCoefs2,
                "x", xs2, tip_base2, zTip, VBias2, T_K, GammaS, GammaT, Wij2, real2, all_params_2, active_sites=real2)

    run_1d_scan(os.path.join(out2, "Vscan"), "2site_Vscan", pSites2, rots2, cs, order, params_pme, hub_params, multipoleCoefs2,
                "V", Vs2, tip_base2, zTip, VBias2, T_K, GammaS, GammaT, Wij2, real2, all_params_2, active_sites=real2)

    xs2d = np.linspace(-10.0, 10.0, 101)
    Vs2d = np.linspace(0.0, 1.5, 76)
    run_2d_scan(os.path.join(out2, "xVmap"), "2site_xV", pSites2, rots2, cs, order, params_pme, hub_params, multipoleCoefs2,
                xs2d, Vs2d, tip_base2, zTip, T_K, GammaS, GammaT, Wij2, real2, all_params_2, active_sites=real2)

    # ==================================================================
    # 4-SITE system (from Ruslan_kite.txt)
    # ==================================================================
    geom4_path = os.path.join(_ocl_dir, "Ruslan_kite.txt")
    pSites4 = load_geometry(geom4_path)  # (4, 4) directly
    nSite4 = 4
    real4 = [0, 1, 2, 3]
    rots4 = np.tile(np.eye(3, dtype=np.float32), (nSite4, 1, 1))
    Wij4 = make_uniform_Wij(nSite4, W)

    multipoleCoefs4 = np.zeros((nSite4, 10), dtype=np.float32)
    multipoleCoefs4[:, 0] = Q0
    multipoleCoefs4[:, 6] = Qzz

    VBias4 = 1.0  # from JSON
    tip_base4 = [0.0, 0.0, zTip]
    out4 = os.path.join(_repo_root, "figs", "parity_4site_ruslan")
    all_params_4 = dict(all_params, geometry=geom4_path, VBias_fixed=VBias4, nSite_real=4)

    print("\n=== 4-SITE (Ruslan_kite.txt) ===")
    xs4 = np.linspace(-10.0, 10.0, 201)
    Vs4 = np.linspace(0.0, 1.5, 151)

    run_1d_scan(os.path.join(out4, "xscan"), "4site_xscan", pSites4, rots4, cs, order, params_pme, hub_params, multipoleCoefs4,
                "x", xs4, tip_base4, zTip, VBias4, T_K, GammaS, GammaT, Wij4, real4, all_params_4, active_sites=real4)

    run_1d_scan(os.path.join(out4, "Vscan"), "4site_Vscan", pSites4, rots4, cs, order, params_pme, hub_params, multipoleCoefs4,
                "V", Vs4, tip_base4, zTip, VBias4, T_K, GammaS, GammaT, Wij4, real4, all_params_4, active_sites=real4)

    xs4_2d = np.linspace(-10.0, 10.0, 101)
    Vs4_2d = np.linspace(0.0, 1.5, 76)
    run_2d_scan(os.path.join(out4, "xVmap"), "4site_xV", pSites4, rots4, cs, order, params_pme, hub_params, multipoleCoefs4,
                xs4_2d, Vs4_2d, tip_base4, zTip, T_K, GammaS, GammaT, Wij4, real4, all_params_4, active_sites=real4)

    print(f"\n=== ALL DONE ===")
    print(f"  2-site results: {out2}")
    print(f"  4-site results: {out4}")


if __name__ == "__main__":
    main()
