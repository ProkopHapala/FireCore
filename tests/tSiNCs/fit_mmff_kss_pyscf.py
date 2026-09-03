#!/usr/bin/env python3
"""Redo MMFF k-fit at frozen QM geometry (bKs + apars[:,2] Kss only), or own-min stretch RMSE vs PySCF PBE L1 NCs.

Molecules: CH4, C2H6, SiH4, Si2H6, adamantane C10H16, sila-adamantane Si10H16.
Frozen-geom path: fails if apars[:,0:2] (θ0) move. Does not FIRE. One molecule per process.
NC path (`--mol pyscf_nc`): FIRE then Hessian with the same switches (own-min). Canonical PBE bank: pyscf_vib_results (tight Si). Not FTIR: cube_C, octahedron_7ring_Si (still 1 imag; mode-follow pending).

Report (C cube vs octa; Morse + split angles; SA did not unify):
  doc/Topics/FTIR_Nanocrystals/MMFF_C_CH_vs_CH2_kfit.md

Run from tests/tMMFF:
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol all_si
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol adamantane
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol SiH4 --ref dftb_pbc-0-3
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref joint
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref nbscan
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref morse2d
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref anneal
  python3 ../tSiNCs/fit_mmff_kss_pyscf.py --mol pyscf_nc --ref octahedron_C
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIRECORE = Path("/home/prokop/git/FireCore")
CCU = Path("/home/prokop/git/CompChemUtils")
REF = CCU / "examples/tSiNCs/results"
CAGES = CCU / "examples/tSiNCs/jobs_dftb_vib_cages"
OUT = FIRECORE / "tests/tSiNCs/OUT_mmff_kss_fit"
FCONV = 1e-4
NSTEP = 8000
sys.path.insert(0, str(FIRECORE))
sys.path.insert(0, str(CCU / "examples/tSiNCs"))
from pyBall.FFfit_utils import PYSCF_VIB_L1_ROOT, resolve_pyscf_vib_case_dir, list_pyscf_vib_cases
PYSCF_NC = PYSCF_VIB_L1_ROOT

MASS = {"H": 1.00784, "C": 12.011, "N": 14.007, "O": 15.999, "Si": 28.085}


def load_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    sym, pos = [], []
    for ln in lines[2:2 + n]:
        t = ln.split()
        sym.append(t[0])
        pos.append([float(t[1]), float(t[2]), float(t[3])])
    return np.array(sym), np.asarray(pos, dtype=float)


def load_ref_freqs(mol, tag):
    meta_path = REF / mol / tag / "modes_meta.json"
    if not meta_path.is_file():
        raise FileNotFoundError(meta_path)
    meta = json.loads(meta_path.read_text())
    return np.asarray(meta["freqs_cm1"], dtype=float)


def hess_input_order(sess, hess):
    inv = sess._perm_in2int
    n = sess.n_atoms
    dof = np.empty(n * 3, dtype=np.int32)
    for ia in range(n):
        ii = int(inv[ia])
        dof[3 * ia:3 * ia + 3] = [3 * ii, 3 * ii + 1, 3 * ii + 2]
    return hess[dof][:, dof]


def mmff_signed(sess, masses):
    from pyBall.FFfit_utils import signed_frequencies_cm1
    sess.restore_input_positions()
    H = hess_input_order(sess, sess.get_hessian(dx=1e-4))
    f = signed_frequencies_cm1(H, masses)
    if not np.all(np.isfinite(f)):
        raise RuntimeError("non-finite MMFF frequencies")
    return f


def mmff_vib(sess, masses, n_ref):
    """FFfit at frozen QM q: take n_ref highest positive freqs (drop ~6 rigid)."""
    f = mmff_signed(sess, masses)
    pos = np.sort(f[f > 0.0])
    n_imag = int(np.sum(f < -10.0))
    if pos.size < n_ref:
        raise RuntimeError(f"only {pos.size} positive freqs, need {n_ref} (n_imag(|ν|>10)={n_imag})")
    return pos[-n_ref:], n_imag, pos[:6]


def rmse_sorted(calc, ref, fmin=10.0):
    a = np.sort(np.asarray(calc, dtype=float))
    b = np.sort(np.asarray(ref, dtype=float))
    a = a[a > fmin]
    b = b[b > fmin]
    if a.size != b.size or a.size == 0:
        return np.inf
    return float(np.sqrt(np.mean((a - b) ** 2)))


def heatmap(err, xs, ys, xlabel, ylabel, title, mark, path):
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    finite = err[np.isfinite(err)]
    plot = err.copy()
    if finite.size:
        plot[~np.isfinite(plot)] = finite.max() * 1.2
    im = ax.imshow(plot, origin="lower", aspect="auto", cmap="viridis_r",
                   extent=[xs[0], xs[-1], ys[0], ys[-1]])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="RMSE (cm⁻¹)")
    if mark is not None:
        ax.plot(mark[0], mark[1], "wo", ms=9, mec="red", mew=1.6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  saved {path}")


def overlay_stick(ref, calc, path, title, ref_label="reference"):
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.vlines(ref, 0.0, 1.0, colors="black", lw=1.6, label=ref_label)
    ax.vlines(calc, 0.0, 0.7, colors="C0", lw=1.2, label="MMFF k-fit")
    ax.set_xlabel("ν (cm⁻¹)")
    ax.set_yticks([])
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  saved {path}")


def real_freq_list(freq):
    f = np.asarray(freq, dtype=float).ravel()
    if np.iscomplexobj(f):
        f = np.where(np.abs(f.imag) > 1e-8, -np.abs(f.imag), f.real)
    return np.sort(f[np.isfinite(f) & (f > 10.0)])


def load_cage_freqs(job):
    p = CAGES / job / "frequencies_cm1.npy"
    if not p.is_file():
        raise FileNotFoundError(p)
    return real_freq_list(np.load(p))


def load_cage_xyz(job):
    p = CAGES / job / "relaxed.xyz"
    if not p.is_file():
        raise FileNotFoundError(p)
    return load_xyz(p)


def outdir_for(mol, ref_tag, pyscf_tag="pyscf_b3lyp_cc-pVDZ"):
    return OUT / mol if ref_tag == pyscf_tag else OUT / f"{mol}_{ref_tag}"


def print_apars(sess, tag):
    ap = np.asarray(sess.MMFF.apars)
    nnode = int(sess.MMFF.nnode)
    cs = np.hypot(ap[:nnode, 0], ap[:nnode, 1])
    print(f"  {tag} nnode={nnode}  Kss={np.unique(np.round(ap[:nnode, 2], 4))}  hypot(c,s)={np.unique(np.round(cs, 6))}")


def stretch_rmse(calc, ref, nH):
    a = np.sort(np.asarray(calc, dtype=float))
    b = np.sort(np.asarray(ref, dtype=float))
    if a.size < nH or b.size < nH:
        return None
    return float(np.sqrt(np.mean((a[-nH:] - b[-nH:]) ** 2)))


def print_freq_table(ref, calc, nmax=20):
    r = np.sort(ref); c = np.sort(calc)
    n = min(len(r), len(c), nmax)
    print(f"  {'mode':<6} {'ref':<10} {'MMFF':<10} {'Δ':<10}")
    for k in range(n):
        print(f"  {k+1:<6} {r[k]:<10.1f} {c[k]:<10.1f} {c[k]-r[k]:<+10.1f}")
    if len(r) > n:
        print(f"  … {len(r)-n} more  νmax ref/MMFF {r[-1]:.1f}/{c[-1] if len(c) else float('nan'):.1f}")


def fit_xh4(sess_cls, mol, ref_tag, n=17):
    xyz = REF / mol / ref_tag / "relaxed.xyz"
    if not xyz.is_file():
        raise FileNotFoundError(xyz)
    sym, pos = load_xyz(xyz)
    masses = np.array([MASS[s] for s in sym], dtype=float)
    ref = np.sort(load_ref_freqs(mol, ref_tag))
    sb = np.linspace(0.8, 2.6, n)
    sa = np.linspace(0.2, 1.8, n)
    err = np.full((n, n), np.inf)
    best = {"rmse": np.inf}
    print(f"\n======== {mol}  vs {ref_tag}  frozen {xyz}  {n}×{n}")
    print(f"  ref n={ref.size}  {np.array2string(ref, precision=1)}")
    sess = sess_cls(pos, list(sym), firecore_path=str(FIRECORE), enable_angles=True)
    print_apars(sess, "init")
    f0, n_imag0, rigid0 = mmff_vib(sess, masses, ref.size)
    print(f"  default MMFF RMSE={rmse_sorted(f0, ref):.1f}  n_imag={n_imag0}  6-lowest+={np.array2string(rigid0, precision=1)}")
    print(f"  default ν={np.array2string(f0, precision=1)}")
    t0 = time.perf_counter()
    for i, b in enumerate(sb):
        for j, a in enumerate(sa):
            sess.set_scales(b, a)
            f, n_imag, _ = mmff_vib(sess, masses, ref.size)
            e = rmse_sorted(f, ref)
            err[i, j] = e
            if e < best["rmse"]:
                best = {"rmse": e, "scale_bond": float(b), "scale_angle": float(a), "freqs": f, "n_imag": n_imag}
                print(f"  [{i*n+j+1}/{n*n}] XH={b:.3f} angK={a:.3f} RMSE={e:.2f} n_imag={n_imag} ***")
    print(f"  grid {time.perf_counter()-t0:.1f}s")
    out = outdir_for(mol, ref_tag)
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "grid.npz", bond_scales=sb, angle_scales=sa, errors=err,
             best_scale_bond=best["scale_bond"], best_scale_angle=best["scale_angle"], best_rmse=best["rmse"],
             freqs_ref=ref, freqs_best=best["freqs"], freqs_default=f0, protocol="k_only_frozen_Kss_col2")
    heatmap(err, sa, sb, "angle Kss scale (apars[:,2])", f"{mol} X–H bond k scale (bKs)",
            f"{mol} vs {ref_tag}  best k={best['scale_bond']:.3f} K∠={best['scale_angle']:.3f}  RMSE={best['rmse']:.1f} cm⁻¹",
            (best["scale_angle"], best["scale_bond"]), out / "heatmap.png")
    overlay_stick(ref, best["freqs"], out / "overlay.png",
                  f"{mol} vs {ref_tag}  k={best['scale_bond']:.3f}  K∠={best['scale_angle']:.3f}  RMSE={best['rmse']:.1f} cm⁻¹",
                  ref_label=ref_tag)
    print(f"  BEST {mol}  k_XH={best['scale_bond']:.4f}  K_ang={best['scale_angle']:.4f}  RMSE={best['rmse']:.2f} cm-1")
    print_freq_table(ref, best["freqs"])
    rec = {"mol": mol, "scale_XH": best["scale_bond"], "scale_angle_K": best["scale_angle"],
           "rmse_all": best["rmse"], "n_imag_best": best["n_imag"], "rmse_default": rmse_sorted(f0, ref), "n": n,
           "reference": f"{ref_tag} frozen relaxed.xyz", "freqs_ref": ref.tolist(),
           "freqs_best": np.sort(best["freqs"]).tolist(), "freqs_default": np.sort(f0).tolist()}
    rec["scale_CH"] = rec["scale_XH"]  # alias for CH4 callers
    rec["rmse_all9"] = rec["rmse_all"]
    (out / "best.json").write_text(json.dumps(rec, indent=2))
    return best


def fit_x2h6(sess_cls, mol, ref_tag, scale_xh, n=15, fmin=500.0):
    xyz = REF / mol / ref_tag / "relaxed.xyz"
    if not xyz.is_file():
        raise FileNotFoundError(xyz)
    sym, pos = load_xyz(xyz)
    masses = np.array([MASS[s] for s in sym], dtype=float)
    ref = np.sort(load_ref_freqs(mol, ref_tag))
    ref_hi = ref[ref > fmin]
    sxx = np.linspace(0.5, 2.6, n)
    sa = np.linspace(0.2, 1.8, n)
    err = np.full((n, n), np.inf)
    best = {"rmse": np.inf}
    print(f"\n======== {mol}  vs {ref_tag}  freeze geom  X–H k fixed={scale_xh:.4f}  {n}×{n}  fmin={fmin}")
    print(f"  ref n={ref.size}  n(>{fmin:g})={ref_hi.size}  {np.array2string(ref, precision=1)}")
    sess = sess_cls(pos, list(sym), firecore_path=str(FIRECORE), enable_angles=True)
    print_apars(sess, "init")
    sess.set_scales_per_bond_type(scale_ch=scale_xh, scale_cc=1.0, scale_angle=1.0)
    f0_all, n_imag0, rigid0 = mmff_vib(sess, masses, ref.size)
    e0 = rmse_sorted(np.sort(f0_all)[-ref_hi.size:], ref_hi, fmin=0.0) if ref_hi.size else np.inf
    print(f"  MMFF (XH×{scale_xh:.3f}, XX=1, K∠=1) RMSE_tail={e0:.1f}  n_imag={n_imag0}  n={f0_all.size}")
    print(f"  6-lowest+={np.array2string(rigid0, precision=1)}")
    t0 = time.perf_counter()
    for i, xx in enumerate(sxx):
        for j, a in enumerate(sa):
            sess.set_scales_per_bond_type(scale_ch=scale_xh, scale_cc=xx, scale_angle=a)
            f_all, n_imag, _ = mmff_vib(sess, masses, ref.size)
            f = f_all  # match n_ref_hi highest vs ref_hi (MMFF torsion may sit above fmin)
            e = rmse_sorted(np.sort(f)[-ref_hi.size:], ref_hi, fmin=0.0) if ref_hi.size else np.inf
            err[i, j] = e
            if e < best["rmse"]:
                best = {"rmse": e, "scale_xx": float(xx), "scale_angle": float(a), "freqs": f_all, "n_imag": n_imag}
                print(f"  [{i*n+j+1}/{n*n}] XX={xx:.3f} angK={a:.3f} RMSE>{fmin:g}={e:.2f} n_imag={n_imag} ***")
    print(f"  grid {time.perf_counter()-t0:.1f}s")
    if not np.isfinite(best.get("rmse", np.inf)):
        raise RuntimeError(f"{mol}: no finite RMSE. MMFF ν={np.array2string(f0_all, precision=1)} ref_hi n={ref_hi.size}")
    out = outdir_for(mol, ref_tag)
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "grid.npz", xx_scales=sxx, angle_scales=sa, errors=err, scale_xh=scale_xh, fmin=fmin,
             best_scale_xx=best["scale_xx"], best_scale_angle=best["scale_angle"], best_rmse=best["rmse"],
             freqs_ref=ref, freqs_best=best["freqs"], protocol="k_only_frozen_Kss_col2_drop_torsion")
    heatmap(err, sa, sxx, "angle Kss scale (apars[:,2])", f"{mol} X–X bond k scale (bKs)",
            f"{mol} vs {ref_tag} (ν>{fmin:g})  XH×{scale_xh:.3f}  XX={best['scale_xx']:.3f} K∠={best['scale_angle']:.3f}  RMSE={best['rmse']:.1f}",
            (best["scale_angle"], best["scale_xx"]), out / "heatmap.png")
    overlay_stick(ref_hi, np.sort(best["freqs"])[-ref_hi.size:], out / "overlay.png",
                  f"{mol} ν>{fmin:g}  XH×{scale_xh:.3f} XX={best['scale_xx']:.3f} K∠={best['scale_angle']:.3f}",
                  ref_label=ref_tag)
    print(f"  BEST {mol}  k_XX={best['scale_xx']:.4f}  K_ang={best['scale_angle']:.4f}  RMSE(ν>{fmin:g})={best['rmse']:.2f}")
    rec = {"mol": mol, "scale_XH_fixed": scale_xh, "scale_XX": best["scale_xx"], "scale_angle_K": best["scale_angle"],
           "rmse_gt_fmin": best["rmse"], "fmin": fmin, "n_imag_best": best["n_imag"], "n": n,
           "reference": f"{ref_tag} frozen; RMSE drops modes <{fmin:g} cm-1",
           "freqs_ref": ref.tolist(), "freqs_best": np.sort(best["freqs"]).tolist(),
           "scale_CH_fixed": scale_xh, "scale_CC": best["scale_xx"], "rmse_gt500": best["rmse"]}
    (out / "best.json").write_text(json.dumps(rec, indent=2))
    return best


def fit_cage(sess_cls, name, job, n=11, mode="xh", k_xx=1.0, k_xh=1.0, sess=None):
    """2D k-grid on a 26-atom cage. mode='xh': k_XH×K with k_XX fixed; mode='xx': k_XX×K with k_XH fixed.
    Pass sess= to reuse the same MMFF init (cannot re-init in one process)."""
    xyz_path = CAGES / job / "relaxed.xyz"
    if sess is None:
        sym, pos = load_xyz(xyz_path)
        masses = np.array([MASS[s] for s in sym], dtype=float)
        sess = sess_cls(pos, list(sym), firecore_path=str(FIRECORE), enable_angles=True)
        print_apars(sess, "init")
    else:
        masses = np.array([MASS[s] for s in sess.symbols], dtype=float)
        nH = int(sum(s == "H" for s in sess.symbols))
    ref = load_cage_freqs(job)
    nH = int(sum(s == "H" for s in sess.symbols))
    sa = np.linspace(0.2, 1.8, n)
    sk = np.linspace(0.8, 3.6, n) if mode == "xx" else np.linspace(0.8, 2.4, n)
    err = np.full((n, n), np.inf)
    best = {"rmse": np.inf}
    print(f"\n======== {name}  vs {job}  mode={mode}  frozen {xyz_path.name}  {n}×{n}  nH={nH}")
    print(f"  ref n={ref.size}  ν={ref.min():.1f}…{ref.max():.1f}  n_stretch_highest={nH}")
    sess.set_scales_per_bond_type(scale_ch=1.0, scale_cc=1.0, scale_angle=1.0)
    f0, n_imag0, rigid0 = mmff_vib(sess, masses, ref.size)
    print(f"  default RMSE_all={rmse_sorted(f0, ref):.1f}  RMSE_stretch={stretch_rmse(f0, ref, nH)}  n_imag={n_imag0}")
    print(f"  6-lowest+={np.array2string(rigid0, precision=1)}  νmax={f0.max():.1f}")
    t0 = time.perf_counter()
    for i, kv in enumerate(sk):
        for j, a in enumerate(sa):
            if mode == "xh":
                sess.set_scales_per_bond_type(scale_ch=kv, scale_cc=k_xx, scale_angle=a)
            elif mode == "xx":
                sess.set_scales_per_bond_type(scale_ch=k_xh, scale_cc=kv, scale_angle=a)
            else:
                raise ValueError(mode)
            f, n_imag, _ = mmff_vib(sess, masses, ref.size)
            e = rmse_sorted(f, ref)
            err[i, j] = e
            if e < best["rmse"]:
                best = {"rmse": e, "scale_k": float(kv), "scale_angle": float(a), "freqs": f, "n_imag": n_imag,
                        "rmse_stretch": stretch_rmse(f, ref, nH)}
                print(f"  [{i*n+j+1}/{n*n}] k={kv:.3f} angK={a:.3f} RMSE_all={e:.2f} RMSE_st={best['rmse_stretch']} n_imag={n_imag} ***")
    print(f"  grid {time.perf_counter()-t0:.1f}s")
    tag = f"{name}_{mode}"
    out = OUT / tag
    out.mkdir(parents=True, exist_ok=True)
    ylab = f"{name} X–H k" if mode == "xh" else f"{name} X–X k"
    heatmap(err, sa, sk, "angle Kss scale (apars[:,2])", ylab,
            f"{name} {mode} vs {job}  k={best['scale_k']:.3f} K∠={best['scale_angle']:.3f}  RMSE={best['rmse']:.1f}",
            (best["scale_angle"], best["scale_k"]), out / "heatmap.png")
    overlay_stick(ref, best["freqs"], out / "overlay.png",
                  f"{name} {mode}  k={best['scale_k']:.3f} K∠={best['scale_angle']:.3f}  RMSE={best['rmse']:.1f}",
                  ref_label=job)
    print(f"  BEST {name}/{mode}  k={best['scale_k']:.4f}  K_ang={best['scale_angle']:.4f}  RMSE_all={best['rmse']:.2f}  RMSE_st={best['rmse_stretch']}")
    rec = {"name": name, "job": job, "mode": mode, "scale_k": best["scale_k"], "scale_angle_K": best["scale_angle"],
           "k_XH_fixed": k_xh if mode == "xx" else best["scale_k"], "k_XX_fixed": k_xx if mode == "xh" else best["scale_k"],
           "rmse_all": best["rmse"], "rmse_stretch_nH": best["rmse_stretch"], "n_imag_best": best["n_imag"],
           "rmse_default": rmse_sorted(f0, ref), "nH": nH, "n": n,
           "nu_max_ref": float(ref.max()), "nu_max_best": float(np.max(best["freqs"])),
           "freqs_ref": ref.tolist(), "freqs_best": np.sort(best["freqs"]).tolist()}
    np.savez(out / "grid.npz", k_scales=sk, angle_scales=sa, errors=err, **{k: rec[k] for k in rec if k not in ("freqs_ref", "freqs_best")})
    (out / "best.json").write_text(json.dumps(rec, indent=2))
    return best, sess


def xh4_best_path(mol, ref_tag):
    p = outdir_for(mol, ref_tag) / "best.json"
    if not p.is_file():
        raise FileNotFoundError(f"{p} — run {mol} vs {ref_tag} first")
    return json.loads(p.read_text())


def load_pyscf_nc(case):
    d = resolve_pyscf_vib_case_dir(PYSCF_NC, case)
    if not (d / "relaxed.xyz").is_file():
        raise FileNotFoundError(d)
    from pyBall.FFfit_utils import as_signed_wavenumbers_cm1
    sym, pos = load_xyz(d / "relaxed.xyz")
    w = as_signed_wavenumbers_cm1(np.load(d / "frequencies_cm1.npy"))
    n_imag = int(np.sum(w < -10.0))
    ref = np.sort(w[w > 10.0])
    nH = int(sum(s == "H" for s in sym))
    masses = np.load(d / "masses.npy").astype(float)
    if masses.shape[0] != len(sym):
        raise ValueError(f"{case}: masses {masses.shape} vs N={len(sym)}")
    st = json.loads((d / "status.json").read_text()) if (d / "status.json").is_file() else {}
    has_si = any(s == "Si" for s in list(sym))
    return {"case": case, "dir": d, "sym": [str(s) for s in list(sym)], "pos": pos, "masses": masses, "nH": nH,
            "ref": ref, "omega_all": w, "n_imag_ref": n_imag, "at_min": n_imag == 0, "status": st, "has_si": has_si,
            "nu_max_ref": float(ref.max()) if ref.size else None,
            "mean_stretch_ref": float(ref[-nH:].mean()) if ref.size >= nH else None}


def pack_nc_eval(sess, rec, k_xh, k_xx, k_ang, tag):
    sess.set_scales_per_bond_type(scale_ch=k_xh, scale_cc=k_xx, scale_angle=k_ang)
    f = mmff_signed(sess, rec["masses"])
    n_imag = int(np.sum(f < -10.0))
    pos = np.sort(f[f > 0.0])
    nH = rec["nH"]
    rst = rec["ref"][-nH:] if rec["ref"].size >= nH else rec["ref"]
    st = pos[-nH:] if pos.size >= nH else pos
    return {"tag": tag, "k_XH": k_xh, "k_XX": k_xx, "Kss": k_ang, "n_imag_mmff": n_imag,
            "nu_max": float(pos.max()) if pos.size else None,
            "mean_stretch": float(st.mean()) if st.size else None,
            "rmse_stretch": stretch_rmse(pos, rec["ref"], nH),
            "n_pos": int(pos.size), "n_ref": int(rec["ref"].size)}


def fit_pyscf_nc(sess_cls, case, n=9):
    rec = load_pyscf_nc(case)
    print(f"\n======== PySCF-NC {case}  PBE/ccecp-cc-pVDZ  N={len(rec['sym'])} nH={rec['nH']}  at_min={rec['at_min']} n_imag_ref={rec['n_imag_ref']}")
    print(f"  ref n>10={rec['ref'].size}  ν={rec['ref'].min():.1f}…{rec['ref'].max():.1f}  <stretch nH>={rec['mean_stretch_ref']:.1f}")
    if not rec["at_min"]:
        print(f"  WARNING: PySCF Hessian is NOT a minimum ({rec['n_imag_ref']} imaginaries). FFfit only.")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    print_apars(sess, "init")
    e0 = pack_nc_eval(sess, rec, 1.0, 1.0, 1.0, "default")
    print(f"  default stretch RMSE={e0['rmse_stretch']}  <st>={e0['mean_stretch']:.1f} vs {rec['mean_stretch_ref']:.1f}  n_imag={e0['n_imag_mmff']}")
    if rec["has_si"]:
        sb = np.linspace(0.6, 1.6, n)
        sa = np.linspace(0.20, 1.2, n)
    else:
        sb = np.linspace(1.0, 2.4, n)
        sa = np.linspace(0.25, 1.2, n)
    err = np.full((n, n), np.inf)
    best = {"rmse": np.inf}
    t0 = time.perf_counter()
    for i, b in enumerate(sb):
        for j, a in enumerate(sa):
            pack = pack_nc_eval(sess, rec, b, 1.0, a, "grid")
            e = pack["rmse_stretch"]
            e = np.inf if e is None else e
            err[i, j] = e
            if e < best["rmse"]:
                best = {"rmse": e, "scale_XH": float(b), "scale_angle": float(a), "pack": pack}
                print(f"  [{i*n+j+1}/{n*n}] XH={b:.3f} K={a:.3f} RMSE_st={e:.2f} <st>={pack['mean_stretch']:.1f} ***")
    print(f"  grid {time.perf_counter()-t0:.1f}s")
    if not np.isfinite(best["rmse"]):
        raise RuntimeError(f"{case}: no finite stretch RMSE")
    kss_freeze = 0.40 if rec["has_si"] else 0.488
    i_k = int(np.argmin(np.min(err, axis=1)))
    kss_span = float(np.nanmax(err[i_k]) - np.nanmin(err[i_k]))
    print(f"  Kss stretch-RMSE span at best k_XH: {kss_span:.2f} cm-1 — freeze Kss={kss_freeze} (XH4 all-mode), 1D-refine k_XH")
    sb1 = np.linspace(max(sb[0], best["scale_XH"] - 0.25), min(sb[-1], best["scale_XH"] + 0.25), 21)
    best1 = {"rmse": np.inf}
    for b in sb1:
        pack = pack_nc_eval(sess, rec, float(b), 1.0, kss_freeze, "k_1d")
        e = pack["rmse_stretch"]
        e = np.inf if e is None else e
        if e < best1["rmse"]:
            best1 = {"rmse": e, "scale_XH": float(b), "scale_angle": kss_freeze, "pack": pack}
    print(f"  1D k_XH={best1['scale_XH']:.4f}  K={kss_freeze:.3f}  RMSE_stretch={best1['rmse']:.2f}  <st>={best1['pack']['mean_stretch']:.1f}")
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "grid.npz", bond_scales=sb, angle_scales=sa, errors=err,
             best_scale_XH=best1["scale_XH"], best_scale_angle=best1["scale_angle"], best_rmse=best1["rmse"],
             grid2d_scale_XH=best["scale_XH"], grid2d_scale_angle=best["scale_angle"], grid2d_rmse=best["rmse"])
    heatmap(err, sa, sb, "angle Kss scale (apars[:,2])", f"{case} X–H k",
            f"{case} vs PBE/cc-pVDZ  k={best1['scale_XH']:.3f} K∠={best1['scale_angle']:.3f}  RMSE_st={best1['rmse']:.1f}",
            (best1["scale_angle"], best1["scale_XH"]), out / "heatmap.png")
    sess.set_scales_per_bond_type(scale_ch=best1["scale_XH"], scale_cc=1.0, scale_angle=best1["scale_angle"])
    fbest = np.sort(mmff_signed(sess, rec["masses"]))
    fbest = fbest[fbest > 0]
    overlay_stick(rec["ref"][-rec["nH"]:], fbest[-rec["nH"]:], out / "overlay_stretch.png",
                  f"{case} stretch  k={best1['scale_XH']:.3f} K∠={best1['scale_angle']:.3f}  RMSE={best1['rmse']:.1f}",
                  ref_label="PySCF PBE")
    js = {"case": case, "reference": "PySCF PBE/ccecp-cc-pVDZ frozen relaxed.xyz",
          "protocol": "stretch RMSE (nH highest); Kss frozen at XH4 all-mode because stretch-only is degenerate",
          "at_min": rec["at_min"], "n_imag_ref": rec["n_imag_ref"], "nH": rec["nH"], "n": n,
          "scale_XH": best1["scale_XH"], "scale_XX": 1.0, "scale_angle_K": best1["scale_angle"],
          "rmse_stretch": best1["rmse"], "kss_stretch_span_cm1": kss_span,
          "grid2d_stretch_min": {"scale_XH": best["scale_XH"], "scale_angle_K": best["scale_angle"], "rmse_stretch": best["rmse"]},
          "default": e0, "best": best1["pack"], "mean_stretch_ref": rec["mean_stretch_ref"]}
    (out / "best.json").write_text(json.dumps(js, indent=2))
    print(f"  BEST {case}  k_XH={best1['scale_XH']:.4f}  K={best1['scale_angle']:.4f}  RMSE_stretch={best1['rmse']:.2f}")
    return best1, sess, rec


def eval_pyscf_nc_scales(sess_cls, case, scale_sets):
    rec = load_pyscf_nc(case)
    print(f"\n======== eval {case}  at_min={rec['at_min']} n_imag_ref={rec['n_imag_ref']}  <st_ref>={rec['mean_stretch_ref']:.1f}")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    rows = []
    for tag, kxh, kxx, ka in scale_sets:
        p = pack_nc_eval(sess, rec, kxh, kxx, ka, tag)
        rows.append(p)
        print(f"  {tag:20s} kXH={kxh:.3f} K={ka:.3f}  RMSE_st={p['rmse_stretch']}  <st>={p['mean_stretch']}  n_imag={p['n_imag_mmff']}")
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    js = {"case": case, "at_min": rec["at_min"], "n_imag_ref": rec["n_imag_ref"], "nH": rec["nH"],
          "mean_stretch_ref": rec["mean_stretch_ref"], "evals": rows}
    (out / "eval.json").write_text(json.dumps(js, indent=2))
    return js


def write_xyz(path, symbols, pos, comment):
    pos = np.asarray(pos, dtype=float)
    with open(path, "w") as f:
        f.write(f"{len(symbols)}\n{comment}\n")
        for s, p in zip(symbols, pos):
            f.write(f"{s}  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}\n")


def rg_of(pos):
    p = np.asarray(pos, dtype=float)
    return float(np.sqrt(np.mean(np.sum((p - p.mean(0)) ** 2, axis=1))))


def fire_then_hessian(sess, masses, k_xh, k_xx, k_ang, tag, pos0, n_kick=3, after_scales=None):
    """Relax with current k, Hessian at that minimum. Do not restore q after FIRE.

    If FIRE stops at |F|≈0 but the Hessian still has |ν|>10 imaginaries, the point is a
    saddle (odd-ring NCs do this). Kick 0.02 Å and re-FIRE. Fail loud if still a saddle.
    """
    from pyBall.FFfit_utils import signed_frequencies_cm1, assert_harmonic_spectrum_at_minimum
    sess.set_scales_per_bond_type(scale_ch=k_xh, scale_cc=k_xx, scale_angle=k_ang)
    if after_scales is not None:
        after_scales(sess)
    sess.restore_input_positions()
    rel = sess.relax(nstep_max=NSTEP, fconv=FCONV)
    pos = sess.positions_input_order()
    rg, rg0 = rg_of(pos), rg_of(pos0)
    print(f"  {tag} FIRE n={rel['n_steps']} fmax={rel['fmax']:.3e} Rg={rg0:.3f}→{rg:.3f} Å")
    if rg < 0.5 * rg0:
        raise RuntimeError(f"{tag} FIRE collapsed (Rg {rg0:.3f}→{rg:.3f} Å). k-scale must not move q0.")
    H = hess_input_order(sess, sess.get_hessian(dx=1e-4))
    if not np.all(np.isfinite(H)):
        raise RuntimeError(f"{tag}: non-finite Hessian")
    f = signed_frequencies_cm1(H, masses)
    for ik in range(n_kick):
        n_imag = int(np.sum(f < -10.0))
        if n_imag == 0:
            break
        print(f"  {tag} saddle: {n_imag} imag (min {float(f.min()):.1f} cm-1) — kick {ik+1}/{n_kick} and re-FIRE")
        rng = np.random.default_rng(1000 + ik)
        pos_in = sess.positions_input_order()
        kicked = pos_in + 0.02 * rng.normal(size=pos_in.shape)
        sess.MMFF.apos[:sess.n_atoms] = np.asarray(kicked, dtype=float)[sess._perm_int2in]
        sess.MMFF.set_opt(dt_max=0.05, dt_min=0.01, damp_max=0.1, cvf_min=-0.1, cvf_max=+0.1)
        rel = sess.relax(nstep_max=NSTEP, fconv=FCONV)
        pos = sess.positions_input_order()
        rg = rg_of(pos)
        print(f"  {tag} re-FIRE n={rel['n_steps']} fmax={rel['fmax']:.3e} Rg={rg:.3f} Å")
        if rg < 0.5 * rg0:
            raise RuntimeError(f"{tag} FIRE collapsed after kick (Rg {rg0:.3f}→{rg:.3f} Å)")
        H = hess_input_order(sess, sess.get_hessian(dx=1e-4))
        f = signed_frequencies_cm1(H, masses)
    st = assert_harmonic_spectrum_at_minimum(f, ctx=f"{tag}: ")
    return {"rel": rel, "pos": pos, "H": H, "freqs": f, "Rg": rg, "Rg0": rg0, "stats": st, "k_XH": k_xh, "Kss": k_ang}


def hydride_hh_stats(sym, pos, r_xh=1.30):
    """CH/CH2 counts and vicinal (non-geminal) H···H. Geminal CH2 is 1–3, skipped by Exclusion2."""
    sym = np.asarray(sym); pos = np.asarray(pos, dtype=float)
    iX = np.where((sym == "C") | (sym == "Si"))[0]; iH = np.where(sym == "H")[0]
    dXH = np.linalg.norm(pos[iX, None, :] - pos[iH][None, :, :], axis=2)
    nH = (dXH < r_xh).sum(axis=1)
    gem = np.zeros((len(iH), len(iH)), dtype=bool)
    for row in dXH:
        hs = np.where(row < r_xh)[0]
        for a, b in ((a, b) for a in range(len(hs)) for b in range(a + 1, len(hs))):
            gem[hs[a], hs[b]] = gem[hs[b], hs[a]] = True
    dHH = np.linalg.norm(pos[iH, None, :] - pos[iH][None, :, :], axis=2); np.fill_diagonal(dHH, np.inf)
    vic = dHH.copy(); vic[gem] = np.inf
    vmin = float(np.min(vic)) if np.isfinite(np.min(vic)) else None
    return {"nXH": int((nH == 1).sum()), "nXH2": int((nH == 2).sum()), "nXH3": int((nH == 3).sum()),
            "hh_min": float(np.min(dHH)), "hh_vic_min": vmin, "n_hh_vic_lt22": int(np.sum(np.triu(vic < 2.2, 1)))}


def general_c_pack():
    p = OUT / "pyscf_pbe_ncs" / "general_ff.json"
    if p.is_file():
        j = json.loads(p.read_text())["C"]
        return float(j["k_XH"]), float(j["Kss"]), float(j["k_XX"])
    return 1.775, 0.488, 1.0


def nbscan_case(sess_cls, case, scales=(0.0, 0.5, 1.0, 2.0, 4.0, 8.0)):
    """Own-min stretch RMSE vs H EvdW scale at frozen general-C k. s=0 is bonded-only (NB off). Cluster NB is clamped LJ, not Morse."""
    rec = load_pyscf_nc(case)
    if rec["has_si"]:
        raise RuntimeError(f"nbscan {case}: Si PBE Hessians are off-minimum; C only")
    k_xh, kss, k_xx = general_c_pack()
    print(f"\n======== nbscan {case}  k_XH={k_xh:.4f} Kss={kss:.3f} k_XX={k_xx:.3f}  at_min_ref={rec['at_min']}")
    hh0 = hydride_hh_stats(rec["sym"], rec["pos"])
    print(f"  DFT geom  XH/XH2/XH3={hh0['nXH']}/{hh0['nXH2']}/{hh0['nXH3']}  HHmin={hh0['hh_min']:.3f}  vic_min={hh0['hh_vic_min']}  n_vic<2.2={hh0['n_hh_vic_lt22']}  <st>_ref={rec['mean_stretch_ref']:.1f}")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    rows = []
    t0 = time.perf_counter()
    for s in scales:
        s = float(s)
        if s == 0.0:
            sess.disable_nonbond(); n_surf = n_bulk = 0; tag = "NB off"
        else:
            n_surf, n_bulk = sess.enable_surface_nonbond(scale_H_EvdW=s)
            tag = f"H-EvdW×{s:g}"
        pack = fire_then_hessian(sess, rec["masses"], k_xh, k_xx, kss, tag, rec["pos"])
        posf = np.sort(pack["freqs"][pack["freqs"] > 0.0])
        e = stretch_rmse(posf, rec["ref"], rec["nH"])
        e = np.inf if e is None else e
        hh = hydride_hh_stats(rec["sym"], pack["pos"])
        mean_st = float(posf[-rec["nH"]:].mean()) if posf.size >= rec["nH"] else None
        print(f"  {tag:16s} RMSE_st={e:.2f}  <st>={mean_st:.1f}  n_surf/bulk={n_surf}/{n_bulk}  HHmin={hh['hh_min']:.3f}  vic={hh['hh_vic_min']}  n_vic<2.2={hh['n_hh_vic_lt22']}")
        rows.append({"s_H_EvdW": s, "nb": s != 0.0, "n_surf": int(n_surf), "n_bulk": int(n_bulk),
                     "rmse_stretch": e, "mean_stretch": mean_st, "nu_max": float(pack["freqs"].max()),
                     "fmax": pack["rel"]["fmax"], "n_steps": pack["rel"]["n_steps"], **hh})
    print(f"  nbscan {time.perf_counter()-t0:.1f}s")
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "nbscan.npz", s_H_EvdW=np.array([r["s_H_EvdW"] for r in rows]),
             rmse_stretch=np.array([r["rmse_stretch"] for r in rows]),
             mean_stretch=np.array([r["mean_stretch"] for r in rows]), k_XH=k_xh, Kss=kss, k_XX=k_xx)
    js = {"case": case, "protocol": "own-min stretch RMSE vs H EvdW; surface LJ Exclusion2; k frozen at general C; s=0 NB off",
          "note": "cluster MMFF nonbond is clamped LJ (evalLJQs_ex2), not Morse/GridFF",
          "k_XH": k_xh, "Kss": kss, "k_XX": k_xx, "dft_geom": hh0, "mean_stretch_ref": rec["mean_stretch_ref"], "rows": rows}
    (out / "nbscan.json").write_text(json.dumps(js, indent=2))
    print(f"  saved {out / 'nbscan.json'}")


def reduce_nbscan(cases):
    root = OUT / "pyscf_pbe_ncs"
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for c in cases:
        z = np.load(root / c / "nbscan.npz")
        s, e = np.asarray(z["s_H_EvdW"], float), np.asarray(z["rmse_stretch"], float)
        ax.plot(s, e, "o-", ms=5, label=c)
        print(f"  {c}  s={s}  RMSE={np.round(e, 1)}")
    ax.set_xlabel("H EvdW scale  (0 = NB off)")
    ax.set_ylabel("own-min stretch RMSE (cm⁻¹)")
    k_xh, kss, k_xx = general_c_pack()
    ax.set_title(f"C surface LJ  k_XH={k_xh:.3f}  Kss={kss:.3f}  k_XX={k_xx:g}")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, None)
    fig.tight_layout()
    path = root / "nbscan_C.png"
    fig.savefig(path, dpi=160); plt.close(fig)
    print(f"  saved {path}")


def c_true_min_cases(cases=None):
    cases = cases or list_pyscf_nc_cases()
    return [c for c in cases if c.endswith("_C") and c != "cube_C"]


def after_morse(k_Morse, s_H=1.0):
    def _fn(sess):
        sess.enable_surface_nonbond(scale_H_EvdW=s_H)
        sess.set_morse_nonbond(k_Morse)
    return _fn


def after_split_morse(pack):
    def _fn(sess):
        sess.enable_each_angle()
        sess.set_angle_type_scales(k_XXX=float(pack["k_XXX"]), k_XXH=float(pack["k_XXH"]), k_HXH=float(pack["k_HXH"]))
        sess.enable_surface_nonbond(scale_H_EvdW=float(pack.get("s_H_EvdW", 1.0)))
        sess.set_morse_nonbond(float(pack["k_Morse"]))
    return _fn


def ownmin_stretch_of(sess, rec, k_xh, kss, tag, after=None):
    pack = fire_then_hessian(sess, rec["masses"], float(k_xh), 1.0, float(kss), tag, rec["pos"], after_scales=after)
    posf = np.sort(pack["freqs"][pack["freqs"] > 0.0])
    e = stretch_rmse(posf, rec["ref"], rec["nH"])
    if e is None:
        raise RuntimeError(f"{tag}: stretch RMSE is None")
    hh = hydride_hh_stats(rec["sym"], pack["pos"])
    mean_st = float(posf[-rec["nH"]:].mean()) if posf.size >= rec["nH"] else None
    print(f"  {tag} RMSE_st={e:.2f}  <st>={mean_st:.1f}  HHmin={hh['hh_min']:.3f}  vic={hh['hh_vic_min']}  n_vic<2.2={hh['n_hh_vic_lt22']}", flush=True)
    return {"rmse_stretch": e, "mean_stretch": mean_st, **hh, "nu_max": float(pack["freqs"].max()), "fmax": pack["rel"]["fmax"]}


def morse2d_case(sess_cls, case, n_xh=4, n_m=4):
    """Own-min stretch RMSE vs (k_XH, Morse α) at frozen Kss, surface Exclusion2 Morse."""
    rec = load_pyscf_nc(case)
    if rec["has_si"]:
        raise RuntimeError(f"morse2d {case}: C only (Si PBE is off-minimum)")
    k_xh0, kss, k_xx = general_c_pack()
    kxh = np.linspace(1.55, 1.95, n_xh)
    km = np.array([0.8, 1.5, 2.2, 3.0]) if n_m == 4 else np.linspace(0.8, 3.0, n_m)
    print(f"\n======== morse2d {case}  Kss={kss}  at_min_ref={rec['at_min']}  <st>_ref={rec['mean_stretch_ref']:.1f}")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    err = np.full((n_xh, km.size), np.inf)
    t0 = time.perf_counter()
    for i, kx in enumerate(kxh):
        for j, kM in enumerate(km):
            row = ownmin_stretch_of(sess, rec, kx, kss, f"kXH={kx:.3f} KM={kM:.2f}", after=after_morse(kM))
            err[i, j] = row["rmse_stretch"]
    print(f"  morse2d {time.perf_counter()-t0:.1f}s")
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    ij = np.unravel_index(int(np.argmin(err)), err.shape)
    np.savez(out / "morse2d.npz", k_XH=kxh, k_Morse=km, errors=err, Kss=kss, k_XX=k_xx,
             best_k_XH=float(kxh[ij[0]]), best_k_Morse=float(km[ij[1]]), best_rmse=float(err[ij]))
    heatmap(err, km, kxh, "K_Morse (Morse α)", f"{case} k_XH",
            f"{case} Morse×k_XH  Kss={kss}  best kXH={kxh[ij[0]]:.3f} α={km[ij[1]]:.2f}  RMSE={err[ij]:.1f}",
            (float(km[ij[1]]), float(kxh[ij[0]])), out / "morse2d.png")
    print(f"  best k_XH={kxh[ij[0]]:.4f}  K_Morse={km[ij[1]]:.3f}  RMSE={err[ij]:.2f}")


def reduce_morse2d(cases):
    root = OUT / "pyscf_pbe_ncs"
    rows, kxh, km = [], None, None
    for c in cases:
        z = np.load(root / c / "morse2d.npz")
        if kxh is None:
            kxh, km = np.asarray(z["k_XH"]), np.asarray(z["k_Morse"])
        rows.append((c, np.asarray(z["errors"])))
        print(f"  {c} minRMSE={float(np.min(z['errors'])):.1f} at kXH={float(z['best_k_XH']):.3f} α={float(z['best_k_Morse']):.2f}")
    R = np.stack([r[1] for r in rows], axis=0)
    mean, mx = R.mean(axis=0), R.max(axis=0)
    ij = np.unravel_index(int(np.argmin(mx)), mx.shape)
    heatmap(mx, km, kxh, "K_Morse (Morse α)", "k_XH",
            f"C Morse×k_XH minimax  kXH={kxh[ij[0]]:.3f} α={km[ij[1]]:.2f}  maxRMSE={mx[ij]:.1f}",
            (float(km[ij[1]]), float(kxh[ij[0]])), root / "morse2d_C.png")
    js = {"k_XH": kxh.tolist(), "k_Morse": km.tolist(), "minimax_k_XH": float(kxh[ij[0]]), "minimax_k_Morse": float(km[ij[1]]),
          "minimax_max_RMSE": float(mx[ij]), "minimax_mean_RMSE": float(mean[ij]),
          "per_case_min": {c: float(np.min(e)) for c, e in rows}}
    (root / "morse2d.json").write_text(json.dumps(js, indent=2))
    print(f"  minimax k_XH={kxh[ij[0]]:.4f}  K_Morse={km[ij[1]]:.3f}  maxRMSE={mx[ij]:.1f}  meanRMSE={mean[ij]:.1f}")


def opteval_case(sess_cls, case):
    rec = load_pyscf_nc(case)
    pack = json.loads((OUT / "pyscf_pbe_ncs" / "opt_trial.json").read_text())
    kss = float(pack.get("Kss", pack.get("k_XXX", 0.488)))
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    row = ownmin_stretch_of(sess, rec, pack["k_XH"], kss, case, after=after_split_morse(pack))
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    js = {"case": case, **pack, **row}
    (out / "opt_eval.json").write_text(json.dumps(js, indent=2))
    print(f"  saved {out / 'opt_eval.json'}")


def anneal_c(cases, nstep=12):
    """SA over (k_XH, k_Morse, k_XXX, k_XXH, k_HXH). Aggregate = mean own-min stretch RMSE. One crystal per subprocess."""
    root = OUT / "pyscf_pbe_ncs"
    k_xh0, kss0, _ = general_c_pack()
    p2 = root / "morse2d.json"
    if p2.is_file():
        j2 = json.loads(p2.read_text())
        k_xh0, kM0 = float(j2["minimax_k_XH"]), float(j2["minimax_k_Morse"])
    else:
        kM0 = 1.5
    names = ["k_XH", "k_Morse", "k_XXX", "k_XXH", "k_HXH"]
    bounds = np.array([[1.45, 2.15], [0.6, 3.5], [0.20, 1.10], [0.20, 1.10], [0.20, 1.10]])
    x = np.array([k_xh0, kM0, kss0, kss0, kss0], dtype=float)
    hist = []
    rng = np.random.default_rng(7)

    def pack_of(v):
        d = {n: float(a) for n, a in zip(names, v)}
        d.update({"k_XX": 1.0, "Kss": float(v[2]), "s_H_EvdW": 1.0, "nb": True, "morse": True, "split_angles": True})
        return d

    def loss(v):
        (root / "opt_trial.json").write_text(json.dumps(pack_of(v), indent=2))
        run_jobs([["--mol", "pyscf_nc", "--ref", f"opteval:{c}"] for c in cases], require_ok=True)
        rms = []
        for c in cases:
            rms.append(float(json.loads((root / c / "opt_eval.json").read_text())["rmse_stretch"]))
        mean, mx = float(np.mean(rms)), float(np.max(rms))
        rec = {**pack_of(v), "mean_RMSE": mean, "max_RMSE": mx, "per_case": dict(zip(cases, rms))}
        hist.append(rec)
        print(f"  LOSS mean={mean:.2f} max={mx:.2f}  { {n: round(float(a),3) for n,a in zip(names,v)} }", flush=True)
        return mean

    fx = loss(x)
    best, fbest = x.copy(), fx
    T0 = 25.0
    for i in range(nstep):
        T = T0 * (1.0 - i / nstep)
        y = x + rng.normal(0.0, 0.07, size=x.size) * (0.25 + T / T0)
        y = np.clip(y, bounds[:, 0], bounds[:, 1])
        fy = loss(y)
        if fy <= fx or rng.random() < np.exp(-(fy - fx) / max(T, 1e-6)):
            x, fx = y, fy
        if fy < fbest:
            best, fbest = y.copy(), fy
    js = {"protocol": "SA mean own-min stretch RMSE; surface Morse Exclusion2; split XXX/XXH/HXH from apars Kss",
          "cases": cases, "nstep": nstep, "best": pack_of(best), "best_mean_RMSE": fbest, "history": hist}
    (root / "anneal.json").write_text(json.dumps(js, indent=2))
    print(f"  SA best meanRMSE={fbest:.2f}  {pack_of(best)}")


ANNEAL_KNOBS = ["k_XH", "k_Morse", "k_XXX", "k_XXH", "k_HXH"]


def plot_anneal_diag():
    """RMSE vs trial and 5-knob trajectories from anneal.json. No MMFF."""
    root = OUT / "pyscf_pbe_ncs"
    js = json.loads((root / "anneal.json").read_text())
    hist = js["history"]
    if not hist:
        raise RuntimeError("plot_anneal_diag: empty history")
    n = len(hist)
    step = np.arange(n)
    mean = np.array([h["mean_RMSE"] for h in hist], float)
    mx = np.array([h["max_RMSE"] for h in hist], float)
    cases = js["cases"]
    per = {c: np.array([h["per_case"][c] for h in hist], float) for c in cases}
    ibest, bi, bm = int(np.argmin(mean)), 0, np.inf
    run_i = []
    for i, m in enumerate(mean):
        if m < bm:
            bm, bi = m, i
        run_i.append(bi)
    run_mean = np.array([mean[i] for i in run_i])
    run_max = np.array([mx[i] for i in run_i])
    print("\n======== SA history  (stretch RMSE cm-1)")
    print(f"  {'i':>3} {'mean':>8} {'max':>8} " + " ".join(f"{c:>16}" for c in cases) + "  " + " ".join(f"{k:>8}" for k in ANNEAL_KNOBS))
    for i, h in enumerate(hist):
        mark = " *" if i == ibest else ""
        pcs = " ".join(f"{h['per_case'][c]:16.2f}" for c in cases)
        ks = " ".join(f"{h[k]:8.3f}" for k in ANNEAL_KNOBS)
        print(f"  {i:3d} {h['mean_RMSE']:8.2f} {h['max_RMSE']:8.2f} {pcs}  {ks}{mark}")
    h0, hb = hist[0], hist[ibest]
    print(f"  start mean={h0['mean_RMSE']:.2f} max={h0['max_RMSE']:.2f}  →  best[i={ibest}] mean={hb['mean_RMSE']:.2f} max={hb['max_RMSE']:.2f}")
    print(f"  Δmean={hb['mean_RMSE']-h0['mean_RMSE']:+.2f}  Δmax={hb['max_RMSE']-h0['max_RMSE']:+.2f}  ntrial={n}  running-best last 4: {np.round(run_mean[-4:], 2)}")
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.plot(step, mean, "k.-", lw=1.4, ms=6, label="mean (loss)")
    ax.plot(step, mx, "k.--", lw=1.0, ms=4, label="max")
    for c, y in per.items():
        ax.plot(step, y, "o-", lw=1.0, ms=4, label=c)
    ax.plot(step, run_mean, color="#dc2626", lw=2.0, label="running-best mean")
    ax.plot(step, run_max, color="#dc2626", lw=1.2, ls=":", label="running-best max")
    ax.axvline(ibest, color="#dc2626", lw=1.0, alpha=0.7)
    ax.set_xlabel("SA trial")
    ax.set_ylabel("own-min stretch RMSE (cm⁻¹)")
    ax.set_title(f"anneal  n={n}  best mean {h0['mean_RMSE']:.1f}→{hb['mean_RMSE']:.1f}  max {h0['max_RMSE']:.1f}→{hb['max_RMSE']:.1f}")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.set_ylim(0, None)
    fig.tight_layout()
    p_rmse = root / "anneal_rmse.png"
    fig.savefig(p_rmse, dpi=160); plt.close(fig)
    print(f"  saved {p_rmse}")
    fig, axes = plt.subplots(len(ANNEAL_KNOBS), 1, figsize=(7.6, 1.35 * len(ANNEAL_KNOBS)), sharex=True)
    for ax, k in zip(axes, ANNEAL_KNOBS):
        y = np.array([h[k] for h in hist], float)
        yb = np.array([hist[i][k] for i in run_i], float)
        ax.plot(step, y, "o", ms=4, color="#64748b", label="trial")
        ax.plot(step, yb, "-", lw=1.8, color="#0f766e", label="running best")
        ax.axvline(ibest, color="#dc2626", lw=0.9, alpha=0.7)
        ax.set_ylabel(k, fontsize=9)
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    axes[0].set_title("SA stiffness knobs  (grey=evaluated, teal=incumbent)")
    axes[-1].set_xlabel("SA trial")
    fig.tight_layout()
    p_par = root / "anneal_params.png"
    fig.savefig(p_par, dpi=160); plt.close(fig)
    print(f"  saved {p_par}")
    return js


def anneal_spectra_case(sess_cls, case):
    """Own-min spectra at SA start pack and best pack. Same crystal, one process."""
    rec = load_pyscf_nc(case)
    js = json.loads((OUT / "pyscf_pbe_ncs" / "anneal.json").read_text())
    start, best = js["history"][0], js["best"]
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    kss_s, kss_b = float(start.get("Kss", start["k_XXX"])), float(best.get("Kss", best["k_XXX"]))
    print(f"\n======== anneal spectra {case}")
    p0 = fire_then_hessian(sess, rec["masses"], start["k_XH"], 1.0, kss_s, "SA start", rec["pos"], after_scales=after_split_morse(start))
    p1 = fire_then_hessian(sess, rec["masses"], best["k_XH"], 1.0, kss_b, "SA best", rec["pos"], after_scales=after_split_morse(best))
    out = OUT / "pyscf_pbe_ncs" / case / "anneal"
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "mmff_start_frequencies_cm1.npy", p0["freqs"])
    np.save(out / "mmff_best_frequencies_cm1.npy", p1["freqs"])
    e0 = stretch_rmse(np.sort(p0["freqs"][p0["freqs"] > 0]), rec["ref"], rec["nH"])
    e1 = stretch_rmse(np.sort(p1["freqs"][p1["freqs"] > 0]), rec["ref"], rec["nH"])
    from pyBall.FFfit_plots import plot_ownmin_method_dos
    plot_ownmin_method_dos(
        [{"label": "PySCF PBE", "omega_cm": rec["omega_all"], "color": "#111827", "require_minimum": False},
         {"label": f"SA start  RMSE={e0:.0f}", "omega_cm": p0["freqs"], "color": "#2563eb"},
         {"label": f"SA best  RMSE={e1:.0f}", "omega_cm": p1["freqs"], "color": "#d97706"}],
        out, title=f"{case}: SA start vs best  (MMFF own min, surface Morse)")
    print(f"  {case} RMSE start={e0:.2f} → best={e1:.2f}  saved {out / 'spectra_overlay.png'}")


def joint1d_case(sess_cls, case, n=17):
    """Own-min stretch RMSE vs k_XH at frozen Kss. One crystal; parent reduces across crystals."""
    rec = load_pyscf_nc(case)
    kss = 0.40 if rec["has_si"] else 0.488
    k_grid = np.linspace(0.90, 1.25, n) if rec["has_si"] else np.linspace(1.50, 2.05, n)
    print(f"\n======== joint1d {case}  n={n}  Kss={kss}  at_min_ref={rec['at_min']}")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    rmses = np.full(n, np.inf)
    t0 = time.perf_counter()
    for i, k in enumerate(k_grid):
        pack = fire_then_hessian(sess, rec["masses"], float(k), 1.0, kss, f"k={k:.3f}", rec["pos"])
        posf = np.sort(pack["freqs"][pack["freqs"] > 0.0])
        e = stretch_rmse(posf, rec["ref"], rec["nH"])
        e = np.inf if e is None else e
        rmses[i] = e
        print(f"  k_XH={k:.4f}  RMSE_st={e:.2f}  νmax={float(pack['freqs'].max()):.1f}")
    print(f"  joint1d {time.perf_counter()-t0:.1f}s")
    out = OUT / "pyscf_pbe_ncs" / case
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / "joint1d.npz", k_XH=k_grid, rmse_stretch=rmses, Kss=kss, k_XX=1.0,
             at_min=rec["at_min"], nH=rec["nH"], mean_stretch_ref=rec["mean_stretch_ref"])
    print(f"  saved {out / 'joint1d.npz'}")


def reduce_joint(cases):
    """One (k_XH, Kss) per element: minimax stretch RMSE across all crystals of that element."""
    root = OUT / "pyscf_pbe_ncs"
    general = {"protocol": "minimax own-min stretch RMSE; Kss frozen at XH4; k_XX=1; FIRE then Hessian",
               "reference": str(PYSCF_NC)}
    for sp, kss in (("C", 0.488), ("Si", 0.40)):
        sp_cases = [c for c in cases if c.endswith("_" + sp)]
        rows = []
        k_ref = None
        for c in sp_cases:
            p = root / c / "joint1d.npz"
            if not p.is_file():
                raise FileNotFoundError(p)
            z = np.load(p)
            k = np.asarray(z["k_XH"], dtype=float)
            if k_ref is None:
                k_ref = k
            elif k.shape != k_ref.shape or np.max(np.abs(k - k_ref)) > 1e-12:
                raise ValueError(f"{c}: k_XH grid mismatch")
            rows.append((c, np.asarray(z["rmse_stretch"], dtype=float), bool(z["at_min"])))
        R = np.vstack([r[1] for r in rows])
        mean, mx = R.mean(axis=0), R.max(axis=0)
        i_max = int(np.argmin(mx))
        i_mean = int(np.argmin(mean))
        rec = {"species": sp, "k_XH": float(k_ref[i_max]), "Kss": kss, "k_XX": 1.0,
               "criterion": "minimax_stretch_RMSE",
               "minimax_k": float(k_ref[i_max]), "minimax_max_RMSE": float(mx[i_max]), "minimax_mean_RMSE": float(mean[i_max]),
               "meanopt_k": float(k_ref[i_mean]), "meanopt_mean_RMSE": float(mean[i_mean]), "meanopt_max_RMSE": float(mx[i_mean]),
               "per_case_at_minimax": {c: float(r[i_max]) for c, r, _ in rows},
               "cases": sp_cases}
        general[sp] = rec
        fig, ax = plt.subplots(figsize=(7.4, 4.6))
        for c, r, at_min in rows:
            ax.plot(k_ref, r, lw=1.1, label=c + ("" if at_min else " *"))
        ax.plot(k_ref, mean, "k-", lw=2.0, label="mean")
        ax.plot(k_ref, mx, "k--", lw=1.6, label="max (minimax)")
        ax.axvline(k_ref[i_max], color="red", lw=1.2, label=f"minimax k={k_ref[i_max]:.3f}")
        ax.set_xlabel("k_XH scale")
        ax.set_ylabel("own-min stretch RMSE (cm⁻¹)")
        ax.set_title(f"{sp} general FF  Kss={kss}  k_XX=1  (* DFT not a min)")
        ax.legend(frameon=False, fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(root / f"joint1d_{sp}.png", dpi=160)
        plt.close(fig)
        print(f"  {sp} minimax k_XH={rec['k_XH']:.4f}  maxRMSE={rec['minimax_max_RMSE']:.1f}  meanRMSE={rec['minimax_mean_RMSE']:.1f}")
        print(f"     mean-opt k_XH={rec['meanopt_k']:.4f}  meanRMSE={rec['meanopt_mean_RMSE']:.1f}  maxRMSE={rec['meanopt_max_RMSE']:.1f}")
        for c, e in rec["per_case_at_minimax"].items():
            print(f"     {c:24s} RMSE_st={e:.1f}")
    (root / "general_ff.json").write_text(json.dumps(general, indent=2))
    print(f"  saved {root / 'general_ff.json'}")
    return general


def ownmin_case(sess_cls, case):
    rec = load_pyscf_nc(case)
    ff = json.loads((OUT / "pyscf_pbe_ncs" / "general_ff.json").read_text())
    sp = "Si" if rec["has_si"] else "C"
    kxh, kss = float(ff[sp]["k_XH"]), float(ff[sp]["Kss"])
    print(f"\n======== ownmin {case}  general {sp}  k_XH={kxh:.4f} Kss={kss:.4f}  DFT_min={rec['at_min']}")
    sess = sess_cls(rec["pos"], rec["sym"], firecore_path=str(FIRECORE), enable_angles=True)
    van = fire_then_hessian(sess, rec["masses"], 1.0, 1.0, 1.0, "vanilla", rec["pos"])
    fit = fire_then_hessian(sess, rec["masses"], kxh, 1.0, kss, "fitted", rec["pos"])
    out = OUT / "pyscf_pbe_ncs" / case / "ownmin"
    out.mkdir(parents=True, exist_ok=True)
    write_xyz(out / "mmff_vanilla_relaxed.xyz", rec["sym"], van["pos"], f"MMFF vanilla own min fmax={van['rel']['fmax']:.3e}")
    write_xyz(out / "mmff_fitted_relaxed.xyz", rec["sym"], fit["pos"], f"MMFF kXH={kxh:.4f} Kss={kss:.4f} own min fmax={fit['rel']['fmax']:.3e}")
    np.save(out / "mmff_vanilla_frequencies_cm1.npy", van["freqs"])
    np.save(out / "mmff_fitted_frequencies_cm1.npy", fit["freqs"])
    np.save(out / "mmff_vanilla_hessian.npy", van["H"])
    np.save(out / "mmff_fitted_hessian.npy", fit["H"])
    def _pack(p, tag):
        posf = np.sort(p["freqs"][p["freqs"] > 0.0])
        return {"tag": tag, "k_XH": p["k_XH"], "Kss": p["Kss"], "n_steps": p["rel"]["n_steps"], "fmax": p["rel"]["fmax"],
                "Rg": p["Rg"], "nu_max": float(p["freqs"].max()), "rmse_stretch": stretch_rmse(posf, rec["ref"], rec["nH"]),
                "mean_stretch": float(posf[-rec["nH"]:].mean()) if posf.size >= rec["nH"] else None,
                "n_imag": int(p["stats"]["n_imag"])}
    js = {"case": case, "protocol": "FIRE then Hessian at that k; start = PySCF xyz",
          "dft_at_min": rec["at_min"], "n_imag_ref": rec["n_imag_ref"], "nH": rec["nH"],
          "mean_stretch_ref": rec["mean_stretch_ref"], "general": {"species": sp, "k_XH": kxh, "Kss": kss, "k_XX": 1.0},
          "vanilla": _pack(van, "vanilla"), "fitted": _pack(fit, "fitted")}
    (out / "ownmin.json").write_text(json.dumps(js, indent=2))
    from pyBall.FFfit_plots import plot_ownmin_method_dos
    dft_lab = "PySCF PBE" if rec["at_min"] else "PySCF PBE (not a min)"
    plot_ownmin_method_dos(
        [{"label": dft_lab, "omega_cm": rec["omega_all"], "color": "#111827", "require_minimum": False},
         {"label": "MMFF vanilla", "omega_cm": van["freqs"], "color": "#2563eb"},
         {"label": f"MMFF kXH={kxh:.3f} K={kss:.3f}", "omega_cm": fit["freqs"], "color": "#d97706"}],
        out, title=f"{case}: vanilla / fitted MMFF / PySCF PBE  (MMFF at own min)")
    print(f"  vanilla RMSE_st={js['vanilla']['rmse_stretch']}  fitted RMSE_st={js['fitted']['rmse_stretch']}  saved {out}")


def replot_ownmin(case):
    rec = load_pyscf_nc(case)
    ff = json.loads((OUT / "pyscf_pbe_ncs" / "general_ff.json").read_text())
    sp = "Si" if rec["has_si"] else "C"
    kxh, kss = float(ff[sp]["k_XH"]), float(ff[sp]["Kss"])
    out = OUT / "pyscf_pbe_ncs" / case / "ownmin"
    van = np.load(out / "mmff_vanilla_frequencies_cm1.npy")
    fit = np.load(out / "mmff_fitted_frequencies_cm1.npy")
    from pyBall.FFfit_plots import plot_ownmin_method_dos
    dft_lab = "PySCF PBE" if rec["at_min"] else "PySCF PBE (not a min)"
    plot_ownmin_method_dos(
        [{"label": dft_lab, "omega_cm": rec["omega_all"], "color": "#111827", "require_minimum": False},
         {"label": "MMFF vanilla", "omega_cm": van, "color": "#2563eb"},
         {"label": f"MMFF kXH={kxh:.3f} K={kss:.3f}", "omega_cm": fit, "color": "#d97706"}],
        out, title=f"{case}: vanilla / fitted MMFF / PySCF PBE  (MMFF at own min)")
    print(f"  replotted {out / 'spectra_overlay.png'}")


def plot_ownmin_gallery(cases):
    from pyBall.FFfit_utils import gaussian_spectrum, as_signed_wavenumbers_cm1
    root = OUT / "pyscf_pbe_ncs"
    for sp, xmax in (("C", 3300.0), ("Si", 2400.0)):
        sp_cases = [c for c in cases if c.endswith("_" + sp)]
        fig, axes = plt.subplots(len(sp_cases), 1, figsize=(10.4, 1.55 * len(sp_cases)), sharex=True)
        if len(sp_cases) == 1:
            axes = [axes]
        grid = np.linspace(0.0, xmax, int(xmax) + 1)
        for ax, case in zip(axes, sp_cases):
            od = root / case / "ownmin"
            dft = as_signed_wavenumbers_cm1(np.load(resolve_pyscf_vib_case_dir(PYSCF_NC, case) / "frequencies_cm1.npy"))
            van = np.load(od / "mmff_vanilla_frequencies_cm1.npy")
            fit = np.load(od / "mmff_fitted_frequencies_cm1.npy")
            for lab, om, col, ls in (("PySCF PBE", dft, "#111827", "-"),
                                     ("MMFF vanilla", van, "#2563eb", "--"),
                                     ("MMFF fitted", fit, "#d97706", "-")):
                v = om[np.isfinite(om) & (om > 10.0)]
                y = gaussian_spectrum(v, grid, 8.0)
                ax.plot(grid, y / max(float(y.max()), 1e-18), color=col, ls=ls, lw=1.15, label=lab)
            ax.set_ylabel(case, fontsize=8)
            ax.set_ylim(0.0, 1.08)
            ax.set_xlim(0.0, xmax)
        axes[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
        axes[0].set_title(f"{sp} L1: whole spectrum  (own-min MMFF vs PySCF PBE)")
        axes[-1].set_xlabel("ω (cm⁻¹)")
        fig.tight_layout()
        fig.savefig(root / f"gallery_ownmin_{sp}.png", dpi=160)
        plt.close(fig)
        print(f"  saved {root / f'gallery_ownmin_{sp}.png'}")


def list_pyscf_nc_cases():
    if not PYSCF_NC.is_dir():
        raise FileNotFoundError(PYSCF_NC)
    return list_pyscf_vib_cases(PYSCF_NC)


def nc_scale_sets(has_si):
    """Transfer packs + any NC-grid best.json already on disk."""
    if has_si:
        packs = [("default", 1.0, 1.0, 1.0),
                 ("SiH4_B3LYP", 1.138, 1.0, 0.40),
                 ("SiH4_pbc", 1.025, 1.0, 0.30)]
        for case in ("octahedron_Si", "cube_5ring_Si"):
            p = OUT / "pyscf_pbe_ncs" / case / "best.json"
            if p.is_file():
                j = json.loads(p.read_text())
                packs.append((f"fit_{case}", float(j["scale_XH"]), 1.0, float(j["scale_angle_K"])))
        return packs
    packs = [("default", 1.0, 1.0, 1.0),
             ("CH4_B3LYP", 2.038, 1.0, 0.488),
             ("CH4_3ob", 1.812, 1.0, 0.488),
             ("adamantane_3ob", 1.76, 1.0, 0.68)]
    for case in ("octahedron_C", "cube_5ring_C"):
        p = OUT / "pyscf_pbe_ncs" / case / "best.json"
        if p.is_file():
            j = json.loads(p.read_text())
            packs.append((f"fit_{case}", float(j["scale_XH"]), 1.0, float(j["scale_angle_K"])))
    return packs


def write_nc_summary():
    root = OUT / "pyscf_pbe_ncs"
    rows = []
    for p in sorted(root.glob("*/eval.json")) + sorted(root.glob("*/best.json")):
        j = json.loads(p.read_text())
        kind = p.name
        if kind == "best.json":
            rows.append({"case": j["case"], "kind": "grid", "at_min": j.get("at_min"), "n_imag_ref": j.get("n_imag_ref"),
                         "tag": "grid_best", "k_XH": j.get("scale_XH"), "Kss": j.get("scale_angle_K"),
                         "rmse_stretch": j.get("rmse_stretch"), "mean_stretch_ref": j.get("mean_stretch_ref")})
            continue
        for e in j.get("evals", []):
            rows.append({"case": j["case"], "kind": "eval", "at_min": j.get("at_min"), "n_imag_ref": j.get("n_imag_ref"),
                         "tag": e["tag"], "k_XH": e["k_XH"], "Kss": e["Kss"],
                         "rmse_stretch": e["rmse_stretch"], "mean_stretch": e.get("mean_stretch"),
                         "mean_stretch_ref": j.get("mean_stretch_ref"), "n_imag_mmff": e.get("n_imag_mmff")})
    (root / "summary.json").write_text(json.dumps({"reference": str(PYSCF_NC), "protocol": "FFfit frozen PySCF xyz; loss=nH highest freqs; k_XX=1", "rows": rows}, indent=2))
    print("\n======== PySCF-NC FFfit (stretch RMSE vs PBE/cc-pVDZ)")
    print(f"  {'case':<24} {'tag':<22} {'k_XH':>6} {'Kss':>6} {'RMSE_st':>8} {'<st>':>8} {'ref<st>':>8} {'min?':>5}")
    for r in rows:
        rmse = r.get("rmse_stretch")
        st = r.get("mean_stretch")
        print(f"  {r['case']:<24} {r.get('tag',''):<22} {r.get('k_XH') or 0:6.3f} {r.get('Kss') or 0:6.3f} "
              f"{rmse if rmse is not None else float('nan'):8.1f} {st if st is not None else float('nan'):8.1f} "
              f"{r.get('mean_stretch_ref') or float('nan'):8.1f} {str(bool(r.get('at_min'))):>5}")
    evals = [r for r in rows if r["kind"] == "eval"]
    if evals:
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), sharey=True)
        for ax, suff, title in ((axes[0], "_C", "C NCs"), (axes[1], "_Si", "Si NCs")):
            cases = sorted({r["case"] for r in evals if r["case"].endswith(suff)})
            tags = sorted({r["tag"] for r in evals if r["case"].endswith(suff)})
            x = np.arange(len(cases))
            w = 0.8 / max(len(tags), 1)
            for i, tag in enumerate(tags):
                ys = []
                for c in cases:
                    hit = [r["rmse_stretch"] for r in evals if r["case"] == c and r["tag"] == tag]
                    ys.append(hit[0] if hit and hit[0] is not None else np.nan)
                ax.bar(x + (i - 0.5 * (len(tags) - 1)) * w, ys, width=w, label=tag)
            ax.set_xticks(x)
            ax.set_xticklabels(cases, rotation=30, ha="right", fontsize=8)
            ax.set_title(title)
            ax.legend(frameon=False, fontsize=7)
            ax.set_ylabel("stretch RMSE (cm⁻¹)" if ax is axes[0] else "")
        fig.suptitle("MMFF k vs PySCF PBE/cc-pVDZ  (FFfit, frozen xyz, nH stretches)")
        fig.tight_layout()
        fig.savefig(root / "eval_stretch_rmse.png", dpi=160)
        plt.close(fig)
        print(f"  saved {root / 'eval_stretch_rmse.png'}")
    print(f"  saved {root / 'summary.json'}")


def write_summary():
    summary = {"angle_buffer": "apars[:,2]=Kss; cols 0,1 = (cos θ0/2, sin θ0/2) untouched",
               "do_not_use": "June 2026 angle×0.30 (that was sin(θ0/2), not Kss)"}
    for p in sorted(OUT.glob("*/best.json")) + sorted(OUT.glob("pyscf_pbe_ncs/*/best.json")):
        summary[str(p.parent.relative_to(OUT))] = json.loads(p.read_text())
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    keys = [k for k in summary if k not in ("angle_buffer", "do_not_use")]
    print("\n======== SUMMARY keys", keys)
    for k in keys:
        rec = summary[k]
        bits = [k]
        for fld in ("scale_XH", "scale_CH", "scale_k", "scale_XX", "scale_CC", "scale_angle_K", "rmse_all", "rmse_all9", "rmse_gt_fmin", "rmse_gt500"):
            if fld in rec:
                bits.append(f"{fld}={rec[fld]}")
        print(" ", "  ".join(str(b) for b in bits))


def run_jobs(jobs, require_ok=True):
    py = sys.executable
    script = str(Path(__file__).resolve())
    failed = []
    for extra in jobs:
        cmd = [py, script, *extra]
        print("\n>>>>", " ".join(cmd), flush=True)
        r = subprocess.run(cmd, cwd=os.getcwd())
        if r.returncode != 0:
            failed.append(extra)
            if require_ok:
                raise SystemExit(r.returncode)
            print(f"  FAILED rc={r.returncode}  {extra}", flush=True)
    return failed


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mol", default="all_si",
                   choices=["CH4", "C2H6", "SiH4", "Si2H6", "adamantane", "Si10H16", "pyscf_nc", "all", "all_si", "all_cages"])
    p.add_argument("--n", type=int, default=0)
    p.add_argument("--ref", default="", help="XH4/X2H6 tag, cage job, PySCF-NC case, or eval[/case]")
    args = p.parse_args()
    os.environ.setdefault("FIRECORE_PATH", str(FIRECORE))
    OUT.mkdir(parents=True, exist_ok=True)
    n_xh4 = args.n or 17
    n_x2 = args.n or 15
    n_cage = args.n or 11
    n_nc = args.n or 9
    if args.mol == "pyscf_nc":
        cases = list_pyscf_nc_cases()
        ref = args.ref or "all"
        if ref == "all":
            fit_cases = [c for c in ("octahedron_C", "cube_5ring_C", "octahedron_Si") if c in cases]
            run_jobs([["--mol", "pyscf_nc", "--ref", c, "--n", str(n_nc)] for c in fit_cases], require_ok=True)
            failed = run_jobs([["--mol", "pyscf_nc", "--ref", f"eval:{c}"] for c in cases], require_ok=False)
            write_nc_summary()
            write_summary()
            if failed:
                raise SystemExit(f"pyscf_nc eval failed for {failed}")
            print("OUT", OUT / "pyscf_pbe_ncs"); return
        if ref == "joint" or ref.startswith("joint1d:") or ref.startswith("ownmin:"):
            n_joint = args.n or 17
            if ref == "joint":
                run_jobs([["--mol", "pyscf_nc", "--ref", f"joint1d:{c}", "--n", str(n_joint)] for c in cases], require_ok=True)
                reduce_joint(cases)
                failed = run_jobs([["--mol", "pyscf_nc", "--ref", f"ownmin:{c}"] for c in cases], require_ok=False)
                plot_ownmin_gallery(cases)
                if failed:
                    raise SystemExit(f"pyscf_nc ownmin failed for {failed}")
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            if ref.startswith("joint1d:"):
                joint1d_case(MMFFMolecularSession, ref.split(":", 1)[1], n=n_joint)
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            ownmin_case(MMFFMolecularSession, ref.split(":", 1)[1])
            print("OUT", OUT / "pyscf_pbe_ncs" / ref.split(":", 1)[1] / "ownmin"); return
        if ref == "nbscan" or ref.startswith("nbscan:"):
            nb_cases = ["octahedron_C", "cube_5ring_C"]
            if ref == "nbscan":
                run_jobs([["--mol", "pyscf_nc", "--ref", f"nbscan:{c}"] for c in nb_cases if c in cases], require_ok=True)
                reduce_nbscan([c for c in nb_cases if c in cases])
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            nbscan_case(MMFFMolecularSession, ref.split(":", 1)[1])
            print("OUT", OUT / "pyscf_pbe_ncs" / ref.split(":", 1)[1]); return
        if ref == "morse2d" or ref.startswith("morse2d:"):
            m_cases = ["octahedron_C", "cube_5ring_C"]
            n_m2 = args.n or 4
            if ref == "morse2d":
                run_jobs([["--mol", "pyscf_nc", "--ref", f"morse2d:{c}", "--n", str(n_m2)] for c in m_cases if c in cases], require_ok=True)
                reduce_morse2d([c for c in m_cases if c in cases])
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            morse2d_case(MMFFMolecularSession, ref.split(":", 1)[1], n_xh=n_m2, n_m=n_m2)
            print("OUT", OUT / "pyscf_pbe_ncs" / ref.split(":", 1)[1]); return
        if ref == "anneal" or ref.startswith("opteval:"):
            a_cases = ["octahedron_C", "cube_5ring_C"]
            if ref == "anneal":
                anneal_c([c for c in a_cases if c in cases], nstep=args.n or 12)
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            opteval_case(MMFFMolecularSession, ref.split(":", 1)[1])
            print("OUT", OUT / "pyscf_pbe_ncs" / ref.split(":", 1)[1]); return
        if ref == "anneal_plot" or ref.startswith("annealspec:"):
            a_cases = ["octahedron_C", "cube_5ring_C"]
            if ref == "anneal_plot":
                plot_anneal_diag()
                run_jobs([["--mol", "pyscf_nc", "--ref", f"annealspec:{c}"] for c in a_cases if c in cases], require_ok=True)
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            anneal_spectra_case(MMFFMolecularSession, ref.split(":", 1)[1])
            print("OUT", OUT / "pyscf_pbe_ncs" / ref.split(":", 1)[1] / "anneal"); return
        if ref == "replot" or ref.startswith("replot:"):
            if ref == "replot":
                for c in cases:
                    replot_ownmin(c)
                plot_ownmin_gallery(cases)
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            replot_ownmin(ref.split(":", 1)[1])
            print("OUT", OUT / "pyscf_pbe_ncs"); return
        if ref == "eval" or ref.startswith("eval:"):
            if ref == "eval":
                failed = run_jobs([["--mol", "pyscf_nc", "--ref", f"eval:{c}"] for c in cases], require_ok=False)
                write_nc_summary()
                if failed:
                    raise SystemExit(f"pyscf_nc eval failed for {failed}")
                print("OUT", OUT / "pyscf_pbe_ncs"); return
            case = ref.split(":", 1)[1]
            from pyBall import MMFF as _MMFF
            _MMFF.setVerbosity(0, 0)
            from mmff_molecular_session import MMFFMolecularSession
            rec = load_pyscf_nc(case)
            eval_pyscf_nc_scales(MMFFMolecularSession, case, nc_scale_sets(rec["has_si"]))
            print("OUT", OUT / "pyscf_pbe_ncs" / case); return
        from pyBall import MMFF as _MMFF
        _MMFF.setVerbosity(0, 0)
        from mmff_molecular_session import MMFFMolecularSession
        fit_pyscf_nc(MMFFMolecularSession, ref, n=n_nc)
        write_summary()
        print("OUT", OUT / "pyscf_pbe_ncs" / ref); return
    if args.mol == "all":
        run_jobs([
            ["--mol", "CH4", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_xh4)],
            ["--mol", "CH4", "--ref", "dftb_3ob-3-1", "--n", str(n_xh4)],
            ["--mol", "C2H6", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_x2)],
            ["--mol", "adamantane", "--n", str(n_cage)],
            ["--mol", "SiH4", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_xh4)],
            ["--mol", "SiH4", "--ref", "dftb_pbc-0-3", "--n", str(n_xh4)],
            ["--mol", "Si2H6", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_x2)],
            ["--mol", "Si2H6", "--ref", "dftb_pbc-0-3", "--n", str(n_x2)],
            ["--mol", "Si10H16", "--n", str(n_cage)],
        ])
        write_summary(); print("OUT", OUT); return
    if args.mol == "all_si":
        run_jobs([
            ["--mol", "SiH4", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_xh4)],
            ["--mol", "SiH4", "--ref", "dftb_pbc-0-3", "--n", str(n_xh4)],
            ["--mol", "Si2H6", "--ref", "pyscf_b3lyp_cc-pVDZ", "--n", str(n_x2)],
            ["--mol", "Si2H6", "--ref", "dftb_pbc-0-3", "--n", str(n_x2)],
            ["--mol", "Si10H16", "--n", str(n_cage)],
        ])
        write_summary(); print("OUT", OUT); return
    if args.mol == "all_cages":
        run_jobs([["--mol", "adamantane", "--n", str(n_cage)], ["--mol", "Si10H16", "--n", str(n_cage)]])
        write_summary(); print("OUT", OUT); return
    from pyBall import MMFF as _MMFF
    _MMFF.setVerbosity(0, 0)
    from mmff_molecular_session import MMFFMolecularSession
    if args.mol == "CH4":
        fit_xh4(MMFFMolecularSession, "CH4", args.ref or "pyscf_b3lyp_cc-pVDZ", n=n_xh4)
    elif args.mol == "SiH4":
        fit_xh4(MMFFMolecularSession, "SiH4", args.ref or "pyscf_b3lyp_cc-pVDZ", n=n_xh4)
    elif args.mol == "C2H6":
        ref = args.ref or "pyscf_b3lyp_cc-pVDZ"
        scale_xh = float(xh4_best_path("CH4", ref if ref != "pyscf_b3lyp_cc-pVDZ" else "pyscf_b3lyp_cc-pVDZ")["scale_XH"])
        if ref == "dftb_3ob-3-1":
            scale_xh = float(xh4_best_path("CH4", "dftb_3ob-3-1")["scale_XH"])
        fit_x2h6(MMFFMolecularSession, "C2H6", ref, scale_xh, n=n_x2, fmin=500.0)
    elif args.mol == "Si2H6":
        ref = args.ref or "pyscf_b3lyp_cc-pVDZ"
        parent = "dftb_pbc-0-3" if ref == "dftb_pbc-0-3" else "pyscf_b3lyp_cc-pVDZ"
        scale_xh = float(xh4_best_path("SiH4", parent)["scale_XH"])
        fit_x2h6(MMFFMolecularSession, "Si2H6", ref, scale_xh, n=n_x2, fmin=200.0)
    elif args.mol == "adamantane":
        job = args.ref or "adamantane_3ob_3_1"
        _, sess = fit_cage(MMFFMolecularSession, "adamantane", job, n=n_cage, mode="xh", k_xx=1.0)
        k_xh_fit = float(json.loads((OUT / "adamantane_xh" / "best.json").read_text())["scale_k"])
        fit_cage(MMFFMolecularSession, "adamantane", job, n=n_cage, mode="xx", k_xh=k_xh_fit, sess=sess)
    elif args.mol == "Si10H16":
        job = args.ref or "Si10H_pbc"
        _, sess = fit_cage(MMFFMolecularSession, "Si10H16", job, n=n_cage, mode="xh", k_xx=1.0)
        k_xh_fit = float(json.loads((OUT / "Si10H16_xh" / "best.json").read_text())["scale_k"])
        fit_cage(MMFFMolecularSession, "Si10H16", job, n=n_cage, mode="xx", k_xh=k_xh_fit, sess=sess)
    write_summary()
    print("OUT", OUT)


if __name__ == "__main__":
    main()
