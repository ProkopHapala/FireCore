#!/usr/bin/env python3
"""L2 C MMFF harmonic spectra: relax with THAT FF, then Hessian. One crystal per process.

Never evaluate MMFF frequencies at a DFTB geometry. See
doc/Topics/FTIR_Nanocrystals/Hessian_at_own_minimum.md

Usage (from tests/tMMFF so MMFF data/ exists):
  python3 ../tSiNCs/run_mmff_ownmin_l2.py rhombic_dodecahedron_C
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path
import numpy as np

FIRECORE = Path("/home/prokop/git/FireCore")
CCU = Path("/home/prokop/git/CompChemUtils")
sys.path.insert(0, str(FIRECORE))
sys.path.insert(0, str(CCU / "examples/tSiNCs"))
sys.path.insert(0, str(CCU))

ATLAS = FIRECORE / "tests/tSiNCs/OUT_chem_atlas/atlas/L2_dftb"
DFTB = Path("/home/prokop/SIMULATIONS/SiNCs/DFTB")
OUT = FIRECORE / "tests/tSiNCs/OUT_dftb_vs_mmff"
FCONV = 1e-4
NSTEP = 8000
SCALE_CH, SCALE_CC, SCALE_ANG = 1.8125, 1.0, 0.4875  # CH4 vs DFTB 3ob-3-1; Kss=apars[:,2]. See MMFF_stiffness_scaling.md


def hess_input_order(sess, hess):
    inv = sess._perm_in2int
    n = sess.n_atoms
    dof = np.empty(n * 3, dtype=np.int32)
    for ia in range(n):
        ii = int(inv[ia])
        dof[3 * ia + 0] = 3 * ii + 0
        dof[3 * ia + 1] = 3 * ii + 1
        dof[3 * ia + 2] = 3 * ii + 2
    return hess[dof][:, dof]


def write_xyz(path, symbols, pos, comment):
    pos = np.asarray(pos, dtype=float)
    with open(path, "w") as f:
        f.write(f"{len(symbols)}\n{comment}\n")
        for s, p in zip(symbols, pos):
            f.write(f"{s}  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}\n")


def freq_report(freqs, tag, nH):
    from pyBall.FFfit_utils import assert_harmonic_spectrum_at_minimum
    f = np.asarray(freqs, dtype=float)
    st = assert_harmonic_spectrum_at_minimum(f, ctx=f"{tag}: ")
    n_gap = int(np.sum((f > 10.0) & (f < 150.0)))
    vib = f[f > 10.0]
    first = float(vib.min()) if vib.size else float("nan")
    n_stretch = int(np.sum((f > 2700.0) & (f < 3300.0)))
    print(f"  {tag}: n_imag(|ν|>10)={st['n_imag']}  n_rigid={st['n_rigid']}  n_in_10_150={n_gap}  first_vib={first:.2f}  νmax={st['nu_max']:.1f}  n_stretch={n_stretch} (nH={nH})")
    return {"n_imag": st["n_imag"], "n_rigid": st["n_rigid"], "n_in_10_150": n_gap, "first_vib_cm1": first, "nu_max_cm1": st["nu_max"], "n_stretch_2700_3300": n_stretch}


def run_kind(sess, masses, nH, kind, outdir, symbols):
    from pyBall.FFfit_utils import signed_frequencies_cm1
    t0 = time.perf_counter()
    if kind == "default":
        sess.set_scales_per_bond_type(scale_ch=1.0, scale_cc=1.0, scale_angle=1.0)
    elif kind == "scaled":
        sess.set_scales_per_bond_type(scale_ch=SCALE_CH, scale_cc=SCALE_CC, scale_angle=SCALE_ANG)
        print(f"  scaled k: CH×{SCALE_CH}  CC×{SCALE_CC}  Kss×{SCALE_ANG}  (CH4 vs 3ob; not June angle×0.30)")
    else:
        raise ValueError(kind)
    rel = sess.relax(nstep_max=NSTEP, fconv=FCONV)
    pos = sess.positions_input_order()
    rg = float(np.sqrt(np.mean(np.sum((pos - pos.mean(0))**2, axis=1))))
    rg0 = float(np.sqrt(np.mean(np.sum((sess.positions - sess.positions.mean(0))**2, axis=1))))
    print(f"  {kind} relax: n_steps={rel['n_steps']}  fmax={rel['fmax']:.3e}  Rg={rg:.3f} Å (init {rg0:.3f})  ({time.perf_counter()-t0:.1f}s)")
    if rg < 0.5 * rg0:
        raise RuntimeError(f"{kind} FIRE collapsed the particle (Rg {rg0:.3f}→{rg:.3f} Å). k-scale must not move q0.")
    t1 = time.perf_counter()
    H = hess_input_order(sess, sess.get_hessian(dx=1e-4))
    print(f"  {kind} Hessian FD: {time.perf_counter()-t1:.1f}s  shape={H.shape}")
    if np.isnan(H).any() or np.isinf(H).any():
        raise ValueError(f"{kind}: NaN/Inf Hessian")
    freqs = signed_frequencies_cm1(H, masses)
    stats = freq_report(freqs, kind, nH)
    pos = sess.positions_input_order()
    write_xyz(outdir / f"mmff_{kind}_relaxed.xyz", symbols, pos, f"MMFF {kind} own minimum fmax={rel['fmax']:.3e}")
    np.save(outdir / f"mmff_{kind}_hessian.npy", H)
    np.save(outdir / f"mmff_{kind}_frequencies_cm1.npy", freqs)
    np.save(outdir / f"mmff_{kind}_pos.npy", pos)
    stats.update(rel)
    stats["kind"] = kind
    stats["protocol"] = "relax_then_hessian_own_minimum"
    stats["fconv"] = FCONV
    stats["Rg_A"] = rg
    stats["Rg_init_A"] = rg0
    return stats, pos


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_mmff_ownmin_l2.py <crystal_id>")
    case = sys.argv[1]
    mol2 = ATLAS / f"{case}.mol2"
    dftb = DFTB / case
    outdir = OUT / case
    outdir.mkdir(parents=True, exist_ok=True)
    if not mol2.is_file():
        raise FileNotFoundError(mol2)
    from pyBall.atomicUtils import loadMol2
    from pyBall import MMFF as _MMFF
    from mmff_molecular_session import MMFFMolecularSession
    _MMFF.setVerbosity(0, 0)
    apos, atypes, enames, qs, bonds, *rest = loadMol2(str(mol2))
    pos0 = np.asarray(apos, dtype=float)
    symbols = [str(s) for s in np.asarray(enames).reshape(-1)]
    masses = np.asarray(np.load(dftb / "masses.npy"), dtype=float)
    if masses.shape[0] != len(symbols):
        raise ValueError(f"masses {masses.shape} vs N={len(symbols)}")
    nH = int(sum(s == "H" for s in symbols))
    def _xyz(path):
        lines = Path(path).read_text().splitlines()
        n = int(lines[0]); rows = []
        for ln in lines[2:2+n]:
            t = ln.split(); rows.append([float(t[1]), float(t[2]), float(t[3])])
        return np.asarray(rows)
    pos_dftb = _xyz(dftb / "relaxed.xyz")
    if pos_dftb.shape != pos0.shape:
        raise ValueError(f"DFTB relaxed {pos_dftb.shape} vs mol2 {pos0.shape}")
    print(f"======== {case}  N={len(symbols)} nH={nH}  mol2={mol2.name}  Fconv={FCONV}")
    print(f"  |F| will be checked after EACH MMFF relax before Hessian")
    sess = MMFFMolecularSession(pos0, symbols, firecore_path=str(FIRECORE), enable_angles=True, mol2_path=str(mol2))
    st_d, pos_d = run_kind(sess, masses, nH, "default", outdir, symbols)
    st_s, pos_s = run_kind(sess, masses, nH, "scaled", outdir, symbols)
    def rmsd(a, b):
        return float(np.sqrt(np.mean((a - b) ** 2)))
    rec = {
        "crystal": case, "natoms": len(symbols), "nH": nH, "fconv": FCONV, "nstep_max": NSTEP,
        "mol2": str(mol2), "protocol": "each method relaxed with its own MMFF then Hessian at that minimum",
        "rmsd_mmff_default_vs_dftb_A": rmsd(pos_d, pos_dftb),
        "rmsd_mmff_scaled_vs_dftb_A": rmsd(pos_s, pos_dftb),
        "rmsd_mmff_scaled_vs_default_A": rmsd(pos_s, pos_d),
        "scaled_k": {"scale_CH": SCALE_CH, "scale_CC": SCALE_CC, "scale_angle_Kss": SCALE_ANG,
                     "source": "CH4 grid vs DFTB 3ob-3-1, apars[:,2]; tests/tSiNCs/OUT_mmff_kss_fit/"},
        "default": st_d, "scaled": st_s,
    }
    (outdir / "mmff_ownmin_protocol.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps({k: rec[k] for k in rec if k not in ("default", "scaled")}, indent=2))
    print(f"  default {st_d}")
    print(f"  scaled  {st_s}")
    print("DONE", case)


if __name__ == "__main__":
    main()
