#!/usr/bin/env python3
"""Comprehensive optimizer parameter scan across force fields (MMFF / UFF / eFF).

For each (force field, molecule):
  1. static grid:  FIRE dt x falpha (eFF: dt x damping),  Anderson eta x m,  BB1
  2. adaptive:     Anderson eta=1e-3 + adapt bits {1,5,7} at m in {5,8}
  3. hybrid:       FIRE_k -> Anderson m in {1,3,5} at best-found eta (quick pass)

Metric: force evaluations to |F|<1e-3 and |F|<1e-6 (computed from outF trace),
plus final energy. Everything runs from the same starting geometry.

Usage (each invocation = separate process, avoids cross-FF state issues):
  python3 scan_all_ff.py --ff mmff --mol propandiol
  python3 scan_all_ff.py --ff uff  --mol xylitol
  python3 scan_all_ff.py --ff eff                      # H2O free nuclei (forked workers)
  python3 scan_all_ff.py --table                       # render combined table from .npz
"""
import os, sys, argparse
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, '../../..'))
sys.path.insert(0, _REPO)                                  # repo root for pyBall
sys.path.insert(0, os.path.join(_REPO, 'tests/tEFF'))      # fire_study_lib

NMAX   = 10000
FCONV  = 1e-6
DTS    = [0.02, 0.05, 0.1, 0.2, 0.5]
FALPHAS= [0.8, 0.95, 0.99]
ETAS   = [1e-3, 3e-3, 1e-2, 3e-2, 0.1]
MS     = [1, 3, 5, 8]
ONLY_ADAPT = os.environ.get('ONLY_ADAPT', '0') == '1'   # re-run only adaptive rows, merge into existing npz
FIRE_MECH  = os.environ.get('FIRE_MECH',  '0') == '1'   # re-run only FIRE-mechanism rows (autoMass/massMetric/modeSupp), merge

# per-molecule tuned FIRE (dt,falpha) from the static grid -> mechanism variants run at the fair baseline
FIRE_BEST = {
    ('mmff','propandiol'): (0.2, 0.8), ('mmff','xylitol'): (0.2, 0.95), ('mmff','nHexadecan_dicarboxylic'): (0.2, 0.8),
    ('uff','propandiol'): (0.1, 0.95), ('uff','xylitol'): (0.1, 0.99), ('uff','guanine'): (0.05, 0.99), ('uff','benzene'): (0.5, 0.95),
}
# name -> (autoMass_mode, massMetric, modeSupp, msMu, msEma); all run via mmff.run_FIRE (DynamicOpt FIRE, same impl for MMFF+UFF)
FIRE_MECHS = [
    ('FIRE* base',       (0, 0, 0, 0.0, 0.7)),
    ('FIRE* aM1',        (1, 0, 0, 0.0, 0.7)),
    ('FIRE* aM2',        (2, 0, 0, 0.0, 0.7)),
    ('FIRE* mMet',       (0, 1, 0, 0.0, 0.7)),
    ('FIRE* aM2+mMet',   (2, 1, 0, 0.0, 0.7)),
    ('FIRE* mSup.7',     (0, 0, 1, 0.7, 0.9)),
    ('FIRE* mSup.9',     (0, 0, 1, 0.9, 0.5)),
    ('FIRE* aM2+mSup.7', (2, 0, 1, 0.7, 0.9)),
    ('FIRE* mMet+mSup.7',(0, 1, 1, 0.7, 0.9)),
    ('FIRE* mSup.5',     (0, 0, 1, 0.5, 0.7)),   # negative-evidence row: mu<=0.5 stalls near basin floor
]

RES = os.path.join(_HERE, 'scan_results')
os.makedirs(RES, exist_ok=True)


def milestones(outF, n, thrs=(1e-3, 1e-6)):
    out = {}
    for t in thrs:
        hit = np.nonzero(outF[:n] < t)[0]
        out[t] = int(hit[0] + 1) if len(hit) else None
    return out


def run_row(rows, name, fn_run, p0, DOFs):
    """reset DOFs, run fn -> (n, outE, outF), record milestones.
    All drivers store |F| in outF. Milestones over finite entries —
    run()/run_FIRE return itr (off-by-one: entries written = itr+1 on break)."""
    DOFs[:] = p0
    outE = np.full(NMAX, np.nan); outF = np.full(NMAX, np.nan)
    try:
        n = fn_run(outE, outF)
    except SystemExit:
        rows.append((name, None, None, np.nan, np.nan, 'CRASH')); print(f'  {name:34s} CRASH', flush=True); return
    n = int(np.isfinite(outE).sum())   # entries actually written (driver-independent)
    E = outE[n-1] if n > 0 else np.nan; F = outF[n-1] if n > 0 else np.nan
    ms = milestones(outF, n)
    conv = F < FCONV
    rows.append((name, ms[1e-3], ms[1e-6], E, F, 'ok' if conv else f'budget'))
    print(f'  {name:34s} n={n:5d} |F|1e-3={str(ms[1e-3]):>5s} |F|1e-6={str(ms[1e-6]):>5s} E={E:10.5f} {"" if conv else "(budget)"}', flush=True)


# ============================= MMFF / UFF =============================

def init_mmff(mol, bUFF):
    from pyBall import MMFF as mmff
    mmff.setVerbosity(verbosity=0, idebug=0)
    kw = dict(xyz_name=f'cpp/common_resources/xyz/{mol}.xyz' if bUFF else f'tests/tMMFF/data/xyz/{mol}.xyz',
              sElementTypes='cpp/common_resources/ElementTypes.dat', sAtomTypes='cpp/common_resources/AtomTypes.dat',
              sBondTypes='cpp/common_resources/BondTypes.dat', sAngleTypes='cpp/common_resources/AngleTypes.dat',
              sDihedralTypes='cpp/common_resources/DihedralTypes.dat',
              bMMFF=not bUFF, bUFF=bUFF, b141=True, bSimple=True, bConj=True, bCumulene=True, nPBC=(0, 0, 0))
    mmff.init(**kw)
    if bUFF:
        mmff.getBuffs_UFF()
        mmff.setSwitches(PBC=-1, NonBonded=1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)
        mmff.setSwitchesUFF_NB(NonBonded=1, NonBondNeighs=-1)
    else:
        mmff.getBuffs()
        mmff.setSwitches(PBC=-1, NonBonded=1, MMFF=1, SurfAtoms=-1, GridFF=-1)
    E0 = mmff.eval()
    return mmff, mmff.DOFs.copy(), E0


def merge_rows(npz_path, new_rows, E0):
    """ONLY_ADAPT/FIRE_MECH mode: update new rows in an existing npz, keep statics.
    FIRE_MECH also purges stale 'FIRE*' rows (mechanism set changed between runs)."""
    old = list(np.load(npz_path, allow_pickle=True)['rows']) if os.path.exists(npz_path) else []
    names = {r[0] for r in new_rows}
    return [r for r in old if r[0] not in names and not (FIRE_MECH and r[0].startswith('FIRE*'))] + new_rows


def scan_mmff_uff(mol, bUFF):
    mmff, p0, E0 = init_mmff(mol, bUFF)
    DOFs = mmff.DOFs
    rows = []
    print(f'=== {"UFF" if bUFF else "MMFF"} {mol}  E0={E0:.5f} nDOFs={len(p0)} ===', flush=True)

    if FIRE_MECH:
        dt, fa = FIRE_BEST[('uff' if bUFF else 'mmff', mol)]
        print(f'-- FIRE mechanisms at tuned dt={dt} fa={fa} (run_FIRE = DynamicOpt FIRE) --', flush=True)
        for name, (am, mm, ms, mu, ema) in FIRE_MECHS:
            def f(oE, oF, am=am, mm=mm, ms=ms, mu=mu, ema=ema, dt=dt, fa=fa):
                mmff.set_opt(dt_max=dt, falpha=fa)
                mmff.setAutoMass(am, 100.0)                          # resets am_init (fresh stiffness state per run)
                mmff.setOptParam('massMetric', mm)
                mmff.setOptParam('modeSupp', ms); mmff.setOptParam('msMu', mu); mmff.setOptParam('msEma', ema); mmff.setOptParam('msInit', 0)
                return mmff.run_FIRE(nstepMax=NMAX, dt=dt, Fconv=FCONV, outE=oE, outF=oF)
            run_row(rows, name, f, p0, DOFs)
    elif not ONLY_ADAPT:
        print('-- FIRE grid (dt x falpha) --', flush=True)
        for dt in DTS:
            for fa in FALPHAS:
                def f(oE, oF, dt=dt, fa=fa):
                    mmff.set_opt(dt_max=dt, falpha=fa)
                    if bUFF:
                        return mmff.run_FIRE(nstepMax=NMAX, dt=dt, Fconv=FCONV, outE=oE, outF=oF)
                    return mmff.run(nstepMax=NMAX, dt=dt, Fconv=FCONV, ialg=2, outE=oE, outF=oF)
                run_row(rows, f'FIRE dt={dt:g} fa={fa:g}', f, p0, DOFs)  # run_no_omp/run_FIRE/run_GH all store |F| (the F^2 convention is only in the unused virtual run())

        print('-- Anderson grid (eta x m) --', flush=True)
        for eta in ETAS:
            for m in MS:
                def f(oE, oF, eta=eta, m=m):
                    return mmff.run_GH(nstepMax=NMAX, Fconv=FCONV, iMethod=2, mHist=m, eta=eta, a_init=-1, outE=oE, outF=oF)
                run_row(rows, f'And m={m} eta={eta:g}', f, p0, DOFs)

        print('-- BB1 --', flush=True)
        run_row(rows, 'BB1 a_init=1e-3', lambda oE, oF: mmff.run_GH(nstepMax=NMAX, Fconv=FCONV, iMethod=1, mHist=1, eta=1e-3, a_init=1e-3, outE=oE, outF=oF), p0, DOFs)

    if not FIRE_MECH:
        print('-- adaptive Anderson (eta=1e-3 + adapt bits) --', flush=True)
        for adapt in (1, 5, 7):
            for m in (5, 8):
                def f(oE, oF, adapt=adapt, m=m):
                    return mmff.run_GH(nstepMax=NMAX, Fconv=FCONV, iMethod=2, mHist=m, eta=1e-3, a_init=-1, adapt=adapt, amax=3.0, outE=oE, outF=oF)
                run_row(rows, f'And m={m} eta=1e-3 adapt={adapt}', f, p0, DOFs)

    tag = 'uff' if bUFF else 'mmff'
    path = os.path.join(RES, f'{tag}_{mol}.npz')
    if ONLY_ADAPT or FIRE_MECH: rows = merge_rows(path, rows, E0)
    np.savez(path, rows=np.array(rows, dtype=object), E0=E0)
    print(f'saved {path}', flush=True)
    return rows


# ============================= eFF =============================

def scan_eff():
    """H2O frozen-core, free nuclei. Forked workers isolate exit(1) crashes."""
    import multiprocessing as mp
    import fire_study_lib as fs
    from pyBall import eFF as eff

    fs.set_atom_params()
    na, ne, _ = fs.load_e('H2O_pairs_fc.xyz', verbosity=0)
    p0 = fs.snapshot()
    print(f'=== eFF H2O_pairs_fc free nuclei na={na} ne={ne} nDOFs={eff.nDOFs} ===', flush=True)

    EFF_DTS   = [0.01, 0.03, 0.05, 0.07, 0.1]
    EFF_DAMPS = [0.5, 0.7, 0.9]

    jobs = []
    if FIRE_MECH:
        for name, (am, mm, ms, mu, ema) in FIRE_MECHS:
            jobs.append((name, ('firex', 0.03, 0.5, am, mm, ms, mu, ema)))   # tuned dt=0.03 damp=0.5
    elif not ONLY_ADAPT:
        for dt in EFF_DTS:
            for d in EFF_DAMPS:
                jobs.append((f'FIRE dt={dt:g} damp={d:g}', ('fire', dt, d)))
        for eta in ETAS:
            for m in MS:
                jobs.append((f'And m={m} eta={eta:g}', ('gh', 2, m, eta, 0, 0.0)))
        jobs.append(('BB1 a_init=1e-3', ('gh', 1, 1, 1e-3, 0, 0.0)))
    if not FIRE_MECH:
        for adapt in (1, 5, 7):
            for m in (5, 8):
                jobs.append((f'And m={m} eta=1e-3 adapt={adapt}', ('gh', 2, m, 1e-3, adapt, 3.0)))

    rows = []
    ctx = mp.get_context('fork')
    for name, spec in jobs:
        q = ctx.Queue()
        p = ctx.Process(target=_eff_child, args=(p0, spec, q))
        p.start(); p.join(300)
        if p.is_alive():
            p.terminate(); res = ('timeout',)
        elif q.empty():
            res = ('crash', p.exitcode)
        else:
            res = q.get()
        if res[0] == 'ok':
            _, n, E, F, m1, m2 = res
            conv = F < FCONV
            rows.append((name, m1, m2, E, F, 'ok' if conv else 'budget'))
            print(f'  {name:34s} n={n:5d} |F|1e-3={str(m1):>5s} |F|1e-6={str(m2):>5s} E={E:10.5f} {"" if conv else "(budget)"}', flush=True)
        else:
            rows.append((name, None, None, np.nan, np.nan, res[0]))
            print(f'  {name:34s} {res[0]}', flush=True)
    path = os.path.join(RES, 'eff_H2O_free.npz')
    if ONLY_ADAPT or FIRE_MECH: rows = merge_rows(path, rows, np.nan)
    np.savez(path, rows=np.array(rows, dtype=object), E0=np.nan)
    print(f'saved {RES}/eff_H2O_free.npz', flush=True)


def _eff_child(p0, spec, q):
    import fire_study_lib as fs
    from pyBall import eFF as eff
    import numpy as np
    try:
        fs.restore(p0)
        eff.setAutoMass(0, 100.0); eff.setMassMetric(0)
        if spec[0] == 'firex':
            _, dt, damp, am, mm, ms, mu, ema = spec
            eff.initOpt(dt=dt, damping=damp, f_limit=100.0, bMass=False)
            eff.setAutoMass(am, 100.0); eff.setMassMetric(mm)
            eff.setOptParam('modeSupp', ms); eff.setOptParam('msMu', mu); eff.setOptParam('msEma', ema); eff.setOptParam('msInit', 0)
            outE, outF = eff.run(NMAX, dt, FCONV, 2)
            n = int(np.isfinite(outE).sum())
        elif spec[0] == 'fire':
            _, dt, damp = spec
            eff.initOpt(dt=dt, damping=damp, f_limit=100.0, bMass=False)
            outE, outF = eff.run(NMAX, dt, FCONV, 2)
            n = int(np.isfinite(outE).sum())
        else:
            _, im, m, eta, adapt, amax = spec
            n, outE, outF = eff.run_GH(NMAX, FCONV, im, m, eta, -1.0, adapt, amax)
        k = min(n, len(outF)); outF = outF[:k]
        m1 = next((i + 1 for i, f in enumerate(outF) if f < 1e-3), None)
        m2 = next((i + 1 for i, f in enumerate(outF) if f < 1e-6), None)
        q.put(('ok', n, outE[k - 1], outF[k - 1], m1, m2))
    except SystemExit as e:
        q.put(('exit', e.code))
    except Exception as e:
        q.put(('err', repr(e)))
    q.close(); q.join_thread()
    os._exit(0)


# ============================= table =============================

def table():
    files = sorted(f for f in os.listdir(RES) if f.endswith('.npz'))
    for fn in files:
        d = np.load(os.path.join(RES, fn), allow_pickle=True)
        rows = d['rows']
        print(f'\n### {fn[:-4]}')
        print(f'{"method":36s} {"|F|<1e-3":>9s} {"|F|<1e-6":>9s} {"E_final":>11s}')
        for r in rows:
            print(f'{r[0]:36s} {str(r[1]):>9s} {str(r[2]):>9s} {r[3]:11.5f}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--ff', choices=['mmff', 'uff', 'eff'])
    ap.add_argument('--mol', default=None)
    ap.add_argument('--table', action='store_true')
    a = ap.parse_args()
    os.chdir(_REPO)   # repo root for relative paths
    if a.table:
        table()
    elif a.ff in ('mmff', 'uff'):
        bUFF = a.ff == 'uff'
        mols = a.mol.split(',') if a.mol else (['propandiol', 'xylitol', 'guanine', 'benzene'] if bUFF else ['propandiol', 'xylitol', 'nHexadecan_dicarboxylic'])
        if len(mols) > 1:   # repeated mmff.init() in one process leaves DOFs aliased to freed buffers (UAF) -> one subprocess per molecule
            import subprocess
            env = dict(os.environ); env['PYTHONUNBUFFERED'] = '1'
            for mol in mols:
                r = subprocess.run([sys.executable, __file__, '--ff', a.ff, '--mol', mol], env=env)
                if r.returncode != 0:
                    print(f'!! {mol} exited {r.returncode}', flush=True)
        else:
            scan_mmff_uff(mols[0], bUFF)
    elif a.ff == 'eff':
        os.chdir(os.path.join(_REPO, 'tests/tEFF'))  # fire_study_lib loads H2O_pairs_fc.xyz relatively
        scan_eff()
