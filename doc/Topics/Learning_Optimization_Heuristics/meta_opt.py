#!/usr/bin/env python3
"""Meta-optimization of DynamicOpt knobs via scipy differential_evolution (§9.1).

Finds the per-system (or global) optimum of the FULL knob vector — FIRE or
gradient-history (GH) family — measured in force evaluations to |F|<1e-6.
This is the control experiment every learned/adaptive method must beat.

Knob vectors (bounded, decoded in _decode):
  fire: log10 dt, falpha, fdec, finc, damp_max, autoMass{0,1,2}, massMetric{0,1},
        modeSupp{0,1}, msMu, msEma
  gh:   log10 eta, mHist{1..16}, adapt bits{0..31}, amax, lsEmaCoef, lsMomBeta,
        lsMomGain, lsTurnCos, iMethod{2,3}, autoMass{0,1,2}

Cost: evals to |F|<1e-6; non-converged -> 2*NMAX + 200*log10(|F|end).
Global mode (--systems with >1 entry): cost = sum over systems of evals/BEST[sys]
(normalized so no system dominates).

Each candidate is evaluated inside a persistent worker process per system
(init once at spawn; a LINE-SEARCH-FAIL exit(1) kills the worker -> parent
detects, penalizes, respawns). Same pattern as scan_all_ff.py.

Usage:
  python3 meta_opt.py --family fire --systems mmff:propandiol
  python3 meta_opt.py --family gh   --systems mmff:propandiol,uff:propandiol,eff:H2O --name global3
  python3 meta_opt.py --family fire --systems eff:H2O
"""
import os, sys, json, time, argparse
import numpy as np
import multiprocessing as mp

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, '../../..'))
sys.path.insert(0, _REPO)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_REPO, 'tests/tEFF'))
os.chdir(_REPO)   # relative paths for params/xyz

NMAX, FCONV = 10000, 1e-6
META = os.path.join(_HERE, 'meta_results'); os.makedirs(META, exist_ok=True)

# best-known evals to 1e-6 across all tested methods (for normalized global cost)
BEST = {('mmff','propandiol'):345, ('mmff','xylitol'):473, ('mmff','nHexadecan_dicarboxylic'):1720,
        ('uff','propandiol'):544, ('uff','xylitol'):382, ('uff','guanine'):188, ('uff','benzene'):29,
        ('eff','H2O'):268}

FIRE_BOUNDS = [(-2.0,-0.2),(0.70,0.999),(0.30,0.90),(1.01,1.60),(0.02,0.50),
               (0,2.999),(0,1.999),(0,1.999),(0.50,0.99),(0.30,0.97)]
GH_BOUNDS   = [(-4.0,-0.3),(1,8.999),(0,31.999),(1.0,6.0),(0.02,0.6),   # m<=8: DynamicOpt hard-caps mMax=8 (m>8 -> worker exit)
               (0.0,0.95),(0.0,4.0),(0.0,0.95),(2,3.999),(0,2.999)]

def decode(family, x):
    if family == 'fire':
        return dict(dt=10**x[0], falpha=x[1], fdec=x[2], finc=x[3], damp=x[4],
                    am=int(x[5]), mm=int(x[6]), ms=int(x[7]), mu=x[8], ema=x[9])
    return dict(eta=10**x[0], m=int(x[1]), adapt=int(x[2]), amax=x[3], emaC=x[4],
                beta=x[5], gain=x[6], tc=x[7], im=int(x[8]), am=int(x[9]))

# ---------------- workers (persistent, one per system) ----------------

def _mmff_worker(mol, bUFF, qin, qout):
    import scan_all_ff as S
    mmff, p0, E0 = S.init_mmff(mol, bUFF)
    DOFs = mmff.DOFs
    qout.put(('ready', E0))
    while True:
        msg = qin.get()
        if msg[0] == 'stop':
            qout.put(('bye',)); return
        _, tag, family, th = msg
        try:
            qout.put((tag, _eval_mmff(mmff, DOFs, p0, family, th)))
        except Exception as e:
            qout.put((tag, ('err', repr(e))))

def _eval_mmff(mmff, DOFs, p0, family, th):
    DOFs[:] = p0
    # reset all opt-in state every evaluation (state must not leak between candidates)
    mmff.setAutoMass(0, 100.0)
    for k, v in [('massMetric',0),('modeSupp',0),('msInit',0),('lsMomBeta',0.5),
                 ('lsMomGain',1.0),('lsTurnCos',0.5),('lsEmaCoef',0.2),('lsEtaLo',1e-4),('lsEtaHi',0.2)]:
        mmff.setOptParam(k, v)
    mmff.set_opt(dt_max=0.1, falpha=0.8, fdec=0.5, finc=1.1, damp_max=0.2)   # defaults first
    oE = np.full(NMAX, np.nan); oF = np.full(NMAX, np.nan)
    if family == 'fire':
        mmff.set_opt(dt_max=th['dt'], falpha=th['falpha'], fdec=th['fdec'], finc=th['finc'], damp_max=th['damp'])
        if th['am']: mmff.setAutoMass(th['am'], 100.0)
        mmff.setOptParam('massMetric', th['mm'])
        if th['ms']:
            mmff.setOptParam('modeSupp', 1); mmff.setOptParam('msMu', th['mu'])
            mmff.setOptParam('msEma', th['ema']); mmff.setOptParam('msInit', 0)
        mmff.run_FIRE(nstepMax=NMAX, dt=th['dt'], Fconv=FCONV, outE=oE, outF=oF)
    else:
        mmff.setOptParam('lsEmaCoef', th['emaC']); mmff.setOptParam('lsMomBeta', th['beta'])
        mmff.setOptParam('lsMomGain', th['gain']); mmff.setOptParam('lsTurnCos', th['tc'])
        if th['am']: mmff.setAutoMass(th['am'], 100.0)
        mmff.run_GH(nstepMax=NMAX, Fconv=FCONV, iMethod=th['im'], mHist=th['m'], eta=th['eta'],
                    a_init=-1.0, adapt=th['adapt'], amax=th['amax'], outE=oE, outF=oF)
    return _cost(oE, oF)

def _eff_worker(qin, qout):
    os.chdir(os.path.join(_REPO, 'tests/tEFF'))   # fire_study_lib loads H2O_pairs_fc.xyz relatively
    import fire_study_lib as fs
    from pyBall import eFF as eff
    fs.set_atom_params()
    na, ne, _ = fs.load_e('H2O_pairs_fc.xyz', verbosity=0)
    p0 = fs.snapshot()
    qout.put(('ready', np.nan))
    while True:
        msg = qin.get()
        if msg[0] == 'stop':
            qout.put(('bye',)); return
        _, tag, family, th = msg
        try:
            qout.put((tag, _eval_eff(fs, eff, p0, family, th)))
        except Exception as e:
            qout.put((tag, ('err', repr(e))))

def _eval_eff(fs, eff, p0, family, th):
    fs.restore(p0)
    eff.setAutoMass(0, 100.0); eff.setMassMetric(0)
    for k, v in [('modeSupp',0),('msInit',0),('lsMomBeta',0.5),('lsMomGain',1.0),
                 ('lsTurnCos',0.5),('lsEmaCoef',0.2)]:
        eff.setOptParam(k, v)
    if family == 'fire':
        eff.initOpt(dt=th['dt'], damping=th['damp'], f_limit=100.0, bMass=False)
        for k, kk in [('fdec','fdec'),('finc','finc'),('falpha','falpha')]:
            eff.setOptParam(kk, th[k])
        if th['am']: eff.setAutoMass(th['am'], 100.0)
        eff.setMassMetric(th['mm'])
        if th['ms']:
            eff.setOptParam('modeSupp', 1); eff.setOptParam('msMu', th['mu'])
            eff.setOptParam('msEma', th['ema']); eff.setOptParam('msInit', 0)
        oE, oF = eff.run(NMAX, th['dt'], FCONV, 2)
    else:
        eff.initOpt(dt=0.03, damping=0.5, f_limit=100.0, bMass=False)
        eff.setOptParam('lsEmaCoef', th['emaC']); eff.setOptParam('lsMomBeta', th['beta'])
        eff.setOptParam('lsMomGain', th['gain']); eff.setOptParam('lsTurnCos', th['tc'])
        if th['am']: eff.setAutoMass(th['am'], 100.0)
        n, oE, oF = eff.run_GH(NMAX, FCONV, th['im'], th['m'], th['eta'], -1.0, th['adapt'], th['amax'])
    return _cost(oE, oF)

def _cost(oE, oF):
    n = int(np.isfinite(oE).sum())
    if n == 0: return ('crash',)
    F = oF[:n]
    hit = np.nonzero(F < 1e-6)[0]
    evals = int(hit[0] + 1) if len(hit) else None
    if evals is None:
        c = 2*NMAX + 200*max(0.0, np.log10(max(F[-1],1e-30)/1e-6))
    else:
        c = float(evals)
    h3 = np.nonzero(F < 1e-3)[0]
    return ('ok', evals, int(h3[0]+1) if len(h3) else None, float(oE[n-1]), float(F[-1]), c)

# ---------------- driver ----------------

class System:
    """Persistent eval worker for one (ff,mol) system."""
    def __init__(self, ff, mol):
        self.key = (ff, mol); self.proc = None; self.qin = None; self.qout = None
    def start(self, ctx):
        self.qin, self.qout = ctx.Queue(), ctx.Queue()
        target, args = ((_eff_worker, (self.qin, self.qout)) if self.key[0] == 'eff'
                        else (_mmff_worker, (self.key[1], self.key[0] == 'uff', self.qin, self.qout)))
        self.proc = ctx.Process(target=target, args=args, daemon=True)
        self.proc.start()
        msg = self.qout.get(timeout=120)
        assert msg[0] == 'ready', msg
        print(f'[worker {self.key}] ready', flush=True)
    def eval(self, tag, family, th, ctx):
        """Send eval; on worker death return crash + mark for respawn."""
        import queue as _q
        self.qin.put(('eval', tag, family, th))
        t0 = time.time()
        while time.time() - t0 < 180:
            try:
                _, res = self.qout.get(timeout=0.05)
                return res
            except _q.Empty:
                if not self.proc.is_alive():
                    return ('dead',)
        return ('timeout',)
    def stop(self):
        if self.proc and self.proc.is_alive():
            self.qin.put(('stop',)); self.proc.join(10)
            if self.proc.is_alive(): self.proc.terminate()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--family', choices=['fire', 'gh'], required=True)
    ap.add_argument('--systems', required=True, help='comma list like mmff:propandiol,uff:xylitol,eff:H2O')
    ap.add_argument('--maxiter', type=int, default=40)
    ap.add_argument('--popsize', type=int, default=8)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--name', default=None)
    a = ap.parse_args()
    from scipy.optimize import differential_evolution

    systems = []
    for s in a.systems.split(','):
        ff, mol = s.split(':'); systems.append(System(ff, mol))
    ctx = mp.get_context('spawn')   # spawn, not fork: forked children deadlock on mp.Queue feeder-thread locks inherited mid-flight
    for s in systems: s.start(ctx)
    multi = len(systems) > 1
    bounds = FIRE_BOUNDS if a.family == 'fire' else GH_BOUNDS
    logf = os.path.join(META, f"{a.name or a.family+'_'+a.systems.replace(':','_').replace(',','+')}.log")
    evals_done = [0]

    def objective(x):
        th = decode(a.family, x)
        tot, det = 0.0, {}
        for s in systems:
            res = s.eval(evals_done[0], a.family, th, ctx)
            if res[0] == 'dead':
                print(f'  !! worker {s.key} died -> penalty + respawn', flush=True)
                s.start(ctx); res = ('dead',); c = 3*NMAX
            elif res[0] in ('err', 'crash', 'timeout'):
                print(f'  !! {s.key} {res[0]} {res[1:] if len(res)>1 else ""} -> penalty', flush=True); c = 3*NMAX
            else:
                _, ev, e3, E, F, c = res
                det[str(s.key)] = ev
            tot += c / (BEST.get(s.key, 1) if multi else 1)
        evals_done[0] += 1
        line = f'eval#{evals_done[0]:4d} cost={tot:9.2f} th={json.dumps(th,default=float)} sys={det}'
        print(line, flush=True)
        with open(logf, 'a') as f: f.write(line + '\n')
        return tot

    print(f'=== meta_opt family={a.family} systems={a.systems} popsize={a.popsize} maxiter={a.maxiter} ===', flush=True)
    t0 = time.time()
    r = differential_evolution(objective, bounds, popsize=a.popsize, maxiter=a.maxiter,
                               seed=a.seed, polish=False, tol=1e-4, updating='immediate', workers=1)
    dt_ = time.time() - t0
    best = decode(a.family, r.x)
    print(f'\n=== DONE {dt_:.0f}s cost={r.fun:.2f} ===\nbest={json.dumps(best, indent=1)}', flush=True)
    out = dict(family=a.family, systems=a.systems, cost=float(r.fun), theta=best, x=list(map(float, r.x)), secs=dt_)
    jf = logf.replace('.log', '.json')
    json.dump(out, open(jf, 'w'), indent=1)
    print(f'saved {jf}', flush=True)
    for s in systems: s.stop()

if __name__ == '__main__':
    main()
