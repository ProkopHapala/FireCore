#!/usr/bin/env python3
"""Branch lab (§13.4): measure learnable headroom with an offline teacher.

From checkpoints along a FIRE trajectory, apply structured memory/velocity/scale perturbations
("actions"), continue each branch for H force evaluations, and record energy outcomes.
Answers two questions before any learned optimizer is trained:
  1) how much better than the baseline continuation can the best branch be? (teacher headroom)
  2) can the permitted dynamical features predict which action wins?      (predictability)

All actions mutate ONLY optimizer state (velocity, masses, FIRE scalars) — positions are never
teleported; a branch is a counterfactual *decision*, not a restart.

Usage: python3 branch_lab.py [--system mmff_propandiol|uff_propandiol|eff_H2O] [--parents K] [--horizon H]
"""
import os, sys, argparse, json
import numpy as np
import ctypes as ct

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
sys.path.insert(0, os.path.join(REPO, 'pyBall'))

c_dbl_p = ct.POINTER(ct.c_double)

def bind_step_api(lib):
    lib.optBind.argtypes=[];            lib.optBind.restype=None
    lib.optInitOpt.argtypes=[ct.c_double,ct.c_double]; lib.optInitOpt.restype=None
    lib.optEval.argtypes=[];            lib.optEval.restype=ct.c_double
    lib.optMoveFIRE.argtypes=[];        lib.optMoveFIRE.restype=ct.c_double
    lib.optNDOFs.argtypes=[];           lib.optNDOFs.restype=ct.c_int
    lib.optBuffer.argtypes=[ct.c_int];  lib.optBuffer.restype=c_dbl_p
    lib.optStateSize.argtypes=[];       lib.optStateSize.restype=ct.c_int
    lib.optSaveState.argtypes=[c_dbl_p];lib.optSaveState.restype=ct.c_int
    lib.optRestoreState.argtypes=[c_dbl_p]; lib.optRestoreState.restype=None
    lib.setOptParam.argtypes=[ct.c_char_p, ct.c_double]; lib.setOptParam.restype=ct.c_double
    lib.setAutoMass.argtypes=[ct.c_int, ct.c_double];    lib.setAutoMass.restype=None


class Driver:
    """Stepwise FIRE driver on top of the opt.* step API. One system, one optimizer."""
    def __init__(self, lib, dt, damping, fire_params):
        self.lib=lib; lib.optBind(); lib.optInitOpt(dt,damping)
        for k,v in fire_params.items(): lib.setOptParam(k.encode(), v)
        self.n = lib.optNDOFs()
        self.pos   = np.ctypeslib.as_array(lib.optBuffer(0),(self.n,))
        self.vel   = np.ctypeslib.as_array(lib.optBuffer(1),(self.n,))
        self.force = np.ctypeslib.as_array(lib.optBuffer(2),(self.n,))
        self.invm  = np.ctypeslib.as_array(lib.optBuffer(3),(self.n,))
        self.bsz   = lib.optStateSize()
        self.dt0, self.damp0 = dt, damping
    def eval(self):  return self.lib.optEval()
    def move(self):  return self.lib.optMoveFIRE()
    def f2(self):    return float(self.force@self.force)
    def setp(self,k,v): self.lib.setOptParam(k.encode(),float(v))
    def save(self):
        sz=self.lib.optStateSize(); b=np.empty(sz); self.lib.optSaveState(b.ctypes.data_as(c_dbl_p)); return b  # size can grow when lazy buffers (autoMass,vSlow) allocate
    def restore(self,b):
        self.lib.optRestoreState(b.ctypes.data_as(c_dbl_p))

# ---------------- branch actions: memory/velocity/scale perturbations only ----------------
def make_actions():
    A=[('base', lambda d,c: None)]
    for a in (0.0,0.5,1.5,2.0):
        A.append((f'v{a}', lambda d,c,a=a: d.vel.__imul__(a)))
    for b in (0.25,0.5,1.0):
        def steer(d,c,b=b):
            f=d.force; nf=np.linalg.norm(f)
            if nf<1e-30: return
            w=d.vel + b*f/nf*np.linalg.norm(d.vel); nw=np.linalg.norm(w)
            if nw>1e-30: d.vel[:] = w*(np.linalg.norm(d.vel)/nw)
        A.append((f'steer{b}', steer))
    for g in (0.5,1.0):
        def killmode(d,c,g=g):
            u=c.get('u')
            if u is None: return
            d.vel -= g*(d.vel@u)*u
        A.append((f'killmode{g}', killmode))
    for s in (0.5,2.0):
        A.append((f'dt{s}', lambda d,c,s=s: dt_scale(d,s)))
    A.append(('dampreset', lambda d,c: (d.setp('damping',d.damp0), d.setp('lastNeg',0))))
    A.append(('fresh',     lambda d,c: (d.vel.__imul__(0.0), d.setp('damping',d.damp0), d.setp('lastNeg',0))))
    for s in (0.5,2.0):
        A.append((f'im{s}', lambda d,c,s=s: d.invm.__imul__(s)))
    def kick(d,c):
        s=c.get('s')
        if s is None: return
        nv=np.linalg.norm(d.vel); d.vel += 0.5*nv*s
    A.append(('kick', kick))
    return A

def dt_scale(d, s):   # read-modify-write through setOptParam's previous-value return
    old = d.lib.setOptParam(b'dt', 0.0); d.lib.setOptParam(b'dt', old*s)

def get_scalar(d, name):
    old=d.lib.setOptParam(name.encode(),0.0); d.lib.setOptParam(name.encode(),old); return old

# ---------------- param-switch actions: mechanism toggles on the STRONG (champion) baseline ----------------
# Each action flips opt params; restore() resets them (all serialized in the state blob).
def make_param_actions():
    A=[('base', lambda d,c: None)]
    for cg in (0.5,0.7,0.9):
        for g in (0.25,0.5):
            A.append((f'g2ts_c{cg}g{g}', lambda d,c,cg=cg,g=g: (d.setp('gate_cos',cg), d.setp('gate_g',g))))
    def _mson(d,c,mu):   # enable modeSupp once; msInit=0 only on first call so the EMA mode accumulates under --persist
        if c.get('ms_on'): return
        d.setp('modeSupp',1); d.setp('msMu',mu); d.setp('msInit',0); c['ms_on']=1
    A.append(('modeSupp',  lambda d,c: _mson(d,c,0.5)))
    A.append(('msKill',    lambda d,c: _mson(d,c,0.0)))
    A.append(('massMet1',  lambda d,c: d.setp('massMetric',1)))
    A.append(('massMet0',  lambda d,c: d.setp('massMetric',0)))
    for s in (0.67,1.5):
        A.append((f'dt{s}', lambda d,c,s=s: dt_scale(d,s)))
    A.append(('fa.99', lambda d,c: d.setp('falpha',0.99)))
    A.append(('fd.3',  lambda d,c: d.setp('fdec',0.3)))
    A.append(('v1.5',  lambda d,c: d.vel.__imul__(1.5)))
    A.append(('v0.5',  lambda d,c: d.vel.__imul__(0.5)))
    A.append(('amOff', lambda d,c: d.lib.setAutoMass(0,100.0)))
    return A

# verified champion configs (meta_results/*.json, clean-replay validated in §12/§15.9)
CHAMP = {
 'mmff_propandiol': dict(dt=0.1929774550017525, falpha=0.9586366974522995, fdec=0.8055255163316105, finc=1.2591175080572874, damp=0.22743413548266378, am=2, mm=0),
 'uff_propandiol':  dict(dt=0.14206022184072467, falpha=0.964287081970814,  fdec=0.8365692820117381, finc=1.0266027038409915, damp=0.26215187288341313, am=0, mm=1),
 'eff_H2O':         dict(dt=0.11106794702650019, falpha=0.8549622176007076, fdec=0.8691954800280984, finc=1.1975211964102057, damp=0.36628736177905297, am=2, mm=1),
}

# ---------------- Python FIRE engine: exact replica of DynamicOpt::move_FIRE ----------------
# Needed for register/memory policies the C++ step can't express (§14.6): policies hook between the
# FIRE velocity update and the leapfrog displacement. Verified trajectory-identical via --verify.
class PyFire:
    def __init__(self, drv, fire_params, f_limit):
        self.d=drv; self.p=dict(fire_params); self.p['f_limit']=f_limit
        self.p.setdefault('cv_kill',0.0); self.p.setdefault('minLastNeg',3)
        self.p.setdefault('fTminmax',0.1); self.p.setdefault('dr_limit',0.2); self.p.setdefault('v_limit',10.0)
        self.st={}
    def sync(self):   # pull mutable scalars from opt after a blob restore
        for k in ('dt','damping','lastNeg','dt_max','dt_min','damp_max'):
            self.st[k]=get_scalar(self.d,k)
    def step(self, hook=None, ctx=None):
        d=self.d; f=d.force; v=d.vel; fi=f*d.invm
        ff=float(fi@fi); vv=float(v@v); vf=float(v@fi)
        fl=np.sqrt(ff); vl=np.sqrt(vv)
        cos_vf=vf/(fl*vl+1e-32); ren=vl/(fl+1e-32)
        cv = self.p['cv_kill'] if cos_vf<0 else 1.0-self.st['damping']
        cf = (0.0 if cos_vf<0 else self.st['damping'])*ren
        v*=cv; v+=fi*cf
        if cos_vf<0.0:
            self.st['dt']=max(self.st['dt']*self.p['fdec'], self.st['dt_min']); self.st['lastNeg']=0; self.st['damping']=self.st['damp_max']
        else:
            if self.st['lastNeg']>self.p['minLastNeg']:
                self.st['dt']=min(self.st['dt']*self.p['finc'], self.st['dt_max']); self.st['damping']*=self.p['falpha']
            self.st['lastNeg']+=1
        dt_=self.st['dt']
        if ff>self.p['f_limit']**2: dt_*=np.sqrt(self.p['f_limit']/np.sqrt(ff))
        dr=(v+fi*dt_)*dt_; drmax=float(np.abs(dr).max())
        if drmax>self.p['dr_limit']:
            s=self.p['dr_limit']/drmax; v[:]=s*(v+fi*dt_)-fi*dt_
        if hook: hook(d, self.st, ctx)
        v+=fi*dt_; d.pos+=v*dt_

# ---------------- register/memory policies (hook = mutate v before displacement) ----------------
def make_policies():
    """Policy = callable(d, st, ctx) with private memory in ctx. Runs BEFORE leapfrog each step."""
    P=[('base', None)]
    for g in (0.5,1.0):   # two-timescale momentum: slow EMA register boosts velocity direction
        def slow(d,st,ctx,g=g):
            vs=ctx.setdefault('vs',np.zeros_like(d.vel))
            vs*=0.9; vs+=0.1*d.vel; d.vel+=g*vs
        P.append((f'2tslow{g}',slow))
    for g in (0.5,1.0):   # online transverse-mode suppression with learned u (EMA, sign-aligned)
        def mkill(d,st,ctx,g=g):
            fo=ctx.get('fold'); xo=ctx.get('xold')
            if fo is not None:
                s=d.pos-xo; ss=s@s
                if ss>1e-24:
                    df=d.force-fo; q=df-(df@s/ss)*s; q2=q@q
                    if q2>1e-12*(df@df):
                        u=ctx.setdefault('u',np.zeros_like(q))
                        sgn=1.0 if u@q>=0 else -1.0
                        u[:]=0.7*u+0.3*sgn*q/np.sqrt(q2); u/=np.linalg.norm(u)+1e-30
            ctx['fold']=d.force.copy(); ctx['xold']=d.pos.copy()
            u=ctx.get('u')
            if u is not None: d.vel-=g*(d.vel@u)*u
        P.append((f'mkill{g}',mkill))
    for b in (0.1,0.3):   # gentle persistent steering of v toward f (steer1.0 was catastrophic)
        def psteer(d,st,ctx,b=b):
            f=d.force; nf=np.linalg.norm(f); nv=np.linalg.norm(d.vel)
            if nf<1e-30 or nv<1e-30: return
            w=d.vel+b*f/nf*nv; nw=np.linalg.norm(w)
            if nw>1e-30: d.vel[:]=w*(nv/nw)
        P.append((f'psteer{b}',psteer))
    for k in (0.5,1.0):   # curved-valley extrapolation: extra displacement along turning normal
        def vcurve(d,st,ctx,k=k):
            s0=ctx.get('s0')
            fo=ctx.get('fold'); xo=ctx.get('xold')
            if xo is not None:
                s=d.pos-xo; ns=np.linalg.norm(s)
                if s0 is not None and ns>1e-30:
                    n=s/ns-s0
                    nn=np.linalg.norm(n)
                    if nn>1e-6: d.pos += k*0.5*ns*ns*n/nn
                if ns>1e-30: ctx['s0']=s/ns
            ctx['fold']=d.force.copy(); ctx['xold']=d.pos.copy()
        P.append((f'vcurve{k}',vcurve))
    # ---- gated composites: apply a register only when a dynamical feature favors it (switching test)
    P+= [('gmkill',make_gmkill(0.7,0.5)),('g2ts',make_g2ts(0.8,0.5))]
    return P

def _upd_u(d,ctx):           # shared: online transverse-mode register update (EMA, sign-aligned)
    fo=ctx.get('fold'); xo=ctx.get('xold')
    if fo is not None:
        s=d.pos-xo; ss=s@s
        if ss>1e-24:
            df=d.force-fo; q=df-(df@s/ss)*s; q2=q@q
            if q2>1e-12*(df@df):
                u=ctx.setdefault('u',np.zeros_like(q))
                sgn=1.0 if u@q>=0 else -1.0
                u[:]=0.7*u+0.3*sgn*q/np.sqrt(q2); u/=np.linalg.norm(u)+1e-30
def make_gmkill(cgate,g):    # transverse-mode kill only while cos(F,v)<cgate
    def gmkill(d,st,ctx):
        fi=d.force*d.invm; nf2=fi@fi; vv=d.vel@d.vel
        cos=(d.vel@fi)/np.sqrt(max(nf2*vv,1e-64))
        _upd_u(d,ctx)
        u=ctx.get('u')
        if cos<cgate and u is not None: d.vel-=g*(d.vel@u)*u
        ctx['fold']=d.force.copy(); ctx['xold']=d.pos.copy()
    return gmkill
def make_g2ts(cgate,g):      # slow-momentum boost only while cos(F,v)>cgate
    def g2ts(d,st,ctx):
        fi=d.force*d.invm; nf2=fi@fi; vv=d.vel@d.vel
        cos=(d.vel@fi)/np.sqrt(max(nf2*vv,1e-64))
        vs=ctx.setdefault('vs',np.zeros_like(d.vel))
        vs*=0.9; vs+=0.1*d.vel
        if cos>cgate: d.vel+=g*vs
    return g2ts

# ---------------- drivers per force field ----------------
def setup_mmff(system):   # same init contract as scan_all_ff.init_mmff
    os.chdir(REPO)
    from pyBall import MMFF as mmff
    mmff.setVerbosity(verbosity=0, idebug=0)
    bUFF=system.startswith('uff'); mol=system.split('_',1)[1]
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
    return mmff

def setup_eff():
    os.chdir(os.path.join(REPO,'tests','tEFF'))
    sys.path.insert(0, os.getcwd())
    import fire_study_lib as fs
    from pyBall import eFF as eff
    fs.set_atom_params(); fs.load_e('H2O_pairs_fc.xyz', verbosity=0)
    return fs, eff

# ---------------- main experiment ----------------
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--system', default='mmff_propandiol')
    ap.add_argument('--parents', type=int, default=20)
    ap.add_argument('--horizon', type=int, default=16)
    ap.add_argument('--budget',  type=int, default=600)
    ap.add_argument('--stride',  type=int, default=1)
    ap.add_argument('--persist', action='store_true', help='re-apply each action every step (policy branches, not kicks)')
    ap.add_argument('--policies', action='store_true', help='branch with register/memory policies (PyFire engine)')
    ap.add_argument('--strong',   action='store_true', help='use champion config (autoMass+massMetric) + param-switch actions')
    ap.add_argument('--qscan',    action='store_true', help='score actions by evals-to-convergence (not just ΔE@H) per §17')
    ap.add_argument('--rollout',  action='store_true', help='run each policy as a full trajectory to convergence; compare evals')
    ap.add_argument('--gatescan', action='store_true', help='scan (cos_gate, gain) grid for gated register policies on rollouts')
    ap.add_argument('--fconv',    type=float, default=1e-4, help='convergence |F| threshold for --rollout')
    ap.add_argument('--maxeval',  type=int,   default=5000, help='eval budget for --rollout')
    ap.add_argument('--verify', action='store_true', help='verify PyFire == C++ move_FIRE trajectory, then exit')
    ap.add_argument('--fire', default='', help='fire params k=v,k=v (e.g. falpha=0.964,fdec=0.837,dt=0.142)')
    ap.add_argument('--out', default='')
    a=ap.parse_args(); print('branch_lab:', vars(a), flush=True)

    ch = CHAMP[a.system] if a.strong else None
    if a.system.startswith('eff'):
        fs,eff=setup_eff(); p0=eff.pDOFs.copy()
        eff.initOpt(dt=0.03, damping=0.5, f_limit=100.0, bMass=False)
        bind_step_api(eff.lib)
        fp={'finc':1.1,'fdec':0.5,'falpha':0.8,'dt_max':0.03}; dt=0.03; damp=0.5; flim=100.0
        if ch: dt=ch['dt']; damp=ch['damp']; fp={'finc':ch['finc'],'fdec':ch['fdec'],'falpha':ch['falpha'],'dt_max':ch['dt']}
        if a.fire:
            for kv in a.fire.split(','):
                k,v=kv.split('='); k=k.strip(); v=float(v)
                if k=='dt': dt=v
                else: fp[k]=v
        drv=Driver(eff.lib, dt, damp, fp)
        if ch: eff.lib.setAutoMass(ch['am'],100.0); drv.setp('massMetric',ch['mm'])
        def reset(): eff.pDOFs[:]=p0; eff.vDOFs[:]=0
    else:
        mmff=setup_mmff(a.system); lib=mmff.lib
        bind_step_api(lib)
        p0=mmff.DOFs.copy()
        fp={'finc':1.1,'fdec':0.5,'falpha':0.8,'dt_max':1.0,'minLastNeg':5}; dt=0.1 if a.system.startswith('uff') else 0.2; damp=0.95; flim=10.0
        if ch: dt=ch['dt']; damp=ch['damp']; fp={'finc':ch['finc'],'fdec':ch['fdec'],'falpha':ch['falpha'],'dt_max':ch['dt'],'minLastNeg':5}   # minLastNeg=5: set_opt() semantics used by every meta_opt/scan eval
        if a.fire:
            for kv in a.fire.split(','):
                k,v=kv.split('='); k=k.strip(); v=float(v)
                if k=='dt': dt=v
                else: fp[k]=v
        drv=Driver(lib, dt, damp, fp)
        if ch: lib.setAutoMass(ch['am'],100.0); drv.setp('massMetric',ch['mm'])
        def reset(): mmff.DOFs[:]=p0; drv.vel[:]=0

    pf=PyFire(drv, dict(fp), flim)
    if a.verify:
        reset(); E0=drv.eval(); b0=drv.save()          # post-init complete state incl. force(p0)
        drv.restore(b0); pos_cpp=[]
        for i in range(30): drv.move(); drv.eval(); pos_cpp.append(drv.pos.copy())
        drv.restore(b0); E1=drv.eval(); pf.sync(); pos_py=[]
        for i in range(30): pf.step(); drv.eval(); pos_py.append(drv.pos.copy())
        dmax=max(float(np.abs(pc-pp).max()) for pc,pp in zip(pos_cpp,pos_py))
        print(f'verify PyFire vs move_FIRE: E0={E0:.10f} E1={E1:.10f} max|dpos| over 30 steps = {dmax:.3e}')
        sys.exit(0 if dmax<1e-12 else 1)

    if a.rollout or a.gatescan:
        reset(); drv.eval(); b0=drv.save()
        pols=list(make_policies())
        if a.gatescan:
            pols=[('base',None)]
            for cg in (0.5,0.7,0.9):
                for g in (0.25,0.5,1.0): pols.append((f'g2ts_c{cg}_g{g}',make_g2ts(cg,g)))
            for cg in (0.3,0.5,0.7):
                for g in (0.25,0.5,1.0): pols.append((f'gmk_c{cg}_g{g}',make_gmkill(cg,g)))
        print(f'=== rollouts: |F|<{a.fconv:g}, budget {a.maxeval}, {len(pols)} policies ===', flush=True)
        res=[]
        for name,fn in pols:
            drv.restore(b0); pf.sync(); ctx={}
            nE=0; f2=drv.f2(); E=drv.eval()   # eval count starts at 1 (the force at p0 is from the restore eval)
            nE=1
            while f2>a.fconv*a.fconv and nE<a.maxeval:
                pf.step(fn,ctx); E=drv.eval(); f2=drv.f2(); nE+=1
                if not np.isfinite(E): break
            res.append((name,nE,np.sqrt(f2),E))
            print(f'  {name:10s} evals={nE:5d} |F|={np.sqrt(f2):.3e} E={E:.6f} {"CONV" if f2<=a.fconv*a.fconv else "BUDGET/DIV"}', flush=True)
        res.append(('ref_cpp',-1,0,0))
        sys.exit(0)

    # ---- collect parent trajectory (stepwise FIRE = same loop as C++ run) ----
    snaps=[]; Es=[]; Fs=[]; poss=[]; fors=[]
    reset()
    E=drv.eval(); F2=drv.f2()
    for it in range(a.budget):
        snaps.append(drv.save()); Es.append(E); Fs.append(np.sqrt(F2))
        poss.append(drv.pos.copy()); fors.append(drv.force.copy())
        if F2 < 1e-12: break
        drv.move(); E=drv.eval(); F2=drv.f2()
    print(f'collected {len(snaps)} states, |F| {Fs[0]:.2e}->{Fs[-1]:.2e}', flush=True)

    idxs=list(range(2,len(snaps),max(1,(len(snaps)-2)//a.parents)))[:a.parents]
    print(f'{len(idxs)} parents, horizon={a.horizon}', flush=True)

    if a.qscan:
        # §17.7 test 1: score each action at H=16, H=64 (ΔE) AND by evals-to-convergence from the parent.
        acts=make_param_actions() if a.strong else (make_policies() if a.policies else make_actions())
        Bmax=a.maxeval; Hs=(16,64)
        qE=np.full((len(idxs),len(acts),2), np.nan)   # ΔE at H16,H64
        qN=np.full((len(idxs),len(acts)), np.nan)     # evals to |F|<fconv (nan=budget/diverged)
        for ip,ix in enumerate(idxs):
            s=poss[ix]-poss[ix-1]; ns=np.linalg.norm(s)
            u=None
            if ix>=2:
                df=fors[ix]-fors[ix-1]; dft=df-(df@s/(ns*ns))*s if ns>1e-30 else df
                nu=np.linalg.norm(dft)
                if nu>1e-30: u=dft/nu
            ctx0={'s': s/ns if ns>1e-30 else None, 'u': u}
            for ia,(name,fn) in enumerate(acts):
                drv.restore(snaps[ix]); ctxa=dict(ctx0)
                if not a.persist and not a.policies: fn(drv,ctxa)
                E=Es[ix]; f2=drv.f2(); conv=False
                if a.policies: pf.sync(); pctx={}
                for st in range(Bmax):
                    if a.policies:
                        pf.step(fn,pctx); E=drv.eval()
                    else:
                        if a.persist and ia>0: fn(drv,ctxa)
                        drv.move(); E=drv.eval()
                    f2=drv.f2()
                    if st+1 in Hs: qE[ip,ia,Hs.index(st+1)]=Es[ix]-E
                    if f2 < a.fconv*a.fconv: conv=True; break
                    if not np.isfinite(E): break
                if conv: qN[ip,ia]=st+1
            print(f'  parent {ip+1}/{len(idxs)} (ix={ix}) done', flush=True)
        base=qN[:,0]
        print(f'\n=== Q-scan: evals-to-convergence from parent (|F|<{a.fconv:g}, budget {Bmax}) ===')
        print(f'baseline remaining evals: mean={np.nanmean(base):.0f} med={np.nanmedian(base):.0f} nconv={np.isfinite(base).sum()}/{len(idxs)}')
        for ia,(name,_) in enumerate(acts):
            sav=base-qN[:,ia]
            both=np.isfinite(base)&np.isfinite(qN[:,ia])
            print(f'  {name:14s} evals mean={np.nanmean(qN[:,ia]):6.0f}  saved mean={np.nanmean(sav):+6.0f}  parents better={np.sum(both&(sav>0)):2d} worse={np.sum(both&(sav<0)):2d}  noconv={np.sum(~np.isfinite(qN[:,ia])):2d}')
        print('\n=== ΔE@16 vs ΔE@64 vs evals-saved (rank consistency) ===')
        for ia,(name,_) in enumerate(acts):
            e16=qE[:,ia,0]; e64=qE[:,ia,1]; sav=base-qN[:,ia]
            print(f'  {name:14s} ΔE16={np.nanmean(e16):+.4f} ΔE64={np.nanmean(e64):+.4f} evals-saved={np.nanmean(sav):+6.0f}')
        out=a.out or f'branch_results/{a.system}_qscan.npz'
        os.makedirs(os.path.dirname(out),exist_ok=True)
        np.savez(out, qE=qE, qN=qN, names=np.array([n for n,_ in acts]), idxs=np.array(idxs))
        print('saved', out, flush=True)
        sys.exit(0)

    acts=make_param_actions() if a.strong else (make_policies() if a.policies else make_actions())
    H=a.horizon
    nF=6
    feat=np.zeros((len(idxs),nF)); adv=np.zeros((len(idxs),len(acts)))
    emaF=0.0
    for ip,ix in enumerate(idxs):
        # context: last accepted step + transverse mode direction
        s=poss[ix]-poss[ix-1]; ns=np.linalg.norm(s)
        u=None
        if ix>=2:
            df=fors[ix]-fors[ix-1]; dft=df-(df@s/(ns*ns))*s if ns>1e-30 else df
            nu=np.linalg.norm(dft)
            if nu>1e-30: u=dft/nu
        ctx={'s': s/ns if ns>1e-30 else None, 'u': u}
        # features (restore fills the bound vel/force views)
        drv.restore(snaps[ix])
        f=fors[ix]; nf=np.linalg.norm(f); nv=np.linalg.norm(drv.vel)
        cosfv = float(f@drv.vel)/(nf*nv+1e-30)
        y=fors[ix]-fors[ix-1]; bb1=float(s@y)/(ns*ns+1e-30)
        emaF = nf if ip==0 else 0.9*emaF+0.1*nf
        feat[ip]=[cosfv, np.log10(nf+1e-30), nv/(nf+1e-30), np.clip(bb1,-1e6,1e6),
                  nf/(emaF+1e-30), (Es[ix-1]-Es[ix])/(nf+1e-30)]
        for ia,(name,fn) in enumerate(acts):
            drv.restore(snaps[ix])
            E=Es[ix]   # actions don't move positions -> same energy; only moves cost evals
            if a.policies:
                pf.sync(); pctx={}
                for st in range(H):
                    pf.step(fn,pctx); E=drv.eval()
                    if not np.isfinite(E): E=1e30; break
                    if drv.f2()<1e-12: break
            else:
                ctxa=dict(ctx)                        # per-action ctx so one-shot init flags (e.g. msInit) fire once
                if not a.persist: fn(drv,ctxa)        # one-shot action
                for st in range(H):
                    if a.persist and ia>0: fn(drv,ctxa)   # persistent policy: apply before each move (incl. first)
                    drv.move(); E=drv.eval()
                    if not np.isfinite(E): E=1e30; break
                    if drv.f2()<1e-12: break
            adv[ip,ia]=Es[ix]-E
        if (ip+1)%5==0: print(f'  parent {ip+1}/{len(idxs)} done', flush=True)

    # note: every action branch costs exactly H force evals (or fewer on early conv/fail)
    # checkpoint-equivalence sanity: 'base' continuation from ix must land on the recorded E[ix+H]
    nchk=min(len(idxs), 3)
    for ip in range(nchk):
        ix=idxs[ip]
        if ix+H < len(Es):
            print(f'ckpt-equiv parent ix={ix}: branch_E={Es[ix]-adv[ip,0]:.10f} recorded E[ix+H]={Es[ix+H]:.10f} diff={abs(Es[ix]-adv[ip,0]-Es[ix+H]):.2e}')
    base=adv[:,0]; oracle=adv.max(1); ibest=adv.argmax(1)
    gap=oracle-base
    print('\n=== teacher headroom (E_base - E_branch, positive=branch better, H=%i evals) ==='%H)
    print(f'baseline mean ΔE={base.mean():.3f} | oracle mean ΔE={oracle.mean():.3f} | headroom mean={gap.mean():.3f} med={np.median(gap):.3f} p90={np.percentile(gap,90):.3f} max={gap.max():.3f}')
    print(f'parents where some action beats baseline: {(gap>1e-6).mean()*100:.0f}%')
    print('\naction win-rate and mean ΔE:')
    for ia,(name,_) in enumerate(acts):
        print(f'  {name:10s} wins {(ibest==ia).mean()*100:4.0f}%  mean ΔE={adv[:,ia].mean():+.3f}')
    # predictability: ridge fit per action -> predicted best vs oracle
    X=np.c_[feat, np.ones(len(idxs))]
    lam=1e-4
    W=np.linalg.solve(X.T@X+lam*np.eye(nF+1), X.T@adv)
    pred=(X@W).argmax(1)
    cap=adv[np.arange(len(idxs)),pred]-base
    print(f'\npredictability: oracle gap mean={gap.mean():.3f} | feature-predicted gap mean={cap.mean():.3f} ({100*cap.mean()/max(gap.mean(),1e-9):.0f}% recovered)')
    out=a.out or f'branch_results/{a.system}.npz'
    os.makedirs(os.path.dirname(out),exist_ok=True)
    np.savez(out, feat=feat, adv=adv, names=np.array([n for n,_ in acts]), idxs=np.array(idxs),
             Es=np.array(Es)[idxs], Fs=np.array(Fs)[idxs])
    print('saved', out, flush=True)

if __name__=='__main__': main()
