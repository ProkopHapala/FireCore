"""Experimental derivative-free optimizers. Only NumPy; no problem-specific features.

API: optimizer = Class(x0, e0, ...); repeatedly call x_trial = optimizer.ask(),
     e_trial = energy(x_trial), optimizer.tell(e_trial).
Only scalar energies enter the optimization algorithms. All accepted configurations
are monotone in E (direction exploration still evaluates uphill trials).

These are RESEARCH STARTERS, not implementations of OMoRF, CMA-ES, or Powell.
All O(N) memory except N*k in RidgeTR. Coordinate metric and constraints should
be handled in the caller's adapter: x is an unconstrained, scaled Euclidean vector.
"""
from dataclasses import dataclass, field
import numpy as np

EPS = 1e-14

def norm(x): return float(np.linalg.norm(x))
def unit(x):
    n = norm(x)
    return x/n if n > EPS else np.zeros_like(x)

def random_orth(rng, shape, basis=()):
    """Random isotropic direction orthogonal to existing vectors; O(N*k)."""
    v = rng.standard_normal(shape)
    for _ in range(2):               # reorthogonalize for numerical stability
        for u in basis: v -= np.dot(v, u)*u
    nv = norm(v)
    if nv < EPS: return None
    return v/nv

@dataclass
class Trace:
    nfev: int = 1
    accepted: int = 0
    best_energy: list = field(default_factory=list)
    events: list = field(default_factory=list)

class Base:
    def __init__(self, x0, e0, step=0.1, seed=0, min_step=1.e-9, max_step=10.):
        self.x = np.array(x0, dtype=float, copy=True).ravel()
        self.e = float(e0)
        self.h = float(step)
        self.hmin = float(min_step)
        self.hmax = float(max_step)
        self.rng = np.random.default_rng(seed)
        self.trace = Trace(best_energy=[self.e])
        self.best_e = self.e
        self.best_x = self.x.copy()
        self._pending = None
        self._displacement = None
        if not (np.isfinite(self.e) and self.hmin > 0 and self.hmax >= self.h >= self.hmin):
            raise ValueError('Invalid initial energy or step bounds')

    def _propose(self, d, kind='step'):
        if self._pending is not None: raise RuntimeError('tell() required before next ask()')
        d = np.asarray(d, dtype=float).ravel()
        if d.shape != self.x.shape or not np.all(np.isfinite(d)) or norm(d) < EPS:
            raise ValueError('Invalid displacement')
        trial = self.x + d
        self._pending = (kind, trial, d, self.e)
        return trial.copy()

    def _finish(self, energy):
        if self._pending is None: raise RuntimeError('ask() required before tell()')
        kind, trial, d, e0 = self._pending
        e1 = float(energy)
        if not np.isfinite(e1): e1 = float('inf')  # invalid sample treated as rejection
        accepted = e1 < e0 and not kind.startswith("probe")
        if e1 < self.best_e:
            self.best_e = e1
            self.best_x = trial.copy()
        self.trace.nfev += 1
        if accepted:
            self.x[:] = trial
            self.e = e1
            self.trace.accepted += 1
        self.trace.best_energy.append(self.best_e)
        self.trace.events.append((kind, bool(accepted), e1, self.h))
        self._pending = None
        return kind, d, e0, e1, accepted

class MomentumWalk(Base):
    """Baseline: successful step momentum + orthogonal exploration + reverse trials.

    Vector state: one direction t. No gradients or local quadratic assumption.
    """
    def __init__(self, x0, e0, *, beta=.70, explore=.25, grow=1.12,
                 shrink=.60, **kw):
        super().__init__(x0, e0, **kw)
        self.beta, self.explore = beta, explore
        self.grow, self.shrink = grow, shrink
        self.t = None
        self.rejected = 0
        self.mirror = None

    def ask(self):
        if self.mirror is not None:
            d = self.mirror; self.mirror = None
            return self._propose(d, 'mirror')
        if self.t is None:
            v = random_orth(self.rng, self.x.shape)
            return self._propose(self.h*v, 'discover')
        v = random_orth(self.rng, self.x.shape, [self.t])
        if v is None: v = np.zeros_like(self.t)
        a = min(2., self.explore*(1 + self.rejected))
        d = self.h*unit(self.t + a*v)
        return self._propose(d, 'momentum')

    def tell(self, energy):
        kind, d, e0, e1, ok = self._finish(energy)
        if ok:
            u = unit(d)
            self.t = u if self.t is None else unit(self.beta*self.t+(1-self.beta)*u)
            self.rejected = 0
            self.h = min(self.hmax, self.h*self.grow)
        else:
            self.rejected += 1
            self.h = max(self.hmin, self.h*self.shrink)
            if self.t is None and kind == 'discover':
                self.mirror = -d                  # same anchor, opposite probe
            if self.rejected >= 6: self.t = None  # restart direction discovery

class SecantSketch(Base):
    """Online finite-difference gradient sketch from *every* trial: normalized LMS.

    Each observation (d, deltaE) constrains one gradient projection, even on a
    rejected move. Vector state: g_hat, previous direction. No history matrix.
    Nonlinearity bias is controlled by shrink/trust radius, not removed.
    """
    def __init__(self, x0, e0, *, gain=.9, forget=.015,
                 explore=.30, grow=1.10, shrink=.65, **kw):
        super().__init__(x0, e0, **kw)
        self.gain, self.forget = gain, forget
        self.explore, self.grow, self.shrink = explore, grow, shrink
        self.g = np.zeros_like(self.x)
        self.rejected = 0
        self.last_d = None
        self.mirror = False

    def ask(self):
        if self.mirror and self.last_d is not None:
            self.mirror = False
            return self._propose(-self.last_d, 'mirror')
        gnorm = norm(self.g)
        if gnorm < 1.e-10 or self.rejected > 5:
            v = random_orth(self.rng, self.x.shape)
            if self.rejected > 5: self.g *= .5
            return self._propose(self.h*v, 'explore')
        t = -self.g/gnorm
        v = random_orth(self.rng, self.x.shape, [t])
        if v is None: v = np.zeros_like(t)
        d = self.h*unit(t + self.explore*(1+.3*self.rejected)*v)
        return self._propose(d, 'sketch')

    def tell(self, energy):
        kind, d, e0, e1, ok = self._finish(energy)
        if np.isfinite(e1):
            # Residual is computed against an estimate local to the OLD center.
            resid = (e1-e0) - np.dot(self.g, d)
            self.g *= (1-self.forget)
            self.g += self.gain*resid/(np.dot(d,d)+EPS)*d
        if ok:
            self.h = min(self.hmax, self.h*self.grow)
            self.rejected = 0
        else:
            self.rejected += 1
            self.h = max(self.hmin, self.h*self.shrink)
            # For exploratory measurements use an opposite probe on same anchor.
            if kind == 'explore' and not self.mirror:
                self.last_d = d.copy(); self.mirror = True

class CurvedWalk(Base):
    """Predict a turning valley from successive accepted displacements.

    Stores tangent t, turning normal n, curvature kappa, and rejection count.
    Kappa is geometrical path curvature, NOT an energy Hessian coefficient.
    """
    def __init__(self, x0, e0, *, smooth=.45, exploration=.30,
                 grow=1.10, shrink=.60, **kw):
        super().__init__(x0, e0, **kw)
        self.t = None; self.normal = None; self.kappa = 0.
        self.smooth, self.exploration = smooth, exploration
        self.grow, self.shrink = grow, shrink
        self.rejected = 0; self.mirror = None

    def ask(self):
        if self.mirror is not None:
            d = self.mirror; self.mirror = None
            return self._propose(d, 'mirror')
        if self.t is None:
            return self._propose(self.h*random_orth(self.rng,self.x.shape), 'discover')
        d = self.h*self.t
        if self.normal is not None:
            d += .5*self.kappa*self.h**2*self.normal
        if self.rejected:
            vs = [self.t] + ([self.normal] if self.normal is not None else [])
            v = random_orth(self.rng, self.x.shape, vs)
            if v is None: v = self.normal if self.normal is not None else -self.t
            d += self.exploration*self.h*min(self.rejected,5)*v
        return self._propose(d, 'extrapolate' if self.rejected == 0 else 'redirect')

    def tell(self, energy):
        kind, d, e0, e1, ok = self._finish(energy)
        if ok:
            tnew = unit(d)
            if self.t is not None:
                c = np.clip(np.dot(self.t,tnew),-1,1)
                angle = float(np.arccos(c))
                turn = tnew - c*self.t
                if norm(turn) > 1.e-8:
                    nnew = unit(turn)
                    self.normal = nnew if self.normal is None else unit((1-self.smooth)*self.normal+self.smooth*nnew)
                    self.normal = unit(self.normal-np.dot(self.normal,tnew)*tnew)
                    if norm(self.normal)<EPS: self.normal=None
                    k = angle/max(norm(d), EPS)
                    self.kappa = (1-self.smooth)*self.kappa+self.smooth*k
                else:
                    self.kappa *= (1-self.smooth)
            self.t = tnew; self.rejected = 0
            self.h = min(self.hmax, self.h*self.grow)
        else:
            self.rejected += 1
            self.h = max(self.hmin, self.h*self.shrink)
            if self.t is None and kind == 'discover': self.mirror = -d
            if self.rejected >= 7:
                self.t = None; self.normal = None; self.kappa = 0.

class RidgeTR(Base):
    """Small k-dimensional quadratic trust-region model in a MOVING subspace.

    This is a deliberately simplified reduced-model prototype, NOT OMoRF:
      * initializes tangent from an accepted exploratory move
      * finite differences 2k+ k(k-1)/2 at one fixed anchor
      * model predicts many successive steps (g <- g+H*a) without rebuilding
      * trust ratio checks model; bad fit restarts probes / changes basis
    State: U [N,k], g[k], H[k,k], k^2 scalars. No long history.
    """
    def __init__(self, x0, e0, *, rank=2, probe_frac=.45,
                 grow=1.25, shrink=.50, ridge=1.e-5, **kw):
        super().__init__(x0, e0, **kw)
        self.rank = min(max(1,int(rank)),self.x.size)
        self.probe_frac = probe_frac
        self.grow,self.shrink,self.ridge = grow,shrink,ridge
        self.t = None
        self.U = None
        self.g = None; self.H = None
        self.tasks = []             # tiny list of (kind, tuple(indices), unit coords)
        self.probes = {}
        self.rejected = 0
        self.rho = np.nan
        self._current_a = None; self._prediction = None

    def _new_basis(self, t=None):
        cols = []
        if t is not None and norm(t)>EPS: cols.append(unit(t.copy()))
        while len(cols)<self.rank:
            v = random_orth(self.rng, self.x.shape, cols)
            if v is None: break
            cols.append(v)
        self.U = np.column_stack(cols)
        self._reset_model()

    def _reset_model(self):
        self.g = None; self.H = None; self.probes.clear()
        self.tasks = []
        for i in range(self.U.shape[1]):
            self.tasks.extend([('plus',(i,)),('minus',(i,))])
        for i in range(self.U.shape[1]):
            for j in range(i+1,self.U.shape[1]): self.tasks.append(('cross',(i,j)))
        self._probe_step = self.probe_frac*self.h

    def _probe_offset(self, task):
        typ, ii = task
        if typ=='plus': return self._probe_step*self.U[:,ii[0]]
        if typ=='minus': return -self._probe_step*self.U[:,ii[0]]
        return self._probe_step*(self.U[:,ii[0]]+self.U[:,ii[1]])

    def _fit(self):
        k=self.U.shape[1]; h=self._probe_step
        g=np.zeros(k); H=np.zeros((k,k)); E=self.e
        for i in range(k):
            ep=self.probes[('plus',(i,))]; em=self.probes[('minus',(i,))]
            g[i]=(ep-em)/(2*h)
            H[i,i]=(ep+em-2*E)/(h*h)
        for i in range(k):
            for j in range(i+1,k):
                epq=self.probes[('cross',(i,j))]
                H[i,j]=H[j,i]=(epq-self.probes[('plus',(i,))]-self.probes[('plus',(j,))]+E)/(h*h)
        self.g, self.H = g,H

    def ask(self):
        if self.U is None:
            return self._propose(self.h*random_orth(self.rng,self.x.shape), 'discover')
        if self.tasks:
            task=self.tasks[0]
            return self._propose(self._probe_offset(task), 'probe')
        # Regularized local Newton/trust-region step in k dimensions.
        eig, V=np.linalg.eigh(self.H)
        # Clamp negative/unreliable eigenvalues (not a global-min guarantee).
        scale=max(float(np.max(np.abs(eig))),1.e-8)
        eig=np.maximum(eig, self.ridge*scale)
        a=-V@((V.T@self.g)/eig)
        an=norm(a)
        if an<self.hmin*0.01 or not np.all(np.isfinite(a)):
            # exhausted known subspace: trigger new exploration basis
            self._new_basis(self.t)
            return self.ask()
        if an>self.h: a*=self.h/an
        self._current_a=a
        self._prediction=float(np.dot(self.g,a)+.5*a@self.H@a)
        return self._propose(self.U@a,'model')

    def tell(self, energy):
        kind,d,e0,e1,ok=self._finish(energy)
        if kind=='discover':
            if ok:
                self.t=unit(d)
                self._new_basis(self.t)
                self.rejected=0
            else:
                # Mirror initial random direction rather than redrawing indefinitely
                self.t=-unit(d)
                self._new_basis(self.t)
                self.h=max(self.hmin,self.h*self.shrink)
        elif kind=='probe':
            task=self.tasks.pop(0)
            if not np.isfinite(e1):
                # Invalid physical DOFs or overflowing forcefield: don't fit
                # a non-finite quadratic. Retry at a smaller trust radius.
                self.h=max(self.hmin,self.h*self.shrink)
                self._reset_model()
                return
            self.probes[task]=e1
            # Anchor remains fixed during probes; best evaluated point is saved.
            if not self.tasks: self._fit()
        elif kind=='model':
            a=self._current_a
            pred=self._prediction
            actual=e1-e0
            rho=actual/pred if np.isfinite(e1) and pred < -EPS else float('-inf')
            self.rho=rho
            if ok:
                self.t=unit(d)
                self.g=self.g+self.H@a
                self.rejected=0
                if rho>.75: self.h=min(self.hmax,self.h*self.grow)
                elif rho<.25: self.h=max(self.hmin,self.h*self.shrink)
                if (rho<.25 or rho>2.5): self._reset_model()
            else:
                self.rejected+=1
                self.h=max(self.hmin,self.h*self.shrink)
                if self.rejected>=2:
                    self._new_basis(self.t)
                else: self._reset_model()


class WallProbeWalk(CurvedWalk):
    """Curved valley extrapolation + ± probes through a candidate canyon wall.

    If the forward move fails, evaluate energy at y +/- q*v where y is the
    rejected point, v is an estimated normal (or a new random orth direction).
    Fit a 1D parabola E(y+z*v) to infer steering correction z; evaluate it once.
    State: t, n, one probe direction v, three energies. No history arrays.
    This is a new experimental heuristic, NOT a reimplementation of a paper.
    """
    def __init__(self, x0, e0, *, probe_ratio=.7, max_wall_trials=2, **kw):
        super().__init__(x0,e0,**kw)
        self.probe_ratio = float(probe_ratio)
        self.max_wall_trials = int(max_wall_trials)
        self.phase = 'forward'
        self.wall_trials = 0
        self.wall_x = None
        self.wall_e = None
        self.probe_dir = None
        self.probe_q = 0.
        self.eplus = None; self.eminus = None
        self.correct_d = None

    def ask(self):
        if self.phase == 'plus':
            return self._propose(self.wall_x+self.probe_q*self.probe_dir-self.x,'probe_plus')
        if self.phase == 'minus':
            return self._propose(self.wall_x-self.probe_q*self.probe_dir-self.x,'probe_minus')
        if self.phase == 'correct':
            return self._propose(self.correct_d,'wall_correct')
        return super().ask()

    def tell(self, energy):
        if self._pending is None: raise RuntimeError('ask() required before tell()')
        kind,trial,d,eold=self._pending
        if kind in ('probe_plus','probe_minus'):
            # These are informative probes; DO NOT move the incumbent while
            # doing finite differences about the previously rejected point.
            _,_,_,e1,_ = self._finish(energy)
            if kind=='probe_plus':
                self.eplus=e1
                self.phase='minus'
            else:
                self.eminus=e1
                q=self.probe_q
                if np.isfinite(self.eplus) and np.isfinite(e1) and np.isfinite(self.wall_e):
                    gv=(self.eplus-e1)/(2*q)
                    kv=(self.eplus+e1-2*self.wall_e)/(q*q)
                    # Fit a parabola only for positive curvature. For negative
                    # curvature prefer the lower-energy endpoint.
                    if kv>1.e-9: z=float(np.clip(-gv/kv,-2*q,2*q))
                    else: z=q if self.eplus<e1 else -q
                    self.correct_d=(self.wall_x + z*self.probe_dir)-self.x
                    # z may round to 0; retry with a new orthogonal direction.
                    self.phase='correct' if norm(self.correct_d)>EPS and abs(z)>q*1.e-5 else 'forward'
                else: self.phase='forward'
                if self.phase=='forward':
                    self.h=max(self.hmin,self.h*.7)
            return
        if kind == 'wall_correct':
            self.phase='forward'
            before=self.e
            super().tell(energy)
            if self.e<before:
                self.wall_trials=0
                return
            self.wall_trials+=1
            # failed wall correction: try another transverse direction or shrink
            if self.wall_trials<self.max_wall_trials and self.t is not None:
                self._setup_wall()
            else:
                self.wall_trials=0
                self.phase='forward'; self.h=max(self.hmin,self.h*.7)
            return
        # ordinary forward trial (includes initial mirror discovery)
        was_forward=kind in ('extrapolate','redirect')
        before=self.e
        # skip parent failed-step redirection until probing gives information
        if was_forward and (not np.isfinite(float(energy)) or float(energy)>=before) and self.t is not None:
            self._finish(energy)     # record the failure without side effects
            self.rejected+=1
            self.wall_x=trial.copy()
            self.wall_e=float(energy)
            self._setup_wall()
            return
        super().tell(energy)
        if self.e<before: self.wall_trials=0

    def _setup_wall(self):
        if self.t is None or self.wall_x is None:
            self.phase='forward'; return
        v=None
        if self.normal is not None:
            w=self.normal-np.dot(self.normal,self.t)*self.t
            if norm(w)>EPS and self.wall_trials==0: v=unit(w)
        if v is None: v=random_orth(self.rng,self.x.shape,[self.t])
        if v is None:
            self.phase='forward'; self.h=max(self.hmin,.7*self.h); return
        self.probe_dir=v
        self.probe_q=max(self.hmin*10.,self.probe_ratio*self.h)
        self.phase='plus'


def minimize(energy, method, x0, max_evals=1000, **kwargs):
    """Convenience sequential driver; an integration can use ask/tell directly.

    Args:
      energy: callable on a COPY of x; must set every DOF each call (no drift).
      method: class (MomentumWalk, SecantSketch, CurvedWalk, RidgeTR).
    """
    x0=np.asarray(x0, dtype=float).ravel()
    solver=method(x0, float(energy(x0.copy())), **kwargs)
    while solver.trace.nfev < max_evals and solver.h > solver.hmin*1.001:
        x_trial=solver.ask()
        try: e=float(energy(x_trial))
        except (FloatingPointError,ValueError,OverflowError): e=float('inf')
        solver.tell(e)
    return solver

METHODS = {
    'momentum': MomentumWalk,
    'secant': SecantSketch,
    'curved': CurvedWalk,
    'ridge': RidgeTR,
    'wall': WallProbeWalk,
}
