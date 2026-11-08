"""Small offline meta-search over interpretable optimizer parameters (NumPy only).

This is NOT gradient-based meta-learning. Replace the outer random search with
CMA-ES / Nevergrad / Bayesian optimization only after the benchmark is stable.
Separate train/validation/test systems and seeds in real use.
"""
import argparse
import numpy as np
from bench import rotated_valley, rosenbrock
from valleyopt import METHODS, minimize

SPACES = {
    'momentum': dict(beta=(.1,.95), explore=(.0,1.),grow=(1.01,1.3),shrink=(.25,.85),step=(.015,.30)),
    'secant': dict(gain=(.2,1.5),forget=(0.,.15),explore=(0.,1.),grow=(1.01,1.30),shrink=(.25,.85),step=(.015,.3)),
    'curved': dict(smooth=(.1,.95),exploration=(.0,1.),grow=(1.01,1.30),shrink=(.25,.85),step=(.015,.3)),
    'wall': dict(smooth=(.1,.95),exploration=(.0,1.),probe_ratio=(.25,1.2),grow=(1.01,1.30),shrink=(.25,.85),step=(.015,.3)),
    'ridge': dict(probe_frac=(.15,.85),grow=(1.01,1.5),shrink=(.30,.85),ridge=(1.e-7,1.e-2),step=(.015,.3)),
}

def sample(rng, space):
    # Mostly linear bounds; use log-space for strictly positive wide-ranging scales.
    d={}
    for k,(a,b) in space.items():
        d[k]=float(np.exp(rng.uniform(np.log(a),np.log(b)))) if a>0 and b/a > 1.e3 else float(rng.uniform(a,b))
    return d

def score(method, params, problems, neval=200, seed=123):
    losses=[]
    for i,(f,x0,target) in enumerate(problems):
        s=minimize(f, method, x0, max_evals=neval, seed=seed+i, **params)
        # Target energy must be known (synthetic benchmark), or from independent
        # high-quality reference calculations. Normalize objective energy scale.
        scale=max(abs(f(x0)-target),1.e-8)
        regret=np.maximum(np.asarray(s.trace.best_energy)-target,0.)/scale
        losses.append(float(np.mean(np.log1p(30*regret)) + 1.0*np.log1p(30*regret[-1])))
    return float(np.mean(losses))

def run(method='ridge', trials=30, neval=200, seed=123):
    rng=np.random.default_rng(seed)
    problems=[(*rotated_valley(n=n,bend=b,seed=100003+seed+17*j),0.)
              for j,(n,b) in enumerate([(2,.2),(3,.3),(6,.15),(12,.25)])]
    problems.append((rosenbrock,np.array([-1.2,1.]),0.))
    best=(float('inf'),None)
    for i in range(trials):
        p=sample(rng,SPACES[method])
        loss=score(METHODS[method],p,problems,neval=neval,seed=seed)
        if loss<best[0]:
            best=(loss,p)
            print('best',i,'loss',round(loss,6),p,flush=True)
    print('FINAL:',best)
    return best

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--method',choices=list(METHODS),default='ridge')
    p.add_argument('--trials',type=int,default=30)
    p.add_argument('--neval',type=int,default=200)
    p.add_argument('--seed',type=int,default=123)
    run(**vars(p.parse_args()))
