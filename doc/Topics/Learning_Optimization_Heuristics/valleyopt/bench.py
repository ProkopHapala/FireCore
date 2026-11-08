"""Small reproducible smoke benchmark (not evidence of superiority)."""
import argparse
import numpy as np
from valleyopt import minimize, METHODS

def rotated_valley(n=12, bend=.20, k_soft=0.15, k_hard=15., seed=123):
    """2D Rosenbrock-like curved centerline embedded in n dims via rotation.
    Energy minimum y=[2, 0, 0, ...], valley y1=bend*(y0**2-4).
    Rotation makes the identity of 'soft coordinates' unavailable to optimizer.
    """
    if n < 2: raise ValueError('n must be at least 2')
    # Two random orthonormal vectors, O(N) memory and arithmetic.
    # Avoid allocating an NxN rotation matrix even for N=1_000_000.
    rng=np.random.default_rng(seed)
    u=rng.standard_normal(n); u/=np.linalg.norm(u)
    v=rng.standard_normal(n); v-=np.dot(v,u)*u; v/=np.linalg.norm(v)
    x0=-2.*u
    def energy(x):
        y0=float(np.dot(u,x)); y1=float(np.dot(v,x))
        bath=max(0.,float(np.dot(x,x))-y0*y0-y1*y1)
        return .5*k_soft*(y0-2.)**2 + .5*k_hard*(y1-bend*(y0**2-4))**2 + .5*k_hard*bath
    return energy, x0

def rosenbrock(x):
    x=np.asarray(x)
    return float(np.sum(100*(x[1:]-x[:-1]**2)**2+(1-x[:-1])**2))

def run(n=12, neval=600, seed=3):
    tests={'rotated_valley':rotated_valley(n=n,seed=seed+10003),
           'rosenbrock':(rosenbrock,np.array([-1.2,1.0]))}
    for test,(f,x0) in tests.items():
        E0=f(x0)
        print(f'\n{test}: N={len(x0)} E0={E0:.8g}')
        for name,method in METHODS.items():
            sol=minimize(f,method,x0,max_evals=neval,step=.10,seed=seed)
            print(f'  {name:9s} E_inc={sol.e:.7g} E_best={sol.best_e:.7g} '
                  f'nfev={sol.trace.nfev:4} accepts={sol.trace.accepted:3} step={sol.h:.2g}')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=12)
    p.add_argument('--neval',type=int,default=600);p.add_argument('--seed',type=int,default=3)
    run(**vars(p.parse_args()))
