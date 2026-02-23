import numpy as np
import os
import sys

sys.path.append("../../")
from pyBall import FitREQ as fit

def parse_dof_file(fname):
    dofs = []
    hyper = {}
    with open(fname, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#HYPER'):
                parts = stripped.split()[1:]
                for p in parts:
                    if '=' in p:
                        k, v = p.split('=', 1)
                        hyper[k] = float(v)
                continue
            if stripped.startswith('#'):
                continue
            parts = stripped.split()
            if len(parts) >= 10:
                dofs.append({
                    'typename': parts[0],
                    'comp': int(parts[1]),
                    'Min': float(parts[2]),
                    'Max': float(parts[3]),
                    'xlo': float(parts[4]),
                    'xhi': float(parts[5]),
                    'Klo': float(parts[6]),
                    'Khi': float(parts[7]),
                    'K0': float(parts[8]),
                    'xstart': float(parts[9]),
                    'invMass': float(parts[10]) if len(parts)>10 else 1.0
                })
    return dofs, hyper

def write_dof_file(fname, dofs, hyper=None):
    with open(fname, 'w') as f:
        if hyper:
            kv = " ".join(f"{k}={v}" for k, v in hyper.items())
            f.write(f"#HYPER {kv}\n")
        f.write("# typename comp Min Max xlo xhi Klo Khi K0 xstart invMass\n")
        for d in dofs:
            f.write(f"{d['typename']} {d['comp']} {d['Min']} {d['Max']} {d['xlo']} {d['xhi']} {d['Klo']} {d['Khi']} {d['K0']} {d['xstart']} {d['invMass']}\n")

def get_dof_values(dofs):
    return np.array([d['xstart'] for d in dofs])

def set_dof_values(dofs, values):
    for i, d in enumerate(dofs):
        d['xstart'] = values[i]

TYPES_LOADED = False

def setup_fitreq(dof_file, xyz_path, kMorse=1.7, Lepairs=0.9, hScale=0.0, softClamp=True, sc_start=4.0, sc_max=6.0, imodel=7, Regularize=1):
    global TYPES_LOADED
    fit.setVerbosity(0)
    fit.setup(imodel=imodel, EvalJ=1, WriteJ=1, Regularize=Regularize, SaveJustElementXYZ=-1)
    if softClamp:
        fit.setGlobalParams(kMorse=kMorse, Lepairs=Lepairs, hScale=hScale, softClamp_start=sc_start, softClamp_max=sc_max)
    else:
        fit.setGlobalParams(kMorse=kMorse, Lepairs=Lepairs, hScale=hScale)
    if not TYPES_LOADED:
        fit.loadTypes()
        TYPES_LOADED = True
    fit.loadDOFSelection(fname=dof_file)
    fit.loadXYZ(xyz_path, bAddEpairs=True)
    fit.getBuffs()
    return fit.loadDOFnames(dof_file)

def evaluate_current(bFs=False):
    # Evaluates energy (and optionally gradients) at current fit.DOFs
    Eerr, Es, Fs = fit.getEs(bOmp=False, bDOFtoTypes=True, bEs=True, bFs=bFs)
    return Eerr, Es, Fs

def relax_tethered(dofs, kMorse=1.7, Lepairs=0.9, nstep=100, Fmax=1e-8, dt=0.05, damping=0.0):
    # Write temp file to update constraints
    tmp_dof = "_tmp_dofs_tether.dat"
    write_dof_file(tmp_dof, dofs)
    fit.loadDOFSelection(tmp_dof)
    fit.getBuffs()
    
    # Run
    fit.run(iparallel=0, ialg=1, nstep=nstep, Fmax=Fmax, dt=dt, damping=damping, max_step=0.1, bClamp=True)
    fit.getBuffs()
    
    # Update dofs with relaxed values
    vals = np.array(fit.DOFs)
    set_dof_values(dofs, vals)
    return vals

def compute_fd_hessian(dofs, h=1e-4):
    n = len(dofs)
    H = np.zeros((n, n))
    base_vals = get_dof_values(dofs)
    
    for i in range(n):
        # +h
        v_plus = base_vals.copy()
        v_plus[i] += h
        for j in range(n): fit.DOFs[j] = v_plus[j]
        _, _, F_plus = evaluate_current(bFs=True)
        
        # -h
        v_minus = base_vals.copy()
        v_minus[i] -= h
        for j in range(n): fit.DOFs[j] = v_minus[j]
        _, _, F_minus = evaluate_current(bFs=True)
        
        # Central difference for gradient
        # Fs are negative gradients (forces) typically, let's assume they are gradients here, 
        # actually C++ evalFitError accumulates fDOFs which are gradients.
        H[i, :] = (F_plus - F_minus) / (2 * h)
        
    # Restore base values
    for j in range(n): fit.DOFs[j] = base_vals[j]
    
    # Symmetrize
    H = (H + H.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(H)
    return H, eigvals, eigvecs
