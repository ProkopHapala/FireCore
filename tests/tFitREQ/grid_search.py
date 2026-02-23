#!/usr/bin/env python3
import os
import sys
import numpy as np
import itertools
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import FitREQ as fit
import barrier_utils as bu

def run_grid_search():
    # Use the Z-file (double-well) as required
    xyz_path = "../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz"
    # Epairs-only DOF selection (no charges, no atom H-corrections)
    base_dof_file = "dofSelection_MorseSR_H2O_epairOnly.dat"
    out_dir = "grid_search_out"
    os.makedirs(out_dir, exist_ok=True)
    
    kMorse_vals       = [1.6, 1.7, 1.8]
    #kMorse_vals       = [1.7]
    Lepairs_vals      = [0.6, 0.8, 1.0, 1.2]
    #Lepairs_vals      = [1.0]
    #weight_alpha_vals = [1.0, 2.0, 3.0, 4.0]
    weight_alpha_vals = [ 4.0]

    hScale_vals      = [0.0,0.5,1.0]
    
    # Parse reference mapping
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(xyz_path)
    
    # Read XYZ data for weights
    Erefs, x0s = fit.read_xyz_data(xyz_path)
    
    base_dofs, _ = bu.parse_dof_file(base_dof_file)
    #bu.write_dof_file("tmp_grid_dofs_base.dat", base_dofs)
    dof_names = bu.setup_fitreq(base_dof_file, xyz_path, imodel=7, Regularize=1)
    #dof_names = bu.setup_fitreq("tmp_grid_dofs_base.dat", xyz_path, imodel=7, Regularize=1)
    # Silence C++/Python verbosity to avoid slow/large logs
    fit.setVerbosity(1, 0, 0, 0, 0, 0, 0)
    
    results = []
    
    print(f"Starting grid search. Total runs: {len(kMorse_vals)*len(Lepairs_vals)*len(weight_alpha_vals)*len(hScale_vals)}")
    
    run_idx = 0
    for kM, Lep, w_alpha, hScale in itertools.product(kMorse_vals, Lepairs_vals, weight_alpha_vals,hScale_vals):
        print(f"Run {run_idx}: kMorse={kM}, Lepairs={Lep}, weight_alpha={w_alpha}, hScale={hScale}")
        
        # Reset DOFs to base
        dofs = [dict(d) for d in base_dofs]
        bu.write_dof_file("tmp_grid_dofs.dat", dofs)
        fit.loadDOFSelection("tmp_grid_dofs.dat")
        
        # Set globals
        fit.setGlobalParams(kMorse=kM, Lepairs=Lep, softClamp_start=4.0, softClamp_max=6.0, hScale=hScale)
        
        # Set weights
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=w_alpha)
        weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
        fit.setWeights(weights0)
        
        # Relax
        fit.setVerbosity(0, 0, 0, 0, 0, 0, 0)
        fit.run(iparallel=0, ialg=1, nstep=150, Fmax=1e-8, dt=0.05, damping=0.0, max_step=0.1, bClamp=True)
        fit.getBuffs()
        
        # Get final DOFs and Error
        final_dofs = np.array(fit.DOFs)
        Eerr, Es, Fs = fit.getEs(bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False)
        
        # Compute model grid for plotting
        run_params = dict(nstep=150, Fmax=1e-8, dt=0.05, max_step=0.1, damping=0.0)
        Gmodel = fit.compute_model_grid(xyz_path, seq, Gref.shape, do_fit=False, bAddEpairs=True, run_params=run_params, bOutXYZ=False)
        
        # Plot (single combined figure only)
        prefix = f"{out_dir}/run_{run_idx}_kM{kM}_L{Lep}_wa{w_alpha}_hS{hScale}"
        title = f"kM={kM} L={Lep} wa={w_alpha} hS={hScale} Err={Eerr:.4e}"
        params_lines = [
            f"kMorse={kM:.3f}",
            f"Lepairs={Lep:.3f}",
            f"hScale={hScale:.3f}",
            f"w_alpha={w_alpha:.3f}",
        ]
        params_lines += [f"{name}={val:.4f}" for name, val in zip(dof_names, final_dofs)]
        params_text = "\n".join(params_lines)

        fit.plot_compare_combined(Gref, Gmodel, angles, distances, title, save_path=f"{prefix}.png", kcal=True, params_text=params_text)
        
        # Save DOFs
        for i, d in enumerate(dofs):
            d['xstart'] = final_dofs[i]
        
        hyper_meta = {'kMorse': kM, 'Lepairs': Lep, 'weight_alpha': w_alpha, 'hScale': hScale}
        bu.write_dof_file(f"{prefix}.dat", dofs, hyper=hyper_meta)
        
        results.append({
            'kMorse': kM,
            'Lepairs': Lep,
            'weight_alpha': w_alpha,
            'error': Eerr,
            'dofs': final_dofs
        })
        
        run_idx += 1
        
    # Save ensemble
    np.savez(f"{out_dir}/ensemble_data.npz", 
             kMorse=np.array([r['kMorse'] for r in results]),
             Lepairs=np.array([r['Lepairs'] for r in results]),
             weight_alpha=np.array([r['weight_alpha'] for r in results]),
             error=np.array([r['error'] for r in results]),
             dofs=np.array([r['dofs'] for r in results]),
             dof_names=dof_names)
             
    print("Grid search complete.")

if __name__ == "__main__":
    run_grid_search()
