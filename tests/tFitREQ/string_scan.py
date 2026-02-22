#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import FitREQ as fit
import barrier_utils as bu

def plot_string_results(alphas, energies_init, energies_relax, chain_relaxed, gradients, param_names, out_prefix):
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # --- Plot 1: The Energy Barrier (Loss Landscape) ---
    ax = axes[0]
    ax.plot(alphas, energies_init, 'k--', label='Linear Interpolation (Rigid)')
    if energies_relax is not None and len(energies_relax) > 0:
        ax.plot(alphas, energies_relax, 'b-o', label='Restrained Relaxation')
    ax.set_ylabel("Fitting Error (Objective)")
    ax.set_title("Energy Barrier Between Solutions")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # --- Plot 2: Parameter Evolution ---
    ax = axes[1]
    # Scale parameters that are too large or too small for [-1, 1] plot
    dofs_plot = np.copy(chain_relaxed)
    scales = np.ones(len(param_names))
    for j in range(len(param_names)):
        max_abs = np.max(np.abs(dofs_plot[:, j]))
        if max_abs > 1.0:
            scales[j] = 0.1
        elif max_abs > 0.0 and max_abs < 0.1:
            scales[j] = 10.0
        dofs_plot[:, j] *= scales[j]
        
    for j in range(len(param_names)):
        label = param_names[j]
        if scales[j] != 1.0:
            label += f" (x{scales[j]:.1f})"
        ax.plot(alphas, dofs_plot[:, j], '.-', label=label)
        
    ax.set_ylabel("Parameter Value")
    ax.set_title("Parameter Evolution along Relaxed Path")
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize='small', ncol=2) 
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1.1, 1.1)
    
    # --- Plot 3: Variational Forces (Gradients) ---
    ax = axes[2]
    for i, name in enumerate(param_names):
        ax.plot(alphas, gradients[:, i], '.-', label=f"dE/d({name})")
    ax.axhline(0, color='k', linestyle=':')
    ax.set_ylabel("Gradient (Force)")
    ax.set_xlabel("Interpolation Progress (α)")
    ax.set_title("Gradients along the path")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_string.png", dpi=150)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="1D String interpolation and Hessian analysis")
    parser.add_argument("--xyz", type=str, default="../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-z.xyz")
    parser.add_argument("--dofA", type=str, default="dofSelection_saved_Fit1.dat")
    parser.add_argument("--dofB", type=str, default="dofSelection_MorseSR_H2O.dat")
    parser.add_argument("--n_points", type=int, default=11)
    parser.add_argument("--kmorse", type=float, default=1.7)
    parser.add_argument("--lepairs", type=float, default=0.9)
    parser.add_argument("--springK", type=float, default=100.0)
    parser.add_argument("--nstep", type=int, default=50)
    parser.add_argument("--out_prefix", type=str, default="path")
    parser.add_argument("--relax", action="store_true", help="Run constrained relaxation")
    parser.add_argument("--hessian", action="store_true", help="Compute Hessian at alpha=0, 0.5, 1")
    args = parser.parse_args()

    # Load DOFs to get arrays A and B
    dofsA, hyperA = bu.parse_dof_file(args.dofA)
    dofsB, hyperB = bu.parse_dof_file(args.dofB)
    valsA = bu.get_dof_values(dofsA)
    valsB = bu.get_dof_values(dofsB)
    
    def hval(h, key, default):
        return h.get(key, default)
        
    kA, kB = hval(hyperA, 'kMorse', args.kmorse), hval(hyperB, 'kMorse', args.kmorse)
    lA, lB = hval(hyperA, 'Lepairs', args.lepairs), hval(hyperB, 'Lepairs', args.lepairs)
    wA, wB = hval(hyperA, 'weight_alpha', 1.0), hval(hyperB, 'weight_alpha', 1.0)
    hSA, hSB = hval(hyperA, 'hScale', 0.0), hval(hyperB, 'hScale', 0.0)
    
    # Setup initially with A to get names
    dof_names = bu.setup_fitreq(args.dofA, args.xyz, kMorse=kA, Lepairs=lA, hScale=hSA)
    Erefs, x0s = fit.read_xyz_data(args.xyz)
    
    alphas = np.linspace(0, 1, args.n_points)
    chain_initial = [(1 - a) * valsA + a * valsB for a in alphas]
    
    chain_relaxed = []
    energies_init = []
    energies_relax = []
    gradients = []
    
    kM_path, Lp_path, wa_path, hS_path = [], [], [], []
    
    # Copy structure for temporary use
    import copy
    dof_tmp = copy.deepcopy(dofsA)

    print(f"Running String Scan ({args.n_points} points). Relax: {args.relax}")
    for i, p_target in enumerate(chain_initial):
        a = alphas[i]
        print(f"Step {i+1}/{args.n_points} (alpha={a:.2f})")
        
        # Interpolate hyperparams
        kM = (1-a)*kA + a*kB
        Lp = (1-a)*lA + a*lB
        wa = (1-a)*wA + a*wB
        hS = (1-a)*hSA + a*hSB
        
        kM_path.append(kM); Lp_path.append(Lp); wa_path.append(wa); hS_path.append(hS)
        
        # Set interpolated physics
        fit.setGlobalParams(kMorse=kM, Lepairs=Lp, hScale=hS, softClamp_start=4.0, softClamp_max=6.0)
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=wa)
        weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
        fit.setWeights(weights0)
        
        # 1. Unrelaxed point
        bu.set_dof_values(dof_tmp, p_target)
        for j in range(len(p_target)): fit.DOFs[j] = p_target[j]
        E_init, _, grads_init = bu.evaluate_current(bFs=True)
        energies_init.append(E_init)
        
        # 2. Relaxed point (optional)
        if args.relax:
            # Set K0 for spring constraints
            for d in dof_tmp: d['K0'] = args.springK
            p_final = bu.relax_tethered(dof_tmp, kMorse=kM, Lepairs=Lp, nstep=args.nstep)
            # Eval relaxed to get final energy & gradients
            E_final, _, grads_final = bu.evaluate_current(bFs=True)
            
            chain_relaxed.append(p_final)
            energies_relax.append(E_final)
            gradients.append(grads_final)
        else:
            chain_relaxed.append(p_target)
            gradients.append(grads_init)
            
    chain_relaxed = np.array(chain_relaxed)
    energies_init = np.array(energies_init)
    energies_relax = np.array(energies_relax) if args.relax else None
    gradients = np.array(gradients)
    
    # Save trajectory data for DR
    np.savez(f"{args.out_prefix}_trajectory.npz", 
             alphas=alphas, 
             chain_initial=np.array(chain_initial),
             chain_relaxed=chain_relaxed, 
             energies_init=energies_init, 
             energies_relax=energies_relax if args.relax else np.array([]),
             param_names=dof_names,
             kMorse_path=np.array(kM_path),
             Lepairs_path=np.array(Lp_path),
             weight_alpha_path=np.array(wa_path),
             hScale_path=np.array(hS_path))
             
    plot_string_results(alphas, energies_init, energies_relax, chain_relaxed, gradients, dof_names, args.out_prefix)
    
    if args.hessian:
        print("\nComputing Hessians...")
        for h_alpha in [0.0, 0.5, 1.0]:
            idx = np.argmin(np.abs(alphas - h_alpha))
            a = alphas[idx]
            print(f"Hessian at alpha ~ {a:.2f}")
            
            # Re-apply interpolated physics for this point
            fit.setGlobalParams(kMorse=kM_path[idx], Lepairs=Lp_path[idx], hScale=hS_path[idx], softClamp_start=4.0, softClamp_max=6.0)
            weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=wa_path[idx])
            weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
            fit.setWeights(weights0)
            
            bu.set_dof_values(dof_tmp, chain_relaxed[idx])
            H, eigvals, eigvecs = bu.compute_fd_hessian(dof_tmp)
            
            # Print minimal summary
            print(f"  Eigenvalues: {eigvals}")
            
            # Plot eigenvalues
            plt.figure(figsize=(6,4))
            plt.plot(eigvals, 'o-')
            plt.axhline(0, color='k', linestyle='--')
            plt.title(f"Hessian Eigenvalues (alpha={alphas[idx]:.2f})")
            plt.ylabel("Eigenvalue")
            plt.xlabel("Mode Index")
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{args.out_prefix}_hessian_eigs_alpha_{alphas[idx]:.2f}.png")
            plt.close()

if __name__ == "__main__":
    main()
