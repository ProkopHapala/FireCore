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
    for i, name in enumerate(param_names):
        ax.plot(alphas, chain_relaxed[:, i], '.-', label=name)
    ax.set_ylabel("Parameter Value")
    ax.set_title("How Parameters Shift Over the Ridge")
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize='small', ncol=2) 
    ax.grid(True, alpha=0.3)
    
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

def plot_hessian_analysis(H, eigvals, eigvecs, param_names, alpha, out_prefix):
    # 1. Plot soft mode composition (Eigenvector for lowest eigenvalue)
    softest_idx = np.argmin(eigvals)
    softest_vec = eigvecs[:, softest_idx]
    
    plt.figure(figsize=(8, 4))
    x = np.arange(len(param_names))
    plt.bar(x, softest_vec, color='coral')
    plt.axhline(0, color='k', linewidth=0.5)
    plt.xticks(x, param_names, rotation=45, ha='right')
    plt.title(f"Softest Mode Composition (alpha={alpha:.2f}, eigval={eigvals[softest_idx]:.2e})")
    plt.ylabel("Component weight")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_hessian_softmode_alpha_{alpha:.2f}.png", dpi=150)
    plt.close()
    
    # 2. Correlation Heatmap (Inverse Hessian / Pseudo-Covariance)
    # Filter out near-zero eigenvalues to avoid huge numbers
    valid_idx = np.abs(eigvals) > 1e-6
    if np.sum(valid_idx) > 0:
        H_inv = (eigvecs[:, valid_idx] / eigvals[valid_idx]) @ eigvecs[:, valid_idx].T
        # Normalize to [-1, 1] correlation matrix style
        d = np.diag(H_inv).copy()
        # Avoid div zero
        d[np.abs(d) < 1e-10] = 1e-10
        corr = H_inv / np.sqrt(np.outer(d, d))
        
        plt.figure(figsize=(8, 6))
        plt.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1)
        plt.colorbar(label="Correlation (Red=Positive, Blue=Negative)")
        plt.xticks(x, param_names, rotation=45, ha='right')
        plt.yticks(x, param_names)
        plt.title(f"Parameter Correlation (alpha={alpha:.2f})")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_hessian_corr_alpha_{alpha:.2f}.png", dpi=150)
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Hyperparameter-Aware String interpolation and Hessian analysis")
    parser.add_argument("--data", type=str, default="grid_search_out/ensemble_data.npz")
    parser.add_argument("--idxA", type=int, default=0, help="Index of Endpoint A in ensemble")
    parser.add_argument("--idxB", type=int, default=47, help="Index of Endpoint B in ensemble")
    parser.add_argument("--n_points", type=int, default=11)
    parser.add_argument("--springK", type=float, default=100.0)
    parser.add_argument("--nstep", type=int, default=50)
    parser.add_argument("--out_prefix", type=str, default="grid_search_out/path")
    args = parser.parse_args()

    xyz_path = "../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-z.xyz"
    base_dof_file = "dofSelection_MorseSR_H2O.dat"
    
    data = np.load(args.data)
    dofs_all = data['dofs']
    kM_all = data['kMorse']
    Lep_all = data['Lepairs']
    wa_all = data['weight_alpha']
    dof_names = data['dof_names']
    
    valA = dofs_all[args.idxA]
    valB = dofs_all[args.idxB]
    kMA, kMB = kM_all[args.idxA], kM_all[args.idxB]
    LepA, LepB = Lep_all[args.idxA], Lep_all[args.idxB]
    waA, waB = wa_all[args.idxA], wa_all[args.idxB]

    # Setup 
    base_dofs = bu.parse_dof_file(base_dof_file)
    # Filter DOFs: remove epair charge Q (comp == 2)
    dofs_struct = []
    for d in base_dofs:
        if int(d['comp']) == 2: continue
        dofs_struct.append(d)
        
    bu.write_dof_file("tmp_string_dofs.dat", dofs_struct)
    bu.setup_fitreq("tmp_string_dofs.dat", xyz_path, imodel=7, Regularize=1)
    Erefs, x0s = fit.read_xyz_data(xyz_path)
    
    alphas = np.linspace(0, 1, args.n_points)
    
    chain_initial = [(1 - a) * valA + a * valB for a in alphas]
    interp_kM = [(1 - a) * kMA + a * kMB for a in alphas]
    interp_Lep = [(1 - a) * LepA + a * LepB for a in alphas]
    interp_wa = [(1 - a) * waA + a * waB for a in alphas]
    
    chain_relaxed = []
    energies_init = []
    energies_relax = []
    gradients = []
    
    print(f"Running Hyperparameter-Aware String Scan ({args.n_points} points).")
    
    for i, a in enumerate(alphas):
        print(f"Step {i+1}/{args.n_points} (alpha={a:.2f}): kM={interp_kM[i]:.2f}, Lep={interp_Lep[i]:.2f}, wa={interp_wa[i]:.2f}")
        
        # Set hyperparameters
        fit.setGlobalParams(kMorse=interp_kM[i], Lepairs=interp_Lep[i], softClamp_start=4.0, softClamp_max=6.0)
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=interp_wa[i])
        weights0, _ = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
        fit.setWeights(weights0)
        
        p_target = chain_initial[i]
        
        # 1. Unrelaxed point
        bu.set_dof_values(dofs_struct, p_target)
        for j in range(len(p_target)): fit.DOFs[j] = p_target[j]
        E_init, _, grads_init = bu.evaluate_current(bFs=True)
        energies_init.append(E_init)
        
        # 2. Relaxed point (always on for this detailed scan)
        for d in dofs_struct: d['K0'] = args.springK
        p_final = bu.relax_tethered(dofs_struct, kMorse=interp_kM[i], Lepairs=interp_Lep[i], nstep=args.nstep)
        E_final, _, grads_final = bu.evaluate_current(bFs=True)
        
        chain_relaxed.append(p_final)
        energies_relax.append(E_final)
        gradients.append(grads_final)
            
    chain_relaxed = np.array(chain_relaxed)
    energies_init = np.array(energies_init)
    energies_relax = np.array(energies_relax)
    gradients = np.array(gradients)
             
    plot_string_results(alphas, energies_init, energies_relax, chain_relaxed, gradients, dof_names, args.out_prefix)
    
    print("\nComputing Hessians...")
    for h_alpha in [0.0, 0.5, 1.0]:
        idx = np.argmin(np.abs(alphas - h_alpha))
        print(f"Hessian at alpha ~ {alphas[idx]:.2f}")
        
        # Set hyperparams again
        fit.setGlobalParams(kMorse=interp_kM[idx], Lepairs=interp_Lep[idx], softClamp_start=4.0, softClamp_max=6.0)
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=interp_wa[idx])
        weights0, _ = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
        fit.setWeights(weights0)
        
        bu.set_dof_values(dofs_struct, chain_relaxed[idx])
        H, eigvals, eigvecs = bu.compute_fd_hessian(dofs_struct)
        
        plot_hessian_analysis(H, eigvals, eigvecs, dof_names, alphas[idx], args.out_prefix)

if __name__ == "__main__":
    main()
