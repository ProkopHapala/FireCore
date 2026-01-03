import argparse
import numpy as np
import pyopencl as cl
import scipy.linalg
import matplotlib.pyplot as plt

import sys
from pathlib import Path
# add repo root (…/FireCore) to sys.path so pyBall is importable
repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.append(str(repo_root))
from pyBall.FireballOCL.CheFSI import FrontierSolver

def generate_1d_system(n_atoms=64, distort_amp=0.1, seed=None):
    if seed is not None: np.random.seed(seed)
    # Positions
    positions = np.arange(n_atoms, dtype=np.float32) * 1.0
    positions += np.random.uniform(-distort_amp, distort_amp, n_atoms).astype(np.float32)
    # Matrices
    H_dense = np.zeros((n_atoms, n_atoms), dtype=np.float32)
    S_dense = np.eye(n_atoms, dtype=np.float32)
    n_max_neighs = 3
    indices = np.full((n_atoms, n_max_neighs), -1, dtype=np.int32)
    H_val = np.zeros((n_atoms, n_max_neighs), dtype=np.float32)
    S_val = np.zeros((n_atoms, n_max_neighs), dtype=np.float32)
    
    for i in range(n_atoms):
        H_dense[i, i] = np.random.uniform(-0.05, 0.05)
        indices[i, 1] = i; H_val[i, 1] = H_dense[i, i]; S_val[i, 1] = 1.0
        
        for idx, j in enumerate([i-1, i+1]):
            if 0 <= j < n_atoms:
                dist = abs(positions[i] - positions[j])
                fac = np.exp(-1.0 * (dist - 1.0)**2)
                h_elem = -1.0 * fac
                s_elem = 0.2 * fac
                
                H_dense[i, j] = h_elem; S_dense[i, j] = s_elem
                col = 0 if idx==0 else 2
                indices[i, col] = j; H_val[i, col] = h_elem; S_val[i, col] = s_elem
    return H_dense, S_dense, indices, H_val, S_val, positions

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CheFSI (variance) debug test with stability diagnostics.")
    parser.add_argument("--n-atoms",      type=int,   default=64)
    parser.add_argument("--n-vecs",       type=int,   default=8)
    parser.add_argument("--n-iter",       type=int,   default=20,  help="Outer iterations (more, but filter_order small).")
    parser.add_argument("--filter-order", type=int,   default=5,   help="Gradient steps between orthogonalizations.")
    parser.add_argument("--target",       type=float, default=0.5)
    parser.add_argument("--tol",          type=float, default=1e-2, help="Convergence on max_non_orth.")
    parser.add_argument("--distort-amp",  type=float, default=0.1,  help="Lattice distortion amplitude.")
    parser.add_argument("--seed",         type=int,   default=None, help="Random seed.")
    parser.add_argument("--plot-orbitals",type=str,   default="0,1,2,3", help="Comma-separated orbital indices to plot.")
    parser.add_argument("--error-scale",  type=float, default=1.0, help="Multiplier for plotted (GPU-Ref) error.")
    args = parser.parse_args()

    # System
    H_dns, S_dns, ell_idx, ell_H, ell_S, pos = generate_1d_system(args.n_atoms, distort_amp=args.distort_amp, seed=args.seed)
    
    # Reference
    ref_evals, ref_evecs = scipy.linalg.eigh(H_dns, S_dns)
    closest_indices = np.argsort(np.abs(ref_evals - args.target))
    print(f"Target Energy: {args.target}")
    print(f"Exact spectrum range: [{ref_evals[0]:.2f}, {ref_evals[-1]:.2f}]")
    print(f"Closest Exact Eigenvalues: {ref_evals[closest_indices[:4]]}")
    
    # Solver (let FrontierSolver pick device/queue internally)
    solver = FrontierSolver(n_vecs=args.n_vecs, tile_size=32, verbose=True)
    solver.prepare(ell_idx, ell_H, ell_S)
    
    print("\n--- Solver Parameters ---")
    print(f"n_atoms={args.n_atoms}, n_vecs={args.n_vecs}, n_iter={args.n_iter}, filter_order={args.filter_order}, tol={args.tol}")
    print(f"target={args.target}, bounds=auto, distort_amp={args.distort_amp}, seed={args.seed}, error_scale={args.error_scale}")
    print("Overlap = |<v_gpu | S | v_ref>| (generalized inner product)")
    
    gpu_evals, gpu_vecs_dev, logs = solver.solve(
        target_energy=args.target, 
        bounds=None, 
        n_iter=args.n_iter, 
        filter_order=args.filter_order,
        tol=args.tol
    )
    gpu_vecs = gpu_vecs_dev.get()
    
    # Analysis
    print("\n--- Comparison Results ---")
    print("Columns: GPU state | Ref idx | GPU eigenvalue | Exact eigenvalue | Abs diff | Residual ||(H - eS)v|| | Overlap | RMSE vs exact eigenvector")
    print(f"{'GPU':<6} | {'Ref':<6} | {'GPU E':<10} | {'Exact E':<10} | {'Diff':<10} | {'Residual':<10} | {'Overlap':<10} | {'RMSE':<10}")
    print("-" * 100)
    
    rmse_list = []
    order = np.argsort(gpu_evals)  # sort by GPU eigenvalue
    for i in order:
        e_gpu = gpu_evals[i]
        idx_ref = np.argmin(np.abs(ref_evals - e_gpu))
        e_ref = ref_evals[idx_ref]
        diff = abs(e_gpu - e_ref)
        
        v = gpu_vecs[:, i]
        Hv = H_dns @ v
        Sv = S_dns @ v
        res_vec = Hv - e_gpu * Sv
        residual = np.linalg.norm(res_vec)
        
        v_ref = ref_evecs[:, idx_ref]
        overlap = abs(np.dot(v.T, np.dot(S_dns, v_ref)))
        rmse = np.sqrt(np.mean((v - v_ref)**2))
        rmse_list.append(rmse)
        
        print(f"{i:<6} | {idx_ref:<6} | {e_gpu:<10.5f} | {e_ref:<10.5f} | {diff:<10.1e} | {residual:<10.1e} | {overlap:<10.4f} | {rmse:<10.3e}")

    print(f"\nMax non-orth logs (len={len(logs['max_non_orth'])}): last={logs['max_non_orth'][-1]:.3e}")
    print(f"Total SpMM calls (matvecs): {logs.get('matvec_count', 'n/a')}")
    if rmse_list:
        print(f"Mean orbital RMSE: {np.mean(rmse_list):.3e}, max RMSE: {np.max(rmse_list):.3e}")
    
    # Plot Convergence
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(logs['max_non_orth'], 'r-o')
    plt.title("Non-Orthogonality (before correction)")
    plt.yscale('log')
    plt.xlabel("Iteration")
    
    plt.subplot(1, 3, 2)
    plt.plot(logs['cond'], 'b-o')
    plt.title("Gram Matrix Condition Number")
    plt.yscale('log')
    plt.axhline(1e12, color='k', linestyle='--', label="Unstable")
    plt.xlabel("Iteration")
    plt.legend()
    
    plt.subplot(1, 3, 3)
    plt.plot(logs['max_non_orth'], label="Non-orth")
    plt.plot(logs['cond'], label="Cond")
    plt.yscale('log')
    plt.title("Combined Diagnostics")
    plt.legend()
    plt.tight_layout()
    
    # Plot selected orbitals along the chain
    plot_indices = [int(x) for x in args.plot_orbitals.split(",") if x.strip() != ""]
    n_plots = len(plot_indices)
    fig, axes = plt.subplots(n_plots, 1, figsize=(8, 2.5*n_plots), sharex=True)
    if n_plots == 1:
        axes = [axes]
    for ax, idx in zip(axes, plot_indices):
        if idx >= args.n_vecs:
            ax.set_title(f"Orbital {idx} (out of range)")
            continue
        e_gpu = gpu_evals[idx]
        ref_idx = np.argmin(np.abs(ref_evals - e_gpu))
        ref_vec = ref_evecs[:, ref_idx]
        gpu_vec = gpu_vecs[:, idx]
        phase = np.sign(np.dot(ref_vec, gpu_vec))
        gpu_vec_aligned = gpu_vec * phase
        error = (gpu_vec_aligned - ref_vec) * args.error_scale
        ax.plot(pos, ref_vec, label=f"Ref {ref_idx}", linewidth=2)
        ax.plot(pos, gpu_vec_aligned, "--", label=f"GPU {idx}", linewidth=2)
        ax.plot(pos, error, ":", label=f"(GPU-Ref)*{args.error_scale:g}", linewidth=1.5)
        ax.set_ylabel("Coeff / Error*scale")
        ax.legend()
        ax.grid(True, linestyle="--", alpha=0.3)
    axes[-1].set_xlabel("Position (1D chain)")
    fig.tight_layout()
    plt.show()