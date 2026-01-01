import numpy as np
import os
import sys
import argparse
import matplotlib.pyplot as plt

# Adjust path to FireCore/pyBall
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL import OMM_ocl as omm

_fb_initialized = False

def init_fireball_water(verbosity=0):
    """Big-picture: Initialize Fireball with a water molecule for testing."""
    global _fb_initialized
    if _fb_initialized:
        return  # Avoid double-initialize -> Fortran double alloc crash

    atomTypes_Z = np.array([8, 1, 1], dtype=np.int32)
    atomPos = np.array([
        [0.00000000, 0.00000000, 0.00000000],
        [-0.75806320, 0.63581013, 0.00000000],
        [0.75806639, 0.63580735, 0.00000000]
    ], dtype=np.float64)

    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    fc.setVerbosity(verbosity)
    fc.assembleH(positions=atomPos, iforce=0, Kscf=1)
    _fb_initialized = True

def check_hamiltonian_diagonal(solver, nAtoms):
    """Optional: Compare sparse Hamiltonian diagonal terms with k-space reference."""
    h_mat, _, atomNeigh, _, nMaxNeigh = solver.get_fireball_data()
    
    kvec_gamma = [0.0, 0.0, 0.0]
    dims = fc.get_HS_dims()
    Hk, _ = fc.get_HS_k(kvec_gamma, dims.norbitals)
    s_indices = [0, 4, 5]
    Hk_s = Hk[np.ix_(s_indices, s_indices)].real
    
    print("\n--- Hamiltonian Diagonal Comparison ---")
    for i in range(nAtoms):
        ineg_self = -1
        for ineg in range(nMaxNeigh):
            if atomNeigh[i, ineg] == i:
                ineg_self = ineg
                break
        if ineg_self >= 0:
            print(f"Atom {i} (s-orb) self: Sparse={h_mat[i,ineg_self]:.6f}, Hk_diag={Hk_s[i,i]:.6f}")
        else:
            print(f"Atom {i} (s-orb) self: NOT FOUND in neighbor list")

def test_omm(args):
    """Main test function for OMM OpenCL implementation."""
    init_fireball_water(verbosity=args.verbosity)
    
    # Define some "Molecular Orbitals" (LMOs)
    # For H2O with 3 atoms, let's make 3 LMOs, each localized on its atom
    nOrbs = 3
    nMaxSupport = 4
    orbSupport = np.full((nOrbs, nMaxSupport), -1, dtype=np.int32)
    orbCoefs   = np.zeros((nOrbs, nMaxSupport), dtype=np.float32)
    nSupport   = np.ones(nOrbs, dtype=np.int32)
    
    for i in range(nOrbs):
        orbSupport[i, 0] = i
        orbCoefs[i, 0]   = 1.0
        
    # Initialize OMM solver and setup data from Fireball
    solver = omm.OMM_OCL()
    solver.setup_from_fireball(nOrbs, orbSupport, nSupport, K_stiff=args.K)
    # MUST upload initial coefficients, otherwise they are zero in GPU!
    solver.upload_data(orbCoefs=orbCoefs)
    
    # Optional debugging
    check_hamiltonian_diagonal(solver, nOrbs)
    
    # Use built-in optimizer (simple gradient descent)
    history, last_s, last_h, grads_final, orbCoefs = solver.optimize(
        orbCoefs, steps=args.steps, lr=args.lr, log=bool(args.verbosity)
    )

    print("\nFinal Orbital Gradients (last step):")
    print(grads_final)

    if args.plot:
        # Plot histories
        fig, axes = plt.subplots(2, 1, figsize=(8, 8))
        
        axes[0].plot(history['hist_maxE'], label='maxE (diag H)')
        axes[0].plot(history['hist_totE'], label='totalE (sum diag H)')
        axes[0].set_ylabel('Energy (eV)')
        axes[0].set_title('OMM Energy History')
        axes[0].legend()
        axes[0].grid(True)

        axes[1].plot(history['hist_offdiag'], label='max|S_off|')
        axes[1].plot(history['hist_diagdev'], label='max|Sii-1|')
        if 'hist_gnorm' in history:
            axes[1].plot(history['hist_gnorm'], label='|g|')
        axes[1].set_yscale('log')
        axes[1].set_ylabel('Error / Gradient Norm')
        axes[1].set_xlabel('Iteration')
        axes[1].set_title('OMM Convergence History')
        axes[1].legend()
        axes[1].grid(True)

        plt.tight_layout()
        plt.savefig(args.metrics_png)

        # Plot final overlap matrix heatmap
        plt.figure(figsize=(4,3))
        plt.imshow(last_s, cmap='viridis')
        plt.colorbar(label='S_ij')
        plt.title('Final Overlap Matrix S_ij')
        plt.tight_layout()
        plt.savefig(args.overlap_png)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--steps",       type=int,   default=10,    help="Number of OMM iterations.")
    p.add_argument("--lr",          type=float, default=1e-3,  help="Gradient descent step size.")
    p.add_argument("--K",           type=float, default=100.0, help="Penalty stiffness for overlap (Lambda).")
    p.add_argument("--verbosity",   type=int,   default=1,     help="Fireball verbosity (0=quiet).")
    p.add_argument("--plot",        type=int,   default=1,     help="Enable plotting (0/1).")
    p.add_argument("--metrics-png", type=str,   default="omm_metrics.png", help="Output path for metrics plot.")
    p.add_argument("--overlap-png", type=str,   default="omm_overlap.png", help="Output path for overlap heatmap.")
    p.add_argument("--show",        type=int,   default=1,     help="Show plots interactively (0/1).")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()
    test_omm(args)
    if args.show:
        plt.show()
