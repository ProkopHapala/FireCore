"""
KekuleCLI.py - Command line interface for Kekule Topology Solver
================================================================
Usage: python KekuleCLI.py --delta 0.5 --pbc_x --pbc_y --save
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
from KekuleSolver import KekuleSolver

def main():
    parser = argparse.ArgumentParser(description="Kekule Solver CLI")
    parser.add_argument("--delta", type=float, default=0.4)
    parser.add_argument("--phi2", type=float, default=2.094)
    parser.add_argument("--Lx", type=float, default=20.0)
    parser.add_argument("--Ly", type=float, default=10.0)
    parser.add_argument("--pbc_x", action="store_true")
    parser.add_argument("--pbc_y", action="store_true")
    parser.add_argument("--onsite", type=float, default=0.0)
    parser.add_argument("--window", type=float, nargs=2, default=[-0.1, 0.1])
    parser.add_argument("--save", action="store_true")

    args = parser.parse_args()

    solver = KekuleSolver(delta=args.delta)
    solver.phi2 = args.phi2
    solver.Lx = args.Lx
    solver.Ly = args.Ly
    solver.pbc = (args.pbc_x, args.pbc_y)
    
    print("Generating lattice...")
    solver.generate_lattice()
    
    print(f"Applying onsite energy {args.onsite} at boundary...")
    solver.set_onsite_at_boundary(energy=args.onsite)
    
    print("Building Hamiltonian and solving...")
    solver.build_hamiltonian()
    solver.solve()
    
    mask = (solver.evals >= args.window[0]) & (solver.evals <= args.window[1])
    print(f"Found {np.sum(mask)} states in window {args.window}")
    if np.sum(mask) > 0:
        print(f"Eigenvalues in window: {solver.evals[mask]}")
    else:
        # Print the 4 states closest to zero
        closest_idx = np.argsort(np.abs(solver.evals))[:4]
        print(f"Closest eigenvalues to zero: {solver.evals[closest_idx]}")
    
    if args.save:
        plt.figure(figsize=(10, 6))
        if np.sum(mask) > 0:
            ldos = solver.get_ldos(args.window[0], args.window[1])
        else:
            # Average LDOS of 4 closest states
            closest_idx = np.argsort(np.abs(solver.evals))[:4]
            ldos = np.sum(np.abs(solver.evecs[:, closest_idx])**2, axis=1)
            print(f"Plotting average LDOS of {len(closest_idx)} closest states.")
        
        plt.scatter(solver.pos[:,0], solver.pos[:,1], c=ldos, cmap='magma', s=50)
        plt.colorbar(label='LDOS')
        plt.title(f"Kekule Topology - PBC: {solver.pbc}")
        plt.axis('equal')
        plt.savefig("Kekule_CLI_Output.png")
        print("Plot saved to Kekule_CLI_Output.png")

if __name__ == "__main__":
    main()
