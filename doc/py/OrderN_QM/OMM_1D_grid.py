
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

# unlimited line lenght in numpy when printing
np.set_printoptions(linewidth=np.inf)

class OMM1DSolver:
    def __init__(self, n_grid=200, l_max=20.0):
        self.n = n_grid
        self.x = np.linspace(0, l_max, n_grid)
        self.dx = self.x[1] - self.x[0]
        
        # 1. Setup Hamiltonian: Kinetic (Finite Difference) + Potential
        # T = -0.5 * d^2/dx^2
        self.main_diag = np.ones(self.n) / (self.dx**2)
        self.off_diag = -0.5 * np.ones(self.n - 1) / (self.dx**2)
        self.H = np.diag(self.main_diag) + np.diag(self.off_diag, k=1) + np.diag(self.off_diag, k=-1)
        
        # V = harmonic well in the middle
        self.V = 0.5 * 0.5 * (self.x - l_max/2)**2
        np.fill_diagonal(self.H, np.diag(self.H) + self.V)

    def setup_orbitals(self, n_orb=4, support_size=40):
        """Initialize localized orbitals with fixed masks."""
        self.n_orb = n_orb
        self.support_size = support_size
        self.orbitals = []
        self.masks = []
        
        # Space centers evenly
        centers = np.linspace(self.n*0.3, self.n*0.7, n_orb).astype(int)
        
        for c in centers:
            mask = np.arange(max(0, c - support_size//2), min(self.n, c + support_size//2))
            self.masks.append(mask)
            
            # Initial guess: a normalized blob
            phi = np.zeros(self.n)
            phi[mask] = np.exp(-0.5 * ((self.x[mask] - self.x[c])/(self.dx*5))**2)
            phi /= np.linalg.norm(phi) 
            self.orbitals.append(phi)

    def recenter_supports(self, support_size=None):
        """Shift each support to center-of-mass of |phi|^2 and trim outside region."""
        if support_size is None:
            support_size = self.support_size
        half = support_size // 2
        for i in range(self.n_orb):
            phi = self.orbitals[i]
            rho = phi * phi
            norm_rho = rho.sum()
            if norm_rho < 1e-14:
                print(f"#DEBUG recenter skipped orb {i} (zero norm)")
                continue
            cog = (self.x * rho).sum() / norm_rho
            idx = int(np.clip(np.searchsorted(self.x, cog), 0, self.n - 1))
            new_start = max(0, idx - half)
            new_end = min(self.n, idx + half)
            new_mask = np.arange(new_start, new_end)
            old_center = self.masks[i][len(self.masks[i])//2] if len(self.masks[i]) else -1
            if old_center != idx:
                print(f"#DEBUG recenter orb {i}: center {old_center} -> {idx} (cog {cog:.3f})")
            new_phi = np.zeros_like(phi)
            new_phi[new_mask] = phi[new_mask]
            nrm = np.linalg.norm(new_phi)
            if nrm > 1e-12:
                new_phi /= nrm
            self.orbitals[i] = new_phi
            self.masks[i] = new_mask

    def orthogonalize(self, n_iter=10, damping=0.5):
        """
        The 'Projective' Constraint Solver.
        Forces <phi_i | phi_j> = 0 and <phi_i | phi_i> = 1
        """
        for _ in range(n_iter):
            # We use Gauss-Seidel style (immediate update) for stability
            for i in range(self.n_orb):
                mi = self.masks[i]
                
                # 1. Orthogonalization Step (Push away from neighbors)
                # Correction: d_phi_i = - sum_{j!=i} <phi_i | phi_j> * phi_j
                correction = np.zeros(self.n)
                for j in range(self.n_orb):
                    if i == j: continue
                    
                    overlap = np.dot(self.orbitals[i], self.orbitals[j])
                    # Only project within our own support
                    correction[mi] += overlap * self.orbitals[j][mi]
                
                self.orbitals[i][mi] -= damping * correction[mi]
                
                # 2. Normalization Step (Project onto unit sphere)
                norm = np.linalg.norm(self.orbitals[i])
                if norm > 1e-9:
                    self.orbitals[i] /= norm
                    
    def energy_step(self, step_size=0.01):
        """Gradient descent on the energy functional: E = sum <phi|H|phi>"""
        for i in range(self.n_orb):
            mi = self.masks[i]
            # Gradient of <phi|H|phi> is 2 * H * phi
            grad = 2.0 * (self.H @ self.orbitals[i])
            
            # Update only within support
            self.orbitals[i][mi] -= step_size * grad[mi]

    def get_stats(self):
        """Calculate total energy and max overlap error."""
        total_e = 0
        max_overlap = 0
        for i in range(self.n_orb):
            total_e += np.dot(self.orbitals[i], self.H @ self.orbitals[i])
            for j in range(i + 1, self.n_orb):
                max_overlap = max(max_overlap, abs(np.dot(self.orbitals[i], self.orbitals[j])))
        return total_e, max_overlap

    def density_from_orbitals(self):
        """Return electron density rho(x) from current (normalized) orbitals."""
        psi = np.array(self.orbitals)
        rho = np.sum(psi * psi, axis=0)
        return rho

    def reference_ground_state(self, n_occ=None):
        """
        Solve the banded eigenproblem with SciPy eigh_tridiagonal.
        Returns eigenvalues, eigenvectors, and ground-state density.
        """
        n_occ = n_occ if n_occ is not None else self.n_orb
        main = self.main_diag + self.V
        evals, evecs = eigh_tridiagonal(main, self.off_diag, select="i", select_range=(0, n_occ - 1))
        rho = np.sum(evecs[:, :n_occ] ** 2, axis=1)
        return evals, evecs, rho

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="1D orbital minimization vs banded reference")
    ap.add_argument("--n_grid",       type=int,   default=100,   help="number of grid points")
    ap.add_argument("--l_max",        type=float, default=30.0,  help="box length")
    ap.add_argument("--n_orb",        type=int,   default=3,     help="number of occupied orbitals")
    ap.add_argument("--support_size", type=int,   default=20,    help="mask size for each orbital")
    ap.add_argument("--n_iter",       type=int,   default=1000,   help="OMM outer iterations")
    ap.add_argument("--step_size",    type=float, default=0.01,  help="energy gradient step")
    ap.add_argument("--ortho_iter",   type=int,   default=5,     help="orthogonalization iterations per step")
    ap.add_argument("--ortho_damp",   type=float, default=0.4,   help="orthogonalization damping")
    ap.add_argument("--report_every", type=int,   default=20,    help="print energy/overlap every k steps")
    ap.add_argument("--compare_ref",  type=int,   default=1,     help="solve reference banded eigenproblem and compare densities")
    ap.add_argument("--recenter_every", type=int, default=1,     help="recenter supports every k steps (0 to disable)")
    ap.add_argument("--no_plot",      type=int,   default=0,     help="skip matplotlib plots")
    args = ap.parse_args()

    solver = OMM1DSolver(n_grid=args.n_grid, l_max=args.l_max)
    solver.setup_orbitals(n_orb=args.n_orb, support_size=args.support_size)

    print(f"{'Iter':>5} | {'Energy':>10} | {'Max Overlap':>12}")
    print("-" * 35)

    for i in range(args.n_iter):
        solver.energy_step(step_size=args.step_size)
        solver.orthogonalize(n_iter=args.ortho_iter, damping=args.ortho_damp)
        if args.recenter_every > 0 and (i % args.recenter_every == 0):
            solver.recenter_supports()
        if args.report_every > 0 and (i % args.report_every == 0):
            e, err = solver.get_stats()
            print(f"{i:5d} | {e:10.5f} | {err:12.6e}")

    rho_omm = solver.density_from_orbitals()

    if args.compare_ref:
        evals, evecs, rho_ref = solver.reference_ground_state(n_occ=args.n_orb)
        diff = rho_omm - rho_ref
        l2 = np.linalg.norm(diff)
        max_abs = np.max(np.abs(diff))
        total_charge_diff = np.sum(diff) * solver.dx
        print(f"#INFO ref_eigs_min={evals.min():.6f} ref_eigs_max={evals.max():.6f}")
        print(f"#INFO density L2={l2:.6e} max|diff|={max_abs:.6e} integral_diff={total_charge_diff:.6e}")

    if not args.no_plot:
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Orbital shapes
        ax[0].plot(solver.x, solver.V * 0.1, 'k--', alpha=0.3, label='Potential (scaled)')
        for i in range(solver.n_orb):
            ax[0].fill_between(solver.x, solver.orbitals[i], alpha=0.5, label=f'OMM Orbital {i}')
        if args.compare_ref:
            for i in range(args.n_orb):
                ax[0].plot(solver.x, evecs[:, i], lw=1.2, label=f'Ref Orbital {i}')
        ax[0].set_ylabel("ψ")
        ax[0].legend(ncol=2, fontsize=8)
        ax[0].set_title("Occupied orbitals (OMM filled; reference lines)")

        # Density comparison
        ax[1].plot(solver.x, rho_omm, 'k-', lw=2, label='Density OMM')
        if args.compare_ref:
            ax[1].plot(solver.x, rho_ref, 'r--', lw=1.5, label='Density reference')
        ax[1].set_xlabel("x")
        ax[1].set_ylabel("ρ")
        ax[1].legend()
        ax[1].set_title("Electron density comparison")

        # DEBUG: visualize current supports as light spans
        colors = plt.cm.tab10.colors
        for i, mi in enumerate(solver.masks):
            if len(mi) == 0:
                continue
            x0 = solver.x[mi[0]]
            x1 = solver.x[mi[-1]]
            ax[0].axvspan(x0, x1, color=colors[i % len(colors)], alpha=0.1, lw=0)
            ax[0].axvline(x0, color=colors[i % len(colors)], ls=":", lw=0.8)
            ax[0].axvline(x1, color=colors[i % len(colors)], ls=":", lw=0.8)

        plt.tight_layout()
        plt.show()

