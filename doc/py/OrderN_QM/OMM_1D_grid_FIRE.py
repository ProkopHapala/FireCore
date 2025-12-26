import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy.linalg import eigh_tridiagonal

class FIRE:
    """
    Fast Iterative Relaxation Engine.
    Used here to relax the orbital shapes toward the ground state.
    It adaptively adjusts the time-step (dt) and mixing (alpha) based on 
    the alignment of the force and velocity (Power).
    """
    def __init__(self, size, dt=0.05, dt_max=0.2, N_min=5, f_inc=1.1, f_dec=0.5, alpha_start=0.1, f_alpha=0.99):
        self.v = [np.zeros(s) for s in size]
        self.dt = dt
        self.dt_max = dt_max
        self.N_min = N_min
        self.f_inc = f_inc
        self.f_dec = f_dec
        self.alpha_start = alpha_start
        self.alpha = alpha_start
        self.f_alpha = f_alpha
        self.n_pos = 0

    def step(self, x_list, f_list):
        for i in range(len(x_list)):
            P = np.dot(f_list[i], self.v[i])
            v_norm = np.linalg.norm(self.v[i])
            f_norm = np.linalg.norm(f_list[i])
            
            if f_norm > 1e-12:
                self.v[i] = (1.0 - self.alpha) * self.v[i] + self.alpha * (f_list[i] / f_norm) * v_norm
            
            if P > 0:
                if self.n_pos > self.N_min:
                    self.dt = min(self.dt * self.f_inc, self.dt_max)
                    self.alpha *= self.f_alpha
                self.n_pos += 1
            else:
                self.dt *= self.f_dec
                self.alpha = self.alpha_start
                self.v[i][:] = 0.0
                self.n_pos = 0
            
            self.v[i] += f_list[i] * self.dt
            x_list[i] += self.v[i] * self.dt

class OMM_Solver:
    def __init__(self, n_grid=150, l_max=10.0, n_orb=4, support_size=30):
        self.n = n_grid
        self.x = np.linspace(0, l_max, n_grid)
        self.dx = self.x[1] - self.x[0]
        self.n_orb = n_orb
        self.support_size = support_size
        
        # Potential: Harmonic Well
        self.V = 2.0 * (self.x - l_max/2)**2 
        
        self.orbitals = []
        self.masks = []      # [start_idx, end_idx]
        self.prev_delta = [] # Momentum buffer for Orthogonalizer
        self.add_orbitals(n_orb=n_orb, support_size=support_size)

    def add_orbitals(self, n_orb=4, support_size=30):
        self.support_size = support_size
        centers = np.linspace(self.n*0.2, self.n*0.8, n_orb).astype(int)
        self.orbitals = []
        self.masks = []
        self.prev_delta = []
        for c in centers:
            s = c - support_size//2
            self.masks.append([s, s + support_size])
            # Start with a simple blob
            phi = np.exp(-0.5 * (np.arange(support_size)-support_size/2)**2 / 4.0)
            phi /= np.linalg.norm(phi)
            self.orbitals.append(phi)
            self.prev_delta.append(np.zeros(support_size))

    def get_forces(self):
        """
        Calculates the Residual Force: F = -(H*phi - epsilon*phi)
        This ensures the force is always orthogonal to the orbital.
        """
        forces = []
        total_E = 0
        for i in range(len(self.orbitals)):
            s, e = self.masks[i]
            phi_loc = self.orbitals[i]
            
            # 1. 3-point Laplacian with 1-pixel zero padding
            # This forces the boundary values and derivatives to zero.
            padded = np.zeros(self.support_size + 2)
            padded[1:-1] = phi_loc
            lap = (padded[:-2] - 2*padded[1:-1] + padded[2:]) / (self.dx**2)
            
            h_phi_kin = -0.5 * lap
            
            # 2. Potential Energy
            grid_idx = np.arange(s, e) % self.n
            h_phi_pot = self.V[grid_idx] * phi_loc
            
            h_phi = h_phi_kin + h_phi_pot
            
            # 3. Local Energy (Rayleigh Quotient)
            epsilon = np.dot(phi_loc, h_phi)
            total_E += epsilon
            
            # 4. Residual Force (Negative Gradient on the unit sphere)
            force = -(h_phi - epsilon * phi_loc)
            forces.append(force)
            
        return forces, total_E

    def energy_terms(self):
        """Return total, kinetic, potential energies."""
        E_kin = 0.0
        E_pot = 0.0
        for i in range(len(self.orbitals)):
            s, e = self.masks[i]
            phi = self.orbitals[i]
            padded = np.zeros(self.support_size + 2)
            padded[1:-1] = phi
            lap = (padded[:-2] - 2*padded[1:-1] + padded[2:]) / (self.dx**2)
            E_kin += -0.5 * np.dot(phi, lap) * self.dx
            grid_idx = np.arange(s, e) % self.n
            E_pot += np.dot(phi * phi, self.V[grid_idx]) * self.dx
        return E_kin + E_pot, E_kin, E_pot

    def overlap_error(self):
        """Max |<phi_i|phi_j>|."""
        max_ov = 0.0
        for i in range(len(self.orbitals)):
            si, ei = self.masks[i]
            mi = np.arange(si, ei) % self.n
            for j in range(i+1, len(self.orbitals)):
                sj, ej = self.masks[j]
                mj = np.arange(sj, ej) % self.n
                shared, idx_i, idx_j = np.intersect1d(mi, mj, return_indices=True)
                if len(shared) == 0:
                    continue
                ov = abs(np.dot(self.orbitals[i][idx_i], self.orbitals[j][idx_j]))
                if ov > max_ov:
                    max_ov = ov
        return max_ov

    def density(self):
        """Return density on full grid."""
        rho = np.zeros(self.n)
        for phi, (s, e) in zip(self.orbitals, self.masks):
            idx = np.arange(s, e) % self.n
            rho[idx] += phi * phi
        return rho

    def reference_density(self):
        """Canonical density from tridiagonal eigen solve."""
        main = np.ones(self.n)/(self.dx**2) + self.V
        off = -0.5 * np.ones(self.n-1)/(self.dx**2)
        evals, evecs = eigh_tridiagonal(main, off, select="i", select_range=(0, self.n_orb-1))
        rho = np.sum(evecs[:, :self.n_orb]**2, axis=1)
        return rho, evals

    def orthogonalize_momentum(self, n_iter=8, damping=0.4, bmix=0.5, track_err=False):
        """
        Constraint Solver: Projects orbitals onto the orthogonal manifold.
        Uses a Jacobi update with Momentum (bmix) to accelerate exclusion.
        """
        err_history = [] if track_err else None
        for _ in range(n_iter):
            new_coeffs = []
            for i in range(len(self.orbitals)):
                si, ei = self.masks[i]
                mi = np.arange(si, ei) % self.n
                
                # Push away from neighbors
                correction = np.zeros(self.support_size)
                for j in range(len(self.orbitals)):
                    if i == j: continue
                    sj, ej = self.masks[j]
                    mj = np.arange(sj, ej) % self.n
                    
                    shared, idx_i, idx_j = np.intersect1d(mi, mj, return_indices=True)
                    if len(shared) > 0:
                        overlap = np.dot(self.orbitals[i][idx_i], self.orbitals[j][idx_j])
                        correction[idx_i] += overlap * self.orbitals[j][idx_j]
                
                # Update and Normalize
                phi_new = self.orbitals[i] - damping * correction
                phi_new /= (np.linalg.norm(phi_new) + 1e-12)
                
                # Apply Momentum Mixing
                delta = phi_new - self.orbitals[i]
                final_phi = phi_new + bmix * self.prev_delta[i]
                
                # Update buffers
                self.prev_delta[i] = delta
                new_coeffs.append(final_phi)
            
            self.orbitals = new_coeffs
            if track_err:
                err_history.append(self.overlap_error())
        return err_history

    def shift_windows(self):
        """Move the support window to track the orbital center of mass."""
        for i in range(len(self.orbitals)):
            phi = self.orbitals[i]
            # Weights based on density (phi^2)
            com = np.sum(np.arange(self.support_size) * phi**2) / (np.sum(phi**2) + 1e-12)
            shift = int(np.round(com - self.support_size/2))
            
            if shift != 0:
                self.masks[i][0] += shift
                self.masks[i][1] += shift
                self.orbitals[i] = np.roll(self.orbitals[i], -shift)
                if shift > 0: self.orbitals[i][-shift:] = 0
                else: self.orbitals[i][:-shift] = 0

# --- Runner Script ---
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters",    type=int,    default=1000  )
    ap.add_argument("--support",  type=int,    default=40   )
    ap.add_argument("--bmix",     type=float,  default=0.5  )
    ap.add_argument("--n_grid",   type=int,    default=100  )
    ap.add_argument("--l_max",    type=float,  default=10.0 )
    ap.add_argument("--n_orb",    type=int,    default=4    )
    ap.add_argument("--dt",       type=float,  default=0.05 )
    ap.add_argument("--dt_max",   type=float,  default=0.3  )
    ap.add_argument("--orth_iter",type=int,    default=4    )
    ap.add_argument("--orth_damp",type=float,  default=0.5  )
    ap.add_argument("--shift_every",type=int,  default=5    )
    ap.add_argument("--report_every",type=int, default=1   )
    ap.add_argument("--verbosity", type=int, default=3)
    ap.add_argument("--no_plot", action="store_true")
    args = ap.parse_args()

    solver = OMM_Solver(n_grid=args.n_grid, l_max=args.l_max, n_orb=args.n_orb, support_size=args.support)
    fire = FIRE([len(o) for o in solver.orbitals], dt=args.dt, dt_max=args.dt_max)

    plt.ion()
    fig, ax = plt.subplots(figsize=(10, 5))
    ov_track = []

    rho_ref, evals_ref = solver.reference_density()

    for it in range(args.iters):
        # 1. Relax Energy
        forces, energy = solver.get_forces()
        fire.step(solver.orbitals, forces)
        
        # 2. Hard Orthogonality Constraint (Projective Dynamics)
        err_hist = solver.orthogonalize_momentum(
            n_iter=args.orth_iter,
            damping=args.orth_damp,
            bmix=args.bmix,
            track_err=args.verbosity >= 3,
        )
        if args.verbosity >= 3 and err_hist:
            print("#ORTH step_err", " ".join(f"{v:.3e}" for v in err_hist))
        
        # 3. Adapt Geometry
        if args.shift_every > 0 and it % args.shift_every == 0:
            solver.shift_windows()

        if args.verbosity >= 3:
            ov_track.append((it, solver.overlap_error()))

        if args.report_every > 0 and it % args.report_every == 0:
            Etot, Ekin, Epot = solver.energy_terms()
            ov_err = solver.overlap_error()
            rho_omm = solver.density()
            err_rho = np.linalg.norm(rho_omm - rho_ref)
            if args.verbosity >= 3 and len(ov_track) > 0:
                xs, ys = zip(*ov_track)
                print(f"#ORTH window it[{xs[0]}..{xs[-1]}] min={min(ys):.3e} max={max(ys):.3e} last={ys[-1]:.3e}")
                ov_track.clear()
                if err_hist:
                    print("#ORTH step_err", " ".join(f"{v:.3e}" for v in err_hist))
            print(f"#ITER {it:04d}  E={Etot:.6f}  T={Ekin:.6f}  V={Epot:.6f}  ov_err={ov_err:.3e}  |rho-rho_ref|={err_rho:.3e}  FIRE_dt={fire.dt:.4f}")
            if not args.no_plot:
                ax.clear()
                ax.axhline(0.0, color="gray", lw=0.6, ls="--")
                ax.plot(solver.x, solver.V * 0.05, 'k--', alpha=0.2, label="V(x)")
                ax.plot(solver.x, rho_omm, 'k-', lw=1.5, label="Density OMM")
                ax.plot(solver.x, rho_ref, 'r--', lw=1.0, label="Density ref")
                for i, (phi, (s, e)) in enumerate(zip(solver.orbitals, solver.masks)):
                    grid_idx = np.arange(s, e) % solver.n
                    ax.plot(solver.x[grid_idx], phi, label=f"Orb {i}")
                ax.set_ylim(-0.5, 1.5)
                ax.set_title(f"FIRE+PD Iter {it}")
                ax.legend(loc='upper right', ncol=4, fontsize=8)
                plt.pause(0.01)

    if not args.no_plot:
        plt.ioff()
        plt.show()