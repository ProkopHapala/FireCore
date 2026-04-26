"""
KekuleSolver.py - General Tight-Binding Solver for Kekulé Topology
==================================================================
This module provides a general engine for solving hexagonal lattices 
with Kekulé distortions, periodic boundary conditions, and onsite energies.
"""

import numpy as np
from scipy.linalg import eigh

class KekuleSolver:
    def __init__(self, t0=2.8, delta=0.4, a_CC=1.0):
        # Physical parameters
        self.t0 = t0
        self.delta = delta
        self.a_CC = a_CC
        self.K_mag = 4 * np.pi / (3 * np.sqrt(3) * a_CC)
        self.K_vec = np.array([self.K_mag, 0.0])
        
        # Lattice and Periodicity
        self.Lx = 10.0
        self.Ly = 10.0
        self.pbc = (False, False) # (x, y)
        
        # Domain Settings
        self.domain_type = "Step" # Step, Circle, Rect
        self.phi1 = 0.0
        self.phi2 = 2*np.pi/3
        self.delta2 = delta
        self.patch_size = [5.0, 5.0]
        
        # Data
        self.pos = None
        self.sub = None
        self.H = None
        self.onsite = None
        self.bonds = []
        self.evals = None
        self.evecs = None

    def get_phi_delta(self, r):
        """Returns (phi, delta) for a given position r."""
        if self.domain_type == "Step":
            # Use patch_size[0] as the x-position of the boundary
            # Shifted so that 0 means left edge and Lx means right edge
            x_bound = self.patch_size[0] - self.Lx/2
            if r[0] < x_bound: return self.phi1, self.delta
            else: return self.phi2, self.delta2
        elif self.domain_type == "Circle":
            if np.linalg.norm(r) < self.patch_size[0]/2: return self.phi2, self.delta2
            else: return self.phi1, self.delta
        elif self.domain_type == "Rect":
            w, h = self.patch_size
            if abs(r[0]) < w/2 and abs(r[1]) < h/2:
                return self.phi2, self.delta2
            else:
                return self.phi1, self.delta
        return self.phi1, self.delta

    def set_onsite_at_boundary(self, width=1.0, energy=1.0):
        """Adds onsite energy (e.g. Nitrogen) near the grain boundary at x=0."""
        self.onsite = np.zeros(len(self.pos))
        for i in range(len(self.pos)):
            if abs(self.pos[i, 0]) < width:
                self.onsite[i] = energy

    def generate_lattice(self, Nx=None, Ny=None):
        """Generates a rectangular honeycomb lattice with exact snapping for PBC."""
        ax = np.sqrt(3) * self.a_CC
        ay = 3.0 * self.a_CC
        
        if Nx is None: Nx = int(round(self.Lx / ax))
        if Ny is None: Ny = int(round(self.Ly / ay))
        
        Nx = max(1, Nx); Ny = max(1, Ny)
        self.Lx = Nx * ax
        self.Ly = Ny * ay
        
        a1 = np.array([ax, 0.0])
        a2 = np.array([ax/2, 1.5 * self.a_CC])
        
        pos_A = []
        for n1 in range(Nx):
            for n2 in range(Ny * 2):
                r = n1*a1 + n2*a2
                r[0] = r[0] % self.Lx
                r[1] = r[1] % self.Ly
                pos_A.append(r)
        
        pos_A = np.unique(np.round(pos_A, 6), axis=0)
        pos_B = pos_A + np.array([0.0, self.a_CC])
        pos_B[:, 0] = pos_B[:, 0] % self.Lx
        pos_B[:, 1] = pos_B[:, 1] % self.Ly
        pos_B = np.unique(np.round(pos_B, 6), axis=0)
        
        self.pos = np.vstack([pos_A, pos_B])
        self.pos[:,0] -= self.Lx/2
        self.pos[:,1] -= self.Ly/2
        
        self.sub = np.array([0]*len(pos_A) + [1]*len(pos_B))
        self.onsite = np.zeros(len(self.pos))
        return self.pos, self.sub

    def build_hamiltonian(self):
        N = len(self.pos)
        self.H = np.zeros((N, N))
        self.bonds = []
        
        delta_vecs = [
            np.array([0.0, 1.0]) * self.a_CC,
            np.array([-np.sqrt(3)/2, -0.5]) * self.a_CC,
            np.array([ np.sqrt(3)/2, -0.5]) * self.a_CC
        ]
        Kdot_delta = [np.dot(self.K_vec, d) for d in delta_vecs]
        
        shifts = [np.array([0,0])]
        if self.pbc[0]: shifts += [np.array([self.Lx, 0]), np.array([-self.Lx, 0])]
        if self.pbc[1]: shifts += [np.array([0, self.Ly]), np.array([0, -self.Ly])]
        if self.pbc[0] and self.pbc[1]:
            shifts += [np.array([self.Lx, self.Ly]), np.array([self.Lx, -self.Ly]),
                       np.array([-self.Lx, self.Ly]), np.array([-self.Lx, -self.Ly])]

        tol = 0.1 * self.a_CC
        for i in range(N):
            if self.sub[i] != 0: continue
            p_i = self.pos[i]
            for j in range(N):
                if self.sub[j] != 1: continue
                p_j_base = self.pos[j]
                for s in shifts:
                    p_j = p_j_base + s
                    vec = p_j - p_i
                    dist = np.linalg.norm(vec)
                    if abs(dist - self.a_CC) < tol:
                        best_k = np.argmax([np.dot(vec, d)/(dist*self.a_CC) for d in delta_vecs])
                        r_mid = (p_i + p_j) / 2
                        phi, d_loc = self.get_phi_delta(r_mid)
                        dt = d_loc * np.cos(Kdot_delta[best_k] + phi)
                        t = self.t0 + dt
                        self.H[i, j] = -t
                        self.H[j, i] = -t
                        self.bonds.append((i, j, t, r_mid, best_k, s))

        for i in range(N): self.H[i, i] = self.onsite[i]

    def solve(self):
        self.evals, self.evecs = eigh(self.H)
        return self.evals, self.evecs

    def get_ldos(self, emin, emax):
        mask = (self.evals >= emin) & (self.evals <= emax)
        if not np.any(mask): return np.zeros(len(self.pos))
        return np.sum(np.abs(self.evecs[:, mask])**2, axis=1)

    def get_dos(self, energy_range, sigma=0.05):
        dos = np.zeros_like(energy_range)
        for e in self.evals:
            dos += np.exp(-0.5 * ((energy_range - e) / sigma)**2)
        return dos / (len(self.pos) * np.sqrt(2*np.pi) * sigma)
