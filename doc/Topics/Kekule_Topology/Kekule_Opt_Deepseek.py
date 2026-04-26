#!/usr/bin/env python3
"""
Self‑consistent Peierls (SSH) model for zigzag graphene nanoribbons
with chemically functionalised edges.

Usage examples:
  python ribbon_pcet.py --width 4 --length 10 --gamma 1.2 --delta 0.4 --top "O,OH,O,C" --bottom "C,C,C,C"
  python ribbon_pcet.py --width 6 --length 8 --gamma 1.5 --t0 2.8 --top "N,N,N,N,N,N,N,N" --bottom "O,O,O,O,O,O,O,O"
"""

import argparse
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# =============================================================================
# 1. Lattice generation for a zigzag nanoribbon
# =============================================================================
a_CC = 1.0                     # carbon‑carbon distance (length unit)
# Graphene lattice vectors (orthogonal basis)
a1 = np.array([np.sqrt(3), 0.0])
a2 = np.array([np.sqrt(3)/2, 1.5])

def build_zigzag_ribbon(width, length):
    """
    Build a zigzag graphene ribbon periodic along a2, finite along a1.

    Parameters
    ----------
    width : int        number of unit cells along a1 (defines ribbon width)
    length : int       number of unit cells along a2 (periodic direction)

    Returns
    -------
    positions : ndarray (N_atoms, 2)   Cartesian coordinates of all atoms
    sublattice : ndarray (N_atoms,)     0 for A, 1 for B
    edges : dict                         indices of atoms on top/bottom edge
           'top':    n1 == width-1 (the outermost zigzag row)
           'bottom': n1 == 0
    """
    A_pos, B_pos = [], []
    for n1 in range(width):          # finite direction
        for n2 in range(length):     # periodic direction
            rA = n1*a1 + n2*a2
            rB = rA + np.array([0.0, a_CC])
            A_pos.append(rA)
            B_pos.append(rB)
    A_pos = np.array(A_pos)
    B_pos = np.array(B_pos)
    positions = np.vstack([A_pos, B_pos])
    sublattice = np.array([0]*len(A_pos) + [1]*len(B_pos))

    # Identify edge atoms by their n1 index
    edge_top = []
    edge_bottom = []
    for i, r in enumerate(positions):
        # find n1 from the position: n1 = round( (r dot a1) / (a1.a1) )? Better:
        # Since we built from loops, we can record n1 per atom when generating.
        # More robust: determine n1 from coordinates by projecting onto a1.
        # a1 = (sqrt(3), 0), so n1 = x / sqrt(3)
        n1 = round(r[0] / a1[0])   # x coordinate integer multiple of sqrt(3)
        if n1 == 0:
            edge_bottom.append(i)
        elif n1 == width-1:
            edge_top.append(i)
    # The edges should each have exactly length atoms of one sublattice.
    return positions, sublattice, {'top': edge_top, 'bottom': edge_bottom}

# =============================================================================
# 2. Build bond list with periodic boundary conditions along a2
# =============================================================================
def build_bonds(positions, sublattice, length, width):
    """
    Find nearest‑neighbour (A–B) bonds.  Bonds that cross the periodic
    boundary are wrapped using the periodicity vector a2*length.
    """
    # Bond vectors from A to B
    delta1 = np.array([0.0, a_CC])
    delta2 = np.array([-np.sqrt(3)/2*a_CC, -0.5*a_CC])
    delta3 = np.array([ np.sqrt(3)/2*a_CC, -0.5*a_CC])
    delta_all = [delta1, delta2, delta3]
    # Bravais periodicity vector (only along a2 direction)
    periodic_vec = a2 * length
    bonds = []   # each: (i, j, bond_type)
    N = len(positions)
    for i in range(N):
        if sublattice[i] != 0:   # only from A atoms to avoid duplicates
            continue
        pos_i = positions[i]
        for d_idx, d in enumerate(delta_all):
            pos_j_ideal = pos_i + d
            # Find the real B atom closest to pos_j_ideal (with possible wrapping)
            best_j = -1
            best_dist = 1e9
            for j in range(N):
                if sublattice[j] != 1:
                    continue
                pos_j = positions[j]
                # try all possible shifts by multiples of periodic_vec
                for shift in [np.array([0.,0.]), periodic_vec, -periodic_vec]:
                    vec = pos_j + shift - pos_i
                    dist = np.linalg.norm(vec - d)
                    if dist < best_dist:
                        best_dist = dist
                        best_j = j
            if best_dist < 0.2 * a_CC:
                bonds.append((i, best_j, d_idx))
    # Now we have unique bonds.  Check that every A atom has 3 bonds (except edges maybe).
    return bonds

# =============================================================================
# 3. Parse edge functionalisation strings into pinning rules
# =============================================================================
def parse_edge_string(s, edge_atoms, length, default_onsite=0.0,
                      onsite_C=0.0, onsite_N=-0.5, onsite_O=-1.0,
                      onsite_OH=0.0, onsite_NH=-0.2):
    """
    Convert a string like "O,OH,N,C" into a dict mapping edge atom indices
    to (pinned_bond_strength_type, onsite_energy).
    """
    parts = [x.strip() for x in s.split(',')]
    if len(parts) == 1:
        # Broadcast single code to all edge atoms
        parts = parts * len(edge_atoms)
    elif len(parts) != len(edge_atoms):
        raise ValueError(f"Edge string length ({len(parts)}) does not match number of edge atoms ({len(edge_atoms)})")
    
    mapping = {}
    for i, code in enumerate(parts):
        atom_idx = edge_atoms[i]
        onsite = default_onsite
        pin = None
        if code == 'C':
            onsite = onsite_C
            pin = None
        elif code == 'O':
            onsite = onsite_O
            pin = 'single'   # makes adjacent C-C bonds weak
        elif code == 'OH':
            onsite = onsite_OH
            pin = 'double'   # makes adjacent C-C bonds strong
        elif code == 'N':
            onsite = onsite_N
            pin = None       # pyridinic N doesn't force alternation strongly
        elif code == 'NH':
            onsite = onsite_NH
            pin = 'double'   # amino nitrogen donates, strengthens adjacent bonds
        else:
            raise ValueError(f"Unknown edge code: {code}")
        mapping[atom_idx] = (pin, onsite)
    return mapping

# =============================================================================
# 4. SCF iteration: t_ij = t0 + gamma * P_ij
# =============================================================================
def scf_relaxation(positions, sublattice, bonds, width, length,
                   t0, gamma, delta, pinning_dict, onsite_dict,
                   max_iter=500, mix=0.3, tol=1e-6):
    """
    Self‑consistent bond order relaxation for the ribbon.
    
    pinning_dict: dict {atom_idx: (pin_type, onsite)}
       pin_type is None, 'single', or 'double'.
    onsite_dict: applies on-site energy shifts to particular atoms.
    """
    N = len(positions)
    M = len(bonds)
    # initialise hopping
    t = np.full(M, t0)
    # apply initial pinning to bonds that involve a pinned atom
    for idx, (i, j, btype) in enumerate(bonds):
        # check if either endpoint is pinned with a specific type
        pin_i = pinning_dict.get(i, (None,0))[0]
        pin_j = pinning_dict.get(j, (None,0))[0]
        pin_type = None
        if pin_i is not None:
            pin_type = pin_i
        if pin_j is not None:
            pin_type = pin_j   # if both, later define precedence
        if pin_type == 'single':
            t[idx] = t0 - delta
        elif pin_type == 'double':
            t[idx] = t0 + delta
        # else keep t0

    # list of bond indices that are pinned permanently
    pinned_bonds = set()
    for idx, (i, j, btype) in enumerate(bonds):
        for atom in (i, j):
            pin_type = pinning_dict.get(atom, (None,0))[0]
            if pin_type is not None:
                pinned_bonds.add(idx)
                break

    for iteration in range(max_iter):
        # Build Hamiltonian
        H = np.zeros((N,N))
        for idx, (i, j, _) in enumerate(bonds):
            H[i, j] = t[idx]
            H[j, i] = t[idx]      # real Hamiltonian
        # Apply on-site energies
        for atom, (_, eps) in pinning_dict.items():
            H[atom, atom] += eps
        # Recognise that some edge atoms may not be in pinning_dict but default.

        evals, evecs = eigh(H)
        n_occ = N // 2
        # Compute bond orders P_ij for all bonds
        P = np.zeros((N,N))
        for n in range(n_occ):
            c = evecs[:, n]
            P += 2.0 * np.outer(c, c)   # factor 2 for spin

        # Target hopping from SSH equation
        t_new = np.zeros(M)
        for idx, (i, j, _) in enumerate(bonds):
            P_ij = P[i, j]   # off-diagonal
            t_new[idx] = t0 + gamma * P_ij

        # Re-apply pinning (restore pinned bonds to their target values)
        for idx in pinned_bonds:
            i, j, btype = bonds[idx]
            # determine target value based on atom's pin type
            for atom in (i, j):
                pin_type = pinning_dict.get(atom, (None,0))[0]
                if pin_type is not None:
                    if pin_type == 'single':
                        t_new[idx] = t0 - delta
                    elif pin_type == 'double':
                        t_new[idx] = t0 + delta
                    break   # only one atom's pinning matters

        # Mixing
        t = (1 - mix) * t + mix * t_new
        max_diff = np.max(np.abs(t_new - t))
        if max_diff < tol:
            break
    else:
        print("  Warning: SCF not converged.")

    # Compute total energy
    E_pi = 2.0 * np.sum(evals[:n_occ])
    E_sigma = 0.0
    for idx in range(M):
        E_sigma += (t0 - t[idx])**2
    E_sigma *= 0.5 / gamma
    E_total = E_pi + E_sigma
    return t, E_total, evals, evecs, n_occ

# =============================================================================
# 5. LDOS and plotting utilities
# =============================================================================
def compute_ldos(evals, evecs, n_occ, sigma_broad=0.05):
    """Local density of states at E=0 (Fermi level)."""
    N_states = len(evals)
    weights = np.exp(- (evals - 0.0)**2 / (2*sigma_broad**2))
    weights /= (np.sqrt(2*np.pi)*sigma_broad)
    ldos = np.sum(evecs**2 * weights, axis=1)
    return ldos / ldos.max()

def plot_results(positions, sublattice, bonds, t_vals, ldos, title, width, length, pinning_dict=None):
    fig, axes = plt.subplots(1, 2, figsize=(14,6))
    periodic_vec = a2 * length
    
    # Color Scheme
    from matplotlib.colors import LinearSegmentedColormap
    custom_cmap = LinearSegmentedColormap.from_list("kekule", ["#f0f0d0", "#000080"])
    
    # 1. Bond strength map
    ax = axes[0]
    ax.set_facecolor('white')
    vmin, vmax = np.min(t_vals), np.max(t_vals)
    
    for idx, (i, j, _) in enumerate(bonds):
        p1 = positions[i]
        p2_orig = positions[j]
        # Find which shift was used (same logic as in build_bonds)
        best_shift = np.array([0., 0.])
        best_d = 1e9
        delta_all = [np.array([0.0, a_CC]), np.array([-np.sqrt(3)/2*a_CC, -0.5*a_CC]), np.array([ np.sqrt(3)/2*a_CC, -0.5*a_CC])]
        d_target = delta_all[bonds[idx][2]]
        for shift in [np.array([0.,0.]), periodic_vec, -periodic_vec]:
            d = np.linalg.norm(p2_orig + shift - p1 - d_target)
            if d < best_d:
                best_d = d
                best_shift = shift
        
        p2 = p2_orig + best_shift
        col = plt.cm.RdBu( (t_vals[idx] - vmin) / (vmax - vmin + 1e-12) )
        lw = 1.0 + 3.0 * (t_vals[idx] - vmin) / (vmax - vmin + 1e-12)
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=col, linewidth=lw, alpha=0.6, zorder=1)
        
        if np.any(best_shift): # Draw the other side of PBC bond
            p1_alt = p1 - best_shift
            p2_alt = p2_orig
            ax.plot([p1_alt[0], p2_alt[0]], [p1_alt[1], p2_alt[1]], color=col, linewidth=lw, alpha=0.6, zorder=1)

    # Plot atoms
    ax.scatter(positions[:,0], positions[:,1], c='lightgrey', s=30, zorder=5, edgecolors='k', linewidth=0.2)
    
    # Highlight functionalized atoms
    if pinning_dict:
        for idx, (pin, eps) in pinning_dict.items():
            color = 'red' if pin == 'single' else ('blue' if pin == 'double' else 'grey')
            if eps != 0: # Chemical marker
                ax.scatter(positions[idx,0], positions[idx,1], c=color, s=80, zorder=6, edgecolors='black')

    # Strict limits for clipping
    xmin, xmax = np.min(positions[:,0])-1, np.max(positions[:,0])+1
    ymin, ymax = np.min(positions[:,1])-1, np.max(positions[:,1])+1
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(f'Bond Alternation Pattern\n{title}')
    ax.set_aspect('equal')

    # 2. LDOS
    ax2 = axes[1]
    sc = ax2.scatter(positions[:,0], positions[:,1], c=ldos, cmap=custom_cmap, s=50, edgecolors='k', linewidth=0.3)
    plt.colorbar(sc, ax=ax2, label="LDOS @ Fermi")
    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(ymin, ymax)
    ax2.set_title('LDOS at E=0')
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()

# =============================================================================
# 6. Main entry point with CLI
# =============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SSH Kekulé ribbon solver")
    parser.add_argument('--width', type=int, default=4, help='number of unit cells across ribbon (zigzag rows)')
    parser.add_argument('--length', type=int, default=10, help='length of periodic unit cell (repetitions along edge)')
    parser.add_argument('--t0', type=float, default=2.8, help='average hopping (eV)')
    parser.add_argument('--gamma', type=float, default=1.2, help='electron-phonon coupling γ = 2α²/K (eV)')
    parser.add_argument('--delta', type=float, default=0.4, help='bond strength pinning amplitude (eV)')
    parser.add_argument('--top_edge', type=str, default="C", help='comma-separated edge codes for top edge (C,O,OH,N,NH)')
    parser.add_argument('--bottom_edge', type=str, default="C", help='edge codes for bottom edge')
    parser.add_argument('--onsite_C', type=float, default=0.0, help='on-site energy for carbon (eV)')
    parser.add_argument('--onsite_N', type=float, default=-0.5, help='on-site for pyridinic N')
    parser.add_argument('--onsite_O', type=float, default=-1.0, help='on-site for quinonic O')
    parser.add_argument('--onsite_OH', type=float, default=0.0, help='on-site for hydroxy carbon')
    parser.add_argument('--onsite_NH', type=float, default=-0.2, help='on-site for amine N')
    parser.add_argument('--mix', type=float, default=0.3, help='SCF mixing parameter')
    parser.add_argument('--tol', type=float, default=1e-6, help='SCF convergence tolerance')
    args = parser.parse_args()
    positions, sublattice, edges = build_zigzag_ribbon(args.width, args.length)
    N = len(positions)
    print(f"Ribbon width = {args.width}, length = {args.length}, total atoms = {N}")
    bonds = build_bonds(positions, sublattice, args.length, args.width)
    M = len(bonds)
    print(f"Number of bonds = {M}")

    # Parse edge configurations
    top_mapping = parse_edge_string(args.top_edge, edges['top'], args.length, default_onsite=0,
                                    onsite_C=args.onsite_C, onsite_N=args.onsite_N,
                                    onsite_O=args.onsite_O, onsite_OH=args.onsite_OH,
                                    onsite_NH=args.onsite_NH)
    bottom_mapping = parse_edge_string(args.bottom_edge, edges['bottom'], args.length, default_onsite=0,
                                       onsite_C=args.onsite_C, onsite_N=args.onsite_N,
                                       onsite_O=args.onsite_O, onsite_OH=args.onsite_OH,
                                       onsite_NH=args.onsite_NH)
    # Combine into a global pinning dict
    pinning_dict = {}
    pinning_dict.update(top_mapping)
    pinning_dict.update(bottom_mapping)
    print("Edge functionalisation loaded.")

    # Run SCF
    print("Running SCF...")
    t_final, E_total, evals, evecs, n_occ = scf_relaxation(
        positions, sublattice, bonds, args.width, args.length,
        t0=args.t0, gamma=args.gamma, delta=args.delta,
        pinning_dict=pinning_dict, onsite_dict=None,
        max_iter=500, mix=args.mix, tol=args.tol)
    print(f"Total energy = {E_total:.4f} eV")

    # Compute LDOS
    ldos = compute_ldos(evals, evecs, n_occ)
    # Plot
    plot_results(positions, sublattice, bonds, t_final, ldos,
                 f"width={args.width}, len={args.length}, E={E_total:.2f} eV",
                 args.width, args.length, pinning_dict=pinning_dict)
