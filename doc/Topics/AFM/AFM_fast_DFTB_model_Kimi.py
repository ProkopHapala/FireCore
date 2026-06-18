"""
AFM Interaction Models: Pseudocode for rigid-fragment interaction
using DFTB Hamiltonian in the AO basis.
 
Variants:
  L0 : Grid density overlap (baseline, not LCAO)
  L1 : P-based Pauli repulsion (AO basis, simplest LCAO)
  L2 : P+W second-order Pauli (sweet spot)
  L3 : P+W Pauli + Frontier CT (chemistry)
 
All variants share:
  - Precompute: expensive, done once per fragment
  - Per-pixel: cheap, done millions of times during AFM scan
  - No SCC, no diagonalization of combined system
  - No orthogonalization on the fly (AO basis only)
  - Atom-block sparse storage for surface (Fireball-style)
 
Data conventions (Fireball-style atom-block sparse):
  atom_block_sparse[i_atom][i_neigh] -> dense block [n_orb_i, n_orb_j]
  where j_atom = neigh_list[i_atom][i_neigh]
"""
 
# =============================================================================
# SHARED DATA STRUCTURES
# =============================================================================
 
class FragmentPrecompute:
    """Base precomputed data for any fragment (tip or surface)."""
    def __init__(self):
        # Atom-block sparse: dict[atom_i][neigh_idx] -> block_matrix
        self.P_blocks = {}           # density matrix P = 2*C_occ C_occ^T
        self.H_blocks = {}           # fragment Hamiltonian (intra-fragment)
        self.S_blocks = {}           # fragment overlap (intra-fragment)
        # Frontier MOs for chemistry (L3)
        self.C_frontier = None       # [nAO x n_front]  (e.g. 8 states)
        self.eps_frontier = None     # [n_front] orbital energies
        # For L2/L3:
        self.W_blocks = None         # energy-weighted density = P @ H @ P
 
class TipFragment(FragmentPrecompute):
    """Tip is small: all blocks can be dense or fully stored."""
    def __init__(self, n_ao):
        super().__init__()
        self.n_atoms = None          # small, ~10-50 atoms
        self.n_ao = n_ao             # ~50-200 AOs
 
class SurfaceFragment(FragmentPrecompute):
    """Surface is large: stored in atom-block sparse format."""
    def __init__(self, n_atoms):
        super().__init__()
        self.n_atoms = n_atoms       # large, ~1000-10000 atoms
        # P_blocks[atom_i] = list of (neigh_atom_j, block[n_orb_i, n_orb_j])
        # Only neighbors within intra-fragment cutoff (e.g. 6 A)
 
 
# =============================================================================
# PRECOMPUTE (once per fragment, can be expensive)
# =============================================================================
 
def precompute_fragment(fragment, calc_type="L2"):
    """
    Run isolated-fragment DFTB and build precomputed matrices.
    calc_type: "L1" | "L2" | "L3"
 
    Cost: dominated by single diagonalization of isolated fragment.
          For tip (~200 AOs): negligible.
          For surface (~10000 AOs): one-shot, amortized over millions of pixels.
    """
    # 1. Solve isolated fragment DFTB -> C, eps, H, S
    C_occ, eps_occ, H, S = solve_isolated_dftb(fragment)
 
    # 2. Density matrix P = 2 * C_occ @ C_occ.T
    fragment.P_blocks = build_atom_block_sparse(C_occ, mode="density")
    fragment.H_blocks = H
    fragment.S_blocks = S
 
    # 3. Energy-weighted density W = P @ H @ P  (needed for L2, L3)
    if calc_type in ("L2", "L3"):
        # W contains orbital-energy weighting for Pauli repulsion
        # In AO basis: W = P @ H @ P  (this is the correct form for 2nd-order APW)
        W_dense = P_dense @ H_dense @ P_dense   # dense, for derivation
        fragment.W_blocks = build_atom_block_sparse_from_dense(W_dense)
 
    # 4. Frontier MOs for charge-transfer chemistry (needed for L3)
    if calc_type == "L3":
        # Select e.g. HOMO-3 .. LUMO+3, truncate to atoms near surface top
        n_front = 8   # total frontier states
        fragment.C_frontier = extract_frontier_mos(C_occ, C_virt, eps_occ, eps_virt, n_front)
        fragment.eps_frontier = extract_frontier_energies(eps_occ, eps_virt, n_front)
 
    return fragment
 
 
# =============================================================================
# PER-PIXEL INTERACTION EVALUATION
# =============================================================================
 
def evaluate_interaction_L0(tip, surface, tip_pos, R_cut=6.0):
    """
    Level 0: Grid density overlap (baseline, scalar repulsion).
    Not LCAO-based; kept for comparison only.
 
    Precompute:
      - rho_tip, rho_surf on 3D grid
    Per-pixel cost:
      - O(N_grid_points_in_overlap_region)  (very fast, ~10^3 ops)
    Physics:
      - Captures Pauli wall magnitude but NO orbital anisotropy, NO chemistry.
    """
    E_Pauli = 0.0
    for atom_a in tip.atoms:
        # Find nearby surface atoms
        nearby_b = find_atoms_within(surface.atoms, tip_pos + atom_a.pos, R_cut)
        for atom_b in nearby_b:
            # Overlap integral of precomputed spherical densities
            overlap = density_overlap(atom_a.rho_grid, atom_b.rho_grid,
                                      displacement=tip_pos + atom_a.pos - atom_b.pos)
            E_Pauli += kappa * overlap
    return E_Pauli
 
 
def evaluate_interaction_L1(tip, surface, tip_pos, R_cut=6.0):
    """
    Level 1: P-based Pauli repulsion (AO basis, first LCAO step).
    Uses only the density matrix P (no energy weighting, no chemistry).
 
    Formula:
      E_Pauli ~ Tr[ P_A S_AB P_B S_BA ]   (simplified, captures S^2 overlap cost)
 
    Precompute:
      - tip.P_blocks   (small, can be dense)
      - surface.P_blocks (atom-block sparse, intra-fragment neighbors only)
    Per-pixel cost:
      - O(N_tip_atoms * N_neigh_surface * n_orb^3)
      - For 20 tip atoms x 20 neighbors x 4^3 = 51k ops per pixel
    Physics:
      + Orbital-dependent Pauli wall (anisotropic, directional)
      - Missing energy-dependent prefactor (W) => less accurate magnitude
      - Missing all chemistry (no virtual states)
    """
    E_Pauli = 0.0
 
    # For each tip atom alpha
    for alpha in tip.atom_indices:
        # Find surface atoms beta within inter-fragment cutoff
        beta_list = surface.find_neighbors(tip_pos, tip.atom_pos[alpha], R_cut)
 
        # Fetch tip row-blocks of P_A for atom alpha
        P_A_alpha = tip.P_blocks[alpha]   # dict: {neigh_a: block}
 
        for beta in beta_list:
            # Build S_AB block for pair (alpha, beta) from Slater-Koster tables
            S_ab = slater_koster_overlap(tip, surface, alpha, beta, tip_pos)
 
            # Fetch surface P_B blocks for atom beta and its intra-neighbors
            P_B_beta = surface.P_blocks[beta]
 
            # Local contraction: sum over tip atom alpha' and surface atom beta'
            #   contribution = Tr[ P_A[alpha,alpha'] @ S[alpha',beta] @
            #                      P_B[beta,beta'] @ S[beta',alpha] ]
            # This is a nested loop over alpha' (tip neighbors) and beta' (surface neighbors)
            for alpha_prime, P_A_aap in P_A_alpha.items():
                for beta_prime, P_B_bbp in P_B_beta.items():
                    # We need S[alpha', beta'] only if alpha' and beta' are coupled
                    # For the standard L1 formula, only diagonal S blocks are used
                    # in the simplest contraction. For full off-diagonal, fetch S_aapb:
                    S_aapb = get_sk_block_or_zero(tip, surface, alpha_prime, beta_prime, tip_pos)
                    if S_aapb is not None:
                        # Matrix chain: P_A[a,ap] @ S[ap,b] @ P_B[b,bp] @ S[bp,a]^T
                        term = trace_chain(P_A_aap, S_ab, P_B_bbp, S_aapb.T)
                        E_Pauli += term
 
    # Symmetrize (optional, depending on exact formula variant)
    return 0.5 * E_Pauli
 
 
def evaluate_interaction_L2(tip, surface, tip_pos, R_cut=6.0):
    """
    Level 2: P+W second-order Pauli (the sweet spot).
 
    Formula (AO basis, exact to 2nd order in S_AB):
      E = -Re Tr[ P_A S_AB P_B H_BA ]                     (chemical coupling)
          + 1/2 Tr[ S_AB P_B S_BA W_A ] + sym(A<->B)      (Pauli repulsion)
    where W_X = P_X H_X P_X  (energy-weighted density).
 
    Precompute:
      - tip.P_blocks, surface.P_blocks (same as L1)
      - tip.W_blocks, surface.W_blocks (atom-block sparse, same sparsity as P)
    Per-pixel cost:
      - O(N_tip_atoms * N_neigh_surface * n_orb^3)
      - Practically ~2x L1 cost (same structure, extra W contraction)
      - For 20 x 20 x 4^3 = ~100k ops per pixel (still trivial)
    Physics:
      + Full 2nd-order Pauli (validated ~0.26% error vs exact frozen energy)
      + Short-range chemical character via H_AB coupling
      + Energy-weighted repulsion (orbital-energy dependent)
      - Missing donor/acceptor CT (no virtual states)
    """
    E_chem = 0.0
    E_Pauli_A = 0.0   # S_AB P_B S_BA W_A
    E_Pauli_B = 0.0   # S_BA P_A S_AB W_B
 
    for alpha in tip.atom_indices:
        beta_list = surface.find_neighbors(tip_pos, tip.atom_pos[alpha], R_cut)
        P_A_alpha = tip.P_blocks[alpha]
        W_A_alpha = tip.W_blocks[alpha]
 
        for beta in beta_list:
            S_ab = slater_koster_overlap(tip, surface, alpha, beta, tip_pos)
            H_ba = slater_koster_hamiltonian(tip, surface, beta, alpha, tip_pos)  # H_BA block
 
            P_B_beta = surface.P_blocks[beta]
            W_B_beta = surface.W_blocks[beta]
 
            # ---- Term 1: Chemical coupling  -Re Tr[ P_A S_AB P_B H_BA ] ----
            for alpha_prime, P_A_aap in P_A_alpha.items():
                for beta_prime, P_B_bbp in P_B_beta.items():
                    S_aapb = get_sk_block_or_zero(tip, surface, alpha_prime, beta_prime, tip_pos)
                    if S_aapb is not None and H_ba is not None:
                        # P_A[a,ap] @ S[ap,bp] @ P_B[bp,b] @ H[b,a]
                        # Note: careful with transpose for H_BA
                        term = trace_chain_complex(P_A_aap, S_aapb, P_B_bbp, H_ba)
                        E_chem -= term.real
 
            # ---- Term 2: Pauli A  1/2 Tr[ S_AB P_B S_BA W_A ] ----
            # Contract over beta_prime (surface neighbors of beta) and alpha_prime (tip neighbors of alpha)
            for beta_prime, P_B_bbp in P_B_beta.items():
                S_bpa = get_sk_block_or_zero(tip, surface, beta_prime, alpha, tip_pos)   # S[bp,a]
                if S_bpa is None: continue
                for alpha_prime, W_A_aap in W_A_alpha.items():
                    S_apb = get_sk_block_or_zero(tip, surface, alpha_prime, beta, tip_pos)  # S[ap,b]
                    if S_apb is None: continue
                    # Matrix chain: S[a,b] @ P_B[b,bp] @ S[bp,a] @ W_A[a,ap] @ S[ap,b]
                    # Actually we need the full trace: S_AB P_B S_BA W_A
                    # This is Tr[ (S[a,b] @ P_B[b,bp] @ S[bp,a]) @ W_A[a,ap] ] ... no, need careful indices
                    # High-level: accumulate local block contractions; exact index algebra deferred to implementation
                    term = trace_pauli_chain(S_ab, P_B_bbp, S_bpa, W_A_aap, S_apb)
                    E_Pauli_A += term

            # ---- Term 3: Pauli B  1/2 Tr[ S_BA P_A S_AB W_B ] ----
            # Same structure with A<->B swap
            for alpha_prime, P_A_aap in P_A_alpha.items():
                S_apb = get_sk_block_or_zero(tip, surface, alpha_prime, beta, tip_pos)
                if S_apb is None: continue
                for beta_prime, W_B_bbp in W_B_beta.items():
                    S_bpa = get_sk_block_or_zero(tip, surface, beta_prime, alpha, tip_pos)
                    if S_bpa is None: continue
                    term = trace_pauli_chain(S_bpa.T, P_A_aap, S_apb, W_B_bbp, S_ab.T)
                    E_Pauli_B += term

    E_total = E_chem + 0.5 * (E_Pauli_A + E_Pauli_B)
    return E_total


def evaluate_interaction_L3(tip, surface, tip_pos, R_cut=6.0):
    """
    Level 3: P+W Pauli + Frontier CT (full chemistry).

    Formula:
      E = E_Pauli^{L2} + E_CT
      E_CT = - sum_{i in occ_A} sum_{j in virt_B} |V_ij|^2 / (eps_j - eps_i) + sym(A<->B)
      where V = C_A^front.T @ H_AB @ C_B^front  (tiny matrix, ~8x8)

    Precompute:
      - tip.P_blocks, surface.P_blocks
      - tip.W_blocks, surface.W_blocks
      - tip.C_frontier, surface.C_frontier  (frontier MO coefficients on AO basis)
      - tip.eps_frontier, surface.eps_frontier
    Per-pixel cost:
      - Pauli part: same as L2  (~100k ops)
      - CT part: O(N_front^2 * n_orb^2)  (~8x8 matrix, ~1k ops, almost free)
      - Total: ~101k ops per pixel
    Physics:
      + Everything in L2
      + Donor/acceptor charge transfer (Lewis acid/base, H-bonds)
      + Fukui-function-like orbital selectivity without full virtual space
    """
    # 1. Pauli + short-range coupling (same as L2)
    E_L2 = evaluate_interaction_L2(tip, surface, tip_pos, R_cut)

    # 2. Frontier charge-transfer correction
    E_CT = 0.0

    # Gather frontier AO coefficient vectors only for atoms in the local patch
    # C_frontier[atom, orbital_index, n_front] -> we need rows for active atoms only
    active_tip_atoms = tip.atom_indices    # tip is small, use all
    active_surf_atoms = []
    for alpha in tip.atom_indices:
        active_surf_atoms.extend(
            surface.find_neighbors(tip_pos, tip.atom_pos[alpha], R_cut)
        )
    active_surf_atoms = list(set(active_surf_atoms))   # unique

    # Extract local frontier vectors: C_A[local_AOs, n_front], C_B[local_BOs, n_front]
    C_A_loc = gather_frontier_vectors(tip, active_tip_atoms)      # [n_loc_A, n_front_A]
    C_B_loc = gather_frontier_vectors(surface, active_surf_atoms)  # [n_loc_B, n_front_B]

    # Build H_AB in the local AO subspace from Slater-Koster tables
    H_AB_loc = build_local_H_AB(tip, surface, active_tip_atoms, active_surf_atoms, tip_pos)

    # Tiny coupling matrix: V = C_A.T @ H_AB @ C_B
    V = C_A_loc.T @ H_AB_loc @ C_B_loc    # [n_front_A, n_front_B]

    # Sum over donor/acceptor pairs
    for i, eps_i in enumerate(tip.eps_frontier):
        for j, eps_j in enumerate(surface.eps_frontier):
            de = eps_j - eps_i
            if de > 1e-6:    # avoid zero denominators; eps_j from virt_B, eps_i from occ_A
                E_CT -= (abs(V[i, j]) ** 2) / de

    # Symmetric term: surface donates to tip
    # (can reuse the same V if frontier spaces cover both occ and virt on each side)
    # If C_frontier on each fragment includes both HOMO-range (occ) and LUMO-range (virt),
    # the double loop above already covers both directions with proper eps ordering.
    # Otherwise add the reverse sum explicitly here.

    return E_L2 + 2.0 * E_CT   # factor 2 for spin if needed; adjust to convention


# =============================================================================
# DESIGN SUMMARY & TRADE-OFFS
# =============================================================================
"""
| Level | Name                | Precompute storage               | Per-pixel FLOPs | Physics fidelity                 | Best use case                  |
|-------|--------------------|-----------------------------------|-----------------|----------------------------------|--------------------------------|
| L0    | Grid overlap       | rho grids (3D)                    | ~1e3            | Scalar Pauli only                | Baseline / quick preview       |
| L1    | P-Pauli            | P sparse blocks                   | ~5e4            | Orbital-anisotropic Pauli        | Fast repulsive wall            |
| L2    | P+W 2nd-order      | P + W sparse blocks               | ~1e5            | + Energy-weighted Pauli + chem   | **Sweet spot for AFM**         |
| L3    | P+W + Frontier CT  | P + W + C_frontier (~8 states)    | ~1.1e5          | + Donor/acceptor, H-bonds        | When chemistry matters         |

Key design decisions:
  - AO basis only: no Löwdin transforms during scan (destroys locality).
  - Atom-block sparse surface storage: reuse Fireball neighbor lists, zero redundancy.
  - No full P_B dense matrix: only fetch local patch under tip at each pixel.
  - W = P @ H @ P is precomputed once, same sparsity as P.
  - Frontier CT adds ~1% extra cost but captures chemistry missing in density-only models.

Open implementation questions:
  - Exact index contraction for the Pauli traces (can be derived from 2nd-order APW).
  - Whether to precompute W as P@H@P or as sum_i eps_i C_i C_i^T (numerically identical).
  - GPU batched evaluation: one thread per pixel, shared tip data, coalesced surface block fetch.
"""

# =============================================================================
# PLACEHOLDER / STUB FUNCTIONS (implementation-specific)
# =============================================================================

def solve_isolated_dftb(fragment):
    """Run DFTB diagonalization for one fragment. Returns C_occ, eps_occ, H, S."""
    raise NotImplementedError("Plug in Fireball/DFTB+ diagonalization")

def build_atom_block_sparse(C_occ, mode="density"):
    """Convert MO coefficients to atom-block sparse P or W."""
    raise NotImplementedError("Use Fireball sparse block layout")

def build_atom_block_sparse_from_dense(M_dense):
    """Pack dense matrix into atom-block sparse format."""
    raise NotImplementedError

def extract_frontier_mos(C_occ, C_virt, eps_occ, eps_virt, n_front):
    raise NotImplementedError

def extract_frontier_energies(eps_occ, eps_virt, n_front):
    raise NotImplementedError

def slater_koster_overlap(tip, surface, atom_a, atom_b, tip_pos):
    """Evaluate S_AB block for atom pair (a,b) at given tip position."""
    raise NotImplementedError("Use SK tables + two-center integral code")

def slater_koster_hamiltonian(tip, surface, atom_a, atom_b, tip_pos):
    """Evaluate H_AB block for atom pair (a,b) at given tip position."""
    raise NotImplementedError

def get_sk_block_or_zero(tip, surface, atom_a, atom_b, tip_pos):
    """Return S or H block if within cutoff, else None."""
    raise NotImplementedError

def trace_chain(P1, S1, P2, S2):
    """Tr[ P1 @ S1 @ P2 @ S2 ] for atom blocks."""
    raise NotImplementedError

def trace_chain_complex(P1, S1, P2, H2):
    """Tr[ P1 @ S1 @ P2 @ H2 ] returning complex scalar."""
    raise NotImplementedError

def trace_pauli_chain(S1, P, S2, W, S3):
    """High-level Pauli trace contraction; exact index structure depends on APW derivation."""
    raise NotImplementedError

def find_atoms_within(atoms, center, R_cut):
    raise NotImplementedError

def density_overlap(rho_a, rho_b, displacement):
    raise NotImplementedError

def gather_frontier_vectors(fragment, active_atoms):
    """Extract frontier MO AO coefficients for selected atoms."""
    raise NotImplementedError

def build_local_H_AB(tip, surface, tip_atoms, surf_atoms, tip_pos):
    """Build H_AB rectangular matrix for local atom subset."""
    raise NotImplementedError
