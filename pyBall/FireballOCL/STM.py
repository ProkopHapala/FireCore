"""
STM Simulator Module — Shared functions for Phase 1.

Provides:
- DOS computation via Fireball SCF (spectral function with Gamma broadening)
- Wolfsberg-Helmholtz hopping model
- STM current calculation
- 1D/2D scanning
- PDOS real-space projection (GridProjector)
- Plotting utilities

Usage:
  from pyBall.FireballOCL.STM import compute_dos, compute_hopping, stm_current, ...
"""

"""
Shared STM transport module for Fireball-based molecular junction simulation.

PURPOSE:
  Provides common functions for STM (Scanning Tunneling Microscopy) analysis including:
  - DOS computation with Γ broadening
  - Hamiltonian assembly from Fireball sparse format
  - NEGF transport (Caroli formula, iterative solvers, response metrics)
  - Inter-system coupling (vacuum exponential + Slater-Koster angular)
  - MO overlap diagnostics for sum-over-states validation
  - PDOS real-space projection

PHYSICAL BACKGROUND:

1) Spectral Function with Γ Broadening:
   The retarded Green's function for a system coupled to leads is:
     G^r(E) = [(E + iη)S - H - Σ]^{-1}
   where Σ = -iΓ is the wide-band limit self-energy on contact atoms.
   The spectral function (projected DOS) is:
     A(E) = (1/π) Im[G^r(E)]
   For each molecular orbital n with energy ε_n and coupling Γ_n:
     A_n(E) = (Γ_n/π) / [(E - ε_n)² + Γ_n²]

2) NEGF Caroli Formula (Landauer-Büttiker):
   For a two-terminal junction with leads L and R:
     G^r(E) = [(E + iη)S - H - Σ_L - Σ_R]^{-1}
     T(E) = Tr[Γ_L(E) G^r(E) Γ_R(E) G^a(E)]
     I = (2e/h) ∫ T(E) [f(E-μ_L) - f(E-μ_R)] dE
   where Γ_{L/R} = i(Σ_{L/R} - Σ_{L/R}†) are the broadening matrices.

3) Iterative Solver (GMRES + Hutchinson Trace):
   To avoid O(N³) inversion, use iterative linear solver:
     A x = b, where A = (E+iη)S - H - Σ_L - Σ_R
   Hutchinson stochastic trace estimator for transmission:
     T(E) ≈ (1/N_rand) Σ_{k=1}^{N_rand} z_k† Γ_L G^r Γ_R G^a z_k
   where z_k are random Rademacher vectors (±1 entries).

4) Deterministic Response Metric:
   Single-solve probe of coupling strength:
     A = (E+iη)S - H - Σ_L - Σ_R
     b = source_vector (e.g., unit vector on lead L)
     x = solve(A, b)
     resp = x† Γ_measure x
   This is not a transmission coefficient but correlates with it.

5) MO Overlap (Sum Over States):
   For a "featureless metallic tip" approximation:
     v_smp(s) = Σ_{i∈tip} H_{TS}(i,:)  # effective coupling vector
     t_n(s) = v_smp(s) · C_{:,n}        # overlap with MO n
     I(s) = |t_n(s)|²
   This is the "poles of Green's function" picture:
     T(E) ≈ Σ_n |t_n(E)|² L_n(E), where L_n(E) is Lorentzian lineshape.

6) Inter-System Coupling (Vacuum Exponential + SK Angular):
   Tunneling matrix elements between tip and sample:
     H_μν(r) = A_μν * exp(-β(r - r₀)) * SK_angular(l,m,n)
   where:
     - exp(-β(r - r₀)): vacuum radial decay (β ~1.0 Å⁻¹, r₀ ~3.0 Å)
     - SK_angular: Slater-Koster factors (l,m,n direction cosines)
     - A_μν: orbital-pair strengths (A_ss, A_sp, A_ppσ, A_ppπ)
   This replaces Fireball's covalent-bond radial functions with physically
   appropriate vacuum tunneling decay while retaining correct angular symmetry.

7) Float4 Per Atom Structure:
   For compatibility with GridFF/OCL_Hamiltonian, we pad to exactly 4 orbitals
   per atom: (s, px, py, pz). For atoms with fewer orbitals (e.g., H with only s),
   padded orbitals are decoupled by setting H_ii = E_pad (large) and S_ii = 1.

KEY FUNCTIONS:
  - compute_dos(): Fireball SCF + spectral function with Γ broadening
  - build_dense_HS(): Convert sparse Fireball H/S to dense matrices
  - build_inter_system_blocks_exp_sk(): Vacuum exponential + SK coupling
  - assemble_combined_HS(): Build tip+sample Hamiltonian for NEGF
  - negf_current(): Direct Caroli transmission
  - negf_current_iterative(): GMRES + Hutchinson trace
  - negf_response_at_energy(): Deterministic response metric
  - coupling_vec_tip_to_sample(): Featureless tip coupling vector
  - mo_overlap_amplitude(): MO overlap intensity
  - find_homo_from_eps(): Identify HOMO from eigenvalues
  - pad_HS_to_float4(): Convert to float4-per-atom structure
  - project_pdos_to_grid(): Real-space PDOS projection

USAGE EXAMPLE:
  from pyBall.FireballOCL.STM import compute_dos, negf_response_at_energy

  # Compute DOS for tip and sample
  tip = compute_dos(tip_types, tip_pos, gamma=0.1, out_path="dos_tip.npz")
  smp = compute_dos(smp_types, smp_pos, gamma=0.1, out_path="dos_smp.npz")

  # Build inter-system coupling
  pairs, Hb, Sb = build_inter_system_blocks_exp_sk(
      tip_types, tip_pos, tip['norb_per'],
      smp_types, smp_pos, smp['norb_per'],
      beta=1.0, r0=3.0, A_ss=-1.0, A_sp=-1.0
  )

  # Assemble combined system
  Hc, Sc = assemble_combined_HS(tip['H'], tip['S'], ..., smp['H'], smp['S'], ...)

  # Compute NEGF response at Fermi level
  resp = negf_response_at_energy(E_F, Hc, Sc, leadL, leadR, gammaL=2.0, gammaR=2.0)

NOTES:
  - All energies in eV
  - All positions in Å
  - Default export directory: tests/pyFireball/export/stm_phase1/
"""

import os, sys, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from numpy.linalg import solve

# ── Constants ─────────────────────────────────────────────────────────────────────
ELEM_Z  = {'H':1, 'C':6, 'N':7, 'O':8, 'S':16}
Z_ELEM  = {v:k for k,v in ELEM_Z.items()}
RC_ANG  = {1:2.87, 6:3.15, 7:2.87, 8:2.82, 16:3.70}  # Fdata_HCNOS rcutoffs (Å)

# ── Default output directory ───────────────────────────────────────────────────────
DEFAULT_EXPORT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                                   "tests", "pyFireball", "export", "stm_phase1")

def set_export_dir(path):
    """Set global output directory for plots and npz files."""
    global DEFAULT_EXPORT_DIR
    DEFAULT_EXPORT_DIR = path
    os.makedirs(DEFAULT_EXPORT_DIR, exist_ok=True)

# ── Geometry loading ───────────────────────────────────────────────────────────────
def load_xyz(path):
    """Returns atomTypes (int32 Z), atomPos (float64 [N,3]) from XYZ file."""
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0]); types, pos = [], []
    for ln in lines[2:2+n]:
        p = ln.split(); types.append(ELEM_Z[p[0]]); pos.append([float(p[1]),float(p[2]),float(p[3])])
    return np.array(types,dtype=np.int32), np.array(pos,dtype=np.float64)

def make_mol(name):
    """Built-in molecules (tip models). Returns atomTypes, atomPos."""
    if name=="H2":  return np.array([1,1],dtype=np.int32), np.array([[0,0,-0.37],[0,0,0.37]],dtype=np.float64)
    if name=="CO":  return np.array([6,8],dtype=np.int32), np.array([[0,0,-0.565],[0,0,0.565]],dtype=np.float64)
    raise ValueError(f"Unknown built-in molecule: {name}")

# ── Fireball SCF + DOS computation ────────────────────────────────────────────────
def build_dense_HS(dims, data, atomTypes):
    """
    Build dense H[norb,norb] and S[norb,norb] from Fireball sparse output.
    Uses nzx[ispec]=Z mapping to get norb per atom.
    Returns H, S, orb2atom, norb_per_atom, starts.
    """
    # nzx[ispec]=Z, num_orb[ispec]=norb (ispec 0-based, nspecies entries valid)
    z_to_norb = {int(data.nzx[isp]): int(data.num_orb[isp]) for isp in range(dims.nspecies)}
    norb_per = [z_to_norb[int(data.iatyp[ia])] for ia in range(dims.natoms)]
    starts = np.zeros(dims.natoms+1, dtype=int)
    for ia in range(dims.natoms): starts[ia+1] = starts[ia]+norb_per[ia]
    norb = starts[-1]
    assert norb == dims.norbitals, f"{norb} vs {dims.norbitals}"
    H = np.zeros((norb,norb)); S = np.zeros((norb,norb))
    for ia in range(dims.natoms):
        ni = norb_per[ia]; i0 = starts[ia]
        for k in range(dims.neigh_max):
            j1 = int(data.neigh_j[ia,k])
            if j1<=0: break
            ja=j1-1; nj=norb_per[ja]; j0=starts[ja]
            H[i0:i0+ni,j0:j0+nj] += data.h_mat[ia,k,:ni,:nj]
            S[i0:i0+ni,j0:j0+nj] += data.s_mat[ia,k,:ni,:nj]
    H=0.5*(H+H.T); S=0.5*(S+S.T)
    orb2atom = np.array([ia for ia in range(dims.natoms) for _ in range(norb_per[ia])])
    return H, S, orb2atom, norb_per, starts

def compute_spectral(H, S, eigen_fc, contact_orbs, gamma=0.1, E_range=6.0, dE=0.02):
    """
    PDOS A_mu(E) using generalized eigenbasis H C = S C eps.
    A_mu(E) = sum_n |C_mu_n|^2 * Gamma_n/pi / [(E-eps_n)^2 + Gamma_n^2]
    where Gamma_n = gamma * sum_{mu in contact} |C_mu_n|^2  (WBL coupling).
    Returns E_grid [nE], A [nE,norb], E_fermi.
    """
    from scipy.linalg import eigh
    norb = H.shape[0]
    eps, C = eigh(H, S)          # C[:,n] = eigenvector n in AO basis, S-orthonormal
    n_occ = max(1, norb//2)
    E_F = eps[n_occ-1]           # HOMO as Fermi level (isolated molecule)
    E_LUMO = eps[n_occ] if n_occ<norb else float('nan')
    print(f"  Eigvals (eV): {np.round(eps,4)}")
    print(f"  E_F=HOMO={E_F:.4f}  LUMO={E_LUMO:.4f}  gap={E_LUMO-E_F:.3f} eV")

    E_grid = np.arange(E_F-E_range, E_F+E_range+dE*0.5, dE)
    contact_mask = np.zeros(norb); contact_mask[contact_orbs]=1.0
    Gamma_n = gamma * np.sum(contact_mask[:,None] * C**2, axis=0)
    Gamma_n = np.maximum(Gamma_n, 1e-4)

    A = np.zeros((len(E_grid), norb))
    for n in range(norb):
        lor = (Gamma_n[n]/np.pi) / ((E_grid - eps[n])**2 + Gamma_n[n]**2)
        A += np.outer(lor, C[:,n]**2)
    return E_grid, A, E_F, eps, C

def compute_dos(atomTypes, atomPos, gamma=0.1, contact="all", nmax_scf=200,
                E_range=6.0, dE=0.02, fdata_dir=None, out_path=None):
    """
    Run Fireball SCF, compute spectral function PDOS, save to .npz.
    This is intended to be run as a subprocess to avoid Fireball reinit issues.
    Returns dict with all data (also saves to out_path if given).
    """
    from pyBall import FireCore as fc
    if fdata_dir is None:
        _REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
        fdata_dir = os.path.join(_REPO_ROOT, "tests", "Fireball", "Fdata_HCNOS")

    fc.setVerbosity(0); fc.preinit(); fc.init(atomTypes, atomPos)
    print(f"SCF (nmax={nmax_scf})..."); t0=time.time()
    fc.SCF(atomPos, nmax_scf=nmax_scf); print(f"SCF done {time.time()-t0:.2f}s")
    dims=fc.get_HS_dims(); data=fc.get_HS_neighs(dims); data=fc.get_HS_sparse(dims,data)
    eigen=fc.get_eigen()
    print(f"norb={dims.norbitals}  eigen={np.round(eigen,4)}")

    H, S, orb2atom, norb_per, starts = build_dense_HS(dims, data, atomTypes)
    norb = H.shape[0]
    print(f"H diag: [{np.diag(H).min():.3f},{np.diag(H).max():.3f}]  S diag: [{np.diag(S).min():.3f},{np.diag(S).max():.3f}]")

    if contact=="all":    contact_orbs = np.arange(norb)
    elif contact=="lowest_z": contact_orbs = np.where(orb2atom==int(np.argmin(atomPos[:,2])))[0]
    else: contact_orbs = np.array([int(x) for x in contact.split(",")])

    E_grid, A, E_F, eps, C = compute_spectral(H, S, eigen, contact_orbs, gamma, E_range, dE)
    dE_actual = E_grid[1]-E_grid[0]
    sr = (A*dE_actual).sum(0)
    print(f"Sum-rule (A*dE per orb): mean={sr.mean():.3f} min={sr.min():.3f} max={sr.max():.3f}")

    result = dict(E_grid=E_grid, A=A, E_fermi=E_F, eigen=eigen,
                  atomTypes=atomTypes, atomPos=atomPos, H=H, S=S, orb2atom=orb2atom,
                  norb_per=np.array(norb_per), contact_orbs=contact_orbs, gamma=gamma,
                  eps=eps, C=C)
    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        np.savez(out_path, **result)
        print(f"Saved → {out_path}")
    return result

# ── Hopping (Wolfsberg-Helmholtz) ──────────────────────────────────────────────────
def compute_hopping(H_tip_diag, orb2atom_tip, atomTypes_tip, atomPos_tip,
                    H_smp_diag, orb2atom_smp, atomTypes_smp, atomPos_smp,
                    tip_offset):
    """
    Wolfsberg-Helmholtz hopping model.
    Returns Vsq [nt,ns], R [nt,ns] (inter-orbital distances).
    V_mn = k_WH * (H_mm + H_nn)/2 * exp(-0.5*(z_m+z_n)*r_mn)
    zeta_m = 1/rc_m  (STO decay length ~ orbital radius).
    """
    pos_tip = atomPos_tip + tip_offset
    nt = len(H_tip_diag); ns = len(H_smp_diag)
    zeta_t = np.array([1.0/RC_ANG.get(int(atomTypes_tip[orb2atom_tip[m]]),3.0) for m in range(nt)])
    zeta_s = np.array([1.0/RC_ANG.get(int(atomTypes_smp[orb2atom_smp[n]]),3.0) for n in range(ns)])
    pt = np.array([pos_tip[orb2atom_tip[m]] for m in range(nt)])
    ps = np.array([atomPos_smp[orb2atom_smp[n]] for n in range(ns)])
    diff = pt[:,None,:] - ps[None,:,:]
    R    = np.linalg.norm(diff, axis=2)
    S_inter = np.exp(-0.5*(zeta_t[:,None]+zeta_s[None,:])*R)
    E_avg   = 0.5*(H_tip_diag[:,None] + H_smp_diag[None,:])
    V = 1.75 * E_avg * S_inter
    return V**2, R

# ── STM current ────────────────────────────────────────────────────────────────────
def stm_current(A_tip, A_smp, Vsq, E_grid, E_F, bias=0.1):
    """
    T(E) = einsum(A_tip, Vsq, A_smp); I = integral T dE over [E_F, E_F+bias].
    Returns I (scalar), T (array [nE]).
    """
    T = np.einsum('ei,ij,ej->e', A_tip, Vsq, A_smp)
    mask = (E_grid >= E_F) & (E_grid <= E_F+bias)
    I = float(np.trapz(T[mask], E_grid[mask])) if mask.sum()>1 else 0.0
    return I, T

# ── Plotting utilities ────────────────────────────────────────────────────────────
def plot_dos(ax, E_grid, A, eigen, title, E_F):
    """Plot PDOS grouped by orbital type (s,px,py,pz)."""
    OL = ['s','px','py','pz']; norb = A.shape[1]
    for oi, lab in enumerate(OL):
        cols = np.arange(oi, norb, 4)
        if len(cols)==0: continue
        ax.plot(E_grid, A[:,cols].sum(1), label=lab, lw=1.2)
    ax.plot(E_grid, A.sum(1), 'k--', lw=0.8, label='total')
    for e in eigen: ax.axvline(e, color='gray', lw=0.3, alpha=0.4)
    ax.axvline(E_F, color='red', lw=1.0, ls='--', label='E_F')
    ax.set_xlabel("E (eV)"); ax.set_ylabel("PDOS (1/eV)"); ax.set_title(title); ax.legend(fontsize=7)

def save_plot(fig, name, export_dir=None, dpi=130):
    """Save figure to export_dir (default: DEFAULT_EXPORT_DIR)."""
    if export_dir is None: export_dir = DEFAULT_EXPORT_DIR
    os.makedirs(export_dir, exist_ok=True)
    path = os.path.join(export_dir, f"{name}.png")
    plt.savefig(path, dpi=dpi); print(f"  → {path}")

# ── PDOS real-space projection ─────────────────────────────────────────────────────
def project_pdos_to_grid(A_smp, E_grid, E_F, atomTypes, atomPos, orb2atom,
                         fdata_basis_dir, dE_window=0.3, step=0.2, margin=3.5):
    """
    Project PDOS at E_F ± dE_window onto 3D grid using GridProjector.
    Returns dens3d [nx,ny,nz], grid_spec dict.
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    from pyBall.FireballOCL.FdataParser import FdataParser

    fparser = FdataParser(os.path.dirname(fdata_basis_dir))
    fparser.parse_info()
    z_list = sorted(set(int(z) for z in atomTypes))
    norb_for = {nz: sum(2*l+1 for l in fparser.species_info[nz]['lssh']) for nz in z_list}
    norb_per = [norb_for[int(z)] for z in atomTypes]
    starts   = np.zeros(len(atomTypes)+1,dtype=int)
    for ia,no in enumerate(norb_per): starts[ia+1]=starts[ia]+no
    norb_total = starts[-1]

    mask_e = (E_grid>=E_F-dE_window) & (E_grid<=E_F+dE_window)
    A_win  = A_smp[mask_e].sum(0) * (E_grid[1]-E_grid[0])

    numorb_max = max(norb_per)
    rho = np.zeros((len(atomTypes), 1, numorb_max, numorb_max))
    neigh_j  = np.ones((len(atomTypes), 1), dtype=np.int32)
    for ia in range(len(atomTypes)):
        neigh_j[ia,0] = ia+1
        no = norb_per[ia]; i0 = starts[ia]
        for io in range(no):
            val = A_win[i0+io] if (i0+io)<len(A_win) else 0.0
            rho[ia,0,io,io] = float(val)
    neighn = np.ones(len(atomTypes), dtype=np.int32)

    class Neighs:
        pass
    neighs = Neighs()
    neighs.neigh_j   = neigh_j
    neighs.neighn    = neighn
    neighs.neigh_max = 1
    neighs.iatyp     = atomTypes.copy()
    num_orb_arr = np.zeros(20,dtype=np.int32)
    for nz in z_list: num_orb_arr[nz] = norb_for[nz]
    neighs.num_orb   = num_orb_arr
    neighs.nzx       = np.array(z_list, dtype=np.int32)

    pmin = atomPos.min(0)-margin; pmax = atomPos.max(0)+margin
    ngrid = ((np.ceil((pmax-pmin)/step).astype(int)+7)//8)*8
    grid_spec = {'origin':pmin,'dA':[step,0,0],'dB':[0,step,0],'dC':[0,0,step],'ngrid':ngrid}

    projector = ocl_grid.GridProjector(fdata_basis_dir)
    projector.load_basis(z_list)
    atoms_dict = {'pos':atomPos,'Rcut':np.array([RC_ANG.get(int(z),3.0) for z in atomTypes]),'type':atomTypes}
    dens3d = projector.project(rho, neighs, atoms_dict, grid_spec, nMaxAtom=64)
    return dens3d, grid_spec

def find_homo_from_eps(eps, E_fermi=None):
    """Return HOMO index/energy from eigenvalues.

    Theory/intent:
      - For our STM tests we often want to evaluate NEGF at HOMO energy.
      - We define HOMO as the highest eigenvalue <= E_fermi (if provided).
      - If E_fermi is None, fall back to the maximum eigenvalue (debug use only).
    """
    eps = np.asarray(eps, dtype=float)
    if E_fermi is None:
        i = int(np.argmax(eps))
        return i, float(eps[i])
    m = eps <= float(E_fermi) + 1e-12
    if not np.any(m):
        i = int(np.argmin(np.abs(eps - float(E_fermi))))
        return i, float(eps[i])
    i = int(np.where(m)[0][np.argmax(eps[m])])
    return i, float(eps[i])

def pad_HS_to_float4(atomTypes, atomPos, H, S, norb_per, E_pad=1e3):
    """Pad Fireball AO basis to fixed 4 orbitals per atom: (s,px,py,pz).

    This is a compatibility layer for interfacing with GridFF-like `float4` layout.

    Important numerical note:
      - Pure zero-padding of missing orbitals would make S singular.
      - We therefore set S_ii=1 for padded orbitals, and set H_ii=E_pad (large)
        so padded orbitals are far outside the energy window and effectively decouple.

    Returns:
      H4,S4 : dense matrices of shape (4*natoms, 4*natoms)
      mask  : bool array (natoms,4) True where orbital exists in original basis
      map_o : int array (natoms,4) global AO index in original basis, or -1 if padded
    """
    nat = len(atomTypes)
    H = np.asarray(H, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    norb_per = np.asarray(norb_per, dtype=int)
    starts = _atom_orb_starts(norb_per)
    n4 = 4
    n = nat * n4

    H4 = np.zeros((n, n), dtype=np.float64)
    S4 = np.eye(n, dtype=np.float64)
    mask = np.zeros((nat, n4), dtype=bool)
    map_o = -np.ones((nat, n4), dtype=int)

    for ia in range(nat):
        ni = int(min(norb_per[ia], n4))
        mask[ia, :ni] = True
        i0o = int(starts[ia])
        for io in range(ni):
            map_o[ia, io] = i0o + io
        # set large onsite for padded orbitals (decouple)
        for io in range(ni, n4):
            H4[ia*n4+io, ia*n4+io] = float(E_pad)

    # copy blocks for existing orbitals only
    for ia in range(nat):
        ni = int(min(norb_per[ia], n4))
        i0o = int(starts[ia]); i0 = ia*n4
        for ja in range(nat):
            nj = int(min(norb_per[ja], n4))
            j0o = int(starts[ja]); j0 = ja*n4
            if (ni == 0) or (nj == 0):
                continue
            H4[i0:i0+ni, j0:j0+nj] = H[i0o:i0o+ni, j0o:j0o+nj]
            S4[i0:i0+ni, j0:j0+nj] = S[i0o:i0o+ni, j0o:j0o+nj]

    return H4, S4, mask, map_o

# ── NEGF / Caroli transport on combined system ───────────────────────────────────

_ORT_SPP_TO_STD = np.array([0, 3, 1, 2], dtype=int)  # Ortega (s,py,pz,px) -> (s,px,py,pz)

def _reorder_sp_block_ortega_to_std(b4):
    """Reorder 4x4 block from Ortega order (s,py,pz,px) to (s,px,py,pz)."""
    p = _ORT_SPP_TO_STD
    return b4[np.ix_(p, p)]

def _atom_orb_starts(norb_per_atom):
    starts = np.zeros(len(norb_per_atom) + 1, dtype=int)
    for i, n in enumerate(norb_per_atom):
        starts[i + 1] = starts[i] + int(n)
    return starts

def select_lead_orbitals_by_atoms(orb2atom, atom_ids):
    atom_ids = set(int(i) for i in atom_ids)
    return np.where(np.isin(orb2atom.astype(int), list(atom_ids)))[0].astype(int)

def select_atoms_by_xy_radius(atomPos, rmin=None, rmax=None, center_xy=None):
    if center_xy is None:
        center_xy = atomPos[:, :2].mean(axis=0)
    d = atomPos[:, :2] - center_xy[None, :]
    r = np.sqrt((d * d).sum(axis=1))
    m = np.ones(len(atomPos), dtype=bool)
    if rmin is not None: m &= (r >= float(rmin))
    if rmax is not None: m &= (r <= float(rmax))
    return np.where(m)[0].astype(int)

def build_inter_system_blocks_fdata(atomTypes_T, atomPos_T, norb_per_T, atomTypes_S, atomPos_S, norb_per_S,
                                    fdata_dir, rcut=6.0, root_H='kinetic', root_S='overlap', exp_beta=0.0, exp_r0=0.0, scale_H=1.0, scale_S=1.0):
    """
    Build inter-system coupling blocks using FireballOCL `OCL_Hamiltonian.assemble_2c`.

    Notes:
      - This uses SK-rotated sp-only 4x4 blocks from Fdata tables for given `root`.
      - We slice blocks down to (ni,nj) based on actual Fireball orbital counts (e.g. H has 1 orbital).
      - Optional extra exponential factor: exp(-beta*(r-r0)).

    Returns:
      pairs: (iT,iS) atom index pairs
      Hb: list of H blocks (ni x nj)
      Sb: list of S blocks (ni x nj)
    """
    from pyBall.FireballOCL.OCL_Hamiltonian import OCL_Hamiltonian

    rcut = float(rcut)
    ham = OCL_Hamiltonian(fdata_dir)
    species = sorted(set(int(z) for z in np.concatenate([atomTypes_T, atomTypes_S])))
    ham.prepare_splines(species)

    starts_T = _atom_orb_starts(norb_per_T)
    starts_S = _atom_orb_starts(norb_per_S)

    pairs = []
    dRs = []
    ptH = []
    ptS = []
    for iT, rT in enumerate(atomPos_T):
        for iS, rS in enumerate(atomPos_S):
            dR = (rS - rT).astype(np.float32)
            r = float(np.linalg.norm(dR))
            if r > rcut: continue
            nz1 = int(atomTypes_T[iT])
            nz2 = int(atomTypes_S[iS])
            pth = ham._resolve_pair_type(root_H, nz1, nz2)
            pts = ham._resolve_pair_type(root_S, nz1, nz2)
            if pth is None or pts is None:
                continue
            pairs.append((iT, iS))
            dRs.append(dR)
            ptH.append(int(pth))
            ptS.append(int(pts))

    if len(pairs) == 0:
        return pairs, [], []

    ratoms = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float32)
    neighbors = np.array([[0, 1]], dtype=np.int32)

    Hb = []
    Sb = []
    for (dR, pth, pts, (iT, iS)) in zip(dRs, ptH, ptS, pairs):
        ratoms[1, :] = dR
        bH4 = ham.assemble_2c(ratoms, neighbors, np.array([pth], dtype=np.int32))[0]
        bS4 = ham.assemble_2c(ratoms, neighbors, np.array([pts], dtype=np.int32))[0]
        bH4 = _reorder_sp_block_ortega_to_std(bH4)
        bS4 = _reorder_sp_block_ortega_to_std(bS4)
        ni = int(norb_per_T[iT]); nj = int(norb_per_S[iS])
        bH = bH4[:ni, :nj].astype(np.float64) * float(scale_H)
        bS = bS4[:ni, :nj].astype(np.float64) * float(scale_S)
        if exp_beta != 0.0:
            r = float(np.linalg.norm(dR))
            fac = float(np.exp(-float(exp_beta) * (r - float(exp_r0))))
            bH *= fac
            bS *= fac
        Hb.append(bH)
        Sb.append(bS)
    return pairs, Hb, Sb

def coupling_vec_tip_to_sample(nT, nS, norb_per_T, norb_per_S, pairs, Hb, reduce_tip='sum'):
    """Build effective coupling vector v_smp (len=nS) from inter-system blocks.

    We first assemble the full rectangular H_TS (nT x nS) into a dense buffer, then reduce over tip orbitals.
    This is intended for the "featureless metallic tip" MO-overlap diagnostic.

    reduce_tip:
      - 'sum'  : v = sum_i H_TS[i,:]
      - 'norm' : v = sqrt(sum_i |H_TS[i,:]|^2)
    """
    starts_T = _atom_orb_starts(norb_per_T)
    starts_S = _atom_orb_starts(norb_per_S)
    HTS = np.zeros((int(nT), int(nS)), dtype=np.float64)
    for (iT, iS), bH in zip(pairs, Hb):
        i0 = int(starts_T[iT]); ni = int(norb_per_T[iT])
        j0 = int(starts_S[iS]); nj = int(norb_per_S[iS])
        HTS[i0:i0+ni, j0:j0+nj] += bH
    if reduce_tip == 'sum':
        v = HTS.sum(axis=0)
    elif reduce_tip == 'norm':
        v = np.sqrt((HTS*HTS).sum(axis=0))
    else:
        raise ValueError(f"reduce_tip must be 'sum' or 'norm', got {reduce_tip}")
    return v

def mo_overlap_amplitude(v_smp, C_smp, mo_index):
    """Compute overlap amplitude t = v·C[:,mo] (signed, real) and intensity |t|^2."""
    v = np.asarray(v_smp, dtype=np.float64)
    C = np.asarray(C_smp, dtype=np.float64)
    mo = int(mo_index)
    assert C.shape[0] == v.shape[0]
    t = float(np.dot(v, C[:, mo]))
    return t, t*t


def _blocks_to_dense_vector(pairs, blocks, starts_tip, starts_smp, norb_per_tip, norb_per_smp, ntip, nsmp, out=None):
    """Convert list of (tip_atom, smp_atom) blocks to dense ntip×nsmp matrix (or vector if ntip=1)."""
    if out is None:
        out = np.zeros((ntip, nsmp), dtype=np.float64)
    for (iT, iS), blk in zip(pairs, blocks):
        i0 = int(starts_tip[iT]); ni = int(norb_per_tip[iT])
        j0 = int(starts_smp[iS]); nj = int(norb_per_smp[iS])
        out[i0:i0+ni, j0:j0+nj] += np.asarray(blk, dtype=np.float64)[:ni, :nj]
    return out


def _get_mo_vec(C_s, mo_index, norb):
    C_s = np.asarray(C_s)
    mo = int(mo_index)
    if C_s.ndim != 2:
        raise ValueError(f"C_s must be 2D, got shape={C_s.shape}")
    # Ambiguous case: square matrix (norb,norb). This can appear for small systems.
    # In that case, either rows or columns could be MO vectors depending on the wrapper.
    # Choose the orientation whose selected vector is closer to unit-norm.
    if C_s.shape[0] == norb and C_s.shape[1] == norb and 0 <= mo < norb:
        c_row = np.asarray(C_s[mo, :], dtype=np.complex128)
        c_col = np.asarray(C_s[:, mo], dtype=np.complex128)
        nr = float(np.linalg.norm(c_row))
        nc = float(np.linalg.norm(c_col))
        if abs(nc - 1.0) <= abs(nr - 1.0):
            return c_col
        return c_row
    # fc.get_wfcoef() is typically returned as (nmo,norb) in this repo; tested code
    # (STM_utils.project_orbital_to_points) uses row indexing first.
    if C_s.shape[1] == norb and 0 <= mo < C_s.shape[0]:
        return np.asarray(C_s[mo, :], dtype=np.complex128)
    if C_s.shape[0] == norb and 0 <= mo < C_s.shape[1]:
        return np.asarray(C_s[:, mo], dtype=np.complex128)
    raise ValueError(f"Cannot extract MO vector: C_s.shape={C_s.shape} norb={norb} mo={mo}")


def response_amplitude_map(
    points_tip, H_s, S_s, C_s, mo_index, E,
    tip_types, tip_pos_base, tip_norb_per,
    smp_types, smp_pos, smp_norb_per,
    coupling_builder,
    eta=1e-6, E_tip=0.0
):
    """Compute STM response amplitude map t_resp(r) = |C^T x_s|^2 via block elimination.

    Solves A(E,r) x = b where A = (E+iη)S(r) - H(r), b = [1,0,...]^T on tip.
    Uses precomputed G0 = A_ss^{-1} and block elimination for O(Ns) per point.

    Args:
        points_tip: (npts, 3) tip positions
        H_s, S_s: (ns, ns) sample Hamiltonian and overlap (dense)
        C_s: (ns, nmo) sample MO coefficients
        mo_index: int, which MO to project
        E: float, energy (typically MO eigenvalue)
        tip_types, tip_pos_base, tip_norb_per: tip atom metadata (usually 1 H atom)
        smp_types, smp_pos, smp_norb_per: sample atom metadata
        coupling_builder: callable(tip_pos) -> (pairs, Hb, Sb)
                          e.g. lambda pos: build_inter_system_blocks_exp_sk(..., tip_pos=pos)
        eta: broadening
        E_tip: tip onsite energy

    Returns:
        resp: (npts,) real response amplitudes |C^T x_s|^2
    """
    import numpy as np
    from numpy.linalg import solve

    ns = H_s.shape[0]
    npts = len(points_tip)
    z = float(E) + 1j * float(eta)

    # Precompute sample Green's function G0 = A_ss^{-1}
    A_ss = z * S_s.astype(np.complex128) - H_s.astype(np.complex128)
    G0 = np.linalg.inv(A_ss)

    # MO coefficient vector for this MO
    c = _get_mo_vec(C_s, mo_index, ns)

    # Precompute v = c^T G0 (row vector)
    v = c @ G0  # shape (ns,)

    # Precompute vv = G0^H c*  (for |v·a|^2 = a^H vv vv^H a)
    # Actually we need |v·a|^2 = |sum_i v_i a_i|^2, so we just need v

    starts_tip = _atom_orb_starts(tip_norb_per)
    starts_smp = _atom_orb_starts(smp_norb_per)
    ntip = int(starts_tip[-1])

    resp = np.zeros(npts, dtype=np.float64)

    for ip, r_tip in enumerate(points_tip):
        tip_pos = np.array([r_tip], dtype=np.float64)
        pairs, Hb, Sb = coupling_builder(tip_pos)
        if len(pairs) == 0:
            resp[ip] = 0.0
            continue

        # Build dense H_ts, S_ts (tip×sample)
        H_ts = np.zeros((ntip, ns), dtype=np.float64)
        S_ts = np.zeros((ntip, ns), dtype=np.float64)
        _blocks_to_dense_vector(pairs, Hb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=H_ts)
        _blocks_to_dense_vector(pairs, Sb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=S_ts)

        # a_st = z*S_ts - H_ts  (tip row, sample cols)
        a_st = z * S_ts.astype(np.complex128) - H_ts.astype(np.complex128)  # (ntip, ns)

        # For single tip orbital (ntip=1):
        # x_tip = 1 / (a_tt - a_st G0 a_st^H)
        # x_s = -G0 a_st^H x_tip
        # resp = |c^T x_s|^2 = |x_tip|^2 |v a_st^H|^2

        if ntip == 1:
            a_vec = a_st[0, :]  # (ns,)
            s1 = np.dot(v, a_vec.conj())  # v · a_st^H (complex)
            s2 = np.dot(a_vec, G0 @ a_vec.conj())  # a_st G0 a_st^H (complex)
            a_tt = z * 1.0 - E_tip  # tip has 1 orbital with S=1, H=E_tip
            denom = a_tt - s2
            if abs(denom) < 1e-30:
                resp[ip] = 0.0
                continue
            x_tip = 1.0 / denom
            resp[ip] = (abs(x_tip) ** 2) * (abs(s1) ** 2)
        else:
            # Multi-orbital tip: solve full system (rare)
            A = np.zeros((ntip + ns, ntip + ns), dtype=np.complex128)
            A[:ntip, :ntip] = z * np.eye(ntip) - E_tip * np.eye(ntip)
            A[:ntip, ntip:] = a_st
            A[ntip:, :ntip] = a_st.conj().T
            A[ntip:, ntip:] = A_ss
            b = np.zeros(ntip + ns, dtype=np.complex128)
            b[0] = 1.0
            x = solve(A, b)
            x_s = x[ntip:]
            resp[ip] = abs(np.dot(c, x_s)) ** 2

    return resp


def response_amplitude_simple_lu(
    points_tip, H_s, S_s, C_s, mo_index, E,
    tip_norb_per, smp_norb_per,
    coupling_builder,
    tip_src=None,
    eta=1e-6
):
    from scipy.linalg import lu_factor, lu_solve

    ns = H_s.shape[0]
    npts = len(points_tip)
    z = float(E) + 1j * float(eta)

    starts_tip = _atom_orb_starts(tip_norb_per)
    starts_smp = _atom_orb_starts(smp_norb_per)
    ntip = int(starts_tip[-1])

    if tip_src is None:
        tip_src = np.zeros(ntip, dtype=np.complex128)
        tip_src[0] = 1.0
    else:
        tip_src = np.asarray(tip_src, dtype=np.complex128)
        if tip_src.size != ntip:
            raise ValueError(f"tip_src size mismatch: ntip={ntip} tip_src.size={tip_src.size}")

    A_ss = z * S_s.astype(np.complex128) - H_s.astype(np.complex128)
    c = _get_mo_vec(C_s, mo_index, ns)
    lu, piv = lu_factor(A_ss)

    resp = np.zeros(npts, dtype=np.float64)
    for ip, r_tip in enumerate(points_tip):
        tip_pos = np.array([r_tip], dtype=np.float64)
        pairs, Hb, Sb = coupling_builder(tip_pos)
        if len(pairs) == 0:
            resp[ip] = 0.0
            continue
        H_ts = np.zeros((ntip, ns), dtype=np.float64)
        S_ts = np.zeros((ntip, ns), dtype=np.float64)
        _blocks_to_dense_vector(pairs, Hb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=H_ts)
        _blocks_to_dense_vector(pairs, Sb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=S_ts)
        a_st = z * S_ts.astype(np.complex128) - H_ts.astype(np.complex128)  # (ntip, ns)
        rhs = a_st.conj().T @ tip_src  # (ns,)
        y = lu_solve((lu, piv), rhs)
        amp = np.vdot(c, y)
        resp[ip] = float(abs(amp) ** 2)

    return resp

def _sk_sp_block(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi):
    """Return 4x4 SK block in standard (s,px,py,pz) basis with given direction cosines.

    Row order = tip orbitals, col order = sample orbitals, both in (s,px,py,pz).
    IMPORTANT: caller must reorder columns to match Fireball Ortega order (s,py,pz,px)
    before inserting into H_ts that is multiplied against Fireball's H_s/S_s matrices.
    Use _sk_sp_block_ortega() for that purpose.
    """
    # row: orb on atom1 (tip), col: orb on atom2 (sample), both (s,px,py,pz)
    b = np.zeros((4, 4), dtype=np.float64)
    b[0, 0] = Vss
    b[0, 1] = l * Vsp; b[0, 2] = m * Vsp; b[0, 3] = n * Vsp
    b[1, 0] = l * Vps; b[2, 0] = m * Vps; b[3, 0] = n * Vps
    d = Vpp_sig - Vpp_pi
    b[1, 1] = l*l*d + Vpp_pi; b[2, 2] = m*m*d + Vpp_pi; b[3, 3] = n*n*d + Vpp_pi
    b[1, 2] = l*m*d; b[1, 3] = l*n*d
    b[2, 1] = m*l*d; b[2, 3] = m*n*d
    b[3, 1] = n*l*d; b[3, 2] = n*m*d
    return b


# Fireball/Ortega order: [s, py, pz, px]  (indices 0,1,2,3 in std = s,px,py,pz)
# Ortega col permutation: std col 0->0(s), 1->2(py), 2->3(pz), 3->1(px)
_STD_TO_ORT_COL = np.array([0, 2, 3, 1], dtype=int)  # (s,px,py,pz) -> (s,py,pz,px)
_STD_TO_ORT_ROW = np.array([0, 2, 3, 1], dtype=int)  # same for rows


def _sk_sp_block_ortega(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi):
    """Return 4x4 SK block with rows in (s,py,pz,px) Fireball/Ortega order and
    columns also in Ortega order.

    This matches Fireball's internal H/S matrix convention so that H_ts blocks
    can be directly placed into an augmented H matrix beside H_s (Fireball dense).
    """
    b_std = _sk_sp_block(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi)
    return b_std[np.ix_(_STD_TO_ORT_ROW, _STD_TO_ORT_COL)]

def build_inter_system_blocks_exp_sk(atomTypes_T, atomPos_T, norb_per_T, atomTypes_S, atomPos_S, norb_per_S,
                                     rcut=10.0, beta=1.0, r0=3.0,
                                     A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0,
                                     A_ps=None, overlap_scale=0.0):
    """
    Inter-system tunneling blocks with *vacuum exponential radial* and *SK angular* (sp-only).

    Radial replacement: Fireball/Fdata radial is NOT used.
      V(r) = A * exp(-beta*(r-r0))

    Angular: SK sp formulas for all s/p combinations.
    Orbital order: blocks are produced in Fireball/Ortega order [s,py,pz,px] for both
    tip rows and sample columns, so they can be directly assembled beside Fireball H_s/S_s.

    overlap_scale:
      If >0, S_TS = overlap_scale * H_TS.  If 0, S_TS = 0 (tunneling approximation).
    """
    if A_ps is None: A_ps = A_sp
    rcut = float(rcut); beta = float(beta); r0 = float(r0)
    pairs = []
    Hb = []
    Sb = []
    for iT, rT in enumerate(atomPos_T):
        ni = int(norb_per_T[iT])
        for iS, rS in enumerate(atomPos_S):
            nj = int(norb_per_S[iS])
            dR = (rT - rS).astype(np.float64)
            r = float(np.linalg.norm(dR))
            if r > rcut or r < 1e-8: continue
            l, m, n = (dR / r).tolist()
            f = float(np.exp(-beta * (r - r0)))
            Vss = float(A_ss) * f
            Vsp = float(A_sp) * f
            Vps = float(A_ps) * f
            Vpp_sig = float(A_pp_sig) * f
            Vpp_pi  = float(A_pp_pi)  * f
            # Use Ortega-ordered block so columns match Fireball H_s/S_s orbital basis
            bH4 = _sk_sp_block_ortega(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi)
            if nj == 1:
                bH4[1:4, 0] = 0.0
            bS4 = (overlap_scale * bH4) if overlap_scale != 0.0 else np.zeros((4, 4), dtype=np.float64)
            bH = bH4[:ni, :nj].copy()
            bS = bS4[:ni, :nj].copy()
            pairs.append((iT, iS))
            Hb.append(bH)
            Sb.append(bS)
    return pairs, Hb, Sb

def assemble_combined_HS(H_T, S_T, orb2atom_T, norb_per_T, H_S, S_S, orb2atom_S, norb_per_S,
                         atomPos_T, atomPos_S, pairs, Hb, Sb, tip_offset):
    """Assemble dense combined H,S for tip+sample with provided inter-system blocks."""
    H_T = np.asarray(H_T, dtype=np.float64)
    S_T = np.asarray(S_T, dtype=np.float64)
    H_S = np.asarray(H_S, dtype=np.float64)
    S_S = np.asarray(S_S, dtype=np.float64)
    nT = H_T.shape[0]; nS = H_S.shape[0]
    n = nT + nS
    H = np.zeros((n, n), dtype=np.float64)
    S = np.zeros((n, n), dtype=np.float64)
    H[:nT, :nT] = H_T
    S[:nT, :nT] = S_T
    H[nT:, nT:] = H_S
    S[nT:, nT:] = S_S

    starts_T = _atom_orb_starts(norb_per_T)
    starts_S = _atom_orb_starts(norb_per_S)
    for (iT, iS), bH, bS in zip(pairs, Hb, Sb):
        i0 = int(starts_T[iT]); ni = int(norb_per_T[iT])
        j0 = int(starts_S[iS]); nj = int(norb_per_S[iS])
        # tip->sample blocks (upper-right)
        H[i0:i0+ni, nT+j0:nT+j0+nj] += bH
        S[i0:i0+ni, nT+j0:nT+j0+nj] += bS
        # sample->tip (lower-left) as transpose (real-valued)
        H[nT+j0:nT+j0+nj, i0:i0+ni] += bH.T
        S[nT+j0:nT+j0+nj, i0:i0+ni] += bS.T
    # Symmetrize
    H = 0.5 * (H + H.T)
    S = 0.5 * (S + S.T)
    return H, S

def response_amplitude_full_matrix(
    points_tip, H_s, S_s, C_s, mo_index, E,
    tip_norb_per, smp_norb_per,
    coupling_builder,
    eta=1e-6, E_tip=0.0
):
    """Phase 1: Compute STM response by explicit inversion of full augmented (tip+sample) matrix.

    Builds A_full = (E+iη)S_full - H_full for each tip position and solves A x = b directly.
    Slowest but most explicit; used for validating block elimination approaches.

    A_full = [ a_tt    a_ts^T ]    b = [ 1 ]
             [ a_st    A_ss   ]        [ 0 ]

    resp = |C^T x_s|^2 where x_s = x[n_t:] from the full solve.
    """
    ns = H_s.shape[0]
    npts = len(points_tip)
    z = float(E) + 1j * float(eta)

    starts_tip = _atom_orb_starts(tip_norb_per)
    starts_smp = _atom_orb_starts(smp_norb_per)
    ntip = int(starts_tip[-1])
    n_full = ntip + ns

    A_ss = z * S_s.astype(np.complex128) - H_s.astype(np.complex128)
    c = _get_mo_vec(C_s, mo_index, ns)
    a_tt = np.eye(ntip, dtype=np.complex128) * (z - float(E_tip))

    resp = np.zeros(npts, dtype=np.float64)
    b = np.zeros(n_full, dtype=np.complex128)
    b[0] = 1.0  # excitation on first (only) tip orbital

    for ip, r_tip in enumerate(points_tip):
        tip_pos = np.array([r_tip], dtype=np.float64)
        pairs, Hb, Sb = coupling_builder(tip_pos)
        H_ts = np.zeros((ntip, ns), dtype=np.float64)
        S_ts = np.zeros((ntip, ns), dtype=np.float64)
        _blocks_to_dense_vector(pairs, Hb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=H_ts)
        _blocks_to_dense_vector(pairs, Sb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=S_ts)
        a_st = z * S_ts.astype(np.complex128) - H_ts.astype(np.complex128)  # (ntip, ns)

        A = np.zeros((n_full, n_full), dtype=np.complex128)
        A[:ntip, :ntip] = a_tt
        A[:ntip, ntip:] = a_st
        A[ntip:, :ntip] = a_st.conj().T
        A[ntip:, ntip:] = A_ss

        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            resp[ip] = 0.0
            continue
        x_s = x[ntip:]
        resp[ip] = abs(np.dot(c, x_s)) ** 2

    return resp


def response_amplitude_lu(
    points_tip, H_s, S_s, C_s, mo_index, E,
    tip_norb_per, smp_norb_per,
    coupling_builder,
    eta=1e-6, E_tip=0.0
):
    """Phase 3: Compute STM response via LU precomputation of A_ss (no explicit matrix inverse).

    Precomputes LU factorization of A_ss = (E+iη)S_s - H_s once.
    Per grid point: solves A_ss · y = a_st via forward/back substitution (O(n_s^2)).
    Then applies Schur complement / block elimination.

    Algorithm per tip position:
      y = LU_solve(A_ss, a_st)         # O(n_s^2) — LU already factorized
      s2 = a_st^T · y                   # O(n_s)
      x_tip = 1 / (a_tt - s2)
      resp = |x_tip|^2 · |c^T G0 a_st^H|^2

    Note: we also need v = c^T G0 = c^T LU^{-1}.
    Precompute once: v = LU_solve(A_ss^H, c)^H (complex adjoint solve).
    """
    from scipy.linalg import lu_factor, lu_solve

    ns = H_s.shape[0]
    npts = len(points_tip)
    z = float(E) + 1j * float(eta)

    starts_tip = _atom_orb_starts(tip_norb_per)
    starts_smp = _atom_orb_starts(smp_norb_per)
    ntip = int(starts_tip[-1])

    A_ss = z * S_s.astype(np.complex128) - H_s.astype(np.complex128)
    c = _get_mo_vec(C_s, mo_index, ns)
    a_tt = z * 1.0 - float(E_tip)   # scalar (single tip orbital)

    # Precompute LU factorization once
    lu, piv = lu_factor(A_ss)

    # Precompute v = c^T G0 via adjoint solve: G0^H c* = (A_ss^H)^{-1} c*
    # Then v = conj result conjugated back
    # v_i = sum_j c_j [A_ss^{-1}]_{ji}  =>  v = (A_ss^{-T} c)^*  [for complex symmetric]
    # More directly: v = c @ inv(A_ss) = (inv(A_ss)^H @ c^*)^H
    # Use: A_ss^H y = c  =>  y = (A_ss^H)^{-1} c,  then v = y^H
    v = lu_solve((lu, piv), c.conj(), trans=2).conj()   # v = c^T A_ss^{-1}, shape (ns,)

    resp = np.zeros(npts, dtype=np.float64)

    for ip, r_tip in enumerate(points_tip):
        tip_pos = np.array([r_tip], dtype=np.float64)
        pairs, Hb, Sb = coupling_builder(tip_pos)
        if len(pairs) == 0:
            resp[ip] = 0.0
            continue
        H_ts = np.zeros((ntip, ns), dtype=np.float64)
        S_ts = np.zeros((ntip, ns), dtype=np.float64)
        _blocks_to_dense_vector(pairs, Hb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=H_ts)
        _blocks_to_dense_vector(pairs, Sb, starts_tip, starts_smp, tip_norb_per, smp_norb_per, ntip, ns, out=S_ts)
        a_vec = (z * S_ts.astype(np.complex128) - H_ts.astype(np.complex128))[0, :]  # (ns,)

        # y = A_ss^{-1} · a_vec  (via LU)
        y = lu_solve((lu, piv), a_vec)

        s1 = np.dot(v, a_vec.conj())          # v · a_st^H
        s2 = np.dot(a_vec, y)                  # a_st · G0 · a_st^H
        denom = a_tt - s2
        if abs(denom) < 1e-30:
            resp[ip] = 0.0
            continue
        x_tip = 1.0 / denom
        resp[ip] = abs(x_tip) ** 2 * abs(s1) ** 2

    return resp


# ──────────────────────────────────────────────────────────────────────────
#  Rigid-body rotation utilities  (STM with moving molecules)
# ──────────────────────────────────────────────────────────────────────────

def build_atom_rotation_matrix(norb, R):
    """Build orbital rotation matrix U for a single atom under 3×3 rotation R.

    For s+p basis (norb = 1 or 4):
      - s orbital is invariant   → U[0,0] = 1
      - p orbitals rotate as vectors → U[1:4, 1:4] = R_ortega

    IMPORTANT: Uses Fireball / Ortega orbital ordering [s, py, pz, px].
    The standard Cartesian p-ordering [px, py, pz] must be permuted:
        Ortega [py, pz, px] = standard [py, pz, px] at indices [1,2,0].

    If norb == 1 (H atom) only s is present.

    Args:
        norb: int, number of orbitals on this atom (1 or 4)
        R:    (3,3) float rotation matrix (orthogonal, det=+1)
    Returns:
        U: (norb, norb) float rotation matrix
    """
    norb = int(norb)
    if norb == 1:
        return np.array([[1.0]], dtype=np.float64)
    elif norb == 4:
        # Fireball uses real p orbitals in Ortega order: [py, pz, px].
        # Rotation of p orbitals must match Fortran ROTATIONS/twister.f90:
        #   pmat(1,1)=eps(2,2); pmat(1,2)=eps(2,3); pmat(1,3)=eps(2,1)
        #   pmat(2,1)=eps(3,2); pmat(2,2)=eps(3,3); pmat(2,3)=eps(3,1)
        #   pmat(3,1)=eps(1,2); pmat(3,2)=eps(1,3); pmat(3,3)=eps(1,1)
        # Here we take eps := R for rigid-body rotation.
        eps = np.asarray(R, dtype=np.float64)
        pmat = np.array([
            [eps[1, 1], eps[1, 2], eps[1, 0]],
            [eps[2, 1], eps[2, 2], eps[2, 0]],
            [eps[0, 1], eps[0, 2], eps[0, 0]],
        ], dtype=np.float64)
        U = np.eye(4, dtype=np.float64)
        U[1:4, 1:4] = pmat
        return U
    else:
        raise ValueError(f"build_atom_rotation_matrix: unsupported norb={norb} (only 1 or 4)")


def rotate_dense_hs(H, S, norb_per_atom, R):
    """Apply a global rigid-body rotation R to dense H and S matrices.

    For each atom i with orbital rotation U_i = build_atom_rotation_matrix(norb_i, R):
        H' = U_full @ H @ U_full^T
        S' = U_full @ S @ U_full^T
    where U_full = block_diag(U_0, U_1, ..., U_{N-1}).

    This is the *global* rotation (molecule rotates as a rigid body).
    It is distinct from the bond-direction Slater–Koster rotation used
    during Hamiltonian assembly.

    Args:
        H, S: (norb, norb) dense Hamiltonian and overlap
        norb_per_atom: (natoms,) int, orbitals per atom
        R: (3,3) rotation matrix
    Returns:
        H_rot, S_rot: rotated dense matrices (same dtype)
    """
    H = np.asarray(H)
    S = np.asarray(S)
    natoms = len(norb_per_atom)
    norb = H.shape[0]

    # Build full U = block_diag(U_0, U_1, ...)
    U = np.zeros((norb, norb), dtype=np.float64)
    i0 = 0
    for ia in range(natoms):
        ni = int(norb_per_atom[ia])
        Ui = build_atom_rotation_matrix(ni, R)
        U[i0:i0+ni, i0:i0+ni] = Ui
        i0 += ni
    if i0 != norb:
        raise ValueError(f"rotate_dense_hs: orbital count mismatch {i0} != {norb}")

    # Keep complex dtype if input is complex
    Uc = U.astype(np.complex128) if (np.iscomplexobj(H) or np.iscomplexobj(S)) else U
    # Active rotation acting on orbital basis vectors: H' = U H U^T
    H_rot = Uc @ H @ Uc.T
    S_rot = Uc @ S @ Uc.T
    return H_rot, S_rot


def rotate_mo_coefficients(C, norb_per_atom, R):
    """Rotate MO coefficients C (norb, nmo) by global rotation R.

    Each column of C is an MO vector in the AO basis.
    Under rotation of the molecule, the AO basis rotates, so
    the coefficients transform as  C' = U @ C.

    Args:
        C: (norb, nmo) MO coefficients
        norb_per_atom: (natoms,) int
        R: (3,3) rotation matrix
    Returns:
        C_rot: (norb, nmo) rotated coefficients
    """
    C = np.asarray(C)
    natoms = len(norb_per_atom)
    norb = C.shape[0]

    U = np.zeros((norb, norb), dtype=np.float64)
    i0 = 0
    for ia in range(natoms):
        ni = int(norb_per_atom[ia])
        Ui = build_atom_rotation_matrix(ni, R)
        U[i0:i0+ni, i0:i0+ni] = Ui
        i0 += ni

    Uc = U.astype(np.complex128) if np.iscomplexobj(C) else U
    return Uc @ C


def rotate_points(points, R, center=None):
    """Rotate 3D points by rotation matrix R around optional center.

    Args:
        points: (npts, 3) array of positions
        R:      (3,3) rotation matrix
        center: (3,) center of rotation (default = centroid of points)
    Returns:
        points_rot: (npts, 3) rotated positions
    """
    points = np.asarray(points, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)
    if center is None:
        center = points.mean(axis=0)
    else:
        center = np.asarray(center, dtype=np.float64)
    return (points - center) @ R.T + center


def build_rotated_plane_grid(mol_pos, z, size=20.0, n=40, R=None):
    """Build a 2D imaging grid in a plane that may be rotated.

    Original (unrotated) grid lies in the xy-plane at height z above
    the molecule centroid.  After rotation R the grid normal is
    R @ [0,0,1] and the in-plane x/y axes are R @ [1,0,0] and R @ [0,1,0].

    This means: if the molecule is rotated by R, we rotate the imaging
    plane by the *same* R so that the STM image is rotationally invariant.

    Args:
        mol_pos: (natoms, 3) atomic positions
        z:       float, perpendicular distance from centroid (Å)
        size:    float, lateral extent of grid (Å)
        n:       int, grid resolution (n×n)
        R:       (3,3) rotation matrix or None (no rotation)
    Returns:
        X, Y:    (n,n) meshgrid coordinates (in rotated frame, for plotting)
        points:  (n*n, 3) 3D points in the rotated plane
        extent:  [xmin, xmax, ymin, ymax] for matplotlib
    """
    mol_pos = np.asarray(mol_pos, dtype=np.float64)
    origin = mol_pos.mean(axis=0)

    if R is None:
        R = np.eye(3, dtype=np.float64)
    else:
        R = np.asarray(R, dtype=np.float64)

    # Match STM_utils.build_plane_grid() convention exactly when R=I
    xs = np.linspace(origin[0] - size * 0.5, origin[0] + size * 0.5, int(n))
    ys = np.linspace(origin[1] - size * 0.5, origin[1] + size * 0.5, int(n))
    X, Y = np.meshgrid(xs, ys, indexing='ij')

    # Active rotation for row-vector points: p_rot = p @ R^T
    # Plane basis vectors are columns of R: ex=R[:,0], ey=R[:,1], ez=R[:,2].
    ex = R[:, 0]
    ey = R[:, 1]
    ez = R[:, 2]

    points = np.zeros((int(n)*int(n), 3), dtype=np.float64)
    for i in range(int(n)):
        for j in range(int(n)):
            idx = i * int(n) + j
            dx = xs[i] - origin[0]
            dy = ys[j] - origin[1]
            points[idx] = origin + dx * ex + dy * ey + float(z) * ez

    extent = [xs[0], xs[-1], ys[0], ys[-1]]

    return X, Y, points, extent


def wbl_self_energy(n, lead_orbs, gamma):
    """Wide-band limit self-energy Σ = -i Γ/2 on selected orbitals."""
    lead_orbs = np.array(lead_orbs, dtype=int)
    Sig = np.zeros((n, n), dtype=np.complex128)
    if lead_orbs.size == 0: return Sig
    g = float(gamma)
    Sig[lead_orbs, lead_orbs] = -0.5j * g
    return Sig

def caroli_T(E, H, S, SigmaL, SigmaR, eta=1e-6):
    """Compute transmission T(E) by Caroli formula for dense matrices."""
    n = H.shape[0]
    z = (E + 1j * float(eta))
    A = z * S.astype(np.complex128) - H.astype(np.complex128) - SigmaL - SigmaR
    # Full inverse via solve(A, I)
    G = solve(A, np.eye(n, dtype=np.complex128))
    # Γ = i(Σ - Σ†)
    GammaL = 1j * (SigmaL - SigmaL.conj().T)
    GammaR = 1j * (SigmaR - SigmaR.conj().T)
    T = np.real(np.trace(GammaL @ G @ GammaR @ G.conj().T))
    return float(T)

def _make_A_linop(E, H, S, SigmaL, SigmaR, eta=1e-6, herm=False):
    """Build scipy LinearOperator for A(E) or A(E)^H if herm=True."""
    from scipy.sparse.linalg import LinearOperator
    n = H.shape[0]
    z = (float(E) + 1j * float(eta))
    A = z * S.astype(np.complex128) - H.astype(np.complex128) - SigmaL - SigmaR
    if herm:
        AH = A.conj().T
        def mv(x):
            return AH @ x
        return LinearOperator((n, n), matvec=mv, dtype=np.complex128), AH
    else:
        def mv(x):
            return A @ x
        return LinearOperator((n, n), matvec=mv, dtype=np.complex128), A

def caroli_T_iterative(E, H, S, SigmaL, SigmaR, eta=1e-6, nrand=8, tol=1e-6, maxiter=500, x0_cache=None, seed=0):
    """
    Estimate T(E)=Tr[GammaL G GammaR G†] without forming G.

    Uses Hutchinson estimator with Rademacher random vectors z:
      T = E_z [ z† GammaL G GammaR G† z ]

    Each sample requires two linear solves:
      A† u = z
      A  x = GammaR u
    then contribution = (x† GammaL x)

    x0_cache: optional dict used to warm-start GMRES across energies/scan points.
              Keys: ('u',E,i), ('x',E,i) if you want; we keep it lightweight here.
    """
    from scipy.sparse.linalg import gmres
    n = H.shape[0]
    Aop, _ = _make_A_linop(E, H, S, SigmaL, SigmaR, eta=eta, herm=False)
    AHop, _ = _make_A_linop(E, H, S, SigmaL, SigmaR, eta=eta, herm=True)
    GammaL = 1j * (SigmaL - SigmaL.conj().T)
    GammaR = 1j * (SigmaR - SigmaR.conj().T)

    rng = np.random.default_rng(int(seed))
    Tacc = 0.0
    for ir in range(int(nrand)):
        z = rng.choice(np.array([-1.0, 1.0], dtype=np.float64), size=n).astype(np.complex128)

        # Solve A† u = z
        u0 = None
        if x0_cache is not None:
            u0 = x0_cache.get(('u', float(E), ir), None)
        u, info = gmres(AHop, z, x0=u0, tol=float(tol), maxiter=int(maxiter))
        if info != 0:
            raise RuntimeError(f"GMRES failed for A^H u=z at E={E} ir={ir} info={info}")

        # Solve A x = GammaR u
        b = GammaR @ u
        x0 = None
        if x0_cache is not None:
            x0 = x0_cache.get(('x', float(E), ir), None)
        x, info = gmres(Aop, b, x0=x0, tol=float(tol), maxiter=int(maxiter))
        if info != 0:
            raise RuntimeError(f"GMRES failed for A x=b at E={E} ir={ir} info={info}")

        if x0_cache is not None:
            x0_cache[('u', float(E), ir)] = u
            x0_cache[('x', float(E), ir)] = x

        Tacc += float(np.real(np.vdot(x, GammaL @ x)))
    return Tacc / float(nrand)

def negf_current(E_grid, H, S, leadL_orbs, leadR_orbs, gammaL=1.0, gammaR=1.0, eta=1e-6, E0=None, E1=None):
    """Compute T(E) on grid and integrate current over [E0,E1] (defaults full range)."""
    n = H.shape[0]
    SigmaL = wbl_self_energy(n, leadL_orbs, gammaL)
    SigmaR = wbl_self_energy(n, leadR_orbs, gammaR)
    T = np.array([caroli_T(float(E), H, S, SigmaL, SigmaR, eta=eta) for E in E_grid], dtype=np.float64)
    if E0 is None: E0 = float(E_grid[0])
    if E1 is None: E1 = float(E_grid[-1])
    m = (E_grid >= float(E0)) & (E_grid <= float(E1))
    I = float(np.trapz(T[m], E_grid[m])) if m.sum() > 1 else 0.0
    return I, T

def negf_response_at_energy(E, H, S, leadL_orbs, leadR_orbs, gammaL=1.0, gammaR=1.0, eta=1e-6,
                            source='L', measure='R', normalize=True):
    """
    Deterministic single-solve "response" metric between leads.

    Construct a source vector b on one lead and measure response amplitude on the other lead:
      A(E) x = b
      resp = Re[ x† Γ_meas x ]

    This is *not* the exact Caroli trace, but is a practical probe-like metric with no stochastic sampling.
    """
    n = H.shape[0]
    SigmaL = wbl_self_energy(n, leadL_orbs, gammaL)
    SigmaR = wbl_self_energy(n, leadR_orbs, gammaR)
    GammaL = 1j * (SigmaL - SigmaL.conj().T)
    GammaR = 1j * (SigmaR - SigmaR.conj().T)

    # A(E)
    z = (float(E) + 1j * float(eta))
    A = z * S.astype(np.complex128) - H.astype(np.complex128) - SigmaL - SigmaR

    # Source vector (inject)
    b = np.zeros(n, dtype=np.complex128)
    if source.upper() == 'L':
        idx = np.array(leadL_orbs, dtype=int)
        g = np.sqrt(np.maximum(np.real(np.diag(GammaL)[idx]), 0.0))
        b[idx] = g
    elif source.upper() == 'R':
        idx = np.array(leadR_orbs, dtype=int)
        g = np.sqrt(np.maximum(np.real(np.diag(GammaR)[idx]), 0.0))
        b[idx] = g
    else:
        raise ValueError(f"source must be 'L' or 'R', got {source}")

    # Measurement projector
    if measure.upper() == 'L':
        Gmeas = GammaL
    elif measure.upper() == 'R':
        Gmeas = GammaR
    elif measure.lower() == 'mol':
        # "molecular susceptibility" proxy: measure total |x|^2 on the molecule subspace
        # (here: all orbitals) -> identity
        Gmeas = np.eye(n, dtype=np.complex128)
    else:
        raise ValueError(f"measure must be 'L','R','mol', got {measure}")

    if normalize:
        nb = float(np.linalg.norm(b))
        assert nb > 0.0
        b /= nb

    x = solve(A, b)
    resp = float(np.real(np.vdot(x, Gmeas @ x)))
    return resp

def negf_current_iterative(E_grid, H, S, leadL_orbs, leadR_orbs, gammaL=1.0, gammaR=1.0, eta=1e-6,
                           E0=None, E1=None, nrand=8, tol=1e-6, maxiter=500, seed=0, reuse=True):
    """Iterative version of negf_current() using caroli_T_iterative()."""
    n = H.shape[0]
    SigmaL = wbl_self_energy(n, leadL_orbs, gammaL)
    SigmaR = wbl_self_energy(n, leadR_orbs, gammaR)
    cache = {} if reuse else None
    T = np.array([caroli_T_iterative(float(E), H, S, SigmaL, SigmaR, eta=eta, nrand=nrand, tol=tol, maxiter=maxiter, x0_cache=cache, seed=seed)
                  for E in E_grid], dtype=np.float64)
    if E0 is None: E0 = float(E_grid[0])
    if E1 is None: E1 = float(E_grid[-1])
    m = (E_grid >= float(E0)) & (E_grid <= float(E1))
    I = float(np.trapz(T[m], E_grid[m])) if m.sum() > 1 else 0.0
    return I, T
