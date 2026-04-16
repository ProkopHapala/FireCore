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

def _sk_sp_block(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi):
    """Return 4x4 SK block for (s,px,py,pz) basis with given direction cosines."""
    # row: orb on atom1 (tip), col: orb on atom2 (sample)
    b = np.zeros((4, 4), dtype=np.float64)
    b[0, 0] = Vss
    # s-p (atom2)
    b[0, 1] = l * Vsp
    b[0, 2] = m * Vsp
    b[0, 3] = n * Vsp
    # p-s (atom1)
    b[1, 0] = l * Vps
    b[2, 0] = m * Vps
    b[3, 0] = n * Vps
    # p-p
    d = Vpp_sig - Vpp_pi
    b[1, 1] = l*l*d + Vpp_pi
    b[2, 2] = m*m*d + Vpp_pi
    b[3, 3] = n*n*d + Vpp_pi
    b[1, 2] = l*m*d
    b[1, 3] = l*n*d
    b[2, 1] = m*l*d
    b[2, 3] = m*n*d
    b[3, 1] = n*l*d
    b[3, 2] = n*m*d
    return b

def build_inter_system_blocks_exp_sk(atomTypes_T, atomPos_T, norb_per_T, atomTypes_S, atomPos_S, norb_per_S,
                                     rcut=10.0, beta=1.0, r0=3.0,
                                     A_ss=-1.0, A_sp=-1.0, A_pp_sig=-1.0, A_pp_pi=+1.0,
                                     A_ps=None, overlap_scale=0.0):
    """
    Inter-system tunneling blocks with *vacuum exponential radial* and *SK angular* (sp-only).

    Radial replacement (your requirement): Fireball/Fdata radial is NOT used.
      V(r) = A * exp(-beta*(r-r0))

    Angular:
      SK sp formulas give all s/px/py/pz combinations.

    overlap_scale:
      If >0, we set S_TS = overlap_scale * exp(-beta*(r-r0)) * angular_block with same channel amps.
      If 0, S_TS = 0 (common tunneling approximation).
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
            dR = (rS - rT).astype(np.float64)
            r = float(np.linalg.norm(dR))
            if r > rcut or r < 1e-8: continue
            l, m, n = (dR / r).tolist()
            f = float(np.exp(-beta * (r - r0)))
            Vss = float(A_ss) * f
            Vsp = float(A_sp) * f
            Vps = float(A_ps) * f
            Vpp_sig = float(A_pp_sig) * f
            Vpp_pi  = float(A_pp_pi)  * f
            bH4 = _sk_sp_block(l, m, n, Vss, Vsp, Vps, Vpp_sig, Vpp_pi)
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
