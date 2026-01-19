import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, Tuple, List

@dataclass
class ComparisonResult:
    """Simple result container for comparison operations.
    
    Replaces complex dict returns with a clean interface.
    Detailed data can be accessed via .details dict if needed.
    """
    ok: bool
    err: float
    details: Dict[str, Any] = field(default_factory=dict)
    
    def __bool__(self) -> bool:
        return self.ok
    
    def __repr__(self) -> str:
        return f"ComparisonResult(ok={self.ok}, err={self.err:.3e})"

import pyopencl as cl
import os
from .FdataParser import FdataParser
from .OCL_Hamiltonian import OCL_Hamiltonian
from pyBall import FireCore as fc

## =============================================================
##   Utilities for testing and checking wrt Fotran reference
## ==============================================================

# NOTE: these helpers are intentionally module-level so they can be imported
#       from testing scripts like tests/pyFireball/verify_C3.py without having to
#       instantiate OCL_Hamiltonian. Keep them lightweight and pure numpy.

def compare_matrices(name, fortran_mat, ocl_mat, tol=1e-5, require_nonzero=False):
    print(f"\n--- Comparing {name} ---")
    diff = np.abs(fortran_mat - ocl_mat)
    max_diff = np.max(diff) if diff.size else 0.0
    max_val = max(np.max(np.abs(fortran_mat)) if fortran_mat.size else 0.0, np.max(np.abs(ocl_mat)) if ocl_mat.size else 0.0)
    print(f"Max difference: {max_diff:.2e}")

    if require_nonzero and max_val < tol:
        print(f"WARNING: {name} is (near) zero for both; treating as failure.")
        return False, max_diff
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        print("Fortran Matrix:")
        print(fortran_mat)
        print("PyOpenCL Matrix:")
        print(ocl_mat)
        print("Abs Diff:")
        print(diff)
        return False, max_diff
    print(f"SUCCESS: {name} matches.")
    return True, max_diff


def compare_matrices_brief(name, A, B, tol=1e-5, require_nonzero=False):
    diff = np.abs(A - B)
    max_diff = float(np.max(diff)) if diff.size else 0.0
    max_val = max(float(np.max(np.abs(A))) if A.size else 0.0, float(np.max(np.abs(B))) if B.size else 0.0)
    print(f"\n--- Comparing {name} (brief) ---")
    print(f"Max difference: {max_diff:.2e}")
    if diff.size:
        ij = np.unravel_index(int(np.argmax(diff)), diff.shape)
        print(f"Argmax index: {ij}, A={float(A[ij]):.6g}, B={float(B[ij]):.6g}")
    if require_nonzero and max_val < tol:
        print(f"WARNING: {name} is (near) zero for both; treating as failure.")
        return False, max_diff
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        return False, max_diff
    print(f"SUCCESS: {name} matches.")
    return True, max_diff


def _orbital_layout(sparse_data, natoms):
    nzx = np.array(sparse_data.nzx, dtype=np.int32)
    num_orb_species = np.array(sparse_data.num_orb, dtype=np.int32)
    iatyp = np.array(sparse_data.iatyp, dtype=np.int32)
    n_orb_atom = np.zeros(natoms, dtype=np.int32)
    for ia in range(natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec = int(w[0])
        n_orb_atom[ia] = int(num_orb_species[ispec])
    if np.any(n_orb_atom == 0):
        raise RuntimeError(f"Zero-orbital atom encountered: n_orb_atom={n_orb_atom}, iatyp={iatyp}, num_orb_species={num_orb_species}, nzx={nzx}")
    offs = np.zeros(natoms + 1, dtype=np.int32)
    offs[1:] = np.cumsum(n_orb_atom)
    return n_orb_atom, offs


def _blocked_to_dense(sparse_data, H_blocks, natoms):
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    norb = int(offs[-1])
    M = np.zeros((norb, norb), dtype=np.float64)
    neighn = np.array(sparse_data.neighn, dtype=np.int32)
    neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
    neigh_b = np.array(sparse_data.neigh_b, dtype=np.int32)

    # Fortran uses neigh_self(iatom) = ineigh index where (jatom==iatom && mbeta==0).
    # Do NOT assume it is the last neighbor slot.
    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1  # Fortran 1-based atom index stored in neigh_j
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break

    for i in range(natoms):
        ni = int(n_orb_atom[i])
        i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                j = i
            else:
                j = int(neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            nj = int(n_orb_atom[j])
            j0 = int(offs[j])
            blk = H_blocks[i, ineigh, :nj, :ni]
            M[i0:i0+ni, j0:j0+nj] += blk.T
    return M


def _blocked_to_dense_vca_atom(sparse_data, A_blocks, natoms):
    # VCA atom term is a diagonal update on atom i accumulated over non-self neighbors.
    # For the core debug export we therefore interpret A_blocks[i,ineigh] as contributing to (i,i),
    # NOT to (i,j). This matches OpenCL assemble_vca which returns diag updates separately.
    n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
    norb = int(offs[-1])
    M = np.zeros((norb, norb), dtype=np.float64)
    neighn = np.array(sparse_data.neighn, dtype=np.int32)
    neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
    neigh_b = np.array(sparse_data.neigh_b, dtype=np.int32)

    # Detect whether A_blocks is stored predominantly in the self-slot or in neighbor slots.
    # We must not guess from interface functions; infer from exported array content only.
    max_self = 0.0
    max_nons = 0.0
    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(neighn[i])):
            if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break
        if neigh_self[i] >= 0:
            ni = int(n_orb_atom[i])
            ms = float(np.max(np.abs(A_blocks[i, neigh_self[i], :ni, :ni])))
            if ms > max_self:
                max_self = ms
        ni = int(n_orb_atom[i])
        for ineigh in range(int(neighn[i])):
            if ineigh == int(neigh_self[i]):
                continue
            mn = float(np.max(np.abs(A_blocks[i, ineigh, :ni, :ni])))
            if mn > max_nons:
                max_nons = mn

    use_self_slot = (max_self > (10.0 * max_nons))
    print(f"  [VCA_ATOM export layout] max_self={max_self:.3e} max_nonself={max_nons:.3e} => use_self_slot={int(use_self_slot)}")

    for i in range(natoms):
        ni = int(n_orb_atom[i])
        i0 = int(offs[i])
        if use_self_slot and neigh_self[i] >= 0:
            blk = A_blocks[i, neigh_self[i], :ni, :ni]
            M[i0:i0+ni, i0:i0+ni] += blk.T
        else:
            ii = i + 1
            for ineigh in range(int(neighn[i])):
                if int(neigh_j[i, ineigh]) == ii and int(neigh_b[i, ineigh]) == 0:
                    continue
                j = int(neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                blk = A_blocks[i, ineigh, :ni, :ni]
                M[i0:i0+ni, i0:i0+ni] += blk.T
    return M


def compare_blocks(name, A, B, tol=1e-5, require_nonzero=False):
    if not np.all(np.isfinite(A)):
        raise RuntimeError(f"compare_blocks: non-finite values in A for {name}")
    if not np.all(np.isfinite(B)):
        raise RuntimeError(f"compare_blocks: non-finite values in B for {name}")
    diff = np.abs(A - B)
    max_diff = float(np.max(diff)) if diff.size else 0.0
    max_val = max(float(np.max(np.abs(A))) if A.size else 0.0, float(np.max(np.abs(B))) if B.size else 0.0)
    print(f"\n--- Comparing {name} ---")
    print(f"Max difference: {max_diff:.2e}")
    if require_nonzero and max_val < tol:
        print(f"WARNING: {name} is (near) zero for both; treating as failure.")
        return False, max_diff
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        return False, max_diff
    print(f"SUCCESS: {name} matches.")
    return True, max_diff


def dense_from_pair_blocks(B00, B01, B10, B11):
    """
    Assemble a dense matrix from four 4x4 blocks representing a 2-atom system.
    Matches the layout used in verify_C2.py for C2 dimer tests.
    """
    M = np.zeros((8, 8), dtype=np.float64)
    M[0:4, 0:4] = B00
    M[0:4, 4:8] = B01
    M[4:8, 0:4] = B10
    M[4:8, 4:8] = B11
    return M


def dense_from_neighbor_list(neighs, blocks, n_orb_atom, offs):
    """
    Convert neighbor list + blocks into dense matrix.
    neighs: list of (i, j) pairs
    blocks: array of shape (len(neighs), nj, ni) in (nu, mu) order
    n_orb_atom: per-atom orbital counts
    offs: orbital offsets (cumulative sum of n_orb_atom)
    """
    norb = int(offs[-1])
    M = np.zeros((norb, norb), dtype=np.float64)
    for idx, (i, j) in enumerate(neighs):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        # blocks are (nu,mu); dense expects (mu,nu)
        M[i0:i0+ni, j0:j0+nj] = blocks[idx, :nj, :ni].T
    return M


def scan2c_fortran(fc, interaction, dR, in3=None, isub=0, applyRotation=True):
    """
    Thin wrapper around FireCore.scanHamPiece2c for 2-center interactions.
    Mirrors usage in verify_C2.py and verify_C3.py.
    Reference: pyBall/FireCore.py scanHamPiece2c
    """
    in1 = 1
    in2 = 1
    if in3 is None:
        in3 = in2
    return fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation=applyRotation)


def scan2c_ocl(ham, root, nz1, nz2, dR, applyRotation=True):
    """
    Thin wrapper around OCL_Hamiltonian.scanHamPiece2c for 2-center interactions.
    Mirrors usage in verify_C2.py and verify_C3.py.
    """
    return ham.scanHamPiece2c(root, int(nz1), int(nz2), dR, applyRotation=applyRotation)


def firecore_sparse_to_dense(fc, export_mode=0, natoms=None):
    """
    Retrieve full H and S matrices from Fortran FireCore as dense matrices.
    Automates fc.set_options/export_mode + fc.get_HS_* calls.
    Reference: libFireCore.f90 exports (h_mat, s_mat, neighn, neigh_j, etc.)
    """
    fc.set_options(1, 1, 1, 1, 1, 1, 1)
    fc.set_export_mode(export_mode)
    dims_ = fc.get_HS_dims()
    sd_ = fc.get_HS_neighs(dims_)
    sd_ = fc.get_HS_sparse(dims_, sd_)
    if natoms is None:
        natoms = int(dims_.natoms)
    H = _blocked_to_dense(sd_, sd_.h_mat, natoms)
    S = _blocked_to_dense(sd_, sd_.s_mat, natoms)
    neighbors = []
    for i in range(natoms):
        for ineigh in range(int(sd_.neighn[i])):
            j = int(sd_.neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            neighbors.append((i, j))
    return H, S, neighbors, sd_


def firecore_sparse_H_with_options(fc, ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, ioff_Ewald=1, export_mode=1, natoms=None):
    """
    Retrieve H matrix from Fortran FireCore with specific component toggles.
    Reference: libFireCore.f90 set_options / get_HS_sparse
    """
    fc.set_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, ioff_Ewald)
    fc.set_export_mode(export_mode)
    dims_ = fc.get_HS_dims()
    sd_ = fc.get_HS_neighs(dims_)
    sd_ = fc.get_HS_sparse(dims_, sd_)
    if natoms is None:
        natoms = int(dims_.natoms)
    H = _blocked_to_dense(sd_, sd_.h_mat, natoms)
    neighbors = []
    for i in range(natoms):
        for ineigh in range(int(sd_.neighn[i])):
            j = int(sd_.neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            neighbors.append((i, j))
    return H, neighbors, sd_


def cl_sp_from_ham(ham, atomTypes_Z, atom_idx):
    """
    Extract sp-projector CL vector for an atom from OCL_Hamiltonian parser cache.
    Enforces sp-only (len=4) layout for VNL contraction tests.
    Reference: assemble_VNL.f90 projector construction
    """
    Z = int(atomTypes_Z[atom_idx])
    cl_shell = ham.cl_pseudo.get(Z, None)
    if cl_shell is None:
        raise RuntimeError(f"Missing ham.cl_pseudo for Z={Z}")
    info = ham.parser.species_info.get(Z, None)
    if info is None:
        raise RuntimeError(f"Missing species_info for Z={Z}")
    lsshPP = info.get('lsshPP', [])
    cl_full = []
    for ish, l in enumerate(lsshPP):
        c = float(cl_shell[ish])
        cl_full += [c] * (2 * int(l) + 1)
    cl_full = np.array(cl_full, dtype=np.float64)
    if cl_full.shape[0] != 4:
        raise RuntimeError(f"Expected sp projectors (len=4), got len={cl_full.shape[0]} for Z={Z}")
    return cl_full


def contract_vnl_blocks(A, B, clv):
    """
    Contract VNL blocks with projector coefficients.
    A,B: (4,4) blocks (nu,mu) => nu=cc, mu=basis index
    clv: projector coefficients (len=4)
    Output: (nu,mu)
    Reference: OpenCL CPU contraction / assemble_vnl Fortran
    """
    out = np.zeros((4, 4), dtype=np.float64)
    for nu in range(4):
        for mu in range(4):
            v = 0.0
            for cc in range(4):
                v += A[cc, mu] * clv[cc] * B[cc, nu]
            out[nu, mu] = v
    return out


def compare_ewaldsr(ham, atomPos, atomTypes_Z, sd, dims, neigh_self, neigh_back, dq_atom, ewaldsr4_f, ewaldsr3c4_f, tol=1e-5, verbose=False):
    res = ham.assemble_ewaldsr_blocks(atomPos, atomTypes_Z, sd, neigh_self, neigh_back, dq_atom, ewaldsr3c_mask=ewaldsr3c4_f)
    ewaldsr4_o = res['ewaldsr4']
    EwaldSR_o = _blocked_to_dense(sd, ewaldsr4_o, int(dims.natoms))
    EwaldSR_f = _blocked_to_dense(sd, ewaldsr4_f, int(dims.natoms))
    ok, err = compare_matrices_brief("EwaldSR (OCL helper)", EwaldSR_f, EwaldSR_o, tol=tol, require_nonzero=True)
    
    if verbose:
        print(f"  EwaldSR comparison: ok={ok}, err={err:.3e}")
    
    return ComparisonResult(
        ok=ok,
        err=err,
        details={
            'EwaldSR_f': EwaldSR_f,
            'EwaldSR_o': EwaldSR_o,
            'ewaldsr4_o': ewaldsr4_o,
            'meta': res,
        }
    )


def compare_dip(ham, atomPos, atomTypes_Z, sd, dims, dip4_f, tol=1e-6, verbose=False):
    res = ham.assemble_dip_blocks_per_slot(atomPos, atomTypes_Z, sd)
    dip4_o = res['dip4']
    # Compare only on slots we actually computed (others remain zero in dip4_o)
    dip4_f_use = dip4_f.copy()
    natoms = int(dims.natoms)
    for ia in range(natoms):
        nn = int(sd.neighn[ia])
        for ineigh in range(nn):
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = int(sd.neigh_b[ia, ineigh])
            if ja < 0 or ja >= natoms:
                continue
            if ia == ja and mb == 0:
                dip4_f_use[ia, ineigh, :, :] = 0.0
    Dip_f = _blocked_to_dense(sd, dip4_f_use, natoms)
    Dip_o = _blocked_to_dense(sd, dip4_o, natoms)
    ok, err = compare_matrices_brief("Dip (OCL helper)", Dip_f, Dip_o, tol=tol, require_nonzero=True)
    
    if verbose:
        print(f"  Dip comparison: ok={ok}, err={err:.3e}")
    
    return ComparisonResult(
        ok=ok,
        err=err,
        details={
            'Dip_f': Dip_f,
            'Dip_o': Dip_o,
            'dip4_o': dip4_o,
            'dip4_f_use': dip4_f_use,
            'meta': res,
        }
    )


def compare_vnl_ref(ham, atomPos, atomTypes_Z, sd, s_map, cls, n_orb_atom, offs, use_gpu_vnl=False, verbose=False):
    """Compare VNL reference blocks with OpenCL implementation.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It computes VNL reference blocks by contracting sVNL blocks with projector coefficients,
    then compares with OpenCL CPU/GPU implementations.

    Parameters
    ----------
    ham : OCL_Hamiltonian instance
    atomPos : (natoms,3) float
    atomTypes_Z : (natoms,) int
    sd : FireCore sparse-data struct
    s_map : dict mapping (i,k) -> sVNL block (4,4)
    cls : list of projector coefficients per atom
    n_orb_atom : (natoms,) int
    offs : (natoms+1,) int
    use_gpu_vnl : bool, whether to use GPU VNL contraction
    verbose : bool, whether to print detailed comparison info

    Returns
    -------
    ComparisonResult with ok, err, and details dict containing Vnl_ref, Vnl_ocl, neighs_vnl
    """
    natoms = int(atomPos.shape[0])

    # Build VNL neighbor list exactly like Fortran export_mode>=2 does:
    # iterate neighPP list and map each (jatom,mbeta) to an index in the normal neigh list.
    neighs_vnl = []
    npp_tot = 0
    npp_mapped = 0
    for i in range(natoms):
        nn = int(sd.neighn[i])
        npp = int(sd.neighPPn[i])
        npp_tot += npp
        for ipp in range(npp):
            j = int(sd.neighPP_j[ipp, i]) - 1
            b = int(sd.neighPP_b[ipp, i])
            ineigh0 = -1
            for ineigh in range(nn):
                jj = int(sd.neigh_j[i, ineigh]) - 1
                bb = int(sd.neigh_b[i, ineigh])
                if (jj == j) and (bb == b):
                    ineigh0 = ineigh
                    break
            if ineigh0 >= 0:
                npp_mapped += 1
                neighs_vnl.append((i, j))
    neighs_vnl = list(dict.fromkeys(neighs_vnl))
    if len(neighs_vnl) == 0:
        print(f"[VNL_DIAG] neighPPn={sd.neighPPn.tolist()}  total_PP_edges={int(npp_tot)}  mapped_PP_edges={int(npp_mapped)}")
        for i in range(natoms):
            nn = int(sd.neighn[i])
            npp = int(sd.neighPPn[i])
            if npp <= 0:
                continue
            pp_list = [(int(sd.neighPP_j[ipp, i]) - 1, int(sd.neighPP_b[ipp, i])) for ipp in range(npp)]
            nb_list = [(int(sd.neigh_j[i, ineigh]) - 1, int(sd.neigh_b[i, ineigh])) for ineigh in range(nn)]
            print(f"[VNL_DIAG] i={i} pp_list={pp_list} neigh_list={nb_list}")
        raise RuntimeError("VNL: neighs_vnl is empty (no PP neighbors mapped into neigh list)")

    Vnl_ref_blocks = []
    for (i, j) in neighs_vnl:
        acc = np.zeros((4, 4), dtype=np.float64)
        for k in range(natoms):
            acc += contract_vnl_blocks(s_map[(i, k)], s_map[(j, k)], cls[k])
        Vnl_ref_blocks.append(acc)
    Vnl_ref_blocks = np.array(Vnl_ref_blocks, dtype=np.float64)
    Vnl_ref = dense_from_neighbor_list(neighs_vnl, Vnl_ref_blocks, n_orb_atom, offs)

    H_vnl_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighs_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=use_gpu_vnl)
    Vnl_ocl = dense_from_neighbor_list(neighs_vnl, np.array(H_vnl_blocks, dtype=np.float64), n_orb_atom, offs)

    ok, err = compare_matrices_brief(f"VNL ({'GPU' if use_gpu_vnl else 'CPU'} contraction)", Vnl_ref, Vnl_ocl, tol=1e-5, require_nonzero=True)
    
    if verbose:
        print(f"  VNL comparison ({'GPU' if use_gpu_vnl else 'CPU'}): ok={ok}, err={err:.3e}")

    return ComparisonResult(
        ok=ok,
        err=err,
        details={
            'Vnl_ref': Vnl_ref,
            'Vnl_ocl': Vnl_ocl,
            'neighs_vnl': neighs_vnl,
            'Vnl_ref_blocks': Vnl_ref_blocks,
        }
    )


def compare_avgrho(ham, atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, sd, dims, neigh_index, cn_counts, pair_2c_types=None, mu3c_map=None, nu3c_map=None, mv3c_map=None, Qin_shell=None, DEBUG_QIN_TEST=False, verbose=False):
    """Compare AvgRho blocks from OpenCL with Fortran reference.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It calls ham.compute_avg_rho_3c() and compares the results with Fortran reference blocks from sd.rho_off.

    Parameters
    ----------
    ham : OCL_Hamiltonian instance
    atomPos : (natoms,3) float
    pairs : (n_pairs,2) int32
    pair_triplet_types : (n_triplets,) int32
    cn_offsets : (n_triplets,) int32
    cn_indices : (n_triplets,3) int32
    S_blocks : (n_pairs,4,4) float64
    rho_blocks : (n_pairs,4,4) float64
    Qneutral_sh : (nsh_max, nspecies) float
    ispec_of_atom : (natoms,) int32
    nssh_species : (nspecies,) int32
    sd : FireCore sparse-data struct
    dims : FireCore dims struct
    neigh_index : dict mapping (i,j) -> ineigh
    cn_counts : (n_pairs,) int
    pair_2c_types : optional (n_pairs,) int32
    mu3c_map : optional dict
    nu3c_map : optional dict
    mv3c_map : optional dict
    Qin_shell : optional (nsh_max, natoms) float
    DEBUG_QIN_TEST : bool, whether to test with Qin weights
    verbose : bool, whether to print detailed comparison info

    Returns
    -------
    ComparisonResult with ok, err, and details dict containing rho_avg_blocks, ref_blocks, mask
    """
    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, sd.lssh[:dims.nspecies], dims.nsh_max, pair_2c_types=pair_2c_types, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)
    rho_avg_blocks_qin = None
    if DEBUG_QIN_TEST and Qin_shell is not None:
        rho_avg_blocks_qin = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qin_shell, ispec_of_atom, nssh_species, sd.lssh[:dims.nspecies], dims.nsh_max, pair_2c_types=pair_2c_types, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)
        # For itheory=1 SCF path the Fortran reference (average_ca_rho) uses Qin weights.
        # Use the Qin-weighted GPU result as the primary comparison target.
        rho_avg_blocks = rho_avg_blocks_qin

    # Build Fortran reference blocks: rho_off is stored in sd.rho_off in the same blocked layout as sd.rho.
    ref_blocks = np.zeros_like(rho_avg_blocks)
    mask = np.zeros(pairs.shape[0], dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        ineigh = neigh_index.get((int(ia), int(ja)), None)
        if ineigh is None:
            continue
        ref_blocks[ip, :, :] = sd.rho_off[int(ia), ineigh, :4, :4].astype(np.float64)
        mask[ip] = 1

    res_AvgRho = False
    err_AvgRho = float('nan')
    if np.any(mask == 1):
        ok_avg, max_diff_avg = compare_blocks("AvgRho_off (Fortran rho_off vs OpenCL compute_avg_rho)", ref_blocks, rho_avg_blocks, tol=1e-4)
        res_AvgRho = ok_avg
        err_AvgRho = max_diff_avg
        if DEBUG_QIN_TEST and rho_avg_blocks_qin is not None:
            ok_avg_qin, max_diff_avg_qin = compare_blocks("AvgRho_off (Qin weights test)", ref_blocks, rho_avg_blocks_qin, tol=1e-4)
            if verbose:
                print(f"  DEBUG AvgRho Qin test: ok={ok_avg_qin} err={max_diff_avg_qin:.2e}")
    else:
        if verbose:
            print("WARNING: No comparable neighbor blocks found (mask empty)")
    
    if verbose:
        print(f"  AvgRho comparison: ok={res_AvgRho}, err={err_AvgRho:.3e}")

    return ComparisonResult(
        ok=res_AvgRho,
        err=err_AvgRho,
        details={
            'rho_avg_blocks': rho_avg_blocks,
            'ref_blocks': ref_blocks,
            'mask': mask,
            'rho_avg_blocks_qin': rho_avg_blocks_qin,
        }
    )


def check_ewaldsr_decomposition(sd, ewaldsr2cA4_f, ewaldsr2cO4_f, ewaldsr3c4_f, ewaldsr4_f, natoms, verbose=False):
    """Check EwaldSR decomposition: 2c_atom + 2c_ontop + 3c should equal ewaldsr4_f.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It verifies that the EwaldSR components sum correctly.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    ewaldsr2cA4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_2c_atom
    ewaldsr2cO4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_2c_ontop
    ewaldsr3c4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_3c
    ewaldsr4_f : (natoms,neigh_max,4,4) float64 - ewaldsr total
    natoms : int
    verbose : bool, whether to print detailed comparison info

    Returns
    -------
    ComparisonResult with ok, err, and details dict containing component matrices
    """
    # NOTE: 2c_atom is accumulated in self-slot (matom=neigh_self(iatom)), so it needs the same dense reconstruction as vca_atom.
    EwaldSR_2c_atom_f = _blocked_to_dense_vca_atom(sd, ewaldsr2cA4_f, natoms)
    EwaldSR_2c_ontop_f = _blocked_to_dense(sd, ewaldsr2cO4_f, natoms)
    EwaldSR_3c_f = _blocked_to_dense(sd, ewaldsr3c4_f, natoms)
    EwaldSR_sum_f = EwaldSR_2c_atom_f + EwaldSR_2c_ontop_f + EwaldSR_3c_f
    EwaldSR_f = _blocked_to_dense(sd, ewaldsr4_f, natoms)
    max_diff = float(np.max(np.abs(EwaldSR_sum_f - EwaldSR_f)))
    ok = max_diff < 1e-12
    
    if verbose:
        print(f"  EwaldSR decomposition check: ok={ok}, max_diff={max_diff:.3e}")

    return ComparisonResult(
        ok=ok,
        err=max_diff,
        details={
            'EwaldSR_2c_atom_f': EwaldSR_2c_atom_f,
            'EwaldSR_2c_ontop_f': EwaldSR_2c_ontop_f,
            'EwaldSR_3c_f': EwaldSR_3c_f,
            'EwaldSR_sum_f': EwaldSR_sum_f,
            'EwaldSR_f': EwaldSR_f,
        }
    )


def compare_raw3c_rotation(fc, ham, atomPos, atomTypes_Z, sd, dims, ispec_of_atom, Qneutral_sh, n_nz_max, mu3c_map, nu3c_map, mv3c_map, verbose=False):
    """Compare raw 3c interpolation and rotation between Fortran and OpenCL.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It checks the raw 3c output, hlist construction, SP recovery, and rotation.

    Parameters
    ----------
    fc : FireCore instance
    ham : OCL_Hamiltonian instance
    atomPos : (natoms,3) float
    atomTypes_Z : (natoms,) int
    sd : FireCore sparse-data struct
    dims : FireCore dims struct
    ispec_of_atom : (natoms,) int32
    Qneutral_sh : (nsh_max, nspecies) float
    n_nz_max : int
    mu3c_map : (1,n_nz_max) int16
    nu3c_map : (1,n_nz_max) int16
    mv3c_map : (1,n_nz_max) int8
    verbose : bool, whether to print detailed debug info

    Returns
    -------
    ComparisonResult with ok, err, and details dict containing rotation errors
    """
    from .OCL_Hamiltonian import _epsilon_fb_py, _recover_sp, _recover_rotate_sp
    
    natoms = int(atomPos.shape[0])
    if natoms < 3:
        return ComparisonResult(ok=True, err=0.0, details={'skipped': 'natoms < 3'})
    
    # Select a single triplet for debugging
    i = 0; j = 2; k = 1
    dRj = (atomPos[j] - atomPos[i]).copy()
    dRk = (atomPos[k] - atomPos[i]).copy()
    in1 = int(ispec_of_atom[i]) + 1
    in2 = int(ispec_of_atom[j]) + 1
    indna = int(ispec_of_atom[k]) + 1
    nz1 = int(atomTypes_Z[i]); nz2 = int(atomTypes_Z[j]); nz3 = int(atomTypes_Z[k])
    nssh_k = int(sd.nssh[indna - 1])
    
    if verbose:
        print(f"\n[DEBUG] 3c raw check pair({i},{j}) cn={k} nssh_k={nssh_k}")
        dims_dbg = fc.get_HS_dims(force_refresh=True)
        print(f"[DEBUG] max_mu_dim3(dims)={dims_dbg.max_mu_dim3}  n_nz_3c_max(ham)={ham.n_nz_3c_max}  n_nz_max(local)={n_nz_max}")
    
    # Compute geometry
    r1 = atomPos[i]
    r2 = atomPos[j]
    r21 = r2 - r1
    y = np.linalg.norm(r21)
    sighat = r21 / y if y > 1e-12 else np.array([0.0, 0.0, 1.0])
    mid = 0.5 * (r1 + r2)
    rnabc = atomPos[k] - mid
    x = np.linalg.norm(rnabc)
    rhat = rnabc / x if x > 1e-12 else np.array([0.0, 0.0, 0.0])
    cost = float(np.dot(sighat, rhat))
    cost2 = cost * cost
    argument = 1.0 - cost2
    if argument < 1e-5:
        argument = 1e-5
    sint = np.sqrt(argument)
    P = np.zeros(5, dtype=np.float64)
    P[0] = 1.0
    P[1] = cost
    P[2] = (3.0 * cost2 - 1.0) * 0.5
    P[3] = (5.0 * cost2 * cost - 3.0 * cost) * 0.5
    P[4] = (35.0 * cost2 * cost2 - 30.0 * cost2 + 3.0) * 0.125
    eps = _epsilon_fb_py(rhat, sighat)
    
    if verbose:
        print("  [PY] eps(3x3):")
        print(eps)
        nssh_in1 = int(sd.nssh[in1-1])
        print(f"  [PY] nssh in1={in1} nssh[in1]={nssh_in1} lssh_row={sd.lssh[in1-1, :nssh_in1]}")
        nssh_in2 = int(sd.nssh[in2-1])
        print(f"  [PY] nssh in2={in2} nssh[in2]={nssh_in2} lssh_row={sd.lssh[in2-1, :nssh_in2]}")
    
    block_f_sum = np.zeros((4, 4), dtype=np.float64)
    block_o_sum = np.zeros((4, 4), dtype=np.float64)
    err_rot_isorp = 0.0
    
    for isorp in range(1, nssh_k + 1):
        out_f = np.zeros((5, n_nz_max), dtype=np.float64, order='F')
        raw_f = fc.scanHamPiece3c_raw(3, isorp, in1, in2, indna, dRj, dRk, out=out_f)
        raw_o = ham.scanHamPiece3c_raw_batch('den3', nz1, nz2, nz3, np.array([dRj], dtype=np.float32), np.array([dRk], dtype=np.float32), isorp=isorp)
        
        if raw_o is None:
            if verbose:
                print(f"  raw3c isorp={isorp}: missing OCL table")
            continue
        
        raw_o0 = raw_o[0]
        
        if verbose:
            print(f"  [PY] raw_f (theta x ME) isorp={isorp}:")
            print(raw_f[:, :min(10, raw_f.shape[1])])
            print(f"  [PY] raw_o (theta x ME) isorp={isorp}:")
            print(raw_o0[:, :min(10, raw_o0.shape[1])])
        
        n_me = min(raw_f.shape[1], raw_o0.shape[1], n_nz_max)
        diff_raw = np.max(np.abs(raw_f[:, :n_me] - raw_o0[:, :n_me]))
        if verbose:
            print(f"  raw3c isorp={isorp}: max|F-OCL|={diff_raw:.2e}")
        
        # build hlist and rotate
        h_f = np.zeros(n_nz_max, dtype=np.float64)
        h_o = np.zeros(n_nz_max, dtype=np.float64)
        for iME in range(n_me):
            hf = float(np.dot(P, raw_f[:, iME]))
            ho = float(np.dot(P, raw_o0[:, iME]))
            if int(mv3c_map[0, iME]) == 1:
                hf *= sint
                ho *= sint
            h_f[iME] = hf
            h_o[iME] = ho
        
        if verbose:
            print(f"  [PY] hlist_f (first 10) isorp={isorp}:", h_f[:10])
            print(f"  [PY] hlist_o (first 10) isorp={isorp}:", h_o[:10])
        
        m_f = _recover_sp(h_f, mu3c_map[0], nu3c_map[0])
        m_o = _recover_sp(h_o, mu3c_map[0], nu3c_map[0])
        
        if verbose:
            print("  [PY] bcnam_f pre-rot (4x4):")
            print(m_f)
            print("  [PY] bcnam_o pre-rot (4x4):")
            print(m_o)
        
        wf = float(Qneutral_sh[isorp - 1, indna - 1])
        bcnax_f = _recover_rotate_sp(h_f, mu3c_map[0], nu3c_map[0], eps, in1, in2, sd.lssh, sd.nssh)
        bcnax_o = _recover_rotate_sp(h_o, mu3c_map[0], nu3c_map[0], eps, in1, in2, sd.lssh, sd.nssh)
        bcnax_f_sd = _recover_rotate_sp(h_f, sd.mu[in2-1, in1-1, :n_nz_max].astype(np.int16), sd.nu[in2-1, in1-1, :n_nz_max].astype(np.int16), eps, in1, in2, sd.lssh, sd.nssh)
        
        if verbose:
            print(f"  [PY] bcnax_f post-rot isorp={isorp}:")
            print(bcnax_f)
            print(f"  [PY] bcnax_o post-rot isorp={isorp}:")
            print(bcnax_o)
            print(f"  [PY] bcnax_f_sdmap post-rot isorp={isorp}:")
            print(bcnax_f_sd)
        
        ref_is = fc.scanHamPiece3c(3, isorp, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4)
        ref_is_T = ref_is.T.copy()
        
        if verbose:
            print(f"  [PY] scanHamPiece3c block isorp={isorp} (raw, col-major):")
            print(ref_is)
            print(f"  [PY] scanHamPiece3c block isorp={isorp} (transposed):")
            print(ref_is_T)
        
        err_is = float(np.max(np.abs(bcnax_f - ref_is_T)))
        err_rot_isorp = max(err_rot_isorp, err_is)
        if verbose:
            print(f"  [PY] max|bcnax_f - scanHamPiece3c(T)| isorp={isorp} = {err_is:.3e}")
            err_sdmap = float(np.max(np.abs(bcnax_f_sd - ref_is_T)))
            print(f"  [PY] max|bcnax_f_sdmap - scanHamPiece3c(T)| isorp={isorp} = {err_sdmap:.3e}")
        
        block_f_sum += bcnax_f * wf
        block_o_sum += bcnax_o * wf
    
    block_ref = fc.scanHamPiece3c(3, 1, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4).T * float(Qneutral_sh[0, indna - 1])
    if nssh_k > 1:
        block_ref += fc.scanHamPiece3c(3, 2, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4).T * float(Qneutral_sh[1, indna - 1])
    
    err_rawrot_f = float(np.max(np.abs(block_f_sum - block_ref)))
    err_rawrot_o = float(np.max(np.abs(block_o_sum - block_ref)))
    
    if verbose:
        print(f"  raw->rot (Fraw) vs Fortran scan: max|diff|={err_rawrot_f:.2e}")
        print(f"  raw->rot (OCLraw) vs Fortran scan: max|diff|={err_rawrot_o:.2e}")
    
    ok = (err_rawrot_f < 1e-6) and (err_rawrot_o < 1e-6)
    
    return ComparisonResult(
        ok=ok,
        err=max(err_rawrot_f, err_rawrot_o),
        details={
            'err_rawrot_f': err_rawrot_f,
            'err_rawrot_o': err_rawrot_o,
            'err_rot_isorp': err_rot_isorp,
            'block_ref': block_ref,
            'block_f_sum': block_f_sum,
            'block_o_sum': block_o_sum,
        }
    )

def _scan_table_2c(interaction_fortran, root_ocl, applyRotation_offdiag=True, applyRotation_self=False, in3=1):
    blocks_f = []
    blocks_o = []
    for (i, j) in neighs_all:
        dR = (atomPos[j] - atomPos[i]).copy()
        if i == j:
            dR[:] = 0.0
            Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_self)
            Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_self)
        else:
            Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_offdiag)
            Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_offdiag)
        blocks_f.append(Af)
        blocks_o.append(Ao)
    blocks_f = np.array(blocks_f, dtype=np.float64)
    blocks_o = np.array(blocks_o, dtype=np.float64)
    Mf = dense_from_neighbor_list(neighs_all, blocks_f, n_orb_atom, offs)
    Mo = dense_from_neighbor_list(neighs_all, blocks_o, n_orb_atom, offs)
    return Mf, Mo

def _scan_table_2c(interaction_fortran, root_ocl, applyRotation_offdiag=True, applyRotation_self=False, in3=1):
    blocks_f = []
    blocks_o = []
    for (i, j) in neighs_all:
        dR = (atomPos[j] - atomPos[i]).copy()
        if i == j:
            dR[:] = 0.0
            Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_self)
            Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_self)
        else:
            Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_offdiag)
            Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_offdiag)
        blocks_f.append(Af)
        blocks_o.append(Ao)
    blocks_f = np.array(blocks_f, dtype=np.float64)
    blocks_o = np.array(blocks_o, dtype=np.float64)
    Mf = dense_from_neighbor_list(neighs_all, blocks_f, n_orb_atom, offs)
    Mo = dense_from_neighbor_list(neighs_all, blocks_o, n_orb_atom, offs)
    return Mf, Mo

def print_ewald_debug(fc, sd, atomPos, neigh_self, dip4_f, eq2=14.39975, verbose=False):
    """Print detailed EwaldSR debug information matching Fortran [EW2C_A]/[EW2C_O]/[EW3C] lines.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It prints detailed debug information for EwaldSR calculations.

    Parameters
    ----------
    fc : FireCore instance
    sd : FireCore sparse-data struct
    atomPos : (natoms,3) float
    neigh_self : (natoms,) int
    dip4_f : (natoms,neigh_max,4,4) float64
    eq2 : float, Coulomb constant
    verbose : bool, whether to print detailed debug info
    """
    if not verbose:
        return

    from .OCL_Hamiltonian import emit_ewald_2c_debug, emit_ewald_3c_debug

    natoms = int(atomPos.shape[0])

    dims = fc.get_HS_dims(force_refresh=True)
    Qin_shell = fc.get_Qin_shell(dims)
    Qneutral_sh = fc.get_Qneutral_shell(dims)
    ispec_of_atom = np.zeros(natoms, dtype=np.int32)
    for ia in range(natoms):
        Z = int(sd.iatyp[ia])
        w = np.where(np.array(sd.nzx, dtype=np.int32) == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list")
        ispec_of_atom[ia] = int(w[0])

    def _dq_atom(ia):
        ispec = int(ispec_of_atom[ia])
        s = 0.0
        for issh in range(int(sd.nssh[ispec])):
            s += float(Qin_shell[issh, ia] - Qneutral_sh[issh, ispec])
        return s

    n_orb_atom, offs = _orbital_layout(sd, natoms)

    lines_2c = emit_ewald_2c_debug(sd, atomPos, neigh_self, dip4_f, _dq_atom, natom_cap=3, eq2=eq2)
    lines_3c = emit_ewald_3c_debug(sd, atomPos, dip4_f, _dq_atom, n_orb_atom, offs, natom_cap=3, eq2=eq2)

    for line in lines_2c + lines_3c:
        print(line)


def apply_sfire_overwrite(ewaldsr4_sum_f, ewaldsr2cA4_f, ewaldsr2cO4_f, ewaldsr3c4_f, sd, neigh_back, n_orb_atom, verbose=False):
    """Apply SFIRE overwrite semantics to EwaldSR 4D blocks.

    This is a verification-oriented routine extracted from tests/pyFireball/verify_C3.py.
    It applies the SFIRE overwrite: ewaldsr(jneigh,jatom) = ewaldsr(ineigh,iatom).T

    Parameters
    ----------
    ewaldsr4_sum_f : (natoms,neigh_max,4,4) float64 - sum of 2c_atom + 2c_ontop + 3c
    ewaldsr2cA4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_2c_atom
    ewaldsr2cO4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_2c_ontop
    ewaldsr3c4_f : (natoms,neigh_max,4,4) float64 - ewaldsr_3c
    sd : FireCore sparse-data struct
    neigh_back : (natoms,neigh_max) int
    n_orb_atom : (natoms,) int
    verbose : bool, whether to print debug info

    Returns
    -------
    ComparisonResult with ok, err, and details dict containing overwritten ewaldsr4_sum_f
    """
    from .OCL_Hamiltonian import enforce_sfire_symmetry

    n_sfire_apply, n_sfire_skip = enforce_sfire_symmetry(ewaldsr4_sum_f, ewaldsr3c4_f, sd, neigh_back, n_orb_atom)

    if verbose:
        print(f"  [EWALD_SPLIT_4D] SFIRE overwrite applied={n_sfire_apply} skipped_bad_map={n_sfire_skip}")

    ok = (float(np.max(np.abs(ewaldsr4_sum_f - (ewaldsr2cA4_f + ewaldsr2cO4_f + ewaldsr3c4_f))) ) < 1e-12)
    err = float(np.max(np.abs(ewaldsr4_sum_f - (ewaldsr2cA4_f + ewaldsr2cO4_f + ewaldsr3c4_f))))

    return ComparisonResult(
        ok=ok,
        err=err,
        details={
            'n_sfire_apply': n_sfire_apply,
            'n_sfire_skip': n_sfire_skip,
            'ewaldsr4_sum_f': ewaldsr4_sum_f,
        }
    )

def validate_neigh_back(sd, neigh_back, dims):
    """Validate neigh_back mapping consistency.

    This helper checks that neigh_back correctly maps (iatom,ineigh)->jneigh
    in the neighbor list of jatom such that neigh_j[jatom, jneigh] == iatom+1.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    neigh_back : (natoms,neigh_max) int
    dims : FireCore dimensions struct

    Returns
    -------
    nb_bad : int, number of bad mappings found
    """
    nb_min = int(np.min(neigh_back))
    nb_max = int(np.max(neigh_back))
    print(f"[NEIGH_BACK] min={nb_min} max={nb_max} (expect 0 or 1..neigh_max)")
    nb_bad = 0
    for i in range(int(dims.natoms)):
        nn = int(sd.neighn[i])
        for ineigh in range(nn):
            j = int(sd.neigh_j[i, ineigh]) - 1
            jb = int(neigh_back[i, ineigh])
            if j < 0 or j >= int(dims.natoms):
                continue
            if jb <= 0:
                nb_bad += 1
                continue
            jneigh = jb - 1
            if jneigh < 0 or jneigh >= int(sd.neighn[j]):
                nb_bad += 1
                continue
            if int(sd.neigh_j[j, jneigh]) - 1 != i:
                nb_bad += 1
    return nb_bad

def summarize_neigh_lists(sd, dims):
    """Build neighbor lists from Fortran export.

    This helper constructs neighbor lists for each atom from the Fortran
    sparse data structure.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    dims : FireCore dimensions struct

    Returns
    -------
    neigh_lists : list of list, neighbor lists (excluding self)
    neigh_lists_self : list of list, neighbor lists (including self)
    """
    neigh_lists = []
    neigh_lists_self = []
    for ia in range(dims.natoms):
        nn = int(sd.neighn[ia])
        js = []
        js_self = []
        for ineigh in range(nn):
            j = int(sd.neigh_j[ia, ineigh]) - 1
            if (j >= 0) and (j < dims.natoms):
                js_self.append(j)
                if j != ia:
                    js.append(j)
        neigh_lists.append(sorted(set(js)))
        neigh_lists_self.append(sorted(set(js_self)))
    return neigh_lists, neigh_lists_self

def build_overlap_pairs(sd, atomPos, ham, pairs, neigh_index, ispec_of_atom, atomTypes_Z):
    """Build overlap pair data structures for verification.

    This helper extracts S blocks, pair types, and geometry information
    for all directed pairs.

    Parameters
    ----------
    sd : FireCore sparse-data struct
    atomPos : (natoms,3) float
    ham : OCL_Hamiltonian instance
    pairs : (npairs,2) int, directed atom pairs
    neigh_index : dict, mapping (i,j) -> ineigh
    ispec_of_atom : (natoms,) int, species index per atom
    atomTypes_Z : (natoms,) int, atomic number per atom

    Returns
    -------
    S_blocks : (npairs,4,4) float64, overlap blocks
    pair_2c_types : (npairs,) int, pair type indices
    pair_dR : (npairs,3) float, displacement vectors
    pair_in_i : (npairs,) int, species index for atom i
    pair_in_j : (npairs,) int, species index for atom j
    pair_valid : (npairs,) int, validity flags
    pair_mbeta : (npairs,) int, mbeta values
    """
    npairs = pairs.shape[0]
    S_blocks = np.zeros((npairs, 4, 4), dtype=np.float64)
    pair_2c_types = np.zeros(npairs, dtype=np.int32)
    pair_dR = np.zeros((npairs, 3), dtype=np.float64)
    pair_in_i = np.zeros(npairs, dtype=np.int32)
    pair_in_j = np.zeros(npairs, dtype=np.int32)
    pair_valid = np.zeros(npairs, dtype=np.int32)
    pair_mbeta = np.zeros(npairs, dtype=np.int32)

    for ip, (ia, ja) in enumerate(pairs):
        ineigh = neigh_index.get((int(ia), int(ja)), None)
        if ineigh is None:
            continue
        S_blocks[ip, :, :] = sd.s_mat[int(ia), ineigh, :4, :4]

        ispec_i = int(ispec_of_atom[int(ia)])
        ispec_j = int(ispec_of_atom[int(ja)])
        nz1 = int(atomTypes_Z[int(ia)])
        nz2 = int(atomTypes_Z[int(ja)])
        pt = ham._resolve_pair_type('overlap', nz1, nz2)
        if pt is None:
            raise RuntimeError(f"Missing 2c pair type for overlap ({nz1},{nz2})")
        pair_2c_types[ip] = int(pt)
        mbeta = int(sd.neigh_b[int(ia), ineigh])
        r1 = atomPos[int(ia)]
        r2 = atomPos[int(ja)] + sd.xl[mbeta]
        dR = (r2 - r1).copy()
        in_i = ispec_i + 1
        in_j = ispec_j + 1
        pair_dR[ip, :] = dR
        pair_in_i[ip] = in_i
        pair_in_j[ip] = in_j
        pair_valid[ip] = 1
        pair_mbeta[ip] = mbeta

    return S_blocks, pair_2c_types, pair_dR, pair_in_i, pair_in_j, pair_valid, pair_mbeta

def assign_triplet_types(ham, pairs, cn_offsets, cn_indices, atomTypes_Z, root='den3'):
    """Assign triplet type indices to pairs.

    This helper determines the triplet type for each pair based on
    common neighbor information.

    Parameters
    ----------
    ham : OCL_Hamiltonian instance
    pairs : (npairs,2) int, directed atom pairs
    cn_offsets : (npairs+1,) int, offsets into cn_indices
    cn_indices : (ncn_total,) int, common neighbor indices
    atomTypes_Z : (natoms,) int, atomic number per atom
    root : str, root key for triplet tables (default 'den3')

    Returns
    -------
    pair_triplet_types : (npairs,) int, triplet type indices
    """
    npairs = pairs.shape[0]
    pair_triplet_types = np.zeros(npairs, dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        nz1 = int(atomTypes_Z[int(ia)]); nz2 = int(atomTypes_Z[int(ja)])
        cn0 = int(cn_offsets[ip]); cn1 = int(cn_offsets[ip + 1])
        if cn1 > cn0:
            k = int(cn_indices[cn0])
            nz3 = int(atomTypes_Z[k])
        else:
            nz3 = nz1
        key = (root, nz1, nz2, nz3)
        if key not in ham.species_triplet_map:
            found = None
            for kk in ham.species_triplet_map.keys():
                if (kk[0] == root) and (kk[1] == nz1) and (kk[2] == nz2):
                    found = kk
                    break
            if found is None:
                raise RuntimeError(f"[PLUMBING] Missing 3c triplet table for {key}")
            key = found
        pair_triplet_types[ip] = int(ham.species_triplet_map[key])
    return pair_triplet_types
