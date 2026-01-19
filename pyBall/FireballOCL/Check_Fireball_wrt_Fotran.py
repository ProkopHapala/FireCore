import numpy as np
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


if __name__ == "__main__":
    import os
    fdata_dir = "/home/prokophapala/git/FireCore/tests/pyFireball/Fdata"
    ocl = OCL_Hamiltonian(fdata_dir)
    
    # H2 molecule test
    species_nz = [1, 1]
    ocl.prepare_splines(species_nz)
    
    ratoms = np.array([[0, 0, 0], [0, 0, 0.74]], dtype=np.float32)
    neighbors = [(0, 1), (1, 0)]
    
    H, S = ocl.assemble_full(ratoms, species_nz, neighbors)
    print("H2 Overlap S[0,1]:\n", S[0,0,0])
    print("H2 Kinetic T[0,1]:\n", H[0,0,0])
