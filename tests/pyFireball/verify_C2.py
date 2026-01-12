import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import OCL_Hamiltonian

np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)


def compare_matrices(name, fortran_mat, ocl_mat, tol=1e-5, require_nonzero=False):
    print(f"\n--- Comparing {name} ---")
    diff = np.abs(fortran_mat - ocl_mat)
    max_diff = np.max(diff)
    max_val = max(np.max(np.abs(fortran_mat)), np.max(np.abs(ocl_mat)))
    print(f"Max difference: {max_diff:.2e}")
    print("Fortran Matrix:")
    print(fortran_mat)
    print("PyOpenCL Matrix:")
    print(ocl_mat)
    print("Abs Diff:")
    print(diff)

    if require_nonzero and max_val < tol:
        print(f"WARNING: {name} is (near) zero for both; treating as failure.")
        return False
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        return False
    print(f"SUCCESS: {name} matches.")
    return True


def run_verification():
    # C2 molecule (simple dimer). Distance ~1.25 A (approx CC bond order between single/double)
    atomTypes_Z = np.array([6, 6], dtype=np.int32)
    atomPos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.25],
        ],
        dtype=np.float64,
    )

    fdata_dir = "./Fdata"

    print("Initializing Fortran FireCore...")
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    fc.setVerbosity(0)

    # Protocol: run SCF once (full Hamiltonian, no gating), then keep density frozen.
    # Do NOT call assembleH(Kscf>1) afterwards: assemble_mcweda skips rebuilding 2c pieces for Kscf!=1,
    # which can desynchronize component arrays from h_mat and break reconstruction checks.
    fc.set_export_mode(0)
    fc.set_options(1, 1, 1, 1, 1, 1, 1)
    fc.SCF(positions=atomPos, iforce=0, nmax_scf=50)

    print("Initializing PyOpenCL Hamiltonian...")
    ham = OCL_Hamiltonian(fdata_dir)
    ham.prepare_splines(atomTypes_Z)

    def _dense_from_pair_blocks(B00, B01, B10, B11):
        M = np.zeros((8, 8), dtype=np.float64)
        M[0:4, 0:4] = B00
        M[0:4, 4:8] = B01
        M[4:8, 0:4] = B10
        M[4:8, 4:8] = B11
        return M

    def _scan2c_fortran(interaction, dR, in3=None, isub=0, applyRotation=True):
        in1 = 1
        in2 = 1
        if in3 is None:
            in3 = in2
        return fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation=applyRotation)

    def _scan2c_ocl(root, dR, applyRotation=True):
        return ham.scanHamPiece2c(root, int(atomTypes_Z[0]), int(atomTypes_Z[1]), dR, applyRotation=applyRotation)

    def _orbital_layout(sparse_data, natoms):
        """Return per-atom n_orb and offsets using sparse_data.nzx/num_orb mapping."""
        # nzx: list of Z numbers in fdata species order; num_orb indexed by that species order
        nzx = np.array(sparse_data.nzx, dtype=np.int32)
        num_orb_species = np.array(sparse_data.num_orb, dtype=np.int32)
        iatyp = np.array(sparse_data.iatyp, dtype=np.int32)
        n_orb_atom = np.zeros(natoms, dtype=np.int32)
        for ia in range(natoms):
            Z = iatyp[ia]
            w = np.where(nzx == Z)[0]
            if w.size == 0:
                raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
            ispec = int(w[0])
            n_orb_atom[ia] = int(num_orb_species[ispec])
        if np.any(n_orb_atom == 0):
            raise RuntimeError(
                f"Zero-orbital atom encountered: n_orb_atom={n_orb_atom}, iatyp={iatyp}, "
                f"num_orb_species={num_orb_species}, nzx={nzx}"
            )
        offs = np.zeros(natoms + 1, dtype=np.int32)
        offs[1:] = np.cumsum(n_orb_atom)
        return n_orb_atom, offs

    def _blocked_to_dense(sparse_data, H_blocks, natoms):
        """Map FireCore blocked sparse (iatom, ineigh, inu, imu) into dense matrix."""
        n_orb_atom, offs = _orbital_layout(sparse_data, natoms)
        norb = int(offs[-1])
        M = np.zeros((norb, norb), dtype=np.float64)
        neighn = np.array(sparse_data.neighn, dtype=np.int32)
        neigh_j = np.array(sparse_data.neigh_j, dtype=np.int32)
        for i in range(natoms):
            ni = int(n_orb_atom[i])
            i0 = int(offs[i])
            for ineigh in range(int(neighn[i])):
                j = int(neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                nj = int(n_orb_atom[j])
                j0 = int(offs[j])
                # H_blocks is [iatom, ineigh, inu, imu] (Fortran order preserved in Python wrapper)
                blk = H_blocks[i, ineigh, :nj, :ni]
                # rows=imu (i), cols=inu (j)
                M[i0:i0+ni, j0:j0+nj] = blk.T
        return M

    def get_fortran_sparse_full(export_mode=0):
        fc.set_options(1, 1, 1, 1, 1, 1, 1)
        fc.set_export_mode(export_mode)
        dims = fc.get_HS_dims()
        sparse_data = fc.get_HS_neighs(dims)
        sparse_data = fc.get_HS_sparse(dims, sparse_data)
        natoms = 2
        H = _blocked_to_dense(sparse_data, sparse_data.h_mat, natoms)
        S = _blocked_to_dense(sparse_data, sparse_data.s_mat, natoms)
        neighbors = []
        for i in range(natoms):
            for ineigh in range(int(sparse_data.neighn[i])):
                j = int(sparse_data.neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                neighbors.append((i, j))
        return H, S, neighbors, sparse_data

    def get_fortran_sparse_H_with_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, export_mode=1):
        fc.set_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca)
        fc.set_export_mode(export_mode)
        dims = fc.get_HS_dims()
        sparse_data = fc.get_HS_neighs(dims)
        sparse_data = fc.get_HS_sparse(dims, sparse_data)
        natoms = 2
        H = _blocked_to_dense(sparse_data, sparse_data.h_mat, natoms)
        neighbors = []
        for i in range(natoms):
            for ineigh in range(int(sparse_data.neighn[i])):
                j = int(sparse_data.neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                neighbors.append((i, j))
        return H, neighbors, sparse_data

    print("\nTesting Overlap S...")
    dR01 = (atomPos[1] - atomPos[0]).copy()
    dR10 = (atomPos[0] - atomPos[1]).copy()
    S00_f = _scan2c_fortran(1, np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    S11_f = S00_f.copy()
    S01_f = _scan2c_fortran(1, dR01, applyRotation=True)
    S10_f = _scan2c_fortran(1, dR10, applyRotation=True)
    S_f = _dense_from_pair_blocks(S00_f, S01_f, S10_f, S11_f)
    S00_o = _scan2c_ocl('overlap', np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    S11_o = S00_o.copy()
    S01_o = _scan2c_ocl('overlap', dR01, applyRotation=True)
    S10_o = _scan2c_ocl('overlap', dR10, applyRotation=True)
    S_o = _dense_from_pair_blocks(S00_o, S01_o, S10_o, S11_o)
    res_S = compare_matrices("Overlap S", S_f, S_o)

    print("\nTesting Kinetic T...")
    T00_f = _scan2c_fortran(13, np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    T11_f = T00_f.copy()
    T01_f = _scan2c_fortran(13, dR01, applyRotation=True)
    T10_f = _scan2c_fortran(13, dR10, applyRotation=True)
    T_f = _dense_from_pair_blocks(T00_f, T01_f, T10_f, T11_f)
    T00_o = _scan2c_ocl('kinetic', np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    T11_o = T00_o.copy()
    T01_o = _scan2c_ocl('kinetic', dR01, applyRotation=True)
    T10_o = _scan2c_ocl('kinetic', dR10, applyRotation=True)
    T_o = _dense_from_pair_blocks(T00_o, T01_o, T10_o, T11_o)
    res_T = compare_matrices("Kinetic T", T_f, T_o)

    print("\nTesting Vna...")
    Vna01_f = _scan2c_fortran(2, dR01, in3=1, applyRotation=True) + _scan2c_fortran(3, dR01, in3=1, applyRotation=True)
    Vna10_f = _scan2c_fortran(2, dR10, in3=1, applyRotation=True) + _scan2c_fortran(3, dR10, in3=1, applyRotation=True)
    Vna00_f = _scan2c_fortran(4, dR01, in3=1, applyRotation=True)
    Vna11_f = _scan2c_fortran(4, dR10, in3=1, applyRotation=True)
    Vna_f = _dense_from_pair_blocks(Vna00_f, Vna01_f, Vna10_f, Vna11_f)
    Vna01_o = _scan2c_ocl('vna_ontopl_00', dR01, applyRotation=True) + _scan2c_ocl('vna_ontopr_00', dR01, applyRotation=True)
    Vna10_o = _scan2c_ocl('vna_ontopl_00', dR10, applyRotation=True) + _scan2c_ocl('vna_ontopr_00', dR10, applyRotation=True)
    Vna00_o = _scan2c_ocl('vna_atom_00', dR01, applyRotation=True)
    Vna11_o = _scan2c_ocl('vna_atom_00', dR10, applyRotation=True)
    Vna_o = _dense_from_pair_blocks(Vna00_o, Vna01_o, Vna10_o, Vna11_o)
    res_Vna = compare_matrices("Vna", Vna_f, Vna_o)

    print("\nTesting sVNL (PP overlap, interaction=5)...")
    sVNL01_f = _scan2c_fortran(5, dR01, in3=1, applyRotation=True)
    sVNL10_f = _scan2c_fortran(5, dR10, in3=1, applyRotation=True)
    sVNL00_f = _scan2c_fortran(5, dR01, in3=1, applyRotation=True)
    sVNL11_f = _scan2c_fortran(5, dR10, in3=1, applyRotation=True)
    sVNL_f = _dense_from_pair_blocks(sVNL00_f, sVNL01_f, sVNL10_f, sVNL11_f)
    sVNL01_o = _scan2c_ocl('vnl', dR01, applyRotation=True)
    sVNL10_o = _scan2c_ocl('vnl', dR10, applyRotation=True)
    sVNL00_o = _scan2c_ocl('vnl', dR01, applyRotation=True)
    sVNL11_o = _scan2c_ocl('vnl', dR10, applyRotation=True)
    sVNL_o = _dense_from_pair_blocks(sVNL00_o, sVNL01_o, sVNL10_o, sVNL11_o)
    require_nz = np.max(np.abs(sVNL_f)) > 1e-5
    res_Vnl = compare_matrices("sVNL", sVNL_f, sVNL_o, require_nonzero=require_nz)

    print("\nTesting sVNL self (dR=0, applyRotation=False)...")
    z0 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    sVNL_self_f = _scan2c_fortran(5, z0, in3=1, applyRotation=False)
    sVNL_self_o = _scan2c_ocl('vnl', z0, applyRotation=False)
    _ = compare_matrices("sVNL self (0)", sVNL_self_f, sVNL_self_o, tol=1e-6, require_nonzero=True)

    print("\nTesting contracted VNL (scan-based reference vs OpenCL CPU/GPU contraction)...")

    # Neighbor list for a 2-atom system (include self edges)
    neighbors_vnl = [(0, 0), (0, 1), (1, 0), (1, 1)]

    # Build cl vectors from OpenCL-side parser cache (must match Fortran cl_value expansion)
    # For now, enforce sp projectors only (len=4).
    def _cl_sp_from_ham(atom_idx):
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

    cl0 = _cl_sp_from_ham(0)
    cl1 = _cl_sp_from_ham(1)
    cls = [cl0, cl1]

    # Fortran sVNL blocks for all (phi_atom, pp_atom) pairs
    # Use applyRotation=True for off-diagonal geometry; for self (i,i) use dR=0 and applyRotation=False.
    s_map = {}
    for i in (0, 1):
        for k in (0, 1):
            if i == k:
                dR = np.array([0.0, 0.0, 0.0], dtype=np.float64)
                s_map[(i, k)] = _scan2c_fortran(5, dR, in3=1, applyRotation=False)
            else:
                dR = (atomPos[k] - atomPos[i]).copy()
                s_map[(i, k)] = _scan2c_fortran(5, dR, in3=1, applyRotation=True)

    # Contract reference in the same block convention as OpenCL blocks (nu,mu)
    def _contract_blocks(A, B, clv):
        # A,B: (4,4) blocks (nu,mu) => nu=cc, mu=basis index
        # Output: (nu,mu)
        out = np.zeros((4, 4), dtype=np.float64)
        for nu in range(4):
            for mu in range(4):
                v = 0.0
                for cc in range(4):
                    v += A[cc, mu] * clv[cc] * B[cc, nu]
                out[nu, mu] = v
        return out

    # Full VNL reference per (i,j): sum over k
    Vnl_ref_blocks = []
    for (i, j) in neighbors_vnl:
        acc = np.zeros((4, 4), dtype=np.float64)
        for k in (0, 1):
            acc += _contract_blocks(s_map[(i, k)], s_map[(j, k)], cls[k])
        Vnl_ref_blocks.append(acc)
    Vnl_ref_blocks = np.array(Vnl_ref_blocks, dtype=np.float64)

    # Convert blocks to dense 8x8 for comparison
    def _dense_from_neighbor_list(neighs, blocks):
        # neighs: [(i,j)], blocks: [n,4,4] (nu,mu)
        M = np.zeros((8, 8), dtype=np.float64)
        for idx, (i, j) in enumerate(neighs):
            i0 = i * 4
            j0 = j * 4
            # dense uses (mu,nu), so transpose
            M[i0:i0+4, j0:j0+4] = blocks[idx].T
        return M

    Vnl_ref = _dense_from_neighbor_list(neighbors_vnl, Vnl_ref_blocks)
    print("  Reference VNL (scan+contract):")
    print(Vnl_ref)

    # OpenCL: assemble_full returns blocks for the neighbor list in the same (nu,mu) convention
    H_vnl_cpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=False)
    H_vnl_gpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=True)

    Vnl_cpu = _dense_from_neighbor_list(neighbors_vnl, H_vnl_cpu_blocks)
    Vnl_gpu = _dense_from_neighbor_list(neighbors_vnl, H_vnl_gpu_blocks)

    res_Vnl_cpu = compare_matrices("VNL (CPU contraction)", Vnl_ref, Vnl_cpu)
    res_Vnl_gpu = compare_matrices("VNL (GPU contraction)", Vnl_ref, Vnl_gpu)

    print("\nTesting contracted VNL (Fortran gated export vs OpenCL CPU/GPU)...")
    # export_mode=1 reconstructs h_mat like buildh.f90 (does NOT include VNL)
    # export_mode>=2 provides explicit H+VNL export; with T/Vna/Vxc/Vca off we get VNL-only reference.
    Vnl_fortran, _neighs_vnl_f, _sd_vnl_f = get_fortran_sparse_H_with_options(0, 0, 0, 1, 0, 0, 0, export_mode=2)
    res_Vnl_cpu_fortran = compare_matrices("VNL (CPU) vs Fortran export", Vnl_fortran, Vnl_cpu)
    res_Vnl_gpu_fortran = compare_matrices("VNL (GPU) vs Fortran export", Vnl_fortran, Vnl_gpu)

    print("\nTesting Fortran export reconstruction (raw full H vs reconstructed full H)...")
    H_full_raw, S_full_raw, neighbors, sparse_data = get_fortran_sparse_full(export_mode=0)
    H_full_rec, S_full_rec, _neighbors2, _sparse2 = get_fortran_sparse_full(export_mode=1)
    res_H_recon = compare_matrices("Fortran full H: raw vs reconstructed", H_full_raw, H_full_rec, tol=1e-6, require_nonzero=True)
    _ = compare_matrices("Fortran full S: raw vs reconstructed(export)", S_full_raw, S_full_rec, tol=1e-12, require_nonzero=True)

    H_full_f = H_full_raw

    print("\nTesting Vxc (off-site)...")
    H_f = H_full_f
    # Note: OpenCL Vxc not yet implemented; placeholder for future
    # H_o_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=False, include_Vnl=False, include_Vxc=True)
    # n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    # norb = int(offs[-1])
    # Vxc_o = np.zeros((norb, norb), dtype=np.float64)
    # for idx, (i, j) in enumerate(neighbors):
    #     ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
    #     i0 = int(offs[i]); j0 = int(offs[j])
    #     Vxc_o[i0:i0+ni, j0:j0+nj] = H_o_blocks[idx, :nj, :ni].T
    # res_Vxc = compare_matrices("Vxc", H_f, Vxc_o)
    print("  Fortran Vxc matrix:")
    print(H_f)
    print("  SKIPPED: OpenCL Vxc not yet implemented")
    res_Vxc = False  # Placeholder

    print("\nTesting Vxc_1c (on-site)...")
    # Vxc_1c is part of Vxc in Fortran; we can extract diagonal blocks
    # For now, just print the Fortran diagonal
    n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    norb = int(offs[-1])
    diag_Vxc_1c_f = np.zeros(norb)
    for i in range(2):
        ni = int(n_orb_atom[i])
        i0 = int(offs[i])
        diag_Vxc_1c_f[i0:i0+ni] = np.diag(H_f[i0:i0+ni, i0:i0+ni])
    print("  Fortran Vxc_1c diagonal:")
    print(diag_Vxc_1c_f)
    print("  SKIPPED: OpenCL Vxc_1c not yet implemented")
    res_Vxc_1c = False  # Placeholder

    print("\nTesting Vca (charge-dependent)...")
    H_f = H_full_f
    print("  Fortran Vca matrix:")
    print(H_f)
    print("  SKIPPED: OpenCL Vca not yet implemented")
    res_Vca = False  # Placeholder

    print("\nTesting Vxc_ca (charge-dependent XC)...")
    H_f = H_full_f
    print("  Fortran Vxc_ca matrix:")
    print(H_f)
    print("  SKIPPED: OpenCL Vxc_ca not yet implemented")
    res_Vxc_ca = False  # Placeholder

    print("\nTesting combined 2-center H = T + Vna (from scan API)...")
    H2c_f = T_f + Vna_f
    H2c_o = T_o + Vna_o
    res_H2c = compare_matrices("H2c (T+Vna)", H2c_f, H2c_o)

    print("\nTesting Full Hamiltonian H (Fortran full vs OpenCL partial) ...")
    print("  SKIPPED: Fortran H contains additional terms (XC/Vca/Vxc_ca and full Vnl contraction) not yet implemented/comparable in OpenCL")
    res_H_partial = False
    res_H_full = False

    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY")
    print(f"Overlap S:   {'PASSED' if res_S else 'FAILED'}")
    print(f"Kinetic T:   {'PASSED' if res_T else 'FAILED'}")
    print(f"Vna:         {'PASSED' if res_Vna else 'FAILED'}")
    print(f"Vnl:         {'PASSED' if res_Vnl else 'FAILED'}")
    print(f"Vxc:         {'PASSED' if res_Vxc else 'NOT IMPLEMENTED'}")
    print(f"Vxc_1c:      {'PASSED' if res_Vxc_1c else 'NOT IMPLEMENTED'}")
    print(f"Vca:         {'PASSED' if res_Vca else 'NOT IMPLEMENTED'}")
    print(f"Vxc_ca:      {'PASSED' if res_Vxc_ca else 'NOT IMPLEMENTED'}")
    print(f"H2c (T+Vna): {'PASSED' if res_H2c else 'FAILED'}")
    print(f"H raw==recon:{'PASSED' if res_H_recon else 'FAILED'}")
    print(f"VNL vs F:    {'PASSED' if (res_Vnl_cpu_fortran and res_Vnl_gpu_fortran) else 'FAILED'}")
    print(f"Full H:      {'PASSED' if res_H_full else 'NOT IMPLEMENTED'}")
    print("=" * 40)


if __name__ == "__main__":
    run_verification()
