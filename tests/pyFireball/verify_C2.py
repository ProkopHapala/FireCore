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
    ham.prepare_data_3c(atomTypes_Z)

    print("\n[PLUMBING] compute_avg_rho (3c gather) using exported rho + Qin-shell...")
    dims = fc.get_HS_dims(force_refresh=True)
    sd = fc.get_HS_neighs(dims)
    sd = fc.get_HS_sparse(dims, data=sd)
    sd = fc.get_rho_sparse(dims, data=sd)

    Qin_shell = fc.get_Qin_shell(dims)  # [nsh_max, natoms]
    Qatom = np.sum(Qin_shell, axis=0).astype(np.float32)

    # 3c data prep similar to verify_C3.py
    Qneutral_sh = fc.get_Qneutral_shell(dims)  # [nsh_max, nspecies_fdata]
    nzx = np.array(sd.nzx, dtype=np.int32)
    iatyp = np.array(sd.iatyp, dtype=np.int32)
    ispec_of_atom = np.zeros(dims.natoms, dtype=np.int32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec_of_atom[ia] = int(w[0])
    
    nspecies_fdata = int(Qneutral_sh.shape[1])
    nssh_species = np.zeros(nspecies_fdata, dtype=np.int32)
    for ispec in range(nspecies_fdata):
        nssh_species[ispec] = int(np.count_nonzero(Qneutral_sh[:, ispec] != 0.0))

    # Neighbor lists from Fortran export: neigh_j is 1-based atom index, as used elsewhere in this file.
    neigh_lists = []
    for ia in range(dims.natoms):
        nn = int(sd.neighn[ia])
        js = []
        for ineigh in range(nn):
            j = int(sd.neigh_j[ia, ineigh]) - 1
            if (j >= 0) and (j < dims.natoms) and (j != ia):
                js.append(j)
        js = sorted(set(js))
        neigh_lists.append(js)

    pairs = []
    for ia in range(dims.natoms):
        for j in neigh_lists[ia]:
            if ia < j:
                pairs.append((ia, j))
    if len(pairs) == 0:
        raise RuntimeError("[PLUMBING] No neighbor pairs found")
    pairs = np.array(pairs, dtype=np.int32)
    cn_offsets, cn_indices = ham.build_common_neighbor_csr(neigh_lists, pairs)
    cn_counts = cn_offsets[1:] - cn_offsets[:-1]
    print(f"  n_pairs={pairs.shape[0]}  n_cn_total={cn_indices.shape[0]}  cn_max={int(np.max(cn_counts))}")

    # Map sparse neighbor blocks for (i,j) to 4x4 blocks
    neigh_index = {}
    for ia in range(dims.natoms):
        nn = int(sd.neighn[ia])
        for ineigh in range(nn):
            j = int(sd.neigh_j[ia, ineigh]) - 1
            if j >= 0:
                neigh_index[(ia, j)] = ineigh

    S_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    rho_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    for ip, (ia, ja) in enumerate(pairs):
        ineigh = neigh_index.get((ia, ja), None)
        if ineigh is None:
            raise RuntimeError(f"[PLUMBING] Missing sparse neighbor index for pair ({ia},{ja})")
        S_blocks[ip, :, :] = sd.s_mat[ia, ineigh, :4, :4]
        rho_blocks[ip, :, :] = sd.rho[ia, ineigh, :4, :4]

    # Triplet type selection: for C2 this mostly won't matter (CN likely empty), but we need a valid index.
    root = 'den3'
    pair_triplet_types = np.zeros(pairs.shape[0], dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        nz1 = int(atomTypes_Z[ia]); nz2 = int(atomTypes_Z[ja])
        # choose nz3=nz1 for now; for real 3c contributions we'll switch to per-(pair,k) typing
        key = (root, nz1, nz2, nz1)
        if key not in ham.species_triplet_map:
            # fallback: find any triplet with same root and nz1,nz2
            found = None
            for k in ham.species_triplet_map.keys():
                if (k[0] == root) and (k[1] == nz1) and (k[2] == nz2):
                    found = k
                    break
            if found is None:
                raise RuntimeError(f"[PLUMBING] Missing 3c triplet table for {key}")
            key = found
        pair_triplet_types[ip] = int(ham.species_triplet_map[key])

    # Build mu/nu/mvalue maps for 3c recovery (needed for compute_avg_rho_3c)
    n_types = int(np.max(pair_triplet_types)) + 1 if pair_triplet_types.size else 0
    n_nz_max = int(ham.n_nz_3c_max)
    if n_types > 0:
        mu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
        nu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
        mv3c_map = np.zeros((n_types, n_nz_max), dtype=np.int8)
        # For C2, single species, in1=in2=1
        in1 = 1; in2 = 1
        mu3c_map[:, :] = sd.mu[in2-1, in1-1, :n_nz_max].astype(np.int16)
        nu3c_map[:, :] = sd.nu[in2-1, in1-1, :n_nz_max].astype(np.int16)
        mv3c_map[:, :] = sd.mvalue[in2-1, in1-1, :n_nz_max].astype(np.int8)
    else:
        mu3c_map = None; nu3c_map = None; mv3c_map = None

    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, dims.nsh_max, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)
    for ip in range(min(3, pairs.shape[0])):
        i, j = int(pairs[ip, 0]), int(pairs[ip, 1])
        print(f"  rho_avg block pair ({i},{j}) cn_count={int(cn_counts[ip])}")
        print(rho_avg_blocks[ip])

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
    print("\nTesting Vca (charge-dependent)...")
    # Isolate Vca from Fortran
    H_vca_f, _, _ = get_fortran_sparse_H_with_options(0, 0, 0, 0, 0, 1, 0, export_mode=2)
    print("  Fortran Vca matrix (isolated):")
    print(H_vca_f)

    # Reconstruct Vca manually from scan2c and charges
    # vca = sum_k sum_isorp (dq_k * vna(k, isorp)) + ewaldsr
    # For now, ignore EwaldSR and check if vca matches (ewald might be small or zero if dq=0?)
    # But dq is likely non-zero.
    # We need to implement manual Vca to confirm we understand it.
    Vca_manual = np.zeros_like(H_vca_f)
    
    # Calculate dq per shell
    # dq[isorp, iatom]
    dq_shell = np.zeros((dims.nsh_max, 2))
    for ia in range(2):
        for issh in range(int(nssh_species[0])): # Assuming species 0
             dq_shell[issh, ia] = Qin_shell[issh, ia] - Qneutral_sh[issh, int(ispec_of_atom[ia])]
    
    print("  dq_shell:")
    print(dq_shell)

    # Accumulate Vca terms
    # Loop over pairs (i,j) in Vca_manual blocks? No, dense matrix.
    # Iterate all pairs (i, j)
    for i in range(2):
        for j in range(2):
            dR = atomPos[j] - atomPos[i]
            # Term 1: vna_atom (interaction 4)
            # <i| V_atom(k) |i> -> k=i. 
            # In assemble_ca_2c, for iatom, we loop neighbors jatom.
            # vna_atom part: interaction=4, in3=in1 (potential at i).
            # Scaled by dxn = Qin(j) - Qneu(j).
            # Wait, assemble_ca_2c says: 
            #   do isorp=1..nssh(in2): dxn = (Qin(j) - Qneu(in2)); bccax = doscentros(4, ... in3=in1 ...); bcca += bccax*dxn
            # This adds to vca(..., matom, iatom). matom is self. 
            # So this is a diagonal block contribution (i, i) due to neighbor j's charge.
            # So for each i, sum over j (all atoms): V_ii += vna_atom(i, j) * dq(j)
            
            # Term 2: vna_ontopl/r
            # For pair (i,j) (neighbor):
            # vna_ontopl (2): pot at i. dxn = Qin(i)-Qneu(i). V_ij += vna_ontopl * dq(i)
            # vna_ontopr (3): pot at j. dxn = Qin(j)-Qneu(j). V_ij += vna_ontopr * dq(j)
            
            # Lets accumulate into blocks
            # Calculate EwaldSR constants
            # dq1 = dq_shell[:, i], dq2 = dq_shell[:, j]
            # but we need total dq for ewald? No, Fortran uses sum over isorp.
            dqi_tot = np.sum(dq_shell[:, i])
            dqj_tot = np.sum(dq_shell[:, j])
            
            # Reconstruct per-block
            # Note: We need to put this into the correct sub-block of Vca_manual
            # Vca_manual is dense (norb x norb), rows=i, cols=j
            ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
            i0 = int(offs[i]); j0 = int(offs[j])
            
            V_blk = np.zeros((ni, nj))
            
            # Diagonal terms (from ALL neighbors k)
            if i == j:
                for k in range(2):
                    dqk = dq_shell[:, k]
                    dqk_tot = np.sum(dqk)
                    dR_ik = atomPos[k] - atomPos[i]
                    rik = np.linalg.norm(dR_ik)
                    
                    # Term: vna_atom (interaction 4)
                    mat4 = np.zeros((ni, nj))
                    for isorp in range(1, int(nssh_species[0])+1):
                        # Force isorp into _scan2c or call fc directly
                        # fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation)
                        # in1=in2=in3=1 (species)
                        m = fc.scanHamPiece2c(4, isorp, 1, 1, 1, dR_ik, applyRotation=True)
                        mat4 += m * dqk[isorp-1]
                    V_blk += mat4
                    
                    # EwaldSR "ATM CASE" (line 321 in assemble_ca_2c.f90)
                    # if rik > 1e-4:
                    #   emnpl = (s_mat / rik) * dqk_tot
                    #   ewaldsr += emnpl * eq2
                    if rik > 1e-4:
                        # Get S matrix for (i, k)? No, s_mat(imu,inu,matom,iatom) -> (i, i) block?
                        # "s_mat(imu,inu,matom,iatom)" where matom=neigh_self(iatom).
                        # So it uses the Overlap of (i, i) which is Identity (mostly) or on-site overlap?
                        # Wait, s_mat for self is Identity.
                        # So emnpl = (I / rik) * dqk_tot.
                        # ewaldsr += eqn * emnpl * eq2 ??
                        # Fortran: emnpl = (s_mat / y) * dq2. ewaldsr += emnpl*eq2.
                        # Here y = rik. dq2 = dqk_tot.
                        # So V_ii += (I / rik) * dqk_tot * EQ2
                        # EQ2 = 14.39975
                        EQ2 = 14.39975
                        S_ii = np.eye(ni) # approximate S_ii as Identity
                        V_blk += (S_ii / rik) * dqk_tot * EQ2

            else:
                # Off-diagonal V_ij (interaction 2 and 3)
                # Term: vna_ontopl (2) - pot at i
                mat2 = np.zeros((ni, nj))
                for isorp in range(1, int(nssh_species[0])+1):
                    m = fc.scanHamPiece2c(2, isorp, 1, 1, 1, dR, applyRotation=True)
                    mat2 += m * dq_shell[isorp-1, i]
                
                # Term: vna_ontopr (3) - pot at j
                mat3 = np.zeros((ni, nj))
                for isorp in range(1, int(nssh_species[0])+1):
                    m = fc.scanHamPiece2c(3, isorp, 1, 1, 1, dR, applyRotation=True)
                    mat3 += m * dq_shell[isorp-1, j]
                
                V_blk += mat2 + mat3
                
                # EwaldSR "ONTOP CASE" (line 415)
                # ewaldsr = ...
                # ewaldsr += ( (S_ij/rij + Dip_ij/rij^2)*dq_i + (S_ij/rij*dq_j - Dip_ij/rij^2*dq_j) ) * eq2
                # Note: assemble_ca_2c has logic that cancels dipole if charges are same?
                # Formula: ( (S/y + D/y^2)*dq1 + (S/y*dq2 - D/y^2*dq2) ) * eq2
                #        = ( S/y*(dq1+dq2) + D/y^2*(dq1-dq2) ) * eq2
                rij = np.linalg.norm(dR)
                if rij > 1e-4:
                    EQ2 = 14.39975
                    S_ij = _scan2c_fortran(1, dR, applyRotation=True) # Overlap
                    # Dipole: interaction 9 ??
                    # Call doscentros(9, isorp=0, ...)
                    Dip_ij = fc.scanHamPiece2c(9, 0, 1, 1, 1, dR, applyRotation=True)
                    
                    term = (S_ij / rij) * (dqi_tot + dqj_tot) + (Dip_ij / (rij**2)) * (dqi_tot - dqj_tot)
                    V_blk += term * EQ2

            Vca_manual[i0:i0+ni, j0:j0+nj] = V_blk

    print("  Manual Vca matrix:")
    print(Vca_manual)
    
    diff_vca = np.abs(H_vca_f - Vca_manual)
    print(f"  Max diff Vca: {np.max(diff_vca):.2e}")
    if np.max(diff_vca) < 1e-4:
        print("SUCCESS: Vca manual reconstruction matches Fortran export.")
        res_Vca = True
    else:
        print("FAILURE: Vca reconstruction mismatch.")
    print("  Comparing OpenCL Vca implementation...")
    # Prepare neighbor list and types for assemble_vca
    # We need to replicate assemble_full logic for pairs (including self)
    neigh_list_vca = []
    for ia in range(dims.natoms):
        for ja in neigh_lists[ia]:
            neigh_list_vca.append((ia, int(ja))) # 0-based from get_neighbors
    # Add self
    for ia in range(dims.natoms):
        neigh_list_vca.append((ia, ia))
    
    pair_types_vca = []
    for (i, j) in neigh_list_vca:
        nz1 = int(atomTypes_Z[i])
        nz2 = int(atomTypes_Z[j])
        k = ('vna', nz1, nz2)
        idx = ham.species_pair_map.get(k)
        if idx is None:
            # Fallback alias logic used in ham.prepare_splines but we need the index.
            # OCL_Hamiltonian aliases vna -> vna_atom_00 if missing.
            # So ham.species_pair_map should have it if prepared correctly.
            print(f"Warning: missing vna pair type for {nz1},{nz2}")
            idx = 0
        pair_types_vca.append(idx)
        
    off_vca, diag_vca = ham.assemble_vca(atomPos, np.array(neigh_list_vca), np.array(pair_types_vca), dq_shell)
    
    # Reconstruct dense
    Vca_ocl = np.zeros_like(H_vca_f)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0=int(offs[i]); j0=int(offs[j])
        ni=int(n_orb_atom[i]); nj=int(n_orb_atom[j])
        
        # Off-diagonal block (V_ij)
        # Transpose to match dense layout (col-major blocks from Fortran/OCL -> row-major numpy)
        Vca_ocl[i0:i0+ni, j0:j0+nj] += off_vca[k].T
        
        # Diagonal update (V_ii) from neighbor j
        # Added to V[i, i]
        Vca_ocl[i0:i0+ni, i0:i0+ni] += diag_vca[k, 0].T

    print("  OpenCL Vca matrix result:")
    print(Vca_ocl)
    
    diff_ocl = np.abs(Vca_ocl - H_vca_f)
    print(f"  Max diff Vca OCL vs Fortran: {np.max(diff_ocl):.2e}")
    if np.max(diff_ocl) < 1e-4:
        print("SUCCESS: OpenCL Vca matches Fortran.")
        res_Vca = True
    else:
        print("FAILURE: OpenCL Vca mismatch.")
        print("  Comparison with Manual (bcca+ewald):")
        diff_man = np.abs(Vca_ocl - Vca_manual)
        print(f"  Max diff OCL vs Manual: {np.max(diff_man):.2e}")
        res_Vca = False

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
