import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import (
    OCL_Hamiltonian,
    compare_matrices as _compare_matrices_shared,
    dense_from_pair_blocks,
    scan2c_fortran,
    scan2c_ocl,
    _orbital_layout as _orbital_layout_shared,
    _blocked_to_dense as _blocked_to_dense_shared,
    firecore_sparse_to_dense,
    firecore_sparse_H_with_options,
    cl_sp_from_ham,
    contract_vnl_blocks,
)

# Reuse testing helpers relocated into OCL_Hamiltonian
compare_matrices = _compare_matrices_shared
_orbital_layout = _orbital_layout_shared
_blocked_to_dense = _blocked_to_dense_shared

np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)


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
    with open("param.dat", "w") as f:
        f.write("&OPTION\n")
        f.write(" itheory = 1\n")
        f.write("&END\n")
    
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    fc.setVerbosity(1, 1)

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
    Qneutral_sh_full = fc.get_Qneutral_shell(dims)  # [nsh_max, nspecies_fdata]
    Qneutral_sh = Qneutral_sh_full[:, :dims.nspecies]
    nzx_full = np.array(sd.nzx, dtype=np.int32)
    nzx = nzx_full[:dims.nspecies]
    iatyp = np.array(sd.iatyp, dtype=np.int32)
    ispec_of_atom = np.zeros(dims.natoms, dtype=np.int32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec_of_atom[ia] = int(w[0])
    
    nssh_species = np.ascontiguousarray(sd.nssh[:dims.nspecies], dtype=np.int32)

    # Neighbor lists from Fortran export: neigh_j is 1-based atom index, as used elsewhere in this file.
    # Build both: with self (Fortran diag Vca uses self neighbor) and without self (off-diagonal pairs).
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

    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, sd.lssh[:dims.nspecies], dims.nsh_max, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)
    for ip in range(min(3, pairs.shape[0])):
        i, j = int(pairs[ip, 0]), int(pairs[ip, 1])
        print(f"  rho_avg block pair ({i},{j}) cn_count={int(cn_counts[ip])}")
        print(rho_avg_blocks[ip])

    print("\nTesting Overlap S...")
    dR01 = (atomPos[1] - atomPos[0]).copy()
    dR10 = (atomPos[0] - atomPos[1]).copy()
    S00_f = _scan2c_fortran(1, np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    S11_f = S00_f.copy()
    S01_f = scan2c_fortran(fc, 1, dR01, applyRotation=True)
    S10_f = scan2c_fortran(fc, 1, dR10, applyRotation=True)
    S_f = dense_from_pair_blocks(S00_f, S01_f, S10_f, S11_f)
    S00_o = scan2c_ocl(ham, 'overlap', np.array([0.0, 0.0, 0.0], dtype=np.float64), applyRotation=False)
    S11_o = S00_o.copy()
    S01_o = scan2c_ocl(ham, 'overlap', dR01, applyRotation=True)
    S10_o = scan2c_ocl(ham, 'overlap', dR10, applyRotation=True)
    S_o = dense_from_pair_blocks(S00_o, S01_o, S10_o, S11_o)
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
    cl0 = cl_sp_from_ham(ham, atomTypes_Z, 0)
    cl1 = cl_sp_from_ham(ham, atomTypes_Z, 1)
    cls = [cl0, cl1]

    # Fortran sVNL blocks for all (phi_atom, pp_atom) pairs
    # Use applyRotation=True for off-diagonal geometry; for self (i,i) use dR=0 and applyRotation=False.
    s_map = {}
    for i in (0, 1):
        for k in (0, 1):
            if i == k:
                dR = np.array([0.0, 0.0, 0.0], dtype=np.float64)
                s_map[(i, k)] = scan2c_fortran(fc, 5, dR, in3=1, applyRotation=False)
            else:
                dR = (atomPos[k] - atomPos[i]).copy()
                s_map[(i, k)] = scan2c_fortran(fc, 5, dR, in3=1, applyRotation=True)

    # Full VNL reference per (i,j): sum over k
    Vnl_ref_blocks = []
    for (i, j) in neighbors_vnl:
        acc = np.zeros((4, 4), dtype=np.float64)
        for k in (0, 1):
            acc += contract_vnl_blocks(s_map[(i, k)], s_map[(j, k)], cls[k])
        Vnl_ref_blocks.append(acc)
    Vnl_ref_blocks = np.array(Vnl_ref_blocks, dtype=np.float64)

    # Convert blocks to dense 8x8 for comparison
    Vnl_ref = dense_from_neighbor_list(neighbors_vnl, Vnl_ref_blocks, np.array([4, 4], dtype=np.int32), np.array([0, 4, 8], dtype=np.int32))
    print("  Reference VNL (scan+contract):")
    print(Vnl_ref)

    # OpenCL: assemble_full returns blocks for the neighbor list in the same (nu,mu) convention
    H_vnl_cpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=False)
    H_vnl_gpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=True)

    Vnl_cpu = dense_from_neighbor_list(neighbors_vnl, H_vnl_cpu_blocks, np.array([4, 4], dtype=np.int32), np.array([0, 4, 8], dtype=np.int32))
    Vnl_gpu = dense_from_neighbor_list(neighbors_vnl, H_vnl_gpu_blocks, np.array([4, 4], dtype=np.int32), np.array([0, 4, 8], dtype=np.int32))

    res_Vnl_cpu = compare_matrices("VNL (CPU contraction)", Vnl_ref, Vnl_cpu)
    res_Vnl_gpu = compare_matrices("VNL (GPU contraction)", Vnl_ref, Vnl_gpu)

    print("\nTesting contracted VNL (Fortran gated export vs OpenCL CPU/GPU)...")
    # export_mode=1 reconstructs h_mat like buildh.f90 (does NOT include VNL)
    # export_mode>=2 provides explicit H+VNL export; with T/Vna/Vxc/Vca off we get VNL-only reference.
    Vnl_fortran, _neighs_vnl_f, _sd_vnl_f = firecore_sparse_H_with_options(fc, 0, 0, 0, 1, 0, 0, 0, export_mode=2)
    res_Vnl_cpu_fortran = compare_matrices("VNL (CPU) vs Fortran export", Vnl_fortran, Vnl_cpu)
    res_Vnl_gpu_fortran = compare_matrices("VNL (GPU) vs Fortran export", Vnl_fortran, Vnl_gpu)

    print("\nTesting Fortran export reconstruction (raw full H vs reconstructed full H)...")
    H_full_raw, S_full_raw, neighbors, sparse_data = firecore_sparse_to_dense(fc, export_mode=0)
    H_full_rec, S_full_rec, _neighbors2, _sparse2 = firecore_sparse_to_dense(fc, export_mode=1)
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
    # Isolate Vca from Fortran (Ewald DISABLED)
    H_vca_f, _, _ = get_fortran_sparse_H_with_options(0, 0, 0, 0, 0, 1, 0, ioff_Ewald=0, export_mode=2)
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
    EQ2 = 14.39975
    # Diagonal: <i|V(j)|i> from neighbor charges (interaction 4)
    for i in range(2):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        in_i = int(ispec_of_atom[i]) + 1
        V_ii = np.zeros((ni, ni))
        for k in neigh_lists_self[i]:
            in_k = int(ispec_of_atom[k]) + 1
            dR_ik = atomPos[k] - atomPos[i]
            for isorp in range(1, int(nssh_species[in_k-1]) + 1):
                dqk = dq_shell[isorp-1, k]
                m = fc.scanHamPiece2c(4, isorp, in_i, in_k, in_i, dR_ik, applyRotation=True)
                V_ii += m * dqk
        Vca_manual[i0:i0+ni, i0:i0+ni] += V_ii * EQ2

    # Off-diagonal: ontopl/ontopr (interaction 2/3) for neighbor pairs
    for i in range(2):
        in_i = int(ispec_of_atom[i]) + 1
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        for j in neigh_lists[i]:
            in_j = int(ispec_of_atom[j]) + 1
            nj = int(n_orb_atom[j]); j0 = int(offs[j])
            dR = atomPos[j] - atomPos[i]
            mat2 = np.zeros((ni, nj))
            for isorp in range(1, int(nssh_species[in_i-1]) + 1):
                mat2 += fc.scanHamPiece2c(2, isorp, in_i, in_i, in_j, dR, applyRotation=True) * dq_shell[isorp-1, i]
            mat3 = np.zeros((ni, nj))
            for isorp in range(1, int(nssh_species[in_j-1]) + 1):
                mat3 += fc.scanHamPiece2c(3, isorp, in_i, in_j, in_j, dR, applyRotation=True) * dq_shell[isorp-1, j]
            Vca_manual[i0:i0+ni, j0:j0+nj] += (mat2 + mat3) * EQ2

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
    # Prepare neighbor list and types for assemble_vca (use Fortran neighbor list including self)
    neigh_list_vca = []
    for ia in range(dims.natoms):
        for ja in neigh_lists_self[ia]:
            neigh_list_vca.append((ia, int(ja))) # 0-based from get_neighbors
    
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

    if len(pair_types_vca) > 0 and hasattr(ham, 'vca_map'):
        t0 = int(pair_types_vca[0])
        print(f"  DEBUG vca_map for pair_type={t0}:\n{ham.vca_map[t0]}")
        
    off_vca, diag_vca = ham.assemble_vca(atomPos, np.array(neigh_list_vca), np.array(pair_types_vca), dq_shell)
    
    # Reconstruct dense
    Vca_ocl = np.zeros_like(H_vca_f)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0=int(offs[i]); j0=int(offs[j])
        ni=int(n_orb_atom[i]); nj=int(n_orb_atom[j])

        # Off-diagonal block (V_ij) only for i!=j (Fortran skips ontop for self)
        if i != j:
            Vca_ocl[i0:i0+ni, j0:j0+nj] += off_vca[k].T

        # Diagonal update (V_ii) from neighbor j (includes self neighbor)
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
    print(f"Vxc:         {'PASSED' if res_Vxc else 'SKIPPED'}")
    print(f"Vxc_1c:      {'PASSED' if res_Vxc_1c else 'SKIPPED'}")
    print(f"Vca:         {'PASSED' if res_Vca else 'FAILED (charge export issue)'}")
    print(f"Vxc_ca:      {'PASSED' if res_Vxc_ca else 'SKIPPED'}")
    print(f"H2c (T+Vna): {'PASSED' if res_H2c else 'FAILED'}")
    print(f"H raw==recon:{'PASSED' if res_H_recon else 'FAILED'}")
    print(f"VNL vs F:    {'PASSED' if (res_Vnl_cpu_fortran and res_Vnl_gpu_fortran) else 'FAILED'}")
    print(f"Full H:      {'PASSED' if res_H_full else 'SKIPPED'}")
    print("=" * 40)


if __name__ == "__main__":
    run_verification()
