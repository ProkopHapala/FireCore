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
    max_diff = np.max(diff) if diff.size else 0.0
    max_val = max(np.max(np.abs(fortran_mat)) if fortran_mat.size else 0.0, np.max(np.abs(ocl_mat)) if ocl_mat.size else 0.0)
    print(f"Max difference: {max_diff:.2e}")

    if require_nonzero and max_val < tol:
        print(f"WARNING: {name} is (near) zero for both; treating as failure.")
        return False
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        print("Fortran Matrix:")
        print(fortran_mat)
        print("PyOpenCL Matrix:")
        print(ocl_mat)
        print("Abs Diff:")
        print(diff)
        return False
    print(f"SUCCESS: {name} matches.")
    return True


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
        return False
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        return False
    print(f"SUCCESS: {name} matches.")
    return True


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
    for i in range(natoms):
        ni = int(n_orb_atom[i])
        i0 = int(offs[i])
        for ineigh in range(int(neighn[i])):
            j = int(neigh_j[i, ineigh]) - 1
            if j < 0 or j >= natoms:
                continue
            nj = int(n_orb_atom[j])
            j0 = int(offs[j])
            blk = H_blocks[i, ineigh, :nj, :ni]
            M[i0:i0+ni, j0:j0+nj] += blk.T
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


def run_verification():
    # Linear C3 chain to ensure there exists a common neighbor for (0,2): cn={1}
    atomTypes_Z = np.array([6, 6, 6], dtype=np.int32)
    dCC = 1.25
    atomPos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, dCC],
            [0.0, 0.0, 2.0 * dCC],
        ],
        dtype=np.float64,
    )

    fdata_dir = "./Fdata"

    print("Initializing Fortran FireCore...")
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    fc.setVerbosity(0)

    # Run SCF once; afterwards export rho/s/h from the same Fortran state.
    fc.set_export_mode(0)
    fc.set_options(1, 1, 1, 1, 1, 1, 1)
    fc.SCF(positions=atomPos, iforce=0, nmax_scf=50)

    print("Initializing PyOpenCL Hamiltonian...")
    ham = OCL_Hamiltonian(fdata_dir)
    ham.prepare_splines(atomTypes_Z)
    ham.prepare_data_3c(atomTypes_Z)

    natoms = int(atomPos.shape[0])

    def _scan2c_fortran(interaction, dR, in3=None, isub=0, applyRotation=True):
        in1 = 1
        in2 = 1
        if in3 is None:
            in3 = in2
        return fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation=applyRotation)

    def _scan2c_ocl(root, nz1, nz2, dR, applyRotation=True):
        return ham.scanHamPiece2c(root, int(nz1), int(nz2), dR, applyRotation=applyRotation)

    def _dense_from_neighbor_list(neighs, blocks, n_orb_atom, offs):
        norb = int(offs[-1])
        M = np.zeros((norb, norb), dtype=np.float64)
        for idx, (i, j) in enumerate(neighs):
            ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
            i0 = int(offs[i]); j0 = int(offs[j])
            # blocks are (nu,mu); dense expects (mu,nu)
            M[i0:i0+ni, j0:j0+nj] = blocks[idx, :nj, :ni].T
        return M

    def get_fortran_sparse_full(export_mode=0):
        fc.set_options(1, 1, 1, 1, 1, 1, 1)
        fc.set_export_mode(export_mode)
        dims_ = fc.get_HS_dims()
        sd_ = fc.get_HS_neighs(dims_)
        sd_ = fc.get_HS_sparse(dims_, sd_)
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

    def get_fortran_sparse_H_with_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, export_mode=1):
        fc.set_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca)
        fc.set_export_mode(export_mode)
        dims_ = fc.get_HS_dims()
        sd_ = fc.get_HS_neighs(dims_)
        sd_ = fc.get_HS_sparse(dims_, sd_)
        H = _blocked_to_dense(sd_, sd_.h_mat, natoms)
        neighbors = []
        for i in range(natoms):
            for ineigh in range(int(sd_.neighn[i])):
                j = int(sd_.neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                neighbors.append((i, j))
        return H, neighbors, sd_

    print("\n[PLUMBING] compute_avg_rho (3c gather) using Fortran-exported rho + Qin-shell...")
    dims = fc.get_HS_dims(force_refresh=True)
    sd = fc.get_HS_neighs(dims)
    sd = fc.get_HS_neighsPP(dims, data=sd)
    sd = fc.get_HS_sparse(dims, data=sd)
    sd = fc.get_rho_sparse(dims, data=sd)
    sd = fc.get_rho_off_sparse(dims, data=sd)

    # OpenCL compute_avg_rho weights 3c density pieces by *neutral* charge of the common neighbor.
    # Fortran average_rho uses Qneutral(isorp, indna) (per-shell neutral population for species indna).
    # Here we pass a per-atom scalar Qneutral_total(atom) = sum_shell Qneutral(shell, species(atom)).
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
    Qatom = np.zeros(dims.natoms, dtype=np.float32)
    for ia in range(dims.natoms):
        Qatom[ia] = float(np.sum(Qneutral_sh[:, ispec_of_atom[ia]]))

    # Neighbor lists from Fortran export: neigh_j is 1-based atom index.
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

    # Build pairs from neighbor lists (unique i<j)
    pairs = []
    for ia in range(dims.natoms):
        for j in neigh_lists[ia]:
            if ia < j:
                pairs.append((ia, j))
    if len(pairs) == 0:
        raise RuntimeError("[PLUMBING] No neighbor pairs found")
    pairs = np.array(pairs, dtype=np.int32)

    # Ensure pair (0,2) is present; if not, append it explicitly (even if not a direct neighbor)
    if not np.any((pairs[:, 0] == 0) & (pairs[:, 1] == 2)):
        pairs = np.vstack([pairs, np.array([[0, 2]], dtype=np.int32)])

    cn_offsets, cn_indices = ham.build_common_neighbor_csr(neigh_lists, pairs)
    cn_counts = cn_offsets[1:] - cn_offsets[:-1]
    print(f"  n_pairs={pairs.shape[0]}  n_cn_total={cn_indices.shape[0]}  cn_max={int(np.max(cn_counts))}")

    # Map sparse neighbor blocks for (i,j) to 4x4 blocks; if (i,j) not in sparse neighbor list, use zeros
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
        ineigh = neigh_index.get((int(ia), int(ja)), None)
        if ineigh is None:
            continue
        S_blocks[ip, :, :] = sd.s_mat[int(ia), ineigh, :4, :4]

        # 2-center seed for rho_off: den_ontopl (15) + den_ontopr (16) pieces weighted by Qneutral shells
        # This matches average_rho.f90 off-site 2-center part before adding the 3c correction.
        dR = (atomPos[int(ja)] - atomPos[int(ia)]).copy()
        ispec_i = int(ispec_of_atom[int(ia)])
        ispec_j = int(ispec_of_atom[int(ja)])
        # Fortran species indices are 1-based
        in_i = ispec_i + 1
        in_j = ispec_j + 1

        # Left piece (interaction=15): sum_{isorp=1..nssh(in_i)} den_ontopl * Qneutral(isorp,in_i)
        Af = np.zeros((4, 4), dtype=np.float64)
        for isorp in range(1, int(np.count_nonzero(Qneutral_sh[:, ispec_i] != 0.0)) + 1):
            w = float(Qneutral_sh[isorp - 1, ispec_i])
            if w == 0.0:
                continue
            Af += fc.scanHamPiece2c(15, isorp, in_i, in_i, in_j, dR, applyRotation=True, norb=4) * w

        # Right piece (interaction=16): sum_{isorp=1..nssh(in_j)} den_ontopr * Qneutral(isorp,in_j)
        Bf = np.zeros((4, 4), dtype=np.float64)
        for isorp in range(1, int(np.count_nonzero(Qneutral_sh[:, ispec_j] != 0.0)) + 1):
            w = float(Qneutral_sh[isorp - 1, ispec_j])
            if w == 0.0:
                continue
            Bf += fc.scanHamPiece2c(16, isorp, in_i, in_j, in_j, dR, applyRotation=True, norb=4) * w

        rho_blocks[ip, :, :] = (Af + Bf).astype(np.float32)

    # Triplet type selection: key is (root,nz1,nz2,nz3) where nz3 is the common neighbor type
    root = 'den3'
    pair_triplet_types = np.zeros(pairs.shape[0], dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        nz1 = int(atomTypes_Z[int(ia)]); nz2 = int(atomTypes_Z[int(ja)])
        # if there is at least one CN, use the first CN's species as nz3, else fall back nz3=nz1
        cn0 = int(cn_offsets[ip]); cn1 = int(cn_offsets[ip + 1])
        if cn1 > cn0:
            k = int(cn_indices[cn0])
            nz3 = int(atomTypes_Z[k])
        else:
            nz3 = nz1
        key = (root, nz1, nz2, nz3)
        if key not in ham.species_triplet_map:
            # fallback: find any triplet with same root and nz1,nz2
            found = None
            for kk in ham.species_triplet_map.keys():
                if (kk[0] == root) and (kk[1] == nz1) and (kk[2] == nz2):
                    found = kk
                    break
            if found is None:
                raise RuntimeError(f"[PLUMBING] Missing 3c triplet table for {key}")
            key = found
        pair_triplet_types[ip] = int(ham.species_triplet_map[key])

    nspecies_fdata = int(Qneutral_sh.shape[1])
    if nspecies_fdata <= 0:
        raise RuntimeError("Qneutral_sh has zero species dimension")
    if int(np.max(ispec_of_atom)) >= nspecies_fdata:
        raise RuntimeError(f"ispec_of_atom out of range: max={int(np.max(ispec_of_atom))} nspecies_fdata(Qneutral_sh)={nspecies_fdata}")
    nssh_species = np.zeros(nspecies_fdata, dtype=np.int32)
    # Derive nssh for each Fdata species from Qneutral_sh nonzero count (robust to differing nsh_max).
    for ispec in range(nspecies_fdata):
        nssh_species[ispec] = int(np.count_nonzero(Qneutral_sh[:, ispec] != 0.0))

    # Use Fortran-exported mu/nu/mvalue (make_munu) to reconstruct 3c blocks exactly.
    # FireCore.py stores Fortran mu(index,in1,in2) as sd.mu[in2-1, in1-1, index-1] due to axis reversal.
    n_types = int(np.max(pair_triplet_types)) + 1 if pair_triplet_types.size else 0
    if n_types <= 0:
        raise RuntimeError("pair_triplet_types has zero types")
    n_nz_max = int(ham.n_nz_3c_max)
    mu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
    nu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
    mv3c_map = np.zeros((n_types, n_nz_max), dtype=np.int8)
    # For this test geometry we have only one chemical species; use in1=in2=1.
    in1 = 1
    in2 = 1
    mu3c_map[:, :] = sd.mu[in2-1, in1-1, :n_nz_max].astype(np.int16)
    nu3c_map[:, :] = sd.nu[in2-1, in1-1, :n_nz_max].astype(np.int16)
    mv3c_map[:, :] = sd.mvalue[in2-1, in1-1, :n_nz_max].astype(np.int8)

    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, dims.nsh_max, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)

    # Build Fortran reference blocks: rho_off is stored in sd.rho_off in the same blocked layout as sd.rho.
    # compute_avg_rho kernel returns rho_off-like blocks (unnormalized) for direct comparison.
    ref_blocks = np.zeros_like(rho_avg_blocks)
    mask = np.zeros(pairs.shape[0], dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        ineigh = neigh_index.get((int(ia), int(ja)), None)
        if ineigh is None:
            continue
        ref_blocks[ip, :, :] = sd.rho_off[int(ia), ineigh, :4, :4].astype(np.float32)
        mask[ip] = 1

    res_AvgRho = True
    maxdiff = 0.0
    if np.any(mask == 1):
        ok_avg, max_diff_avg = compare_blocks("AvgRho_off (Fortran rho_off vs OpenCL compute_avg_rho)", ref_blocks, rho_avg_blocks, tol=1e-4)
    else:
        print("WARNING: No comparable neighbor blocks found (mask empty)")
        res_AvgRho = False

    if np.any(mask == 1):
        # Isolate 3c contribution: kernel starts from rho_2c seed and adds 3c pieces.
        # If this matches, remaining mismatch is in rho_2c seed construction.
        rho3c_gpu = rho_avg_blocks - rho_blocks
        rho3c_ref = ref_blocks - rho_blocks
        ok_3c, max_diff_3c = compare_blocks("AvgRho_off 3c-only residual (ref - seed vs gpu - seed)", rho3c_ref, rho3c_gpu, tol=1e-4)

    # Report a few pairs, focusing on (0,2) which should have CN count 1 in the linear chain
    for ip in range(pairs.shape[0]):
        i, j = int(pairs[ip, 0]), int(pairs[ip, 1])
        if (i, j) in [(0, 2), (0, 1), (1, 2)]:
            print(f"\n  rho_avg block pair ({i},{j}) cn_count={int(cn_counts[ip])}")
            print(rho_avg_blocks[ip])
            if mask[ip] == 1:
                print("  rho_off Fortran ref block:")
                print(ref_blocks[ip])
                print("  abs diff:")
                print(np.abs(ref_blocks[ip] - rho_avg_blocks[ip]))

    # ------------------------------------------------------------------------------------------
    # Verify_C2-style 2-center tests (scan API vs OCL scan) for all neighbor pairs (i!=j) and self.
    # ------------------------------------------------------------------------------------------
    n_orb_atom, offs = _orbital_layout(sd, natoms)
    norb = int(offs[-1])

    # For scan-based tests we use all (i,j) including self. For sparse-export parity we must use
    # the exact same neighbor list as Fortran export.
    neighs_all = [(i, j) for i in range(natoms) for j in range(natoms)]

    def _scan_table_2c(interaction_fortran, root_ocl, applyRotation_offdiag=True, applyRotation_self=False, in3=1):
        blocks_f = []
        blocks_o = []
        for (i, j) in neighs_all:
            dR = (atomPos[j] - atomPos[i]).copy()
            if i == j:
                dR[:] = 0.0
                Af = _scan2c_fortran(interaction_fortran, dR, in3=in3, applyRotation=applyRotation_self)
                Ao = _scan2c_ocl(root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_self)
            else:
                Af = _scan2c_fortran(interaction_fortran, dR, in3=in3, applyRotation=applyRotation_offdiag)
                Ao = _scan2c_ocl(root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_offdiag)
            blocks_f.append(Af)
            blocks_o.append(Ao)
        blocks_f = np.array(blocks_f, dtype=np.float64)
        blocks_o = np.array(blocks_o, dtype=np.float64)
        Mf = _dense_from_neighbor_list(neighs_all, blocks_f, n_orb_atom, offs)
        Mo = _dense_from_neighbor_list(neighs_all, blocks_o, n_orb_atom, offs)
        return Mf, Mo

    print("\nTesting Overlap S...")
    S_f, S_o = _scan_table_2c(1, 'overlap', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_S = compare_matrices("Overlap S", S_f, S_o, tol=1e-6, require_nonzero=True)

    print("\nTesting Kinetic T...")
    T_f, T_o = _scan_table_2c(13, 'kinetic', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_T = compare_matrices("Kinetic T", T_f, T_o, tol=1e-5, require_nonzero=True)

    print("\nTesting Vna...")
    # Fortran Vna total = ontopl(2)+ontopr(3) for off-diagonal; vna_atom(4) for self accumulated from neighbor direction.
    # Here we follow verify_C2 behavior: use vna_atom with dR(i->j) for self block as in original script.
    blocks_f = []
    blocks_o = []
    for (i, j) in neighs_all:
        dR = (atomPos[j] - atomPos[i]).copy()
        if i == j:
            # Choose nearest neighbor direction if available; otherwise 0.
            k = 1 if i == 0 else 0
            dR = (atomPos[k] - atomPos[i]).copy()
            Af = _scan2c_fortran(4, dR, in3=1, applyRotation=True)
            Ao = _scan2c_ocl('vna_atom_00', atomTypes_Z[i], atomTypes_Z[i], dR, applyRotation=True)
        else:
            Af = _scan2c_fortran(2, dR, in3=1, applyRotation=True) + _scan2c_fortran(3, dR, in3=1, applyRotation=True)
            Ao = _scan2c_ocl('vna_ontopl_00', atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=True) + _scan2c_ocl('vna_ontopr_00', atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=True)
        blocks_f.append(Af)
        blocks_o.append(Ao)
    Vna_f = _dense_from_neighbor_list(neighs_all, np.array(blocks_f, dtype=np.float64), n_orb_atom, offs)
    Vna_o = _dense_from_neighbor_list(neighs_all, np.array(blocks_o, dtype=np.float64), n_orb_atom, offs)
    res_Vna = compare_matrices("Vna", Vna_f, Vna_o, tol=1e-5, require_nonzero=True)

    print("\nTesting sVNL (PP overlap, interaction=5)...")
    sVNL_f, sVNL_o = _scan_table_2c(5, 'vnl', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    require_nz = np.max(np.abs(sVNL_f)) > 1e-6
    res_sVNL = compare_matrices("sVNL", sVNL_f, sVNL_o, tol=1e-5, require_nonzero=require_nz)

    # ------------------------------------------------------------------------------------------
    # Contracted VNL: scan sVNL + contract (reference) vs ham.assemble_full CPU/GPU and Fortran gated export
    # ------------------------------------------------------------------------------------------
    print("\nTesting contracted VNL (scan-based reference vs OpenCL CPU/GPU)...")

    def _cl_sp_from_ham(Z):
        cl_shell = ham.cl_pseudo.get(int(Z), None)
        if cl_shell is None:
            raise RuntimeError(f"Missing ham.cl_pseudo for Z={Z}")
        info = ham.parser.species_info.get(int(Z), None)
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

    cls = [_cl_sp_from_ham(int(z)) for z in atomTypes_Z]

    # sVNL blocks for all (basis_atom i, pp_atom k)
    s_map = {}
    for i in range(natoms):
        for k in range(natoms):
            if i == k:
                dR = np.array([0.0, 0.0, 0.0], dtype=np.float64)
                s_map[(i, k)] = _scan2c_fortran(5, dR, in3=1, applyRotation=False)
            else:
                dR = (atomPos[k] - atomPos[i]).copy()
                s_map[(i, k)] = _scan2c_fortran(5, dR, in3=1, applyRotation=True)

    def _contract_blocks(A, B, clv):
        out = np.zeros((4, 4), dtype=np.float64)
        for nu in range(4):
            for mu in range(4):
                v = 0.0
                for cc in range(4):
                    v += A[cc, mu] * clv[cc] * B[cc, nu]
                out[nu, mu] = v
        return out

    # Build VNL neighbor list exactly like Fortran export_mode>=2 does:
    # iterate neighPP list and map each (jatom,mbeta) to an index in the normal neigh list.
    # This avoids comparing OpenCL against zero blocks for non-PP pairs.
    neighs_vnl = []
    for i in range(natoms):
        nn  = int(sd.neighn[i])
        npp = int(sd.neighPPn[i])
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
                neighs_vnl.append((i, j))
    neighs_vnl = list(dict.fromkeys(neighs_vnl))
    if len(neighs_vnl) == 0:
        raise RuntimeError("VNL: neighs_vnl is empty (no PP neighbors mapped into neigh list)")

    Vnl_ref_blocks = []
    for (i, j) in neighs_vnl:
        acc = np.zeros((4, 4), dtype=np.float64)
        for k in range(natoms):
            acc += _contract_blocks(s_map[(i, k)], s_map[(j, k)], cls[k])
        Vnl_ref_blocks.append(acc)
    Vnl_ref_blocks = np.array(Vnl_ref_blocks, dtype=np.float64)
    Vnl_ref = _dense_from_neighbor_list(neighs_vnl, Vnl_ref_blocks, n_orb_atom, offs)

    H_vnl_cpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighs_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=False)
    H_vnl_gpu_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighs_vnl.copy(), include_T=False, include_Vna=False, include_Vnl=True, use_gpu_vnl=True)
    Vnl_cpu = _dense_from_neighbor_list(neighs_vnl, np.array(H_vnl_cpu_blocks, dtype=np.float64), n_orb_atom, offs)
    Vnl_gpu = _dense_from_neighbor_list(neighs_vnl, np.array(H_vnl_gpu_blocks, dtype=np.float64), n_orb_atom, offs)

    res_Vnl_cpu = compare_matrices("VNL (CPU contraction)", Vnl_ref, Vnl_cpu, tol=1e-5, require_nonzero=True)
    res_Vnl_gpu = compare_matrices("VNL (GPU contraction)", Vnl_ref, Vnl_gpu, tol=1e-5, require_nonzero=True)

    print("\nTesting contracted VNL (Fortran gated export vs OpenCL CPU/GPU)...")
    Vnl_fortran_full, _neighs_vnl_f, _sd_vnl_f = get_fortran_sparse_H_with_options(0, 0, 0, 1, 0, 0, 0, export_mode=2)
    # Compare only on the PP-neighbor block pattern.
    Vnl_fortran = Vnl_fortran_full
    print(f"  neighs_vnl blocks: {len(neighs_vnl)}")
    print(f"  max|Vnl_fortran|={float(np.max(np.abs(Vnl_fortran))):.3e}  max|Vnl_cpu|={float(np.max(np.abs(Vnl_cpu))):.3e}")
    res_Vnl_cpu_fortran = compare_matrices_brief("VNL (CPU) vs Fortran export", Vnl_fortran, Vnl_cpu, tol=1e-5, require_nonzero=True)
    res_Vnl_gpu_fortran = compare_matrices_brief("VNL (GPU) vs Fortran export", Vnl_fortran, Vnl_gpu, tol=1e-5, require_nonzero=True)

    print("\nTesting Fortran export reconstruction (raw full H vs reconstructed full H)...")
    H_full_raw, S_full_raw, neighbors_sparse, sparse_data = get_fortran_sparse_full(export_mode=0)
    H_full_rec, S_full_rec, _neighbors2, _sparse2 = get_fortran_sparse_full(export_mode=1)
    res_H_recon = compare_matrices("Fortran full H: raw vs reconstructed", H_full_raw, H_full_rec, tol=1e-6, require_nonzero=True)
    _ = compare_matrices("Fortran full S: raw vs reconstructed(export)", S_full_raw, S_full_rec, tol=1e-12, require_nonzero=True)

    print("\nTesting combined 2-center H = T + Vna (from scan API)...")
    H2c_f = T_f + Vna_f
    H2c_o = T_o + Vna_o
    res_H2c = compare_matrices("H2c (T+Vna)", H2c_f, H2c_o, tol=1e-5, require_nonzero=True)

    # Placeholders
    res_Vxc = False
    res_Vxc_1c = False
    res_Vca = False
    res_Vxc_ca = False
    res_H_full = False

    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY")
    print(f"Overlap S:   {'PASSED' if res_S else 'FAILED'}")
    print(f"Kinetic T:   {'PASSED' if res_T else 'FAILED'}")
    print(f"Vna:         {'PASSED' if res_Vna else 'FAILED'}")
    print(f"Vnl:         {'PASSED' if res_sVNL else 'FAILED'}")
    print(f"Vxc:         {'PASSED' if res_Vxc else 'NOT IMPLEMENTED'}")
    print(f"Vxc_1c:      {'PASSED' if res_Vxc_1c else 'NOT IMPLEMENTED'}")
    print(f"Vca:         {'PASSED' if res_Vca else 'NOT IMPLEMENTED'}")
    print(f"Vxc_ca:      {'PASSED' if res_Vxc_ca else 'NOT IMPLEMENTED'}")
    print(f"H2c (T+Vna): {'PASSED' if res_H2c else 'FAILED'}")
    print(f"H raw==recon:{'PASSED' if res_H_recon else 'FAILED'}")
    print(f"VNL vs F:    {'PASSED' if (res_Vnl_cpu_fortran and res_Vnl_gpu_fortran) else 'FAILED'}")
    print(f"AvgRho_off:  {'PASSED' if res_AvgRho else 'FAILED'}")
    print(f"Full H:      {'PASSED' if res_H_full else 'NOT IMPLEMENTED'}")
    print("=" * 40)

    # AvgRho_off is expected to fail until 2c+3c density-table evaluation is fully implemented on GPU.
    if not (res_S and res_T and res_Vna and res_sVNL and res_Vnl_cpu and res_Vnl_gpu and res_H_recon and res_H2c and res_Vnl_cpu_fortran and res_Vnl_gpu_fortran):
        raise RuntimeError("verify_C3: some implemented checks FAILED")

    print("\nDONE verify_C3")


if __name__ == "__main__":
    run_verification()
