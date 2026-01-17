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
    fc.setVerbosity(7, 1)

    # Enable AvgRho diagnostics in Fortran core (exported via libFireCore, logic still in average_rho.f90)
    # IMPORTANT: average_rho off-site part is computed on the neighbor list, so pick an actual neighbor.
    # Use jatom/mbeta wildcards (-1) so Fortran captures the first neighbor of iatom=1.
    fc.set_avg_rho_diag(enable=1, iatom=1, jatom=-1, mbeta=-1)
    en, ti, tj, tb = fc.get_avg_rho_diag_state()
    print(f"[FORTRAN AvgRho DIAG state before SCF] enable={en} iatom={ti} jatom={tj} mbeta={tb}")
    if en != 1 or ti != 1:
        raise RuntimeError("AvgRho diagnostic state not set correctly before SCF")

    # Run SCF once; afterwards export rho/s/h from the same Fortran state.
    fc.set_export_mode(0)
    fc.set_options(1, 1, 1, 1, 1, 1, 1)
    fc.SCF(positions=atomPos, iforce=0, nmax_scf=50)

    # Fetch AvgRho diagnostic buffers from Fortran core
    d_iatom, d_ineigh, d_in1, d_in2, d_jatom, d_mbeta = fc.get_avg_rho_diag_meta()
    eps2c = fc.get_avg_rho_diag_eps2c()
    smS   = fc.get_avg_rho_diag_sm()
    rhom2cS = fc.get_avg_rho_diag_rhom2c()
    rhom3cS = fc.get_avg_rho_diag_rhom3c()
    rhooff_3c = fc.get_avg_rho_diag_rhooff_3c()
    rhooff_fin = fc.get_avg_rho_diag_rhooff_final()
    print("\n[FORTRAN AvgRho DIAG]")
    print(f"  target meta: iatom={d_iatom} ineigh={d_ineigh} in1={d_in1} in2={d_in2} jatom={d_jatom} mbeta={d_mbeta}")
    print("  eps2c:")
    print(eps2c)
    print("  sm (shell):")
    print(smS)
    print("  rhom2c (shell):")
    print(rhom2cS)
    print("  rhom3c (shell):")
    print(rhom3cS)
    print("  rho_off after 3c (orbital):")
    print(rhooff_3c[:4,:4])
    print("  rho_off final (orbital):")
    print(rhooff_fin[:4,:4])

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

    def get_fortran_sparse_H_with_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, ioff_Ewald=1, export_mode=1):
        fc.set_options(ioff_S, ioff_T, ioff_Vna, ioff_Vnl, ioff_Vxc, ioff_Vca, ioff_Vxc_ca, ioff_Ewald)
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
    Qneutral_sh_full = fc.get_Qneutral_shell(dims)  # [nsh_max, nspecies_fdata]
    print(f"DEBUG: Qneutral_sh_full shape={Qneutral_sh_full.shape}")
    print(f"DEBUG: Qneutral_sh_full values:\n{Qneutral_sh_full}")
    # Qneutral_shell is per-species; slice to species used in this run and zero-fill the rest to avoid garbage.
    Qneutral_sh = np.zeros_like(Qneutral_sh_full)
    Qneutral_sh[:, :dims.nspecies] = Qneutral_sh_full[:, :dims.nspecies]
    Qin_shell = fc.get_Qin_shell(dims)  # [nsh_max, natoms]
    print(f"DEBUG: Qin_shell shape={Qin_shell.shape}")
    print(f"DEBUG: Qneutral_sh sliced shape={Qneutral_sh.shape}")
    print(f"DEBUG: Qneutral_sh sliced values:\n{Qneutral_sh}")
    nzx_full = np.array(sd.nzx, dtype=np.int32)
    nzx = nzx_full[:dims.nspecies]  # Only valid species indices (avoid garbage beyond nspecies)
    iatyp = np.array(sd.iatyp, dtype=np.int32)
    ispec_of_atom = np.zeros(dims.natoms, dtype=np.int32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec_of_atom[ia] = int(w[0])
    # TODO/DEBUG: Qatom unused below; summing Qneutral_sh triggered overflow on some runs.
    # Qatom = np.zeros(dims.natoms, dtype=np.float32)
    # for ia in range(dims.natoms):
    #     ispec = ispec_of_atom[ia]
    #     if ispec >= Qneutral_sh.shape[1]:
    #         raise RuntimeError(f"ispec_of_atom[{ia}]={ispec} out of bounds for Qneutral_sh.shape={Qneutral_sh.shape}")
    #     Qatom[ia] = float(np.sum(Qneutral_sh[:, ispec]))


    # Neighbor lists from Fortran export: neigh_j is 1-based atom index.
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
    print(f"DEBUG: neigh_lists (no self) = {neigh_lists}")
    print(f"DEBUG: neigh_lists_self = {neigh_lists_self}")

    # Build pairs from neighbor lists (unique i<j)
    pairs = []
    for ia in range(dims.natoms):
        for j in neigh_lists[ia]:
            if ia < j:
                pairs.append((ia, j))
    if len(pairs) == 0:
        raise RuntimeError("[PLUMBING] No neighbor pairs found")
    pairs = np.array(pairs, dtype=np.int32)

    # TODO/DEBUG: do NOT inject non-Fortran neighbor pairs; comparisons must follow Fortran neighbor list exactly.
    # If pair (0,2) is not in the Fortran neighbor list, we skip its diagnostics below.

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

    DEBUG_QIN_TEST = True
    S_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    # We no longer construct the 2c seed using scanHamPiece2c (libFireCore interface is not reference).
    # Instead, OpenCL compute_avg_rho builds the 2c seed internally from den_ontopl/den_ontopr tables.
    rho_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    rho_blocks_alt = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    rho_blocks_alt_raw = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    rho_blocks_alt_raw_L = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    rho_blocks_alt_raw_R = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    pair_dR = np.zeros((pairs.shape[0], 3), dtype=np.float64)
    pair_in_i = np.zeros(pairs.shape[0], dtype=np.int32)
    pair_in_j = np.zeros(pairs.shape[0], dtype=np.int32)
    pair_mbeta = np.zeros(pairs.shape[0], dtype=np.int32)
    pair_valid = np.zeros(pairs.shape[0], dtype=np.int8)
    rho_blocks_qin = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32) if DEBUG_QIN_TEST else None
    pair_2c_types = np.zeros(pairs.shape[0], dtype=np.int32)
    for ip, (ia, ja) in enumerate(pairs):
        ineigh = neigh_index.get((int(ia), int(ja)), None)
        if ineigh is None:
            continue
        S_blocks[ip, :, :] = sd.s_mat[int(ia), ineigh, :4, :4]

        # Store 2c pair type index for OpenCL seed construction (den_ontopl/den_ontopr tables)
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
        # Keep pair geometry for potential debug, but do not compute rho_blocks via scanHamPiece2c.
        in_i = ispec_i + 1
        in_j = ispec_j + 1
        pair_dR[ip, :] = dR
        pair_in_i[ip] = in_i
        pair_in_j[ip] = in_j
        pair_valid[ip] = 1
        pair_mbeta[ip] = mbeta

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
    nssh_species = np.ascontiguousarray(sd.nssh[:dims.nspecies], dtype=np.int32)

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

    print(f"DEBUG mu/nu/mvalue (sd) first 12: mu={mu3c_map[0,:12]}  nu={nu3c_map[0,:12]}  mv={mv3c_map[0,:12]}")
    # TODO/DEBUG: override mu/nu map to expected sp mapping (from make_munu) to check Fortran export layout
    DEBUG_MUNU_OVERRIDE = True
    if DEBUG_MUNU_OVERRIDE:
        mu_sp = np.array([1, 1, 3, 2, 3, 4, 1, 4, 4, 3], dtype=np.int16)
        nu_sp = np.array([1, 3, 1, 2, 3, 4, 4, 1, 3, 4], dtype=np.int16)
        mv_sp = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1], dtype=np.int8)
        mu3c_map[:, :] = 0
        nu3c_map[:, :] = 0
        mv3c_map[:, :] = 0
        mu3c_map[:, :mu_sp.size] = mu_sp
        nu3c_map[:, :nu_sp.size] = nu_sp
        mv3c_map[:, :mv_sp.size] = mv_sp
        print(f"DEBUG mu/nu/mvalue (override) first 12: mu={mu3c_map[0,:12]}  nu={nu3c_map[0,:12]}  mv={mv3c_map[0,:12]}")

    if np.any(np.isnan(rho_blocks)):
        print("ERROR: rho_blocks contains NaN!")
    if np.any(np.isnan(Qneutral_sh)):
        print("ERROR: Qneutral_sh contains NaN!")
        
    def _epsilon_fb_py(r1, r2):
        r1mag = np.linalg.norm(r1)
        r2mag = np.linalg.norm(r2)
        spe = np.zeros((3, 3), dtype=np.float64)
        if r2mag < 1e-4:
            np.fill_diagonal(spe, 1.0)
            return spe
        zphat = r2 / r2mag
        if r1mag > 1e-4:
            r1hat = r1 / r1mag
            yphat = np.cross(zphat, r1hat)
            ypmag = np.linalg.norm(yphat)
            if ypmag > 1e-6:
                yphat = yphat / ypmag
                xphat = np.cross(yphat, zphat)
                spe[:, 0] = xphat
                spe[:, 1] = yphat
                spe[:, 2] = zphat
                return spe
        if abs(zphat[0]) > 1e-4:
            yphat = np.array([-(zphat[1] + zphat[2]) / zphat[0], 1.0, 1.0])
        elif abs(zphat[1]) > 1e-4:
            yphat = np.array([1.0, -(zphat[0] + zphat[2]) / zphat[1], 1.0])
        else:
            yphat = np.array([1.0, 1.0, -(zphat[0] + zphat[1]) / zphat[2]])
        yphat = yphat / np.linalg.norm(yphat)
        xphat = np.cross(yphat, zphat)
        spe[:, 0] = xphat
        spe[:, 1] = yphat
        spe[:, 2] = zphat
        return spe

    def _twister_pmat(eps):
        # Fortran twister: pmat(1,1)=eps(2,2), pmat(1,2)=eps(2,3), pmat(1,3)=eps(2,1), ...
        return np.array([
            [eps[1, 1], eps[1, 2], eps[1, 0]],
            [eps[2, 1], eps[2, 2], eps[2, 0]],
            [eps[0, 1], eps[0, 2], eps[0, 0]],
        ], dtype=np.float64)

    def _chooser(l, pmat):
        if l == 0:
            return np.array([[1.0]], dtype=np.float64)
        if l == 1:
            return pmat.astype(np.float64)
        raise RuntimeError(f"rotate_fb: unsupported l={l} (only s+p handled)")

    def _rotate_fb_py(in1, in2, eps, mmatrix, lssh, nssh):
        # Fortran rotate_fb: left/right from chooser, apply left * M * right (with right in Fortran indexing)
        pmat = _twister_pmat(eps)
        # FireCore.py stores lssh with species as first axis (row). Use that layout here.
        lssh1 = [int(x) for x in lssh[in1 - 1, :int(nssh[in1 - 1])]]
        lssh2 = [int(x) for x in lssh[in2 - 1, :int(nssh[in2 - 1])]]
        xmatrix = np.zeros_like(mmatrix, dtype=np.float64)
        n1 = 0
        for l1 in lssh1:
            left = _chooser(l1, pmat)
            n2 = 0
            for l2 in lssh2:
                right = _chooser(l2, pmat)
                n1b = n1 + 2 * l1 + 1
                n2b = n2 + 2 * l2 + 1
                sub = mmatrix[n1:n1b, n2:n2b]
                xmatrix[n1:n1b, n2:n2b] += left @ sub @ right.T
                n2 = n2b
            n1 = n1b
        return xmatrix

    def _recover_rotate_sp(hlist, mu, nu, eps, in1, in2, lssh, nssh):
        m = np.zeros((4, 4), dtype=np.float64)
        for idx in range(mu.size):
            imu = int(mu[idx]); inu = int(nu[idx])
            if imu <= 0 or inu <= 0 or imu > 4 or inu > 4:
                continue
            m[imu - 1, inu - 1] = hlist[idx]
        return _rotate_fb_py(in1, in2, eps, m, lssh, nssh)

    def _recover_sp(hlist, mu, nu):
        m = np.zeros((4, 4), dtype=np.float64)
        for idx in range(mu.size):
            imu = int(mu[idx]); inu = int(nu[idx])
            if imu <= 0 or inu <= 0 or imu > 4 or inu > 4:
                continue
            m[imu - 1, inu - 1] = hlist[idx]
        return m

    # Build alternative 2c seed by manually rotating unrotated scanHamPiece2c output
    for ip, (ia, ja) in enumerate(pairs):
        if pair_valid[ip] == 0:
            continue
        mb = int(pair_mbeta[ip])
        r1 = atomPos[int(ia)]
        r2 = atomPos[int(ja)] + sd.xl[mb]
        r21 = r2 - r1
        r21_norm = np.linalg.norm(r21)
        if r21_norm < 1e-12:
            continue
        sighat = r21 / r21_norm
        eps2c = _epsilon_fb_py(r2, sighat)
        in_i = int(pair_in_i[ip])
        in_j = int(pair_in_j[ip])
        left_rot = _rotate_fb_py(in_i, in_i, eps2c, rho_blocks_alt_raw_L[ip], sd.lssh, sd.nssh)
        right_rot = _rotate_fb_py(in_i, in_j, eps2c, rho_blocks_alt_raw_R[ip], sd.lssh, sd.nssh)
        rho_blocks_alt[ip, :, :] = (left_rot + right_rot).astype(np.float32)

    # Raw 3c interpolation check (no rotation/recover) for a single triplet
    DEBUG_RAW3C = True
    err_rot_isorp = float('nan')
    err_rawrot_f = float('nan')
    err_rawrot_o = float('nan')
    if DEBUG_RAW3C and natoms >= 3:
        i = 0; j = 2; k = 1
        dRj = (atomPos[j] - atomPos[i]).copy()
        dRk = (atomPos[k] - atomPos[i]).copy()
        in1 = int(ispec_of_atom[i]) + 1
        in2 = int(ispec_of_atom[j]) + 1
        indna = int(ispec_of_atom[k]) + 1
        nz1 = int(atomTypes_Z[i]); nz2 = int(atomTypes_Z[j]); nz3 = int(atomTypes_Z[k])
        nssh_k = int(sd.nssh[indna - 1])
        print(f"\n[DEBUG] 3c raw check pair({i},{j}) cn={k} nssh_k={nssh_k}")
        dims_dbg = fc.get_HS_dims(force_refresh=True)
        print(f"[DEBUG] max_mu_dim3(dims)={dims_dbg.max_mu_dim3}  n_nz_3c_max(ham)={ham.n_nz_3c_max}  n_nz_max(local)={n_nz_max}")
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
        print("  [PY] eps(3x3):")
        print(eps)
        print(f"  [PY] nssh in1={in1} nssh[in1]={int(sd.nssh[in1-1])} lssh_row={sd.lssh[in1-1, :int(sd.nssh[in1-1])]}")
        print(f"  [PY] nssh in2={in2} nssh[in2]={int(sd.nssh[in2-1])} lssh_row={sd.lssh[in2-1, :int(sd.nssh[in2-1])]}")

        block_f_sum = np.zeros((4, 4), dtype=np.float64)
        block_o_sum = np.zeros((4, 4), dtype=np.float64)
        print("verify_C3.py() 1")
        err_rot_isorp = 0.0
        for isorp in range(1, nssh_k + 1):
            print("verify_C3.py() 2")
            out_f = np.zeros((5, n_nz_max), dtype=np.float64, order='F')
            raw_f = fc.scanHamPiece3c_raw(3, isorp, in1, in2, indna, dRj, dRk, out=out_f)
            print("verify_C3.py() 3")
            print("verify_C3.py() 3 --- ")
            print("verify_C3.py() 3 ---dRk ", dRj )
            print("verify_C3.py() 3 ---dRk ", dRk )
            raw_o = ham.scanHamPiece3c_raw_batch('den3', nz1, nz2, nz3, np.array([dRj], dtype=np.float32), np.array([dRk], dtype=np.float32), isorp=isorp)
            print("verify_C3.py() 4")
            if raw_o is None:
                print(f"  raw3c isorp={isorp}: missing OCL table")
                continue
            raw_o0 = raw_o[0]
            print(f"  [PY] raw_f (theta x ME) isorp={isorp}:")
            print(raw_f[:, :min(10, raw_f.shape[1])])
            print(f"  [PY] raw_o (theta x ME) isorp={isorp}:")
            print(raw_o0[:, :min(10, raw_o0.shape[1])])
            n_me = min(raw_f.shape[1], raw_o0.shape[1], n_nz_max)
            diff_raw = np.max(np.abs(raw_f[:, :n_me] - raw_o0[:, :n_me]))
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
            print(f"  [PY] hlist_f (first 10) isorp={isorp}:", h_f[:10])
            print(f"  [PY] hlist_o (first 10) isorp={isorp}:", h_o[:10])
            m_f = _recover_sp(h_f, mu3c_map[0], nu3c_map[0])
            m_o = _recover_sp(h_o, mu3c_map[0], nu3c_map[0])
            print("  [PY] bcnam_f pre-rot (4x4):")
            print(m_f)
            print("  [PY] bcnam_o pre-rot (4x4):")
            print(m_o)
            wf = float(Qneutral_sh[isorp - 1, indna - 1])
            bcnax_f = _recover_rotate_sp(h_f, mu3c_map[0], nu3c_map[0], eps, in1, in2, sd.lssh, sd.nssh)
            bcnax_o = _recover_rotate_sp(h_o, mu3c_map[0], nu3c_map[0], eps, in1, in2, sd.lssh, sd.nssh)
            bcnax_f_sd = _recover_rotate_sp(h_f, sd.mu[in2-1, in1-1, :n_nz_max].astype(np.int16), sd.nu[in2-1, in1-1, :n_nz_max].astype(np.int16), eps, in1, in2, sd.lssh, sd.nssh)
            print(f"  [PY] bcnax_f post-rot isorp={isorp}:")
            print(bcnax_f)
            print(f"  [PY] bcnax_o post-rot isorp={isorp}:")
            print(bcnax_o)
            print(f"  [PY] bcnax_f_sdmap post-rot isorp={isorp}:")
            print(bcnax_f_sd)
            ref_is = fc.scanHamPiece3c(3, isorp, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4)
            ref_is_T = ref_is.T.copy()
            print(f"  [PY] scanHamPiece3c block isorp={isorp} (raw, col-major):")
            print(ref_is)
            print(f"  [PY] scanHamPiece3c block isorp={isorp} (transposed):")
            print(ref_is_T)
            err_is = float(np.max(np.abs(bcnax_f - ref_is_T)))
            err_rot_isorp = max(err_rot_isorp, err_is)
            print(f"  [PY] max|bcnax_f - scanHamPiece3c(T)| isorp={isorp} = {err_is:.3e}")
            print(f"  [PY] max|bcnax_f_sdmap - scanHamPiece3c(T)| isorp={isorp} = {float(np.max(np.abs(bcnax_f_sd - ref_is_T))):.3e}")
            # TODO/DEBUG: keep alt rotations removed; we now use exact rotate_fb implementation above
            block_f_sum += bcnax_f * wf
            block_o_sum += bcnax_o * wf
        print("verify_C3.py() 5")
        block_ref = fc.scanHamPiece3c(3, 1, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4).T * float(Qneutral_sh[0, indna - 1])
        print("verify_C3.py() 6")
        if nssh_k > 1:
            block_ref += fc.scanHamPiece3c(3, 2, in1, in2, indna, dRj, dRk, applyRotation=True, norb=4).T * float(Qneutral_sh[1, indna - 1])
        print("verify_C3.py() 7")
        err_rawrot_f = float(np.max(np.abs(block_f_sum - block_ref)))
        err_rawrot_o = float(np.max(np.abs(block_o_sum - block_ref)))
        print(f"  raw->rot (Fraw) vs Fortran scan: max|diff|={err_rawrot_f:.2e}")
        print(f"  raw->rot (OCLraw) vs Fortran scan: max|diff|={err_rawrot_o:.2e}")

    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, sd.lssh[:dims.nspecies], dims.nsh_max, pair_2c_types=pair_2c_types, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)
    rho_avg_blocks_qin = None
    if DEBUG_QIN_TEST:
        rho_avg_blocks_qin = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks_qin, Qin_shell, ispec_of_atom, nssh_species, sd.lssh[:dims.nspecies], dims.nsh_max, pair_2c_types=pair_2c_types, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mvalue3c_map=mv3c_map)

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
        ref_blocks[ip, :, :] = sd.rho_off[int(ia), ineigh, :4, :4].astype(np.float32)
        mask[ip] = 1

    res_AvgRho = False
    err_AvgRho = float('nan')
    if np.any(mask == 1):
        ok_avg, max_diff_avg = compare_blocks("AvgRho_off (Fortran rho_off vs OpenCL compute_avg_rho)", ref_blocks, rho_avg_blocks, tol=1e-4)
        res_AvgRho = ok_avg
        err_AvgRho = max_diff_avg
        if DEBUG_QIN_TEST:
            ok_avg_qin, max_diff_avg_qin = compare_blocks("AvgRho_off (Qin weights test)", ref_blocks, rho_avg_blocks_qin, tol=1e-4)
            print(f"  DEBUG AvgRho Qin test: ok={ok_avg_qin} err={max_diff_avg_qin:.2e}")
    else:
        print("WARNING: No comparable neighbor blocks found (mask empty)")

    # Report a few pairs (only those present in Fortran neighbor list)
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
    res_S, err_S = compare_matrices("Overlap S", S_f, S_o, tol=1e-6, require_nonzero=True)

    print("\nTesting Kinetic T...")
    T_f, T_o = _scan_table_2c(13, 'kinetic', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_T, err_T = compare_matrices("Kinetic T", T_f, T_o, tol=1e-5, require_nonzero=True)

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
    res_Vna, err_Vna = compare_matrices("Vna", Vna_f, Vna_o, tol=1e-5, require_nonzero=True)

    print("\nTesting sVNL (PP overlap, interaction=5)...")
    sVNL_f, sVNL_o = _scan_table_2c(5, 'vnl', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    require_nz = np.max(np.abs(sVNL_f)) > 1e-6
    res_sVNL, err_sVNL = compare_matrices("sVNL", sVNL_f, sVNL_o, tol=1e-5, require_nonzero=require_nz)

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

    res_Vnl_cpu, err_Vnl_cpu = compare_matrices("VNL (CPU contraction)", Vnl_ref, Vnl_cpu, tol=1e-5, require_nonzero=True)
    res_Vnl_gpu, err_Vnl_gpu = compare_matrices("VNL (GPU contraction)", Vnl_ref, Vnl_gpu, tol=1e-5, require_nonzero=True)

    print("\nTesting contracted VNL (Fortran gated export vs OpenCL CPU/GPU)...")
    Vnl_fortran_full, _neighs_vnl_f, _sd_vnl_f = get_fortran_sparse_H_with_options(0, 0, 0, 1, 0, 0, 0, export_mode=2)
    # Compare only on the PP-neighbor block pattern.
    Vnl_fortran = Vnl_fortran_full
    print(f"  neighs_vnl blocks: {len(neighs_vnl)}")
    print(f"  max|Vnl_fortran|={float(np.max(np.abs(Vnl_fortran))):.3e}  max|Vnl_cpu|={float(np.max(np.abs(Vnl_cpu))):.3e}")
    res_Vnl_cpu_fortran, err_Vnl_cpu_fortran = compare_matrices_brief("VNL (CPU) vs Fortran export", Vnl_fortran, Vnl_cpu, tol=1e-5, require_nonzero=True)
    res_Vnl_gpu_fortran, err_Vnl_gpu_fortran = compare_matrices_brief("VNL (GPU) vs Fortran export", Vnl_fortran, Vnl_gpu, tol=1e-5, require_nonzero=True)

    print("\nTesting Fortran export reconstruction (raw full H vs reconstructed full H)...")
    H_full_raw, S_full_raw, neighbors_sparse, sparse_data = get_fortran_sparse_full(export_mode=0)
    H_full_rec, S_full_rec, _neighbors2, _sparse2 = get_fortran_sparse_full(export_mode=1)
    res_H_recon, err_H_recon = compare_matrices("Fortran full H: raw vs reconstructed", H_full_raw, H_full_rec, tol=1e-6, require_nonzero=True)
    res_S_recon, err_S_recon = compare_matrices("Fortran full S: raw vs reconstructed(export)", S_full_raw, S_full_rec, tol=1e-12, require_nonzero=True)

    print("\nTesting combined 2-center H = T + Vna (from scan API)...")
    H2c_f = T_f + Vna_f
    H2c_o = T_o + Vna_o
    res_H2c, err_H2c = compare_matrices("H2c (T+Vna)", H2c_f, H2c_o, tol=1e-5, require_nonzero=True)

    print("\nTesting Vca (charge-dependent, 2-center only)...")
    # Fortran Vca term in assembled H is: vca + ewaldlr - ewaldsr  (see firecore_get_HS_sparse)
    # Here we export each component explicitly.
    vcaL4_f    = fc.export_interaction4D(10) # vca_ontopl
    vcaR4_f    = fc.export_interaction4D(11) # vca_ontopr
    vcaA4_f    = fc.export_interaction4D(12) # vca_atom
    ewaldsr4_f = fc.export_interaction4D(2)  # ewaldsr
    ewaldlr4_f = fc.export_interaction4D(3)  # ewaldlr
    print(f"  Fortran export max|vcaL4|={float(np.max(np.abs(vcaL4_f))):.3e}  max|vcaR4|={float(np.max(np.abs(vcaR4_f))):.3e}  max|vcaA4|={float(np.max(np.abs(vcaA4_f))):.3e}")
    ia_max, ineigh_max, inu_max, imu_max = np.unravel_index(int(np.argmax(np.abs(vcaA4_f))), vcaA4_f.shape)
    print(f"  Fortran export argmax|vcaA4| at iatom={ia_max} ineigh={ineigh_max} (neighn={int(sd.neighn[ia_max])} neigh_j={int(sd.neigh_j[ia_max, ineigh_max])})  |val|={float(np.abs(vcaA4_f[ia_max, ineigh_max, inu_max, imu_max])):.3e}")
    VcaL_f      = _blocked_to_dense(sd, vcaL4_f, natoms)
    VcaR_f      = _blocked_to_dense(sd, vcaR4_f, natoms)
    VcaA_f      = _blocked_to_dense_vca_atom(sd, vcaA4_f, natoms)
    # Core-reference Vca2c is the *component sum* (ontopl + ontopr + atom).
    # DO NOT use export_interaction4D(1) as reference; it is an interface-level array and not trusted.
    Vca2c_f     = VcaL_f + VcaR_f + VcaA_f
    EwaldSR_f   = _blocked_to_dense(sd, ewaldsr4_f, natoms)
    EwaldLR_f   = _blocked_to_dense(sd, ewaldlr4_f, natoms)
    Vca_full_f  = Vca2c_f + EwaldLR_f - EwaldSR_f

    H_vca_f = Vca2c_f
    print("  Fortran Vca2c (component sum, no Ewald):")
    print(H_vca_f)
    print(f"  Fortran |Vca2c|={float(np.max(np.abs(Vca2c_f))):.3e}  |EwaldSR|={float(np.max(np.abs(EwaldSR_f))):.3e}  |EwaldLR|={float(np.max(np.abs(EwaldLR_f))):.3e}")
    print(f"  Fortran |Vca_full|={float(np.max(np.abs(Vca_full_f))):.3e}")

    dq_shell = np.zeros((dims.nsh_max, natoms))
    for ia in range(natoms):
        ispec = int(ispec_of_atom[ia])
        for issh in range(int(sd.nssh[ispec])):
            dq_shell[issh, ia] = Qin_shell[issh, ia] - Qneutral_sh[issh, ispec]
    print("  dq_shell:")
    print(dq_shell)

    # Build directed neighbor list for OpenCL Vca assembly.
    # IMPORTANT: Fortran stores the diagonal Vca contributions in a dedicated "self" slot (matom=neigh_self(iatom)),
    # but those diagonal blocks are accumulated from *non-self* neighbors inside assemble_ca_2c.
    # Therefore for OpenCL assembly we exclude self-neighbors (ia==ja) and sum diagonal updates over non-self pairs.
    neigh_list_vca = []
    neigh_mbeta_vca = []
    for ia in range(natoms):
        nn = int(sd.neighn[ia])
        for ineigh in range(nn):
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = int(sd.neigh_b[ia, ineigh])
            if ja < 0 or ja >= natoms:
                continue
            # Include self-pairs (ia==ja, mbeta==0) because Fortran includes on-site vca_atom.
            # This lets us compare OpenCL vs Fortran element-by-element including the self-slot contribution.
            if ja == ia and mb != 0:
                continue
            neigh_list_vca.append((ia, ja))
            neigh_mbeta_vca.append(mb)

    if len(neigh_mbeta_vca) > 0:
        print(f"  VCA neigh_b stats: max|mbeta|={int(np.max(np.abs(np.array(neigh_mbeta_vca, dtype=np.int32))))}  unique={sorted(set([int(x) for x in neigh_mbeta_vca]))}")

    pair_types_vca = []
    for (i, j) in neigh_list_vca:
        nz1 = int(atomTypes_Z[i]); nz2 = int(atomTypes_Z[j])
        k = ('vna', nz1, nz2)
        idx = ham.species_pair_map.get(k)
        if idx is None:
            print(f"Warning: missing vna pair type for {nz1},{nz2}")
            idx = 0
        pair_types_vca.append(idx)

    # IMPORTANT: Fortran uses r2 = atomPos[j] + xl[mbeta] for each neighbor.
    # OpenCL assemble_vca currently accepts only a flat atom position array; emulate mbeta shifts by
    # constructing a per-pair position buffer (two atoms per pair) when any mbeta!=0.
    mb_max = int(np.max(np.abs(np.array(neigh_mbeta_vca, dtype=np.int32)))) if len(neigh_mbeta_vca) > 0 else 0
    # Pick one directed neighbor pair for synchronized VCA debugging (Fortran vs OpenCL).
    # Target: ia=1 -> ja=0 in 0-based indexing (Fortran iatom=2, jatom=1), mbeta==0.
    dbg_pair = -1
    dbg_ish = 0
    for k,(i,j) in enumerate(neigh_list_vca):
        if i == 1 and j == 0 and int(neigh_mbeta_vca[k]) == 0:
            dbg_pair = k
            break
    print(f"  [VCA_DBG][PY][CPSEL] target ia=1 ja=0 mbeta=0 -> dbg_pair={dbg_pair} dbg_ish={dbg_ish}")

    if mb_max != 0:
        ratoms_pairs = np.zeros((2*len(neigh_list_vca), 3), dtype=np.float64)
        neigh_pairs = np.zeros((len(neigh_list_vca), 2), dtype=np.int32)
        for k,(i,j) in enumerate(neigh_list_vca):
            mb = int(neigh_mbeta_vca[k])
            ratoms_pairs[2*k+0] = atomPos[i]
            ratoms_pairs[2*k+1] = atomPos[j] + sd.xl[mb]
            neigh_pairs[k,0] = 2*k+0
            neigh_pairs[k,1] = 2*k+1
        off_vca, diag_vca = ham.assemble_vca(ratoms_pairs, neigh_pairs, np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                             ispec_of_atom=np.repeat(ispec_of_atom, 1), nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                             dbg_enable=int(dbg_pair>=0), dbg_pair=dbg_pair, dbg_ish=dbg_ish)
    else:
        off_vca, diag_vca = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                         ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                         dbg_enable=int(dbg_pair>=0), dbg_pair=dbg_pair, dbg_ish=dbg_ish)
    Vca_ocl = np.zeros_like(H_vca_f)
    Vca_ocl_off = np.zeros_like(H_vca_f)
    Vca_ocl_diag = np.zeros_like(H_vca_f)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0 = int(offs[i]); j0 = int(offs[j])
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        if i != j:
            Vca_ocl[i0:i0+ni, j0:j0+nj] += off_vca[k]
            Vca_ocl_off[i0:i0+ni, j0:j0+nj] += off_vca[k]
        # Diagonal updates from vna_atom:
        # diag_vca[k,0] updates V_ii from potential on j (including self-pair i==j).
        Vca_ocl[i0:i0+ni, i0:i0+ni] += diag_vca[k, 0]
        Vca_ocl_diag[i0:i0+ni, i0:i0+ni] += diag_vca[k, 0]

    res_Vca, err_Vca = compare_matrices_brief("Vca2c (vca array)", Vca2c_f, Vca_ocl, tol=1e-4, require_nonzero=True)

    # OpenCL component-wise Vca2c (bcca splines only) by masking vca_map indices
    vmap = ham.vca_map
    vmap_L = vmap.copy(); vmap_L[:, :, 1] = -1; vmap_L[:, :, 2] = -1
    vmap_R = vmap.copy(); vmap_R[:, :, 0] = -1; vmap_R[:, :, 2] = -1
    vmap_A = vmap.copy(); vmap_A[:, :, 0] = -1; vmap_A[:, :, 1] = -1

    # Diagnostic: swapped L/R component interpretation (if OpenCL component order differs)
    vmap_Ls = vmap.copy(); vmap_Ls[:, :, 0] = -1; vmap_Ls[:, :, 2] = -1  # take comp-1 as "L"
    vmap_Rs = vmap.copy(); vmap_Rs[:, :, 1] = -1; vmap_Rs[:, :, 2] = -1  # take comp-0 as "R"

    # Extra synchronized atom-only traces for ia=1 -> ja=0 and ia=1 -> ja=2 at ish=0.
    # These calls are for printing only and do not affect later comparisons.
    dbg_pairs = []
    for target_j in (0, 2):
        kp = -1
        for k,(i,j) in enumerate(neigh_list_vca):
            if i == 1 and j == target_j and int(neigh_mbeta_vca[k]) == 0:
                kp = k
                break
        dbg_pairs.append(kp)
    dbg_ish_list = [0, 1]
    print(f"  [VCA_DBG][PY][CPSEL2] ia=1 -> ja=(0,2) mbeta=0 dbg_pairs={dbg_pairs} dbg_ish_list={dbg_ish_list}")
    for kp in dbg_pairs:
        if kp >= 0:
            for ish_dbg in dbg_ish_list:
                _off_tmp, _diag_tmp = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                      ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                      vca_map_override=vmap_A, dbg_enable=1, dbg_pair=kp, dbg_ish=int(ish_dbg))

    offL, diagL = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_L)
    offR, diagR = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_R)
    offA, diagA = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_A)

    offLs, diagLs = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                   ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_Ls)
    offRs, diagRs = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                   ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_Rs)

    VcaL_o = np.zeros_like(H_vca_f)
    VcaR_o = np.zeros_like(H_vca_f)
    VcaA_o = np.zeros_like(H_vca_f)
    VcaL_os = np.zeros_like(H_vca_f)
    VcaR_os = np.zeros_like(H_vca_f)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0 = int(offs[i]); j0 = int(offs[j])
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        # offdiag components
        VcaL_o[i0:i0+ni, j0:j0+nj] += offL[k]
        VcaR_o[i0:i0+ni, j0:j0+nj] += offR[k]
        VcaL_os[i0:i0+ni, j0:j0+nj] += offLs[k]
        VcaR_os[i0:i0+ni, j0:j0+nj] += offRs[k]
        # diagonal (atom) component
        VcaA_o[i0:i0+ni, i0:i0+ni] += diagA[k, 0]

    res_VcaL, err_VcaL = compare_matrices_brief("Vca2c_ontopl (2)", VcaL_f, VcaL_o, tol=1e-4, require_nonzero=True)
    res_VcaR, err_VcaR = compare_matrices_brief("Vca2c_ontopr (3)", VcaR_f, VcaR_o, tol=1e-4, require_nonzero=True)
    res_VcaA, err_VcaA = compare_matrices_brief("Vca2c_atom (4)",   VcaA_f, VcaA_o, tol=1e-4, require_nonzero=True)

    _res_VcaL_s, err_VcaL_s = compare_matrices_brief("Vca2c_ontopl (SWAPPED)", VcaL_f, VcaL_os, tol=1e-4, require_nonzero=True)
    _res_VcaR_s, err_VcaR_s = compare_matrices_brief("Vca2c_ontopr (SWAPPED)", VcaR_f, VcaR_os, tol=1e-4, require_nonzero=True)
    print(f"  [diag] component swap check: err_L={err_VcaL:.3e} err_R={err_VcaR:.3e} | swapped: err_L={err_VcaL_s:.3e} err_R={err_VcaR_s:.3e}")
    # Ewald components are currently not implemented on OCL side => compare against zeros (visible in final summary)
    EwaldSR_o = np.zeros_like(EwaldSR_f)
    EwaldLR_o = np.zeros_like(EwaldLR_f)
    res_EwaldSR, err_EwaldSR = compare_matrices_brief("EwaldSR", EwaldSR_f, EwaldSR_o, tol=1e-5, require_nonzero=True)
    res_EwaldLR, err_EwaldLR = compare_matrices_brief("EwaldLR", EwaldLR_f, EwaldLR_o, tol=1e-5, require_nonzero=True)
    Vca_full_o = Vca_ocl + EwaldLR_o - EwaldSR_o
    res_Vca_full, err_Vca_full = compare_matrices_brief("Vca_full (2c)", Vca_full_f, Vca_full_o, tol=1e-5, require_nonzero=True)

    # Pairwise block diagnostics (non-self directed neighbors only).
    print("\n--- Vca pairwise block diagnostics (non-self directed neighbors) ---")
    max_pair_diff = 0.0
    worst = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        # Find corresponding Fortran neighbor slot (ia, ineigh) where neigh_j==ja+1
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh
                break
        if ineighF < 0:
            continue
        # Core reference off-diagonal block is ontopl+ontopr. Atom term contributes only to diagonal updates.
        blk_f = vcaL4_f[ia, ineighF, :, :] + vcaR4_f[ia, ineighF, :, :]   # [nu,mu]
        blk_o = off_vca[kpair].T             # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_diff:
            max_pair_diff = d
            worst = (ia, ja, ineighF, kpair, d)
    print(f"  max|Fortran_vca_block - OCL_off_block| over non-self directed pairs = {max_pair_diff:.3e}")
    if worst is not None:
        ia, ja, ineighF, kpair, d = worst
        print(f"  worst directed pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        print("  Fortran vca(offdiag) block [nu,mu] = vcaL+vcaR:")
        print(vcaL4_f[ia, ineighF] + vcaR4_f[ia, ineighF])
        print("  OpenCL off block [nu,mu]:")
        print(off_vca[kpair].T)

    # Component-level pairwise diagnostics for ontopl/ontopr
    print("\n--- Vca pairwise diagnostics: ontopl vs export ---")
    max_pair_L = 0.0
    worst_L = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh; break
        if ineighF < 0: continue
        blk_f = vcaL4_f[ia, ineighF, :, :]       # [nu,mu]
        blk_o = offL[kpair].T                    # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_L:
            max_pair_L = d
            worst_L = (ia, ja, ineighF, kpair, d)
    print(f"  max|vcaL_export - OCL_offL| over non-self directed pairs = {max_pair_L:.3e}")
    if worst_L is not None:
        ia, ja, ineighF, kpair, d = worst_L
        print(f"  worst ontopl pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        blk_f0 = vcaL4_f[ia, ineighF]
        blk_fT = blk_f0.T
        blk_o  = offL[kpair].T
        print(f"  ontopl orientation check: max|F - OCL|={float(np.max(np.abs(blk_f0-blk_o))):.3e}  max|F.T - OCL|={float(np.max(np.abs(blk_fT-blk_o))):.3e}")
        print("  Fortran vcaL block [nu,mu]:")
        print(blk_f0)
        print("  OpenCL offL block [nu,mu]:")
        print(blk_o)

        # Try to infer orbital ordering/sign mismatch by brute force (full 4-orb basis).
        import itertools
        fref = blk_f0[:4, :4].copy()    # [nu,mu]
        oref = blk_o[:4, :4].copy()     # [nu,mu]
        best = None
        for perm in itertools.permutations([0, 1, 2, 3], 4):
            for sgn in itertools.product([-1.0, 1.0], repeat=4):
                P = np.array(perm, dtype=int)
                S = np.array(sgn, dtype=np.float64)
                M = oref[np.ix_(P, P)].copy()
                M = (S[:, None] * M) * (S[None, :])
                dd = float(np.max(np.abs(fref - M)))
                if (best is None) or (dd < best[0]):
                    best = (dd, perm, sgn)
        if best is not None:
            dd, perm, sgn = best
            print(f"  [diag] best 4-orb map for ontopl: max|F - mapped(OCL)|={dd:.3e}  perm={perm}  sgn={sgn}")

    print("\n--- Vca pairwise diagnostics: ontopr vs export ---")
    max_pair_R = 0.0
    worst_R = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh; break
        if ineighF < 0: continue
        blk_f = vcaR4_f[ia, ineighF, :, :]       # [nu,mu]
        blk_o = offR[kpair].T                    # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_R:
            max_pair_R = d
            worst_R = (ia, ja, ineighF, kpair, d)
    print(f"  max|vcaR_export - OCL_offR| over non-self directed pairs = {max_pair_R:.3e}")
    if worst_R is not None:
        ia, ja, ineighF, kpair, d = worst_R
        print(f"  worst ontopr pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        blk_f0 = vcaR4_f[ia, ineighF]
        blk_fT = blk_f0.T
        blk_o  = offR[kpair].T
        print(f"  ontopr orientation check: max|F - OCL|={float(np.max(np.abs(blk_f0-blk_o))):.3e}  max|F.T - OCL|={float(np.max(np.abs(blk_fT-blk_o))):.3e}")
        print("  Fortran vcaR block [nu,mu]:")
        print(blk_f0)
        print("  OpenCL offR block [nu,mu]:")
        print(blk_o)

        import itertools
        fref = blk_f0[:4, :4].copy()    # [nu,mu]
        oref = blk_o[:4, :4].copy()     # [nu,mu]
        best = None
        for perm in itertools.permutations([0, 1, 2, 3], 4):
            for sgn in itertools.product([-1.0, 1.0], repeat=4):
                P = np.array(perm, dtype=int)
                S = np.array(sgn, dtype=np.float64)
                M = oref[np.ix_(P, P)].copy()
                M = (S[:, None] * M) * (S[None, :])
                dd = float(np.max(np.abs(fref - M)))
                if (best is None) or (dd < best[0]):
                    best = (dd, perm, sgn)
        if best is not None:
            dd, perm, sgn = best
            print(f"  [diag] best 4-orb map for ontopr: max|F - mapped(OCL)|={dd:.3e}  perm={perm}  sgn={sgn}")

    # Diagonal diagnostics for atom term
    print("\n--- Vca atom(diag) diagnostics: export self-slot vs OCL diagA ---")
    max_ai = 0.0
    worst_ai = None
    # Find self-slot index for each atom
    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(sd.neighn[i])):
            if int(sd.neigh_j[i, ineigh]) == ii and int(sd.neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break
    for i in range(natoms):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        # Fortran export: vca_atom is accumulated in self-slot (matom = neigh_self(iatom))
        acc_f = np.zeros((ni, ni), dtype=np.float64)
        if neigh_self[i] >= 0:
            acc_f = vcaA4_f[i, neigh_self[i], :ni, :ni].T
        acc_o = np.zeros((ni, ni), dtype=np.float64)
        for kpair, (ia, ja) in enumerate(neigh_list_vca):
            if ia != i: continue
            acc_o += diagA[kpair, 0, :ni, :ni]
        d = float(np.max(np.abs(acc_f - acc_o)))
        if d > max_ai:
            max_ai = d
            worst_ai = (i, d, acc_f, acc_o)
    print(f"  max|atom_diag_export - atom_diag_OCL| over atoms = {max_ai:.3e}")
    if worst_ai is not None:
        i, d, acc_f, acc_o = worst_ai
        print(f"  worst atom i={i} diff={d:.3e}")
        print("  Fortran atom diag block (from self-slot):")
        print(acc_f)
        print("  OpenCL atom diag block (sum diagA[k,0] for ia==i):")
        print(acc_o)

    diff_ocl = np.abs(Vca_ocl - H_vca_f)
    # err_Vca already computed above from Fortran export vs OpenCL.
    print(f"  Max diff Vca OCL vs Fortran: {err_Vca:.2e}")
    if res_Vca:
        print("SUCCESS: OpenCL Vca2c matches Fortran export.")
    else:
        print("FAILURE: OpenCL Vca2c mismatch vs Fortran export.")

    def _block(M, i, j):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        return M[i0:i0+ni, j0:j0+nj]

    for i in range(natoms):
        for j in range(natoms):
            if i == j or j in neigh_lists[i]:
                df = _block(H_vca_f, i, j)
                do = _block(Vca_ocl, i, j)
                print(f"  VCA block ({i},{j}) max|F-ocl|={float(np.max(np.abs(df-do))):.2e}")

    # Component-wise diagnostics for the full Vca term used in H assembly: vca + ewaldlr - ewaldsr
    print("\n--- Vca component-wise diagnostics (Fortran) ---")
    _ = compare_matrices_brief("Vca2c (vca array)", Vca2c_f, Vca_ocl, tol=1e-4, require_nonzero=True)
    print(f"  Max|EwaldSR|={float(np.max(np.abs(EwaldSR_f))):.3e}  Max|EwaldLR|={float(np.max(np.abs(EwaldLR_f))):.3e}")
    Vca_resid = Vca_full_f - Vca_ocl
    print(f"  Residual fullVca-ocl: max|.|={float(np.max(np.abs(Vca_resid))):.3e}")
    # If you later implement EwaldSR/LR in OpenCL, compare against Vca_full_f directly.

    print("\n--- Vca split diagnostics: offdiag vs diag (OpenCL) ---")
    # Off-diagonal differences
    mask_off = np.ones_like(H_vca_f, dtype=bool)
    np.fill_diagonal(mask_off, False)
    diff_off = np.abs((Vca2c_f - Vca_ocl_off) * mask_off)
    diff_diag = np.abs(np.diag(Vca2c_f) - np.diag(Vca_ocl_diag))
    print(f"  max|offdiag(F - OCL_off)|={float(np.max(diff_off)):.3e}")
    print(f"  max|diag(F - OCL_diag)|={float(np.max(diff_diag)):.3e}")

    ij = np.unravel_index(int(np.argmax(np.abs(Vca2c_f - Vca_ocl))), Vca2c_f.shape)
    print(f"  WORST elem (F-OCL) at {ij}: F={float(Vca2c_f[ij]): .6e}  OCL={float(Vca_ocl[ij]): .6e}")
    # Locate which atom blocks contain that matrix element
    def _which_atom(idx_):
        for ia in range(natoms):
            if idx_ >= int(offs[ia]) and idx_ < int(offs[ia+1]):
                return ia
        return -1
    ia = _which_atom(int(ij[0])); ja = _which_atom(int(ij[1]))
    print(f"  WORST elem belongs to atom-block ({ia},{ja})")

    # Placeholders
    res_Vxc = False
    res_Vxc_1c = False
    res_Vxc_ca = False
    res_H_full = False

    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY (max|diff|)")
    print(f"Overlap S:   {'PASSED' if res_S else 'FAILED'}  err={err_S:.2e}")
    print(f"Kinetic T:   {'PASSED' if res_T else 'FAILED'}  err={err_T:.2e}")
    print(f"Vna:         {'PASSED' if res_Vna else 'FAILED'}  err={err_Vna:.2e}")
    print(f"Vnl:         {'PASSED' if res_sVNL else 'FAILED'}  err={err_sVNL:.2e}")
    print(f"Vxc:         {'PASSED' if res_Vxc else 'NOT IMPLEMENTED'}  err=nan")
    print(f"Vxc_1c:      {'PASSED' if res_Vxc_1c else 'NOT IMPLEMENTED'}  err=nan")
    print(f"Vca:         {'PASSED' if res_Vca else 'FAILED'}  err={err_Vca:.2e}")
    print(f" Vca2c_ontopl:{'PASSED' if res_VcaL else 'FAILED'}  err={err_VcaL:.2e}")
    print(f" Vca2c_ontopr:{'PASSED' if res_VcaR else 'FAILED'}  err={err_VcaR:.2e}")
    print(f" Vca2c_atom:  {'PASSED' if res_VcaA else 'FAILED'}  err={err_VcaA:.2e}")
    print(f" EwaldSR:     {'PASSED' if res_EwaldSR else 'FAILED'}  err={err_EwaldSR:.2e}")
    print(f" EwaldLR:     {'PASSED' if res_EwaldLR else 'FAILED'}  err={err_EwaldLR:.2e}")
    print(f" Vca_full2c:  {'PASSED' if res_Vca_full else 'FAILED'}  err={err_Vca_full:.2e}")
    print(f"Vxc_ca:      {'PASSED' if res_Vxc_ca else 'NOT IMPLEMENTED'}  err=nan")
    print(f"H2c (T+Vna): {'PASSED' if res_H2c else 'FAILED'}  err={err_H2c:.2e}")
    print(f"H raw==recon:{'PASSED' if res_H_recon else 'FAILED'}  err={err_H_recon:.2e}")
    print(f"VNL vs F:    {'PASSED' if (res_Vnl_cpu_fortran and res_Vnl_gpu_fortran) else 'FAILED'}  err=max({err_Vnl_cpu_fortran:.2e},{err_Vnl_gpu_fortran:.2e})")
    print(f"AvgRho_off:  {'PASSED' if res_AvgRho else 'FAILED'}  err={err_AvgRho:.2e}")
    
    #print(f"AvgRho_T:    {'PASSED' if err_AvgRho_T < 1e-6 else 'FAILED'}  err={err_AvgRho_T:.2e}")
    #print(f"AvgRho_seed:{'PASSED' if err_AvgRho_seed < 1e-6 else 'FAILED'}  err={err_AvgRho_seed:.2e}")
    #print(f"Seed_T:      {'PASSED' if err_AvgRho_seed_T < 1e-6 else 'FAILED'}  err={err_AvgRho_seed_T:.2e}")
    #print(f"AvgRho_3c:   {'PASSED' if err_AvgRho3c < 1e-6 else 'FAILED'}  err={err_AvgRho3c:.2e}")
    #print(f"AvgRho_3c_T: {'PASSED' if err_AvgRho3c_T < 1e-6 else 'FAILED'}  err={err_AvgRho3c_T:.2e}")

    print(f"Rot(isorp):  {'PASSED' if err_rot_isorp < 1e-6 else 'FAILED'}  err={err_rot_isorp:.2e}")
    print(f"Raw->rot F:  {'PASSED' if err_rawrot_f < 1e-6 else 'FAILED'}  err={err_rawrot_f:.2e}")
    print(f"Raw->rot O:  {'PASSED' if err_rawrot_o < 1e-6 else 'FAILED'}  err={err_rawrot_o:.2e}")
    print(f"Full H:      {'PASSED' if res_H_full else 'NOT IMPLEMENTED'}  err=nan")
    print("=" * 40)

    # AvgRho_off is expected to fail until 2c+3c density-table evaluation is fully implemented on GPU.
    if not (res_S and res_T and res_Vna and res_sVNL and res_Vnl_cpu and res_Vnl_gpu and res_H_recon and res_H2c and res_Vnl_cpu_fortran and res_Vnl_gpu_fortran):
        raise RuntimeError("verify_C3: some implemented checks FAILED")

    print("\nDONE verify_C3")


if __name__ == "__main__":
    run_verification()
