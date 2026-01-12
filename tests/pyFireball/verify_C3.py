import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import OCL_Hamiltonian

np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)


def compare_blocks(name, A, B, tol=1e-5, require_nonzero=False):
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

    print("\n[PLUMBING] compute_avg_rho (3c gather) using Fortran-exported rho + Qin-shell...")
    dims = fc.get_HS_dims(force_refresh=True)
    sd = fc.get_HS_neighs(dims)
    sd = fc.get_HS_sparse(dims, data=sd)
    sd = fc.get_rho_sparse(dims, data=sd)
    sd = fc.get_rho_off_sparse(dims, data=sd)

    # OpenCL compute_avg_rho weights 3c density pieces by *neutral* charge of the common neighbor.
    # Fortran average_rho uses Qneutral(isorp, indna) (per-shell neutral population for species indna).
    # Here we pass a per-atom scalar Qneutral_total(atom) = sum_shell Qneutral(shell, species(atom)).
    Qneutral_sh = fc.get_Qneutral_shell(dims)  # [nsh_max, nspecies]
    nzx = np.array(sd.nzx, dtype=np.int32)
    iatyp = np.array(sd.iatyp, dtype=np.int32)
    Qatom = np.zeros(dims.natoms, dtype=np.float32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec = int(w[0])
        Qatom[ia] = float(np.sum(Qneutral_sh[:, ispec]))

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
        rho_blocks[ip, :, :] = sd.rho[int(ia), ineigh, :4, :4]

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

    rho_avg_blocks = ham.compute_avg_rho_3c(atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qatom)

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

    ok = True
    maxdiff = 0.0
    if np.any(mask == 1):
        ok, maxdiff = compare_blocks("AvgRho_off (Fortran rho_off vs OpenCL compute_avg_rho)", ref_blocks[mask == 1], rho_avg_blocks[mask == 1], tol=1e-4, require_nonzero=True)
    else:
        print("WARNING: No comparable neighbor blocks found (mask empty)")
        ok = False

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

    print("\n========================================")
    print("VERIFICATION SUMMARY")
    print(f"AvgRho_off:   {'PASSED' if ok else 'FAILED'}   (max diff {maxdiff:.2e})")
    print("========================================")
    if not ok:
        raise RuntimeError("verify_C3: AvgRho_off comparison FAILED")

    print("\nDONE verify_C3")


if __name__ == "__main__":
    run_verification()
