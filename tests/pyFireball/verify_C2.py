import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL import OCL_Hamiltonian as ocl

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

    print("Initializing PyOpenCL Hamiltonian...")
    ham = ocl(fdata_dir)
    ham.prepare_splines(atomTypes_Z)

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

    def get_fortran_HS(ioffs):
        # ioffs: [S, T, Vna, Vnl, Vxc, Vca, Vxc_ca]
        fc.set_options(*ioffs)
        fc.assembleH(positions=atomPos, iforce=0, Kscf=1)
        dims = fc.get_HS_dims()
        sparse_data = fc.get_HS_sparse(dims)
        natoms = 2
        H = _blocked_to_dense(sparse_data, sparse_data.h_mat, natoms)
        S = _blocked_to_dense(sparse_data, sparse_data.s_mat, natoms)
        # Build neighbor list (iatom,jatom) in the same order as sparse export
        neighbors = []
        for i in range(natoms):
            for ineigh in range(int(sparse_data.neighn[i])):
                j = int(sparse_data.neigh_j[i, ineigh]) - 1
                if j < 0 or j >= natoms:
                    continue
                neighbors.append((i, j))
        return H, S, neighbors, sparse_data

    print("\nTesting Overlap S...")
    H_f, S_f, neighbors, sparse_data = get_fortran_HS([1, 0, 0, 0, 0, 0, 0])
    _, S_o_blocks = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=False, include_Vnl=False)
    n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    norb = int(offs[-1])
    S_o = np.zeros((norb, norb), dtype=np.float64)
    for idx, (i, j) in enumerate(neighbors):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        # blocks are (inu,imu); map to dense (imu,inu)
        S_o[i0:i0+ni, j0:j0+nj] = S_o_blocks[idx, :nj, :ni].T
    res_S = compare_matrices("Overlap S", S_f, S_o)

    print("\nTesting Kinetic T...")
    H_f, S_f, neighbors, sparse_data = get_fortran_HS([0, 1, 0, 0, 0, 0, 0])
    H_o_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=True, include_Vna=False, include_Vnl=False)
    n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    norb = int(offs[-1])
    T_o = np.zeros((norb, norb), dtype=np.float64)
    for idx, (i, j) in enumerate(neighbors):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        T_o[i0:i0+ni, j0:j0+nj] = H_o_blocks[idx, :nj, :ni].T
    res_T = compare_matrices("Kinetic T", H_f, T_o)

    print("\nTesting Vna...")
    H_f, S_f, neighbors, sparse_data = get_fortran_HS([0, 0, 1, 0, 0, 0, 0])
    H_o_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=True, include_Vnl=False)
    n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    norb = int(offs[-1])
    Vna_o = np.zeros((norb, norb), dtype=np.float64)
    for idx, (i, j) in enumerate(neighbors):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        Vna_o[i0:i0+ni, j0:j0+nj] = H_o_blocks[idx, :nj, :ni].T
    res_Vna = compare_matrices("Vna", H_f, Vna_o)

    print("\nTesting Vnl...")
    H_f, S_f, neighbors, sparse_data = get_fortran_HS([0, 0, 0, 1, 0, 0, 0])
    H_o_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=False, include_Vnl=True)
    n_orb_atom, offs = _orbital_layout(sparse_data, 2)
    norb = int(offs[-1])
    Vnl_o = np.zeros((norb, norb), dtype=np.float64)
    for idx, (i, j) in enumerate(neighbors):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        Vnl_o[i0:i0+ni, j0:j0+nj] = H_o_blocks[idx, :nj, :ni].T
    require_nz = np.max(np.abs(H_f)) > 1e-5
    res_Vnl = compare_matrices("Vnl", H_f, Vnl_o, require_nonzero=require_nz)

    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY")
    print(f"Overlap S: {'PASSED' if res_S else 'FAILED'}")
    print(f"Kinetic T: {'PASSED' if res_T else 'FAILED'}")
    print(f"Vna:       {'PASSED' if res_Vna else 'FAILED'}")
    print(f"Vnl:       {'PASSED' if res_Vnl else 'FAILED'}")
    print("=" * 40)


if __name__ == "__main__":
    run_verification()
