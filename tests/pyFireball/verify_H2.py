import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall import elements
from pyBall.FireballOCL import OCL_Hamiltonian as ocl

# Set print options to show full matrices
np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)

def compare_matrices(name, fortran_mat, ocl_mat, tol=1e-5):
    print(f"\n--- Comparing {name} ---")
    diff = np.abs(fortran_mat - ocl_mat)
    max_diff = np.max(diff)
    print(f"Max difference: {max_diff:.2e}")
    print("Fortran Matrix:")
    print(fortran_mat)
    print("PyOpenCL Matrix:")
    print(ocl_mat)
    print("Abs Diff:")
    print(diff)
    
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        return False
    else:
        print(f"SUCCESS: {name} matches.")
        return True

def run_verification():
    # 1. Setup H2 molecule
    atomTypes_Z = np.array([1, 1], dtype=np.int32)
    atomPos = np.array([
        [0.0, 0.0, 0.0 ],
        [0.0, 0.0, 0.74]
    ], dtype=np.float64)
    
    fdata_dir = "./Fdata"
    
    # 2. Initialize Fortran FireCore
    print("Initializing Fortran FireCore...")
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    fc.setVerbosity(0)
    
    # 3. Initialize PyOpenCL Hamiltonian
    print("Initializing PyOpenCL Hamiltonian...")
    ham = ocl(fdata_dir)
    ham.prepare_splines(atomTypes_Z)
    
    # Helper to get Fortran H/S and neighbor list from sparse output
    def get_fortran_HS(ioffs):
        # ioffs: [S, T, Vna, Vnl, Vxc, Vca, Vxc_ca]
        fc.set_options(*ioffs)
        fc.assembleH(positions=atomPos, iforce=0, Kscf=1)
        dims = fc.get_HS_dims()
        sparse_data = fc.get_HS_sparse(dims)
        # For H2 (1s orbital), matrices are 2x2
        # However, sparse_data returns neighbors. We need to reconstruct.
        S = np.zeros((2, 2))
        H = np.zeros((2, 2))
        
        # Reconstrcut from sparse
        # neigh_j is [natoms, max_neigh]
        neighbors = []
        for i in range(2):
            for inigh in range(sparse_data.neighn[i]):
                j = sparse_data.neigh_j[i, inigh] - 1
                if j < 0 or j >= 2:
                    continue
                neighbors.append((i, j))
                S[i, j] = sparse_data.s_mat[i, inigh, 0, 0]
                H[i, j] = sparse_data.h_mat[i, inigh, 0, 0]
        return H, S, neighbors

    # --- Test 1: Overlap S ---
    print("\nTesting Overlap S...")
    H_f, S_f, neighbors = get_fortran_HS([1, 0, 0, 0, 0, 0, 0])
    _, S_o_blocks = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=False)
    S_o = np.zeros((2, 2))
    for idx, (i, j) in enumerate(neighbors):
        S_o[i, j] = S_o_blocks[idx, 0, 0]
    res_S = compare_matrices("Overlap S", S_f, S_o)

    # --- Test 2: Kinetic T ---
    print("\nTesting Kinetic T...")
    H_f, S_f, neighbors = get_fortran_HS([0, 1, 0, 0, 0, 0, 0])
    pairs_T = [(i, j, ham.species_pair_map[('kinetic', 1, 1)]) for i, j in neighbors]
    neigh_arr = np.array([p[:2] for p in pairs_T], dtype=np.int32)
    type_arr = np.array([p[2] for p in pairs_T], dtype=np.int32)
    T_o_blocks = ham.assemble_2c(atomPos, neigh_arr, type_arr)
    T_o = np.zeros((2, 2))
    for idx, (i, j) in enumerate(neighbors):
        T_o[i, j] = T_o_blocks[idx, 0, 0]
    res_T = compare_matrices("Kinetic T", H_f, T_o)

    # --- Test 3: Vna ---
    print("\nTesting Vna...")
    H_f, S_f, neighbors = get_fortran_HS([0, 0, 1, 0, 0, 0, 0])
    H_o_blocks, _ = ham.assemble_full(atomPos, atomTypes_Z, neighbors, include_T=False, include_Vna=True)
    V_o = np.zeros((2, 2))
    for idx, (i, j) in enumerate(neighbors):
        V_o[i, j] = H_o_blocks[idx, 0, 0]
    
    res_Vna = compare_matrices("Vna", H_f, V_o)

    print("\n" + "="*40)
    print("VERIFICATION SUMMARY")
    print(f"Overlap S: {'PASSED' if res_S else 'FAILED'}")
    print(f"Kinetic T: {'PASSED' if res_T else 'FAILED'}")
    print(f"Vna:       {'PASSED' if res_Vna else 'FAILED'}")
    print("="*40)

if __name__ == "__main__":
    run_verification()
