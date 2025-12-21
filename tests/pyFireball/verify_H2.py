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
    
    if max_diff > tol:
        print(f"WARNING: {name} discrepancy exceeds tolerance!")
        print("Fortran Matrix:")
        print(fortran_mat)
        print("PyOpenCL Matrix:")
        print(ocl_mat)
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
    
    # Define neighbors for H2
    neighbors = [(0, 1), (1, 0), (0, 0), (1, 1)] # include on-site
    
    # Helper to get Fortran H/S
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
        for i in range(2):
            for inigh in range(sparse_data.neighn[i]):
                j = sparse_data.neigh_j[i, inigh] - 1
                if j < 2:
                    S[i, j] = sparse_data.s_mat[i, inigh, 0, 0]
                    H[i, j] = sparse_data.h_mat[i, inigh, 0, 0]
        return H, S

    # --- Test 1: Overlap S ---
    print("\nTesting Overlap S...")
    H_f, S_f = get_fortran_HS([1, 0, 0, 0, 0, 0, 0])
    _, S_o_blocks = ham.assemble_full(atomPos, atomTypes_Z, neighbors)
    S_o = np.zeros((2, 2))
    for idx, (i, j) in enumerate(neighbors):
        S_o[i, j] = S_o_blocks[idx, 0, 0]
    res_S = compare_matrices("Overlap S", S_f, S_o)

    # --- Test 2: Kinetic T ---
    print("\nTesting Kinetic T...")
    H_f, S_f = get_fortran_HS([0, 1, 0, 0, 0, 0, 0])
    pairs_T = [(i, j, ham.species_pair_map[('kinetic', 1, 1)]) for i, j in neighbors]
    T_o_blocks = ham.assemble_2c(atomPos, np.array([p[:2] for p in pairs_T]), np.array([p[2:] for p in pairs_T]))
    T_o = np.zeros((2, 2))
    for idx, (i, j) in enumerate(neighbors):
        T_o[i, j] = T_o_blocks[idx, 0, 0]
    res_T = compare_matrices("Kinetic T", H_f, T_o)

    # --- Test 3: Vna ---
    print("\nTesting Vna...")
    H_f, S_f = get_fortran_HS([0, 0, 1, 0, 0, 0, 0])
    # For Vna, the Fortran side sums ontopl, ontopr, and atom contributions.
    # Our OCL currently only has one 'vna' interaction type per pair.
    # Let's try to match by summing OCL components if they exist.
    V_o = np.zeros((2, 2))
    for root in ['vna_atom_00', 'vna_ontopl_00', 'vna_ontopr_00', 'vna']:
        t = ham.species_pair_map.get((root, 1, 1))
        if t is not None:
            print(f"  Adding OCL component: {root}")
            pairs = [(i, j, t) for i, j in neighbors]
            blocks = ham.assemble_2c(atomPos, np.array([p[:2] for p in pairs]), np.array([p[2:] for p in pairs]))
            for idx, (i, j) in enumerate(neighbors):
                V_o[i, j] += blocks[idx, 0, 0]
    
    # Firecore might have a huge diagonal due to sum_j <i | v_j | i>
    # but verify against what it exports.
    res_Vna = compare_matrices("Vna", H_f, V_o)

    print("\n" + "="*40)
    print("VERIFICATION SUMMARY")
    print(f"Overlap S: {'PASSED' if res_S else 'FAILED'}")
    print(f"Kinetic T: {'PASSED' if res_T else 'FAILED'}")
    print(f"Vna:       {'PASSED' if res_Vna else 'FAILED'}")
    print("="*40)

if __name__ == "__main__":
    run_verification()
