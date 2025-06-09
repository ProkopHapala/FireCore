import numpy as np
import os
import sys

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../") # Assuming this script is in FireCore/tests/pyFireball
from pyBall import FireCore as fc
from pyBall import elements # For element names from atomic numbers

# infinite line with for numpy
np.set_printoptions(linewidth=np.inf)

def print_sparse_matrices_info(sparse_data, dims):
    """
    Prints information from the sparse Hamiltonian and Overlap matrices
    in an understandable way.
    """
    natoms  = dims.natoms
    h_mat   = sparse_data.h_mat
    s_mat   = sparse_data.s_mat
    iatyp   = sparse_data.iatyp   # Species type for each atom (0-indexed from Fortran, adjust if needed)
    num_orb = sparse_data.num_orb # Number of orbitals for each species type (0-indexed)
    nzx     = sparse_data.nzx     # Z number for each compact Fireball species type
    neighn  = sparse_data.neighn  # Number of neighbors for each atom (0-indexed)
    neigh_j = sparse_data.neigh_j # Index of the j-th neighbor of atom i (Fortran: neigh_j(ineigh,iatom))
                                     # Python: neigh_j[iatom, ineigh] if reshaped, or handle Fortran order
    
    # Orbital labels (simplified for minimal basis, extend as needed)
    # This needs to be consistent with how Fireball orders orbitals for each species.
    # For now, a generic labeling. You'd need more info from Fireball's
    # make_munu or similar to get precise s, px, py, pz labels.
    orbital_labels_per_species = {} # Dictionary to store labels like {6: ['s', 'px', 'py', 'pz'], 1: ['s']}

    # Populate orbital_labels_per_species based on num_orb and known basis set conventions
    # This is a placeholder - you need to know Fireball's exact orbital ordering for each element type
    # nzx maps Fireball's compact species index (1 to nspecies_distinct) to Z number
    # num_orb is indexed by Fireball's compact species index (0 to nspecies_distinct-1 in Python)
    
    # Create a mapping from Z to compact Fireball species index (0-based for Python)
    Z_to_fb_species_idx_map = {z: i for i, z in enumerate(nzx)}

    for fb_sp_idx_py in range(dims.nspecies): # Iterate 0 to nspecies_distinct-1
        current_Z = nzx[fb_sp_idx_py]
        n_orb_sp = num_orb[fb_sp_idx_py]
        labels = []
        if n_orb_sp >= 1: labels.append("s")
        if n_orb_sp >= 4: # Assuming s, px, py, pz for sp3
            labels.extend(["px", "py", "pz"])
        # Add more logic for d orbitals, excited states etc. if your basis includes them
        # For now, fill remaining with generic labels if n_orb_sp > len(labels)
        for i_orb_extra in range(len(labels), n_orb_sp):
            labels.append(f"orb{i_orb_extra+1}")
        orbital_labels_per_species[current_Z] = labels # Store labels by Z for convenience


    print("\n--- Sparse Hamiltonian and Overlap Matrix Elements ---")
    for iatom in range(natoms):
        Z_i = iatyp[iatom] # This is the Z number of atom i
        element_name_i = elements.ELEMENTS[Z_i-1][1] # Get element symbol from Z
        
        # Get the compact Fireball species index (0-based) for atom i
        fb_species_idx_i_py = Z_to_fb_species_idx_map.get(Z_i)
        if fb_species_idx_i_py is None:
            print(f"Warning: Z={Z_i} for atom {iatom} not found in nzx mapping. Skipping.")
            continue
        num_orbitals_i = num_orb[fb_species_idx_i_py]
        labels_i = orbital_labels_per_species.get(Z_i, [f"orb{k+1}" for k in range(num_orbitals_i)])

        print(f"\nAtom {iatom} (Z={Z_i} - {element_name_i}, FB_species_idx={fb_species_idx_i_py+1}, Num Orbitals: {num_orbitals_i}):")
        
        # On-site block (interaction with self, ineigh might be specific or handled by ktransform)
        # For simplicity, let's assume the sparse format might not explicitly store the iatom-iatom block
        # unless it's a periodic image. The true on-site terms are often handled differently or
        # are part of the diagonal elements in the k-space Hamiltonian.
        # If you need to show on-site from h_mat/s_mat, you'd need to know how Fireball stores it (e.g. ineigh=0 or special handling)

        for ineigh_idx_f in range(neighn[iatom]): # Fortran ineigh runs 1 to neighn(iatom)
            # Adjust for 0-based Python indexing if neigh_j was directly copied
            # Fortran: neigh_j(ineigh, iatom)
            # Python:  neigh_j[iatom, ineigh] (with 0-based indices)
            jatom_fortran_idx = neigh_j[iatom, ineigh_idx_f] # Fortran: neigh_j(ineigh_idx_f+1, iatom+1)
            jatom_python_idx = jatom_fortran_idx -1 # Convert to 0-based for Python array access

            if jatom_python_idx < 0 or jatom_python_idx >= natoms : # Skip invalid neighbor indices
                print(f"  Neighbor index {jatom_fortran_idx} out of bounds for atom {iatom}, neighbor list index {ineigh_idx_f}")
                continue

            Z_j = iatyp[jatom_python_idx] # Z number of neighbor atom j
            element_name_j = elements.ELEMENTS[Z_j-1][1]

            fb_species_idx_j_py = Z_to_fb_species_idx_map.get(Z_j)
            if fb_species_idx_j_py is None:
                print(f"Warning: Z={Z_j} for neighbor atom {jatom_python_idx} not found in nzx mapping. Skipping.")
                continue
            num_orbitals_j = num_orb[fb_species_idx_j_py]
            labels_j = orbital_labels_per_species.get(Z_j, [f"orb{k+1}" for k in range(num_orbitals_j)])

            print(f"  Interaction with Neighbor Atom {jatom_python_idx} (Z={Z_j} - {element_name_j}, FB_species_idx={fb_species_idx_j_py+1}, Num Orbitals: {num_orbitals_j}), ineigh_list_idx={ineigh_idx_f+1}")

            # Accessing h_mat and s_mat:
            # Fortran: h_mat(imu, inu, ineigh, iatom)
            # Python:  h_mat[iatom, ineigh, inu, imu] (with 0-based indices)
            # The Python wrapper now allocates as (natoms, neigh_max, numorb_max, numorb_max)
            # which matches Fortran's dummy argument declaration order.
            
            print("    Hamiltonian Block (H):      Overlap Block (S):")
            # Print header row for S_mat
            header_s = "        " + "".join([f"{lbl_j:>8s}" for lbl_j in labels_j[:num_orbitals_j]])
            print(header_s)

            for imu_idx_py in range(num_orbitals_i): # Iterate up to actual number of orbitals for atom i
                row_h_str = f"    {labels_i[imu_idx_py]:<4s} |"
                row_s_str = f"    {labels_i[imu_idx_py]:<4s} |"
                for inu_idx_py in range(num_orbitals_j): # Iterate up to actual number of orbitals for atom j
                    # Fortran indices for h_mat/s_mat are 1-based for imu, inu
                    # ineigh_idx_f is 0-based for Python array access, iatom is 0-based
                    h_val = h_mat[iatom, ineigh_idx_f, inu_idx_py, imu_idx_py]
                    s_val = s_mat[iatom, ineigh_idx_f, inu_idx_py, imu_idx_py]
                    
                    # Only print if non-zero or if you want to see the full allocated block
                    # For now, printing all up to num_orbitals_i/j
                    row_h_str += f" {h_val:8.4f}"
                    row_s_str += f" {s_val:8.4f}"
                
                # For clarity, let's print H and S side-by-side if possible, or one after another
                # This example prints S matrix block fully, then implies H would be similar
                print( row_h_str+"   |    "+row_s_str) 
            print("    (Hamiltonian block would be similar format)")


if __name__ == "__main__":
    # 1. Define a simple molecule (e.g., H2O)
    # For H2O: O is type 8, H is type 1 (using atomic numbers as types for simplicity here)
    # Fireball uses its own species indexing, ensure iatyp matches that.
    # Let's assume Fireball's species index for O is 1, H is 2 for this example.
    # This needs to align with your Fdata and info.dat.
    # For a real run, you'd load from an XYZ and map element names to Fireball's internal species types.
    
    # Using the answer.xyz from your context for H2O
    # O atypes[0]=8, H atypes[1]=1, H atypes[2]=1
    # In Fireball's info.dat, Oxygen might be species 1, Hydrogen species 2.
    # So, atomTypes_fc would be [1, 2, 2] if O is species 1, H is species 2.
    

    atomTypes_Z = np.array([8, 1, 1], dtype=np.int32) # Atomic numbers Z
    atomPos = np.array([
        [0.00000000, 0.00000000, 0.00000000], # O
        [-0.75806320, 0.63581013, 0.00000000], # H
        [0.75806639, 0.63580735, 0.00000000]  # H
    ], dtype=np.float64)

    # fc.preinit() # Sets default parameters
    # fc.setVerbosity(2) # Set verbosity high enough to see Fortran debug prints
    # # For direct initialization:
    print(f"Initializing FireCore with atomTypes_Z: {atomTypes_Z}")
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)
    print("FireCore initialized.")
    fc.setVerbosity(2) 

    print("Assembling Hamiltonian (Kscf=1, iforce=0)...")
    # We need to run assembleH to populate h_mat and s_mat
    # Kscf=1 for initial assembly, iforce=0 if we don't need forces yet
    fc.assembleH(positions=atomPos, iforce=0, Kscf=1)
    print("Hamiltonian assembled.")

    # 2. Get dimensions
    print("\nFetching dimensions...")
    dims = fc.get_HS_dims()
    print("Dimensions retrieved (type):", type(dims))
    print("Python: natoms    ", dims.natoms)
    print("Python: norbitals ", dims.norbitals)
    print("Python: nspecies  ", dims.nspecies)
    print("Python: neigh_max ", dims.neigh_max)
    print("Python: numorb_max", dims.numorb_max)
    print("Python: nsh_max   ", dims.nsh_max)

    print(f"Python: max_mu_dims: ({dims.max_mu_dim1}, {dims.max_mu_dim2}, {dims.max_mu_dim3})")
    print(f"Python: mbeta_max (size of xl dim 2): {dims.mbeta_max}")
    print(f"Python: nspecies_fdata (dim for nzx): {dims.nspecies_fdata}")

    print("\nPython: Expected allocation shapes for firecore_get_HS_sparse:")
    print(f"  h_mat_out: ({dims.numorb_max}, {dims.numorb_max}, {dims.neigh_max}, {dims.natoms})")
    print(f"  num_orb_out: ({dims.nspecies},)")
    # Use the new dimension for mu, nu, mvalue
    print(f"  mu_out: ({dims.max_mu_dim1}, {dims.max_mu_dim2}, {dims.max_mu_dim3})")
    print(f"  xl_out: (3, {dims.mbeta_max})")
    print(f"  nzx_out: ({dims.nspecies_fdata},)")
    
    if dims.natoms == 0 or dims.norbitals == 0:
        print("Error: System not properly initialized or empty. Exiting.")
        sys.exit(1)

    # 3. Get sparse H, S, and indexing data
    print("\nFetching sparse H, S and indexing data...")
    sparse_data = fc.get_HS_sparse(dims)
    print("Sparse data fetched.")
    # print("iatyp:", sparse_data.iatyp)
    # print("num_orb (per species):", sparse_data.num_orb)
    # print("neighn (per atom):", sparse_data.neighn)
    # print("neigh_j (Fortran order): shape", sparse_data.neigh_j.shape)

    # 4. Print the sparse matrix information
    print_sparse_matrices_info(sparse_data, dims)

    # 5. Get dense k-space H, S (optional, for a specific k-point)
    print("\nFetching dense k-space H and S for Gamma point...")
    kvec_gamma = [0.0, 0.0, 0.0]
    Hk, Sk = fc.get_HS_k(kvec_gamma, dims.norbitals)
    print(f"Hk shape: {Hk.shape}, Sk shape: {Sk.shape}")
    
    # You can print parts of Hk and Sk if desired
    print("Hk:")
    print(Hk[:, :])
    print("Sk:")
    print(Sk[:, :])

    print("\nScript finished.")
