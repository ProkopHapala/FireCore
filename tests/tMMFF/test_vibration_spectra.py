import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

sys.path.append("../../")
from pyBall import MMFF as mmff
from pyBall import FTIR

def main():
    parser = argparse.ArgumentParser(description="Linear response vibration spectroscopy using MMFF.")
    parser.add_argument("-i", "--input", type=str, default="../../web/common_resources/mol/adamantane.mol2", help="Input molecule file (.mol2 or .xyz)")
    parser.add_argument("--fmin", type=float, default=0.01, help="Minimum frequency for scan")
    parser.add_argument("--fmax", type=float, default=10.0, help="Maximum frequency for scan")
    parser.add_argument("--nfreq", type=int, default=1000, help="Number of frequency samples")
    parser.add_argument("--eta", type=float, default=0.05, help="Damping parameter (broadening)")
    parser.add_argument("--dx", type=float, default=1e-4, help="Finite difference step for Hessian")
    parser.add_argument("--shift", type=float, default=1e6, help="Rigid body mode penalty shift")
    parser.add_argument("--out_spec", type=str, default="vibration_spectrum.png", help="Output spectrum image")
    parser.add_argument("--save_modes", action="store_true", default=True, help="Save eigenmodes to XYZ")
    parser.add_argument("--save_responses", action="store_true", default=True, help="Save response vectors to XYZ")
    parser.add_argument("--save_excitation", action="store_true", default=True, help="Save excitation vectors to XYZ")
    parser.add_argument("--verbosity", type=int, default=0, help="MMFF verbosity")
    
    args = parser.parse_args()

    # 1. Initialize MMFF
    mmff.setVerbosity(verbosity=args.verbosity, idebug=0)
    mmff.init(xyz_name=args.input, bEpairs=False, bMMFF=True)

    # 2. Configure forcefield
    mmff.setSwitches(
        NonBonded=-1,
        MMFF=+1,
        Angles=+1,
        SurfAtoms=-1,
        GridFF=-1,
        PiSigma=-1,
        PiPiI=-1
    )

    mmff.getBuffs()
    n_atoms = mmff.natoms
    print(f"Loaded {n_atoms} atoms from {args.input}")

    # 3. Compute Hessian
    inds = np.arange(n_atoms, dtype=np.int32)
    H = mmff.getHessian3Nx3N(inds, dx=args.dx)
    H = 0.5 * (H + H.T)

    # 4. Compute Mass matrix
    M = FTIR.get_mass_matrix(mmff, n_atoms)
    M_inv_sqrt = np.diag(1.0 / np.sqrt(np.diag(M)))

    # 4.5 Project out rigid body modes
    pos = mmff.apos[:n_atoms].copy()
    H = FTIR.project_rigid_modes(H, M, pos, shift=args.shift)

    # 5. Define elements and excitation vectors
    elements = []
    for i in range(n_atoms):
        type_name = mmff.getTypeName(i)
        elem_str = type_name.split('_')[0]
        if len(elem_str) == 0: elem_str = "C"
        elements.append(elem_str)

    charges = np.ones(n_atoms)
    com = np.mean(pos, axis=0)
    
    # Symmetric stretching excitations
    dir_x = np.zeros((n_atoms, 3)); dir_x[:, 0] = pos[:, 0] - com[0]
    dir_y = np.zeros((n_atoms, 3)); dir_y[:, 1] = pos[:, 1] - com[1]
    dir_z = np.zeros((n_atoms, 3)); dir_z[:, 2] = pos[:, 2] - com[2]

    if args.save_excitation:
        FTIR.save_xyz_vib("excitation_x.xyz", elements, pos, dir_x, comments="Symmetric X Excitation")
        print("Saved excitation_x.xyz")

    # 6. Probing the vibration spectrum
    omegas = np.linspace(args.fmin, args.fmax, args.nfreq)
    print(f"Computing vibration spectrum (eta={args.eta})...")
    res_x = FTIR.mechanical_greens_probing(H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)
    res_y = FTIR.mechanical_greens_probing(H, M, omegas, eta=args.eta, direction_vec=dir_y, charges=charges)
    res_z = FTIR.mechanical_greens_probing(H, M, omegas, eta=args.eta, direction_vec=dir_z, charges=charges)

    mag_x = np.linalg.norm(res_x["dipole"], axis=1)
    mag_y = np.linalg.norm(res_y["dipole"], axis=1)
    mag_z = np.linalg.norm(res_z["dipole"], axis=1)

    # 7. Diagonalization for Eigenmodes
    print("Diagonalizing mass-weighted Hessian...")
    Hm = M_inv_sqrt @ H @ M_inv_sqrt
    eigvals, eigvecs = np.linalg.eigh(Hm)
    eig_freqs = np.sqrt(np.maximum(eigvals, 0.0))

    # 8. Save Eigenmodes and Response Vectors to XYZ
    if args.save_modes or args.save_responses:
        print("Saving vectors to XYZ...")
        u_modes = []
        comments_modes = []
        u_responses = []
        comments_responses = []

        valid_indices = [i for i, f in enumerate(eig_freqs) if args.fmin < f < args.fmax]
        for i in valid_indices:
            f = eig_freqs[i]
            # Eigenmode
            u = M_inv_sqrt @ eigvecs[:, i]
            u = u / np.linalg.norm(u)
            u_modes.append(u.reshape(n_atoms, 3))
            comments_modes.append(f"Mode {i} Freq {f:.4f}")
            
            # Response at resonance
            u_resp = FTIR.solve_response(H, M, f, eta=args.eta, direction_vec=dir_x, charges=charges)
            u_responses.append(np.imag(u_resp))
            comments_responses.append(f"Response (Imag) to X-excitation at Freq {f:.4f}")

        if args.save_modes and u_modes:
            FTIR.save_xyz_vib("eigenmodes.xyz", elements, pos, np.array(u_modes), comments=comments_modes)
            print("Saved eigenmodes.xyz")
        if args.save_responses and u_responses:
            FTIR.save_xyz_vib("response_x.xyz", elements, pos, np.array(u_responses), comments=comments_responses)
            print("Saved response_x.xyz")

    # 9. Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(omegas, mag_x, label="X-excitation")
    plt.plot(omegas, mag_y, label="Y-excitation")
    plt.plot(omegas, mag_z, label="Z-excitation")

    valid_freqs = eig_freqs[(eig_freqs > args.fmin) & (eig_freqs < args.fmax)]
    for f in valid_freqs:
        plt.axvline(x=f, color='gray', linestyle='--', alpha=0.3)

    plt.xlabel("Frequency $\omega$")
    plt.ylabel("|Dipole Response|")
    plt.title(f"FTIR Vibration Spectrum: {os.path.basename(args.input)}\n(Bonding Only, Symmetric Excitation)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(args.out_spec)
    print(f"Spectrum saved to {args.out_spec}")

    mmff.clear()

if __name__ == "__main__":
    main()
