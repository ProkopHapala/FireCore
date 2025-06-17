"""basis.py
Creation and orthogonalisation of analytical basis functions.
All functions are light-weight and purely functional.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import qr
from typing import Tuple, List, Dict
from basis_utils import poly_basis, cutoff_poly_basis, gram_schmidt_weighted

if __name__ == "__main__":
    z = np.linspace(0, 5, 100)

    print("--- Polynomial Basis Test ---")
    Phi, scale_info, labels = poly_basis(z, 4)
    Q = gram_schmidt_weighted(Phi)
    print("Phi shape", Phi.shape, "Q shape", Q.shape)
    print("Labels:", labels)
    print("Scale Info:", scale_info)

    from matplotlib import pyplot as plt
    from plot_utils import plot1D
    plot1D(z, Phi, "Raw polynomial basis")
    plot1D(z, Q, "Orthogonalised basis")

    print("\n--- Cutoff Polynomial Basis Test ---")
    z_cut_test = 3.0
    max_n_test = 3
    Phi_cutoff, scale_info_cutoff, labels_cutoff = cutoff_poly_basis(z, z_cut_test, max_n_test)
    Q_cutoff = gram_schmidt_weighted(Phi_cutoff)
    print("Phi_cutoff shape", Phi_cutoff.shape)
    print("Scale Info:", scale_info_cutoff)
    print("Labels:", labels_cutoff)
    plot1D(z, Phi_cutoff, f"Cutoff polynomial basis (z_cut={z_cut_test}, max_2n={2*max_n_test})")
    plot1D(z, Q_cutoff, f"Orthogonalised cutoff polynomial basis (z_cut={z_cut_test}, max_2n={2*max_n_test})")
    
    plt.show()

