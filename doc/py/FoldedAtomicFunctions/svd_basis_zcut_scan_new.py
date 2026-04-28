#!/usr/bin/env python3
"""svd_basis_zcut_scan_new.py

Optimized version to study singular value decay and reconstruction error
of Morse samples with cutoff-polynomial basis, varying z_cut.
Uses utility functions for core analysis to keep this script concise.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

# project helpers -------------------------------------------------------------
from basis_utils import gen_morse_prms, gen_morse_curves, construct_composite_cutoff_basis
from optimize import fit_coefficients
import plot_utils as pu # custom plotting helpers

# ----------------------------------------------------------------------------
# 0. Parameters
# ----------------------------------------------------------------------------

# z grid
z = np.linspace(0.0, 12.0, 250)

# Morse sample family  – vary both r0 and a
n_samples = 20
rng = np.random.default_rng(0)
prms_list = gen_morse_prms(n_samples, (1.5, 2.0), (2.0, 4.0), 1.0)
samples, _ = gen_morse_curves(z, prms=prms_list)

# Composite cutoff-poly basis parameters
# Each tuple is (z_cut, list_of_n_factors)
# where n_factor corresponds to (z_cut - z)^(2*n_factor)
# basis_definitions_custom = [
#     (8.0, [2, 3, 4]),    # 3 functions with z_cut=8.0 and powers 4, 6, 8
#     (6.0, [1, 2]),       # 2 functions with z_cut=6.0 and powers 2, 4
#     (10.0, [3, 5])       # 2 functions with z_cut=10.0 and powers 6, 10
# ]

basis_definitions_custom = [
    #(10.0, [1,2,3,4,5,6,7,8]),  
    #(10.0, [2,3,4,5,6,7,8]),  
    #(10.0, [4,6,8]),  
    #(12.0, [4,5,6,7,8]),  
    #(6.0, [1,2,3,4])   
    #(8.0, [2,3,4,5])   
    (8.0, [2,4,6,8])   
    #(4.0, [1,2])      
]


# ----------------------------------------------------------------------------
# 1. Generate reference sample functions (columns of Y)
# ----------------------------------------------------------------------------

Y = np.vstack(samples)  # (M, Nz)
# plot1D(z, Y, "Example Morse samples") # Optional: plot original samples

# Define per-sample mask: 1.0 for relevant regions, 0.0 for ignored regions
V_repulsive_thresh = 0.5   # eV, values AT or ABOVE this threshold are ignored for each sample

# weights_Y_mask will be (M, Nz), with 1.0 where Y < V_repulsive_thresh, and 0.0 otherwise.
weights_Y_mask = (Y < V_repulsive_thresh).astype(float)

# ----------------------------------------------------------------------------
# 2. Construct the composite basis and evaluate its error
# ----------------------------------------------------------------------------
print(f"Constructing composite basis from definitions: {basis_definitions_custom}")
Phi_composite, composite_labels = construct_composite_cutoff_basis(z, basis_definitions_custom, z0=z.min())

if Phi_composite.shape[0] == 0:
    print("Error: The defined custom basis is empty. Exiting.")
    exit()

print(f"Constructed composite basis with {Phi_composite.shape[0]} functions.")
print("Labels of composite basis functions:")
for i, label in enumerate(composite_labels):
    print(f"  {i}: {label}")

# Plot the constructed basis functions
pu.plot1D(z, Phi_composite, title=f"Constructed Composite Basis ({Phi_composite.shape[0]} functions)", labels=composite_labels, ylims=(0,1.0))

# Fit coefficients using the composite basis and per-sample masked regions
S_composite = fit_coefficients(Y, Phi_composite, weights=weights_Y_mask)

# Reconstruct using the composite basis
Y_reconstructed_composite = S_composite.T @ Phi_composite

# Calculate total error for the composite basis
total_masked_error_composite = np.linalg.norm((Y - Y_reconstructed_composite) * weights_Y_mask)

print(f"\nTotal masked error for the composite basis: {total_masked_error_composite:.2e}")

# ----------------------------------------------------------------------------
# 3. Show fits for example samples using the composite basis
# ----------------------------------------------------------------------------
rng_plot        = np.random.default_rng(1)
idx_samples_plot = rng_plot.choice(Y.shape[0], size=min(5, Y.shape[0]), replace=False)

data_pairs_plot = []
for idx in idx_samples_plot:
    y_samp = Y[idx]
    # For individual sample fitting, use its specific coefficients from S_composite
    y_fit_samp  = S_composite[:, idx].T @ Phi_composite

    # Error for this specific sample, considering its mask
    error_samp_masked  = np.linalg.norm((y_samp - y_fit_samp) * weights_Y_mask[idx, :])
    print(f"Absolute masked reconstruction error for sample {idx} using composite basis: {error_samp_masked:.2e}")
    data_pairs_plot.append((y_samp, y_fit_samp))

# Use the new function to plot multiple samples and their approximations
pu.plotMultiFunctionApprox(z, data_pairs_plot, bError=False, errMax=0.01, scMin=1.1, title=f"Reconstruction using Composite Basis ({Phi_composite.shape[0]} functions)")

plt.show()