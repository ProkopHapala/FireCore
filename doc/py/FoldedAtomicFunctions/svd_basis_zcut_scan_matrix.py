#!/usr/bin/env python3
"""svd_basis_zcut_scan_matrix.py

Scans a 2D parameter space (z_cut, n0_start_power) for cutoff_poly_basis,
evaluates total masked fitting error, and plots it as a 2D matrix.
The basis for each point (z_cut, n0) is [(z_cut, [n0, n0+1, ..., n0+nn-1])].
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# project helpers -------------------------------------------------------------
from basis_utils import gen_morse_prms, gen_morse_curves, construct_composite_cutoff_basis
from optimize import fit_coefficients
import plot_utils

# ----------------------------------------------------------------------------
# 0. Parameters
# ----------------------------------------------------------------------------

# z grid
z_grid = np.linspace(1.5, 12.0, 250)

# Morse sample family
n_samples = 20
# For reproducibility with gen_morse_prms (which uses global np.random):
np.random.seed(0) 
prms_list = gen_morse_prms(n_samples, (1.0, 2.0), (2.0, 4.0), 1.0)
samples, _ = gen_morse_curves(z_grid, prms=prms_list)

# Basis construction parameters for the scan
z_cut_scan_values = np.linspace(4.0, 16.0, 13)  # Axis 1 for the error matrix
n0_start_scan_values = np.arange(1, 7)       # Axis 2 for the error matrix (n0 must be >= 1)
# dns defines the offsets from n0_start to form the power factors.
#dns = [0, 1, 2, 3]  # e.g., if n0=1, powers are [1,2,3,4]. If dns=[0,2], powers are [n0, n0+2]
#dns = [0, 1] 
#dns = [0,4] 
dns = [0] 


V_repulsive_thresh = 0.5  # eV, for masking

# ----------------------------------------------------------------------------
# 1. Generate reference sample functions and mask
# ----------------------------------------------------------------------------
Y_samples = np.vstack(samples)  # (M, Nz)
weights_Y_mask = (Y_samples < V_repulsive_thresh).astype(float) # (M, Nz)

# ----------------------------------------------------------------------------
# 2. Perform the 2D scan and calculate errors
# ----------------------------------------------------------------------------
error_matrix_rmse = np.zeros((len(z_cut_scan_values), len(n0_start_scan_values)))
error_matrix_max = np.zeros((len(z_cut_scan_values), len(n0_start_scan_values)))

print(f"Starting 2D scan: {len(z_cut_scan_values)} z_cuts, {len(n0_start_scan_values)} n0_starts, with dns={dns}")

for i, z_c_val in enumerate(z_cut_scan_values):
    for j, n0_val in enumerate(n0_start_scan_values):
        power_factors = [int(d + n0_val) for d in dns] # Ensure Python ints
        if not power_factors or power_factors[0] <= 0: # Ensure n0 is valid
            print(f"Skipping z_cut={z_c_val:.2f}, n0={n0_val} (invalid power_factors: {power_factors})")
            error_matrix_rmse[i, j] = np.inf
            error_matrix_max[i, j] = np.inf
            continue

        current_basis_def = [(z_c_val, power_factors)]
        
        Phi_current, _ = construct_composite_cutoff_basis(z_grid, current_basis_def)

        if Phi_current.shape[0] == 0:
            print(f"Warning: Empty basis for z_cut={z_c_val:.2f}, n0={n0_val}. Assigning infinite error.")
            error_matrix_rmse[i, j] = np.inf
            error_matrix_max[i, j] = np.inf
            continue

        S_fit = fit_coefficients(Y_samples, Phi_current, weights=weights_Y_mask)
        Y_reconstructed = S_fit.T @ Phi_current
        
        # Calculate masked errors for each sample, then aggregate
        errors_per_sample_sq = np.sum(((Y_samples - Y_reconstructed) * weights_Y_mask)**2, axis=1)
        
        # RMSE over all *masked* points across all samples
        total_masked_points = np.sum(weights_Y_mask) # Count of points where weight is 1.0
        if total_masked_points > 0:
            total_sum_sq_error = np.sum(errors_per_sample_sq) # This is sum of (y-y_fit)^2 over all masked points
            error_matrix_rmse[i, j] = np.sqrt(total_sum_sq_error / total_masked_points)
        else:
            error_matrix_rmse[i, j] = np.inf # Or 0 if no points to fit

        # Max error over all *masked* points (absolute difference)
        abs_diff_masked = np.abs((Y_samples - Y_reconstructed) * weights_Y_mask)
        error_matrix_max[i, j] = np.max(abs_diff_masked) if total_masked_points > 0 else np.inf

        print(f"  z_cut={z_c_val:.2f}, n0={n0_val}, powers={power_factors}: RMSE={error_matrix_rmse[i,j]:.2e}, MaxErr={error_matrix_max[i,j]:.2e}")

print("Scan completed.")

# ----------------------------------------------------------------------------
# 3. Plot the error matrix
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

im_rmse = ax[0].imshow(error_matrix_rmse.T, origin='lower', aspect='auto',
                   extent=[z_cut_scan_values[0], z_cut_scan_values[-1], n0_start_scan_values[0]-0.5, n0_start_scan_values[-1]+0.5],
                   norm=LogNorm()) # Log scale for errors
fig.colorbar(im_rmse, ax=ax[0], label='Total Masked RMSE')
ax[0].set_xlabel('z_cut (Å)')
ax[0].set_ylabel('n0_start_power_offset')
ax[0].set_title(f'Total Masked RMSE (dns={dns})')

im_max = ax[1].imshow(error_matrix_max.T, origin='lower', aspect='auto',
                   extent=[z_cut_scan_values[0], z_cut_scan_values[-1], n0_start_scan_values[0]-0.5, n0_start_scan_values[-1]+0.5],
                   norm=LogNorm())
fig.colorbar(im_max, ax=ax[1], label='Max Masked Absolute Error')
ax[1].set_xlabel('z_cut (Å)')
ax[1].set_ylabel('n0_start_power_offset')
ax[1].set_title(f'Max Masked Absolute Error (dns={dns})')

plt.tight_layout()
plt.show()