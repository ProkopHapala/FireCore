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
from basis_utils import cutoff_poly_basis, gen_morse_prms, gen_morse_curves, calc_fit_svd, eval_multi_basis_recon_err, calc_svd_reconstruction_errors
from optimize import fit_coefficients  # functional least-squares helper
from plot_utils import plot1D, plotFunctionApprox  # custom plotting helpers

# ----------------------------------------------------------------------------
# 0. Parameters
# ----------------------------------------------------------------------------

# z grid
z = np.linspace(1.5, 12.0, 250)

# Morse sample family  – vary both r0 and a
n_samples = 100
rng = np.random.default_rng(0)
prms_list = gen_morse_prms(n_samples, (1.0, 2.0), (2.0, 4.0), 1.0)
samples, _ = gen_morse_curves(z, prms=prms_list)

# cutoff-poly basis parameters
max_n       = 8                # highest n in (z_cut - z)^(2 n)
zcut_values = np.linspace(4.0, 16.0, 25)

# ----------------------------------------------------------------------------
# 1. Generate reference sample functions (columns of Y)
# ----------------------------------------------------------------------------

Y = np.vstack(samples)  # (M, Nz)

# plot1D(z, np.array(samples)[:5], "Example Morse samples (first 5)")

# ----------------------------------------------------------------------------
# 2. Prepare list of basis sets and evaluate errors
# ----------------------------------------------------------------------------

phi_list_scan = []
for zc in zcut_values:
    Phi, *_ = cutoff_poly_basis(z, zc, max_n, scale=False)
    phi_list_scan.append(Phi)

ks_to_plot = [1,2,3,4,5,6] 

err_full, err_ks = eval_multi_basis_recon_err(Y, phi_list_scan, ks_to_plot)

print(f"Initial full basis fit errors for max_n={max_n} vs z_cut (relative): {['%.2e' % e for e in err_full]}")

# ----------------------------------------------------------------------------
# 3. Plot error vs z_cut for several nBasisOpt values
# ----------------------------------------------------------------------------
err_k = np.array([err_ks[k] for k in ks_to_plot])

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(zcut_values, err_full, 'o:', lw=1.5, color='black', label='Rel. Initial Full Basis Fit Error')
for i, k in enumerate(ks_to_plot):
    ax.plot(zcut_values, err_k[i], '.-', lw=0.5, label=f'Actual Recon Error k={k}')
ax.set_yscale('log')
ax.set_xlabel('z_cut (Å)')
ax.set_ylabel('Relative Error')
ax.set_title( f'Error vs z_cut: Initial Fit and SVD Compression (max_n={max_n})')
ax.legend()
ax.grid(True, which="both", ls=":")

# ----------------------------------------------------------------------------
# 4. Choose optimal z_cut for chosen nBasisOpt and show fits
# ----------------------------------------------------------------------------
nBasisOpt = 6  # user-selected principal components
errors_k_actual_recon = err_ks[nBasisOpt]
zc_best_idx = int(np.argmin(errors_k_actual_recon))
zc_best = zcut_values[zc_best_idx]

print(f"Best z_cut for k={nBasisOpt} (Actual Recon Error): {zc_best:.2f} Å (Actual Recon rel. error {min(errors_k_actual_recon):.2e})")
print(f"Initial full basis fit rel. error at z_cut={zc_best:.2f}: {err_full[zc_best_idx]:.2e}")

Phi_best = phi_list_scan[zc_best_idx]
S_best, U_S_best, _, _ = calc_fit_svd(Y, Phi_best) # s_vals and initial_err not needed here
U          = U_S_best                              # This is U from SVD of S_best
B_opt_rows = (U[:, :nBasisOpt].T @ Phi_best)       # shape (k, Nz)

rng_plot        = np.random.default_rng(1)
idx_samples_plot = rng_plot.choice(Y.shape[0], size=1, replace=False)
ys_approx_plot   = []
for idx in idx_samples_plot:
    y_samp = Y[idx]
    coeffs_samp = np.linalg.lstsq(B_opt_rows.T, y_samp, rcond=None)[0]
    y_fit_samp  = coeffs_samp @ B_opt_rows
    error_samp  = np.linalg.norm(y_samp - y_fit_samp)
    print(f"Reconstruction error for sample {idx} using {nBasisOpt} optimal funcs from zc_best={zc_best:.2f}: {error_samp:.2e}")
    ys_approx_plot.append((y_fit_samp, zc_best, f'sample {idx}'))

plotFunctionApprox(z, Y[idx_samples_plot[0]], ys_approx_plot, bError=True, errMax=0.01, scMin=1.1)
plt.title(f"Reconstruction Error for Sample {idx_samples_plot[0]} using {nBasisOpt} Optimal Funcs from zc_best={zc_best:.2f}")

plt.show()