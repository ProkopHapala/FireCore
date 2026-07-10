#!/usr/bin/env python3
"""
Demo: Fitting angularly-dependent 2D potential with radial basis functions
on two sites (host atom at origin + off-center electron pair).

=== PHYSICS MOTIVATION ===

In molecular modeling, the interaction between an atom and a localized electron
density (e.g. a lone pair, or a bond electron pair) is often described by a
potential V(r, theta) that depends on both the distance r and the angle theta
relative to some symmetry axis. This angular dependence arises from the
multipole structure of the electron density:

  - Monopole (charge):       isotropic, ~1/r or constant
  - Dipole (polar):          ~cos(theta), breaks forward/backward symmetry
  - Quadrupole:              ~P2(cos(theta)) = (3cos^2(theta)-1)/2,
                              breaks symmetry between axial and equatorial

A common model potential combines a radial envelope (e.g. Morse potential for
bonding interactions) with angular modulation from these multipoles:

  V(r, theta) = V_morse(r) * [1 + c1 * cos(theta) * f1(r) + c2 * P2(cos(theta)) * f2(r)]

where f1(r), f2(r) are radial decay functions (here: exponentials).

=== THE FITTING PROBLEM ===

We want to approximate V(r, theta) as a sum of RADIAL basis functions centered
on two sites:
  - "Host" atom at the origin (r=0)
  - "Electron pair" (epair) at position (L, 0) along the x-axis

The approximation is:

  V_fit(x, y) = sum_i  c_i^host  * phi_i^host( |r - r_host| )
              + sum_j  c_j^epair * phi_j^epair( |r - r_epair| )

Each basis function phi_i(r) depends ONLY on the radial distance from its site.
The angular dependence of V emerges from the SUPERPOSITION of two radially-
symmetric basis sets centered at different points. This is analogous to how
two overlapping s-orbital-like functions can produce a p-orbital-like angular
pattern through linear combination (LCAO idea).

The coefficients c_i are found by LINEAR LEAST SQUARES: we sample V on a 2D
grid, build a design matrix A where each column is a basis function evaluated
at all grid points, and solve  A @ c = V  in the least-squares sense.

=== BASIS TYPES ===

Four radial basis families are supported (selectable per site):

  1. Exponential:  phi_i(r) = exp(-alpha_i * r)
     - Decay rate controlled by alpha_i. Long tails for small alpha.
     - Simple but not compactly supported; can be ill-conditioned if
       exponents are too close together (near-linear dependence).

  2. Gaussian:     phi_i(r) = exp(-r^2 / sigma_i^2)
     - Faster decay than exponential. Very smooth, infinitely differentiable.
     - Natural choice for quantum chemistry (Gaussian basis sets).
     - Width sigma_i controls locality: small sigma = very localized.

  3. B-spline:     piecewise polynomial of degree p on a knot sequence.
     - Compact support: each basis function is nonzero only on p+1 knot spans.
     - Local control: changing one coefficient only affects a local region.
     - Used in density fitting, nuclear potentials, etc.
     - Parameters: knot spacing (resolution), Rcut (domain), degree (smoothness).

  4. Jacobi:       P_n^{alpha,beta}(2r/Rcut - 1), n = 0, 1, ..., N
     - Orthogonal polynomials on [0, Rcut] (mapped from [-1,1]).
     - Global support: each polynomial spans the entire domain.
     - Good for smooth functions; spectral convergence for analytic V.
     - alpha=beta=0 gives Legendre polynomials (simplest case).

=== EXCLUDED ZONE AND HALF-PLANE ===

Near the host atom (r < exclude_radius), the potential is singular or
physically meaningless for our basis expansion (the atom core region). We
exclude this zone from the fit by setting a boolean mask.

We also restrict to the half-plane x >= 0, since the epair sits on the +x
axis and we are interested in the interaction region on that side.

=== LEPAIR SCAN ===

We scan the epair distance L (called "Lepair") from 1.0 to 3.0 Angstrom.
For each L, we:
  1. Fit the basis coefficients on a LARGE grid (10 A half-plane)
  2. Evaluate and plot the fit on a SMALL grid (4 A zoomed view)
  3. Compute the residual = V_true - V_fit and its RMS

This shows how the quality of the radial-basis approximation depends on the
geometry of the two-site decomposition.

Usage examples:

  # Default: exp on host, gaussian on epair
  python3 demo_epair_basis.py

  # B-spline on host, Jacobi on epair
  python3 demo_epair_basis.py --host-basis bspline --epair-basis jacobi

  # Compare all 4 basis types, save to separate folders + summary plot
  python3 demo_epair_basis.py --compare-all --outdir output_compare

  # Custom exponents and widths
  python3 demo_epair_basis.py --host-basis exp --host-exponents 0.3 0.8 1.5 3.0  --epair-basis gaussian --epair-widths 0.2 0.5 1.0 2.0

"""

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os
from scipy.interpolate import BSpline
from scipy.special import eval_jacobi

# ============================================================================
# Basis function generators
# Each returns a callable: basis(r) -> array shape (len(r), n_basis)
#
# The design matrix A for the least-squares fit is built by evaluating each
# basis function at all grid points. A has shape (n_points, n_basis_total)
# where n_basis_total = n_host_basis + n_epair_basis.
#
# A good basis set should:
#   - Span the function space where V has significant amplitude
#   - Be well-conditioned (columns not nearly linearly dependent)
#   - Have enough functions to capture the angular structure emerging
#     from the two-site superposition
# ============================================================================

# --- Exponential basis ---
# phi_i(r) = exp(-alpha_i * r)
#
# These are Slater-type functions (STO-like). They have a cusp at r=0
# (nonzero derivative) which matches the behavior of Coulombic potentials.
# However, they lack compact support and the design matrix can become
# ill-conditioned if exponents are too similar (columns nearly proportional).
#
# Choosing a geometric progression of exponents (e.g. 0.5, 1.0, 1.5, 2.0, 3.0)
# gives a multi-resolution representation: small alpha captures long-range
# tails, large alpha captures the near-core region.

def make_exp_basis(exponents):
	"""exp(-alpha*r) basis with given exponents."""
	exponents = np.atleast_1d(np.asarray(exponents, dtype=float))
	def basis(r):
		return np.stack([np.exp(-a * r) for a in exponents], axis=-1)
	return basis

# --- Gaussian basis ---
# phi_i(r) = exp(-r^2 / sigma_i^2)
#
# Gaussian-type functions (GTO-like) are the workhorse of quantum chemistry.
# They are smooth (C^infinity), decay faster than exponentials, and have the
# key computational advantage that products of Gaussians are also Gaussians
# (the product of two Gaussians centered at different points is another
# Gaussian centered at a weighted average point).
#
# Unlike exponentials, Gaussians have zero derivative at r=0 (flat top),
# which means they don't naturally capture Coulomb cusps. But for fitting
# smooth effective potentials, this is usually not a problem.
#
# Width sigma controls the spatial extent: sigma=0.3 A is very localized
# (captures near-site structure), sigma=3.0 A is broad (captures long-range).

def make_gaussian_basis(widths):
	"""exp(-r^2/sigma^2) basis with given widths."""
	widths = np.atleast_1d(np.asarray(widths, dtype=float))
	def basis(r):
		return np.stack([np.exp(-r**2 / s**2) for s in widths], axis=-1)
	return basis

# --- B-spline basis ---
# B-splines are piecewise polynomials defined by a knot sequence and degree.
#
# Key properties:
#   - COMPACT SUPPORT: each basis function B_i(r) is nonzero only on the
#     interval [t_i, t_{i+p+1}] where t_i are knots and p is the degree.
#     This means the design matrix is BANDED (sparse), and changing one
#     coefficient only affects the fit locally.
#   - PARTITION OF UNITY: sum_i B_i(r) = 1 for all r in [0, Rcut].
#     This means a constant function can be represented exactly.
#   - SMOOTHNESS: B-splines of degree p have C^{p-1} continuity at knots.
#     Cubic (p=3) gives C^2 continuity, which is sufficient for smooth
#     potentials and ensures continuous forces (derivative of V).
#
# The knot sequence is "clamped" (open): the first and last knots are
# repeated p+1 times so that the spline interpolates the endpoints.
# This is standard for bounded-domain fitting.
#
# Parameters:
#   knot_spacing = 0.5 A  -> resolution of local features
#   Rcut = 6.0 A         -> domain [0, Rcut]
#   degree = 3 (cubic)   -> smoothness C^2
#
# With knot_spacing=0.5 and Rcut=6.0, we get 13 interior knots, giving
# n_basis = 13 + 2*3 - 3 - 1 = 15 basis functions (for clamped cubic).

def make_bspline_basis(knot_spacing=0.5, Rcut=6.0, degree=3):
	"""Uniform B-spline basis on [0, Rcut] with clamped boundary knots."""
	interior = np.arange(0, Rcut + 1e-10, knot_spacing)
	knots = np.r_[[interior[0]]*degree, interior, [interior[-1]]*degree]
	n_basis = len(knots) - degree - 1
	def basis(r):
		r_c = np.clip(r, 0, Rcut)
		cols = []
		for i in range(n_basis):
			c = np.zeros(n_basis)
			c[i] = 1.0
			sp = BSpline(knots, c, degree, extrapolate=False)
			vals = np.nan_to_num(sp(r_c), nan=0.0)
			cols.append(vals)
		return np.stack(cols, axis=-1)
	return basis

# --- Jacobi polynomial basis ---
# P_n^{alpha,beta}(x) where x = 2r/Rcut - 1 maps [0, Rcut] -> [-1, 1]
#
# Jacobi polynomials are a family of ORTHOGONAL polynomials on [-1, 1] with
# weight function w(x) = (1-x)^alpha * (1+x)^beta.
#
# With alpha=beta=0, the weight is uniform and we get LEGENDRE polynomials:
#   P_0(x) = 1
#   P_1(x) = x
#   P_2(x) = (3x^2 - 1)/2
#   P_3(x) = (5x^3 - 3x)/2
#   ...
#
# Key properties:
#   - GLOBAL SUPPORT: each polynomial is nonzero over the entire [0, Rcut].
#     This means the design matrix is DENSE (no sparsity).
#   - ORTHOGONALITY: <P_m, P_n> = 0 for m != n (with appropriate weight).
#     This gives excellent conditioning of the normal equations.
#   - SPECTRAL CONVERGENCE: for smooth (analytic) V, the coefficients decay
#     EXPONENTIALLY with degree n. This means a small number of polynomials
#     can represent very smooth functions with high accuracy.
#   - However, for functions with localized features or discontinuities,
#     spectral methods suffer from GIBBS oscillations (ringing artifacts).
#
# max_degree=6 gives 7 basis functions (P_0 through P_6), which can represent
# polynomials up to degree 6 exactly. For our Morse x dipole/quadrupole
# potential, this should capture most of the structure if Rcut is large enough.

def make_jacobi_basis(alpha=0.0, beta=0.0, max_degree=6, Rcut=6.0):
	"""Jacobi polynomial basis on [0, Rcut] mapped to [-1, 1]."""
	def basis(r):
		x = 2.0 * np.clip(r, 0, Rcut) / Rcut - 1.0
		cols = [eval_jacobi(n, alpha, beta, x) for n in range(max_degree + 1)]
		return np.stack(cols, axis=-1)
	return basis

BASIS_FACTORIES = {
	'exp':      make_exp_basis,
	'gaussian': make_gaussian_basis,
	'bspline':  make_bspline_basis,
	'jacobi':   make_jacobi_basis,
}

def build_basis(kind, params):
	if kind not in BASIS_FACTORIES:
		raise ValueError(f"Unknown basis '{kind}'. Choose from {list(BASIS_FACTORIES.keys())}")
	return BASIS_FACTORIES[kind](**params)

# ============================================================================
# Test potential: Morse × (dipole + quadrupole) modulation in x-y plane
#
# The MORSE potential is a standard model for diatomic interactions:
#   V_morse(r) = D * (1 - exp(-a*(r - r0)))^2 - D
#
#   D  = well depth (eV)        -- how strongly the atoms bind
#   a  = width parameter (1/A)  -- controls how narrow the well is
#   r0 = equilibrium distance   -- where V_morse = -D (minimum)
#
# As r -> infinity: V_morse -> 0 (dissociation limit)
# At r = r0:        V_morse = -D (bond minimum)
# As r -> 0:        V_morse -> +infinity (Pauli repulsion, steep wall)
#
# The ANGULAR MODULATION introduces directionality via multipole expansion:
#
#   cos(theta) = (r_vec . d_hat) / |r_vec|
#
# where d_hat is the (normalized) dipole direction. This is the angle between
# the position vector r and the dipole axis.
#
#   - Dipole term:    c1 * cos(theta) * exp(-alpha1 * r)
#     This breaks the forward/backward symmetry. When the dipole points along
#     +x (default), the potential is stronger on the +x side and weaker on -x.
#     The exponential exp(-alpha1*r) makes this modulation decay with distance,
#     so the potential becomes isotropic far from the origin.
#
#   - Quadrupole term: c2 * P2(cos(theta)) * exp(-alpha2 * r)
#     P2(cos(theta)) = (3*cos^2(theta) - 1) / 2
#     This is the second Legendre polynomial. It distinguishes:
#       - Axial directions (theta=0 or pi):  P2 = +1 (enhanced)
#       - Equatorial directions (theta=pi/2): P2 = -1/2 (suppressed)
#     The quadrupole creates a "dumbbell" shaped anisotropy.
#
# The full potential: V = V_morse(r) * [1 + dipole_mod + quadrupole_mod]
# The factor [1 + ...] ensures the potential reduces to V_morse when the
# modulation strengths c1, c2 -> 0.
# ============================================================================

def test_potential(X, Y, D=1.0, a=1.0, r0=2.0,
                   c1=0.3, c2=0.15, alpha1=1.0, alpha2=0.8,
                   dip_dir=(1.0, 0.0)):
	"""
	V(r,theta) = V_morse(r) * (1 + c1*cos_theta*exp(-alpha1*r) + c2*P2(cos_theta)*exp(-alpha2*r))
	cos_theta = (r . d_hat) / r  using normalized direction (x/r, y/r) dot dip_dir.
	"""
	r = np.sqrt(X**2 + Y**2 + 1e-12)
	nx, ny = dip_dir
	norm = np.sqrt(nx**2 + ny**2) + 1e-12
	nx, ny = nx / norm, ny / norm
	cos_theta = (X * nx + Y * ny) / r  # normalized direction dot dipole_dir
	P2 = 0.5 * (3 * cos_theta**2 - 1)
	V_morse = D * (1 - np.exp(-a * (r - r0)))**2 - D
	mod = 1.0 + c1 * cos_theta * np.exp(-alpha1 * r) + c2 * P2 * np.exp(-alpha2 * r)
	return V_morse * mod

# ============================================================================
# Fitting: linear least squares with radial basis on multiple sites
#
# The model is LINEAR in the coefficients c_i:
#
#   V_fit(x, y) = sum_i c_i * phi_i(|r - r_site_i|)
#
# This can be written as a matrix equation:  A @ c = V
#
# where:
#   A = design matrix, shape (n_points, n_basis_total)
#     - Each column j of A is the j-th basis function evaluated at all grid points
#     - Columns are grouped by site: first n_host columns for host basis,
#       then n_epair columns for epair basis
#   c = coefficient vector, shape (n_basis_total,)
#   V = potential values at grid points, shape (n_points,)
#
# Since we have more grid points than basis functions (overdetermined system),
# we solve in the LEAST SQUARES sense: find c that minimizes ||A@c - V||^2.
#
# This is solved via SVD (Singular Value Decomposition) internally by
# np.linalg.lstsq. The SVD approach is numerically stable even if A is
# rank-deficient (some basis functions nearly linearly dependent).
#
# The RANK of A tells us how many independent basis functions we effectively
# have. If rank < n_basis_total, some functions are redundant.
#
# MASKING: grid points inside the excluded zone (r < exclude_radius) are
# removed from both A and V before fitting. This prevents the singular
# near-core region from biasing the fit. We use NaN in the output arrays
# for masked-out points so they appear blank in plots.
# ============================================================================

def fit_potential(X, Y, V, site_positions, basis_funcs, mask=None):
	"""
	Fit V(x,y) ~= sum_sites sum_basis c_i * phi_i(|r - r_site|)
	via linear least squares. If mask is given (boolean array same shape as V),
	only masked-in points are used for fitting.
	Returns: coeffs_per_site (list of arrays), V_fit (full grid, NaN outside mask), residual (NaN outside mask).
	"""
	if mask is None:
		mask = np.ones_like(V, dtype=bool)
	cols = []
	n_per_site = []
	# Build design matrix: for each site, evaluate all its basis functions
	# at distances |r - r_site| from that site. Each site contributes
	# n_basis_site columns to A.
	for pos, basis_fn in zip(site_positions, basis_funcs):
		r_site = np.sqrt((X - pos[0])**2 + (Y - pos[1])**2)  # distance from this site
		B = basis_fn(r_site.ravel())  # shape (n_points, n_basis_site)
		cols.append(B)
		n_per_site.append(B.shape[1])
	A = np.hstack(cols)  # shape (n_points, n_basis_total) = [host_cols | epair_cols]
	# Apply mask: remove excluded-zone points from both A and V
	A_fit = A[mask.ravel()]
	V_fit_flat = V.ravel()[mask.ravel()]
	# Solve least-squares: minimize ||A_fit @ c - V_fit_flat||^2
	# Uses SVD internally for numerical stability
	coeffs_all, _, rank, sv = np.linalg.lstsq(A_fit, V_fit_flat, rcond=None)
	# Evaluate fit and residual on masked points only; NaN elsewhere for plotting
	V_fit = np.full_like(V, np.nan)
	residual = np.full_like(V, np.nan)
	V_fit[mask] = (A_fit @ coeffs_all)
	residual[mask] = V_fit_flat - (A_fit @ coeffs_all)
	# Split coefficient vector back into per-site arrays
	coeffs_per_site = []
	idx = 0
	for n in n_per_site:
		coeffs_per_site.append(coeffs_all[idx:idx + n])
		idx += n
	return coeffs_per_site, V_fit, residual

# ============================================================================
# Plotting helpers
# ============================================================================

def plot_2d_map(ax, x, y, Z, title, cmap='RdBu_r', vmin=None, vmax=None):
	im = ax.imshow(Z, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap=cmap, aspect='equal',  vmin=vmin, vmax=vmax)
	ax.set_title(title)
	ax.set_xlabel('x (A)')
	ax.set_ylabel('y (A)')
	plt.colorbar(im, ax=ax, shrink=0.8)

def plot_radial_profile(ax, basis_fn, coeffs, r_max, title, color='C0'):
	r = np.linspace(0, r_max, 300)
	B = basis_fn(r)
	f = B @ coeffs
	ax.plot(r, f, color=color, lw=2)
	ax.set_title(title)
	ax.set_xlabel('r (A)')
	ax.set_ylabel('f(r) (eV)')
	ax.axhline(0, color='gray', ls='--', lw=0.5)
	ax.grid(True, alpha=0.3)

# ============================================================================
# CLI parameter extraction
# ============================================================================

def get_basis_config(args, prefix):
	kind = getattr(args, f'{prefix}_basis')
	if kind == 'exp':
		return 'exp', {'exponents': getattr(args, f'{prefix}_exponents')}
	elif kind == 'gaussian':
		return 'gaussian', {'widths': getattr(args, f'{prefix}_widths')}
	elif kind == 'bspline':
		return 'bspline', {'knot_spacing': getattr(args, f'{prefix}_bspline_knot'), 'Rcut': getattr(args, f'{prefix}_bspline_rcut')}
	elif kind == 'jacobi':
		return 'jacobi', {'max_degree': getattr(args, f'{prefix}_jacobi_maxdeg'), 'Rcut': getattr(args, f'{prefix}_jacobi_rcut')}
	else:
		raise ValueError(f"Unknown basis: {kind}")

# ============================================================================
# Single-basis scan (host+epair pair)
# ============================================================================

def get_default_params_for_kind(kind, prefix):
	"""Return default parameters dict for a basis kind."""
	if kind == 'exp':
		return {'exponents': [0.5, 1.0, 1.5, 2.0, 3.0]}
	elif kind == 'gaussian':
		if prefix == 'host':
			return {'widths': [0.5, 1.0, 1.5, 2.0, 3.0]}
		else:
			return {'widths': [0.3, 0.6, 1.0, 1.5, 2.5]}
	elif kind == 'bspline':
		return {'knot_spacing': 0.5, 'Rcut': 6.0}
	elif kind == 'jacobi':
		return {'max_degree': 6, 'Rcut': 6.0}


def run_basis_scan(outdir, host_kind, host_params, epair_kind, epair_params,
                   X_fit, Y_fit, V_fit_grid, mask_fit,
                   X_plot, Y_plot, V_plot, mask_plot, V_plot_masked,
                   Lepairs, potential_params, plot_range, show=False):
	"""Run the full Lepair scan for a single (host, epair) basis pair and save figures."""
	os.makedirs(outdir, exist_ok=True)
	host_basis = build_basis(host_kind, host_params)
	epair_basis = build_basis(epair_kind, epair_params)

	x_plot = X_plot[0, :]
	y_plot = Y_plot[:, 0]
	rms_residuals = []

	print(f"Host basis:  {host_kind} ({host_params})")
	print(f"Epair basis: {epair_kind} ({epair_params})")

	for iL, L in enumerate(Lepairs):
		epair_pos = np.array([L, 0.0])
		site_positions = [np.array([0.0, 0.0]), epair_pos]
		basis_funcs = [host_basis, epair_basis]

		coeffs_list, _, _ = fit_potential(X_fit, Y_fit, V_fit_grid, site_positions, basis_funcs, mask=mask_fit)

		# Evaluate the fitted model on the PLOT grid (not the fit grid).
		# We rebuild the design matrix A_plot using the plot-grid coordinates,
		# then multiply by the coefficients found from the fit grid.
		# This shows how well the fit generalizes to the zoomed-in region.
		# The residual = V_true(plot) - V_fit(plot) measures the quality
		# of the approximation in the region of physical interest.
		cols_plot = []
		for pos, basis_fn in zip(site_positions, basis_funcs):
			r_s = np.sqrt((X_plot - pos[0])**2 + (Y_plot - pos[1])**2)
			cols_plot.append(basis_fn(r_s.ravel()))
		A_plot = np.hstack(cols_plot)
		coeffs_flat = np.concatenate(coeffs_list)
		V_fit_plot = (A_plot @ coeffs_flat).reshape(V_plot.shape)
		V_fit_plot[~mask_plot] = np.nan
		residual_plot = V_plot - V_fit_plot
		residual_plot[~mask_plot] = np.nan
		rms = np.sqrt(np.nanmean(residual_plot**2))
		rms_residuals.append(rms)

		# Figure: 2x2 layout
		fig, axes = plt.subplots(2, 2, figsize=(14, 12))

		# Top-left: Fitted potential
		plot_2d_map(axes[0, 0], x_plot, y_plot, V_fit_plot, f'Fitted V (L={L:.2f} A)',cmap='viridis', vmin=np.nanmin(V_plot_masked), vmax=np.nanmax(V_plot_masked))
		axes[0, 0].plot(0, 0, 'k+', markersize=12, markeredgewidth=2, label='host')
		axes[0, 0].plot(L, 0, 'rx', markersize=10, markeredgewidth=2, label='epair')
		axes[0, 0].legend()

		# Top-right: Residual
		r_max_res = np.nanmax(np.abs(residual_plot))
		plot_2d_map(axes[0, 1], x_plot, y_plot, residual_plot,  f'Residual (RMS={rms:.4f} eV)',cmap='RdBu_r', vmin=-r_max_res, vmax=r_max_res)
		axes[0, 1].plot(0, 0, 'k+', markersize=12, markeredgewidth=2)
		axes[0, 1].plot(L, 0, 'rx', markersize=10, markeredgewidth=2)

		# Bottom-left: Host radial profile
		plot_radial_profile(axes[1, 0], host_basis, coeffs_list[0], plot_range,  f'Host radial ({host_kind})', color='C0')

		# Bottom-right: Epair radial profile
		plot_radial_profile(axes[1, 1], epair_basis, coeffs_list[1], plot_range,  f'Epair radial ({epair_kind})', color='C1')

		plt.suptitle(f'Lepair = {L:.2f} A', fontsize=14, fontweight='bold')
		plt.tight_layout()
		plt.savefig(os.path.join(outdir, f'Lepair_{L:.2f}.png'), dpi=150)
		if show:
			plt.show()
		plt.close(fig)

		print(f"  L={L:.2f} A: RMS residual = {rms:.6f} eV")

	# --- Basis functions figure ---
	r_plot = np.linspace(0, plot_range, 300)
	fig, axes = plt.subplots(1, 2, figsize=(14, 5))
	for ax, basis_fn, kind, label, color in [
		(axes[0], host_basis, host_kind, 'Host', 'C0'),
		(axes[1], epair_basis, epair_kind, 'Epair', 'C1'),
	]:
		B = basis_fn(r_plot)
		for i in range(B.shape[1]):
			ax.plot(r_plot, B[:, i], lw=1.5, alpha=0.8, label=f'phi_{i}')
		ax.set_title(f'{label} basis ({kind}), {B.shape[1]} functions')
		ax.set_xlabel('r (A)')
		ax.set_ylabel('phi(r)')
		ax.axhline(0, color='gray', ls='--', lw=0.5)
		ax.grid(True, alpha=0.3)
		ax.legend(fontsize=8, ncol=2)
	plt.tight_layout()
	plt.savefig(os.path.join(outdir, 'basis_functions.png'), dpi=150)
	if show:
		plt.show()
	plt.close(fig)

	# --- Original potential figure ---
	fig, ax = plt.subplots(1, 1, figsize=(7, 6))
	plot_2d_map(ax, x_plot, y_plot, V_plot_masked, 'Test potential (Morse x dipole/quadrupole)')
	ax.plot(0, 0, 'k+', markersize=12, markeredgewidth=2, label='host')
	ax.legend()
	plt.tight_layout()
	plt.savefig(os.path.join(outdir, 'potential_original.png'), dpi=150)
	if show:
		plt.show()
	plt.close(fig)

	# --- Residual vs Lepair ---
	fig, ax = plt.subplots(1, 1, figsize=(8, 5))
	ax.plot(Lepairs, rms_residuals, 'o-', lw=2, markersize=8, color='C2')
	ax.set_xlabel('Lepair (A)')
	ax.set_ylabel('RMS residual (eV)')
	ax.set_title('Total residual vs Lepair')
	ax.grid(True, alpha=0.3)
	plt.tight_layout()
	plt.savefig(os.path.join(outdir, 'residual_vs_Lepair.png'), dpi=150)
	if show:
		plt.show()
	plt.close(fig)

	# Save data
	np.savez(os.path.join(outdir, 'results.npz'),   Lepairs=Lepairs, rms_residuals=np.array(rms_residuals))

	print(f"\nDone. {len(Lepairs)} figures + 2 summary plots saved to '{outdir}/'")
	return np.array(rms_residuals)


# ============================================================================
# Main
# ============================================================================

def parse_args():
	parser = ArgumentParser(description="Demo: fit 2D angular potential with radial basis on two sites")

	# Fit grid (large sampling box)
	parser.add_argument('--fit-range', type=float, default=10.0, help='Fit grid half-extent (A); sampled from -fit-range to +fit-range in both x and y')
	parser.add_argument('--fit-n', type=int, default=200, help='Fit grid points per axis')

	# Plot grid (zoomed view)
	parser.add_argument('--plot-range', type=float, default=4.0, help='Plot extent (A); x from 0 to range, y from -range to range')
	parser.add_argument('--plot-n', type=int, default=120, help='Plot grid points per axis')

	parser.add_argument('--exclude-radius', type=float, default=0.5, help='Excluded zone radius around origin (A), diameter=2*radius')

	# Morse potential
	parser.add_argument('--D', type=float, default=1.0, help='Morse well depth (eV)')
	parser.add_argument('--a', type=float, default=1.0, help='Morse width parameter')
	parser.add_argument('--r0', type=float, default=2.0, help='Morse equilibrium distance (A)')

	# Modulation
	parser.add_argument('--c1', type=float, default=0.3, help='Dipole modulation strength')
	parser.add_argument('--c2', type=float, default=0.15, help='Quadrupole modulation strength')
	parser.add_argument('--alpha1', type=float, default=1.0, help='Dipole exponential decay (1/A)')
	parser.add_argument('--alpha2', type=float, default=0.8, help='Quadrupole exponential decay (1/A)')
	parser.add_argument('--dip-x', type=float, default=1.0, help='Dipole direction x-component')
	parser.add_argument('--dip-y', type=float, default=0.0, help='Dipole direction y-component')

	# Lepair scan
	parser.add_argument('--Lepair-min', type=float, default=1.0, help='Min epair distance (A)')
	parser.add_argument('--Lepair-max', type=float, default=3.0, help='Max epair distance (A)')
	parser.add_argument('--Lepair-n', type=int, default=8, help='Number of Lepair values')

	# Host basis
	parser.add_argument('--host-basis', type=str, default='exp', choices=list(BASIS_FACTORIES.keys()))
	parser.add_argument('--host-exponents', type=float, nargs='+', default=[0.5, 1.0, 1.5, 2.0, 3.0])
	parser.add_argument('--host-widths', type=float, nargs='+', default=[0.5, 1.0, 1.5, 2.0, 3.0])
	parser.add_argument('--host-bspline-knot', type=float, default=0.5)
	parser.add_argument('--host-bspline-rcut', type=float, default=6.0)
	parser.add_argument('--host-jacobi-maxdeg', type=int, default=6)
	parser.add_argument('--host-jacobi-rcut', type=float, default=6.0)

	# Epair basis
	parser.add_argument('--epair-basis', type=str, default='gaussian', choices=list(BASIS_FACTORIES.keys()))
	parser.add_argument('--epair-exponents', type=float, nargs='+', default=[0.5, 1.0, 1.5, 2.0, 3.0])
	parser.add_argument('--epair-widths', type=float, nargs='+', default=[0.3, 0.6, 1.0, 1.5, 2.5])
	parser.add_argument('--epair-bspline-knot', type=float, default=0.5)
	parser.add_argument('--epair-bspline-rcut', type=float, default=6.0)
	parser.add_argument('--epair-jacobi-maxdeg', type=int, default=6)
	parser.add_argument('--epair-jacobi-rcut', type=float, default=6.0)

	# Output / modes
	parser.add_argument('--outdir', type=str, default='output')
	parser.add_argument('--show', action='store_true', help='Show plots interactively')
	parser.add_argument('--compare-all', action='store_true', help='Run all 4 basis types (same on host+epair) into separate subfolders and make a comparison plot')

	return parser.parse_args()


def main():
	args = parse_args()

	# === TWO-GRID STRATEGY ===
	# We use separate grids for FITTING and PLOTTING:
	#
	# - FIT grid: large (10 A half-plane, 200x200) to capture the full
	#   spatial extent of the potential. The Morse potential decays slowly
	#   (~exp(-2*a*r)), so we need a large domain to avoid boundary effects
	#   in the least-squares fit. If we fit only on a small grid, the basis
	#   functions may "overfit" the local region and produce artifacts
	#   outside the fitting domain.
	#
	# - PLOT grid: small (4 A half-plane, 120x120) for visualization.
	#   This zooms into the physically interesting region near the host
	#   atom and the epair site.
	#
	# The fit coefficients are solved on the fit grid, then the model is
	# EVALUATED on the plot grid for display. This is analogous to fitting
	# a function on a fine mesh but only displaying a cross-section.

	# --- Fit grid: large box, half-plane x>=0 ---
	# Half-plane x>=0 because the epair is on the +x axis; the interesting
	# physics (dipole asymmetry, two-site interference) is on this side.
	x_fit = np.linspace(0, args.fit_range, args.fit_n)
	y_fit = np.linspace(-args.fit_range, args.fit_range, args.fit_n)
	X_fit, Y_fit = np.meshgrid(x_fit, y_fit)

	potential_params = dict(D=args.D, a=args.a, r0=args.r0, c1=args.c1, c2=args.c2, alpha1=args.alpha1, alpha2=args.alpha2,  dip_dir=(args.dip_x, args.dip_y))

	V_fit_grid = test_potential(X_fit, Y_fit, **potential_params)

	r_fit = np.sqrt(X_fit**2 + Y_fit**2)
	# Mask: exclude points within exclude_radius of the origin.
	# The Morse potential diverges as r->0 (Pauli repulsion wall),
	# and the radial basis functions may not have enough flexibility
	# there. Excluding this zone prevents the steep wall from dominating
	# the least-squares fit at the expense of the softer angular region.
	mask_fit = r_fit > args.exclude_radius

	# --- Plot grid: zoomed half-plane x>=0 ---
	x_plot = np.linspace(0, args.plot_range, args.plot_n)
	y_plot = np.linspace(-args.plot_range, args.plot_range, args.plot_n)
	X_plot, Y_plot = np.meshgrid(x_plot, y_plot)

	V_plot = test_potential(X_plot, Y_plot, **potential_params)
	r_plot_grid = np.sqrt(X_plot**2 + Y_plot**2)
	mask_plot = r_plot_grid > args.exclude_radius
	V_plot_masked = V_plot.copy()
	V_plot_masked[~mask_plot] = np.nan  # NaN renders as blank in imshow

	# === LEPAIR SCAN ===
	# We scan the epair distance L from Lepair_min to Lepair_max.
	# For each L, the epair site is at (L, 0) on the x-axis.
	#
	# The physics: as L increases, the two basis centers (host at origin,
	# epair at (L,0)) move apart. When they are close (small L), their
	# basis functions overlap strongly, giving the fit more "resolution"
	# to capture angular structure. When they are far apart (large L),
	# the overlap decreases and the angular resolution degrades.
	#
	# The RMS residual vs Lepair curve shows this tradeoff:
	#   - Small L: good angular resolution but the epair basis may
	#     "compete" with the host basis in the core region
	#   - Large L: poor angular resolution, basis functions don't overlap
	#     enough to represent the angular modulation
	Lepairs = np.linspace(args.Lepair_min, args.Lepair_max, args.Lepair_n)

	if args.compare_all:
		results = {}
		for kind in ['exp', 'gaussian', 'bspline', 'jacobi']:
			outdir = os.path.join(args.outdir, f'basis_{kind}')
			host_params = get_default_params_for_kind(kind, 'host')
			epair_params = get_default_params_for_kind(kind, 'epair')
			print(f"\n===== Running {kind} basis in {outdir} =====")
			rms = run_basis_scan(outdir, kind, host_params, kind, epair_params,
			                     X_fit, Y_fit, V_fit_grid, mask_fit,
			                     X_plot, Y_plot, V_plot, mask_plot, V_plot_masked,
			                     Lepairs, potential_params, args.plot_range, args.show)
			results[kind] = rms

		# --- Comparison summary plot ---
		fig, ax = plt.subplots(1, 1, figsize=(8, 5))
		for kind, rms in results.items():
			ax.plot(Lepairs, rms, 'o-', lw=2, markersize=8, label=kind)
		ax.set_xlabel('Lepair (A)')
		ax.set_ylabel('RMS residual (eV)')
		ax.set_title('RMS residual vs Lepair: basis comparison')
		ax.legend(title='basis')
		ax.grid(True, alpha=0.3)
		plt.tight_layout()
		comp_path = os.path.join(args.outdir, 'comparison_Lepair.png')
		plt.savefig(comp_path, dpi=150)
		if args.show:
			plt.show()
		plt.close(fig)
		print(f"\nComparison plot saved to {comp_path}")

		np.savez(os.path.join(args.outdir, 'comparison_results.npz'),
		         Lepairs=Lepairs, results=results)
	else:
		host_kind, host_params = get_basis_config(args, 'host')
		epair_kind, epair_params = get_basis_config(args, 'epair')
		run_basis_scan(args.outdir, host_kind, host_params, epair_kind, epair_params,
		               X_fit, Y_fit, V_fit_grid, mask_fit,
		               X_plot, Y_plot, V_plot, mask_plot, V_plot_masked,
		               Lepairs, potential_params, args.plot_range, args.show)

if __name__ == '__main__':
	main()
