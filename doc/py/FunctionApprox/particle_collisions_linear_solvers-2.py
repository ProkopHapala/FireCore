import numpy as np
import matplotlib.pyplot as plt

def harmonic_potential(r, R_min, R_cut, k, E_min):
    """
    Computes the unsmoothed harmonic potential and force.

    For r < R_cut:
      E(r) = k*(r - R_min)**2 + E_min
      F(r) = -2*k*(r - R_min)

    Parameters:
      r     : numpy array of distances.
      R_min : position of potential minimum.
      R_cut : cutoff distance for the harmonic region.
      k     : stiffness.
      E_min : depth of the potential at R_min.

    Returns:
      (E, F) : Tuple of arrays for potential and force.
    """
    E = np.zeros_like(r, dtype=float)
    F = np.zeros_like(r, dtype=float)
    mask = r < R_cut
    E[mask] = k * (r[mask] - R_min)**2 + E_min
    F[mask] = -2 * k * (r[mask] - R_min)
    return E, F

def smoothing_potential(r, R_cut, R_cut2, k, R_min, E_min):
    """
    Computes the smoothing (quadratic) potential and force in the interval
    [R_cut, R_cut2] with the form:

       E(r) = k2*(r - R_cut2)**2 + delta
       F(r) = -2*k2*(r - R_cut2)

    where k2 is tuned to match the derivative of the harmonic potential
    at r = R_cut, and delta is chosen to match the energy at r = R_cut.

    Parameters:
      r      : numpy array of distances (assumed within [R_cut, R_cut2]).
      R_cut  : end of the harmonic region.
      R_cut2 : cutoff for the smoothing region.
      k      : stiffness of the harmonic potential.
      R_min  : location of the potential minimum.
      E_min  : potential value at R_min.

    Returns:
      (E_smooth, F_smooth): arrays of the smoothing potential and force.
    """
    # Derivative matching:
    # Harmonic derivative at R_cut: 2*k*(R_cut - R_min)
    # Smoothing derivative at R_cut: 2*k2*(R_cut - R_cut2)
    # --> k2 = k*(R_cut-R_min)/(R_cut-R_cut2)
    k2 = k * (R_cut - R_min) / (R_cut - R_cut2)

    # Energy matching at R_cut:
    # E_harm(R_cut) = k*(R_cut-R_min)**2 + E_min
    # E_smooth(R_cut) = k2*(R_cut-R_cut2)**2 + delta
    # --> delta = E_harm(R_cut) - k2*(R_cut-R_cut2)**2
    delta = k * (R_cut - R_min)**2 + E_min - k2 * (R_cut - R_cut2)**2

    E_smooth = k2 * (r - R_cut2)**2 + delta
    F_smooth = -2 * k2 * (r - R_cut2)
    return E_smooth, F_smooth

def composite_potential(r, R_min, R_cut, R_cut2, k, E_min):
    """
    Computes the composite potential and force.

    For r < R_cut: Uses the harmonic potential.
    For R_cut <= r < R_cut2: Uses the smoothing quadratic.
    For r >= R_cut2: Sets the potential to a constant (delta) and force to zero.

    Parameters:
      r      : numpy array of distances.
      R_min  : position of potential minimum.
      R_cut  : cutoff for the harmonic region.
      R_cut2 : cutoff for the smoothing region.
      k      : stiffness.
      E_min  : potential value at R_min.

    Returns:
      (E, F): arrays for the composite potential and force.
    """
    E = np.zeros_like(r, dtype=float)
    F = np.zeros_like(r, dtype=float)

    # Harmonic region: r < R_cut.
    mask1 = r < R_cut
    E_harm, F_harm = harmonic_potential(r, R_min, R_cut, k, E_min)
    E[mask1] = E_harm[mask1]
    F[mask1] = F_harm[mask1]

    # Smoothing region: R_cut <= r < R_cut2.
    mask2 = np.logical_and(r >= R_cut, r < R_cut2)
    if np.any(mask2):
        E_smooth, F_smooth = smoothing_potential(r[mask2], R_cut, R_cut2, k, R_min, E_min)
        E[mask2] = E_smooth
        F[mask2] = F_smooth

    # Beyond R_cut2: potential is constant and force is zero.
    mask3 = r >= R_cut2
    if np.any(mask3):
        # Calculate delta as above.
        k2 = k * (R_cut - R_min) / (R_cut - R_cut2)
        delta = k * (R_cut - R_min)**2 + E_min - k2 * (R_cut - R_cut2)**2
        E[mask3] = delta
        F[mask3] = 0.0

    return E, F

# ======== Main

# Set the parameters.
R_min  = 2.0   # Distance at which the potential is minimum.
R_cut  = 2.5   # Cutoff for the harmonic potential.
R_cut2 = 3.5   # Cutoff for the smoothing (anharmonic) region.
k      = 5.0   # Stiffness of the harmonic potential.
E_min  = -1.0  # Depth of the potential at R_min.

# Create an array of distances.
r = np.linspace(0, 5.0, 1000)

# Calculate the potentials and forces.
E_harm, F_harm = harmonic_potential(r, R_min, R_cut, k, E_min)
E_comp, F_comp = composite_potential(r, R_min, R_cut, R_cut2, k, E_min)

# Plot the results.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ax1.plot(r, E_harm, label="Harmonic (unsmoothed)", linestyle="--")
ax1.plot(r, E_comp, label="Composite (smoothened)", lw=3.0)
ax1.axvline(x=R_min,  color='b', linestyle='--', label="R_min")
ax1.axvline(x=R_cut,  color='r', linestyle='--', label="R_cut")
ax1.axvline(x=R_cut2, color='g', linestyle='--', label="R_cut2")
ax1.set_xlabel('Distance [A]')
ax1.set_ylabel('Energy [eV]')
ax1.set_title('Short-Range Potential')
ax1.legend()
ax1.grid(True)

ax2.plot(r, F_harm, label="Harmonic (unsmoothed)", linestyle="--")
ax2.plot(r, F_comp, label="Composite (smoothened)", lw=3.0)
ax2.axvline(x=R_min,  color='b', linestyle='--')
ax2.axvline(x=R_cut,  color='r', linestyle='--')
ax2.axvline(x=R_cut2, color='g', linestyle='--')
ax2.set_xlabel('Distance [A]')
ax2.set_ylabel('Force [eV/A]')
ax2.set_title('Force')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()