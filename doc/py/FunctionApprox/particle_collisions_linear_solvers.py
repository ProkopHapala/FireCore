import numpy as np
import matplotlib.pyplot as plt

def getSR_r2(r, R_min, R_cut):
  """
  Calculates the Lennard-Jones potential and force between spherical atoms.

  Args:
    r: A NumPy array of distances between atoms.
    R_min: The distance at which the potential is minimum.
    R_cut: The cutoff distance beyond which the potential is zero.

  Returns:
    A tuple containing the potential and force arrays.
  """
  # Calculate the potential using a shifted parabola.
  x       = r - R_min
  mask    = r < R_cut
  E       = np.zeros_like(r, dtype=float)
  E[mask] = (x**2)[mask] - (R_cut - R_min)**2

  # Calculate the force as the negative derivative of the potential.
  F       = np.zeros_like(r, dtype=float)
  F[mask] = -2 * x[mask]

  return E, F

def getSR_r2_smooth(r, R_min, R_cut, R_cut2):
  """
  Calculates the short-range potential and force with a smoother transition.

  Args:
    r: A NumPy array of distances between atoms.
    R_min: The distance at which the potential is minimum.
    R_cut: The distance at which the first parabola ends.
    R_cut2: The distance at which the potential becomes zero.

  Returns:
    A tuple containing the potential and force arrays.
  """
  # Calculate the potential using two shifted parabolas.
  x     = r - R_min
  mask1 = r < R_cut
  mask2 = np.logical_and(r >= R_cut, r < R_cut2)
  E = np.zeros_like(r, dtype=float)
  F = np.zeros_like(r, dtype=float)

  # First parabola (r < R_cut)
  E[mask1] = (x**2)[mask1]  # No shift yet

  # Second parabola (R_cut <= r < R_cut2)
  # Match derivative (force) at R_cut
  k1 = 2 * (R_cut - R_min)
  E[mask2]        = k1 * (r[mask2] - R_cut2)

  # Make the second parabola with negative curvature and go to zero at R_cut2
  k2        = -k1 / (2 * (R_cut2-R_cut))                                 # curvature must be negative
  E[mask2] = E[mask2] + k2*((r[mask2]-R_cut)**2) - k2*((R_cut2-R_cut)**2)
  # Enforce to be zero at R_cut2

  # Calculate the force as the negative derivative of the potential.
  F[mask1] = -2 * x[mask1]
  F[mask2] = -k1 - 2*k2*(r[mask2] - R_cut)
  # Shift the first parabola to match the value of the second parabola at R_cut
  C = (R_cut - R_min)**2 - k2*((R_cut2-R_cut)**2)                              # Find the constant C
  E[mask1] = E[mask1] - C                                               # Apply the shift

  return E, F

def getSR_r2_smoothing(r, R_min, R_cut, R_cut2):
  """
    Computes the smoothening potential and force for the transition region.

    Args:
        r (numpy.ndarray): Array of distances between atoms.
        R_min (float): Distance at which the potential is minimum.
        R_cut (float): End of the first parabola and start of the smoothing region.
        R_cut2 (float): Distance where the potential becomes zero.

    Returns:
        tuple: Arrays of the smoothening potential and force.
  """
  x = r - R_min
  mask1 = r < R_cut
  mask2 = np.logical_and(r >= R_cut, r < R_cut2)
  E = np.zeros_like(r, dtype=float)
  F = np.zeros_like(r, dtype=float)

  # Match derivative (force) at R_cut
  k1 = 2 * (R_cut - R_min)

  k2        = -k1 / (2 * (R_cut2-R_cut))                                 # curvature must be negative
  E[mask2] =  k1 * (r[mask2] - R_cut2)  + k2*((r[mask2]-R_cut)**2) - k2*((R_cut2-R_cut)**2)
  F[mask2] = -k1 - 2*k2*(r[mask2] - R_cut)
  # Shift the first parabola to match the value of the second parabola at R_cut
  C = k2*((R_cut2-R_cut)**2)                              # Find the constant C
  F[mask1] = 0
  E[mask1] = C
  return E,F

# ======== Main

# Set the parameters.
R_min  = 2.0  # Distance at which the potential is minimum.
R_cut  = 2.5  # Cutoff distance.
R_cut2 = 3.5

# Create an array of distances.
r = np.linspace(0, 5.0, 1000)

# Calculate the potential and force.
E,F   = getSR_r2          (r, R_min, R_cut)
E2,F2 = getSR_r2_smooth   (r, R_min, R_cut, R_cut2)
E3,F3 = getSR_r2_smoothing(r, R_min, R_cut, R_cut2)

# Plot the potential and force.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

ax1.plot(r, E , label="SR_x2" )
ax1.plot(r, E2, label="SR_x2_smooth", lw=3.0)
ax1.plot(r, E3, label="SR_x2_smoothening")
ax1.axvline(x=R_min,  color='b', linestyle='--')
ax1.axvline(x=R_cut,  color='r', linestyle='--')
ax1.axvline(x=R_cut2, color='g', linestyle='--')
ax1.set_xlabel('Distance [A]')
ax1.set_ylabel('Energy [eV]')
ax1.set_title('Short-Range Potential and Force')
ax1.legend()
ax1.grid(True)

ax2.plot(r, F , label="SR_x2")
ax2.plot(r, F2, label="SR_x2_smooth", lw=3.0)
ax2.plot(r, F3, label="SR_x2_smoothening")
ax2.axvline(x=R_min,  color='b', linestyle='--')
ax2.axvline(x=R_cut,  color='r', linestyle='--')
ax2.axvline(x=R_cut2, color='g', linestyle='--')
ax2.set_ylabel
ax2.set_xlabel('Distance [A]')
ax2.set_ylabel('Force [eV/A]')
ax2.grid(True)

plt.tight_layout()
plt.show()