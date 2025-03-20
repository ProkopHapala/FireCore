import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff

# Create a 2D grid of r_ij (distance between electrons) and s_i (electron cloud radius)
# We'll fix s_j (radius of second electron) to 1.0 Å
nr, ns = 100, 100
r_range = np.linspace(0.1, 5.0, nr)  # Distance range in Bohr
s_range = np.linspace(0.1, 2.0, ns)  # Electron cloud radius range in Bohr
R, S = np.meshgrid(r_range, s_range)

# Fixed parameters
s_j = 1.0  # Fixed radius of second electron in Angstrom
spin = 1   # Same spin electrons (1 for parallel, -1 for antiparallel)

# Prepare input array for sample_ee (each point needs r_ij, s_i, and s_j)
RSs = np.zeros((nr * ns, 3), dtype=np.float64, order='C')
RSs[:, 0] = R.flatten()  # r_ij values (distance between electrons)
RSs[:, 1] = S.flatten()  # s_i values (radius of first electron)
RSs[:, 2] = s_j        # s_j value (radius of second electron) - fixed

# Prepare output array
FEout = np.zeros((nr * ns, 4), dtype=np.float64, order='C')  # Will contain [fx, fy, fz, E] for each point

# Default parameters
KRSrho = np.array([1.125, 0.9, -0.2], dtype=np.float64)

# Call sample_ee
eff.sample_ee(RSs, spin=spin, FEout=FEout, KRSrho=KRSrho, bEvalCoulomb=True, bEvalPauli=True)

# Extract energy values and reshape to 2D grid
E_grid = FEout[:, 3].reshape(nr, ns)  # Energy is in the 4th component (w) of Quat4d

# Create the plot
plt.figure(figsize=(10, 8))
plt.imshow(E_grid, origin='lower', extent=[r_range[0], r_range[-1], s_range[0], s_range[-1]],  aspect='auto')
plt.colorbar(label='Energy (Hartree)')
plt.xlabel('r_ij (Bohr)')
plt.ylabel('s_i (Bohr)')
plt.title(f'Electron-Electron Interaction Energy (s_j = {s_j} Å, spin = {spin})')

# Add contour lines for better visualization
plt.contour(R, S, E_grid, colors='white', alpha=0.3, levels=10)

plt.tight_layout()
plt.show()
