import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import eFF as eff

# Create a 2D grid of r_ia (distance) and s_i (electron cloud radius)
nr, ns = 100, 100
r_range = np.linspace(0.1, 5.0, nr)  # Distance range in Bohr
s_range = np.linspace(0.1, 5.0, ns)  # Electron cloud radius range in Bohr
R, S = np.meshgrid(r_range, s_range)

# Prepare input array for sample_EA (each point needs r_ia and s_i)
RSs = np.zeros((nr * ns, 2))
RSs[:, 0] = R.flatten()  # r_ia values
RSs[:, 1] = S.flatten()  # s_i values

# Prepare output array - must be 2D array for ctypes interface
FEout = np.zeros((nr * ns, 3), dtype=np.float64, order='C')  # Will contain [fx, fy, E] for each point

# Default parameters
KRSrho = np.array([1.125, 0.9, -0.2], dtype=np.float64)
aPar = np.array([4.0, 0.1, 0.1, 2.0], dtype=np.float64)

# Call sample_EA with explicit array passing
#eff.sample_EA(RSs, FEout=FEout, KRSrho=KRSrho, aPar=aPar, bEvalAECoulomb=True, bCoreCoul=True, bEvalAEPauli=True)
eff.sample_EA(RSs, FEout=FEout, KRSrho=KRSrho, aPar=aPar, bEvalAECoulomb=False, bCoreCoul=False, bEvalAEPauli=True)
#eff.sample_EA(RSs, FEout=FEout, KRSrho=KRSrho, aPar=aPar, bEvalAECoulomb=True, bCoreCoul=True, bEvalAEPauli=True)

# Extract energy values and reshape to 2D grid
E_grid = FEout[:, 2].reshape(nr, ns)

# Create the plot
plt.figure(figsize=(10, 8))
plt.imshow(E_grid, origin='lower', extent=[r_range[0], r_range[-1], s_range[0], s_range[-1]], vmax=100.0)
plt.colorbar(label='Energy (Hartree)')
plt.xlabel('r_ia (Bohr)')
plt.ylabel('s_i (Bohr)')
plt.title('Energy Landscape of Electron-Nuclei System')

# Add contour lines for better visualization
plt.contour(R, S, E_grid, colors='white', alpha=0.3, levels=10)



plt.tight_layout()
plt.show()