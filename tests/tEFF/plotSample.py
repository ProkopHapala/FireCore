import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff
from pyBall import eFF_terms as pyeff
# const_bohr = 0.5291772105638411

spin1 = 1.0
spin2 = -1.0

# Define a vector function (example: E, fx, fy, fz)
def example_vector_function(points):
    Fout = eff.sample_ee(points, 1, bEvalCoulomb=True, bEvalPauli=True)
    return Fout

nx = 100
X = np.linspace(-5, 5, nx)
RSs = np.array([[x, spin1, spin2] for x in X])

points = np.zeros( (len(X), 3) )
points[:,0] = np.linspace(-5, 5, 100)
points[:,1] = spin1
points[:,2] = spin2

FEout = example_vector_function( points)
print(FEout)

plt.plot(X, FEout[:,0], label="E")
plt.plot(X, FEout[:,1], label="fx")
plt.plot(X, FEout[:,2], label="fy")
plt.plot(X, FEout[:,3], label="fz")

plt.legend()
plt.title("Sample ee")
plt.show()
# # Deserialize the results back into 2D grids
# fx_grid = FEout[:,0]
# fy_grid = FEout[:,1]
# fz_grid = FEout[:,2]
# E_grid  = FEout[:,3]

# # ========== Plotting ==========

# # Plot the components using imshow
# fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# # Plot E
# ax = axes[0, 0]
# im = ax.imshow(E_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='viridis')
# ax.set_title('E')
# fig.colorbar(im, ax=ax)

# # Plot fx
# ax = axes[0, 1]
# im = ax.imshow(fx_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
# ax.set_title('fx')
# fig.colorbar(im, ax=ax)

# # Plot fy
# ax = axes[1, 0]
# im = ax.imshow(fy_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
# ax.set_title('fy')
# fig.colorbar(im, ax=ax)

# # Plot fz
# ax = axes[1, 1]
# im = ax.imshow(fz_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
# ax.set_title('fz')
# fig.colorbar(im, ax=ax)

# plt.tight_layout()
# plt.show()

# distance = np.linspace(0, 10, 1000)

# force = Fout[:, 1]  # Extract y values
# force2 = Fout[:, 2]
# force3 = Fout[:, 3]
# print(Fout)

# plt.plot(distance, force)
# plt.plot(distance, force2)
# plt.plot(distance, force3)

# plt.show()
