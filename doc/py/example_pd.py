import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from projective_dynamics import build_grid_2d, solve_pd

# Create a simple grid
nx, ny = 5, 5
bonds, points, masses, ks, fixed = build_grid_2d(nx, ny, m=1.0, m_end=1000.0, l=1.0, k=100.0)
velocity = np.zeros_like(points)

# Run simulation
new_points, new_velocity = solve_pd(points, velocity, bonds, masses, ks, 
                                  dt=0.01, n_iter=100, fixed_points=fixed)

# Visualize results
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot initial configuration
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='blue', label='Initial')
for bond in bonds:
    i, j = bond
    ax.plot([points[i, 0], points[j, 0]], 
            [points[i, 1], points[j, 1]], 
            [points[i, 2], points[j, 2]], 'b-', alpha=0.3)

# Plot final configuration
ax.scatter(new_points[:, 0], new_points[:, 1], new_points[:, 2], c='red', label='Final')
for bond in bonds:
    i, j = bond
    ax.plot([new_points[i, 0], new_points[j, 0]], 
            [new_points[i, 1], new_points[j, 1]], 
            [new_points[i, 2], new_points[j, 2]], 'r-', alpha=0.3)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
