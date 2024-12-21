import numpy as np
import matplotlib.pyplot as plt
from projective_dynamics import build_grid_2d, solve_pd

def plot_truss(points1, bonds, c='k', label='', alpha=1.0  ):
    """Plot two configurations of the grid for comparison"""
    plt.plot(points1[:, 0], points1[:, 1], 'o', c=c, label=label)
    for bond in bonds:
        i, j = bond
        plt.plot([points1[i, 0], points1[j, 0]], [points1[i, 1], points1[j, 1]], '-', c=c, alpha=alpha)
    
    plt.grid(True)
    plt.axis('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

# Create a simple grid
nx, ny = 5, 5
bonds, points, masses, ks, fixed = build_grid_2d(nx, ny, m=1.0, m_end=1000.0, l=1.0, k=10000.0)
velocity = np.zeros_like(points)


# Visualize results
plt.figure(figsize=(10, 10))
plot_truss(points,     bonds, c='b', label='Initial Points' )
new_points, new_velocity = solve_pd(points, velocity, bonds, masses, ks, dt=0.001, n_iter=50, fixed_points=fixed, call_back=lambda x: plot_truss(x,bonds, c='k', alpha=0.1 ) )
plot_truss(new_points, bonds, c='r', label='Final Points'   )
plt.show()
