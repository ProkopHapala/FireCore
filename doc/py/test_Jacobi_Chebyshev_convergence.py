import numpy as np
import matplotlib.pyplot as plt
from projective_dynamics import solve_pd, build_grid_2d, make_pd_matrix, make_pd_rhs
from projective_dynamics_iterative import jacobi_iteration, calculate_omega, estimate_spectral_radius

def solve_direct(A, b):
    """Solve system directly using numpy"""
    return np.linalg.solve(A, b)

def solve_jacobi(A, b, x0, n_iter, fixed_points=None):
    """Solve system using Jacobi iteration"""
    x = x0.copy()
    errors = []
    x_direct = solve_direct(A, b)
    
    for _ in range(n_iter):
        x = jacobi_iteration(A, b, x, fixed_points)
        error = np.linalg.norm(x - x_direct) / np.linalg.norm(x_direct)
        errors.append(error)
    
    return x, errors

def solve_chebyshev(A, b, x0, n_iter, fixed_points=None):
    """Solve system using Chebyshev-accelerated Jacobi iteration"""
    x = x0.copy()
    x_prev = x0.copy()
    errors = []
    x_direct = solve_direct(A, b)
    
    # Estimate spectral radius
    rho = estimate_spectral_radius(A, np.random.rand(len(x0)))
    omega = 1.0
    
    for k in range(n_iter):
        # Perform Jacobi iteration
        x_new = jacobi_iteration(A, b, x, fixed_points)
        
        # Apply Chebyshev acceleration
        omega_next = calculate_omega(k, rho, omega)
        x_accel = omega_next * (x_new - x) + x
        
        # Update for next iteration
        x_prev = x.copy()
        x = x_accel
        omega = omega_next
        
        # Calculate error
        error = np.linalg.norm(x - x_direct) / np.linalg.norm(x_direct)
        errors.append(error)
    
    return x, errors

def test_convergence():
    # Create test problem
    nx, ny = 5, 5
    bonds, points, masses, ks, fixed = build_grid_2d(nx, ny)
    dt = 0.1
    
    # Build system matrix and RHS
    neighbs = [[] for _ in range(len(points))]
    for i, (i_, j_) in enumerate(bonds):
        neighbs[i_].append(i)
        neighbs[j_].append(i)
    
    A, Mt = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Create a random initial state and velocity for testing
    velocity = np.random.randn(len(points), 3) * 0.1
    points_pred = points + velocity * dt
    
    # Get RHS for one coordinate (x-coordinate)
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    b = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, points_pred)
    
    # Test both methods for x-coordinate
    n_iter = 100
    x0 = points_pred[:, 0].copy()  # Initial guess
    
    # Solve using both methods
    _, errors_jacobi = solve_jacobi(A, b[:, 0], x0, n_iter, fixed)
    _, errors_cheby = solve_chebyshev(A, b[:, 0], x0, n_iter, fixed)
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.semilogy(range(1, n_iter + 1), errors_jacobi, 'b-', label='Jacobi')
    plt.semilogy(range(1, n_iter + 1), errors_cheby, 'r-', label='Chebyshev')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Relative Error')
    plt.title('Convergence Comparison: Jacobi vs Chebyshev-accelerated Jacobi')
    plt.legend()
    plt.savefig('convergence_comparison.png')
    plt.close()

if __name__ == "__main__":
    test_convergence()
    print("Test completed. Results saved in 'convergence_comparison.png'")
