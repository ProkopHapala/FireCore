# Chebyshev-Accelerated Jacobi Solver for Projective Dynamics

## 1. Introduction

This document provides a complete implementation guide for the Chebyshev-accelerated Jacobi solver as described in [Wang 2015] (Section 2 of the source document). The method combines Jacobi iterations with Chebyshev semi-iterative acceleration for fast convergence in projective dynamics simulations.

## 2. Core Algorithm

## 2.1 Chebyshev Acceleration Formula

From Section 2.1, the definitive Chebyshev acceleration formula is:



lets express this in terms of displacement vectors $d_k$ by which the positions $q_k$ are changed in each iteration:

$ q_{k+1} = ω_{k+1} ( γ {\hat d}_{k+1} + d_k ) + q_{k−1} $

$\hat d_{k+1} = q̂_{k+1}- q_k$
$d_{k+1} = q_{k+1}- q_k$
$d_k     = q_k- q_{k-1}$


Elimineate $q_{k-1}$ by using 
$ q_{k-1} = q_k - d_k $

$ q_{k+1} = ω ( γ ( q̂_{k+1} - q_k ) + q_k - q_k - d_k ) + q_k - d_k $


$ q_{k+1} = ω  γ  q̂_{k+1} - ω  γ q_k  + ω q_k - ω q_k - ω d_k  + q_k - d_k $

$ q_{k+1} = ω  γ  q̂_{k+1} + ( 1 - ω  γ ) q_k         - (1+ω) d_k $





$ q_{k+1} = ω_{k+1} ( γ ( q̂_{k+1} - q_k ) + q_k - q_{k−1} ) + q_{k−1} $
$ q_{k+1} = ω_{k+1}  γ ( q̂_{k+1} - q_k ) + ω_{k+1} q_k - ω_{k+1} q_{k−1} + q_{k−1} $
$ q_{k+1} = ω_{k+1}    q̂_{k+1} -  ω_{k+1} ( 1  - γ ) q_k + (1 - ω_{k+1}) q_{k−1} $

$ q_{k+1} = ω ( γ ( q̂_{k+1} - q_k ) + q_k - q_{k−1} ) + q_{k−1} $
$ q_{k+1} = ω   γ ( q̂_{k+1} - q_k ) + ω q_k - ω q_{k−1} + q_{k−1} $
$ q_{k+1} = ω   q̂_{k+1}     -  ω ( 1  - γ ) q_k + (1 - ω) q_{k−1} $

where:
 - $x_{new}$​ is the result of one Jacobi iteration
 - $x_{k}$​ is the previous position
 - $ω_{k+1}$​ is the Chebyshev weight

### 2.2 Chebyshev Weights Calculation

The weights are calculated as:

```python
def calculate_omega(k, rho, prev_omega):
    if k < DELAY_START:  # typically 10 iterations
        return 1.0
    elif k == DELAY_START:
        return 2.0 / (2.0 - rho**2)
    else:
        return 4.0 / (4.0 - rho**2 * prev_omega)
```

To estimate the spectral radius (ρ), we use power iteration:


```python
def estimate_spectral_radius(self, A, x0):
        """Estimate spectral radius using power iteration"""
        x = x0.copy()
        x_prev = x0.copy()
        
        # Apply matrix twice
        x = A @ x
        x_next = A @ x
        
        # Estimate using Rayleigh quotient
        rho = np.sqrt(abs(np.dot(x_next, x) / np.dot(x, x)))
        return min(rho, 0.9992)  # Cap at 0.9992 as suggested in Section 2.2.2
```


### 2.3 Complete Implementation

```python
        
def solve(self, A, b, x0, max_iter=400, tol=1e-6):
    """
    Solve Ax = b using Chebyshev-accelerated Jacobi
    Args:
        A: System matrix
        b: Right-hand side vector
        x0: Initial guess
    """
    # Initialize
    x_prev = x0.copy()
    x_curr = x0.copy()
    
    # Estimate spectral radius
    rho = estimate_spectral_radius(A, x0)
    
    # Initialize omega
    omega = 1.0
    
    for k in range(max_iter):
        # Regular Jacobi iteration
        D_inv = 1.0 / np.diag(A)
        x_new = D_inv * (b - (A @ x_curr - np.diag(A) * x_curr))
        
        # Calculate new Chebyshev weight
        omega_next = self.calculate_omega(k, rho, omega)
        
        # Chebyshev acceleration
        x_accel = omega_next * (x_new - x_prev) + x_prev
        
        # Update positions
        x_prev = x_curr
        x_curr = x_accel
        omega = omega_next
        
        # Check convergence
        residual = np.linalg.norm(b - A @ x_curr)
        if residual < tol:
            break
            
    return x_curr
```

## 3. Important Implementation Notes

### 3.1 Spectral Radius (ρ)

From Section 2.2.2:

 - Use ρ ≈ 0.9992-0.9999 for most cases
 - Underestimating ρ is better than overestimating
 - Can be estimated using power iteration as shown above

### 3.2 Delayed Start

From Section 2.2.2 of the paper:

 - Delay Chebyshev acceleration for first S iterations (typically 10)
 - Set ω₁ = ω₂ = ... = ωₛ = 1 for these iterations
 - Helps prevent oscillation/divergence in early iterations

### 3.3 Convergence Properties

From Section 2.2.2:

 - Method converges faster in first few iterations
 - Convergence rate drops after initial phase
 - Oscillation indicates ρ is too large
 - Divergence is rare if ρ is properly tuned

4. Usage Example

```python
# Create solver
solver = ChebyshevAcceleratedSolver(delay_start=10)

# Setup system
A = create_system_matrix()  # Your system matrix
b = create_rhs_vector()     # Your right-hand side
x0 = np.zeros(n)           # Initial guess

# Solve system
x = solver.solve(A, b, x0, max_iter=400, tol=1e-6)
```

## 5. Advantages and Limitations

### 5.1 Advantages (Section 2.4.2)

 - Highly parallelizable
 - Small memory footprint
 - No need for expensive linear algebra libraries
 - More efficient than CG on GPU
 - Robust to small system changes

### 5.2 Limitations

 - May require underrelaxation for stability in some cases
 - Convergence affected by mesh quality
 - May need smaller timesteps for highly stiff systems

## References

Wang, H. (2015). A Chebyshev Semi-iterative Approach for Accelerating Projective and Position-based Dynamics. ACM Trans. Graph. (SIGGRAPH Asia) 34, 6, Article 246.