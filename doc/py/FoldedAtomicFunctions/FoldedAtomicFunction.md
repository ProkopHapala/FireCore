# Surface Potential Fitting Framework: Design Document

## Executive Summary

This document describes a novel computational framework for efficiently representing non-covalent interactions between molecules and crystal surfaces. The method addresses the computational bottleneck of repeatedly calculating long-range periodic interactions by pre-computing and compressing the interaction potential into a compact, localized basis set representation.

## 1. Background and Motivation

### 1.1 Problem Statement

Non-covalent interactions (electrostatic, van der Waals) between molecules and crystal surfaces are fundamental to many processes including:
- Molecular adsorption and desorption
- Catalytic reactions at surfaces  
- Crystal growth and dissolution
- Drug-surface interactions

The computational challenge arises from the need to account for the periodic nature of crystal surfaces, requiring summation over many periodic images to achieve convergence. Direct calculation scales poorly with system size and becomes prohibitively expensive for large-scale molecular dynamics or Monte Carlo simulations.

### 1.2 Current Limitations

**Direct Summation Approach:**
- Requires evaluation of interactions with hundreds/thousands of periodic images
- Computational cost scales as O(N_molecules × N_images × N_atoms_per_cell)
- Becomes inefficient for systems with fewer than ~1000 atoms
- Memory intensive for storing interaction matrices

**Existing Methods:**
- Ewald summation: Complex implementation, still computationally expensive
- Cubic B-spline interpolation: Overhead not justified for smaller systems
- Lookup tables: Memory intensive, limited accuracy

## 2. Proposed Method: Basis Set Compression

### 2.1 Core Concept

The method implements a two-stage approach:

1. **Pre-computation Stage** (expensive, done once):
   - Calculate full periodic potential on a dense grid
   - Include contributions from many periodic images
   - Store as reference potential field

2. **Compression Stage** (moderate cost, done once):
   - Fit the grid-based potential to a compact set of basis functions
   - Basis functions are localized within the unit cell
   - Store only the fitted coefficients

3. **Evaluation Stage** (very fast, done many times):
   - Evaluate potential at any point using linear combination of basis functions
   - No periodic image summation required

### 2.2 Mathematical Framework

The periodic potential at position **r** is:
```
V(r) = Σ_n Σ_i V_i(|r - r_i - n*L|)
```
where:
- `V_i` is the potential of atom i
- `r_i` is the position of atom i in the unit cell
- `n` is the periodic image vector
- `L` is the lattice vector

This is approximated as:
```
V(r) ≈ Σ_j c_j φ_j(r)
```
where:
- `φ_j(r)` are basis functions localized in the unit cell
- `c_j` are fitted coefficients

### 2.3 Advantages

- **Speed**: O(N_basis) evaluation vs O(N_images × N_atoms)
- **Memory**: Store only coefficients vs full interaction matrices
- **Accuracy**: Controlled by basis set size and grid resolution
- **Flexibility**: Easy to add new potential types or basis functions
- **Reusability**: Coefficients computed once, used many times

## 3. Technical Design

### 3.1 System Architecture

```
┌─────────────────┐    ┌──────────────────┐    ┌─────────────────┐
│   Grid Manager  │    │ Potential Calc   │    │ Basis Functions │
│                 │    │                  │    │                 │
│ • 2D Grid Setup │    │ • Coulomb 1/r    │    │ • Plane Waves   │
│ • Coordinates   │ ───│ • Morse Potential│───▶│ • Exponentials  │
│ • Boundaries    │    │ • Periodic Images│    │ • PBC Handling  │
└─────────────────┘    └──────────────────┘    └─────────────────┘
         │                        │                        │
         └────────────────────────┼────────────────────────┘
                                  ▼
                      ┌─────────────────┐    ┌─────────────────┐
                      │ Potential Fitter│    │   Visualizer    │
                      │                 │    │                 │
                      │ • Least Squares │    │ • Comparisons   │
                      │ • Regularization│    │ • Basis Plots   │
                      │ • RMSE Analysis │    │ • Error Analysis│
                      └─────────────────┘    └─────────────────┘
```

### 3.2 Module Specifications

#### 3.2.1 GridManager
**Purpose**: Manage 2D spatial grid and coordinate transformations

**Key Features**:
- Rectangular grid spanning unit cell (10Å × 10Å)
- Grid spacing: 0.2Å (50×50 points)
- Periodic boundary conditions in x-direction
- Non-periodic in z-direction (surface → vacuum)

**Interface**:
```python
class GridManager:
    def __init__(cell_x, cell_z, grid_step)
    def get_coordinates() -> (X, Z)  # 2D meshgrids
    def flatten_grid(grid_2d) -> array_1d
    def reshape_grid(array_1d) -> grid_2d
```

#### 3.2.2 PotentialCalculator
**Purpose**: Compute reference potentials with periodic boundary conditions

**Supported Potentials**:
- Coulomb: `V(r) = k_e * q / r`
- Morse: `V(r) = D * (exp(-2a(r-r0)) - 2*exp(-a(r-r0)))`

**Periodic Image Handling**:
- Sum over periodic images in x-direction: `x_image = x_atom + n*L_x`
- Configurable cutoff for convergence
- Efficient numpy broadcasting for vectorized calculation

**Interface**:
```python
class PotentialCalculator:
    def add_atom(x, z, charge, morse_params)
    def calculate_coulomb_periodic(n_images) -> potential_2d
    def calculate_morse_periodic(n_images) -> potential_2d
```

#### 3.2.3 BasisFunctions
**Purpose**: Generate and manage basis function sets

**Current Implementation**:
- **X-direction**: Cosine plane waves `cos(2πnx/L)` for n=1,2,...,8
  - Automatically satisfies periodic boundary conditions
  - Captures spatial oscillations due to atomic periodicity
- **Z-direction**: Exponential decay `exp(-αz)` with α=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0] Å⁻¹
  - Models decay of surface interactions into vacuum
  - Different decay rates capture near-field vs far-field behavior

**Total**: 8×4 = 32 basis functions

**Design Considerations**:
- **Modularity**: Easy to swap basis sets for experimentation
- **Pre-computation**: Store basis functions as 2D grids
- **Linear Operations**: All fitting done via linear combinations

**Interface**:
```python
class BasisFunctions:
    def generate_plane_wave_exponential_basis(n_harmonics, n_z_functions)
    def get_basis_matrix() -> (n_points, n_basis)
    def reconstruct_potential(coefficients) -> potential_2d
```

#### 3.2.4 PotentialFitter
**Purpose**: Fit grid potentials to basis function expansion

**Method**: Regularized Least Squares
```
minimize: ||A*c - b||² + λ||c||²
```
where:
- `A`: basis matrix (n_grid_points × n_basis_functions)
- `c`: coefficients to determine
- `b`: target potential values on grid
- `λ`: regularization strength

**Regularization Benefits**:
- Prevents overfitting to grid noise
- Enables suppression of unnecessary basis functions
- Improves numerical stability

**Interface**:
```python
class PotentialFitter:
    def set_regularization(strength)
    def fit(target_potential) -> (coefficients, rmse)
    def get_fitted_potential() -> potential_2d
```

## 4. Implementation Strategy

### 4.1 Phase 1: Core Framework (Current)
- [x] Basic module structure
- [x] 2D grid management
- [x] Coulomb potential with PBC
- [x] Plane wave × exponential basis
- [x] Regularized least squares fitting
- [x] Visualization tools

### 4.2 Phase 2: Optimization and Extension
- [ ] Morse potential implementation
- [ ] Basis function optimization (genetic algorithms, bayesian optimization)
- [ ] Adaptive grid refinement
- [ ] Performance benchmarking
- [ ] Memory optimization

### 4.3 Phase 3: Advanced Features
- [ ] Multiple atom types and potentials
- [ ] 3D extension for bulk crystals
- [ ] Integration with molecular dynamics codes
- [ ] Machine learning basis function discovery

## 5. Key Design Decisions

### 5.1 2D vs 3D Grid
**Decision**: Start with 2D (x,z) grid
**Rationale**:
- Surface problems are naturally 2D
- Faster computation and visualization
- Easier debugging and validation
- Can extend to 3D later

### 5.2 Basis Function Choice
**Decision**: Plane waves (x) × Exponentials (z)
**Rationale**:
- Plane waves naturally satisfy PBC in x
- Exponentials model surface decay in z
- Separable functions simplify analysis
- Well-understood mathematical properties

### 5.3 Regularization Strategy
**Decision**: L2 regularization (Ridge regression)
**Rationale**:
- Simple to implement and understand
- Provides smooth coefficient shrinkage
- Allows tuning for basis function suppression
- Numerically stable

### 5.4 Grid Storage
**Decision**: Store basis functions as pre-computed grids
**Rationale**:
- Avoids repeated function evaluations
- Enables fast linear algebra operations
- Memory cost is acceptable for 2D grids
- Simplifies basis function swapping

## 6. Performance Considerations

### 6.1 Computational Complexity
- **Grid Setup**: O(N_x × N_z) - done once
- **Potential Calculation**: O(N_atoms × N_images × N_x × N_z) - done once
- **Basis Generation**: O(N_basis × N_x × N_z) - done once
- **Fitting**: O(N_basis³ + N_basis² × N_points) - done once
- **Evaluation**: O(N_basis) - done many times

### 6.2 Memory Requirements
- Grid storage: ~10MB for 50×50×32 basis functions
- Potential storage: ~10KB for 50×50 grid
- Coefficient storage: ~1KB for 32 coefficients
- Total: Manageable for current problem size

### 6.3 Accuracy Control
- **Grid resolution**: Denser grid → higher accuracy, more memory
- **Basis set size**: More functions → better fit, higher cost
- **Periodic images**: More images → better convergence, slower pre-computation
- **Regularization**: Stronger → smoother fit, potential underfitting

## 7. Validation Strategy

### 7.1 Unit Tests
- Grid coordinate calculations
- Periodic image summation convergence
- Basis function orthogonality and completeness
- Fitting algorithm correctness

### 7.2 Physical Validation
- Compare with analytical solutions (simple cases)
- Energy conservation in molecular dynamics
- Convergence with respect to grid resolution
- Convergence with respect to basis set size

### 7.3 Performance Validation
- Timing comparisons with direct summation
- Memory usage profiling
- Scaling behavior with system size

