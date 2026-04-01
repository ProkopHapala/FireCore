# Particle on Spring: Hysteresis and Manipulation Trajectory Demonstrations

This directory contains educational Python scripts that demonstrate the physics of particle manipulation using spring-based force microscopy, particularly focusing on hysteresis phenomena in double-well potentials. These simulations are designed to help students understand the fundamental principles of atomic force microscopy (AFM) and molecular manipulation.

## Overview

The scripts demonstrate how a particle (representing an atom or molecule) behaves when manipulated by a spring-connected tip, similar to how AFM manipulates individual atoms on surfaces. The key phenomenon illustrated is **hysteresis** - where the system's response depends on its history, leading to different paths for forward and backward manipulation.

## Physical Background

### Double-Well Potential
The surface is modeled using a double-well potential, which has two stable equilibrium positions separated by an energy barrier. This is a common model for:
- Atomic positions on crystalline surfaces
- Molecular conformations with two stable states  
- Bistable systems in nanomechanics

Mathematical form:
```
V_surf(x) = A * cos(k * x)
```
where:
- `A` is the amplitude (depth of the wells)
- `k` is the wave vector (determines the spacing between wells)
- The minima occur at x = ±π/k

### Spring-Particle System
The particle is connected to a movable tip via a harmonic spring:
```
V_spring = 0.5 * K_tip * (x - x_tip)²
```
where:
- `K_tip` is the spring stiffness
- `x_tip` is the controlled tip position
- `x` is the particle position

### Total Energy
```
V_total = V_surf(x) + V_spring
```

The system finds equilibrium by minimizing `V_total` for each tip position.

### Hysteresis Origin
Hysteresis occurs when the spring stiffness is comparable to but less than the curvature of the surface potential. As the tip moves:
1. **Forward path**: Particle stays in the left well until the tip pulls it over the barrier
2. **Backward path**: Particle stays in the right well until the tip pushes it back over the barrier
3. **Different paths**: The transition points differ, creating a hysteresis loop

## Scripts Description

### 1. `histeresis.py` - 1D Hysteresis Demonstration

**Purpose**: Demonstrates classic hysteresis in a 1D double-well system with clear visualization of the energy landscape and force relationships.

**Key Features**:
- **Three-panel visualization**:
  1. Energy landscape with tip-molecule connectors
  2. Force vs tip position curves  
  3. Hysteresis loop (tip position vs molecule position)

**Physics Illustrated**:
- Energy minimization at each tip position
- Force transmission through the spring
- History-dependent behavior (hysteresis)
- Bifurcation and sudden transitions

**Usage**:
```bash
python histeresis.py [options]
```

**Key Parameters**:
- `--amplitude A`: Surface potential amplitude (default: 1.0)
- `--stiffness K`: Tip spring stiffness (default: 0.5)
- `--tip-max x`: Maximum tip displacement (default: 4.5)
- `--points n`: Number of sample points (default: 300)

**Educational Value**:
Students can observe how:
- Changing stiffness affects hysteresis loop size
- The system "remembers" which well it came from
- Forces change abruptly at transition points
- Energy barriers determine manipulation difficulty

### 2. `serpent.py` - 2D Manipulation Paths

**Purpose**: Extends the concept to 2D, showing how different manipulation paths (straight vs. serpentine) affect the required forces and can be used to overcome energy barriers.

**Key Features**:
- **2D surface potential**: Double-well in x + harmonic confinement in y + Gaussian hill
- **Two path types**: Straight line vs. serpentine (curved) trajectory
- **Force limits**: Demonstrates bond breaking when force exceeds threshold
- **Real-time force monitoring**: Shows when manipulation fails

**Physics Illustrated**:
- 2D energy landscapes
- Path-dependent manipulation strategies
- Force constraints and bond breaking
- Strategic path planning to avoid high forces

**Usage**:
```bash
python serpent.py [options]
```

**Key Parameters**:
- `--stiffness K`: Tip spring stiffness (default: 3.3)
- `--break-force F`: Maximum sustainable force (default: 4.0)
- `--serp-amp A`: Amplitude of serpentine path (default: 1.5)
- `--hill-amp H`: Central hill height (default: 10.0)

**Educational Value**:
Students learn:
- How 2D paths can reduce required forces
- The importance of path planning in manipulation
- When and why bonds break during manipulation
- Trade-offs between direct vs. indirect paths

### 3. `serpent_bakF.py` - Backup Version

**Purpose**: Backup version of the serpentine path simulation with potentially different parameterizations or experimental features.

## Mathematical Details

### Equilibrium Conditions
The system finds equilibrium where forces balance:
```
F_spring + F_surface = 0
K_tip * (x - x_tip) + dV_surf/dx = 0
```

### Stability Analysis
At equilibrium, stability requires:
```
d²V_total/dx² > 0
K_tip + d²V_surf/dx² > 0
```

### Critical Spring Stiffness
For hysteresis to occur:
```
K_tip < A * k²
```
If `K_tip` is too large, the particle follows the tip perfectly (no hysteresis).
If `K_tip` is too small, the particle doesn't move with the tip (no manipulation).

## Installation and Dependencies

**Required packages**:
```bash
pip install numpy matplotlib scipy
```

**Python version**: 3.7+ recommended

## Teaching Applications

### Physics Courses:
- **Classical Mechanics**: Energy minimization, stability, phase transitions
- **Statistical Mechanics**: Double-well potentials, metastable states
- **Computational Physics**: Numerical optimization, force calculations

### Materials Science/Nanotechnology:
- **AFM Principles**: Force microscopy, atomic manipulation
- **Surface Science**: Adsorption sites, energy barriers
- **Nanomechanics**: Bistable systems, switching mechanisms

### Problem Sets for Students:

1. **Parameter Exploration**: How does changing `K_tip` affect the hysteresis loop width?
2. **Critical Points**: Find the exact spring stiffness where hysteresis disappears
3. **Energy Calculations**: Calculate the work done during forward vs. backward paths
4. **Path Optimization**: In `serpent.py`, find the optimal serpentine amplitude to minimize maximum force
5. **Force Analysis**: Derive the analytical expressions for the forces and verify numerically

## Advanced Topics

### Extensions for Further Study:
- **Thermal Effects**: Add temperature-dependent transitions
- **Damping**: Include viscous damping for dynamic simulations
- **Multiple Particles**: Extend to interacting particle systems
- **Real Potentials**: Use actual DFT-calculated surface potentials

### Connection to Real Experiments:
- **AFM Manipulation**: Direct analogy to real atomic manipulation experiments
- **Molecular Switches**: Models for bistable molecular devices
- **Nanomechanical Memory**: Principles for mechanical information storage

## Running the Scripts

### Quick Start:
```bash
# Basic hysteresis demonstration
python histeresis.py

# Custom parameters
python histeresis.py --stiffness 0.3 --amplitude 2.0 --output my_hysteresis

# 2D manipulation paths
python serpent.py --stiffness 2.0 --serp-amp 2.0

# Force-limited manipulation
python serpent.py --break-force 3.0 --hill-amp 15.0
```

### Output:
- Scripts display interactive plots by default
- Use `--output filename` to save figures
- Supported formats: PNG, SVG
- High-resolution output (300 DPI) for publications

## Troubleshooting

**Common Issues**:
1. **No hysteresis observed**: Reduce spring stiffness `K_tip`
2. **No transitions**: Increase tip displacement range `--tip-max`
3. **Optimization failures**: Check initial conditions and potential parameters
4. **Force too high**: Reduce spring stiffness or increase surface potential amplitude

**Parameter Guidelines**:
- For clear hysteresis: `K_tip < A * k²`
- For smooth transitions: Use sufficient points (`--points > 200`)
- For force-limited manipulation: Set `--break-force` appropriately

## References and Further Reading

1. **AFM Manipulation**: 
   - Meyer, G. et al. "Manipulation of single atoms on insulating surfaces"
   - Hla, S.W. "Atomic manipulation and its applications"

2. **Double-Well Physics**:
   - Landau, L.D. & Lifshitz, E.M. "Statistical Physics" (phase transitions)
   - Strogatz, S.H. "Nonlinear Dynamics and Chaos" (bifurcation theory)

3. **Computational Methods**:
   - Press, W.H. et al. "Numerical Recipes" (optimization algorithms)
   - SciPy documentation for `minimize` function
