CODEMAP: FDBM Fitting Infrastructure: GridFF Sampling, Ewald Electrostatics, and Folded Basis Fitting
https://windsurf.com/codemaps/a8b020b7-8af8-46db-bb73-e6a0627021f6-fe86ab10a43f3d18

---

# USER

I want to fit GridFF potential based on Full density based model.
Full density based model is basically sum of Electrostatic Pential, Pauli repulsion, and Van Der Waals for each atom i
so we have total energy being sum over atoms
E =Sum_i P_i*V_Pauli(r_i) + L_i*V_London(r_i) + Q_i*V_Coulomb(r_i)

In fact Pauli repulsion in first approximation is just electron density to some power, Lenard jones can be basically computed as C6/r^6, co L_i is basically C6 coefficient.
In second phase we should ad the C6/r^6 poential and fit C6 (L_i) for each atom if it imporves the fit

Besid e that I have DFT  ineraction energy of small molecules with the flat NaCl substrate for multiple positions x,y,z and rotation (rigid body) so I can hopefully fit these coefficiets from that

my molecules are these
H2CO-O  H2O-H  H2O-O  HCONH2-H1  HCONH2-H2  HCONH2-O  HCOOH-H  HCOOH-O1  HCOOH-O2  NH3-H  NH3-N


In the first phase we may perhaps even neglect vdW which is typically small and consider just electrostatics and Pauli, then basically we have repulsion proportional to charge density multiplied by P_i which is different for each atoms and corespond to size of that atom (if we assume that density decay exponentilly from factor P_i=exp(-beta*R_i)), and then we have just charge for each atom which we can take as RESP charge basically

we have from DFT there sinputs
charge_density.cube
electrostatic_potential.cube

NaCl.xyz
```
32
Lattice="     11.2826 0.0 0.0 0.0      11.2826 0.0 0.0 0.0      40.0000" Properties=species:S:1:pos:R:3 #comment=T pbc="T T F"
Na    0.00000    0.00000    10.00000
Cl    2.82065    0.00000    10.00000
Na    5.64130    0.00000    10.00000
Cl    8.46195    0.00000    10.00000
Cl    0.00000    2.82065    10.00000
Na    2.82065    2.82065    10.00000
Cl    5.64130    2.82065    10.00000
Na    8.46195    2.82065    10.00000
Na    0.00000    5.64130    10.00000
Cl    2.82065    5.64130    10.00000
Na    5.64130    5.64130    10.00000
Cl    8.46195    5.64130    10.00000
Cl    0.00000    8.46195    10.00000
Na    2.82065    8.46195    10.00000
Cl    5.64130    8.46195    10.00000
Na    8.46195    8.46195    10.00000
Cl    0.00000    0.00000    12.82065
Na    2.82065    0.00000    12.82065
Cl    5.64130    0.00000    12.82065
Na    8.46195    0.00000    12.82065
Na    0.00000    2.82065    12.82065
Cl    2.82065    2.82065    12.82065
Na    5.64130    2.82065    12.82065
Cl    8.46195    2.82065    12.82065
Cl    0.00000    5.64130    12.82065
Na    2.82065    5.64130    12.82065
Cl    5.64130    5.64130    12.82065
Na    8.46195    5.64130    12.82065
Na    0.00000    8.46195    12.82065
Cl    2.82065    8.46195    12.82065
Na    5.64130    8.46195    12.82065
Cl    8.46195    8.46195    12.82065
```

based on this information can you design efficient way how to fit it.
I was thinking in the first phase we simply sample the potential using GridFF atom position of atom for all configurations from the training examples, store them separately for Pauli, Electrostatics and vdW, and then I will randomly vary coefficietsn of P_i, L_i, Q_i (Q_i is actually fixed by RESP) for the atom types , note that multiple atoms can have same type, also for different molecule, type assigment is quite importaint here,

---

# Gemini-3.5-flash

This is a very elegant approach to a Grid Force Field (GridFF) parametrization. By leveraging the spatial grid grids of the substrate ($V_{\text{Coulomb}}$, $\rho_{\text{Pauli}}$) obtained from DFT, you can decouple the expensive 3D spatial evaluation from the coefficient optimization loop.

Because your total interaction energy is linear with respect to the scaling parameters ($P_i$, $L_i$), you don't even need to randomly vary them in a slow Monte Carlo loop. Once the grid potentials are sampled at the atomic coordinates for all training configurations, you can solve for the optimal coefficients using **Linear Least Squares (with regularized/constrained bounds)**.

Here is a comprehensive strategy to structure and execute this fitting procedure efficiently.

## Phase 1: Atom Type Classification

Before collecting the grid samples, you must establish a consistent mapping from your molecular configurations to unique atom types. Since you have specific active sites/orientations labeled (like H2CO-O, H2O-H, HCOOH-O1), your types should reflect chemical environment, not just element.

A robust type map for your dataset would look like this:

| Atom Type Index | Representative Label | Elements Included | Notes / Constraints |
| --- | --- | --- | --- |
| **1** | H_water | H from H2O-H | Fixed RESP charge |
| **2** | O_water | O from H2O-O | Fixed RESP charge |
| **3** | O_carbonyl | O from H2CO-O, HCONH2-O | Can group or keep separate |
| **4** | H_amine | H1, H2 from HCONH2 |  |
| **5** | O_acid_hydroxyl | O from HCOOH-O1 |  |
| **6** | O_acid_carbonyl | O from HCOOH-O2 |  |
| **...** | ... | ... | ... |

## Phase 2: Grid Sampling & Feature Matrix Generation

For every molecular configuration $k$ (which specifies the rigid translation $x,y,z$ and rotation $R$ of the molecule over the NaCl surface), you calculate the absolute position $\vec{r}_{i,k}$ of each atom $i$ in the system.

Instead of re-evaluating grids during optimization, you pass through the training data **exactly once** to build a feature design matrix.

### 1. Prerequisite Cubes from Substrate DFT

From your isolated NaCl slab DFT calculations, you need:

* **$V_{\text{electrostatic}}(\vec{r})$**: The raw electrostatic potential grid.
* **$\rho_{\text{electron}}(\vec{r})$**: The electron density grid. You can approximate the Pauli potential locally as a power of the density:

$$V_{\text{Pauli}}(\vec{r}) = [\rho_{\text{electron}}(\vec{r})]^\alpha$$


*(where $\alpha \approx 1.0$ to $2.0$; typically $\alpha = 1.0$ for a pure density overlay, or $\alpha = 4/3$ or $2$ for kinetic energy density approximations).*
* **$V_{\text{London}}(\vec{r})$**: A pre-computed geometric $1/r^6$ grid generated by placing a dummy $C_6=1$ probe over the NaCl lattice positions under periodic boundary conditions (PBC).

### 2. Math Formulation for the Linear System

For a single configuration $k$, the model energy is:

$$E_k = \sum_{t \in \text{types}} \left( P_t \sum_{i \in t} V_{\text{Pauli}}(\vec{r}_{i,k}) + L_t \sum_{i \in t} V_{\text{London}}(\vec{r}_{i,k}) \right) + E_{\text{Coulomb}, k}$$

Where the Coulomb contribution is completely fixed beforehand since both $Q_i$ (RESP) and $V_{\text{Coulomb}}$ are known:

$$E_{\text{Coulomb}, k} = \sum_{i} Q_i V_{\text{Coulomb}}(\vec{r}_{i,k})$$

We move the known term to the target side to get our residual $\Delta E_k$:


$$\Delta E_k = E_{\text{DFT}, k} - E_{\text{Coulomb}, k}$$

Now we can express this as a matrix equation $\mathbf{A}\mathbf{c} = \mathbf{b}$:

* **$\mathbf{b}$**: A vector of length $M$ (total number of configurations), where $b_k = \Delta E_k$.
* **$\mathbf{c}$**: The vector of unknown coefficients we want to fit: $\mathbf{c} = [P_1, P_2, \dots, L_1, L_2, \dots]^T$.
* **$\mathbf{A}$**: A matrix of shape $M \times (2 \times N_{\text{types}})$.
* The column corresponding to $P_t$ for row $k$ contains $\sum_{i \in t} V_{\text{Pauli}}(\vec{r}_{i,k})$.
* The column corresponding to $L_t$ for row $k$ contains $\sum_{i \in t} V_{\text{London}}(\vec{r}_{i,k})$.

## Phase 3: Step-by-Step Execution Pipeline

To implement this efficiently, follow this structural pipeline:

1. **Trilinear Interpolation Engine:** Highly Optimized.
Load charge_density.cube and electrostatic_potential.cube. Write a fast trilinear interpolation function in Python (using scipy.ndimage.map_coordinates or custom Numba code) that handles the periodic boundary conditions (PBC) along X and Y based on your lattice vector 11.2826.


2. **Extract Raw Coordinates and DFT Targets:** Data Pre-processing.
Parse your rigid-body configuration files. For each configuration $k$, calculate the transformed coordinates $\vec{r}_{i,k}$ for every atom in the molecule. Match these with their corresponding $E_{\text{DFT}, k}$.


3. **Populate the A Matrix and b Vector:** Grid Sampling Loop.
Loop through all configurations. For each atom, sample $V_{\text{electrostatic}}$, $V_{\text{Pauli}}$, and $V_{\text{London}}$ at $\vec{r}_{i,k}$. Accumulate these into their respective atom-type slots in row $k$ of matrix $\mathbf{A}$. Compute $b_k = E_{\text{DFT},k} - E_{\text{Coulomb},k}$.


4. **Bounded Solver Optimization:** Matrix Inversion.
Instead of using standard unconstrained pseudo-inverses, pass $\mathbf{A}$ and $\mathbf{b}$ into a bounded least-squares solver (scipy.optimize.lsq_linear). Set constraints: $P_t > 0$ (repulsion must be positive) and $L_t > 0$ (if $V_{\text{London}}$ is defined as a negative attraction grid, or vice versa).


## Phase 4: Refinement Strategies

Once your infrastructure is built, you can easily tweak the model parameters without re-running the heavy grid sampling loop:

### 1. Multi-Stage Fitting (As you planned)

* **Stage 1:** Zero out the $L_t$ columns entirely. Fit only $P_t$ using your Electrostatic + Pauli model. Check the root-mean-square error (RMSE).
* **Stage 2:** Include the $L_t$ columns. Check if the training residual decreases significantly and ensure that the fitted $L_t$ values correlate well with empirical atomic $C_6$ trends (e.g., $C_6(\text{O}) < C_6(\text{N})$).

### 2. Tuning the Pauli Exponent ($\alpha$)

Because $V_{\text{Pauli}} = \rho^\alpha$ is non-linear with respect to $\alpha$, you can optimize $\alpha$ using a simple outer 1D scalar optimization loop:

```python
import numpy as np
from scipy.optimize import minimize_scalar, lsq_linear

def objective_function(alpha, rho_sampled, matrix_other_columns, b_vector):
    # Recompute only the Pauli columns of matrix A using the new alpha power
    A_Pauli = np.sum(rho_sampled ** alpha, axis=axis_indices) 
    A = np.hstack([A_Pauli, matrix_other_columns])
    
    # Solve linear least squares with lower bounds at 0
    res = lsq_linear(A, b_vector, bounds=(0, np.inf))
    return res.cost # Returns the sum of squared residuals
    
# Find the global optimal alpha value 
optimal_result = minimize_scalar(objective_function, bounds=(1.0, 2.5), method='bounded')

```

### 3. Handling Outliers (Weighting Scheme)

Configurations where the molecule is too close to the surface will have massive, exponential DFT energies that can completely dominate the least-squares loss function. To prevent the fit from ignoring the critical equilibrium well regions, apply a configuration weight factor $W_k$:

* Set $W_k = 1.0$ for configurations near the equilibrium distance and long-range attraction regimes.
* Set $W_k \rightarrow 0.0$ (or a small threshold) for configurations where $E_{\text{DFT}}$ exceeds a high energy cutoff (e.g., $> 1.0\text{ eV}$ above the vacuum energy).
* Implement this by multiplying both row $A_{k,\cdot}$ and $b_k$ by $\sqrt{W_k}$.

This strategy guarantees a fast setup, physically meaningful coefficients, and zero danger of getting stuck in local minima during parameter optimization.

---

# USER

OK, I currently do not have Data on this computer, so I want to generate some mockup

We generate both reference configurations of molecules over the surface and their energies by sampling another GridFF, and then we check if we can reproduce it.

Similarly for electron desnity - we simply sum up exponential functions of the morese potential (just the repulsive part) and call it electron density (it will be like rho_NA sum of atoms, no bond charges)

these modules to geenrate mockup reference energies and mockup density should be separated as we can easily replace them by true data when we get them.

For generating molecular configurations, we should take two apporahcs

1) random confs ( random position and rotation)
2) systematic scan in z and in xor y with given rotation. Plese make sure we can generate homogenous rotation samples.

The second systematic scans are importaint for diagnostics as they can be easily plotted (human understandable trends).

look how I sample rotation here:
/home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/sample_driver.py
/home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/sampling_lib.py

now analyse the problem, and write notes how to do it, key code sample, but do not program the whole things yet, we need to know where we have functions which we can reuse

---

# Gemmini-3.5-flash

## Analysis and Architectural Plan

To test the fitting framework without real DFT data, we will build a self-contained **Mock Generation & Validation Pipeline**. By using an analytic potential to generate both the target energies and the mock substrate grids, we can guarantee that our linear fitting approach works: if the pipeline is correct, the fitting code must recover the exact parameters used to generate the mock data (a self-consistency check).

Based on the provided codebase (`sample_driver.py` and `sampling_lib.py`), we already have excellent infrastructure for handling geometry, curvilinear grids ($u, v$), $y$-level blending, and uniform angular hemisphere distributions (`generate_tilts`). We should leverage this rather than rewriting it.

## 1. Reusing the Existing Architecture

We can reuse several components from your current library:

* **Systematic Sampling Loops:** The `do_scan` and `do_random` paths in `sample_driver.py` already generate 5D coordinate frames ($\vec{r}_{i,k}$) for molecular geometries over the substrate. We can intercept these frames right before or after they write to `.xyz` to evaluate our mock energy functions.
* **Angular Distribution:** `generate_tilts` uses a Fibonacci-like fan disk mapping (`map_to_disk` and `lift_to_dome`) to ensure homogeneous orientation sampling over the hemisphere. This fully satisfies the request for uniform rotation distributions.
* **Grid Interpolation:** The curvilinear $u, v$ mapping handles the tricky job of tracking the undulating landscape of the NaCl slab.

## 2. Mock Data Generator Components

We will separate the mock generator into three decoupled functions so they can be replaced by real `.cube` grids and DFT energy vectors later.

```
       +---------------------------------------------+
       |             Analytical Models               |
       |  (Morse Repulsion / Coulomb / C6 Dispersion)|
       +---------------------------------------------+
                              |
               +--------------+--------------+
               |                             |
               v                             v
   +-----------------------+     +-----------------------+
   |  Mock Grid Generator  |     | Mock Target Generator |
   | (Populates 3D Arrays) |     |  (Computes Exact E)   |
   +-----------------------+     +-----------------------+

```

### Component A: The Mock Grid Generator

Instead of parsing a `.cube` file from VASP or CP2K, we evaluate analytical fields over a standard 3D Cartesian grid box covering the NaCl substrate space:

* **Mock Electron Density ($\rho_{\text{mock}}$):** Sum of atomic repulsions from substrate atoms $j$:

$$\rho_{\text{mock}}(\vec{r}) = \sum_{j \in \text{NaCl}} A_j \exp(-\beta_j |\vec{r} - \vec{R}_j|)$$

* **Mock Electrostatic Potential ($V_{\text{electrostatic, mock}}$):** Ewald or direct sum of point charges:

$$V_{\text{Coulomb, mock}}(\vec{r}) = \sum_{j \in \text{NaCl}} \frac{Q_j}{|\vec{r} - \vec{R}_j|}$$

* **Mock London Grid ($V_{\text{London, mock}}$):** Standard attractive dispersion grid assuming a probe with $C_6 = 1$:

$$V_{\text{London, mock}}(\vec{r}) = \sum_{j \in \text{NaCl}} \frac{-1}{|\vec{r} - \vec{R}_j|^6}$$


### Component B: The Mock Target Energy Function

For any given molecular configuration $k$, this function computes the *exact* analytical interaction energy by calculating the pairwise distance between every molecular atom $i$ and every substrate atom $j$:

$$E_{\text{Target}, k} = \sum_{i \in \text{mol}} \sum_{j \in \text{NaCl}} \left[ A_i A_j e^{-\beta (r_{ij})} - \frac{C_{6,i} C_{6,j}}{r_{ij}^6} + \frac{q_i Q_j}{r_{ij}} \right]$$


## 3. Key Code Design Strategy

Here is how to implement these components cleanly without bloating your production pipeline.

### Core Mock Engine (`mock_engine.py`)

```python
import numpy as np

class MockSubstrateGridFF:
    def __init__(self, nacl_atoms, lattice_dims):
        """
        nacl_atoms: list of dicts with 'species', 'x', 'y', 'z'
        lattice_dims: (abc) for managing periodic boundary conditions
        """
        self.atoms = nacl_atoms
        self.lat = lattice_dims
        
        # Ground-truth parameters for the substrate atoms
        self.param_map = {
            'Na': {'Q': +1.0, 'A': 1.5, 'beta': 1.8, 'C6': 1.0},
            'Cl': {'Q': -1.0, 'A': 2.5, 'beta': 1.6, 'C6': 5.0}
        }

    def generate_cartesian_cube(self, nx, ny, nz, bounds):
        """Generates mock 3D arrays to simulate loading .cube files."""
        x = np.linspace(bounds[0][0], bounds[0][1], nx)
        y = np.linspace(bounds[1][0], bounds[1][1], ny)
        z = np.linspace(bounds[2][0], bounds[2][1], nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        rho = np.zeros_like(X)
        v_coul = np.zeros_like(X)
        v_vdw = np.zeros_like(X)
        
        # Vectorized or simple loop accumulation over substrate atoms with PBC
        for at in self.atoms:
            p = self.param_map[at['species']]
            dx = X - at['x']  # Apply minimum image convention along X/Y if needed
            dy = Y - at['y']
            dz = Z - at['z']
            r = np.sqrt(dx**2 + dy**2 + dz**2) + 1e-8
            
            rho += p['A'] * np.exp(-p['beta'] * r)
            v_coul += p['Q'] / r
            v_vdw += -1.0 / (r**6)
            
        return rho, v_coul, v_vdw

    def compute_exact_interaction(self, mol_coords, mol_types, mol_charges):
        """
        Computes the target ground-truth energy vector directly from coordinates.
        This provides the baseline 'DFT energy' for your fitting loop.
        """
        # Ground truth values for molecular types that we expect to recover
        true_mol_params = {
            'O_water': {'P': 1.2, 'L': 2.3},
            'H_water': {'P': 0.1, 'L': 0.05},
        }
        
        total_energy = 0.0
        for i, (coord, t, q) in enumerate(zip(mol_coords, mol_types, mol_charges)):
            p_mol = true_mol_params[t]
            
            for at in self.atoms:
                p_sub = self.param_map[at['species']]
                r = np.linalg.norm(coord - np.array([at['x'], at['y'], at['z']])) + 1e-8
                
                # Pairwise interactions
                rep = p_mol['P'] * p_sub['A'] * np.exp(-p_sub['beta'] * r)
                disp = - (p_mol['L'] * p_sub['C6']) / (r**6)
                coul = (q * p_sub['Q']) / r
                
                total_energy += (rep + disp + coul)
                
        return total_energy

```

### Linking to Your Configuration Scans

You can intercept `_write_xyz_for_indices` (random mode) and `_write_xyz_line` (systematic scans) within `sample_driver.py` to generate the matching dataset arrays:

```python
def generate_mock_training_set(indices_iter, sys, grid_uv, y_vals, tilts, phi_grid, mock_ff):
    """
    Passes over configurations, building the training vector y and feature matrix X.
    """
    E_target = []
    configurations_coords = []
    
    # Assign mock atom types and RESP charges to your molecule
    # For example, H2O: [O_water, H_water, H_water]
    mol_types = ['O_water', 'H_water', 'H_water'] 
    mol_charges = [-0.834, 0.417, 0.417] 

    for (iu, iv, iy, it, ip) in indices_iter:
        # 1. Generate geometry using your existing framework
        xz = grid_uv[iy, iu, iv]
        y = y_vals[iy]
        up = local_up_from_streamline(grid_uv[iy], iu, iv)
        tvec = tilts[it]
        phi = phi_grid[ip]
        pos = np.array([xz[0], y, xz[1]])
        
        coords = make_pose_coords_local(sys, pos, up, tvec, phi, center=True)
        configurations_coords.append(coords)
        
        # 2. Evaluate target energy function
        energy = mock_ff.compute_exact_interaction(coords, mol_types, mol_charges)
        E_target.append(energy)
        
    return np.array(E_target), configurations_coords

```

## 4. Diagnostic Scans Strategy

Systematic 1D scans are valuable for code validation. Because we can isolate axes like $u$ (sliding along a streamline curve) or $y$ (lifting height above the surface), plotting the mock analytical energy against the grid-sampled fitted energy yields clean curves:

1. **Z-height / U-axis curves:** Show clear Morse-like exponential behavior. If the trilinear grid interpolation is implemented correctly, the grid-sampled values will perfectly track the analytical smooth curve down to the grid resolution limit.
2. **Rotation sweeps ($\phi$):** Create smooth periodic wave profiles. If there is a coordinate alignment bug in `make_pose_coords_local`, it will immediately break the symmetry of these curves, making it obvious during visual diagnostics.

Once you verify that the linear least-squares matrix solver recovers the exact parameter values (`true_mol_params`) from this mock data with an error near zero, you can safely drop in your real DFT `.cube` files and begin production fitting.

---

# USER

We are interested in problem of Full density based model (FDBM) fitting using mockup data (density, and molecular condiguration energy pairs) which we later replace by DFT reference. We have notes about how to do it here
@GridFF_FDBM_Fitting.md plese review our code most importantly @Surface.cl @Surface_utils.py @SurfaceEwald.py @test_electrostatics_comparison.py @test_folded_fit_nacl1x1.py @test_gridff_alignment.py 

and focus on how to integrate the general ideas from @GridFF_FDBM_Fitting.md with this syste, which is also desxribed in @RigidSurfPotential_GridFF.md 


Note the fitting should be done efficiently, which means we do not sample the potential on the fly, considering that filt is linar we somaple te potentil for all atoms for all atoms in all configurations and store them, and then we do just the dot products. For starter we ignore fitting of the the non-linear param alpha E_pauli=P_i * rho^alpha exactly because we want to have simple linear sytem where matrix product can be precomputed nad cashed efficiently. In the end we fit just few coefficients for few atoms types.

We start with H2O, which means we have just two types H and O. @H2O.xyz 
notice that charges are in the last colum 

do not program yet, make detailed notes how to implement it with maximal reuse of what is already doen in our system... 

---

# SWE-1.6

I'll review the documentation and existing code to understand how to implement FDBM fitting efficiently. Let me start by reading all the relevant files.

# FDBM Fitting Implementation Notes

## Overview

Implement Full Density Based Model (FDBM) fitting for H2O molecule on NaCl substrate using mockup data. The goal is to fit atom-type coefficients (P_i, L_i) for Pauli repulsion and London dispersion, while Coulomb charges Q_i are fixed from RESP charges.

## Mathematical Model

**Total energy per configuration k:**
```
E_k = Sum_i [P_i * V_Pauli(r_i) + L_i * V_London(r_i) + Q_i * V_Coulomb(r_i)]
```

**Linear system formulation:**
- For each atom type t, we fit coefficients P_t (Pauli) and L_t (London)
- Coulomb term is precomputed using fixed Q_i (RESP charges)
- Matrix equation: A * c = b, where:
  - b_k = E_DFT,k - E_Coulomb,k (residual after subtracting electrostatics)
  - A[k, :] = [Sum_{i in type1} V_Pauli(r_i), Sum_{i in type1} V_London(r_i), 
                Sum_{i in type2} V_Pauli(r_i), Sum_{i in type2} V_London(r_i), ...]
  - c = [P_1, L_1, P_2, L_2, ...]^T

## Key Design Principles

1. **Efficiency**: Pre-sample all potentials once, then reuse via dot products
2. **Modularity**: Separate mock data generation from fitting (easy DFT replacement)
3. **Reuse**: Leverage existing GridFF, Ewald2D, and fitting infrastructure
4. **Simplicity**: Ignore non-linear alpha parameter initially (E_pauli = P_i * rho, not rho^alpha)

## Implementation Architecture

### 1. Mock Data Generation Module (`pyBall/OCL/FDBM_mock.py`)

**Purpose**: Generate mock substrate grids and reference energies using analytical potentials

**Key Functions to Implement:**

```python
class MockSubstrateGridFF:
    """Generate mock substrate grids (rho, V_coulomb, V_london)"""
    
    def __init__(self, substrate_xyz, lattice_dims):
        # Load substrate atoms from XYZ
        # Define mock parameters for substrate atoms (Na, Cl)
        # param_map = {'Na': {'Q': +1.0, 'A': 1.5, 'beta': 1.8, 'C6': 1.0},
        #              'Cl': {'Q': -1.0, 'A': 2.5, 'beta': 1.6, 'C6': 5.0}}
    
    def generate_density_grid(self, nx, ny, nz, bounds):
        # rho(r) = Sum_j A_j * exp(-beta_j * |r - R_j|)
        # Reuse: exponential decay pattern from getMorsePLQH in Surface.cl
    
    def generate_coulomb_grid(self, nx, ny, nz, bounds):
        # V_coulomb(r) = Sum_j Q_j / |r - R_j|
        # Reuse: SurfaceEwaldCL for accurate periodic electrostatics
        # OR: Use simple Coulomb sum for mock (faster, less accurate)
    
    def generate_london_grid(self, nx, ny, nz, bounds):
        # V_london(r) = Sum_j -C6_j / |r - R_j|^6
        # Reuse: London dispersion pattern from existing GridFF generation
    
    def compute_exact_energy(self, mol_coords, mol_types, mol_charges, mol_params):
        # E = Sum_i Sum_j [P_i*A_j*exp(-beta*r_ij) - L_i*C6_j/r_ij^6 + q_i*Q_j/r_ij]
        # This provides the "DFT reference" energy
```

**Reuse Opportunities:**
- **SurfaceEwaldCL** ([pyBall/OCL/SurfaceEwald.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/SurfaceEwald.py:0:0-0:0)): Use for accurate Coulomb grid generation
- **getMorsePLQH** (`Surface.cl:172`): Reuse exponential repulsion formula
- **GridFF generation pipeline** (`pyBall/tests/ocl_GridFF_new.py`): Adapt grid sampling logic

### 2. Configuration Generation Module (`pyBall/OCL/FDBM_config.py`)

**Purpose**: Generate molecular configurations (positions + rotations) over substrate

**Key Functions to Implement:**

```python
def generate_random_configs(n_configs, substrate_lvec, z_range, mol_template):
    """Random positions and rotations"""
    # Reuse: Rotation generation from test_assembly.py:generate_rotations()
    # Reuse: PBC tiling logic from test_assembly.py
    # z_range: (z_min, z_max) above substrate
    
def generate_systematic_configs(scan_axes, substrate_lvec, mol_template):
    """Systematic scans for diagnostics (z-scan, x-scan, rotation-scan)"""
    # Reuse: sampling_lib.py for homogeneous rotation sampling
    # scan_axes: e.g., {'z': np.linspace(2.0, 6.0, 50), 'phi': np.linspace(0, 2pi, 12)}
    
def apply_rigid_transform(mol_coords, pos, quat):
    """Apply translation and rotation to molecule"""
    # Reuse: Transform logic from MolecularDynamics.eval_rigid_getSurfMorse()
    # Reuse: pack_transforms() from test_folded_fit_nacl1x1.py:34
```

**Reuse Opportunities:**
- **generate_rotations()** (`test_assembly.py:44`): Full3D/inplane/tilt rotation grids
- **sampling_lib.py** (mentioned in docs): Homogeneous hemisphere sampling
- **pack_transforms()** (`test_folded_fit_nacl1x1.py:34`): Transform packing for GPU

### 3. Potential Sampling Module (`pyBall/OCL/FDBM_sampling.py`)

**Purpose**: Pre-sample substrate potentials at all atom positions for all configurations

**Key Functions to Implement:**

```python
class PotentialSampler:
    """Sample substrate grids at atom positions efficiently"""
    
    def __init__(self, grid_data, g0, dg):
        # Load precomputed grids (rho, V_coulomb, V_london)
        # grid_data: (nx, ny, nz, 3) array
        # Reuse: load_gridff_array() from Surface_utils.py:67
    
    def sample_configs(self, configs_coords, configs_types):
        """
        Sample potentials for all atoms in all configurations
        
        Returns:
            sampled_pauli: (n_configs, n_atom_types) - Sum V_pauli per type
            sampled_london: (n_configs, n_atom_types) - Sum V_london per type
            sampled_coulomb: (n_configs, n_atom_types) - Sum V_coulomb per type
        """
        # For each config:
        #   For each atom:
        #     Trilinear interpolate grid at atom position
        #     Accumulate to atom-type bucket
        # Reuse: Trilinear interpolation from Surface_utils.py:sample_gridff_opencl()
        # Reuse: OpenCL B-spline sampling from RigidBodyAFM.sample_gridff_single_atom()
    
    def sample_batch_opencl(self, configs_coords, configs_types):
        """GPU-accelerated batch sampling using OpenCL"""
        # Reuse: RigidBodyDynamics from Surface_utils.py:650
        # Reuse: sample_grid_at_atoms_opencl() from Surface_utils.py:650
```

**Reuse Opportunities:**
- **load_gridff_array()** (`Surface_utils.py:67`): Grid loading with validation
- **sample_gridff_opencl()** (`Surface_utils.py:468`): OpenCL B-spline sampling
- **sample_grid_at_atoms_opencl()** (`Surface_utils.py:650`): Batch sampling interface
- **RigidBodyDynamics** (`pyBall/OCL/RigidBodyDynamics.py`): GPU harness for batch evaluation

### 4. Linear Fitting Module (`pyBall/OCL/FDBM_fit.py`)

**Purpose**: Solve linear least squares to fit P_i, L_i coefficients

**Key Functions to Implement:**

```python
def fit_linear_coefficients(sampled_pauli, sampled_london, target_energies, 
                           atom_types, fit_mask=None, bounds=(0, np.inf)):
    """
    Solve A * c = b using bounded least squares
    
    Args:
        sampled_pauli: (n_configs, n_types) - Pre-sampled Pauli potentials
        sampled_london: (n_configs, n_types) - Pre-sampled London potentials
        target_energies: (n_configs,) - E_DFT - E_Coulomb (residuals)
        atom_types: List of atom type names
        fit_mask: (n_configs,) boolean mask for filtering configs
        bounds: Lower/upper bounds for coefficients (P_i > 0, L_i > 0)
    
    Returns:
        coefficients: dict {type: {'P': value, 'L': value}}
        fit_stats: dict with RMSE, R^2, etc.
    """
    # Build design matrix A: (n_configs, 2*n_types)
    #   Columns: [P_1, L_1, P_2, L_2, ...]
    # Reuse: scipy.optimize.lsq_linear for bounded least squares
    # Reuse: Fitting framework from test_folded_fit_nacl1x1.py:245
    
def validate_fit(coefficients, sampled_pauli, sampled_london, target_energies):
    """Compute validation statistics (RMSE, correlation, etc.)"""
    # Reuse: rmse() from test_folded_fit_nacl1x1.py:14
    # Reuse: Correlation analysis from Surface_utils.py:630
```

**Reuse Opportunities:**
- **fit_folded_surface_basis()** (`MolecularDynamics.py:1029`): Existing fitting framework
- **lsq_linear** (`scipy.optimize`): Bounded least squares solver
- **rmse()** (`test_folded_fit_nacl1x1.py:14`): RMSE calculation utility

### 5. Main Test Script (`tests/tMMFF/test_fdbm_fit_h2o.py`)

**Purpose**: End-to-end FDBM fitting pipeline for H2O

**Workflow:**

```python
def main():
    # 1. Load H2O molecule
    mol_xyz = "cpp/common_resources/xyz/H2O.xyz"
    # Charges: O=-0.4, H=+0.2 (from last column)
    # Atom types: ['O_water', 'H_water', 'H_water']
    
    # 2. Load NaCl substrate
    sub_xyz = "cpp/common_resources/xyz/NaCl_1x1_L3.xyz"
    # Reuse: load_substrate_xyz_with_lvec() from Surface_utils.py:106
    
    # 3. Generate mock substrate grids
    mock_ff = MockSubstrateGridFF(sub_xyz, substrate_lvec)
    rho_grid, coul_grid, london_grid = mock_ff.generate_all_grids(...)
    
    # 4. Generate configurations
    configs_random = generate_random_configs(n=1000, ...)
    configs_systematic = generate_systematic_configs(scan_axes={'z': ...}, ...)
    
    # 5. Sample potentials (precompute and cache)
    sampler = PotentialSampler(rho_grid, coul_grid, london_grid, g0, dg)
    sampled = sampler.sample_configs(all_configs, atom_types)
    
    # 6. Compute reference energies (mock "DFT")
    ref_energies = mock_ff.compute_exact_energy(configs_coords, ...)
    
    # 7. Subtract Coulomb contribution (fixed Q_i)
    coulomb_contrib = sampled['coulomb'] @ charges_per_type
    target_residuals = ref_energies - coulomb_contrib
    
    # 8. Fit P_i, L_i coefficients
    coefficients = fit_linear_coefficients(
        sampled['pauli'], sampled['london'], target_residuals, atom_types
    )
    
    # 9. Validate and visualize
    plot_1d_scans(configs_systematic, coefficients, ...)
    plot_parity(ref_energies, predicted_energies, ...)
    
    # 10. Save results
    save_coefficients(coefficients, output_path)
```

**Reuse Opportunities:**
- **load_substrate_xyz_with_lvec()** (`Surface_utils.py:106`): Substrate loading
- **plot_1d_lines()** (`test_folded_fit_nacl1x1.py:54`): 1D scan visualization
- **compare_sampling_methods()** (`Surface_utils.py:506`): Method comparison utilities

## Data Structures

### Atom Type Mapping (H2O example)
```python
ATOM_TYPES = {
    'O_water': {'element': 'O', 'charge': -0.4, 'indices': [0]},
    'H_water': {'element': 'H', 'charge': +0.2, 'indices': [1, 2]}
}
```

### Configuration Storage
```python
configs = {
    'coords': (n_configs, n_atoms, 3) array,  # Cartesian positions
    'types': (n_configs, n_atoms) int array,   # Atom type indices
    'charges': (n_configs, n_atoms) array,    # RESP charges
    'quaternions': (n_configs, 4) array,      # Orientation quaternions
    'label': str  # 'random' or 'systematic_z_scan' etc.
}
```

### Sampled Potentials (Cached)
```python
sampled = {
    'pauli': (n_configs, n_types) array,   # Sum V_pauli per type per config
    'london': (n_configs, n_types) array,  # Sum V_london per type per config
    'coulomb': (n_configs, n_types) array  # Sum V_coulomb per type per config
}
```

## Implementation Phases

### Phase 1: Mock Data Generation (Standalone)
- Implement `MockSubstrateGridFF` class
- Test grid generation and visualization
- Verify analytical energy computation

### Phase 2: Configuration Generation
- Implement random and systematic config generators
- Reuse rotation sampling from existing code
- Test with H2O template

### Phase 3: Potential Sampling
- Implement trilinear interpolation (Python fallback)
- Implement OpenCL batch sampling (GPU acceleration)
- Benchmark performance

### Phase 4: Linear Fitting
- Implement bounded least squares fitting
- Add fit mask for energy filtering
- Validate with mock data (should recover exact parameters)

### Phase 5: Integration and Testing
- Assemble end-to-end pipeline
- Generate diagnostic plots (1D scans, parity plots)
- Test systematic vs random configurations

## Key Reusable Components Summary

| Component | Location | Reuse For |
|-----------|----------|-----------|
| GridFF loading | `Surface_utils.py:67` | Load mock grids |
| OpenCL sampling | `Surface_utils.py:468` | GPU potential sampling |
| Ewald2D electrostatics | [SurfaceEwald.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/SurfaceEwald.py:0:0-0:0) | Accurate Coulomb grids |
| Rotation generation | `test_assembly.py:44` | Config orientations |
| Transform packing | `test_folded_fit_nacl1x1.py:34` | GPU transform upload |
| Fitting framework | `MolecularDynamics.py:1029` | Linear least squares |
| RMSE utilities | `test_folded_fit_nacl1x1.py:14` | Validation metrics |
| 1D scan plotting | `test_folded_fit_nacl1x1.py:54` | Diagnostics |
| Substrate loading | `Surface_utils.py:106` | Parse NaCl.xyz |

## File Structure

```
pyBall/OCL/
├── FDBM_mock.py          # Mock substrate grid generation
├── FDBM_config.py        # Configuration generation
├── FDBM_sampling.py      # Potential sampling (CPU/GPU)
├── FDBM_fit.py           # Linear least squares fitting
└── FDBM_utils.py         # Shared utilities

tests/tMMFF/
└── test_fdbm_fit_h2o.py # End-to-end test script
```

## Next Steps

1. **Implement Phase 1**: Create `FDBM_mock.py` with analytical substrate grids
2. **Test mock generation**: Visualize grids, verify energy computation
3. **Implement Phase 2**: Add configuration generation with rotation reuse
4. **Implement Phase 3**: Add potential sampling with OpenCL reuse
5. **Implement Phase 4**: Add linear fitting with bounded LS
6. **Integrate**: Assemble pipeline and test with H2O

This approach maximally reuses existing GridFF, Ewald2D, and fitting infrastructure while keeping mock data generation modular for easy DFT replacement.

I've completed the detailed implementation notes for FDBM fitting. Here's a summary of the key points:

## Implementation Strategy

**Core approach**: Pre-sample substrate potentials (Pauli, London, Coulomb) at all atom positions for all configurations, then solve linear least squares to fit atom-type coefficients P_i and L_i.

**Architecture**: 5 modular Python files in `pyBall/OCL/`:
- `FDBM_mock.py` - Mock substrate grid generation (replaced by DFT later)
- `FDBM_config.py` - Configuration generation (random + systematic)
- `FDBM_sampling.py` - Potential sampling (CPU/GPU)
- `FDBM_fit.py` - Linear least squares fitting
- `FDBM_utils.py` - Shared utilities

**Key reuse opportunities**:
- **SurfaceEwaldCL** - Accurate Coulomb grid generation
- **sample_gridff_opencl()** - GPU B-spline sampling
- **RigidBodyDynamics** - Batch evaluation harness
- **fit_folded_surface_basis()** - Existing fitting framework
- **generate_rotations()** - Rotation sampling from assembly code
- **load_gridff_array()** - Grid loading utilities

**H2O-specific setup**:
- 2 atom types: O_water (charge -0.4), H_water (charge +0.2)
- Fit 4 coefficients: P_O, L_O, P_H, L_H
- Coulomb term precomputed using fixed RESP charges
- Linear system: A * [P_O, L_O, P_H, L_H]^T = b

**Implementation phases**:
1. Mock data generation (analytical potentials)
2. Configuration generation (random + systematic scans)
3. Potential sampling (trilinear + OpenCL)
4. Linear fitting (bounded least squares)
5. Integration and testing

The notes provide detailed function signatures, data structures, and specific line references for reusable code. Ready to proceed with implementation when you approve.

---

# What is already implemented (high-value reuse)

## 1) FDBM grid pipeline already exists (PyOpenCL, step-by-step, cached, with plots)
- **File:** [tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py:0:0-0:0)
- Implements stages:
  1. **Density projection to grid** ([step1_density_projection](cci:1://file:///home/prokophapala/git/FireCore/tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py:185:0-321:100)) using `pyBall.FireballOCL.GridProjector`
  2. **Electrostatics** via FFT Poisson ([step2_electrostatics](cci:1://file:///home/prokophapala/git/FireCore/tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py:408:0-491:15))
  3. **Pauli** via FFT convolution ([step3_pauli](cci:1://file:///home/prokophapala/git/FireCore/tests/tAFM/pyocl_fdbm/run_pyocl_fdbm.py:497:0-585:39))
  4. ES convolution, dispersion, composition, etc.
- Writes `.npy` grids + *lots of diagnostics plots* into `tests/tAFM/pyocl_fdbm/debug/step*/...`

So for “mock data / pipeline validation”, this already covers the “generate rho / V_es / pauli field / vdw field” part.

## 2) FDBM Pauli fitting scripts already exist (with good 1D scan visualizations)
You literally already have the “nice 1D scan reference vs fit” pattern you asked me to implement:
- **File:** `tests/tAFM/pyocl_fdbm/fit_fdbm_pauli.py`
- Fits model `E_ref(z) = A * overlap(z)^beta` using extracted z-profiles
- Generates plots (`fit_pauli.png`, summary plots, etc.)

Also a more integrated/utility-heavy version exists in:
- **File:** `pyBall/OCL/AFM_utils.py`
  - `_load_fdbm_grids`, `_fit_pauli_powerlaw`, `_plot_pauli_fit`, `fit_pauli_parameters(...)`

This means: plotting + fitting + extraction utilities are already there; we should reuse them, not rewrite.

## 3) GPU GridFF sampling kernel for many atoms/poses already exists
In [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:0:0-0:0):
- **[initGridFF(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1548:4-1619:53)** uploads the B-spline GridFF data to GPU
- **[run_sampleGridFF_Bspline_points(nSystems, PLQH=...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:953:4-967:27)** runs kernel `sampleGridFF_Bspline_points`
- You already have a mechanism to select components by **`PLQH`** (see [_component_plqh](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:984:4-990:41), [_resolve_plqh](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:992:4-996:46), `COMPONENT_PLQH`)

And in [pyBall/OCL/Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0) you already have:
- [sample_gridff_opencl(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:467:0-502:51)
- [sample_grid_at_atoms_opencl(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:649:0-723:67)
- [compare_sampling_methods(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:505:0-646:18)

So sampling “potential at all atoms for all configurations” can be done via existing OpenCL kernels — no Python hot loops needed.

## 4) Existing least-squares fitting pattern exists (folded basis)
- **File:** [pyBall/OCL/MolecularDynamics.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:0:0-0:0)
- **Method:** [fit_folded_surface_basis(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1156:4-1291:33)
- It builds a basis matrix and uses `np.linalg.lstsq` over weighted samples.

Even if folded basis isn’t the same model, the *fitting style + data flow* is already in-house.

# What is *missing* for your requested FDBM fitting (and should be minimal)

Your requested fitting is different from the existing AFM Pauli fitting in one key way:

- The AFM scripts fit *Pauli parameters* (A,beta) **against a 1D z-scan** extracted from a *field*.
- You want **linear fitting of atom-type coefficients** `P_t`, `L_t` from **many rigid-body configurations**:
  - build features `Φ_k,t = Σ_{i∈type t} ρ(r_{k,i})` (or Pauli field / overlap field)
  - `Ψ_k,t = Σ_{i∈type t} V_London(r_{k,i})`
  - plus Coulomb term using fixed charges.

What’s missing is basically a **small “feature accumulation” glue** that:
1. takes a batch of rigid transforms (configs)
2. samples the grid channels at all atom positions (GPU)
3. reduces per-atom sampled values into per-type sums (features)
4. solves linear LSQ and plots 1D scans

That “reduce to per-type features” can be CPU-side (it’s not hot; it’s `O(nconfigs*natoms)` with natoms small), but if you want, we can even keep it mostly on GPU by reusing repeated kernel calls with different `PLQH` masks. The key is: **no new sampling kernels are needed**.

# Concrete reuse-first plan (before coding)

## A) Decide which grid(s) are the “features”
We need clarity which fields you want to fit against:

- **Option A (closest to your FDBM equation):**
  - `V_coulomb(r)` from DFT electrostatic potential cube (or Ewald2D)
  - `rho(r)` from density cube (Pauli uses `rho`, initially α=1)
  - `V_london(r)` from some model grid (or from your existing dispersion stage)

- **Option B (reuse AFM pipeline outputs directly):**
  - use `tests/tAFM/pyocl_fdbm` outputs:
    - `overlap_raw.npy` (raw overlap ~ pauli feature)
    - `E_ES_field.npy` (electrostatics feature)
    - `E_vdw_field.npy` (dispersion feature)

Option B is *maximal reuse* for mock/testing, because those grids already exist and have diagnostics.

## B) Reuse existing OpenCL sampling (no new kernels)
Use [MolecularDynamics.initGridFF(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1548:4-1619:53) + [run_sampleGridFF_Bspline_points(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:953:4-967:27) to sample grid channels for all atoms in all configurations.

Key detail: your GridFF layout already supports multi-channel PLQ-like combination; by choosing `PLQH`, we can isolate:
- Pauli channel
- London channel
- Coulomb channel

(Exactly the pattern already used in [Surface_utils.compare_sampling_methods()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:505:0-646:18).)

## C) Minimal new code should likely be *one test script*, not new library modules
Something like:
- `tests/tMMFF/test_fdbm_fit_h2o_linear.py`

That script would:
1. load H2O ([cpp/common_resources/xyz/H2O.xyz](cci:7://file:///home/prokophapala/git/FireCore/cpp/common_resources/xyz/H2O.xyz:0:0-0:0))
2. load grids from an existing directory (either real `.cube` converted to `.npy`, or from `pyocl_fdbm` debug outputs)
3. generate rigid configurations (you already have config sampling utilities scattered in tests; we should reuse whichever you prefer)
4. sample features via existing OpenCL sampling
5. solve LSQ (numpy/scipy)
6. plot 1D z-scans / rotation scans using plotting code **reused from** `fit_fdbm_pauli.py` / `AFM_utils.py` style

# Two quick clarification questions (to avoid coding the wrong thing)

1. **What is your intended source of grids for the first implementation?**
   - **(a)** use the existing `tests/tAFM/pyocl_fdbm` generated `.npy` fields (best for quick end-to-end validation)
   - **(b)** use your real DFT `.cube` files (`charge_density.cube`, `electrostatic_potential.cube`) and just reuse sampling/fitting parts

2. **Do you want the first fit to include London (`L_t`) or start with only Pauli + Coulomb?**
   - Starting with Pauli-only linear fit is simplest and matches your “ignore alpha” first phase.

If you answer those, I’ll propose the *smallest set of files* to touch and exactly which existing functions/classes we’ll call — **before** writing code again.

# Status
- **Done:** Inventory found strong existing reuse candidates: `tests/tAFM/pyocl_fdbm/*`, `pyBall/OCL/AFM_utils.py`, `pyBall/OCL/MolecularDynamics.py:initGridFF/run_sampleGridFF_Bspline_points`, `pyBall/OCL/Surface_utils.py:sample_gridff_opencl/sample_grid_at_atoms_opencl`.
- **Next:** Confirm grid source + fit terms, then implement *minimal glue* (likely 1 test script + maybe a tiny helper function if absolutely necessary).

---

# USER

1) the pyocl_fdbm is FDBM for AFM of molecule if I remmber? We want to do it for NaCl surface. So we must make new mockup data. Later we will replace them with DFT data (cube files)
2) all the reusable function should be put in @Surface_utils.py (e.g. the one for generating mocups etc.)
3) make general pipeline icluding London, butmaybe intitially we will turn it of (fix it to zero) for simplicuyt (rerducing number of params), the pipeline must be flaxible allowing it
4) Plese be sure you do not somethig stupid or ineffciient like recomputing GridFF again and again!
Perhaps we can already reuse the GridFF as it is in @Bspline_PLQd.npy @Bspline_PLQd_meta.json 
We can consider that Pauli can be rised to power alpha to get different decays, and that it is more-or less already a density (they are all just sum of exponential decays from atom centers more or less..., so this is already our mockup.

5) So what we will fit you may ask? That is the point we want to fit the same thing by hthe same thing to see if the fit works. But we can make some parturbations (noise, rescaling, rsing to power close to 1.0) to see how the fit deals with these pertubations. Put there options to do that. 


---

# GPT-5.2-low

## What I did (reuse-first, no GridFF recomputation)

I implemented the **reusable mock-fitting helpers directly in** [pyBall/OCL/Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0), reusing:

- The existing **NaCl GridFF Bspline grid** `Bspline_PLQd.npy` + [Bspline_PLQd_meta.json](cci:7://file:///home/prokophapala/git/FireCore/cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd_meta.json:0:0-0:0)
- Existing OpenCL kernel **`sampleGridFF_Bspline_points`** (already compiled/used by [MolecularDynamics](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:59:0-2007:21))
- Existing rigid-batch transform upload ([MolecularDynamics.upload_rigid_transforms](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1811:4-1824:19))

So the pipeline will:
1) **Load GridFF once**
2) **Upload Bspline once**
3) For many configs, **only upload transforms + run sampling kernel** in chunks
4) Build feature matrix on CPU using vectorized matrix multiplications (no “stupid hot loops”)

## New reusable functions added to [Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0)

- **[load_bspline_gridff(grid_path)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:91:0-100:21)**
  - loads `.npy` grid + requires `_meta.json`
  - checks `meta['ns']` matches array shape

- **[init_gridff_sampler_md(grid_path, apos0, nSystems, use_texture=False)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:103:0-122:19)**
  - creates a [MolecularDynamics](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:59:0-2007:21) instance
  - initializes rigid molecule batch once
  - calls [initGridFF(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1548:4-1619:53) to upload Bspline data once

- **[sample_gridff_channels_rigid(md, transforms, PLQH_channels=...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:125:0-148:13)**
  - samples **Pauli/London/Coulomb** channels separately by passing different `PLQH`
  - returns `Es[nconf,natoms,nch]` where each entry is **per-atom energy** from the channel
  - uses chunking via `md.nSystems` (avoids realloc / avoids re-uploading grid)

- **[fdbm_build_feature_matrix(Es_PLQ, type_ids, ntypes, use_london=True)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:151:0-167:15)**
  - turns per-atom samples into a **linear feature matrix**:
    - `Σ_{atoms of type t} V_pauli`
    - optionally `Σ_{atoms of type t} V_london`
  - done via `@` with a one-hot mapping matrix (fast, vectorized)

- **[fdbm_make_mock_reference(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:170:0-208:23)**
  - generates mock “DFT reference” energies by **composing the same sampled channels**
  - adds perturbations you requested:
    - **Pauli power** `pauli_alpha` (close to 1.0)
    - per-channel **rescale**
    - additive **Gaussian noise**
    - **London on/off** (if off, contribution is 0)

These are deliberately *thin glue* around existing GridFF/OpenCL infrastructure.

## Important detail (energy sign)
The OpenCL kernel writes `forces.w = -E` (see [relax_multi.cl](cci:7://file:///home/prokophapala/git/FireCore/cpp/common_resources/cl/relax_multi.cl:0:0-0:0) comment and assignment), so I convert back with:
- `E = -aforce[...,3]`

## Next step (pending)
You asked for a **general pipeline including London (but optionally disabled)** + a minimal test with 1D scans and fitting “self to self with perturbations”.

I have not created the test script yet. Next I will add a **single test** under `tests/tMMFF/` that:

- Loads `cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy`
- Builds a set of rigid configurations (z-scan + random)
- Samples `V_P/V_L/V_Q` using the new [Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0) helpers (GPU)
- Creates mock reference `E_ref` using [fdbm_make_mock_reference(...)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:170:0-208:23)
- Fits `(P_type, L_type)` by linear least squares (London optionally disabled)
- Plots **1D scan: reference vs fitted** and residuals

## One clarification before I write the test
For the very first working demo, do you want to fit **per-element types** (e.g. `H`, `O`) or **per-atom** (each atom its own type) for the self-fit?

- **Per-element type** matches your real goal (generalization).
- **Per-atom** is the strongest sanity check (should fit almost perfectly even with limited configs).

## Status
- **Completed:** Reusable helpers added to [pyBall/OCL/Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0) (GridFF upload once; GPU sampling; mock reference knobs).
- **Pending:** Minimal test script + running/validating it end-to-end.

## Finished: test implemented and run successfully

I completed the minimal **FDBM mock-fit pipeline** that **reuses the existing NaCl GridFF** (`Bspline_PLQd.npy` + meta) and samples channels on GPU (OpenCL) without recomputing the grid.

### What I added
- **Test script:**  
  [tests/tMMFF/test_fdbm_fit_gridff_mock.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_fdbm_fit_gridff_mock.py:0:0-0:0)
- **Reusable helpers (already added earlier):** in  
  [pyBall/OCL/Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0)  
  ([init_gridff_sampler_md](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:103:0-122:19), [sample_gridff_channels_rigid](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:125:0-148:13), [fdbm_build_feature_matrix](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:151:0-167:15), [fdbm_make_mock_reference](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:170:0-208:23), …)

# How the test works (high level)
1. Loads `H2O.xyz` (with charges).
2. Loads `NaCl_1x1_L3/Bspline_PLQd.npy` and meta JSON.
3. Initializes [MolecularDynamics.initGridFF()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:1548:4-1619:53) **once** (uploads Bspline grid once).
4. Generates configurations:
   - 800 random rigid transforms
   - 66-point **z-scan** at `(x,y)=(1,1)` for plotting
5. Samples **3 channels** by repeatedly calling `sampleGridFF_Bspline_points` with different `PLQH`:
   - Pauli: `(1,0,0,0)`
   - London: `(0,1,0,0)`
   - Coulomb: `(0,0,1,0)`
6. Generates mock reference energies `E_ref` by composing the same sampled channels with “true” coefficients.
7. Fits linear least squares for `P_type` and `L_type` (Coulomb fixed by charges).
8. Saves plots + `.npz`.

# Where are the results?
All outputs are in:

[/home/prokophapala/git/FireCore/tests/tMMFF/out_fdbm_mock/](cci:9://file:///home/prokophapala/git/FireCore/tests/tMMFF/out_fdbm_mock:0:0-0:0)

Files produced:
- **Data + fitted coeffs:**  
  `fdbm_mock_fit_results.npz`
- **Parity plot:**  
  `parity.png`
- **Z-scan reference vs fit:**  
  `zscan.png`

# Are the results reasonable?
Yes — for the default “self-fit, no perturbation” case they are **exact**, which is what we want as a baseline sanity check.

From the run output:

- **Recovered coefficients exactly**
  - O: `P_fit=1.300000` vs `P_true=1.300000`, `L_fit=0.800000` vs `L_true=0.800000`
  - H: `P_fit=0.600000` vs `P_true=0.600000`, `L_fit=0.200000` vs `L_true=0.200000`

- **Fit error essentially zero**
  - `RMSE=6.36e-17 eV` (numerical precision)

So the **linear algebra + sampling + feature build is correct**.

Also the generated reference energy range printed:
- `E_ref range=[-0.032552,-0.000004] eV`
This is small in magnitude but consistent with “sampling above surface” and with arbitrary scaling in GridFF channels; the important thing is the pipeline correctness and ability to perturb.

# Notes / small warnings seen (not fatal)
- `RepeatedKernelRetrieval` warning from [MolecularDynamics.run_sampleGridFF_Bspline_points()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:953:4-967:27)  
  This is a performance warning (kernel retrieved repeatedly). The test still uses GPU and is correct; later we can optimize by caching the kernel handle inside [MolecularDynamics](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/MolecularDynamics.py:59:0-2007:21) to avoid repeated retrieval overhead.
- Some “skipping getNonBond…” warnings during kernel setup are expected because we’re **only using the GridFF sampling kernel**, not full MMFF/nonbond machinery.

# Next recommended action
Tell me what perturbation sweep you want first, and I’ll extend the test to run a small matrix of cases and save additional plots:

- **`pauli_alpha`**: e.g. `0.95, 1.00, 1.05`
- **Noise**: e.g. `noise_sigma=1e-4, 1e-3 eV`
- **London off**: `use_london=False` both in reference generation and fitting

## Status
- **Complete:** reusable functions added to [Surface_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/Surface_utils.py:0:0-0:0), test created, test run, results saved and validated.

---

# CRITICAL: REQ/PLQ Convention for GridFF

**ElementTypes.dat stores EvdW (energy), but GridFF generation uses sqrt(EvdW) in REQ.y.**

When converting REQ parameters to PLQ coefficients for GridFF sampling, you MUST sqrt the E value from ElementTypes.dat.

## Why this convention?

The GridFF channels contain substrate contributions with sqrt(Ej) baked in:
- Pauli channel: `sum_j sqrt(Ej) * exp(-2a(r-Rj))`
- London channel: `sum_j -2*sqrt(Ej) * exp(-a(r-Rj))`

When multiplied by probe PLQ coefficients, this reconstructs the mixed interaction:
- `P_probe * Pauli = [exp(2aRi) * sqrt(Ei)] * [sum_j sqrt(Ej) * exp(-2a(r-Rj))]`
- `= sum_j sqrt(Ei*Ej) * exp(-2a(r-(Ri+Rj)))`

This gives the correct `Eij = sqrt(Ei * Ej)` mixing rule.

## Implementation

**GridFF generation** (`pyBall/tests/ocl_GridFF_new.py::make_atoms_arrays`):
```python
REQs[:,1] = sqrt(REvdW[:,1])  # Convert EvdW to sqrt(EvdW)
```

**Probe PLQ conversion** (`pyBall/OCL/RigidBodyDynamics::_reqs_to_plq`):
```python
# Expects REQ.y to be sqrt(EvdW), NOT raw EvdW
e = exp(alpha * R)
cL = e * E              # London coefficient
cP = e * cL = e^2 * E   # Pauli coefficient
```

**When reading from ElementTypes.dat for diagnostics:**
```python
req = _get_element_REQ('O', element_types_path)  # Returns (R, EvdW, Q, H)
req = (req[0], np.sqrt(req[1]), req[2], req[3])  # MUST sqrt E
plq = _req_to_plq(req, alpha=alpha_morse)
```

## Files with this convention documented

- `tests/tMMFF/test_fdbm_fit_gridff_mock.py` - _get_element_REQ, _req_to_plq
- `pyBall/OCL/RigidBodyDynamics.py` - _reqs_to_plq
- `pyBall/tests/ocl_GridFF_new.py` - make_atoms_arrays
- `pyBall/OCL/Surface_utils.py` - import section

---

# Session Report: REQ/PLQ Documentation and Diagnostic Plot Improvements

## Summary

This session focused on documenting the critical REQ/PLQ convention for GridFF and improving the diagnostic plotting in the FDBM mock-fit test.

## REQ/PLQ Convention Documentation

Added comprehensive documentation about the `EvdW` vs `sqrt(EvdW)` convention to:

- `tests/tMMFF/test_fdbm_fit_gridff_mock.py` - Added detailed docstrings to `_get_element_REQ` and `_req_to_plq` functions explaining that ElementTypes.dat stores raw EvdW but GridFF expects sqrt(EvdW) in REQ.y
- `pyBall/OCL/RigidBodyDynamics.py` - Added detailed docstring to `_reqs_to_plq` function with formula derivation
- `pyBall/tests/ocl_GridFF_new.py` - Added docstring to `make_atoms_arrays` explaining the `bSqrtEvdw` parameter
- `pyBall/OCL/Surface_utils.py` - Added critical comment in import section
- `tests/tMMFF/test_gridff_alignment.py` - Added documentation to `test_system` function
- `tests/tMMFF/test_folded_fit_nacl1x1.py` - Added comment in import section
- `doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md` - Added dedicated section with formula derivation
- `doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md` - Added dedicated section with cross-references

## Diagnostic Plot Improvements

### Layout Expansion
- Changed from 2x2 to 2x3 layout
- Added XZ slice at origin (iy=0) to show both Na and Cl ions
- Added Z-profile at origin (cx=0, cy=0) for consistency

### Visual Consistency
- Changed all labels to use 3 decimal digits (`.3f` format)
- Fixed XZ slice vmax calculation to use only slice data above zmin, not whole grid
- Added consistent z-range limits to all XZ slices and Z-profiles (`g0[2]` to `g0[2] + min(nz*dz, 20)`)

### Code Refactoring
- Created reusable `plot_xz_slice(ax, iy, title_suffix)` function for XZ slices
- Created reusable `plot_1d_profile(ax, axis, fixed_idx, title_suffix)` function for 1D profiles
- Eliminated code duplication across plotting functions

### Atom Display Consistency
- Precomputed `atoms_colors` once at function start
- Removed all spatial filtering from XZ slices and 1D profiles
- All panels now show the exact same atom arrays (no on-the-fly processing)
- XZ slices: plot all `atoms_xyz[:, 0], atoms_xyz[:, 2]` with `atoms_colors`
- 1D profiles: iterate through all atoms with `atoms_colors`
- Molecule samples: plot all `mol_apos_flat` without y-filtering

## Configuration Generation Improvements

Added CLI options to `test_fdbm_fit_gridff_mock.py`:
- `--xy-range` (default -2.0 2.0 A) - symmetric XY range around origin
- `--z-range` (default -1.0 5.0 A) - includes repulsion region
- Z-scan now uses the same z-range as random configs
- Energy range increased from -0.028 to -0.216 eV, confirming samples now reach repulsion region

## Key Files Modified

- `tests/tMMFF/test_fdbm_fit_gridff_mock.py` - REQ/PLQ documentation, CLI options, config generation
- `pyBall/OCL/Surface_utils.py` - REQ/PLQ comment, plotting refactoring, atom display consistency
- `pyBall/OCL/RigidBodyDynamics.py` - REQ/PLQ documentation
- `pyBall/tests/ocl_GridFF_new.py` - REQ/PLQ documentation
- `tests/tMMFF/test_gridff_alignment.py` - REQ/PLQ documentation
- `tests/tMMFF/test_folded_fit_nacl1x1.py` - REQ/PLQ documentation
- `doc/Topics/OnSurfaceAssembly/GridFF_FDBM_Fitting.md` - REQ/PLQ documentation, session report
- `doc/Topics/OnSurfaceAssembly/GridFF_atoms_alignment.md` - REQ/PLQ documentation

---

# USER

I'm interested in fitting GridFF/FDBM using DFT references.
read these documents to know what we already have
@GridFF_FDBM_Fitting.md @RigidSurfPotential_GridFF.md @test_folded_fit_nacl1x1.py @Surface_utils.py 

I make link to folfer with DFT data here
/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New
plese investigrate it and tell me where is what, how it is structuree, be warry there is huge number of files

Now one importaint thing is we have 1x1 substrate prepared, we want to use it, but the DFT was actually computed with a large cell 
in particular 
/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/4-elec_pot_chr_den/
there is 
charge_density.cube
electrostatic_potential.cube
NaCl.xyz

We need to see how the molecules are sampled on this substrate before we start fitting. Can you do this?
you may also review @test_gridff_alignment.py 

Do not program anything yet. Revie all the relevant files and make detailed notes where it what and how to use it, 

---

# SWE-1.6

# DFT Data Structure Analysis: small_mols_NaCl_New

## Overview

The dataset contains DFT interaction energies for 11 small molecules on NaCl(001) surface, with systematic sampling of positions and orientations. Total dataset: ~52k configurations per molecule-site.

## Folder Structure

### 0-geoms/ - Molecular Geometries
- **Molecule XYZ files**: H2O, NH3, H2CO, HCONH2, HCOOH (neutral geometries)
- **Oriented molecule-site files**: H2O-H, H2O-O, NH3-H, NH3-N, etc. (H-bond direction aligned to +x axis, active atom at origin)
- **NaCl.xyz**: 32-atom 2-layer NaCl slab (11.2826 × 11.2826 × 40.0 Å)
- **molecules.dat**: Metadata for each molecule-site (flag, atom indices for orientation)
- **rotate.f90/sh**: Fortran script to generate oriented molecule-site geometries

### 1-inputs/ - Configuration Generation
- **assembly.f90**: Fortran program to place molecules on substrate
- **run.sh**: Bash driver that loops over:
  - 11 molecule-site combinations
  - 28 xy-positions (triangular grid: ix=0..6, iy=0..ix for first half, symmetric for second)
  - 33 z-distances (1.40-10.00 Å)
  - 3 xy-orientations (Na, Cl, hollow)
  - 3 z-orientations (-45°, 0°, 45°)
- **confs/**: Generated XYZ configurations (52,272 per molecule-site = 575k total)

### 2-sp/ - Single-Point Calculations
- **confs/**: DFT input/output for each configuration (313,687 items)
- **Scripts**: Setup and run DFT jobs (likely ASE/CP2K based on cube file format)

### 3-results/ - Energy Results
- **Organized by molecule-site**: e.g., H2O-H/, H2O-O/, NH3-N/
- **Data files**: ixX_iyY-orientxy_XY-orientz_Z.dat (z-distance vs energy curves)
- **Format**: 2-column text (z [Å], E [eV])
- **Example**: ix0_iy0-orientxy_Na-orientz_0.dat for H2O-H shows dissociation curve from 1.40-10.00 Å
- **plot.py**: Generates HBond.png comparing H2O-H, H2O-O, NH3-N on Na/Cl sites
- **res.sh**: Extracts energies from 2-sp/ results

### 4-elec_pot_chr_den/ - Substrate DFT Grids
- **NaCl.xyz**: Large-cell substrate (same as 0-geoms but with different comment format)
- **charge_density.cube**: 101 MB, 128×128×480 grid, electron density from DFT
- **electrostatic_potential.cube**: 99 MB, 128×128×480 grid, electrostatic potential from DFT
- **Grid parameters**: 
  - Origin: (0, 0, 0) Å
  - Spacing: (0.166570, 0.166570, 0.157477) Å
  - Extent: 21.32 × 21.32 × 75.59 Å
  - 32 atoms listed in cube header
- **Other files**: ASE/CP2K outputs (gs.gpw, gs.out, POSCAR, etc.)

### opt_slab/ - Substrate Optimization
- **input.xyz**: Initial NaCl geometry
- **opt.py**: Substrate optimization script

## Key Observations

### Substrate Cell Mismatch
- **1x1 substrate available**: cpp/common_resources/xyz/NaCl_1x1_L3.xyz (used in existing GridFF tests)
- **DFT computed with large cell**: 11.2826 × 11.2826 × 40.0 Å (same as 1x1 but with more vacuum)
- **Cube files cover larger area**: 21.32 × 21.32 Å (2×2 unit cells in XY)
- **Implication**: Need to either:
  1. Fold/periodize cube data to 1x1 cell, or
  2. Use larger substrate for fitting, or
  3. Extract 1x1 region from cube with proper PBC handling

### Configuration Sampling Strategy
The sampling is highly systematic:
- **XY positions**: Triangular grid covering irreducible wedge of 2D unit cell (ix=0..6, iy=0..ix)
- **Z distances**: Dense sampling near surface (1.40-3.00 Å in 0.05-0.10 Å steps), sparse at long range
- **XY orientations**: 3 high-symmetry sites (Na, Cl, hollow)
- **Z orientations**: 3 tilt angles (-45°, 0°, 45°)

This is ideal for fitting - provides dense coverage of binding wells and long-range asymptotes.

### Energy Data Format
Each .dat file contains:
- 33 z-points × 3 repeats = 99 lines (duplicates for averaging?)
- Energy values are absolute DFT energies (~ -124 eV for H2O-H)
- Need to subtract reference (isolated molecule + slab) to get interaction energies

## Reusable Components from Existing Code

### Surface_utils.py Functions
- **load_gridff_array()**: Load .npy GridFF with validation
- **load_gridff_metadata()**: Load JSON metadata for grid parameters
- **load_substrate_xyz_with_lvec()**: Parse XYZ with lattice vectors
- **sample_gridff_channels_rigid()**: Sample GridFF at atom positions for many configs
- **fdbm_build_feature_matrix()**: Build linear fitting matrix from sampled potentials
- **fdbm_make_mock_reference()**: Generate mock reference energies (for testing)
- **plot_gridff_diagnostics()**: Visualize grids with atom overlay

### test_folded_fit_nacl1x1.py Patterns
- **Rigid transform packing**: pack_transforms() converts XYZ to 3×4 matrices
- **Component evaluation**: eval_rigid_getSurfMorse_components() for Pauli/London/Coulomb
- **Linear fitting**: Uses scipy.optimize.lsq_linear with bounds
- **1D scan visualization**: plot_1d_lines() for diagnostic curves
- **RMSE calculation**: rmse() and rmse_masked() for validation

### test_gridff_alignment.py
- **Convention testing**: Tests multiple grid origin conventions
- **Alignment verification**: Checks grid-atom correspondence
- **Diagnostic plots**: XY/XZ slices with atom overlays

## Critical Issue: Cube File vs 1x1 Substrate

The cube files cover a 2×2 unit cell area (21.32 × 21.32 Å), but the 1x1 substrate is 11.2826 × 11.2826 Å. Options:

1. **Extract 1x1 region**: Take ix=0..63, iy=0..63 from cube (assuming proper alignment)
2. **Use 2×2 substrate**: Create NaCl_2x2_L3.xyz and use larger GridFF
3. **Fold with PBC**: Average over 2×2 periodic images to get 1x1 effective potential

Need to verify cube grid alignment with substrate lattice vectors before proceeding.

## Next Steps for Fitting

1. **Parse cube files**: Write cube loader (or use ASE) to get grid data and metadata
2. **Verify alignment**: Check cube origin/spacing against NaCl lattice
3. **Map configurations**: Load XYZ configs from 1-inputs/confs/ or regenerate using existing code
4. **Sample potentials**: Use Surface_utils.py functions to sample cube at atom positions
5. **Load energies**: Parse 3-results/*.dat files to get DFT interaction energies
6. **Build feature matrix**: Use fdbm_build_feature_matrix() pattern
7. **Fit coefficients**: Linear least squares for P_i, L_i per atom type

## Data Volume Considerations

- **Cube files**: 100 MB each (manageable)
- **Configurations**: 575k XYZ files (huge - don't load all at once)
- **Results**: 144 .dat files per molecule-site (manageable)
- **Strategy**: Stream configurations, pre-sample potentials in batches, cache intermediate results

## Summary: DFT Data Structure Analysis

I've reviewed the documentation and investigated the DFT data folder structure. Here are the key findings:

### Folder Organization

**0-geoms/**: Molecular geometries and substrate
- Neutral molecules (H2O, NH3, H2CO, HCONH2, HCOOH)
- Oriented molecule-site files (H2O-H, H2O-O, etc.) with H-bond aligned to +x
- NaCl.xyz: 32-atom 2-layer slab (11.2826 × 11.2826 × 40.0 Å)
- molecules.dat: Orientation metadata

**1-inputs/**: Configuration generation
- assembly.f90 + run.sh generate 575k configurations
- Systematic sampling: triangular XY grid (ix=0..6), 33 z-distances (1.40-10.00 Å), 3 XY orientations (Na/Cl/hollow), 3 Z tilts (-45°, 0°, 45°)

**2-sp/**: Single-point DFT calculations (313k items)

**3-results/**: Energy results
- Organized by molecule-site (e.g., H2O-H/, H2O-O/)
- .dat files: z-distance vs energy curves (33 points × 3 repeats)
- Example: ix0_iy0-orientxy_Na-orientz_0.dat shows dissociation from -122 to -124 eV

**4-elec_pot_chr_den/**: Substrate grids
- charge_density.cube: 101 MB, 128×128×480 grid
- electrostatic_potential.cube: 99 MB, 128×128×480 grid
- Grid spacing: (0.166570, 0.166570, 0.157477) Å
- Extent: 21.32 × 21.32 × 75.59 Å (2×2 unit cells in XY)

### Critical Issue: Cell Size Mismatch

The cube files cover a **2×2 unit cell area** (21.32 × 21.32 Å), but the available 1x1 substrate is 11.2826 × 11.2826 Å. Options:
1. Extract 1x1 region from cube (ix=0..63, iy=0..63)
2. Create 2×2 substrate (NaCl_2x2_L3.xyz)
3. Fold with PBC to get 1x1 effective potential

### Reusable Components

From Surface_utils.py:
- [load_gridff_array()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/Surface_utils.py:72:0-92:84), [load_gridff_metadata()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/Surface_utils.py:44:0-67:19) for grid I/O
- [sample_gridff_channels_rigid()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/Surface_utils.py:129:0-152:13) for batch sampling
- [fdbm_build_feature_matrix()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/Surface_utils.py:155:0-171:15) for linear fitting matrix
- [fdbm_make_mock_reference()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/Surface_utils.py:174:0-212:23) for mock data generation

From test_folded_fit_nacl1x1.py:
- [pack_transforms()](cci:1://file:///home/prokop/git/FireCore/tests/tMMFF/test_folded_fit_nacl1x1.py:37:0-42:28) for rigid transform packing
- `eval_rigid_getSurfMorse_components()` for component evaluation
- Linear fitting with scipy.optimize.lsq_linear
- 1D scan visualization patterns

### Next Steps

1. Parse cube files (or use ASE) to get grid data
2. Verify cube alignment with NaCl lattice
3. Map configurations from 1-inputs/confs/ or regenerate
4. Sample cube potentials at atom positions
5. Load energies from 3-results/*.dat
6. Build feature matrix and fit P_i, L_i coefficients

The dataset is well-structured for fitting - dense sampling near binding wells, systematic coverage of orientations, and manageable energy file format. The main challenge is handling the 2×2 vs 1×1 cell mismatch.

---

# USER

I think you are wrong. The CUB files cover the provided unit cell 11.2826 × 11.2826 A, maybe wha you read (  21.32 × 21.32 × 75.59 ) is in bohr radius rather than Angstroem? cube files are in bohr I think ?

this is what I get from VESTA when loading electrostatic_potential.cube

```
OpenGL version: 4.6.0 NVIDIA 550.120
Video configuration: NVIDIA GeForce RTX 3090/PCIe/SSE2
Maximum supported width and height of the viewport: 32768 x 32768
OpenGL depth buffer bit: 24


/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/4-elec_pot_chr_den/electrostatic_potential.cube
====================================================================================
Title               Cube file from ASE, written on Mon May 18 15:03:44 2026
Dimensions           128 128 480

Lattice parameters

   a        b        c       alpha    beta     gamma
11.28257 11.28257 39.99995  90.0000  90.0000  90.0000

Unit-cell volume = 5091.846058 Å^3

Total number of polygons and unique vertices on slices;
      (1 0 0):      4136 (     2926),      4136 (     2926)
      (0 1 0):      4136 (     2926),      4136 (     2926)
      (0 0 1):         0 (        0),         0 (        0)
====================================================================================

====================================================================================
Title               Cube file from ASE, written on Mon May 18 15:03:44 2026

Lattice type        P
Space group name    P 1
Space group number  1
Setting number      1

Lattice parameters

   a        b        c       alpha    beta     gamma
 1.00000  1.00000  1.00000  90.0000  90.0000  90.0000

Unit-cell volume = 1.000000 Å^3

Structure parameters

                       x          y          z          Occ.     U    Site      Sym.
   1 Na  Na1         0.00000    0.00000   10.00000    1.000    0.000    1          -
   2 Cl  Cl1         2.82065    0.00000   10.00000    1.000    0.000    1          -
   3 Na  Na2         5.64130    0.00000   10.00000    1.000    0.000    1          -
   4 Cl  Cl2         8.46195    0.00000   10.00000    1.000    0.000    1          -
   5 Cl  Cl3         0.00000    2.82065   10.00000    1.000    0.000    1          -
   6 Na  Na3         2.82065    2.82065   10.00000    1.000    0.000    1          -
   7 Cl  Cl4         5.64130    2.82065   10.00000    1.000    0.000    1          -
   8 Na  Na4         8.46195    2.82065   10.00000    1.000    0.000    1          -
   9 Na  Na5         0.00000    5.64130   10.00000    1.000    0.000    1          -
  10 Cl  Cl5         2.82065    5.64130   10.00000    1.000    0.000    1          -
  11 Na  Na6         5.64130    5.64130   10.00000    1.000    0.000    1          -
  12 Cl  Cl6         8.46195    5.64130   10.00000    1.000    0.000    1          -
  13 Cl  Cl7         0.00000    8.46195   10.00000    1.000    0.000    1          -
  14 Na  Na7         2.82065    8.46195   10.00000    1.000    0.000    1          -
  15 Cl  Cl8         5.64130    8.46195   10.00000    1.000    0.000    1          -
  16 Na  Na8         8.46195    8.46195   10.00000    1.000    0.000    1          -
  17 Cl  Cl9         0.00000    0.00000   12.82065    1.000    0.000    1          -
  18 Na  Na9         2.82065    0.00000   12.82065    1.000    0.000    1          -
  19 Cl  Cl10        5.64130    0.00000   12.82065    1.000    0.000    1          -
  20 Na  Na10        8.46195    0.00000   12.82065    1.000    0.000    1          -
  21 Na  Na11        0.00000    2.82065   12.82065    1.000    0.000    1          -
  22 Cl  Cl11        2.82065    2.82065   12.82065    1.000    0.000    1          -
  23 Na  Na12        5.64130    2.82065   12.82065    1.000    0.000    1          -
  24 Cl  Cl12        8.46195    2.82065   12.82065    1.000    0.000    1          -
  25 Cl  Cl13        0.00000    5.64130   12.82065    1.000    0.000    1          -
  26 Na  Na13        2.82065    5.64130   12.82065    1.000    0.000    1          -
  27 Cl  Cl14        5.64130    5.64130   12.82065    1.000    0.000    1          -
  28 Na  Na14        8.46195    5.64130   12.82065    1.000    0.000    1          -
  29 Na  Na15        0.00000    8.46195   12.82065    1.000    0.000    1          -
  30 Cl  Cl15        2.82065    8.46195   12.82065    1.000    0.000    1          -
  31 Na  Na16        5.64130    8.46195   12.82065    1.000    0.000    1          -
  32 Cl  Cl16        8.46195    8.46195   12.82065    1.000    0.000    1          -
====================================================================================

Number of polygons and unique vertices on isosurface = 84864 (43588)
32 atoms, 64 bonds, 16 polyhedra; CPU time = 80 ms

```

---

# SWE-1.6

You're absolutely right. The cube file spacing values are in Bohr, not Angstrom. Let me recalculate:

- Grid spacing in cube header: (0.166570, 0.166570, 0.157477) Bohr
- Convert to Angstrom: 0.166570 × 0.529177 = 0.0881 Å, 0.157477 × 0.529177 = 0.0833 Å
- Grid extent: 128 × 0.0881 = 11.28 Å, 480 × 0.0833 = 40.0 Å

This matches the 1x1 substrate cell exactly (11.2826 × 11.2826 × 40.0 Å). **No cell mismatch issue** - the cube files are already in the correct 1x1 unit cell.

## Corrected Analysis

### 4-elec_pot_chr_den/ - Substrate Grids (Corrected)
- **charge_density.cube**: 101 MB, 128×128×480 grid
- **electrostatic_potential.cube**: 99 MB, 128×128×480 grid
- **Grid spacing**: (0.166570, 0.166570, 0.157477) Bohr = (0.0881, 0.0881, 0.0833) Å
- **Extent**: 11.28 × 11.28 × 40.0 Å (matches 1x1 substrate exactly)
- **Format**: Cube files use Bohr for grid spacing (standard convention)

### Implications for Fitting

This simplifies the workflow significantly:
1. **No cell conversion needed** - cube grids align with 1x1 substrate
2. **Direct sampling** - can sample cube potentials at atom positions using trilinear interpolation
3. **Existing GridFF infrastructure** - can reuse 1x1 substrate from cpp/common_resources/xyz/NaCl_1x1_L3.xyz

### Updated Next Steps

1. Parse cube files (convert Bohr spacing to Angstrom for internal use)
2. Verify cube origin alignment with substrate (should be at (0,0,0) based on header)
3. Load configurations from 1-inputs/confs/ or regenerate using assembly.f90
4. Sample cube potentials at atom positions (need trilinear interpolation for cube data)
5. Load energies from 3-results/*.dat
6. Build feature matrix and fit P_i, L_i coefficients

The cube file format uses Bohr, but VESTA and other tools handle the conversion automatically. For the fitting code, we'll need to convert the grid spacing to Angstrom when setting up the interpolation.

---

# USER

OK, so can we use @test_folded_fit_nacl1x1.py or @test_gridff_alignment.py with minimal changes to visualize these potential, and later eventually to plot our sampling data for example for H2O in that potential?

/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/3-results/H2O-H/
/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/3-results/H2O-O/

before we do that, we should perhaps convert the data into .npy or .npz so we can work with that easily.

look we have this file
/home/prokop/git/FireCore/cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy
which contains 4 components [nx,ny,nz,4] or maybe [nz,ny,nx,4] not sure which ordering chek it) 

the 4 compoennts are. {Pauli,London,Electrostatics,HydrogenBond}, HydrogenBond is zero for now. London we can compute by C6/r^6 resp by the atractive part of getSurfMorse in @Surface.cl 
what is key i Pauli and Electrostatics. Electrostatics is directly  electrostatic_potential.cube which should be direcly mulitplied by atomic charge Q_i to get interaction energy with proper atom (just in proper units, we work in eV cube file may be in Hartree? But not sure perhaps not)
The pauli is more nuanced. We can take the density from charge_density.cube as some estimate of Pauli repulsion. Just maybe it should be rised to some powere betaPauli, and multiplied by the coefficient for atom P_i, whre i is index of atom of molecule (absorbate), but lets keep betaPauli=1 for the moment.

E_pauli = sum_i[ P_i rho(r_i)^beta ] 

this can be rigorously derived from original formula of FDBM 

E_pauli(R) =  integral [  (rho_substrater) * rho_molecule(r+R)) ^beta ]
by cosidering  that molecule is compsed of atom i 
rho_molecule(r) = sum_i rho_i(r)
and that these atomic densityes are basically exponential  rho_i(r) = A_i * exp( b|r-r_i|)
now the procuct of exponential function can be extracted R_i which is effective radius of an atom i , so P_i = A_i * exp(b*R_i)
this derivation is just for mathermaticall baground, the point is to justify we do not need to do convolution of the exponential density of each atom, as this convolution is basically encoded in the scaling coeficient P_i, so we only need to fit P_i. Is that clear? Plese write this more mathematically reigorously.

The final point is 
1) we build  4 components [nx,ny,nz,4] Bspline_PLQd.npy but putting electron dnesity in slot 0 (Pauli) and electrostatic potential in slot 3, slot 2 (Hbond) is zero, and slot 2 (Lndon) is less importaint but can be computed from  C6/r^6 or using getSurfMorse kernel. 

2) then we will fit coeffecits P_i asuming constant Q_i given by resp charges. Where can we find those reps charges? I will progide you better reps charge once I have them but provisionally lets just use these -0.4e for O and +0.2e for H @H2O.xyz 

so this would be like gradial fitting process. 
1) step one we fully focus on P_i for the types, currently we have sjut H2O and types are simply H and O, the two hydrogens in H2O molecule must have same P
2) then we can fitt also charges Q_i 
3) then we can compute the London part by getSurfMorse and fit L_i coefs

It would be also good to see if Pauli computed by getSurfMorse is roughly similar to Pauli computed from electron density. MAybe good would be to plot it in log-scale to see the decay fastro  as linear slope. This is something we did in AFM fitting of FDBM, here we generalize it form molecule composed of more atom in AFM we assume single proble particle (single atom)
@AFM_FDBM_fitting.md 

These are scattered ideas, please review the referenced files and write this more systematically, rigrously, and in more organized manner

---

# SWE-1.6

# GridFF/FDBM Fitting with DFT References: Systematic Plan

## Overview

Fit GridFF/FDBM parameters using DFT reference data for small molecules on NaCl(001). The approach uses DFT-computed electron density and electrostatic potential as substrate grids, then fits atom-type coefficients (P_i for Pauli, Q_i for electrostatics, L_i for London) to match DFT interaction energies.

## Mathematical Foundation

### FDBM Pauli Repulsion from Density Convolution

The Full Density Based Model (FDBM) expresses Pauli repulsion as a convolution of substrate and adsorbate electron densities:

```
E_pauli(R) = A_pauli * ∫ [ρ_substrate(r) * ρ_adsorbate(r + R)]^β dr
```

where:
- ρ_substrate(r) is the substrate electron density
- ρ_adsorbate(r) is the adsorbate electron density
- R is the adsorbate position relative to substrate
- β is a power law exponent (typically β ≈ 1)
- A_pauli is a global scaling coefficient

### Simplification to Atom-Type Coefficients

Assume the adsorbate is composed of atoms i with atomic densities:

```
ρ_adsorbate(r) = Σ_i ρ_i(r - r_i)
```

where ρ_i(r) is the electron density of atom i centered at r_i. For exponential atomic densities:

```
ρ_i(r) = A_i * exp(-b_i * |r - r_i|)
```

The product of two exponential functions is also exponential:

```
ρ_substrate(r) * ρ_i(r - r_i) = ρ_substrate(r) * A_i * exp(-b_i * |r - r_i|)
```

The convolution integral can be separated into a substrate-dependent part and an atom-dependent part. Define the substrate sampling function:

```
S_i(R) = ∫ ρ_substrate(r) * exp(-b_i * |r - (R + r_i)|) dr
```

Then the Pauli energy becomes:

```
E_pauli(R) = A_pauli * Σ_i [A_i * S_i(R)]^β
```

For β = 1 (linear case), this simplifies to:

```
E_pauli(R) = Σ_i P_i * ρ_substrate(R + r_i)
```

where P_i = A_pauli * A_i is an atom-type coefficient. The key insight is that the convolution of the exponential atomic density with the substrate density is encoded in the coefficient P_i, so we only need to sample the substrate density at atom positions and fit P_i.

For general β, we have:

```
E_pauli(R) = Σ_i P_i * [ρ_substrate(R + r_i)]^β
```

where P_i = A_pauli * A_i^β.

### Electrostatic Interaction

The electrostatic interaction is linear in charge:

```
E_es(R) = Σ_i Q_i * φ_substrate(R + r_i)
```

where:
- Q_i is the partial charge of atom i
- φ_substrate(r) is the substrate electrostatic potential

### London Dispersion

The London (van der Waals) interaction is typically modeled as:

```
E_london(R) = Σ_i L_i * V_london(R + r_i)
```

where V_london can be computed from C6/r^6 or extracted from the attractive part of a Morse potential.

## Data Structure

### DFT Reference Data

**Location:** `/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/`

**Structure:**
- `0-geoms/`: Molecular geometries (H2O, NH3, H2CO, HCONH2, HCOOH) and substrate (NaCl.xyz)
- `1-inputs/confs/`: Generated configurations (575k total, 52k per molecule-site)
- `3-results/`: DFT energies organized as `ixX_iyY-orientxy_XY-orientz_Z.dat` files
  - Format: 2-column text (z [Å], E [eV])
  - Each file contains z-scan for specific (ix, iy, orientxy, orientz)
- `4-elec_pot_chr_den/`: Substrate grids
  - `charge_density.cube`: 128×128×480 grid, electron density
  - `electrostatic_potential.cube`: 128×128×480 grid, electrostatic potential
  - Grid spacing: (0.166570, 0.166570, 0.157477) Bohr = (0.0881, 0.0881, 0.0833) Å
  - Extent: 11.28 × 11.28 × 40.0 Å (matches 1x1 NaCl substrate)

### Existing GridFF Format

**Location:** `/home/prokop/git/FireCore/cpp/common_resources/NaCl_1x1_L3/Bspline_PLQd.npy`

**Format:** [nx, ny, nz, 4] array with 4 channels:
- Channel 0: Pauli repulsion
- Channel 1: London dispersion
- Channel 2: Coulomb electrostatics
- Channel 3: Hydrogen bond (zero for NaCl)

## Implementation Plan

### Phase 1: Data Conversion and Visualization

**Goal:** Convert cube files to GridFF-compatible format and verify alignment.

**Steps:**

1. **Cube file parser**
   - Write function to read cube files (header + volumetric data)
   - Convert Bohr spacing to Angstrom
   - Extract grid metadata (origin, spacing, dimensions)
   - Output: numpy array [nz, ny, nx] + metadata dict

2. **GridFF format conversion**
   - Create 4-channel array [nx, ny, nz, 4]
   - Channel 0 (Pauli): electron density from charge_density.cube
   - Channel 1 (London): initialize to zero (compute later)
   - Channel 2 (Coulomb): electrostatic potential from electrostatic_potential.cube
   - Channel 3 (Hbond): zero
   - Save as .npy with metadata JSON

3. **Alignment verification**
   - Use `test_gridff_alignment.py` pattern
   - Load NaCl substrate from `0-geoms/NaCl.xyz`
   - Visualize grid slices with atom overlay
   - Verify cube origin aligns with substrate atoms
   - Check grid spacing matches lattice vectors

4. **Visualization tools**
   - Adapt `test_folded_fit_nacl1x1.py` plotting functions
   - Plot XY, XZ, YZ slices of density and potential
   - Overlay substrate atoms for visual verification
   - Compare with existing Bspline_PLQd.npy for reference

### Phase 2: Pauli Coefficient Fitting (β = 1)

**Goal:** Fit P_i coefficients for H and O atom types using density-based Pauli.

**Assumptions:**
- β_pauli = 1 (linear in density)
- Fixed charges: Q_O = -0.4e, Q_H = +0.2e (from H2O.xyz)
- London term ignored initially
- Hydrogen bond term zero

**Steps:**

1. **Load DFT reference data**
   - Parse `3-results/H2O-H/*.dat` and `3-results/H2O-O/*.dat`
   - Extract z vs E curves for each configuration
   - Convert to interaction energy: E_int = E_total - E_slab - E_molecule
   - Need reference energies for isolated slab and molecule

2. **Load configurations**
   - Load XYZ files from `1-inputs/confs/` or regenerate using assembly.f90
   - Extract atom positions and types
   - Apply rigid transforms (rotation, translation) as needed

3. **Sample substrate potentials**
   - Use `Surface_utils.py` functions:
     - `load_gridff_array()` to load converted cube data
     - `sample_gridff_channels_rigid()` to sample at atom positions
   - For each configuration:
     - Sample electron density at each atom position
     - Sample electrostatic potential at each atom position
   - Output: per-atom density and potential values

4. **Build feature matrix**
   - Use `fdbm_build_feature_matrix()` pattern
   - For each configuration c and atom i:
     - Feature_i = ρ(r_i) (density)
     - Target = E_int(c)
   - Linear model: E_int = Σ_i P_i * ρ(r_i) + Σ_i Q_i * φ(r_i)
   - Since Q_i fixed, electrostatic term is known baseline
   - Subtract electrostatic contribution: E_pauli_ref = E_int - Σ_i Q_i * φ(r_i)

5. **Linear least squares fit**
   - Solve for P_H and P_O (2 parameters for H2O)
   - Use scipy.optimize.lsq_linear with bounds (P_i > 0)
   - Validate with RMSE and visual comparison

6. **Validation**
   - Plot fitted vs DFT curves for representative configurations
   - Check consistency across different sites (Na, Cl, hollow)
   - Verify P_H and P_O are physically reasonable

### Phase 3: Charge Fitting

**Goal:** Fit Q_i coefficients (RESP charges) with P_i fixed from Phase 2.

**Steps:**

1. **Linear model with variable charges**
   - E_int = Σ_i P_i * ρ(r_i) + Σ_i Q_i * φ(r_i)
   - P_i fixed from Phase 2
   - Fit Q_H and Q_O (2 parameters)

2. **Constraints**
   - Charge neutrality: Σ_i Q_i = 0 (for H2O: 2*Q_H + Q_O = 0)
   - Use constrained least squares

3. **Validation**
   - Compare fitted charges with RESP reference
   - Check physical reasonableness

### Phase 4: London Term and L_i Fitting

**Goal:** Compute London potential and fit L_i coefficients.

**Options for London potential:**

1. **C6/r^6 model**
   - Compute pairwise C6/r^6 between adsorbate and substrate atoms
   - Sum over substrate atoms with PBC
   - Fast but approximate

2. **getSurfMorse attractive part**
   - Use `getMorsePLQH()` from Surface.cl
   - Extract London component: EL = -2.0*Ee*PLQH.y
   - More accurate but requires GridFF-style preprocessing

**Steps:**

1. **Compute London potential grid**
   - Choose method (C6/r^6 or getSurfMorse)
   - Generate grid [nx, ny, nz] for London
   - Add to channel 1 of GridFF array

2. **Sample London at atom positions**
   - Use same sampling pipeline as Phase 2
   - Extract V_london(r_i) for each atom

3. **Fit L_i coefficients**
   - Linear model: E_int = Σ_i P_i * ρ(r_i) + Σ_i Q_i * φ(r_i) + Σ_i L_i * V_london(r_i)
   - P_i, Q_i fixed from Phases 2-3
   - Fit L_H and L_O

4. **Validation**
   - Check improvement in RMSE
   - Verify L_i physically reasonable

### Phase 5: Pauli Power Law (β ≠ 1)

**Goal:** Test β_pauli ≠ 1 for better fit, following AFM FDBM methodology.

**Mathematical background:**
- From AFM fitting, optimal β ≈ 0.8-0.9 for DFTB+ basis sets
- β < 1 softens the density decay to match physical force gradient
- Fit via nonlinear least squares: E_pauli = A_pauli * overlap^β

**Steps:**

1. **Nonlinear fit**
   - Model: E_pauli_ref = Σ_i P_i * [ρ(r_i)]^β
   - Fit P_H, P_O, and β simultaneously
   - Use scipy.optimize.curve_fit with log-linear initialization

2. **Compare β=1 vs β≠1**
   - RMSE comparison
   - Visual comparison of decay curves (log-scale plot)
   - Check if β improves fit significantly

3. **Log-scale decay analysis**
   - Plot ln(E_pauli) vs ln(ρ) to check linearity
   - Slope should equal β if power law holds
   - Compare with AFM results (β ≈ 0.85 for mio-1-1)

### Phase 6: Comparison with getSurfMorse Pauli

**Goal:** Compare density-based Pauli with getSurfMorse Pauli for validation.

**Steps:**

1. **Compute getSurfMorse Pauli grid**
   - Use `getMorsePLQH()` from Surface.cl
   - Extract Pauli component: EP = Ee*e*PLQH.x
   - Generate grid for NaCl substrate

2. **Compare decay profiles**
   - Extract z-column at representative atom positions
   - Plot both Pauli profiles on log scale
   - Compare decay rates (slope in log-log plot)
   - Check if density-based Pauli captures similar physics

3. **Fit both to DFT**
   - Fit P_i using getSurfMorse Pauli
   - Compare fitted coefficients and RMSE
   - Determine which method gives better fit

## Reusable Components

### From Surface_utils.py

- `load_gridff_array()`: Load .npy GridFF with validation
- `load_gridff_metadata()`: Load JSON metadata
- `load_substrate_xyz_with_lvec()`: Parse XYZ with lattice vectors
- `sample_gridff_channels_rigid()`: Sample GridFF at atom positions
- `fdbm_build_feature_matrix()`: Build linear fitting matrix
- `plot_gridff_diagnostics()`: Visualize grids with atom overlay

### From test_folded_fit_nacl1x1.py

- `pack_transforms()`: Convert XYZ to 3×4 rigid transforms
- `eval_rigid_getSurfMorse_components()`: Evaluate Pauli/London/Coulomb
- `plot_1d_lines()`: 1D scan visualization
- `rmse()`, `rmse_masked()`: Validation metrics

### From test_gridff_alignment.py

- `run_alignment_verification()`: Grid-atom alignment testing
- Visualization patterns for grid slices

## Key Technical Considerations

### Unit Conversions

- Cube files: spacing in Bohr, convert to Angstrom (1 Bohr = 0.529177 Å)
- Electrostatic potential: may be in Hartree, convert to eV (1 Hartree = 27.2114 eV)
- Charges: in elementary charge units
- Energies: DFT in eV, ensure consistency

### Grid Alignment

- Cube origin at (0, 0, 0) in Bohr
- Substrate atoms in Angstrom
- Need to verify alignment before sampling
- Use visual verification with atom overlay

### PBC Handling

- Cube files cover 1x1 unit cell (11.28 × 11.28 Å)
- Configurations may span multiple unit cells
- Apply PBC when sampling: use fractional coordinates modulo 1

### Data Volume

- 575k configurations total
- Don't load all at once
- Stream in batches, cache intermediate results
- Use efficient numpy operations

## Output Structure

```
fit_dft_gridff/
├── cube_to_gridff.py          # Cube to .npy conversion
├── verify_alignment.py        # Grid-atom alignment verification
├── load_dft_results.py        # Parse 3-results/*.dat files
├── sample_configs.py          # Sample potentials at atom positions
├── fit_pauli_beta1.py         # Phase 2: P_i fitting (β=1)
├── fit_charges.py             # Phase 3: Q_i fitting
├── compute_london.py          # Phase 4: London potential computation
├── fit_london.py              # Phase 4: L_i fitting
├── fit_pauli_beta.py          # Phase 5: β fitting
├── compare_pauli_methods.py   # Phase 6: Density vs getSurfMorse
├── data/
│   ├── NaCl_density.npy       # Converted density grid
│   ├── NaCl_potential.npy     # Converted potential grid
│   ├── NaCl_gridff_4ch.npy    # 4-channel GridFF array
│   └── metadata.json          # Grid metadata
├── results/
│   ├── phase2_pauli_beta1/    # P_i fit results
│   ├── phase3_charges/        # Q_i fit results
│   ├── phase4_london/         # L_i fit results
│   ├── phase5_beta/           # β fit results
│   └── phase6_comparison/     # Pauli method comparison
└── plots/
    ├── alignment/             # Alignment verification plots
    ├── phase2/                # P_i fit plots
    ├── phase3/                # Q_i fit plots
    ├── phase4/                # L_i fit plots
    ├── phase5/                # β fit plots
    └── phase6/                # Comparison plots
```

## References

- GridFF_FDBM_Fitting.md: Theoretical framework
- RigidSurfPotential_GridFF.md: GridFF infrastructure overview
- AFM_FDBM_fitting.md: AFM single-atom fitting methodology (β fitting)
- Surface.cl: getSurfMorse implementation
- Surface_utils.py: GridFF I/O and sampling utilities
- test_folded_fit_nacl1x1.py: Fitting example
- test_gridff_alignment.py: Alignment verification

---

# Session Report: DFT Cube File GridFF Generation and Plotting Refinement

## Summary

This session focused on refining the GridFF generation process to use DFT cube files (charge_density.cube and electrostatic_potential.cube) instead of GPU-generated Morse/Ewald potentials. The goal was to ensure correct unit conversions, axis ordering, and physically meaningful potentials for FDBM fitting.

## Key Changes

### 1. Cube File Loading (`tests/tMMFF/gen_gridff_nacl_gpu.py`)

**Switched from GPU to DFT cube files:**
- Pauli channel: loaded from `charge_density.cube` (converted e/Bohr³ → e/Å³)
- Coulomb channel: loaded from `electrostatic_potential.cube` (no Hartree→eV conversion - already in eV/V)
- London channel: set to zeros (not used in DFT-based approach)

**Unit conversions:**
- Grid spacing: Bohr → Angstrom (multiply by 0.529177)
- Density: e/Bohr³ → e/Å³ (divide by bohr²ang³)
- Potential: No conversion (cube values already in Volts/eV, NOT Hartree)

**V0 correction:**
- Subtract top-slice average: `V0 = phi[:,:,-1].mean()`
- Removes arbitrary constant from FFT Poisson solver
- Ensures potential decays to zero in vacuum

**Sign inversion:**
- DFT Hartree potential is for electrons (opposite sign to physical electrostatic)
- Invert sign: `phi = -phi`
- After correction: Na⁺ has positive potential, Cl⁻ has negative potential (physically correct)

**Grid origin shift:**
- Shift `g0[2]` by subtracting max substrate atom z
- Aligns z=0 with top substrate atom
- Enables intuitive plotting: z=0 is surface, z>0 is vacuum

### 2. Plotting Improvements (`pyBall/OCL/Surface_utils.py`)

**New parameters in `plot_alignment_summary`:**
- `plq_coeffs2`: Optional second set of (P, L, Q) coefficients for plotting both H and O total potentials in Z-profiles
- `vmin_vmax_xz`: Optional tuple (vmin, vmax) for manual XZ slice color scale limits

**XZ slice positions:**
- First slice: at center Y (iy_center)
- Second slice: at iy=32 (1/4 of cell size, ~2.82 Å) instead of iy=0

**Z-profile improvements:**
- Above Cl⁻: ix=0, iy=0 (Cl at origin)
- Above Na⁺: ix=32, iy=0 (Na at 2.82 Å along x)
- Both Z-profiles show H-Total (black) and O-Total (magenta) curves
- Enables visualization of differential interaction for different atom types

**Aspect ratio:**
- Changed from 'auto' to 'equal' for XZ plots
- Ensures 1 Å in z equals 1 Å in x

**Plot range:**
- XZ plots now display from z=0 (surface) up to 20 Å
- vmin/vmax still computed from z>2.0 Å only (vacuum region)
- Avoids huge substrate bulk values in color scale

**Color scale:**
- XY slice: symmetric range from single slice
- XZ slice: symmetric range from data above z=2.0 Å, or manual via vmin_vmax_xz

## Physics Validation

**Electrostatic potential after corrections:**
- Range: [-0.8, 330.5] eV (after V0 correction and sign inversion)
- At z=2 Å above surface: [-1.3, 0.8] eV (reasonable vacuum values)
- Na⁺: positive potential (repels H⁺)
- Cl⁻: negative potential (attracts H⁺)

**Total potential for H2O:**
- H (Q=+0.4238 e): repelled by Na⁺ (+0.32 eV), attracted to Cl⁻ (-0.54 eV)
- O (Q=-0.8476 e): attracted to Na⁺, repelled by Cl⁻
- Physics is now correct

## Files Modified

1. **`tests/tMMFF/gen_gridff_nacl_gpu.py`**
   - Added `_load_cube()` function for manual cube file parsing
   - Removed GPU Morse/Ewald generation path
   - Added V0 correction and sign inversion
   - Added grid origin shift for z=0 at surface
   - Updated to pass `plq_coeffs2` and `vmin_vmax_xz` to plotting

2. **`pyBall/OCL/Surface_utils.py`**
   - Added `plq_coeffs2` parameter to `plot_alignment_summary()`
   - Added `vmin_vmax_xz` parameter for manual XZ color scale
   - Changed second XZ slice from iy=0 to iy=32 (1/4 cell)
   - Changed Z-profiles to show above Cl⁻ and above Na⁺
   - Added both H-Total and O-Total curves in Z-profiles
   - Changed XZ aspect ratio to 'equal'
   - Changed XZ plot range to start from z=0

## Output Files

**GridFF:**
- `/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/Bspline_PLQd.npy`
- `/home/prokop/git/FireCore/tests/tSurf/small_mols_NaCl_New/Bspline_PLQd_meta.json`

**Diagnostic plots:**
- `/home/prokop/git/FireCore/tests/tMMFF/out_fdbm_dft_gridff/total_potential_H_atom.png`
- `/home/prokop/git/FireCore/tests/tMMFF/out_fdbm_dft_gridff/total_potential_O_atom.png`

## Workflow

To generate GridFF from DFT cube files and plot diagnostics:

```bash
python tests/tMMFF/gen_gridff_nacl_gpu.py
```

This script:
1. Loads charge_density.cube and electrostatic_potential.cube
2. Applies unit conversions, V0 correction, and sign inversion
3. Shifts grid origin so z=0 at top substrate atom
4. Generates GridFF array (nx, ny, nz, 4) with Pauli, London, Coulomb channels
5. Saves GridFF and metadata
6. Generates diagnostic plots with H and O PLQ coefficients

## Coefficient Estimation vs Fitting

**Current diagnostic coefficients:**
The P, L coefficients used in the diagnostic plots are **estimated from fixed atom parameters**, not fitted from DFT data:

```python
# Source data
- ElementTypes.dat: RvdW, EvdW per element
- H2O.xyz: RESP charges (O=-0.8476 e, H=+0.4238 e)

# PLQ conversion (alpha=1.8)
R = element_types['H'].RvdW
E = sqrt(element_types['H'].EvdW)  # Critical: sqrt convention
Q = RESP charge

e = exp(alpha * R)
cL = e * E
cP = e * cL = e^2 * E
```

**Resulting diagnostic coefficients:**
- H: P=7.877, L=0.587, Q=+0.4238
- O: P=27.778, L=1.190, Q=-0.8476

**Important:**
These coefficients are for **visualization and verification only**. They demonstrate that the GridFF data is physically reasonable (correct signs, reasonable magnitudes, expected H⁺/O⁻ behavior).

**For actual FDBM fitting:**
The P and L coefficients must be **refitted from DFT reference energies** using the linear least squares framework:
- Sample GridFF channels at atom positions for all DFT configurations
- Build feature matrix A with Pauli and London contributions per atom type
- Solve A * c = b where b = E_DFT - E_Coulomb (residual after subtracting fixed Coulomb)
- Fit P_t, L_t for each atom type t (H, O, etc.)

The diagnostic coefficients provide a reasonable starting point, but the fitted values will differ based on the actual DFT interaction energies.

## Next Steps

The GridFF generated from DFT cube files is now ready for FDBM fitting. The diagnostic plots show physically meaningful potentials with correct signs and reasonable magnitudes in the vacuum region. The next step is to:
1. Assemble DFT reference data (configurations + energies)
2. Sample GridFF channels at atom positions for all configurations
3. Build feature matrix and fit P, L coefficients using linear least squares

---

# DFT Reference Data Structure Analysis

## Directory Organization

### Geometries (`1-inputs/confs/`)
```
confs/{molecule}/ix{ix}_iy{iy}/orientxy_{site}/orientz_{angle}/z{distance}/input.xyz
```

**Molecules:** H2O-H, H2O-O, H2CO-O, HCONH2-H1, HCONH2-H2, HCONH2-O, HCOOH-H, HCOOH-O1, HCOOH-O2, NH3-H, NH3-N

**Sampling parameters:**
- `ix, iy`: Indices for irreducible wedge of NaCl substrate (triangular sampling)
- `orientxy`: Surface site orientation (Na, Cl, hollow)
- `orientz`: Tilt angle (-45, 0, 45 degrees)
- `z{distance}`: Height above surface (1.40 to 10.00 Å, 33 values)

**XYZ file format:**
- First 32 atoms: NaCl substrate (fixed, always the same)
- Last 3 atoms: H2O molecule (O, H, H)
- Example: `ix0_iy0/orientxy_Na/orientz_0/z02.00/input.xyz`

### DFT Energies (`3-results/`)
```
3-results/{molecule}/ix{ix}_iy{iy}-orientxy_{site}-orientz_{angle}.dat
```

**Format:** Two-column text file
- Column 1: z-distance (e.g., "01.40", "02.00")
- Column 2: DFT energy in eV (e.g., "-124.0881275172401")

**File sizes:**
- Full files: ~2550 bytes (33 lines × 3 repeats = 99 lines)
- Truncated files: ~850 bytes (e.g., ix6_iy0 files have fewer lines)

**Repetition pattern:**
- Most files: 33 z-distances repeated 3 times (99 lines total)
- Some files: Fewer lines (possibly incomplete calculations)

## Sampling Grid (from run.sh)

### Z-distances (33 values)
```
01.40, 01.50, 01.60, 01.70, 01.80, 01.85, 01.90, 01.95, 02.00,
02.05, 02.10, 02.15, 02.20, 02.25, 02.30, 02.35, 02.40, 02.45, 02.50,
02.60, 02.70, 02.80, 02.90, 03.00, 03.50, 04.00, 04.50, 05.00,
06.00, 07.00, 08.00, 09.00, 10.00
```

### Surface sites (3)
- `orientxy_Na`: Aligned over Na⁺ ion
- `orientxy_Cl`: Aligned over Cl⁻ ion
- `orientxy_hollow`: Aligned over hollow site

### Tilt angles (3)
- `orientz_-45`: -45 degrees
- `orientz_0`: 0 degrees (vertical)
- `orientz_45`: +45 degrees

### Irreducible wedge sampling (ix, iy)
**First half triangle:** ix=0..3, iy=0..ix
- Positions: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3)

**Second half triangle:** ix=4..6, iy=0..(6-ix)
- Positions: (4,0), (4,1), (4,2), (5,0), (5,1), (6,0)

**Total:** 16 unique (ix, iy) positions

### Total configurations per molecule
- 16 positions × 3 sites × 3 angles × 33 z-distances = 4752 configurations
- Some positions have fewer z-distances (e.g., ix6_iy0)

## Molecule Orientation (from molecules.dat)

**H2O-H:**
- Flag: `none`
- Atoms: `2 1 3` (i1=H, i2=O, i3=H)
- Meaning: H atom at origin, O aligned along x-direction, other H in -x/+y plane

**H2O-O:**
- Flag: `bisect`
- Atoms: `1 2 3 3` (i1=O, i21=H, i22=H, i3=H)
- Meaning: O at origin, H-bond direction is bisector of H-O-H angle

## Assembly Process (from assembly.f90)

1. **Read substrate (NaCl.xyz)** and molecule
2. **Translate molecule** so atom at origin is at (0,0,0)
3. **Calculate shift:**
   - `shift.x = ix * side.x / (2*ninter)` where ninter=6
   - `shift.y = iy * side.y / (2*ninter)`
   - `shift.z = zmax + dz` (zmax = top substrate atom z)
4. **Rotate molecule** to align with Na/Cl/hollow direction (xy-plane rotation)
5. **Apply tilt rotation** (orientz angle around y-axis)
6. **Write combined geometry** to input.xyz (substrate + molecule)

## Data Consolidation Strategy

### Proposed approach

1. **Parse all XYZ files:**
   - Extract only molecule coordinates (last 3 atoms for H2O)
   - Store as `(n_configs, n_atoms, 3)` numpy array
   - Store metadata: ix, iy, orientxy, orientz, zdist

2. **Parse all .dat files:**
   - Match filename pattern to geometry metadata
   - Extract z-distance and energy pairs
   - Handle repeated entries (take average or first occurrence)

3. **Consolidate into single structure:**
   ```python
   data = {
       'coords': (n_configs, 3, 3),  # H2O coordinates
       'energies': (n_configs,),      # DFT energies
       'metadata': {
           'ix': (n_configs,),
           'iy': (n_configs,),
           'orientxy': (n_configs,),  # 'Na', 'Cl', 'hollow'
           'orientz': (n_configs,),  # -45, 0, 45
           'zdist': (n_configs,)      # 1.40 to 10.00
       }
   }
   ```

4. **Generate XYZ movies:**
   - Fix z, scan orientation: slice by zdist
   - Fix orientation, scan z: slice by (ix, iy, orientxy, orientz)
   - Optional: include/exclude substrate atoms for debugging

### Key challenges

1. **File naming mismatch:**
   - Geometries: `ix0_iy0/orientxy_Na/orientz_0/z02.00/input.xyz`
   - Energies: `ix0_iy0-orientxy_Na-orientz_0.dat`
   - Need to map between directory structure and filename pattern

2. **Incomplete data:**
   - Some .dat files have fewer lines (e.g., ix6_iy0)
   - Need to handle missing configurations gracefully

3. **Repeated entries:**
   - Most .dat files have 3 repeats of the same z-scan
   - Need to decide whether to average or take first occurrence

4. **Substrate atoms:**
   - XYZ files always include 32 substrate atoms
   - For fitting, only need molecule coordinates
   - For debugging, may want substrate atoms in visualization

---

# DFT Data Consolidation Utilities

## Overview

Python utilities for organizing DFT reference data into consolidated formats suitable for FDBM fitting.

## Location

- **Utility functions:** `pyBall/atomicUtils.py` (end of file)
- **CLI script:** `tests/tSurf/small_mols_NaCl_New/consolidate_dft_data.py`

## Functions in atomicUtils.py

### `load_xyz_molecule_only(fname, n_substrate=32)`
Load XYZ file and extract only molecule atoms (last n atoms).

```python
enames, apos = load_xyz_molecule_only('ix0_iy0/orientxy_Na/orientz_0/z02.00/input.xyz', n_substrate=32)
```

### `load_dft_energies(dat_file, average_repeats=True)`
Load DFT energies from .dat file, averaging repeated entries.

```python
z_dists, energies = load_dft_energies('ix0_iy0-orientxy_Na-orientz_0.dat')
```

### `consolidate_xyz_to_npz(base_dir, molecule='H2O-H')`
Consolidate all XYZ files into single NPZ file with coordinates and metadata.

**Sorting:** Data is sorted by `(ix, iy, orientxy, orientz, zdist)` for consistent ordering.

```python
output = consolidate_xyz_to_npz('1-inputs', molecule='H2O-H')
# Output: 1-inputs/H2O-H_consolidated.npz
# Contents: coords, enames, ix, iy, orientxy, orientz, zdist
# Order: ix → iy → site (Cl/Na/hollow) → tilt (-45/0/45) → z (1.40→10.00)
```

### `match_xyz_with_energies(base_dir, results_dir, molecule='H2O-H')`
Match XYZ geometries with DFT energies from .dat files.

```python
matched = match_xyz_with_energies('1-inputs', '3-results', 'H2O-H')
# Returns list of dicts with: path, ix, iy, orientxy, orientz, zdist,
#                              enames, apos, energy, z_dists_all, energies_all
```

### `generate_xyz_movie(matched_data, output_file, ...)`
Generate XYZ movie from matched data with energies in comment line.

**Sorting:** Data is sorted by `zdist` (ascending) for sequential frame ordering.

```python
generate_xyz_movie(matched, 'movie_scan_z.xyz', include_substrate=False)
# Comment line: "E=-124.088128 eV ix=0 iy=0 site=Na tilt=0 z=2.00"
```

## CLI Script Usage

### Consolidate XYZ files
```bash
cd tests/tSurf/small_mols_NaCl_New
python consolidate_dft_data.py consolidate --molecule H2O-H
# Output: 1-inputs/H2O-H_consolidated.npz
```

### Match with DFT energies
```bash
python consolidate_dft_data.py match --molecule H2O-H
# Output: H2O-H_matched.npz
# Contains: coords, energies, enames, ix, iy, orientxy, orientz, zdist
```

### Generate XYZ movies

**Scan z distances (fixed site, tilt, ix, iy):**
```bash
python consolidate_dft_data.py movie --movie-type scan_z --site Na --tilt 0 --ix 0 --iy 0
# Output: H2O-H_scan_z_ix0_iy0_Na_tilt0.xyz
```

**With substrate atoms (for debugging):**
```bash
python consolidate_dft_data.py movie --movie-type scan_z --site Na --tilt 0 --ix 0 --iy 0 --include-substrate
```

**Scan orientations (fixed z, site, ix, iy):**
```bash
python consolidate_dft_data.py movie --movie-type scan_orientation --site Na --zdist 2.0 --ix 0 --iy 0
```

**Scan sites (fixed z, tilt, ix, iy):**
```bash
python consolidate_dft_data.py movie --movie-type scan_site --zdist 2.0 --tilt 0 --ix 0 --iy 0
```

**Scan XY positions (fixed z, site, tilt):**
```bash
python consolidate_dft_data.py movie --movie-type scan_xy --site Na --zdist 2.0 --tilt 0
```

### Print summary statistics
```bash
python consolidate_dft_data.py summary --molecule H2O-H
# Shows: total configs, breakdown by site/tilt/position, energy statistics
```

## Data Flow for FDBM Fitting

1. **Consolidate:** `1-inputs/confs/H2O-H/**/input.xyz` → `1-inputs/H2O-H_consolidated.npz`
2. **Match:** Combine with `3-results/H2O-H/*.dat` → `H2O-H_matched.npz`
3. **Fit:** Load NPZ, sample GridFF channels, build feature matrix, solve for P, L coefficients

## Output Format

### Consolidated NPZ (`{molecule}_consolidated.npz`)
- `coords`: (n_configs, n_atoms, 3) - molecule positions
- `enames`: (n_atoms,) - element symbols
- `ix`, `iy`: (n_configs,) - position indices
- `orientxy`: (n_configs,) - surface sites (Na/Cl/hollow)
- `orientz`: (n_configs,) - tilt angles (-45, 0, 45)
- `zdist`: (n_configs,) - z distances in Å

**Sorting order:** `(ix, iy, orientxy, orientz, zdist)` - ensures consistent indexing

### Matched NPZ (`{molecule}_matched.npz`)
Same as consolidated plus:
- `energies`: (n_configs,) - DFT energies in eV

**Sorting order:** Same as consolidated (inherited from matched data)

### XYZ Movies
Multi-frame XYZ file with energy in comment line:

**Sorting order:** Frames sorted by `zdist` (ascending) for sequential playback

```
3
E=-124.088128 eV ix=0 iy=0 site=Na tilt=0 z=2.00
O 0.6850813252 -0.0000000000 15.0057313252
H 0.0000000000 0.0000000000 14.3206500000
H 0.1677141234 -0.0000000000 15.8248796774
3
E=-124.110191 eV ix=0 iy=0 site=Na tilt=0 z=2.05
...
```

---

# DFT Fitting Report: Pauli Coefficients for H2O on NaCl

## What We Are Doing

We are fitting Pauli repulsion coefficients (P_H, P_O) and charges (q_H, q_O) for H2O molecules on a NaCl substrate using DFT reference data. The model is:

```
E_int = P_H * sum(rho@H) + P_O * sum(rho@O) + q_H * (sum(phi@H) - 2*phi@O)
```

where:
- `rho` = electron density from DFT `charge_density.cube` (Pauli channel)
- `phi` = electrostatic potential from DFT `electrostatic_potential.cube` (Coulomb channel)
- London/vdW channel is currently zero (not included in fit)
- Charges are constrained by neutrality: q_O = -2 * q_H

We use the consolidated DFT data from `H2O-H_matched.npz` and `H2O-O_matched.npz`, which contain:
- Molecular configurations (positions, orientations)
- DFT total energies
- Metadata: ix, iy, orientxy (Na/Cl/hollow), orientz (tilt angles), zdist

## Problems Encountered

### 1. Unphysical Negative Pauli Coefficients
Initial unconstrained least-squares fitting produced:
- P_H = -1.86 (negative - Pauli cannot be attractive)
- P_O = 43.05 (very large)
- q_H = 0.345, q_O = -0.69 (inflated vs RESP: 0.2, -0.4)

This violates the fundamental physical constraint that Pauli repulsion must be repulsive (P >= 0).

### 2. Jagged/Noisy Fit Curves
The initial 1D scan diagnostic plots showed extremely jagged, noisy fitted curves. Investigation revealed:
- We were using **nearest-neighbor grid sampling** (integer voxel indices)
- As the molecule moves continuously in z, the index snaps to the next voxel every ~0.083 Å
- When multiplied by large coefficients and summed over 3 atoms, these jumps look like "noise"

### 3. Model Cannot Reproduce Key Minima with Physical Constraints
Even after fixing the sampling method and applying physical constraints, the model fails to reproduce the two most important bonding configurations:
- **H2O-H over Cl** (hydrogen bond): DFT shows ~-0.25 eV attraction
- **H2O-O over Na** (electrostatic minimum): DFT shows ~-0.4 eV attraction

With constraints P_H >= 0, P_O >= 0 and fixed RESP charges, the optimizer chooses P_H = 0, P_O = 0 (pure Coulomb), which cannot match these minima.

### 4. Possible Z-Frame Misalignment
Visual inspection of the 1D scans suggests the DFT attraction may start at a different z than the model, indicating a possible coordinate frame offset between the DFT data and the GridFF sampling.

## What We Tried

### 1. Fixed Sampling Method
Replaced nearest-neighbor sampling with **trilinear interpolation**:
- Periodic wrap in x/y (slab geometry)
- Clamp in z (non-periodic)
- Smooth interpolation from 8 surrounding grid corners
- Result: fitted curves are now smooth, no more staircase artifacts

### 2. Restricted Fit to Key Profiles
Instead of fitting all 9504 configurations, we restricted to:
- Only tilt=0 configurations
- Only the two key bonding profiles:
  - H2O-H over Cl (ix=0, iy=0)
  - H2O-O over Na (ix=6, iy=0)
- Energy filter: E_int <= 0.5 eV
- Result: 177 configurations used for fit (down from 9504)

### 3. Exponential Weighting
Applied weighting to prioritize near-minimum energies:
- w = exp(-E_int / kT) with kT = 0.2 eV
- Normalized so max weight = 1
- Result: effective weight range [0.02, 1.0]

### 4. Physical Constraints
Enforced hard constraints:
- P_H >= 0 (Pauli cannot be attractive)
- P_O >= 0 (Pauli cannot be attractive)
- q_H fixed to RESP value (0.2)
- q_O = -2 * q_H (charge neutrality)

Implemented via grid scan over P_H (0 to 60, 601 points) and solving scalar LSQ for P_O, clamping P_O >= 0.

### 5. Z-Shift Scan
To test the alignment hypothesis, we scanned dz in [-1.0, +1.0] Å (81 steps):
- For each dz, we shift all atom positions by dz during GridFF sampling
- Choose dz that minimizes weighted SSE on the fit set
- Result: best dz = -0.100 Å (small shift, not the 0.5-1.0 Å suspected)

### 6. Model Baseline Subtraction
We corrected the fitting metric to baseline-subtract the model energy in the same way as DFT:
- For each group (mol, ix, iy, orientxy, orientz), compute E0 at max z
- Compare E_int = E - E0 for both DFT and model
- This ensures we are comparing interaction energies, not absolute energies

### 7. Diagnostic Plots
Generated two separate 1D scan diagnostic figures:
- `scan_DFT_vs_default.png`: Uses default atomic coefficients (P_H=7.877, P_O=27.778)
- `scan_DFT_vs_fitted.png`: Uses fitted coefficients (P_H=0, P_O=0 with current constraints)
- Each figure shows:
  - 4 panels: H2O-H/Na, H2O-H/Cl, H2O-O/Na, H2O-O/Cl
  - DFT curves (dotted)
  - Model curves (solid)
  - Residuals (model - DFT) on secondary y-axis
  - Caption with fitted parameters and fit settings

## What We Are Struggling With

### Core Issue: Linear Density Model is Inadequate
The fundamental problem is that the linear model `E = P * rho + q * phi` cannot simultaneously:
1. Reproduce the repulsive wall position (zmin)
2. Reproduce the well depth (Emin)
3. Maintain physical constraints (P >= 0, q fixed to RESP)

The optimizer's choice of P_H = 0, P_O = 0 is not a numerical artifact—it is the best solution under the constraints, meaning:
- Adding positive Pauli makes the fit worse
- The Coulomb term alone (with RESP charges) is insufficient to match the DFT minima
- The model is missing key physics (likely polarization/induction or charge transfer)

### Why Default Coefficients Look Better
Visually, the default atomic coefficients (P_H=7.877, P_O=27.778) produce curves that more closely follow the DFT reference than the fitted zeros. This suggests:
- The default coefficients were derived from atomic parameters (RvdW, EvdW) that capture some of the missing physics
- The linear density proxy is not a good representation of true Pauli repulsion
- We may need a nonlinear Pauli feature (e.g., rho^α with α > 1) or additional terms

### Possible Next Directions
1. **Nonlinear Pauli feature**: Replace linear density with rho^α (α ≈ 1.5-2.0) to make repulsion steeper at short range
2. **Polarization term**: Add a term proportional to electric field or gradient to capture induction
3. **Fit (zmin, Emin) directly**: Instead of SSE, match the minimum position and depth for the two key profiles
4. **Include London/vdW**: The London channel is currently zero; adding it may improve the fit significantly
5. **Relax charge constraints**: Allow q to deviate from RESP if it improves the fit, then check physical reasonableness

## Current Status
- GridFF generation from DFT cubes: working
- DFT data consolidation and matching: working
- Trilinear interpolation: implemented and verified
- Constrained fitting with P >= 0, q fixed: implemented
- Result: P_H = 0, P_O = 0 (model reduces to pure Coulomb)
- Conclusion: Linear density model is insufficient for this system with physical constraints

## Key Files
- `tests/tMMFF/gen_gridff_nacl_gpu.py`: Main fitting script
- `tests/tMMFF/consolidate_dft_data.py`: DFT data consolidation and matching
- `tests/tMMFF/out_fdbm_dft_gridff/scan_DFT_vs_default.png`: Default coefficients diagnostic
- `tests/tMMFF/out_fdbm_dft_gridff/scan_DFT_vs_fitted.png`: Fitted coefficients diagnostic
- `tests/tMMFF/out_fdbm_dft_gridff/fitted_PLQ_coeffs.json`: Fitted coefficients JSON

---

# Part 2: H-Bond Correction and Regularized Fitting (2026)

## Motivation

The initial linear density model (Pauli repulsion + Coulomb electrostatics) was insufficient to capture the hydrogen bonding physics of H2O on NaCl. The constrained fit with P ≥ 0 and fixed charges reduced to P_H = P_O = 0, meaning the model effectively became pure Coulomb. This indicated a missing physical component: a directional hydrogen bond correction.

## Channel 4 Evolution: From Kinetic Energy Density to Cl Repulsive Morse

### Attempt 1: Kinetic Energy Density (τ) via Laplacian of Gaussian

The first approach used the kinetic energy density as the 4th GridFF channel, approximated via the Laplacian of Gaussian (LoG) filtering of the electron density:

```
τ(r) ≈ LoG[ρ](r) = ∇²G_σ * ρ(r)
```

where G_σ is a Gaussian kernel with σ controlling the length scale.

**Implementation:**
- Computed τ via FFT-based convolution in `pyBall/OCL/Surface_utils.py:compute_kinetic_energy_density_log()`
- Added `--hbond_model tau` CLI option to `gen_gridff_nacl_gpu.py`
- Saved τ as a cube file for visualization

**Result:**
The τ channel did not provide a physically meaningful hydrogen bond correction. The kinetic energy density is a bulk property of the electron cloud and does not capture the directional, short-range repulsion characteristic of hydrogen bonds.

### Attempt 2: Cl Repulsive Morse Potential

The second approach replaced τ with a repulsive Morse potential computed only from substrate Cl atoms. This is physically motivated: hydrogen bonds involve H atoms approaching Cl (electronegative) sites, and the repulsive wall of the Morse potential provides the short-range correction needed.

**Morse Potential Form:**
```
V_Morse(r) = D_e * [exp(-α(r - R_eq)) - 1]^2
```

For the repulsive part (short-range), we use only the first term:
```
V_repulsive(r) ≈ D_e * exp(-2α(r - R_eq))
```

**Implementation:**
- Reused existing `GridFF_cl.make_MorseFF()` GPU kernel from `pyBall/OCL/GridFF.py`
- Added `compute_cl_repulsive_on_grid()` wrapper in `Surface_utils.py` to filter substrate atoms to Cl only
- Default parameters: R_eq = 1.80 Å, D_e = 0.0116 eV, α = 1.5 Å⁻¹
- Added `--hbond_model cl_morse` CLI option with parameter overrides

**Result:**
The Cl repulsive channel provided a physically meaningful correction that could be tuned via H_H and H_O coefficients in the model.

## Extended Model with H-Bond Correction

The final model includes 5 linear parameters:

```
E = P_H * |ρ_H|^β + P_O * |ρ_O|^β + H_H * τ_H + H_O * τ_O + q_H * φ_H + q_O * φ_O
```

where:
- `P_H, P_O`: Pauli repulsion coefficients (density to power β)
- `H_H, H_O`: H-bond correction coefficients (4th channel: Cl repulsive Morse)
- `q_H, q_O`: Coulomb charges (q_O = -2*q_H for charge neutrality)
- `β`: Density exponent (fixed at 1.15 based on manual fitting)
- `τ_H, τ_O`: 4th channel values sampled at H and O positions

## GUI for Manual Intuition Building

To build physical intuition and obtain reasonable starting parameters, we developed an interactive GUI in `tests/tMMFF/gui_fdbm_fit.py`.

**Features:**
- **Live 2×2 plotting**: DFT vs model for H2O-H/Cl and H2O-O/Na scan profiles
- **Spinbox controls**: Mouse-wheel support, expanded ranges (H_H, H_O: -10000..10000)
- **Per-atom decomposition**: Toggle to see H1, H2, O contributions separately
- **Geometry subplot**: Visualizes H2O orientation in xz plane
- **Save/Load parameters**: JSON export/import for reproducibility
- **Weighted RMS display**: Real-time fit quality feedback

**Manual Fit Result:**
Through interactive tuning, the user obtained physically reasonable parameters:
```json
{
  "P_H": 7.88,
  "P_O": 120.78,
  "H_H": -115.0,
  "H_O": 0.0,
  "q_H": 0.35,
  "beta": 1.15,
  "dz": -0.1
}
```

This manual fit provided a target for the automatic fitting procedure to match.

## Regularized Fitting: Solving the Ill-Conditioned Linear System

### The Problem: Ill-Conditioned Design Matrix

When we attempted an unconstrained linear least-squares fit for all 5 parameters (P_H, P_O, H_H, H_O, q_H), the results were unphysical:

```
P_H = -9.98   (negative repulsion - unphysical)
P_O = 553     (unrealistically large)
H_H = -66
H_O = 1290    (unrealistically large)
q_H = 1.07    (far from RESP reference 0.2)
```

The **condition number** of the design matrix was 1.17×10⁵, indicating severe ill-conditioning. This means the columns of the matrix are nearly linearly dependent: the Pauli density, Coulomb potential, and Cl repulsive channels are correlated through the underlying physics of the system.

When a matrix is ill-conditioned, small numerical errors or noise in the data can cause large, unphysical swings in the fitted parameters. The solver finds a mathematical solution that minimizes the residual, but that solution may violate physical constraints.

### Solution 1: Column Scaling

**Concept:**
When columns of a matrix have very different magnitudes (e.g., density values ~10⁴, Coulomb ~10⁰, tau ~10⁻³), the numerical linear algebra becomes unstable. Column scaling normalizes each column to have similar magnitude.

**Mathematical Form:**
```
A_scaled = A / col_scale
where col_scale[j] = sqrt(mean(A[:,j]²))
```

**Effect:**
- Condition number improved from 1.17×10⁵ to 7.7×10¹
- Solver becomes numerically stable
- Solution is equivalent in scaled coordinates, then transformed back

**Why it works:**
Scaling ensures that numerical precision (machine epsilon ~10⁻¹⁶) is not wasted on columns with large magnitudes. All columns contribute equally to the numerical stability.

### Solution 2: Ridge Regression with Priors (Tikhonov Regularization)

**Concept:**
Pure least-squares minimizes the residual error: `min ||A·c - b||²`. This allows parameters to take any value as long as it fits the data. Ridge regression adds a penalty term that keeps parameters close to reasonable "prior" values unless the data strongly demands otherwise.

**Mathematical Form:**
```
min ||A·c - b||² + λ||c - c_prior||²
```

where:
- `c_prior`: Prior values for parameters (physics-based guesses)
- `λ`: Regularization strength (how strongly to enforce the prior)
- The penalty term `λ||c - c_prior||²` discourages large deviations from the prior

**In our implementation:**
We use **parameter-specific regularization scales** to allow different parameters to deviate by different amounts:

```
λ_j = ridge_strength / prior_scale[j]²
min ||A·c - b||² + Σ_j λ_j (c_j - c_prior[j])²
```

**Prior values used:**
- `P_H, P_O`: Diagnostic values from vdW mapping (plq_H[0], plq_O[0])
- `H_H, H_O`: 0 (no prior knowledge, let data drive)
- `q_H`: RESP reference value (plq_H[2])

**Prior scales (characteristic deviations allowed):**
- `P_H`: 30 (can deviate by ~30 from diagnostic)
- `P_O`: 300 (larger uncertainty for O)
- `H_H, H_O`: 300 (allow large corrections if data demands)
- `q_H`: 0.5 (keep close to RESP)

**Didactic Explanation:**
Think of ridge regression as a "tension spring" connecting each parameter to its prior value. The spring constant is λ. If the data provides strong evidence for a different value, the spring stretches. If the data is weak or noisy, the spring pulls the parameter back toward the prior. This prevents overfitting to noise and ensures physically reasonable solutions.

**Combined with column scaling:**
The regularization is applied in scaled coordinates to ensure numerical stability:

```
c_scaled = c * col_scale
c_prior_scaled = c_prior * col_scale
λ_scaled = λ / col_scale²
min ||A_scaled·c_scaled - b||² + Σ_j λ_scaled[j] (c_scaled[j] - c_prior_scaled[j])²
```

### Results with Regularization

The regularized fit produced physically reasonable parameters close to the manual fit:

| Parameter | Manual Fit | Regularized Fit | Prior |
|-----------|------------|-----------------|-------|
| P_H       | 7.88       | 5.96            | 7.88  |
| P_O       | 120.78     | 104.18          | 27.78 |
| H_H       | -115.0     | -89.21          | 0     |
| H_O       | 0.0        | -5.22           | 0     |
| q_H       | 0.35       | 0.327           | 0.20  |

**Fit quality:**
- wRMS residual: 0.063 eV (higher than unconstrained fit, but physically meaningful)
- Condition number: 7.7×10¹ (well-conditioned)
- Parameters are in the same physical regime as manual fit

**Key insight:**
The regularization successfully prevented the solver from taking unphysical shortcuts (e.g., making P_H negative and compensating with huge H_O). Instead, it found a solution that balances data fit with physical plausibility.

## Fukui Function Integration and Ablation Test

### Background

The 4th channel (H-bond correction) was originally implemented using two approaches:
1. **tau (kinetic energy density)**: Computed via LoG-filtered Laplacian of electron density
2. **cl_morse (Cl repulsive Morse)**: GPU-computed Morse repulsion from Cl atoms

A third approach was identified: using the **Fukui f- function** as the 4th channel. The Fukui function represents the difference in electron density between the neutral and anionic states (f- = ρ(N) - ρ(N-1)), which provides a direct measure of electron affinity and reactivity at each point in space.

### Implementation

The Fukui function was integrated as a third option for the `--hbond_model` CLI argument in `gen_gridff_nacl_gpu.py`:

```python
parser.add_argument('--hbond_model', type=str, default='cl_morse', 
                    choices=['tau', 'cl_morse', 'fukui'], 
                    help='Hbond correction model: tau (kinetic energy density), 
                         cl_morse (Cl repulsive Morse), or fukui (Fukui f- function)')
```

**Loading logic:**
- Load `fukui_minus.cube` from DFT results directory
- Verify grid compatibility (shape, origin, spacing) with DFT density/potential cubes
- No resampling needed (grids are identical)
- Assign to `tau` variable for downstream compatibility

**Metadata updates:**
- Channel labels dynamically updated based on selected model
- Plot titles and y-axis labels reflect the chosen H-bond model
- Saved metadata includes `hbond_model` and `hbond_range`

### Fitted Parameters with Fukui

Using the Fukui f- function as the 4th channel, the regularized fit produced:

| Parameter | Value | Physical meaning |
|-----------|-------|-----------------|
| P_H       | 4.90  | Pauli coefficient for H atoms |
| P_O       | 101.56| Pauli coefficient for O atoms |
| H_H       | -34.15| H-bond coefficient for H atoms (Fukui channel) |
| H_O       | -26.42| H-bond coefficient for O atoms (Fukui channel) |
| q_H       | 0.289 | Charge for H atoms |
| beta      | 1.8   | Pauli exponent |

**Fit quality:**
- wRMS residual: 0.0613 eV
- Fukui f- range: [-0.27, 0.43] e/Ang³

### Ablation Test

To quantify the contribution of the Fukui 4th channel, an ablation test was implemented comparing:

1. **With Fukui**: Using fitted H_H and H_O coefficients
2. **Ablation (no 4th channel)**: Setting H_H = H_O = 0
3. **DFT reference**: Original DFT interaction energies

**Implementation:**
- New function `_plot_ablation()` generates 4-panel comparison plot
- Each panel shows DFT vs model with Fukui vs model without 4th channel
- Residuals plotted on twin y-axis for both cases
- RMS residuals computed for both cases

**Results:**
The ablation test (`scan_ablation_fukui.png`) clearly demonstrates:
- The Fukui 4th channel significantly improves fit quality
- RMS residuals are lower with Fukui vs ablation across all 4 panels
- H2O-H over Cl (hydrogen bonding) shows the largest improvement
- The ablation (no 4th channel) deviates significantly from DFT, especially in the bonding region

**Key insight:**
The Fukui function provides a physically meaningful descriptor for H-bond corrections by directly capturing electron affinity variations in space. This is superior to the Morse repulsion approach which only captures geometric repulsion from Cl atoms.

### Diagnostic Plots

**2D slices (XZ, YZ):**
- `total_potential_H_atom_fitted-1xz.png` / `total_potential_O_atom_fitted-1xz.png`
- Cross-sectional views showing spatial distribution of total potential
- Uses fitted coefficients with Fukui channel

**1D profiles:**
- `rho_fukui_profiles.png`: Z-profiles over Na/Cl atoms showing charge density and Fukui f- function
- `scan_DFT_vs_fitted.png`: 1D scan diagnostic (DFT vs fitted with Fukui)
- `scan_ablation_fukui.png`: Ablation test comparing with/without Fukui

### Summary of the Complete Workflow

1. **GridFF Generation**: Compute 4 channels (Pauli density, London=0, Coulomb potential, H-bond correction) with configurable model (tau/cl_morse/fukui)
2. **DFT Data Loading**: Load H2O-H and H2O-O scan data with interaction energies
3. **Sampling**: Pre-sample all channels at atom positions for all configurations
4. **Regularized Fitting**:
   - Scale design matrix columns for numerical stability
   - Apply ridge regression with physics-based priors
   - Solve linear system with `np.linalg.solve()`
5. **Validation**: Generate diagnostic plots comparing DFT vs fitted model
6. **Ablation Testing**: Quantify contribution of 4th channel by comparing with/without
7. **GUI Refinement**: Use interactive GUI for final manual tuning if needed

## Key Implementation Files

- `tests/tMMFF/gen_gridff_nacl_gpu.py`: Main fitting script with regularized solver
- `pyBall/OCL/Surface_utils.py`: 
  - `compute_kinetic_energy_density_log()`: τ computation via LoG
  - `compute_cl_repulsive_on_grid()`: Cl repulsive Morse wrapper
  - `sample_gridff_trilinear()`: Grid sampling with PBC
- `tests/tMMFF/gui_fdbm_fit.py`: Interactive GUI for manual fitting
- `pyBall/OCL/GridFF.py`: `GridFF_cl.makeMorseFF()`: GPU Morse kernel
- `cpp/common_resources/cl/GridFF.cl`: `makeMorseFF` OpenCL kernel

## References

- Tikhonov, A. N. (1963). "Solution of incorrectly posed problems and the regularization method". Soviet Mathematics Doklady.
- Ridge regression: Hoerl, A. E., & Kennard, R. W. (1970). "Ridge regression: Biased estimation for nonorthogonal problems". Technometrics.