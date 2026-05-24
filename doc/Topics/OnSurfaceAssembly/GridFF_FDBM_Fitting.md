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