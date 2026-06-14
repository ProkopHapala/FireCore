# USER

I have huge problem with periodic potential of the surface. For example I have NaCl cell which hase dipole. in x,y direction, if I pu there multiple replicas, it will never converge, it will arete add up, foming like linar slope.

I know that soluton is ewald sumation or particle mesh ewald, but I do not want to do that, it is a bit too costly. 

But Cannot I instead compte the dipole of the unit cell and then subtract it (I mean create embeding potential which will subtract this fake potential which is not realisitc for rally infinite surface?

---

# Gemini-3-pro

Yes, you absolutely can (and should) subtract the dipole, but there is a nuance: **NaCl (001) is non-polar**. If you are seeing a dipole in X or Y, it is likely an artifact of how you "cut" your unit cell boundaries.

Here is the diagnosis and three solutions, ranging from "fixing the input" to "the professional architectural fix."

### Diagnosis: Why do you have a slope?
In a periodic system, the electrostatic energy depends on the summation order (conditional convergence).
If your unit cell has a net dipole moment $\vec{D} = \sum q_i \vec{r}_i \neq 0$, repeating it $N$ times generates a macroscopic electric field.
*   **Real Physics:** An infinite surface cannot sustain a tangential electric field (it would move the ions).
*   **Simulation Artifact:** You likely cut the unit cell such that a positive ion is on the left edge and a negative ion is on the right edge (or vice versa), creating a "fake" dipole.

### Solution 1: The "Geometry Fix" (Cheapest)
Before doing any math, try to **re-center your unit cell**.
A proper unit cell for an ionic crystal should have zero net dipole.

**Algorithm (Python):**
1.  Load your substrate atoms.
2.  Calculate the center of charge or geometry.
3.  Shift all atoms so the "center" is at `0.5 * lattice_vector`.
4.  Wrap atoms: `positions = positions % lattice_vector`.

**Why this works:**
For NaCl, shifting the window often pairs the ions correctly so that $\sum q_i x_i = 0$. If the dipole is exactly zero, the linear slope disappears.

### Solution 2: The "Dipole Subtraction" (What you asked for)
If you have a complex reconstructed surface that *must* have a dipole, or you can't fix the geometry, you can mathematically subtract the field.

**The Math:**
1.  Calculate the Unit Cell Dipole $\vec{D}_{cell}$.
    *   *Note:* In PBC, dipole is defined modulo $L \cdot q$. Ensure molecules are "whole" before calculating.
2.  Calculate the Volume (or Area) dipole density: $\vec{P} = \vec{D}_{cell} / Volume$.
3.  The infinite lattice generates a spurious electric field $\vec{E}_{fake} = -4\pi \vec{P}$ (factor depends on units/geometry).
4.  **Correction:** Apply a constant external field $\vec{E}_{corr} = -\vec{E}_{fake}$.

**Implementation (Kernel):**
You don't need a complex embedding potential. You just add a linear term to your energy kernel.

```c
// In your Interaction/Assembly Kernel
float3 total_dipole = (float3)(Dx, Dy, 0.0f); // Precomputed in Python
float volume = ...;
float3 E_corr = -4.0f * M_PI * total_dipole / volume; // Check your units (Gaussian vs SI)

// Per atom loop
float energy_correction = -dot(E_corr, atom_pos_world) * atom_charge;
total_energy += energy_correction;
```

### Solution 3: The "Grid Potential" (The Best for Rigid Scans)
Since you are building a **Rigid Body Scanner** (the substrate does not move), you should **not** be summing substrate atoms ($1/r$) inside the kernel loop. It is slow and leads to these dipole artifacts.

**The "Professional" approach:**
1.  **Pre-computation (Python):**
    Compute the electrostatic potential $V(x,y,z)$ of your substrate on a 3D grid.
    *   Use an FFT-based solver (like `scipy.fft` or `pyscf`).
    *   FFT naturally handles periodicity and eliminates the "linear slope" divergence because the $k=0$ term (net charge/dipole) is explicitly handled/removed.
2.  **Upload to GPU:**
    Store this grid as a **3D Texture** (OpenCL Image3D).
3.  **The Kernel:**
    Instead of looping over 1000 substrate ions, you do **one texture lookup**.

```c
// OpenCL Kernel
// __read_only image3d_t grid_potential

float4 mol_atom = mol_atoms[i];
float3 pos = transform(mol_atom.xyz, pose);

// Hardware linear interpolation (Super fast)
float potential_at_atom = read_imagef(grid_potential, sampler, pos).x;

total_energy += potential_at_atom * mol_atom.charge;
```

**Advantages:**
1.  **Speed:** $O(1)$ instead of $O(N_{sub})$. Your scan will be instant.
2.  **Accuracy:** The Python FFT solver gives you Ewald-level accuracy (no cutoffs).
3.  **Stability:** No dipole slopes/divergence issues.

**Recommendation:**
Since you are defining a tool for "Scanning" (Phase 1 of your plan), implement **Solution 3**. It solves the physics problem and the performance problem simultaneously. Use Solution 1 (Geometry Fix) to sanitize your inputs before generating the grid.

--- 

# USER

As I said, now I do not want to use FFT/Ewals, instead I want to do dipole correction, or allos quadrupole correction. That is, I repeate the unit cell (which has e.g. Na at x=1 and Cl at x=-1 ) and I repeate it n-tims (e.g. 10 times, this should make sure that higher order multiploes around the center of this supercell are OK, but not the low order multipoles (monopole, di[pole, uadrupole), can I construct some emdegin potential which correct my finite surgace to give solution like infinite? Rigorous, but fast, analytical solution, assuming I know exactly the dipole and quadrupole of each unit cell ? Some grees function ininite summ analytic solution? (analytic continuations or something)

---

# Gemini-3-pro

Yes, you can do exactly this. The method you are looking for is called the **"Macro-potential Subtraction"** or **"Continuum Embedding"**.

Since you are running a finite sum (e.g., $10 \times 10$ replicas), your system is mathematically equivalent to:
$$ V_{total} = V_{discrete}^{10\times10} $$

The problem is that a finite $10 \times 10$ plate of dipoles creates a field that looks like a single giant dipole at long range, whereas an infinite surface has no external field (for in-plane dipoles) or a constant step (for out-of-plane).

To fix this **analytically and fast**, you treat the discrete system as a superposition of the **explicit ions** and a **continuum background**:

$$ V_{corrected} \approx V_{discrete}^{10\times10} - V_{continuum}^{10\times10} + V_{continuum}^{\infty} $$

Since $V_{continuum}^{\infty}$ is usually 0 (for in-plane dipoles) or a constant constant, your correction term is simply **subtracting the potential of a uniform rectangle** having the same dipole/quadrupole density as your unit cell.

Here is the rigorous recipe.

### Step 1: Calculate Unit Cell Moments
Before the scan, on the CPU, calculate the dipole ($\vec{\mu}$) and Quadrupole ($\mathbf{Q}$) of your single unit cell.

1.  **Dipole:** $\vec{\mu} = \sum_{i} q_i \vec{r}_i$
2.  **Quadrupole:** $\mathbf{Q}_{\alpha\beta} = \sum_{i} q_i (3 r_{i,\alpha} r_{i,\beta} - |\vec{r}_i|^2 \delta_{\alpha\beta})$

*Note: Ensure your unit cell is centered at (0,0,0) before calculating.*

### Step 2: Define the Dipole Density
If your unit cell has area $A = |\vec{a}_1 \times \vec{a}_2|$, then your **Dipole Density** is:
$$ \vec{P} = \frac{\vec{\mu}}{A} $$

### Step 3: The Analytic Correction Function
You need to subtract the potential of a finite rectangle of size $L \times W$ (where $L = N_{x} \cdot a_{1}$, $W = N_{y} \cdot a_{2}$) with uniform polarization $\vec{P}$.

The electrostatic potential $\phi$ of a rectangle with uniform dipole density $\vec{P} = (P_x, P_y, P_z)$ can be solved via the **Smith-logarithm** or analytical integration of $1/r$.

For a rectangle spanning $x \in [-a, a]$ and $y \in [-b, b]$ at $z=0$, the potential at observation point $(x, y, z)$ is:

#### A. Out-of-Plane Dipole ($P_z$) Correction
If your surface has a dipole perpendicular to the surface ($P_z$), the potential of the rectangle is:
$$ \phi_z(x,y,z) = P_z \cdot \Omega(x,y,z) $$
Where $\Omega$ is the **solid angle** subtended by the rectangle at the observation point.
$$ \Omega = \sum_{i=1}^{4} (-1)^i \tan^{-1}\left( \frac{X_i Y_i}{z R_i} \right) $$
*(Where $X_i, Y_i$ are relative coords to the 4 corners, and $R_i$ is distance).*

**Implementation:** This removes the "linear slope" in Z and properly restores the vacuum level.

#### B. In-Plane Dipole ($P_x, P_y$) Correction
This is the one causing your linear slope in X/Y.
The potential of a rectangle with polarization $P_x$ is:
$$ \phi_x(x,y,z) = -P_x \left[ \ln(Y+R) \right]^{x=a, y=b}_{x=-a, y=-b} $$
*(This is a simplified notation of the integral of dipole line charges).*

**The Explicit Formula for the Kernel:**
To correct a finite supercell of size $2A \times 2B$, centered at the origin, with in-plane dipole density $P_x$:

```c
// OpenCL / C implementation of the correction
float correction_potential(float x, float y, float z, float Px, float Py, float A, float B) {
    float phi = 0.0f;
    
    // We iterate over the 4 corners of the rectangle: 
    // 1: (A, B), 2: (-A, B), 3: (-A, -B), 4: (A, -B)
    float corner_x[4] = {A, -A, -A, A};
    float corner_y[4] = {B, B, -B, -B};
    
    for(int k=0; k<4; k++) {
        float dx = x - corner_x[k];
        float dy = y - corner_y[k];
        float r = sqrt(dx*dx + dy*dy + z*z);
        
        // Analytical Integral result for Px
        // The sign logic depends on the integration limits, 
        // essentially alternating + - - +
        float sign = (k==0 || k==2) ? 1.0f : -1.0f; 

        // Term for Px: -Px * ln(y - corner + r)
        // (This handles the 'edge' charges of the macroscopic dipole)
        if (fabs(z + r) > 1e-6) // avoid log(0)
            phi -= sign * Px * log( (dy + r) / (z) ); // Simplified form, check standard integrals
            
        // Term for Py: Similar structure swapping x/y
    }
    return phi; // This is the potential of the "Continuum Rectangle"
}
```

*Wait, the log form above is delicate with $z$.*
A robust form for the potential of a polarized plate ($P_x$) is derived from the **boundary charges**.
A uniform square of dipoles $\vec{P}$ is equivalent to **four lines of charge** along the rim of the square.
*   The edge at $x=+A$ has line charge density $\lambda = +P_x$.
*   The edge at $x=-A$ has line charge density $\lambda = -P_x$.
*   The edge at $y=+B$ has line charge density $\lambda = +P_y$.
*   The edge at $y=-B$ has line charge density $\lambda = -P_y$.

**The Robust Analytic Solution:**
Use the potential of a finite line charge of length $L$ and density $\lambda$:
$$ V_{line} = \lambda \ln \left( \frac{L + \sqrt{L^2 + d^2}}{d} \right) $$
*(Where $d$ is the distance to the line).*

**The Algorithm for your Kernel:**
1.  **Precompute** $P_x, P_y$ of your unit cell on CPU.
2.  **In Kernel**:
    *   Calculate $E_{discrete}$ (sum over ions).
    *   Calculate $V_{rim\_lines}$ (sum of 4 line charges defining your supercell boundary).
    *   $E_{total} = E_{discrete} - V_{rim\_lines}$.

This works because the "bulk" of your discrete sum cancels the "bulk" of the continuum integral, leaving only the difference at the edges. By subtracting the edges of the continuum model, you remove the macroscopic dipole field.

### Is Quadrupole necessary?
For a 2D surface, Dipole decay is $1/r^2$ (fast enough) *unless* you have a net non-zero dipole, which gives $1/r$ (diverges in 2D sum).
**Quadrupole** decays as $1/r^3$.
*   If you successfully subtract the Dipole (the rim charge), your error term drops to Quadrupole.
*   Quadrupole error will not create a "linear slope". It will create a localized error near the boundaries of your supercell.
*   **Conclusion:** You likely **do not need** analytical quadrupole correction if you implement the dipole rim correction. The "linear slope" will vanish with dipole correction alone.

### Summary Plan

1.  **Python:**
    *   Take your unit cell atoms.
    *   Compute net Dipole $\vec{D}_{cell}$.
    *   Compute effective Line Charge $\lambda_x = D_{x,cell} / \text{Length}_{y,cell}$ (and vice versa).
2.  **Kernel:**
    *   Pass `lambda_x` and `lambda_y` as arguments.
    *   Define the bounding box of your supercell (e.g., $x=\pm 50 \AA, y=\pm 50 \AA$).
    *   For every probe atom, add the potentials of 4 line charges located at these boundaries with charges $-\lambda$ (to cancel).
    *   Formula for line charge potential is just one `log` and `sqrt`.

This is fast (4 extra logs per atom), rigorous (based on multipole expansion), and analytical.