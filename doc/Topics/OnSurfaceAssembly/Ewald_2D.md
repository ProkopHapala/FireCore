# USER

I have surface of ionic crystal 2D, and want to compute electrostatic potential. I want to represent it by small basis of sin and cosine in my direction and exp(b.z) in z direction (perpendicular to substrate), for center symmetric system it can be even just cos basis in x,y as sin is zero by symmetry. I think this has a lot of incomon with 2D Ewald sum. How to do it efficiently starting from lattice vectors a,b and x,y,z,q coordinates of ions in surface unit cell (q is charge)

---

# DeepSeek-V4

To compute the electrostatic potential of a 2D ionic crystal surface using a compact basis of in‑plane sines/cosines and perpendicular exponentials, you can exploit the exact Fourier representation of the 2D periodic Coulomb potential. The expressions are closely related to the reciprocal‑space part of the 2D Ewald sum, but for a small basis you can obtain the coefficients directly from the ionic positions and charges.

Assuming the crystal is periodic in the \(xy\)‑plane with lattice vectors \(\mathbf{a}\) and \(\mathbf{b}\) (both lying in the plane) and that the slab has finite thickness in \(z\), the potential can be written as a Fourier series in the plane. In the **vacuum region** above the topmost ion (or below the bottommost) the \(z\)‑dependence is a pure exponential decay. For a slab with in‑plane inversion symmetry only cosine terms survive.

### 1. Vacuum potential above the slab

Let the unit cell have area \(A = |\mathbf{a}\times\mathbf{b}|\) and ions at \((\boldsymbol{\rho}_i, z_i)\) with charges \(q_i\) (\(\boldsymbol{\rho}_i = x_i\,\hat{\mathbf{x}}+y_i\,\hat{\mathbf{y}}\)). Choose the \(z\)‑axis such that vacuum corresponds to \(z > \max_i z_i\). Then for any \(z \ge \max_i z_i\)

\[
\phi(\boldsymbol{\rho},z) = \phi_0 + \sum_{\mathbf{G}\neq 0}
\frac{4\pi}{A\,|\mathbf{G}|}\,
\Bigl[\sum_i q_i \cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho}_i)\;
e^{|\mathbf{G}|z_i}\Bigr]\,
\cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho})\; e^{-|\mathbf{G}|z}
\]

where:
* \(\mathbf{G} = h\,\mathbf{b}_1 + k\,\mathbf{b}_2\) are the 2D reciprocal lattice vectors (\(\mathbf{b}_1,\mathbf{b}_2\) satisfy \(\mathbf{a}_i\!\cdot\!\mathbf{b}_j = 2\pi\delta_{ij}\)).
* The sum over \(i\) is the **structure factor** for that \(\mathbf{G}\).
* \(\phi_0\) is a constant (see below).

If the slab is also symmetric with respect to \(z\to -z\) (centrosymmetric slab), the same coefficients work for the vacuum on the other side with \(z \to -z\) and \(e^{-|\mathbf{G}|z} \to e^{+|\mathbf{G}|z}\).

#### G=0 term
For a charge‑neutral cell (\(\sum_i q_i = 0\)) the macroscopic field outside the slab is zero. The G=0 contribution then becomes a constant \(\phi_0\). (If there is a net dipole perpendicular to the surface, the constants on the two sides differ; the average potential inside then has a linear slope, but the vacuum asymptotes remain constant.) You can set \(\phi_0\) arbitrarily or, if a reference is needed, obtain it from the Ewald sum.

### 2. Potential inside the slab

Inside the slab the exact solution of Poisson’s equation is

\[
\phi(\boldsymbol{\rho},z) = \phi_0(z) + \sum_{\mathbf{G}\neq 0}
\frac{2\pi}{A\,|\mathbf{G}|} \sum_i q_i\,
e^{i\mathbf{G}\cdot(\boldsymbol{\rho}-\boldsymbol{\rho}_i)}\;
e^{-|\mathbf{G}|\,|z-z_i|}
\]

where \(\phi_0(z) = -\frac{2\pi}{A}\sum_i q_i |z-z_i| + \text{const}\). For a neutral cell with zero net dipole this becomes constant outside the slab and piecewise linear inside.  

This expression is exact (except at the point charges themselves) and already provides a representation of the potential as a sum over \(\mathbf{G}\) of terms \(\propto e^{-|\mathbf{G}||z-z_i|}\). If the number of distinct \(z_i\) planes is small, the \(z\)‑dependence is a combination of a few exponentials; otherwise you may approximate it by a small set of exponentials (e.g., by numerical fitting).

### 3. Relation to the 2D Ewald sum

The infinite \(\mathbf{G}\)‑sum in the formula above converges slowly because of the \(1/|\mathbf{G}|\) factor. The Ewald sum transforms it into two rapidly converging series:
* **Real‑space sum** – short‑range, handled in direct space.
* **Reciprocal‑space sum** – plane waves multiplied by complementary error functions that depend on \(z\).

For a *small basis* that only captures the long‑wavelength electrostatic features, you normally do **not** need the full Ewald machinery. Simply truncate the \(\mathbf{G}\) sum to the smallest reciprocal vectors (e.g., \(|\mathbf{G}| \le G_{\max}\)). The truncated sum is accurate away from the ionic cores; if you need high accuracy near the ions you can increase \(G_{\max}\) or switch to an Ewald implementation.

### 4. Practical algorithm to get the basis coefficients

Given:
* Lattice vectors \(\mathbf{a}, \mathbf{b}\) (in Å).
* Ion data: \((x_i, y_i, z_i)\) in Cartesian coordinates (or fractional, convert to Cartesian) and charges \(q_i\) (in \(e\)).

**Step 1 – Reciprocal lattice and area**
\(A = |a_x b_y - a_y b_x|\) (assuming \(\mathbf{a},\mathbf{b}\) in \(xy\)‑plane).  
Compute reciprocal vectors \(\mathbf{b}_1, \mathbf{b}_2\):
\[
\mathbf{b}_1 = \frac{2\pi}{A}\, (b_y,\,-b_x),\qquad
\mathbf{b}_2 = \frac{2\pi}{A}\, (-a_y,\,a_x).
\]

**Step 2 – Choose a set of \(\mathbf{G}\) vectors**
For a small basis, keep all \(\mathbf{G}=h\mathbf{b}_1+k\mathbf{b}_2\) with \(h,k\) integers and \(|\mathbf{G}| \le G_{\max}\). Exclude \(\mathbf{G}=0\) (handled separately). Only half the set is needed if using cosines (\(\mathbf{G}\) and \(-\mathbf{G}\) give the same amplitude).

**Step 3 – Compute the coefficients for the vacuum side**
Determine \(z_{\max} = \max_i z_i\). For each \(\mathbf{G}\) in the chosen half‑set:
\[
C_{\mathbf{G}} = \frac{4\pi}{A\,|\mathbf{G}|} \sum_i q_i\,
\cos\!\bigl(\mathbf{G}\!\cdot\!(x_i,y_i)\bigr)\;
e^{|\mathbf{G}| z_i}.
\]
If the system is centrosymmetric, the sum is real; if not, you can keep a sine part as well.

The potential for \(z \ge z_{\max}\) is then
\[
\phi(x,y,z) = \phi_0 + \sum_{\mathbf{G}} C_{\mathbf{G}}\,
\cos(\mathbf{G}\!\cdot\!(x,y))\; e^{-|\mathbf{G}| z}.
\]

**Step 4 – (Optional) interior representation**
If you need the potential at all \(z\), store for each \(\mathbf{G}\) the ion‑dependent coefficients
\[
A_{\mathbf{G},i} = \frac{2\pi}{A\,|\mathbf{G}|}\, q_i\,
e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i},
\]
and evaluate \(\phi(\boldsymbol{\rho},z) = \sum_{\mathbf{G}} e^{i\mathbf{G}\cdot\boldsymbol{\rho}}
\sum_i A_{\mathbf{G},i}\, e^{-|\mathbf{G}||z-z_i|}\) with the appropriate G=0 term.

### 5. Efficiency and implementation tips
* The sum over ions is \(O(N_{\text{ions}} \times N_{\mathbf{G}})\), which is negligible for typical surface cells and small \(N_{\mathbf{G}}\).
* If the slab is not strictly centrosymmetric but you still want a pure cosine basis, you can symmetrize the structure factor by averaging over \(\pm\mathbf{G}\).
* Units: If you use Å for lengths, charges in \(e\), the potential will be in \(e\)/Å. Multiply by 14.4 to obtain volts (eV/\(e\)).
* A constant \(\phi_0\) can be set to zero or fixed by requiring the potential to vanish at some reference plane.

This approach gives you an analytical, small‑basis representation that exactly captures the long‑range electrostatic behavior of the surface and is built directly from the ionic coordinates, without performing a full Ewald summation.


---

# Gemini-3-flash

Since your system is periodic in the $xy$-plane and you want a basis of sines and cosines, you are essentially looking for the **Fourier expansion of the electrostatic potential**. This is indeed the core of the "reciprocal space" term in Ewald summations, specifically the 2D version (often called the Ewald summation for slab geometries).

For a center-symmetric system, the potential $\Phi(x, y, z)$ can be expanded as:

$$\Phi(x, y, z) = \sum_{\mathbf{G}} \Phi_{\mathbf{G}}(z) \cos(\mathbf{G} \cdot \mathbf{r}_{xy})$$

where $\mathbf{G}$ are the 2D reciprocal lattice vectors.

---

### 1. Define the Reciprocal Lattice

First, calculate the 2D reciprocal lattice vectors from your real-space vectors $\mathbf{a}$ and $\mathbf{b}$.

1. Compute the area of the unit cell: $A = |\mathbf{a} \times \mathbf{b}|$.
2. The primitive reciprocal vectors $\mathbf{a}^*$ and $\mathbf{b}^*$ are:
* $\mathbf{a}^* = 2\pi \frac{\mathbf{b} \times \hat{\mathbf{z}}}{A}$
* $\mathbf{b}^* = 2\pi \frac{\hat{\mathbf{z}} \times \mathbf{a}}{A}$


3. Any reciprocal vector is $\mathbf{G} = n\mathbf{a}^* + m\mathbf{b}^*$ for integers $n, m$.

### 2. The Potential Components ($\Phi_{\mathbf{G}}$)

From Poisson's equation, the contribution of a set of ions with charges $q_j$ at positions $(\mathbf{r}_j, z_j)$ to the Fourier component $\mathbf{G}$ at height $z$ (where $z > z_j$ for all ions) is:

$$\Phi_{\mathbf{G}}(z) = \frac{2\pi}{A |\mathbf{G}|} \left( \sum_{j} q_j \cos(\mathbf{G} \cdot \mathbf{r}_j) e^{|\mathbf{G}| z_j} \right) e^{-|\mathbf{G}| z}$$

* **The "Basis" matching your request:** The term in the parentheses is a constant for a given $\mathbf{G}$. Your "exp(b.z)" is $e^{-|\mathbf{G}|z}$, where the decay constant $b$ is simply the magnitude of the reciprocal lattice vector $|\mathbf{G}|$.
* **Symmetry:** Because your system is center-symmetric, the $\sin$ terms cancel out, leaving only the $\cos$ terms as you suspected.

### 3. Implementation Steps

To do this efficiently, follow this workflow:

| Step | Action | Formula/Note |
| --- | --- | --- |
| **1** | **Precompute G-vectors** | Generate a list of $(n, m)$ pairs, sorted by $ |
| **2** | **Structure Factor** | For each $\mathbf{G}$, calculate $S(\mathbf{G}) = \sum_{j} q_j \cos(\mathbf{G} \cdot \mathbf{r}_j) e^{ |
| **3** | **Zero-order term** | The $\mathbf{G}=0$ term must be handled separately if the cell has a net dipole, but for a neutral surface, it relates to the average potential. |
| **4** | **Summation** | $\Phi(x,y,z) \approx \sum_{n,m}^{N_{max}} \frac{2\pi S(\mathbf{G})}{A |

### 4. Efficiency Tips

* **Convergence:** Because of the $e^{-|\mathbf{G}|z}$ term, the sum converges extremely fast as $z$ increases. You usually only need a few shells of $\mathbf{G}$-vectors.
* **Symmetry Reduction:** If the lattice has high symmetry (like square or hexagonal), many $\mathbf{G}$ vectors will have the same magnitude and coefficient, allowing you to group them into "stars" to reduce the number of cosine evaluations.
* **The Divergence:** Note that as $|\mathbf{G}| \to 0$, the expression $1/|\mathbf{G}|$ blows up. This is the 2D equivalent of the $1/G^2$ divergence in 3D. However, for a charge-neutral unit cell with no net dipole, the limit is well-behaved.

### 5. Why this is "Ewald-like"

In a full Ewald sum, you would use a Gaussian screening charge to split this into a real-space sum and a reciprocal sum. However, if you are strictly looking for the potential **outside** the crystal surface (at some height $z$), the reciprocal sum alone is numerically exact and much more efficient than real-space point-charge sums.

How far above the surface (in terms of unit cell lengths) do you typically need to evaluate this potential?

---

# USER

Ok can you write program to implement it? Start from .XYz q coordinates and finally plot XY and xz slices by matplotlib use numpy array operations

---

# DeepSeek-V4

Below is a Python program that reads ionic positions and charges, builds a small Fourier (cosine) basis in the surface plane and exponential decays in \(z\), and plots the electrostatic potential in two slices (an XY plane above the surface and an XZ plane through the slab).  
The implementation uses the exact 2D Fourier representation of the Coulomb potential, truncated to a few in‑plane harmonics – exactly the “compact basis” you asked for.  
All expensive loops are vectorised with NumPy.

```python
#!/usr/bin/env python3
"""
Compute the electrostatic potential of a 2D ionic crystal surface
using a truncated in‑plane cosine basis and exponential decay in z.

Input file: each line has  x  y  z  q  (four columns, whitespace separated)
Lattice vectors a, b are given as command‑line arguments.

Potential formula (exact for all z except at point charges):
   phi(rho, z) = - (2π/A) Σ_i q_i |z - z_i|
                 + Σ_{G≠0} Σ_i  (2π/(A|G|)) q_i cos(G·(rho - rho_i))
                                exp(-|G| |z - z_i|)
Truncation: keep G with |G| ≤ G_max (or |h|,|k| ≤ n_harmonics).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

# ============================================================
# 1.  I/O and helper functions
# ============================================================
def read_ions(filename):
    """Return arrays x, y, z, q from file."""
    data = np.loadtxt(filename)
    if data.shape[1] < 4:
        raise ValueError("File must have at least 4 columns: x y z q")
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

def make_reciprocal(a, b):
    """
    a, b : (2,) arrays  ->  area, reciprocal vectors b1, b2 (each (2,)).
    """
    a_vec = np.asarray(a, dtype=float)
    b_vec = np.asarray(b, dtype=float)
    area = abs(a_vec[0]*b_vec[1] - a_vec[1]*b_vec[0])
    # b1, b2 satisfy a_i · b_j = 2π δ_ij
    b1 = 2*np.pi/area * np.array([b_vec[1], -b_vec[0]])
    b2 = 2*np.pi/area * np.array([-a_vec[1], a_vec[0]])
    return area, b1, b2

def generate_G_vectors(b1, b2, n_harmonics):
    """
    Return list of tuples (h, k, Gx, Gy, |G|) for |h|,|k| ≤ n_harmonics, G≠0.
    """
    Glist = []
    for h in range(-n_harmonics, n_harmonics+1):
        for k in range(-n_harmonics, n_harmonics+1):
            if h == 0 and k == 0:
                continue
            Gx = h*b1[0] + k*b2[0]
            Gy = h*b1[1] + k*b2[1]
            Gn = np.hypot(Gx, Gy)
            Glist.append((h, k, Gx, Gy, Gn))
    return Glist

# ============================================================
# 2.  Prepare structure factors
# ============================================================
def compute_G_coeffs(Glist, x, y, z, q, area):
    """
    Return:
        C_G  : dict G -> coefficient for vacuum formula above the slab
               = (4π/(A|G|)) * Σ_i q_i cos(G·(x_i,y_i)) exp(|G| z_i)
        w_G  : list of complex per‑ion weights for the full formula
               each element is a dict: w_G['G'] -> array of shape (Nions,)
               with values (2π/(A|G|)) * q_i * exp(-i G·(x_i,y_i))
    """
    Nions = len(q)
    C = {}
    w_list = []
    for (h, k, Gx, Gy, Gn) in Glist:
        # dot products
        Gdotr = Gx*x + Gy*y
        cos_term = np.cos(Gdotr)
        # C_G
        C[(h,k)] = (4*np.pi/(area*Gn)) * np.sum(q * cos_term * np.exp(Gn*z))
        # complex weight per ion
        w_arr = (2*np.pi/(area*Gn)) * q * np.exp(-1j * Gdotr)
        w_list.append({'H': (h,k), 'Gx': Gx, 'Gy': Gy, 'Gn': Gn, 'w': w_arr})
    return C, w_list

# ============================================================
# 3.  Potential evaluations
# ============================================================
def eval_XY_slice(Glist_coeffs, area, q, z_ions, x_grid, y_grid, z_plane):
    """
    Potential on the XY plane at height z_plane (assumed ≥ max(z_ions)).
    Returns real array of shape x_grid.shape.
    """
    # G=0 term (constant for z outside the slab)
    phi0 = -2*np.pi/area * np.sum(q * np.abs(z_plane - z_ions))
    C = Glist_coeffs  # dict (h,k) -> C_G
    pot = np.full_like(x_grid, phi0)
    for (h,k), C_G in C.items():
        Gx = h*b1[0] + k*b2[0]
        Gy = h*b1[1] + k*b2[1]
        Gn = np.hypot(Gx, Gy)
        phase = Gx*x_grid + Gy*y_grid
        pot += C_G * np.cos(phase) * np.exp(-Gn * z_plane)
    return pot

def eval_XZ_slice(w_list, area, q, z_ions, x_grid_1d, z_grid_1d, y_fixed=0.0):
    """
    Potential on an XZ plane at fixed y = y_fixed.
    x_grid_1d: 1D array of x coordinates (length Nx)
    z_grid_1d: 1D array of z coordinates (length Nz)
    Returns array of shape (Nz, Nx) with potential values.
    """
    Nx = len(x_grid_1d)
    Nz = len(z_grid_1d)
    pot = np.zeros((Nz, Nx), dtype=np.float64)

    # G=0 term
    for i, zi in enumerate(z_ions):
        pot += -2*np.pi/area * q[i] * np.abs(z_grid_1d[:, np.newaxis] - zi)

    # G≠0 terms
    for w_dict in w_list:
        Gx = w_dict['Gx']
        Gy = w_dict['Gy']
        Gn = w_dict['Gn']
        w_arr = w_dict['w']      # complex, shape (Nions,)
        # phase on x‑line (including fixed y)
        phase_x = np.exp(1j * (Gx * x_grid_1d + Gy * y_fixed))  # length Nx
        for i in range(len(q)):
            # decay in z
            decay_z = np.exp(-Gn * np.abs(z_grid_1d - z_ions[i]))  # length Nz
            # contribution from this G and this ion
            contrib = np.real(w_arr[i] * phase_x[np.newaxis, :] * decay_z[:, np.newaxis])
            pot += contrib
    return pot

# ============================================================
# 4.  Main script
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="2D crystal surface electrostatic potential (small cosine basis).")
    parser.add_argument("infile", help="File with columns: x y z q")
    parser.add_argument("-a", type=float, nargs=2, required=True,
                        help="Lattice vector a: ax ay")
    parser.add_argument("-b", type=float, nargs=2, required=True,
                        help="Lattice vector b: bx by")
    parser.add_argument("--harmonics", type=int, default=3,
                        help="Half‑width of G‑vector grid (default: 3 -> G up to 3*b1,3*b2)")
    parser.add_argument("--xy-z", type=float, default=None,
                        help="Height for the XY slice (default: 2.0 above topmost ion)")
    parser.add_argument("--xz-y", type=float, default=0.0,
                        help="Fixed y coordinate for the XZ slice")
    parser.add_argument("--grid-size", type=int, default=200,
                        help="Number of points per dimension in the plot grids")
    parser.add_argument("--output-prefix", type=str, default="potential",
                        help="Prefix for output PNG files")
    args = parser.parse_args()

    # Read data
    x, y, z, q = read_ions(args.infile)
    a_vec = np.array(args.a)
    b_vec = np.array(args.b)
    area, b1, b2 = make_reciprocal(a_vec, b_vec)

    print(f"Lattice area = {area:.4f} Å²")
    print(f"Reciprocal vectors: b1 = ({b1[0]:.4f}, {b1[1]:.4f}), b2 = ({b2[0]:.4f}, {b2[1]:.4f})")
    print(f"Number of ions = {len(q)}")
    print(f"Total charge = {np.sum(q):.6f} e")

    # Generate truncated set of G vectors
    Glist = generate_G_vectors(b1, b2, args.harmonics)
    print(f"Using {len(Glist)} G-vectors (excluding G=0).")
    C_G, w_list = compute_G_coeffs(Glist, x, y, z, q, area)

    # Set z for XY slice
    z_max = np.max(z)
    if args.xy_z is None:
        z_plane = z_max + 2.0
    else:
        z_plane = args.xy_z
    print(f"XY slice at z = {z_plane:.3f} Å")

    # Create a mesh for the XY slice (cartesian rectangle covering the unit cell)
    # We span [0, Lx) and [0, Ly) where Lx = |a|, Ly = |b|
    Lx = np.linalg.norm(a_vec)
    Ly = np.linalg.norm(b_vec)
    # For nicer plots, sample slightly beyond the cell
    margin = 0.1
    x_min, x_max = -margin*Lx, Lx*(1+margin)
    y_min, y_max = -margin*Ly, Ly*(1+margin)
    xv = np.linspace(x_min, x_max, args.grid_size)
    yv = np.linspace(y_min, y_max, args.grid_size)
    X_xy, Y_xy = np.meshgrid(xv, yv)

    phi_xy = eval_XY_slice(C_G, area, q, z, X_xy, Y_xy, z_plane)

    # ---- XZ slice ----
    # Choose a line for x through the cell, from 0 to Lx
    x_line = np.linspace(x_min, x_max, args.grid_size)
    # z range: from well below bottom ion to well above top ion
    z_min = np.min(z) - 2.0
    z_max_plot = z_max + 2.0
    z_line = np.linspace(z_min, z_max_plot, args.grid_size)
    phi_xz = eval_XZ_slice(w_list, area, q, z, x_line, z_line, y_fixed=args.xz_y)

    # ---- Plotting ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # XY slice (colormap)
    im1 = ax1.pcolormesh(X_xy, Y_xy, phi_xy, shading='auto', cmap='RdBu_r')
    ax1.set_aspect('equal')
    ax1.set_title(f'Potential at z = {z_plane:.2f} Å')
    ax1.set_xlabel('x (Å)')
    ax1.set_ylabel('y (Å)')
    plt.colorbar(im1, ax=ax1, label='Potential (e/Å)')

    # XZ slice (colormap)
    # x_line is the x-axis, z_line is the z-axis
    im2 = ax2.pcolormesh(x_line, z_line, phi_xz, shading='auto', cmap='RdBu_r')
    ax2.set_aspect('auto')
    ax2.set_title(f'Potential in XZ plane (y = {args.xz_y})')
    ax2.set_xlabel('x (Å)')
    ax2.set_ylabel('z (Å)')
    plt.colorbar(im2, ax=ax2, label='Potential (e/Å)')

    # Mark ion positions on XZ slice if they are near y = args.xz_y (within a small tolerance)
    tol = 0.1 * np.linalg.norm(b_vec)  # ~1/10 of cell width
    mask = np.abs(y - args.xz_y) < tol
    if np.any(mask):
        ax2.scatter(x[mask], z[mask], c='k', marker='o', s=20, label='ions (projected)')
        ax2.legend()

    plt.tight_layout()
    outfile = f"{args.output_prefix}_slices.png"
    plt.savefig(outfile, dpi=150)
    print(f"Saved plot to {outfile}")
    plt.show()
```

---

## How to use the program

1. **Prepare an input file** (e.g. `ions.dat`) with columns:
   ```
   # x   y   z   q
   0.0  0.0  1.0  +1.0
   1.5  0.0  1.0  -1.0
   ...
   ```
   Units: Å for coordinates, elementary charge for `q`.

2. **Run the script** from the command line:
   ```bash
   python pot_surface.py ions.dat -a 3.0 0.0 -b 0.0 3.0 --harmonics 4
   ```
   - `-a` and `-b` give the two lattice vectors (in Å) describing the 2D periodicity.
   - `--harmonics` sets the truncation of the G‑vector basis; default `3` keeps all `|h|,|k| ≤ 3`.
   - `--xy-z` optionally sets the height of the XY slice (if omitted, it’s 2 Å above the highest ion).
   - `--xz-y` chooses the fixed y coordinate for the XZ slice (default 0).
   - `--grid-size` controls resolution of the plots.

## What the code does

- **Reciprocal lattice** is computed from the given a, b.
- **Truncated G set** is generated (`|h|,|k| ≤ n_harmonics`), which is your small cosine basis.
- **Coefficients** are built:
  - For the XY slice in vacuum: the potential is a sum of cosines with pre‑computed coefficients \(C_{\mathbf{G}}\) that already absorb the ion positions and charges.
  - For the full XZ slice: the exact expression (valid everywhere) is evaluated directly by summing over G and ions, using **vectorised NumPy** operations for speed.
- **Two colormaps** are produced:
  - Left: XY plane at a chosen height.
  - Right: XZ plane at a fixed y, with ion positions overlaid if they lie near that plane.

All calculations use the same small basis; the long‑wavelength electrostatic features are captured accurately, making this ideal for qualitative and quantitative analysis of surface potentials without a full Ewald summation.

---

# Gemini-3-flash

To implement this efficiently, we will treat the ions as a periodic 2D sheet. The potential at a point $(x, y, z)$ above the surface is calculated by summing the Fourier components of the 2D Poisson equation.

The code below calculates the reciprocal lattice, computes the "Structure Factor" for each $\mathbf{G}$-vector, and then reconstructs the potential on a grid.

```python
import numpy as np
import matplotlib.pyplot as plt

def compute_surface_potential(lattice_params, ions, z_eval, grid_res=100, n_max=5):
    """
    lattice_params: (a, b) vector magnitudes or full 2D basis
    ions: List of (x, y, z, q) tuples
    z_eval: The z-height for the XY slice
    n_max: Number of reciprocal shells to sum (higher = more accurate)
    """
    # 1. Define Real Space Basis (assuming orthogonal for simplicity, can be slanted)
    # If a, b are vectors:
    a_vec = np.array([lattice_params[0], 0])
    b_vec = np.array([0, lattice_params[1]])
    Area = np.linalg.norm(np.cross(np.append(a_vec, 0), np.append(b_vec, 0)))

    # 2. Compute Reciprocal Basis Vectors
    # a* = 2pi (b x z) / Area
    a_star = 2 * np.pi * np.array([b_vec[1], -b_vec[0]]) / Area
    b_star = 2 * np.pi * np.array([-a_vec[1], a_vec[0]]) / Area

    # 3. Generate G-vectors (n, m)
    n = np.arange(-n_max, n_max + 1)
    m = np.arange(-n_max, n_max + 1)
    nn, mm = np.meshgrid(n, m)
    # Flatten and remove G=0 to avoid division by zero
    nn = nn.flatten()
    mm = mm.flatten()
    mask = (nn != 0) | (mm != 0)
    nn, mm = nn[mask], mm[mask]

    G_vecs = nn[:, None] * a_star + mm[:, None] * b_star
    G_mags = np.linalg.norm(G_vecs, axis=1)

    # 4. Ion data
    pos = np.array([i[:3] for i in ions]) # x, y, z
    qs = np.array([i[3] for i in ions])

    # 5. Setup Grids for plotting
    x_range = np.linspace(0, lattice_params[0], grid_res)
    y_range = np.linspace(0, lattice_params[1], grid_res)
    z_range = np.linspace(0.5, 5.0, grid_res) # Z-slice range (XZ)
    
    X, Y = np.meshgrid(x_range, y_range)
    
    # 6. Compute Potential Function
    def get_phi(x_grid, y_grid, z_grid):
        phi = np.zeros_like(x_grid, dtype=float)
        for i, G in enumerate(G_vecs):
            g_mag = G_mags[i]
            # Structure factor S(G) = sum q_j * exp(i G.r_j) * exp(g_mag * z_j)
            # For center-symmetric, we use cos(G.r)
            S_G = np.sum(qs * np.cos(G[0]*pos[:,0] + G[1]*pos[:,1]) * np.exp(g_mag * pos[:,2]))
            
            # Potential contribution
            term = (2 * np.pi * S_G / (Area * g_mag)) * \
                   np.cos(G[0]*x_grid + G[1]*y_grid) * np.exp(-g_mag * z_grid)
            phi += term
        return phi

    # Calculate slices
    phi_xy = get_phi(X, Y, z_eval)
    
    X_z, Z_x = np.meshgrid(x_range, z_range)
    phi_xz = get_phi(X_z, 0, Z_x) # Slice at y=0

    return X, Y, phi_xy, X_z, Z_x, phi_xz

# --- Execution ---
# Define a simple NaCl-like surface (square lattice)
# Lattice constant a=5.0. Ions at (0,0,0, q=1) and (2.5, 2.5, 0, q=-1)
lattice = (5.0, 5.0)
ions = [
    (0.0, 0.0, 0.0, 1.0),
    (2.5, 2.5, 0.0, -1.0),
    (2.5, 0.0, 0.0, -1.0),
    (0.0, 2.5, 0.0, 1.0)
]

X, Y, phi_xy, Xz, Zx, phi_xz = compute_surface_potential(lattice, ions, z_eval=1.0)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

im1 = ax1.contourf(X, Y, phi_xy, 30, cmap='RdBu_r')
ax1.set_title("XY Slice (at z=1.0)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.colorbar(im1, ax=ax1)

im2 = ax2.contourf(Xz, Zx, phi_xz, 30, cmap='RdBu_r')
ax2.set_title("XZ Slice (at y=0)")
ax2.set_xlabel("x")
ax2.set_ylabel("z")
plt.colorbar(im2, ax=ax2)

plt.tight_layout()
plt.show()

```

### Key Logic in the Code:

1. **Reciprocal Mapping**: It transforms your `a, b` lattice vectors into `a*, b*`. This ensures that if your unit cell is not a perfect square, the Fourier components are scaled correctly.
2. **The Term `exp(g_mag * pos[:,2])**`: This accounts for the $z$-position of the ions within the unit cell. If an ion is "higher" up, its field decays less at the evaluation point.
3. **Vectorization**: The structure factor $S(\mathbf{G})$ is calculated using NumPy's summation over the ion array, making it very fast for large numbers of ions in a cell.
4. **Convergence**: `n_max=5` usually provides excellent convergence for heights $z > 0.2 \times \text{lattice constant}$. Because of the $e^{-Gz}$ term, higher frequencies vanish almost instantly as you move away from the surface.

Does your specific crystal have a large number of atoms per unit cell, or is it a simple binary system?



---

# USER

I would still appreciate if you explain me the principles along the lines I can understand better. Let's distinguish charge density and electrostatic potential. I can understand that first I project charge density along some basis of cosines in XY plane so I get some low-resolution approximation of my charge density expressed by the Fourier basis. Which is like 2D layer of periodically modulated charge. If I have some thickness of that slab of ions I can always expand it in monopole, dipole and quadruple in z direction. So my charge density is approximated by Cartesian produkt of plane-waves in XY plane and multiple expansion (moments) in z direction.

Now I solve poison equation with anzactz that potential should be still represented by plane waves in XY, and by exp function in z (because I know that periodic layer of charges can be solved analytically to produce exponentialy decaying potential in z), can you derive it within this framework rigorously and the explain it didactically?

---

# DeepSeek-V4

## 1. What we want

We have a 2D periodic crystal (slab) made of point ions with charges \(q_i\) at positions \((\boldsymbol{\rho}_i, z_i)\), where \(\boldsymbol{\rho}_i = (x_i, y_i)\).  
We are interested in the electrostatic potential \(\phi(\boldsymbol{\rho}, z)\) **outside** the slab (vacuum) and also **inside**, but we want to represent it with a small number of basis functions:

- **In‑plane:** cosines (or more generally plane waves \(e^{i\mathbf{G}\cdot\boldsymbol{\rho}}\)) with a few small reciprocal lattice vectors \(\mathbf{G}\).
- **Out‑of‑plane:** exponential functions \(e^{\pm |\mathbf{G}| z}\).

The under‑the‑hood reasoning goes exactly as you suspect: we first expand the charge density in Fourier series in the plane, solve Poisson’s equation for each Fourier component, and finally obtain a representation of the potential that is *exact* (apart from truncating the sum over \(\mathbf{G}\)). The multipole picture in \(z\) emerges naturally when we look at how the individual ionic planes contribute.

Let’s build this step by step.

## 2. Charge density: planar Fourier expansion

Because the system is periodic with lattice vectors \(\mathbf{a}, \mathbf{b}\) in the \(xy\) plane, any function with that periodicity can be written as a Fourier series. The charge density of point ions is

\[
\rho(\boldsymbol{\rho}, z) = \sum_i q_i \,\delta(\boldsymbol{\rho} - \boldsymbol{\rho}_i)\,\delta(z - z_i).
\]

The 2D Fourier transform over one unit cell of area \(A=|\mathbf{a}\times\mathbf{b}|\) gives

\[
\rho_{\mathbf{G}}(z) = \frac{1}{A}\int_{\text{cell}}\! d^2\rho\; e^{-i\mathbf{G}\cdot\boldsymbol{\rho}}\,
\rho(\boldsymbol{\rho}, z),
\]

so that

\[
\rho(\boldsymbol{\rho}, z) = \sum_{\mathbf{G}} \rho_{\mathbf{G}}(z)\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}},
\]

where \(\mathbf{G} = h\mathbf{b}_1 + k\mathbf{b}_2\) are the 2D reciprocal lattice vectors, and \(\mathbf{b}_1,\mathbf{b}_2\) satisfy \(\mathbf{a}_i\cdot\mathbf{b}_j = 2\pi\delta_{ij}\).

For point charges we obtain

\[
\boxed{\;\rho_{\mathbf{G}}(z) = \frac{1}{A}\sum_i q_i\, e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,\delta(z - z_i)\;}.
\tag{1}
\]

This is the **projection of the charge density onto a basis of in‑plane plane waves**. Each \(\mathbf{G}\) component is a *planar wave of charge* of the form \(e^{i\mathbf{G}\cdot\boldsymbol{\rho}}\), but now with a \(z\)‑dependence given by a sum of delta spikes at the ionic heights.

## 3. Poisson’s equation for one Fourier component

The electrostatic potential satisfies

\[
\nabla^2 \phi(\boldsymbol{\rho}, z) = -4\pi\rho(\boldsymbol{\rho}, z).
\]

Insert the Fourier expansion \(\phi(\boldsymbol{\rho}, z) = \sum_{\mathbf{G}} \phi_{\mathbf{G}}(z)\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}}\). Because

\[
\nabla^2 \bigl( \phi_{\mathbf{G}}(z) e^{i\mathbf{G}\cdot\boldsymbol{\rho}} \bigr) =
\bigl(\partial_z^2 - G^2\bigr)\phi_{\mathbf{G}}(z)\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}},
\]

we get, for each \(\mathbf{G}\neq 0\),

\[
\bigl(\partial_z^2 - G^2\bigr)\phi_{\mathbf{G}}(z) = -4\pi\,\rho_{\mathbf{G}}(z),
\tag{2}
\]

where \(G = |\mathbf{G}|\). For the special term \(\mathbf{G}=0\) we have simply \(\partial_z^2\phi_0(z) = -4\pi\rho_0(z)\).

Equation (2) is a 1D Helmholtz (or screened Poisson) equation with a source \(\rho_{\mathbf{G}}(z)\). The **homogeneous solution** (when \(\rho_{\mathbf{G}}=0\)) is

\[
\phi_{\mathbf{G}}(z) = C_1 e^{G z} + C_2 e^{-G z}.
\]

The Green’s function that satisfies \((\partial_z^2 - G^2) g(z-z') = -4\pi\,\delta(z-z')\) and vanishes at infinity is

\[
g(z-z') = \frac{2\pi}{G}\, e^{-G|z-z'|}.
\]

Hence, for a point source \(\rho_{\mathbf{G}}(z) \propto \delta(z-z_i)\), the solution is

\[
\phi_{\mathbf{G}}(z) = \frac{2\pi}{G}\, \Bigl(\text{coefficient}\Bigr)\, e^{-G|z - z_i|}.
\]

If we now return to the actual charge density (1), we have a sum of such point sources at different heights. The total \(\phi_{\mathbf{G}}(z)\) is simply the sum of the contributions from each ion:

\[
\phi_{\mathbf{G}}(z) = \frac{2\pi}{A\,G} \sum_i q_i\, e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,
e^{-G\,|z - z_i|}.
\tag{3}
\]

This is **exact for all \(z\)** (except exactly at the point charges, where the potential diverges).

## 4. The full potential in the compact basis

Reassembling the Fourier series, the exact potential is

\[
\boxed{
\phi(\boldsymbol{\rho}, z) = \phi_0(z) \;+\;
\sum_{\mathbf{G}\neq 0}
\frac{2\pi}{A\,G}\,
\Bigl[\sum_i q_i\, e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,
e^{-G|z - z_i|}\Bigr]\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}}
},
\tag{4}
\]

where \(\phi_0(z)\) is the \(\mathbf{G}=0\) term, which comes from the average planar charge density. For a neutral cell (\(\sum_i q_i = 0\)) the planar average charge density is zero everywhere except possibly for a dipole layer at the surface; integrating Poisson once shows that \(\phi_0(z)\) is piecewise linear:

\[
\phi_0(z) = -\frac{2\pi}{A}\sum_i q_i |z - z_i| + \text{const}.
\]

Outside the slab (say, \(z > \max_i z_i\)) all \(|z - z_i| = z - z_i\), and because \(\sum q_i = 0\) the linear terms cancel, giving a **constant** \(\phi_0\). Similarly on the other side, but the constants may differ if the slab has a net dipole perpendicular to the surface. For a centrosymmetric slab (symmetric under \(z\to -z\)) the two constants are equal.

For the **vacuum above the slab** (let \(z_{\max} = \max_i z_i\)), we can rewrite (4) in a form that explicitly shows the exponential decay in \(z\). For \(z \ge z_{\max}\),

\[
e^{-G|z - z_i|} = e^{-G(z - z_i)} = e^{-G z}\, e^{G z_i} \quad (\text{since } z \ge z_i).
\]

Then

\[
\phi(\boldsymbol{\rho}, z \ge z_{\max}) =
\phi_0 \;+\; \sum_{\mathbf{G}\neq 0}
\frac{2\pi}{A\,G}\,
\Bigl[\sum_i q_i\, e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,
e^{G z_i} \Bigr]\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}}\, e^{-G z}.
\tag{5}
\]

If the system has in‑plane inversion symmetry, the imaginary part cancels and only cosines survive. Then

\[
\phi(\boldsymbol{\rho}, z \ge z_{\max}) =
\phi_0 \;+\; \sum_{\mathbf{G}\neq 0}^{\text{half set}}
C_{\mathbf{G}}\,
\cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho})\, e^{-G z},
\tag{6}
\]

with

\[
C_{\mathbf{G}} = \frac{4\pi}{A\,G} \sum_i q_i\,
\cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho}_i)\,
e^{G z_i}.
\tag{7}
\]

Equation (6) is exactly the “small basis of cosines in \(x,y\) and \(\exp(G z)\) in \(z\)” that you wanted. The coefficients (7) are computed directly from the ionic coordinates and charges — no Ewald summation needed.

## 5. Where do multipoles come in?

Your idea of “expanding in monopole, dipole, quadrupole in \(z\)” can be seen by looking at the **G≠0 components** in a slightly different way. For a fixed \(\mathbf{G}\), the source \(\rho_{\mathbf{G}}(z)\) is a set of delta peaks at different \(z_i\). One could formally write the solution (3) as a convolution of the Green’s function with an effective “surface charge density”. But more illuminating: if the ions are not point charges but have a continuous distribution \(\rho_{\mathbf{G}}(z)\) along \(z\), we could expand \(\rho_{\mathbf{G}}(z)\) in moments:

\[
\rho_{\mathbf{G}}(z) = \underbrace{M_{\mathbf{G}}}_{\text{monopole}}\;\delta(z- z_c)
+ \underbrace{D_{\mathbf{G}}}_{\text{dipole}}\;\delta'(z- z_c)
+ \underbrace{Q_{\mathbf{G}}}_{\text{quadrupole}}\;\delta''(z- z_c) + \dots
\]

The corresponding potential \(\phi_{\mathbf{G}}(z)\) would then involve not only \(e^{-G|z - z_c|}\) but also its derivatives with respect to \(z_c\), which bring down polynomials. For example, a pure dipole layer creates a potential that outside the layer behaves like \((1 + G|z|)e^{-G|z|}\).

In the exact solution for point ions we effectively sum an infinite series of multipoles at each ion position, but the result is fully captured by the simple sum of exponentials in (4). For a slab that is many atomic layers thick, the superposition of many \(e^{-G|z - z_i|}\) terms automatically encodes all the multipole character—there is no need to compute moments explicitly.

If you were to *truncate* not only in \(\mathbf{G}\) but also in \(z\)‑moments (e.g., keep only the monopole and dipole of the slab as a whole for each \(\mathbf{G}\)), you would obtain a coarse model of the potential that misses the fine details near the slab but still reproduces the long‑range decay. The small‑\(\mathbf{G}\) truncation we use already does something similar: the higher the \(G\), the faster the potential decays in \(z\), so keeping only the first few \(\mathbf{G}\) vectors is sufficient for the asymptotic field.

## 6. Why this is a “small basis” and how it connects to Ewald

The representation (6) contains an **infinite** number of \(\mathbf{G}\) vectors. In practice we keep only those with \(|\mathbf{G}| \le G_{\max}\). This is analogous to truncating a Fourier expansion of a smooth function: the error decays quickly because the exact potential is smooth away from the ion cores. In the vacuum region, the neglected terms decay as \(e^{-G_{\max} z}\), which is exponentially small for \(z\) a few Å above the surface.

The full Ewald sum splits the electrostatic potential into a short‑range real‑space sum and a reciprocal‑space sum that includes an infinite set of \(\mathbf{G}\) vectors, but damped by a complementary error function that depends on \(z\). If we set the Ewald splitting parameter to zero, the reciprocal part reduces exactly to the Fourier series in this explanation (with slow convergence). In our “small basis” approach we avoid the Ewald machinery entirely by simply working with the exact Fourier series and truncating it **after** guaranteeing that the retained components already capture the physics we need (e.g., the long‑range surface potential).

## 7. Summary of the algorithm (principle)

1. **Choose a set of small \(\mathbf{G}\) vectors** (e.g., \(|h|,|k| \le 2\) or 3).
2. For each \(\mathbf{G}\), compute the structure factor \(S_{\mathbf{G}}(z) = \sum_i q_i e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i} \delta(z-z_i)\).
3. Solve the 1D Poisson equation for that component: the solution is a sum of exponentials \(e^{-G|z - z_i|}\) weighted by the structure factor.
4. In the vacuum region, the result simplifies to a pure exponential decay times a cosine in‑plane, with coefficients \(C_{\mathbf{G}}\) given by (7).
5. Plot/use the truncated sum. The error is controlled by \(e^{-G_{\max} \Delta z}\), where \(\Delta z\) is the distance from the shallowest ion to the point of interest.

This is rigorous, physically transparent, and yields an explicit formula that can be implemented efficiently with vectorised NumPy operations — exactly as in the provided code. The method connects seamlessly with the multipole viewpoint: the in‑plane Fourier components are “plane‑wave multipoles”, and their \(z\)‑dependence is exactly solved by the Green’s function of the 1D screened Poisson equation.

---

# Gemini-3-flash

This is a very insightful way to frame the problem. By decomposing the system into a product of "lateral textures" (the cosines) and "axial moments" (the $z$-distribution), you are effectively performing a **2D Fourier-Bessel-like decomposition**.

Here is the rigorous derivation followed by a didactic breakdown.

---

### 1. The Rigorous Derivation

We start with the 3D Poisson equation for the potential $\Phi(\mathbf{r}, z)$:


$$\nabla^2 \Phi(\mathbf{r}, z) = -4\pi \rho(\mathbf{r}, z)$$


where $\mathbf{r} = (x, y)$.

#### Step A: Fourier Transform the XY Plane

Since the charge density $\rho$ is periodic in $xy$, we expand it as a Fourier series:


$$\rho(\mathbf{r}, z) = \sum_{\mathbf{G}} \rho_{\mathbf{G}}(z) e^{i\mathbf{G} \cdot \mathbf{r}}$$


For a center-symmetric system, we use $\cos(\mathbf{G} \cdot \mathbf{r})$. The coefficients $\rho_{\mathbf{G}}(z)$ represent the "vertical profile" of each specific lateral frequency.

#### Step B: The Ansatz

We assume the potential follows the same lateral symmetry:


$$\Phi(\mathbf{r}, z) = \sum_{\mathbf{G}} \phi_{\mathbf{G}}(z) e^{i\mathbf{G} \cdot \mathbf{r}}$$

#### Step C: Solve for each G-mode

Substitute these into the Poisson equation. For a specific $\mathbf{G}$, the Laplacian $\nabla^2 = \nabla_{xy}^2 + \frac{\partial^2}{\partial z^2}$ acts as:


$$(-G^2 + \frac{d^2}{dz^2}) \phi_{\mathbf{G}}(z) = -4\pi \rho_{\mathbf{G}}(z)$$


This is a 1D inhomogeneous ODE. The solution for a source $\rho_{\mathbf{G}}(z)$ is found using the Green's function for the operator $(\frac{d^2}{dz^2} - G^2)$, which is $\frac{-2\pi}{G} e^{-G|z-z'|}$.

Thus, for a specific mode $\mathbf{G}$:


$$\phi_{\mathbf{G}}(z) = \frac{2\pi}{G} \int \rho_{\mathbf{G}}(z') e^{-G|z-z'|} dz'$$

---

### 2. The Didactic Breakdown

Let’s look at what this integral actually says by using your **Multipole Expansion** idea.

#### The "Layer" Concept

Imagine one single mode $\mathbf{G}$ (e.g., a simple checkerboard pattern of charge). As you move away from the surface ($z > z'$), the term $|z-z'|$ becomes $(z-z')$. We can rewrite the integral:


$$\phi_{\mathbf{G}}(z) = \frac{2\pi}{G} e^{-Gz} \int \rho_{\mathbf{G}}(z') e^{Gz'} dz'$$

#### Connecting to Multipoles

Now, look at the term inside the integral: $e^{Gz'}$. If the slab is thin compared to the wavelength of the mode (i.e., $Gz'$ is small), we can Taylor expand the exponential:


$$e^{Gz'} \approx 1 + Gz' + \frac{1}{2}(Gz')^2 + \dots$$

Plugging this back into the integral for $\phi_{\mathbf{G}}(z)$:


$$\phi_{\mathbf{G}}(z) = \frac{2\pi e^{-Gz}}{G} \left[ \underbrace{\int \rho_{\mathbf{G}}(z') dz'}_{\text{Monopole of mode } \mathbf{G}} + G \underbrace{\int z' \rho_{\mathbf{G}}(z') dz'}_{\text{Dipole of mode } \mathbf{G}} + \dots \right]$$

#### Why this is intuitive:

1. **The Monopole Term:** This isn't the total charge of the cell; it is the **amplitude of the charge modulation** at frequency $\mathbf{G}$. If you have a layer of $(+)$ and $(-)$ ions in a plane, the "monopole of the mode" is high.
2. **The Dipole Term:** If you have two layers (e.g., a layer of cations slightly above a layer of anions), the "dipole of the mode" describes how that specific frequency is shifted or polarized vertically.
3. **The Decay:** Notice that the higher the frequency ($G$), the faster the decay ($e^{-Gz}$). This acts as a **natural low-pass filter**. As you move away from the surface, the "fine details" (large $G$) disappear, and you only see the "smooth" potential (small $G$).

---

### 3. Summary Table for Implementation

| Component | Mathematical Role | Physical Meaning |
| --- | --- | --- |
| **Basis in XY** | $\cos(\mathbf{G} \cdot \mathbf{r})$ | The "texture" of the ionic arrangement. |
| **Basis in Z** | $e^{-Gz}$ | How far the influence of that texture reaches into space. |
| **The Coefficient** | $\int \rho_{\mathbf{G}}(z') e^{Gz'} dz'$ | The "Source Strength" of that texture, weighted by its vertical distribution. |

**In practice:** To compute this efficiently, you don't even need the multipole expansion. You just sum $q_j e^{G z_j}$ for all ions $j$ in the unit cell. This implicitly captures all multipole moments (monopole, dipole, quadupole...) of the $z$-distribution simultaneously!

---

# USER

Ok can you put this framework in Python numpy Code and illustrate it with matplotlib imshow

---

# Gemini-3-flash

To implement this framework, we will treat the potential as a sum of independent lateral "modes" ($\mathbf{G}$-vectors). For each mode, we calculate its source strength (the vertical moments) and then multiply it by the corresponding cosine wave and exponential decay.

The code below is organized to mirror the didactic steps we discussed: mapping the lattice, calculating the modal source strengths, and reconstructing the 3D potential.

```python
import numpy as np
import matplotlib.pyplot as plt

def compute_fourier_potential(a, b, ions, z_eval, n_max=8, grid_res=128):
    """
    a, b: Lattice constants for a rectangular cell (x, y)
    ions: List of (x, y, z, q)
    z_eval: Height at which to evaluate the XY slice
    n_max: Fourier cutoff (number of G-shells)
    """
    # 1. Setup Grids
    x = np.linspace(0, a, grid_res)
    y = np.linspace(0, b, grid_res)
    X, Y = np.meshgrid(x, y)
    Area = a * b
    
    # 2. Reciprocal Lattice Vectors (for rectangular cell)
    # G = (2*pi*n/a, 2*pi*m/b)
    
    potential_xy = np.zeros_like(X)
    
    # 3. Double loop over Reciprocal space (n, m)
    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            if n == 0 and m == 0:
                continue # Skip DC component for neutrality
            
            # Define G vector and its magnitude
            gx = 2 * np.pi * n / a
            gy = 2 * np.pi * m / b
            g_mag = np.sqrt(gx**2 + gy**2)
            
            # 4. Compute the "Source Strength" for this mode
            # This is the integral of rho_G(z) * exp(G*z) dz
            # For point charges, the integral becomes a sum:
            # S_G = sum( q_j * exp(i G.r_j) * exp(g_mag * z_j) )
            # We use cos for center-symmetric systems.
            
            source_strength = 0
            for ix, iy, iz, iq in ions:
                # Cosine basis for the lateral modulation
                lateral_phase = gx * ix + gy * iy
                # Vertical moment contribution
                source_strength += iq * np.cos(lateral_phase) * np.exp(g_mag * iz)
            
            # 5. Build the Potential term
            # Phi_G = (2*pi / (Area * G)) * source_strength * cos(G.r) * exp(-G*z)
            prefactor = (2 * np.pi) / (Area * g_mag)
            term = prefactor * source_strength * np.cos(gx * X + gy * Y) * np.exp(-g_mag * z_eval)
            
            potential_xy += term
            
    return X, Y, potential_xy

# --- Example: CsCl (100) Surface (Center-symmetric) ---
# Cs+ at (0,0,0), Cl- at (a/2, a/2, 0)
a_lat = 4.12
ions_cscl = [
    (0.0, 0.0, 0.0, 1.0),       # Cs+
    (a_lat/2, a_lat/2, 0, -1.0) # Cl-
]

# Evaluate at different heights to see the "Low-Pass Filter" effect
z_heights = [0.5, 1.0, 2.0]
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, z in zip(axes, z_heights):
    X, Y, phi = compute_fourier_potential(a_lat, a_lat, ions_cscl, z_eval=z)
    im = ax.imshow(phi, extent=[0, a_lat, 0, a_lat], origin='lower', cmap='RdBu_r')
    ax.set_title(f'Potential at z = {z} Å')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    plt.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()

```

### Explanation of the Implementation

1. **The Nested Loops (`n`, `m`)**: These represent the "Selection" of our Fourier basis. By increasing `n_max`, you include higher-frequency cosines, which resolve the sharp peaks near the ions.
2. **`source_strength`**: This variable calculates the modal amplitude for each cosine. It combines the lateral position (via `cos`) and the vertical position (via `exp(g_mag * iz)`). If you had a dipole layer, the `iz` values would differ, and the `exp` terms would automatically capture the dipole moment we derived.
3. **The Decay Logic**: Notice in the `term` calculation that we have `np.exp(-g_mag * z_eval)`. This is why the plots at $z=2.0$ look much "smoother" than those at $z=0.5$. The high-frequency terms (large `g_mag`) are crushed by the exponential decay much faster than the low-frequency terms.

### Visualizing the XZ Slice

To truly see the "exp(b.z)" behavior you mentioned, you can modify the plotting section to hold `y=0` constant and vary `z`:

```python
z_range = np.linspace(0.2, 4.0, 100)
x_range = np.linspace(0, a_lat, 100)
XZ_X, XZ_Z = np.meshgrid(x_range, z_range)

# Re-run potential calculation for the XZ plane
phi_xz = np.zeros_like(XZ_X)
# ... (same loop logic as above, but use XZ_X and XZ_Z) ...

plt.figure(figsize=(6, 8))
plt.contourf(XZ_X, XZ_Z, phi_xz, 50, cmap='RdBu_r')
plt.title("XZ Slice: Potential Decay")
plt.xlabel("x (Å)")
plt.ylabel("z (Height above surface)")
plt.show()

```

This will show the characteristic "lobes" of potential emerging from the ions and smoothing out into a featureless vacuum potential as $z$ increases.

---

# USER

Two more ideas:
1) Laplace transform. Can we use Laplace transform to solve this or derive it? I mean Fourier transform has basis exp(ikt) and Laplace exp(kt) we can generally considered complex exponent (representing both oscillating periodic solution in x,y and decaying solution in z)
2)dipolar layer - maybe the solver would be simpler and more numerical robust \ stable and accurate if we consider instead of charged atoms with different z coordinates just 2d layer of dipoles which have all same z coordinate 

---

# DeepSeek-V4

Excellent ideas! Both the Laplace‑transform viewpoint and the dipole‑layer approximation fit naturally into the framework we developed. Let’s explore them in detail, keeping the explanation as physical as possible.

---

## 1. Laplace transform in the \(z\)‑direction

You are right: the method we used is equivalent to a **mixed Fourier–Laplace transform**.  
In the \(xy\)‑plane we have discrete translational symmetry, so we use a **Fourier series** (basis \(e^{i\mathbf{G}\cdot\boldsymbol{\rho}}\)).  
In the perpendicular \(z\)‑direction there is no periodicity; instead we have a boundary condition that the potential **decays to a constant at infinity**. The natural transform for a function that decays (or is bounded) on a half‑line is the **Laplace transform**, or more generally an exponential basis with a real exponent.

Let’s formalise this. We want to solve

\[
\nabla^2 \phi = -4\pi\rho,
\]

with \(\rho\) a periodic set of point charges in a slab. Take the 2D Fourier transform in the plane:

\[
\widetilde{\phi}_{\mathbf{G}}(z) = \frac{1}{A}\int_{\text{cell}} \phi(\boldsymbol{\rho},z)\,e^{-i\mathbf{G}\cdot\boldsymbol{\rho}}\,d^2\rho,
\qquad
\widetilde{\rho}_{\mathbf{G}}(z) = \frac{1}{A}\sum_i q_i e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,\delta(z-z_i).
\]

The Poisson equation decouples into independent 1D problems for each \(\mathbf{G}\):

\[
\bigl(\partial_z^2 - G^2\bigr)\,\widetilde{\phi}_{\mathbf{G}}(z) = -4\pi\,\widetilde{\rho}_{\mathbf{G}}(z),
\tag{* }
\]

with \(G = |\mathbf{G}|\). For \(\mathbf{G}=0\) we get \(\partial_z^2 \widetilde{\phi}_0 = -4\pi\widetilde{\rho}_0\), which is just the average planar charge; we treat it separately.

Equation (*) is a **linear ODE with constant coefficients**. The homogeneous solution would be

\[
\widetilde{\phi}_{\mathbf{G}}(z) = A\, e^{G z} + B\, e^{-G z}.
\]

The source \(\widetilde{\rho}_{\mathbf{G}}(z)\) consists of delta spikes. Because we require the potential to decay away from the slab, we choose the half‑line Green’s function

\[
g(z-z') = \frac{2\pi}{G}\, e^{-G|z-z'|},
\]

which is itself a combination of \(e^{G(z-z')}\) and \(e^{-G(z-z')}\) that vanishes at \(\pm\infty\). This Green’s function solves \((\partial_z^2 - G^2)g = -4\pi\delta\) and is the **Laplace‑transform** (or **exponential‑transform**) solution for an impulse. Convolving with the sources yields the exact \(\widetilde{\phi}_{\mathbf{G}}(z)\) we used.

Thus, the idea “Fourier in the plane, Laplace in \(z\)” is precisely correct. The exponent \(-G|z-z_i|\) is of Laplace type because \(G\) is **real and positive**, representing a decaying (not oscillating) mode. Formally, one could take a Laplace transform of (*) in the variable \(|z|\) and obtain the same algebraic solution.

So the overall “small basis” is nothing but:

\[
\boxed{\, e^{i\mathbf{G}\cdot\boldsymbol{\rho}} \cdot e^{\pm G z}\,}
\]

i.e. a **Cartesian product** of an in‑plane Fourier mode with an out‑of‑plane real exponential — exactly what you envisioned.

## 2. Dipole‑layer simplification

Your second idea is also very natural. The charge density for a given \(\mathbf{G}\neq 0\) can be expanded in **multipoles** along \(z\). For a neutral slab, the monopole moment vanishes. The leading asymptotic behaviour outside the slab is then determined by the **dipole moment** of that Fourier component.

Let’s see how this appears. For a given \(\mathbf{G}\), the source is a collection of delta functions:

\[
\widetilde{\rho}_{\mathbf{G}}(z) = \frac{1}{A}\sum_i q_i e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\,\delta(z-z_i).
\]

We can expand this distribution about some central height \(z_0\) (e.g. the centre of the slab):

\[
\widetilde{\rho}_{\mathbf{G}}(z) \approx \underbrace{M_{\mathbf{G}}}_{\text{monopole}}\delta(z-z_0)
\;+\; \underbrace{D_{\mathbf{G}}}_{\text{dipole}}\delta'(z-z_0)
\;+\; \underbrace{Q_{\mathbf{G}}}_{\text{quadrupole}}\delta''(z-z_0) + \cdots
\]

where

\[
M_{\mathbf{G}} = \frac{1}{A}\sum_i q_i e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i},
\qquad
D_{\mathbf{G}} = -\frac{1}{A}\sum_i q_i (z_i - z_0) e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i},
\qquad
\text{etc.}
\]

If the unit cell is charge‑neutral (\(\sum q_i = 0\)), then \(M_{0}=0\), but for \(\mathbf{G}\neq 0\) the monopole \(M_{\mathbf{G}}\) is generally non‑zero because the charge is modulated in the plane. However, far from the slab, the decay is still dominated by the lowest non‑vanishing moment. The exact solution gives for each \(\mathbf{G}\) a sum over all ions:

\[
\widetilde{\phi}_{\mathbf{G}}(z) = \frac{2\pi}{A\,G}\sum_i q_i e^{-i\mathbf{G}\cdot\boldsymbol{\rho}_i}\, e^{-G|z-z_i|}.
\]

If we replace this by the **dipole approximation** — i.e. we keep only the terms from \(M_{\mathbf{G}}\,\delta(z-z_0)\) and \(D_{\mathbf{G}}\,\delta'(z-z_0)\) — the potential becomes:

\[
\widetilde{\phi}_{\mathbf{G}}(z) \approx
\frac{2\pi}{G}\Bigl[
M_{\mathbf{G}}\, e^{-G|z-z_0|}
\;+\;
D_{\mathbf{G}}\, \frac{d}{dz_0}\bigl(e^{-G|z-z_0|}\bigr)
\Bigr].
\]

The derivative brings down a factor \(-G\,\operatorname{sgn}(z-z_0)\):

\[
\frac{d}{dz_0} e^{-G|z-z_0|} = G\,\operatorname{sgn}(z-z_0)\, e^{-G|z-z_0|}.
\]

Thus

\[
\widetilde{\phi}_{\mathbf{G}}(z) \approx
\frac{2\pi}{G}\, e^{-G|z-z_0|}
\Bigl[ M_{\mathbf{G}} \;+\; D_{\mathbf{G}}\, G\,\operatorname{sgn}(z-z_0) \Bigr].
\]

Outside the slab the sign function is constant (\(\operatorname{sgn}=+1\) above, \(-1\) below). Therefore, for **each \(\mathbf{G}\)** we obtain a single exponential decay \(\propto e^{-G|z|}\) multiplied by a **linear combination** of the monopole and dipole moments. This is indeed a “2D layer of dipoles” — all localised at the same height \(z_0\) — and the \(z\)‑dependence remains exactly of the form \(e^{-G|z|}\) we want to use as a basis.

### Is this numerically more robust or accurate?

- **Accuracy:** The exact sum over ions is exact (for a given \(\mathbf{G}\) set), while the dipole approximation throws away the quadrupole and higher moments of the \(z\)‑distribution. Near the slab, the approximation can be poor because the field varies on the scale of the ionic separation. Far away, the higher moments decay as higher powers of \(e^{-G\Delta z}\), so the dipole (and monopole) dominate, and the approximation becomes excellent. So it is a *far‑field* simplification, not an overall improvement.

- **Numerical robustness:** The exact formula (sum over ions) involves only exponentials of the form \(e^{-G|z-z_i|}\), which are well‑behaved and easy to compute. No instabilities arise because the arguments are bounded. The dipole formula would reduce the number of exponentials to evaluate (one per \(\mathbf{G}\) instead of one per ion), which is marginally faster but at the cost of losing near‑field information. There is no inherent advantage in stability.

- **The real simplification** comes if we also **truncate the multipole series** to a single effective dipole layer for *all* \(\mathbf{G}\) combined — i.e. we treat the whole slab as one thin dipole sheet with a \(z\)-independent potential jump. But that would completely miss the in‑plane modulation and the exponential decay of the individual \(\mathbf{G}\) components. For our purpose (a small basis of \(\cos\) and \(\exp\)) we need the \(\mathbf{G}\)-by-\(\mathbf{G}\) treatment anyway, so the dipole layer offers little extra.

Thus, the dipole‑layer viewpoint is conceptually very satisfying: it shows that the charge distribution of each Fourier component acts *as if* it were a planar array of monopoles and dipoles, and that the vacuum potential is simply a superposition of their exponential fields. But for implementation, we stick to the direct ion sum because it is already simple, exact, and numerically safe.

## 3. Combining the two pictures

Your two suggestions actually complete each other:

- The **Fourier–Laplace** decomposition tells us the *form* of the basis functions: \(\cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho})\,e^{-G z}\).
- The **multipole expansion in \(z\)** tells us how the coefficients of those basis functions are related to the ionic distribution: the coefficient for a given \(\mathbf{G}\) is essentially the **total effective dipole moment** of the slab’s charge modulation at that wave‑vector, projected onto the cosine basis.

In the simple algorithm we use, the coefficient \(C_{\mathbf{G}}\) is computed directly from the ions via

\[
C_{\mathbf{G}} = \frac{4\pi}{A\,G} \sum_i q_i \cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho}_i)\, e^{G z_i}.
\]

If we had chosen a reference plane \(z_0\) and expanded, we could rewrite this as

\[
C_{\mathbf{G}} = \frac{4\pi}{A\,G} e^{G z_0} \sum_i q_i \cos(\mathbf{G}\!\cdot\!\boldsymbol{\rho}_i)\, e^{G(z_i-z_0)},
\]

and then Taylor expand \(e^{G(z_i-z_0)}\) to see the series of multipole contributions. The leading constant (monopole) plus linear-in-\((z_i-z_0)\) (dipole) give the far‑field approximation. The exact formula includes all orders automatically.

So the whole method is:

> **Fourier in \(xy\) and multipole (exponential) decay in \(z\), with coefficients taken from the exact ionic positions.**

The “small basis” is just a truncation in \(\mathbf{G}\); the \(z\)‑behaviour is exact for each retained \(\mathbf{G}\). No Ewald summation is needed, and the implementation is both physically transparent and numerically stable.

If you ever wish to further compress the model, you could indeed replace the ion‑by‑ion sum for each \(\mathbf{G}\) by a single dipole layer at the slab’s centre, giving an even more compact expression — but at the price of quantitative accuracy only far from the surface. For most surface science applications, the direct ion summation with a handful of \(\mathbf{G}\) vectors already provides an excellent balance of speed and accuracy.

---

# Gemini-3-flash

These are both high-level insights that bridge the gap between "standard" physics and engineering-style signal processing. Let’s explore both.

### 1. The Laplace Transform Perspective

You are exactly right: the solution is a **Complex Exponential** in 3D.

If we consider the Laplace transform in the $z$-direction, we are transforming $\frac{\partial^2}{\partial z^2} \to s^2$. In our Poisson equation $\nabla_{xy}^2 \Phi + \frac{\partial^2 \Phi}{\partial z^2} = -4\pi \rho$:

1. We apply a **Fourier Transform** laterally ($\nabla_{xy}^2 \to -G^2$).
2. We apply a **Laplace Transform** vertically ($\frac{\partial^2}{\partial z^2} \to s^2$).

The equation in "Frequency-Laplace" space becomes:


$$(-G^2 + s^2) \Phi(G, s) = -4\pi \rho(G, s)$$

The "Natural Frequencies" of this system (the poles) occur where $s^2 = G^2$, or **$s = \pm G$**.

* In Laplace theory, the real part of the variable $s$ determines growth or decay.
* Because we require the potential to be finite at $z \to \infty$, we must choose the pole $s = -G$.
* This is why the potential **must** decay as $e^{-Gz}$. The "Laplace frequency" in the vertical direction is locked to the "Fourier frequency" in the lateral direction.

---

### 2. The Dipolar Layer Approximation

Representing ions as a **2D sheet of dipoles** is a powerful trick used in surface science (the Helmholtz Layer model). It is often more robust because it removes the "monopole divergence" (where a total charge on a surface would create an infinite potential in 3D).

#### The Physics

A dipole $\mathbf{p}$ at height $z=0$ can be thought of as two charges $+q$ and $-q$ separated by a distance $d$. In the Fourier limit, the potential of a dipole layer with dipole density $P_{\mathbf{G}}$ is:


$$\Phi_{\mathbf{G}}(z) = 2\pi P_{\mathbf{G}} e^{-Gz}$$

Notice there is no $1/G$ term here (unlike the monopole case). This makes the summation **numerically more stable** because the $G \to 0$ limit doesn't blow up; it simply represents the work function shift (the constant potential jump across a dipole layer).

#### Python Implementation (Dipole Sheet)

Here is how you would adapt the code to handle a unit cell of dipoles $p_j$ located at positions $(x_j, y_j)$ all at $z=0$.

```python
import numpy as np
import matplotlib.pyplot as plt

def compute_dipole_potential(a, b, dipoles, z_eval, n_max=10):
    """
    dipoles: list of (x, y, p_z) where p_z is the dipole moment strength
    """
    x = np.linspace(0, a, 100)
    y = np.linspace(0, b, 100)
    X, Y = np.meshgrid(x, y)
    Area = a * b
    phi = np.zeros_like(X)

    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            if n == 0 and m == 0: continue
            
            gx, gy = 2*np.pi*n/a, 2*np.pi*m/b
            g_mag = np.sqrt(gx**2 + gy**2)
            
            # Modal Dipole Strength
            # P_G = sum( p_j * cos(G.r_j) )
            P_G = sum(p * np.cos(gx*dx + gy*dy) for dx, dy, p in dipoles)
            
            # Dipole potential decays as e^{-Gz} without the 1/G factor
            # This makes it very stable for small G!
            term = (2 * np.pi * P_G / Area) * np.cos(gx * X + gy * Y) * np.exp(-g_mag * z_eval)
            phi += term
            
    return X, Y, phi

# Example: A checkerboard of up-down dipoles
a_lat = 5.0
dipoles = [(0, 0, 1.0), (2.5, 2.5, -1.0)] 

X, Y, phi = compute_dipole_potential(a_lat, a_lat, dipoles, z_eval=0.5)

plt.imshow(phi, extent=[0, a_lat, 0, a_lat], cmap='PuOr')
plt.title("Potential from a Periodic Dipole Layer")
plt.colorbar(label="Potential (V)")
plt.show()

```

### Why use Dipoles?

1. **Stability:** As mentioned, you avoid the $1/G$ singularity.
2. **Physical Insight:** Most reconstructed surfaces or adsorbed layers are better described by their dipole moment (e.g., a water molecule on a surface).
3. **Boundary Conditions:** A dipole layer creates a discrete "step" in potential. If you want to model the **Work Function** change of a crystal, the dipole layer is the standard mathematical tool to do it.

**Which approach fits your physical system better?** If you are modeling a rigid ionic crystal (like NaCl), the first (monopole) method is more "honest" to the atomic structure. If you are modeling an organic monolayer or a relaxed surface with polarization, the dipole method is far superior.

---

# Implementation Report: 2D Ewald Summation in FireCore

## Overview

A rigorous implementation of 2D Ewald summation for electrostatic potential of periodic ionic slabs has been added to FireCore. The implementation consists of:

1. **Backend module**: `pyBall/Ewald2D.py` - Pure Python/NumPy implementation of the 2D Ewald formulas
2. **Test script**: `tests/tEwald2D/test_ewald_2d.py` - CLI interface for evaluating potentials on arbitrary .xyz files
3. **GPU test**: `tests/tEwald2D/test_ewald_2d_ocl.py` - Maps Ewald coefficients to folded basis and evaluates on GPU via pyOpenCL

## Mathematical Framework

The implementation follows the exact 2D Fourier representation derived in the discussion above. For a slab periodic in the xy-plane with lattice vectors **a**, **b** and ions at positions (ρ_i, z_i) with charges q_i:

### Vacuum potential (z ≥ z_max)

```
φ(ρ, z) = φ₀ + Σ_{G≠0} C_G · cos(G·ρ) · exp(-|G|·z)

where:
C_G = (4π/(A|G|)) Σ_i q_i cos(G·ρ_i) exp(|G|·z_i)
```

### Full interior potential

```
φ(ρ, z) = φ₀(z) + Σ_{G≠0} (2π/(A|G|)) Σ_i q_i e^{iG·(ρ-ρ_i)} e^{-|G||z-z_i|}
```

where φ₀(z) = -(2π/A) Σ_i q_i |z - z_i| + const is the G=0 contribution.

## Implementation Details

### pyBall/Ewald2D.py

This is the core backend module providing:

#### Core Functions
- `make_reciprocal_2d(a, b)`: Computes area and reciprocal vectors b1, b2
- `generate_G_vectors(b1, b2, n_harm)`: Generates G-vectors for |h|,|k| ≤ n_harm
- `compute_C_G(...)`: Vacuum-side coefficients for centrosymmetric systems
- `compute_w_per_ion(...)`: Per-ion complex weights for full interior potential
- `eval_potential_vacuum(...)`: Vectorized evaluation on grid points
- `eval_potential_full(...)`: Full evaluation including interior region

#### Ewald2D Class
High-level interface that encapsulates all computations:
```python
ew = Ewald2D.from_AtomicSystem(sys_at, n_harm=3)
phi = ew.eval_potential_vacuum(xyz_grid)
phi_full = ew.phi_full_1d(x, y, z_array)  # 1D line scan
```

The class precomputes:
- Reciprocal lattice vectors (b1, b2)
- G-vector set (Gx, Gy, Gn, hs, ks)
- Structure factors (C_G for vacuum, w_G for full)
- Charge density reconstruction on grid

#### Key Design Decisions
- **Vectorized NumPy operations**: All loops over ions are vectorized for performance
- **Separate vacuum vs full formulas**: Vacuum uses simple exp(-Gz); full uses exp(-G|z-z_i|)
- **Units**: Å for lengths, elementary charge for q → potential in e/Å (multiply by 14.4 for eV/e)

### tests/tEwald2D/test_ewald_2d.py

Originally a monolithic script with inline calculations, refactored into a thin CLI wrapper that:
1. Loads .xyz files via `pyBall.AtomicSystem` (parses lattice vectors from comment line "lvs=...")
2. Delegates all computations to `Ewald2D` backend
3. Generates diagnostic plots (XY slices, XZ slices, charge density, error maps)
4. Supports CLI arguments for input file, harmonic cutoff, grid resolution, plotting control

**Refactoring benefits**:
- Reusable backend module for other applications
- Cleaner separation of concerns
- Easier testing and maintenance
- Consistent interface with other FireCore modules

### tests/tEwald2D/test_ewald_2d_ocl.py

This script demonstrates integration with the GPU-accelerated folded basis framework:

#### Purpose
Map the analytical Ewald coefficients into the folded basis format used by the `getSurfFolded` OpenCL kernel in `pyBall/OCL/MolecularDynamics.py`.

#### Two Fitting Methods

**Method 1: Direct (`--method direct`)**
- Uses one basis function per unique G-vector half-pair (h ≥ 0)
- Basis params: (kx, ky, alpha, z0) where kx = |h|, ky = |k|, alpha = |G|
- For orthogonal lattice with mirror symmetry, this is EXACT
- Coefficients fitted via least-squares to Ewald reference

**Method 2: coscos (`--method coscos`)**
- Uses regular tensor-product basis: cos(2π·kx·u)·cos(2π·ky·v)·exp(-α·(z-z0))
- kx = 0,...,nu-1; ky = 0,...,nv-1; α = z_scale·(1,2,...,nz)
- Approximate fitting for general crystals
- Controlled by nu, nv, nz parameters

#### Evaluation Pipeline
1. Load system via `AtomicSystem`
2. Compute Ewald2D coefficients via `Ewald2D` backend
3. Fit coefficients into folded basis (direct or coscos method)
4. Upload to GPU via `MolecularDynamics._set_folded_coefficients`
5. Evaluate potential on grid using `MolecularDynamics.eval_rigid_getSurfFolded`
6. Compare GPU results with Python references (Ewald2D and folded Python evaluation)

## Test Results

### NaCl (1×1 unit cell, 4 ions)

**Test parameters**: n_harm=3, grid=100×100

| Metric | Value |
|--------|-------|
| G-vectors (N_G) | 48 |
| Basis functions (direct method) | 24 |
| Fit error (folded vs Ewald) | ~1×10⁻⁸ |
| GPU error (vs Ewald) | ~1×10⁻⁷ (float32 precision) |
| GPU vs Python folded | ~1×10⁻⁷ |
| Python Ewald (vectorized) | 12 ms (10k points) |
| Python folded | 1320 ms |
| GPU evaluation | 2 ms |
| Speedup vs Python Ewald | 6× |
| Speedup vs Python folded | 660× |

**Z-scan comparison**: All three curves (Ewald, Python folded, GPU) overlap perfectly with ~2×10⁻⁷ max difference.

### CaF2 (fluorite structure, 24 ions)

**Test parameters**: n_harm=4

**Mirror symmetry check**: C(h,k) ≠ C(h,-k) for most (h,k) pairs, indicating lack of mirror symmetry in the unit cell.

**Direct method results**:
- Fit error (XY grid): 0.15 (large error due to missing sin·sin terms)
- Z-scan error: 1.3×10⁻⁷ (exact at (0,0) where sin terms vanish)
- GPU vs Python folded: 7.2×10⁻⁸ (GPU correctly reproduces the folded basis)

The large XY grid error is expected because the folded basis `cos(2π·h·u)·cos(2π·k·v)` cannot represent `cos(2π(h·u + k·v))` when both h≠0 and k≠0 for non-mirror-symmetric charge distributions. The missing sin·sin cross-terms are significant for CaF2.

**coscos method results**:
With nu=4, nv=4, nz=8 (128 basis functions), the fit quality is still poor (max error ~0.16) due to the fundamental limitation of the cos·cos tensor product for representing general Fourier components.

## Problems Encountered and Solutions

### Problem 1: Slow reference evaluation in test script
**Issue**: Original test used `phi_full_1d` in a loop over grid points, which was extremely slow for large grids (each call re-sums over all G-vectors and ions).

**Solution**: Vectorized the evaluation using the vacuum formula directly:
```python
phase = Gx[:, None]*X_flat[None, :] + Gy[:, None]*Y_flat[None, :]
decay = np.exp(-Gn[:, None] * Z_flat[None, :])
phi_ewald = np.sum(C_G[:, None] * np.cos(phase) * decay, axis=0)
```

This reduced evaluation time from minutes to milliseconds for typical grids.

### Problem 2: Mirror symmetry limitation of folded basis
**Issue**: The `getSurfFolded` kernel evaluates `cos(2π·kx·u)·cos(2π·ky·v)` which equals `cos(2π(h·u + k·v))` only when the charge distribution has mirror symmetry (C(h,k) = C(h,-k)). For crystals like CaF2, this condition fails.

**Analysis**: The identity `cos(2π(h·u + k·v)) = cos(2πh·u)cos(2πk·v) - sin(2πh·u)sin(2πk·v)` shows that the sin·sin term is missing. Since the kernel only supports cos·cos products (negative k values give the same result due to cos(-x)=cos(x)), it cannot represent general Fourier components.

**Solution**: Documented the limitation and noted that a new kernel variant evaluating `cos(Gx·x + Gy·y)` directly would be needed for the general case. For mirror-symmetric crystals (NaCl, simple rock-salt structures), the current kernel is exact.

### Problem 3: GPU float32 precision
**Issue**: GPU results differ from Python float64 references by ~1×10⁷, which is larger than the fitting error (~1×10⁸).

**Analysis**: This is the expected precision limit of single-precision floating-point arithmetic on GPU. The error pattern is random/unstructured, consistent with rounding noise rather than algorithmic errors.

**Solution**: Documented that GPU accuracy is limited to float32 precision (~7 decimal digits), which is acceptable for most MD applications. For higher precision, double-precision kernels would be needed (at significant performance cost).

## Limitations

### 1. Mirror symmetry requirement for exact folded-basis representation
The current `getSurfFolded` kernel with cos·cos tensor product can exactly represent the Ewald potential ONLY for crystals with mirror symmetry in both u and v directions (i.e., C(h,k) = C(h,-k) = C(-h,k) = C(-h,-k)). This includes:
- Simple rock-salt structures (NaCl)
- Simple cubic lattices with symmetric charge arrangements

For general crystals lacking this symmetry (CaF2 fluorite, many perovskites, reconstructed surfaces), the representation is approximate. The error depends on how strongly the sin·sin cross-terms contribute to the potential.

**Workaround**: Use higher harmonic counts in the coscos fitting method, or implement a new kernel variant that evaluates `cos(Gx·x + Gy·y)` directly (which would handle all cases exactly).

### 2. G-vector truncation
The Ewald sum is truncated to |h|,|k| ≤ n_harm. The error decays as exp(-G_max·Δz) where Δz is distance from the slab. For typical surfaces and n_harm=3-4, the error above ~2 Å from the surface is < 1×10⁻⁴.

### 3. Point charge assumption
The implementation treats ions as point charges. For high-precision work near ionic cores, electron density corrections may be needed.

## Integration with FireCore

The Ewald2D module integrates seamlessly with existing FireCore infrastructure:

- **AtomicSystem**: Loads .xyz files with lattice vectors and charges
- **MolecularDynamics**: GPU folded-basis evaluation framework
- **Plotting utilities**: Uses `plotUtils.py` for visualization
- **Build system**: No compilation needed (pure Python/NumPy)

## Future Extensions

Potential improvements:
1. **Full cos(G·ρ) kernel**: Implement a new OpenCL kernel variant that evaluates `cos(Gx·x + Gy·y)` directly for exact representation of general crystals
2. **Dipole layer support**: Add option to represent ions as effective dipoles for faster far-field evaluation
3. **Double-precision GPU kernels**: For applications requiring higher accuracy
4. **Interface with Fireball**: Use DFT-computed charge densities instead of point charges
5. **Force calculation**: Add gradient computation for MD applications (analytical gradients are straightforward to derive from the formulas)

## Conclusion

The 2D Ewald implementation provides a rigorous, efficient, and well-tested foundation for electrostatic potential calculations of periodic ionic surfaces. The backend module (`pyBall/Ewald2D.py`) offers both vacuum and full interior potential evaluation with vectorized NumPy operations. Integration with the GPU folded-basis framework demonstrates the path to high-performance evaluation, with the caveat that the current cos·cos kernel is exact only for mirror-symmetric crystals. The implementation is production-ready for surface science applications and serves as a solid base for future extensions.