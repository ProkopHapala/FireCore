# FireBall Hamiltonian Assembly Mapping

This document provides a technical mapping between the theoretical energy expressions in FireBall and the Fortran implementation. It is intended to serve as a guide for reimplementing the core Hamiltonian assembly in OpenCL.

## 1. Hamiltonian Decomposition

The total Hamiltonian $H$ is assembled as a sum of several physical contributions:

$$H = T + V_{na} + V_{nl} + V_{xc} + V_{xc\_1c} + V_{ca} + V_{xc\_ca} + V_{Ewald}$$

| Term | Description | Fortran Variable | Main Subroutine | Core Interpolation |
| :--- | :--- | :--- | :--- | :--- |
| $T$ | Kinetic Energy | `t_mat` | `assemble_2c` | `doscentros(int=13)` |
| $V_{na}$ | Neutral Atom Potential | `vna` | `assemble_2c`, `assemble_3c` | `doscentros(int=2,3,4)`, `trescentros(int=1)` |
| $V_{nl}$ | Non-local Pseudopotential | `vnl` | `assemble_2c_PP`, `assemble_3c_PP` | `doscentrosPP`, `trescentrosPP` |
| $V_{xc}$ | Multi-center XC (Off-site) | `vxc` | `assemble_olsxc_off` | `doscentros(int=6)` |
| $V_{xc\_1c}$ | One-center XC (On-site) | `vxc_1c` | `assemble_olsxc_1c` | `assemble_olsxc_on` |
| $V_{ca}$ | Charge-dependent Potential | `vca` | `assemble_ca_2c`, `assemble_ca_3c` | DOGS-specific kernels |
| $V_{xc\_ca}$ | Charge-dependent XC | `vxc_ca` | `assemble_olsxc_off` | DOGS-specific kernels |
| $V_{Ewald}$ | Long-range Electrostatics | `ewaldlr/sr` | `assemble_lr` | Ewald summation |

---

## 2. Assembly Workflow (Pseudocode)

The assembly is managed by `assemble_mcweda.f90`. It distinguishes between geometry-dependent terms (calculated once at $K_{scf}=1$) and charge-dependent terms (recalculated every SCF step).

### High-Level Algorithm

```python
# PRE-CALCULATION (Kscf == 1)
# ---------------------------------------------------------------------------
# These terms depend ONLY on geometry and neutral atom densities.
# They are calculated once and stored in H and S sparse blocks.

for i in atoms:
    # 2-CENTER TERMS (T, Vna, Vnl, S)
    for j in neighbors[i]:
        r_ij = distance(i, j)
        orient = orientation(i, j)
        
        # Kinetic and Overlap
        T[i,j], S[i,j] = doscentros(interaction=[13, 1], r_ij, orient)
        
        # Neutral Atom (2-center parts: on-top and atom)
        Vna[i,j] += doscentros(interaction=[2, 3, 4], r_ij, orient)

    # 3-CENTER TERMS (Vna - "bcna")
    for (j, k) in common_neighbor_pairs[i]:
        # Potential centered at 'i', orbitals at 'j' and 'k'
        x, y, cost = geometry_3c(i, j, k)
        # interaction=1: bcna (Neutral Atom part)
        Vna[j,k] += trescentros(interaction=1, x, y, cost, orient)

# SCF LOOP (Every Step)
# ---------------------------------------------------------------------------
# These terms depend on the self-consistent charges Qin[i].

for i in atoms:
    # 1. DENSITY AVERAGING (average_ca_rho)
    # 3-center density overlap is used here!
    for (j, k) in common_neighbor_pairs[i]:
        x, y, cost = geometry_3c(i, j, k)
        # interaction=3: den3 (Density Overlap)
        # Note: The integral is fixed, but it's scaled by current Qin[i]
        rho_bar[j,k] += trescentros(interaction=3, x, y, cost) * Qin[i]

    # 2. OFF-SITE XC (assemble_olsxc_off)
    for j in neighbors[i]:
        # OLSXC formula application (Ref: PRB 40, 3979 (1989))
        # This uses rho_bar (averaged density) in a 2-center-like potential form.
        Vxc[i,j] = build_olsxc_off(Vxc_2c[i,j], Overlap[i,j], rho_bar[i,j], Qin[i], Qin[j])

# FINAL SUMMATION
H = T + Vna + Vnl + Vxc + Vxc_1c + Vca + Vxc_ca + V_Ewald
```

---

## 3. Data Flow & Interpolation

FireBall avoids on-the-fly integral calculation by reading pre-tabulated data (`Fdata`).

### Data Structures
-   **`Fdata3c.f90`**: Stores 2D grids for 3-center terms (`bcna_01` to `bcna_05`).
-   **`integrals.f90`**: Stores 1D splines for 2-center terms (`splineint_2c`).

### Interpolation Kernels
-   **`interpolate_1d`**: Uses a "superspline" (cubic spline) implementation.
    ```fortran
    ! aaa, bbb, ccc, ddd are spline coefficients
    dfdx = bbb + dx*(ccc*2 + dx*ddd*3)
    yout = aaa + dx*(bbb + dx*(ccc + dx*ddd))
    ```
-   **`interpolate_2d`**: Uses 4x4 sub-grid interpolation for 3-center terms.
    -   Interpolates over $x$ (distance to NA) and $y$ (bond length).
    -   Legendre polynomials are applied in `trescentros` to handle angular dependence ($\cos \theta$).

---

## 4. McWEDA / OLSXC Breakdown

The **McWEDA** (Multi-center Weighted Exchange-correlation Density Approximation) logic is implemented via **OLSXC** routines.

**Key Formula in `build_olsxc_off.f90`**:
Approximates $\langle \mu | V_{xc} | \nu \rangle$ using average density $\bar{\rho}$ and its derivative $V'_{xc}$:

$$ \langle i, \mu | V_{xc}(\rho) | j, \nu \rangle \approx V_{xc}(\bar{\rho}) S_{\mu \nu} + V'_{xc}(\bar{\rho}) [\langle \mu | \rho | \nu \rangle - \bar{\rho} S_{\mu \nu}] + \dots $$

This allows capturing the non-linearity of XC without a global real-space grid, maintaining the efficiency of a tight-binding-like scheme.
