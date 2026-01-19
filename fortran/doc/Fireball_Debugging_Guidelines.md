This document summarizes the methodology, technical insights, and strict guidelines derived from the porting and debugging of the Fireball $V_{CA}$ (Charged Atom) and Ewald summation terms from Fortran to OpenCL.

This guide is the reference for future work (specifically $V_{XC}$).

---

# I. The Core Methodology: Evidence-Based Porting

The primary source of friction in previous sessions was the attempt to "guess" the implementation based on physics intuition or standard formulas. **This is forbidden.** Fireball is a legacy code with specific, non-standard conventions (Ortega ordering, SFIRE symmetrization, specific rotation frames).

### The Golden Rule
**Fortran is the absolute Source of Truth.**
If the OpenCL result differs from Fortran, OpenCL is wrong, even if the Fortran implementation looks "weird" or non-standard.

### The Proven Debugging Workflow
Do not attempt to match the final result (e.g., the Hamiltonian matrix) in one go. Use the **Checkpoint Strategy**:

1.  **Micro-Verification (The Scan):**
    *   Before assembling neighbors, verify the spline interpolation and rotation of a single pair in isolation.
    *   Vary distance ($R$) and angles.
    *   *Tool:* `scanHamPiece` vs. OpenCL Scan Kernels.
2.  **Structural Integrity (Indices):**
    *   Verify neighbor lists match exactly (integer comparison).
    *   Verify `neigh_self` (diagonal) and `neigh_back` (reverse) mappings.
    *   *Note:* A mismatch here guarantees failure later. Do not look at floats until indices match.
3.  **Intermediate Scalar Injection (The "Gated Print"):**
    *   **Do not** build elaborate C-bindings/interfaces to export internal scalars. It creates boilerplate hell.
    *   **Do** inject temporary `print *, ...` statements directly into the Fortran loops (e.g., `assemble_ca_2c.f90`).
    *   **Gate them:** Use `if (idebugWrite > 0 .and. iatom==1 .and. jatom==2) ...` to restrict output to a single pair.
    *   Mirror these prints in the Python/OpenCL test script. Compare the logs line-by-line.
4.  **Component Isolation:**
    *   Split the term: $V_{CA} = V_{2c\_atom} + V_{2c\_ontop} + V_{3c}$.
    *   Verify $V_{2c\_atom}$ (diagonal) first.
    *   Verify $V_{2c\_ontop}$ (off-diagonal) next.
    *   Only then tackle 3-center terms.

---

# II. Technical Insights: Fireball Internals

### 1. Orbital Ordering & Signs
Fireball uses the **Ortega Ordering**, not the standard quantum chemistry order.
*   **Order:** $s, p_y, p_z, p_x$ (Indices: 0, 1, 2, 3).
*   **The "Ontop" Sign Trap:** In 2-center "ontop" integrals (off-diagonal blocks in VCA), the coupling between $s$ and $p_z$ involves a specific sign convention.
    *   Mapping: $comp[1] \to sm[2][0]$ and $comp[2] \to sm[0][2]$.
    *   **Crucial:** These specific terms require a sign flip (negation) in the OpenCL kernel to match Fortran.

### 2. Spline Types (SK vs. Generic)
Not all tables are Slater-Koster (SK) parametrizable.
*   **S, T, Vna:** These are SK-encoded (ss, sp, ps, pp$\pi$, pp$\sigma$). Use the SK reconstruction path.
*   **Dipoles:** These are **not** SK-encoded. They require a "generic" path (`recover_2c` in Fortran).
    *   *Implementation:* OpenCL must detect dipole tables and switch from SK-math to raw $\mu,\nu$ map interpolation.

### 3. Rotations & Frames
*   **$\epsilon$ (Epsilon) Construction:** Rotation matrices are built from the bond vector.
*   **Inconsistency:** Different terms use different definitions of the rotation frame.
    *   Most terms use normalized $\hat{r}$.
    *   $V_{na}$ (atom part) and Dipoles sometimes calculate $\epsilon$ based on absolute $|r|^2$ or specific definitions in `doscentros`.
    *   *Action:* Always check `assemble_*.f90` to see exactly what is passed to `epsilon_fb` or `twister`.

### 4. SFIRE Symmetrization (The 3-Center Overwrite)
Fortran's 3-center assembly is **not purely additive**.
*   It computes the term for the forward neighbor $(i, ineigh)$.
*   It **overwrites** the reverse neighbor slot $(j, jneigh)$ with the transpose/conjugate of the forward term.
*   *Consequence:* You cannot simply sum $2c + 3c$ arrays to check results. You must replicate this overwrite logic in Python verification scripts to match the Fortran output.

### 5. Neighbor Identification
*   A neighbor is **not** defined just by $(i_{atom}, j_{atom})$.
*   It is defined by the triplet **$(i_{atom}, i_{neigh}, m_{beta})$**, where $m_{beta}$ handles periodic boundary conditions (cell shifts).
*   Verification scripts must respect $m_{beta}$; otherwise, periodic images collapse into one, causing errors.

---

# III. Pitfalls & Anti-Patterns (What Makes the User Angry)

1.  **Boilerplate Explosion:**
    *   *Bad:* Creating new `module_*.f90` files, reallocating arrays, and adding C-bindings just to check one scalar value (like `dq`).
    *   *Good:* Inserting a gated `print` inside the loop and reading `stdout`.

2.  **Modifying Logic for Debugging:**
    *   *Bad:* Hardcoding `k=0` or changing loop bounds in Fortran to isolate a case.
    *   *Good:* Using `if` guards (`if (iatom==1)`) so the physics logic remains untouched.

3.  **Ignoring Prerequisites:**
    *   *Bad:* Trying to debug EwaldLR when EwaldSR (Short Range) is still failing.
    *   *Bad:* Trying to debug Hamiltonian assembly when the Dipole or Overlap matrices don't match yet.

4.  **"Fishy" Zeros:**
    *   *Observation:* If a check returns `err=0.00e+00`, be suspicious. It usually means you are comparing a value to itself (tautology) or checking integers/structure, not float precision.
    *   *Action:* Clearly label structural checks vs. numerical checks in the summary output.

5.  **Assuming Standard Physics:**
    *   *Bad:* "The formula for Dipole rotation is usually X."
    *   *Good:* "Line 405 of `doscentros.f90` does Y."

---

# IV. Strategy for Future Terms ($V_{XC}$ & EwaldLR)

For the upcoming work (Exchange Correlation and Long Range Ewald), follow this specific sequence:

### Step 1: Input Verification
$V_{XC}$ depends on Density ($\rho$). EwaldLR depends on partial charges ($dq$).
*   **Action:** Ensure `AvgRho` (Average Density) and `Qin` (Charges) match Fortran 1:1.
*   **Check:** Use the "Qin-weighted" comparison if `itheory=1`.

### Step 2: Fortran Tracing (The "Gated Print")
Before writing a line of OpenCL:
1.  Run Fortran with a small system (e.g., C3).
2.  Locate the exact loop in `assemble_xc.f90` (or equivalent).
3.  Insert gated prints for **Inputs** (density values on grid/matrix) and **Outputs** (matrix elements before accumulation).

### Step 3: Minimal OpenCL Kernel
1.  Write a kernel that only reproduces the "Inputs" from Step 2. Verify.
2.  Implement the math.
3.  Verify the "Outputs" against Step 2.

### Step 4: Integration
1.  Sum into the global Hamiltonian.
2.  Check for SFIRE-like overwrites or self-interaction terms (diagonal) that behave differently.

### Summary of `verify_*.py` Best Practices
*   **Keep structural checks:** `Step1_neigh` and `Step4_BlockRecon` are essential sanity checks.
*   **Separate 4D vs Dense:** Debug in 4D (sparse blocks) first. Convert to Dense only for final full-matrix parity.
*   **Summary Table:** Maintain the high-level Pass/Fail table at the end of the script. It is the only way to track regression.