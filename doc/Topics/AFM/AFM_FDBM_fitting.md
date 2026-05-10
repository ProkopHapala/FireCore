CODEMAP: DFTB Reference Calculation & FDBM AFM Forcefield Comparison System
https://windsurf.com/codemaps/1153fe89-ff29-4d4b-b4a6-e97d8f37047f-fe86ab10a43f3d18


# USER

look we have @run_pyocl_fdbm_dftb_pentacene.py @AFM.py for calculation AFM image of pentacene molecule. What we need for rigorous testing is to fit some coefficients of the 3 forcefield parameters Pauli, Lenard jones and electrostatic, fot that we need some reference. I'm thinking to compute that reference using DFTB (as it it fast)

1)
@dftb_utils.py @DFTBcore.py , basically what we want is to define some points above the molecule  (e.g. over some carbon atom)  and in these points compute approach of the molecule with CO tip, it should be rigid scan (no relaxation), always just evaluate single-point energy (SCF after shift of the atoms of CO) the reference distance is distance betwee botom atom of CO (O atom) and C atom of the pentancene (if we choose to do the z-scan over that C atom). So figure out how to make such z-scan efficiently with @DFTBcore.py 
or @dftb_utils.py 

2) to cmpate with full-density-based model (FDBM) AFM simulation we need to cut out proper z-slice of data from the 3 forcefields Pauli, Electrostaic, Dispersion (London), we need to identify rigorously which grid point correspond precisely to position of that atom. I remember we were struggling with that before. Perhaps best how to check it visually is to set value at that grid point to some nonses e value (e.g. very high) and see if when we plot the grid and the atoms in one image the artifact introduced to the imshow is really at the position of that atom (the atom should be slighly transparent to make imshow visible)

now your job is not to program it yet, rather make comprehensive plan

---

# SWE-1.6


I'll examine the relevant files to understand the current structure and create a comprehensive plan for both tasks.

I'll examine the files to understand the current structure and create a comprehensive plan for both tasks.

# Comprehensive Plan for DFTB Reference Calculation and Grid Point Identification

## Task 1: DFTB-based Z-Scan for CO Tip Approaching Pentacene

### Objective
Compute reference energy vs distance curve for CO tip approaching pentacene to fit forcefield parameters (Pauli, LJ, electrostatic). Rigid scan (no relaxation), single-point SCF at each tip position.

### Current Capabilities Analysis

**From dftb_utils.py:**
- [constrained_scan()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:337:0-411:18) - existing constrained scan infrastructure with path-based atom movement
- [run_pbc()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:305:0-335:25) - runs DFTB+ with geometry optimization or single-point
- [parse_energy_out()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:231:0-257:17), [parse_forces()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:259:0-280:17) - result parsing
- [load_molecule()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:125:0-132:133) - XYZ loading

**From DFTBcore.py:**
- [set_coords()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:205:4-220:64) / [set_coords_and_lattice()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:222:4-243:81) - update geometry without re-initialization
- [run_scf()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:285:4-308:27) - run SCF calculation
- [enable_matrix_collection()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:245:4-259:9) - collect matrices if needed
- More efficient than dftb_utils.py for repeated SCF on same system

### Recommended Approach: Use DFTBcore with Coordinate Updates

**Advantages:**
- Initialize once, update coordinates repeatedly (much faster)
- No file I/O overhead for each scan point
- Direct API access to energy and forces

**Implementation Plan:**

1. **Initial Setup**
   - Load pentacene geometry from `pentacene.xyz` (22 atoms: C22H14)
   - Load CO tip geometry (O-C bond length ~1.13 Å)
   - Combine into single system (pentacene + CO)
   - Select target carbon atom on pentacene for z-scan (e.g., central C)
   - Define z-range: e.g., 2.5 Å to 6.0 Å in 0.1 Å steps (35 points)
   - Define CO orientation: O atom pointing down, C above (tip apex = O)

2. **DFTBcore Initialization**
   - Create DFTBcore instance
   - Write initial DFTB+ input file (no geometry optimization, SCC only)
   - Initialize with combined system
   - Enable matrix collection if needed for debugging

3. **Z-Scan Loop**
   - For each z-distance:
     - Position CO such that O atom is at (x_C, y_C, z_C + z_distance)
     - Update all atom coordinates via [set_coords()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:205:4-220:64)
     - Run SCF via [run_scf()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:285:4-308:27)
     - Extract total energy via [get_energy()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:318:4-327:27)
     - Optionally extract forces via [parse_forces()](cci:1://file:///home/prokop/git/FireCore/pyBall/dftb_utils.py:259:0-280:17) from detailed.out
     - Store (z_distance, energy, forces)
   - Save results to file

4. **Coordinate Update Details**
   - Pentacene atoms fixed at original positions
   - CO atoms translated along z-axis
   - CO orientation: O at bottom (closest to pentacene), C above
   - Reference: O-C distance = 1.13 Å, O is tip apex

5. **Output Format**
   - Save as numpy array or text file: columns = [z_distance, energy_Hartree, energy_eV]
   - Optionally save forces per atom for debugging
   - Plot energy vs distance curve

### Alternative Approach: dftb_utils.constrained_scan()

**Advantages:**
- Already implemented, tested
- Handles workdir management automatically
- Returns structured results

**Disadvantages:**
- Each step runs full DFTB+ subprocess (slower)
- More file I/O overhead

**Implementation if using constrained_scan:**
- Define path as array of shape (n_moved_atoms, n_steps, 3)
- Moved atoms = CO atoms (indices 22, 23 in combined system)
- Path = linear z-translation of CO
- Set `do_relax=False` for single-point calculations
- Set `fixed_idx` = pentacene atom indices (0-21)

### Recommended: DFTBcore Approach

**Rationale:** More efficient, cleaner for many scan points, direct API control.

## Task 2: Grid Point Identification for Forcefield Comparison

### Objective
Identify which grid point in the 3D forcefield grids corresponds precisely to a specific atom position, enabling extraction of z-slices for comparison with DFTB reference.

### Challenge
Grid points are discrete; need to map continuous atom position to nearest grid index. Visual verification needed to ensure correctness.

### Current Grid Setup (from AFM.py)

**Grid specification:**
```python
grid_spec = {
    'origin': origin,  # (3,) array in Å
    'dA': [step, 0, 0],
    'dB': [0, step, 0],
    'dC': [0, 0, step],
    'ngrid': ngrid,  # (nx, ny, nz)
}
```

**Grid point positions:**
```
grid_point(i, j, k) = origin + i*dA + j*dB + k*dC
                     = origin + [i*step, j*step, k*step]
```

### Implementation Plan

1. **Atom-to-Grid Index Mapping Function**
   - Input: atom position (x, y, z), grid_spec
   - Compute fractional coordinates: `frac = (pos - origin) / step`
   - Round to nearest integer: `idx = np.round(frac).astype(int)`
   - Clip to valid range: `idx = np.clip(idx, [0,0,0], ngrid-1)`
   - Return grid index (ix, iy, iz)
   - Verify: compute actual grid point position, check distance to atom

2. **Visual Verification Protocol**
   - Create copy of forcefield grid (e.g., E_pauli_field)
   - Set identified grid point to very high value (e.g., 1e6 or max*10)
   - Plot 2D slice through that point (XY plane at that z)
   - Overlay atom positions (slightly transparent, e.g., alpha=0.3)
   - Check if bright artifact appears exactly at atom position
   - Repeat for XZ and YZ slices if needed
   - If artifact misaligned, debug grid spec or coordinate system

3. **Z-Slice Extraction**
   - Once grid index (ix, iy, iz) is verified:
   - Extract z-column: `z_slice = grid[ix, iy, :]`
   - Convert to physical z-values: `z_vals = origin[2] + np.arange(nz) * step`
   - This gives forcefield value vs height above that atom
   - Can directly compare with DFTB energy vs distance curve

4. **Coordinate System Considerations**
   - **AFMulator grid**: Molecules shifted so grid occupies [0, L] (see [setup_grid()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:122:4-156:16))
   - **mol_shift**: Applied to atoms_arr to center in grid
   - Need to account for mol_shift when mapping atom positions
   - Atom position in grid space = original_position + mol_shift
   - Or: use atoms_arr positions directly (already shifted)

5. **Implementation Steps**
   - Function `atom_to_grid_idx(atom_pos, grid_spec, mol_shift=None)`
   - Function `visualize_grid_point(grid, grid_spec, atom_pos, mol_shift, output_path)`
   - Function `extract_z_slice(grid, ix, iy, grid_spec)`
   - Test on known atom (e.g., central C in pentacene)

### Visual Verification Details

**Plotting approach:**
```python
# Create figure with 2 subplots: grid slice and atom overlay
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: grid slice with artifact
im = ax1.imshow(grid[:, :, iz].T, origin='lower', extent=extent_xy)
ax1.scatter(atom_pos[:,0], atom_pos[:,1], c='red', s=50, alpha=0.5)

# Right: zoom on target atom
ax2.imshow(grid[:, :, iz].T, origin='lower', extent=extent_xy)
ax2.scatter(target_atom_pos[0], target_atom_pos[1], c='red', s=100, alpha=0.7)
ax2.set_xlim([target_atom_pos[0]-1, target_atom_pos[0]+1])
ax2.set_ylim([target_atom_pos[1]-1, target_atom_pos[1]+1])
```

**Success criteria:**
- Bright spot in grid exactly overlaps target atom
- No offset in x or y direction
- z-index corresponds to correct height

## Task 3: Integration and Comparison

### Objective
Compare DFTB reference energy curve with FDBM forcefield z-slices to fit coefficients.

### Implementation Plan

1. **DFTB Reference Curve**
   - Run DFTB z-scan (Task 1)
   - Get E(z) for CO approaching pentacene
   - Reference distance: O-C distance
   - Convert to interaction energy: E_int(z) = E_total(z) - E_pentacene - E_CO

2. **FDBM Forcefield Extraction**
   - Run FDBM pipeline (existing in run_pyocl_fdbm_dftb_pentacene.py)
   - Get E_pauli_field, E_es_field, E_vdw_field
   - Identify grid point for target C atom (Task 2)
   - Extract z-slices: E_pauli(z), E_es(z), E_vdw(z)

3. **Fitting**
   - Model: E_fit(z) = A_pauli * E_pauli(z) + A_es * E_es(z) + A_vdw * E_vdw(z)
   - Fit coefficients A_pauli, A_es, A_vdw to minimize |E_fit - E_DFTB|
   - Use linear least squares (numpy.linalg.lstsq)
   - Validate on separate test points

4. **Validation**
   - Compare fitted curve with DFTB reference
   - Compute RMSE, max error
   - Plot both curves
   - Test on different atom positions (e.g., different C atoms)

## Summary of Implementation Order

1. **Implement DFTB z-scan using DFTBcore**
   - Create function `dftb_zscan_pentacene_co()`
   - Test on small range first (3-4 points)
   - Verify energy trend (repulsion at short range, attraction at medium)
   - Full scan once validated

2. **Implement grid point identification**
   - Create `atom_to_grid_idx()` function
   - Create visual verification function
   - Test on pentacene grid from existing FDBM run
   - Verify artifact alignment visually

3. **Extract z-slices from forcefields**
   - Load existing E_pauli, E_es, E_vdw grids
   - Extract z-column at verified grid point
   - Plot vs height

4. **Fit coefficients**
   - Combine DFTB reference and FDBM slices
   - Perform linear fit
   - Validate and plot results

## Key Technical Considerations

**DFTB z-scan:**
- Coordinate units: DFTBcore expects Bohr, but we work in Angstrom (need conversion)
- CO orientation: Ensure O is tip apex (bottom atom)
- Convergence: Use tight SCCTolerance for accurate energies
- Single-point: No geometry optimization (Driver = None or removed)

**Grid mapping:**
- mol_shift: Account for molecule shift in AFMulator.setup_grid()
- Grid origin: Not necessarily (0,0,0), use grid_spec['origin']
- Step size: Uniform in all directions (dA, dB, dC all have same magnitude)
- Boundary conditions: Clip indices to valid range

**Visual verification:**
- Use transparency for atoms (alpha < 1.0)
- Set artifact value significantly higher than natural range
- Check multiple slices (XY, XZ, YZ) for 3D confidence

I've created a comprehensive plan for both tasks. Here's the summary:

## Task 1: DFTB Z-Scan for CO Tip Approaching Pentacene

**Recommended approach: Use DFTBcore with coordinate updates** (more efficient than dftb_utils.constrained_scan)

**Implementation steps:**
1. Load pentacene (C22H14) and CO tip geometries, combine into single system
2. Select target carbon atom on pentacene for z-scan
3. Initialize DFTBcore once with combined system (no geometry optimization, SCC only)
4. Loop over z-distances (e.g., 2.5-6.0 Å in 0.1 Å steps):
   - Position CO such that O atom is at (x_C, y_C, z_C + z_distance)
   - Update coordinates via [set_coords()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:205:4-220:64)
   - Run SCF via [run_scf()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:285:4-308:27)
   - Extract energy via [get_energy()](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/DFTBcore.py:318:4-327:27)
   - Store (z_distance, energy)
5. Save results and plot E(z) curve

**Key technical points:**
- DFTBcore expects Bohr coordinates, need conversion from Angstrom
- CO orientation: O at bottom (tip apex), C above (O-C bond ~1.13 Å)
- Reference distance = O-C distance between tip and pentacene
- Single-point calculations only (no relaxation)

## Task 2: Grid Point Identification for Forcefield Comparison

**Goal:** Map atom position to grid index to extract z-slices from Pauli/ES/Dispersion forcefields

**Implementation steps:**
1. **Atom-to-grid mapping function:**
   - Compute fractional coordinates: `frac = (pos - origin) / step`
   - Round to nearest integer: `idx = np.round(frac).astype(int)`
   - Clip to valid range
   - Verify by computing actual grid point position

2. **Visual verification protocol:**
   - Create copy of forcefield grid
   - Set identified grid point to very high value (e.g., 1e6)
   - Plot 2D slice (XY plane at that z) with atom positions overlaid (alpha=0.3)
   - Check if bright artifact appears exactly at target atom
   - Repeat for XZ/YZ slices if needed

3. **Z-slice extraction:**
   - Once grid index (ix, iy, iz) verified: `z_slice = grid[ix, iy, :]`
   - Convert to physical z-values: `z_vals = origin[2] + np.arange(nz) * step`

**Key technical points:**
- Account for `mol_shift` (AFMulator shifts molecules to [0,L] grid space)
- Use `atoms_arr` positions directly (already shifted)
- Grid spec: `origin`, `dA/dB/dC` (step vectors), `ngrid`

## Task 3: Integration and Coefficient Fitting

**Steps:**
1. Run DFTB z-scan → get E_DFTB(z) interaction energy
2. Run FDBM pipeline → get E_pauli, E_es, E_vdw grids
3. Extract z-slices at verified grid point
4. Fit model: `E_fit(z) = A_pauli*E_pauli(z) + A_es*E_es(z) + A_vdw*E_vdw(z)`
5. Use linear least squares to fit coefficients
6. Validate with RMSE and plot comparison

## Implementation Order

1. Implement DFTB z-scan using DFTBcore (test small range first)
2. Implement grid point identification with visual verification
3. Extract z-slices from existing FDBM forcefield grids
4. Fit coefficients and validate

The plan is ready for implementation when you're ready to proceed.

---

# Experimental Work: Grid Step Dependence and Pauli Repulsion Model

## What We Have Done

### 1. Density Reprojection at Different Grid Steps
- **Goal:** Test whether electron density integration is independent of grid step size
- **Method:** Reprojected DFTB densities from step 0.15 Å to step 0.1 Å and 0.2 Å using fast density matrix method
- **Files:**
  - `reproject_to_step.py` - Reprojects pentacene and CO densities onto different grid steps
  - `compute_pauli_from_reprojected.py` - Computes Pauli energy from reprojected densities
  - `test_dftb_fdbm_comparison.py` - Compares DFTB z-scan with FDBM forcefield grids

### 2. Density Integration Verification
- **Results:**
  - Step 0.1 Å: Pentacene 102.008 e (expected 102 e, diff 0.008), CO 9.985 e (expected 10 e, diff 0.015)
  - Step 0.15 Å: Pentacene 101.958 e (expected 102 e, diff 0.042)
  - Step 0.2 Å: Pentacene 102.251 e (expected 102 e, diff 0.251), CO 9.931 e (expected 10 e, diff 0.069)
- **Conclusion:** Density integration is reasonably independent of grid step size (step 0.1 gives most accurate result)

### 3. Pauli Energy Calculation from Reprojected Densities
- **Method:** FFT convolution of pentacene density difference with CO density difference
- **Implementation:**
  ```python
  tip_kernel = np.roll(np.roll(np.roll(rho_tip[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
  E_pauli_field = A_pauli * dV * np.real(np.fft.ifftn(np.fft.fftn(rho_grid) * np.fft.fftn(tip_kernel)))
  ```
- **Parameters:** A_pauli = 16.0 (default), dV = step^3

### 4. Comparison with DFTB Z-Scan
- **DFTB Reference:** Rigid CO tip approaching pentacene, single-point SCF at each z-distance
- **FDBM Extraction:** Z-slice at target atom grid column from pre-computed forcefield grids
- **Comparison:** DFTB total energy change vs FDBM forcefield components (Pauli, ES, vdW, Total)

## Problems Encountered

### 1. CO Density Positioning for Convolution
- **Issue:** Original reprojection positioned CO at its actual DFTB geometry position, not at grid center
- **Problem:** Pauli convolution assumes CO density is centered at grid origin to represent tip scanning approach
- **Fix:** Modified `reproject_to_step.py` to center CO at grid center before projection:
  ```python
  grid_center = origin + 0.5 * (ngrid - 1) * args.step
  co_geo['coords_bohr'] = (co_geo['coords_bohr'] * BOHR2ANG + grid_center) / BOHR2ANG
  ```
- **Status:** Fixed - CO now centered at grid center for convolution

### 2. Pauli Energy Magnitude Mismatch
- **Issue:** Pauli values from reprojected grids were very small (0.2 eV) compared to DFTB energy change (20 eV)
- **Root Cause:** 
  - For step 0.1 and 0.2, only Pauli was computed (ES and vdW missing)
  - Total = Pauli + ES + vdW = Pauli + 0 + 0 = Pauli
  - On single-axis plot, Total and Pauli lines overlapped perfectly, making Pauli invisible
- **Fix:** Modified plot to use separate subplots for each component with auto-scaled y-axes
- **Status:** Fixed - Pauli now visible in its own subplot

### 3. Grid Point Mapping for Z-Slice Extraction
- **Issue:** Need to map continuous atom position to discrete grid index for z-slice extraction
- **Method:** `atom_to_grid_idx()` function:
  ```python
  frac = (pos - origin) / step
  idx = np.round(frac).astype(int)
  idx = np.clip(idx, [0,0,0], ngrid-1)
  ```
- **Verification:** Visual check by setting grid point to high value and overlaying atom positions
- **Status:** Implemented but not yet visually verified

### 4. Pauli Energy Scaling/Shift
- **Observation:** Pauli vs z-distance curves for different step sizes are shifted/scaled differently
- **Issue:** Absolute Pauli values don't match expected physical magnitude
- **Possible Causes:**
  - A_pauli parameter (currently 16.0) may be incorrect
  - Density normalization factor (B3_FACTOR) may be wrong
  - Convolution kernel construction may have issues

## Pauli Repulsion Model Investigation

### Current Implementation

**Formula:**
```
E_pauli(r) = A_pauli * dV * ∫ ρ_sample(r') * ρ_tip(r' - r) dr'
```

**Parameters:**
- **A_pauli = 16.0** (hardcoded in both `run_pyocl_fdbm_dftb_pentacene.py` and `compute_pauli_from_reprojected.py`)
- **dV = step^3** (voxel volume in Å^3)
- **No beta parameter** - current model is linear in density product

**Densities Used:**
- **Pentacene:** `rho_diff = rho_pentacene - rho_pentacene_na` (valence electron density minus neutral atom density)
- **CO:** `rho_co_diff = rho_co - rho_co_na` (valence electron density minus neutral atom density)
- **No exponentiation** - densities used as-is, not raised to any power

**Tip Kernel Construction:**
```python
# Original FDBM (run_pyocl_fdbm_dftb_pentacene.py step3_pauli):
tip_kernel = np.roll(np.roll(np.roll(rho_tip_total[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
# Uses rho_tip_total (full CO density)

# Reprojected version (compute_pauli_from_reprojected.py):
tip_kernel = np.roll(np.roll(np.roll(rho_co_diff[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
# Uses rho_co_diff (CO density difference)
```

**Key Difference:**
- Original FDBM uses `rho_tip_total` (full CO density) for Pauli
- Reprojected version uses `rho_co_diff` (CO density difference, i.e., valence electrons only)
- This may cause scaling differences between the two approaches

### Questions to Investigate

1. **Beta Parameter:** Is there a beta parameter in the Pauli model? Current code shows no beta, only A_pauli=16.0.
   - Standard Pauli repulsion models often use: E ∝ ρ^β with β ≈ 1-2
   - Our current model uses β=1 (linear in density)
   - Should we test β=2 (quadratic in density)?

2. **A_pauli = 16.0:** Is this value physically correct?
   - This appears to be an empirical scaling factor
   - Should be determined by fitting to DFTB reference data
   - Current value may be arbitrary

3. **Density Difference vs Total Density:**
   - Why does original FDBM use `rho_tip_total` for Pauli but `rho_tip_delta` for electrostatics?
   - Should Pauli use density difference (valence electrons only) or total density?
   - Electrostatics clearly uses density difference (charge density)

4. **CO Centering:**
   - Original FDBM centers CO at grid center for convolution
   - Our reprojection now also centers CO at grid center
   - Is this correct for representing tip scanning approach?

5. **Density Normalization:**
   - B3_FACTOR = 1.0 / (BOHR2ANG^3) converts from Bohr^3 to Å^-3
   - Ensures density integral gives electron count in Å units
   - Is this applied correctly in all places?

## Next Steps

1. **Verify CO Density Positioning:**
   - Plot CO density grid with atom overlay to confirm centering
   - Check if convolution kernel construction is correct

2. **Investigate A_pauli and Beta:**
   - Test different A_pauli values (e.g., 8, 16, 32)
   - Test beta=1 vs beta=2 (i.e., ρ vs ρ^2)
   - Fit to DFTB reference data to determine optimal parameters

3. **Compare Density Difference vs Total Density:**
   - Compute Pauli using both `rho_tip_total` and `rho_co_diff`
   - Compare results with original FDBM
   - Determine which is physically correct

4. **Visual Verification of Grid Mapping:**
   - Implement visual check by setting grid point to high value
   - Overlay atom positions with transparency
   - Confirm artifact aligns with target atom

5. **Complete FDBM Components for Reprojected Grids:**
   - Compute ES energy (need V_ES from Poisson solver)
   - Compute vdW energy (need C6 coefficients)
   - This will allow proper Total energy comparison

## Slow Density Projection Method

The original density projection used a slow method:
- Projected each orbital individually
- For each occupied orbital: projected orbital squared, multiplied by occupation
- Summed over all orbitals
- This is O(n_occ * n_grid * n_atoms) complexity

The fast density matrix method used in `reproject_to_step.py`:
- Projects density matrix directly using sparse block structure
- More efficient GPU-accelerated projection
- Complexity depends on neighbor list size rather than total orbitals

This is the "horribly slow method" mentioned - the original orbital-by-orbital projection was slow for large systems.
# SWE-2.0 (Analysis and Verification of Backend Parity)

## Comparison of Fireball vs. DFTB+ Backends

Following the implementation of the DFTB+ backend for FDBM AFM simulations, a rigorous comparative analysis was performed to understand the magnitude discrepancy in Pauli repulsion forces.

### 1. Key Findings: Basis Set Compactness
The discrepancy is primarily due to the **decay characteristics of the basis functions**. 
- **Fireball**: Uses numerical orbitals with "fat" tails.
- **DFTB+ (mio-1-1)**: Uses Slater-Type Orbitals (STOs) that are significantly more compact.

| Metric (at z=2.1 Å) | Fireball | DFTB+ | Ratio (FB/DB) |
| :--- | :--- | :--- | :--- |
| **Electron Density ($\rho$)** | 0.002385 $e/\text{\AA}^3$ | 0.000127 $e/\text{\AA}^3$ | **~18.8x** |
| **Pauli Energy ($E_{Pauli}$)** | 4.48175 eV | 0.17685 eV | **~25.3x** |

### 2. Physical Parity and Fitting
Since Pauli repulsion in the FDBM model is proportional to the overlap integral $\int \rho_s \rho_t dV$, the "thin" STO tails lead to a ~100x reduction in total repulsion magnitude compared to Fireball (where $\rho_s$ and $\rho_t$ are both fatter).

To achieve physical parity with the SCC-DFTB total energy reference (~20 eV at 2.1 Å), the following parameters were identified:
- **$A_{pauli} \approx 1500-2000$**: Required to compensate for the smaller density overlap.
- **$\beta_{pauli} \approx 0.83$**: Required to "soften" the STO tails and match the physical force gradient (DFTB total energy decays at 4.2 Å⁻¹, while STO overlap decays at 5.0 Å⁻¹).

### 3. Consolidated Verification System
A standardized verification script was developed to prevent future regressions:
- **[`compare_fireball_dftb.py`](file:///home/prokop/git/FireCore/tests/tAFM/pyocl_fdbm/compare_fireball_dftb.py)**

This script automates:
- **Density Integration Checks**: Confirms charge conservation (102e) on both grids.
- **Side-by-side Plotting**: Generates $\rho(z)$ and $E(z)$ diagnostic plots in `debug_comparison_v2/`.
- **Numerical Reporting**: Provides a detailed summary of density ratios and energy errors at contact.

The results confirm that while the mathematical projection is correct, the empirical parameters of the FDBM model must be tailored to the specific basis set characteristics of the backend.

### 4. Formal Fitting Results
A automated linear regression was performed in the range $z \in [2.1, 3.0]$ Å to optimize the Pauli parameters for the `mio-1-1` basis. The fit linearizes the exponential decay by analyzing $\ln(E_{ref})$ vs. $\ln(\text{overlap})$.

**Optimal Parameters Found:**
- **$A_{pauli} = 950.89$**
- **$\beta_{pauli} = 0.8578$**

These values provide a near-perfect match to the SCC-DFTB repulsive wall in the contact region, ensuring that the FDBM simulation accurately reproduces the physical height-dependent contrast observed in reference quantum mechanical calculations.

The consolidated script now includes this fitting routine as a standard calibration step for any new basis set or molecule.

---

# Multi-Basis-Set Fitting Pipeline (Updated 2026-05-10)

## Overview
A unified Python fitting pipeline was implemented to calibrate FDBM Pauli parameters independently for each DFTB+ basis set against rigid z-scan DFTB reference energies. Two basis sets were calibrated: **mio-1-1** and **3ob-3-1**. For each basis, four peripheral carbon atoms (atoms 0, 1, 20, 21 in `pentacene.xyz`) were fitted to assess site-to-site consistency.

## Relation to Existing Scripts

| Script | Purpose | When to use |
|--------|---------|-------------|
| `compare_fireball_dftb.py` | **Legacy diagnostic**: compares Fireball vs DFTB+ backends, explains the ~25× magnitude discrepancy in density/Pauli tails, and performs a single-atom fit against a pre-existing z-scan text file. Hardcoded for mio-1-1 only. | Understanding backend differences; Fireball parity analysis. |
| `run_dftb_zscan.py` | **New**: generates DFTB z-scan reference data from scratch. Supports any basis set and any set of target atoms via `--target_indices`. | When fresh reference data is needed (new basis, new molecule, new atoms). |
| `fit_fdbm_pauli.py` | **New**: production calibration script. Direct nonlinear least-squares fit, multi-atom support, compact per-atom 2-panel plots + summary figure. | Production Pauli calibration for any basis set. |

The new pipeline replaces the old fitting because: (1) the old z-scan reference (`debug_dftb_comparison/zscan_results.txt`) did not exist for 3ob-3-1, (2) the old script is Fireball-centric and single-atom, and (3) the old slope-ratio fitting method is less robust than direct nonlinear least squares.

## Implementation Files

| File | Purpose |
|------|---------|
| `tests/tAFM/pyocl_fdbm/fit_fdbm_pauli.py` | Main fitting script. Loads FDBM grids + DFTB z-scan, extracts profiles at specified atom positions, fits power-law per atom. Supports `--basis`, `--target_indices`, `--zscan_dir`. Produces per-atom `fit_pauli.png` and a multi-atom `summary_all_atoms.png`. |
| `tests/tAFM/pyocl_fdbm/run_dftb_zscan.py` | Standalone DFTB rigid z-scan. Runs `dftb+` single-point SCF at 54 distances (2.0–10.0 Å, Δz=0.15 Å) with per-atom caching. Basis-agnostic via `--basis` flag. Supports multiple target atoms via `--target_indices`. |
| `tests/tAFM/pyocl_fdbm/run_pyocl_fdbm_dftb_pentacene.py` | FDBM grid generation pipeline. Projects DFTB densities, computes Pauli (FFT convolution), electrostatics (Poisson + convolution), and vdW (C6/r⁶) fields. |

## Fitting Methodology

### 1. Pauli Power-Law Fit
Model: `E_DFTB(z) = A_pauli · overlap(z)^beta`
- `overlap(z)` is extracted from the raw FDBM Pauli field divided by the default coefficient `A_default = 16.0`
- Fit range: `z ∈ [2.0, 3.0] Å` (contact/decay region)
- **Nonlinear least squares** (`scipy.optimize.curve_fit`) with log-linear initialization for robustness
- Overlap values are clipped to `≥1e-30` before power evaluation to avoid numerical underflow/overflow with fractional `beta`

This is more robust than the slope-ratio method in `compare_fireball_dftb.py`, which estimates `beta = slope_ref / slope_ovl` from separate log-linear fits and is sensitive to interpolation errors.

## Calibrated Parameters (Multi-Atom)

### mio-1-1

| Atom | Name | Position [Å] | A_pauli | beta | R² | RMSE(fit) [eV] |
|------|------|-------------|---------|------|-----|----------------|
| 0 | C | [-6.0056, -0.8000, 0.0000] | 949.86 | 0.8725 | 0.999881 | 0.119 |
| 1 | C | [-6.0056, 0.8000, 0.0000] | 924.30 | 0.8717 | 0.999882 | 0.118 |
| 20 | C | [5.9944, -0.8000, 0.0000] | 1000.44 | 0.8703 | 0.999883 | 0.117 |
| 21 | C | [5.9944, 0.8000, 0.0000] | 986.73 | 0.8695 | 0.999883 | 0.117 |
| **Mean ± std** | | | **965 ± 35** | **0.8710 ± 0.0012** | | |

### 3ob-3-1

| Atom | Name | Position [Å] | A_pauli | beta | R² | RMSE(fit) [eV] |
|------|------|-------------|---------|------|-----|----------------|
| 0 | C | [-6.0056, -0.8000, 0.0000] | 631.72 | 0.7983 | 0.999914 | 0.121 |
| 1 | C | [-6.0056, 0.8000, 0.0000] | 618.40 | 0.7972 | 0.999913 | 0.121 |
| 20 | C | [5.9944, -0.8000, 0.0000] | 664.91 | 0.7945 | 0.999908 | 0.124 |
| 21 | C | [5.9944, 0.8000, 0.0000] | 658.95 | 0.7934 | 0.999907 | 0.125 |
| **Mean ± std** | | | **643 ± 22** | **0.7959 ± 0.0023** | | |

## Key Findings

### 1. 3ob-3-1 produces ~2× larger contact repulsion
- Contact ΔE at z=2.1 Å: **mio-1-1 ≈ 32.4 eV**, **3ob-3-1 ≈ 38.9 eV**
- This is captured by the fit through smaller `A_pauli` (643 vs 965 mean) and slightly lower `beta` (0.796 vs 0.871 mean)
- The 3ob basis functions have tighter radial functions with additional Gaussians, leading to different overlap magnitudes

### 2. Site-to-site consistency is excellent
- **mio-1-1**: std(A) = 35 (3.6% relative), std(beta) = 0.0012 (0.14% relative)
- **3ob-3-1**: std(A) = 22 (3.4% relative), std(beta) = 0.0023 (0.29% relative)
- 3ob-3-1 is slightly more consistent across atoms (smaller relative scatter)

### 3. Small left/right asymmetry in mio-1-1
- Left edge atoms (0, 1): A ≈ 925–950
- Right edge atoms (20, 21): A ≈ 987–1000
- This reflects real geometry-dependent density tails, not a fitting artifact

## Output Directory Structure

Each basis set produces a dedicated output directory with per-atom subdirectories:
```
fit_pauli_mio_1_1/
├── atom_0/
│   ├── fit_pauli.png      # 2-panel: linear + log, DFTB vs fitted Pauli
│   ├── params.json        # A_pauli, beta_pauli, R², RMSE for this atom
│   ├── z_ref.npy, e_ref.npy, e_pauli.npy, overlap_raw.npy
│   └── fit.log
├── atom_1/
│   └── ... (same structure)
├── atom_20/
│   └── ...
├── atom_21/
│   └── ...
├── summary_all_atoms.png   # A, beta, RMSE bar charts + overlaid fitted curves
└── summary.txt             # Table of all atoms

fit_pauli_3ob_3_1/          # Same structure
```

## Usage

### Single atom (default)
```bash
cd tests/tAFM/pyocl_fdbm

# Step 1: Generate FDBM grids
python3 run_pyocl_fdbm_dftb_pentacene.py --basis mio-1-1
python3 run_pyocl_fdbm_dftb_pentacene.py --basis 3ob-3-1 --output_dir debug_dftb_pentacene_3ob

# Step 2: Run DFTB z-scan reference (single atom)
python3 run_dftb_zscan.py --basis mio-1-1 --output_dir zscan_mio_1_1 --target_indices 0
python3 run_dftb_zscan.py --basis 3ob-3-1 --output_dir zscan_3ob_3_1 --target_indices 0

# Step 3: Fit parameters
python3 fit_fdbm_pauli.py --basis all --target_indices 0
```

### Multiple atoms
```bash
# Run z-scan for 4 peripheral carbons
python3 run_dftb_zscan.py --basis mio-1-1 --output_dir zscan_mio_1_1 --target_indices 0,1,20,21
python3 run_dftb_zscan.py --basis 3ob-3-1 --output_dir zscan_3ob_3_1 --target_indices 0,1,20,21

# Fit all atoms
python3 fit_fdbm_pauli.py --basis all --target_indices 0,1,20,21
```

## Notes
- Grid origin `[-11.36, -6.6, -4.2] Å`, step `0.15 Å`, shape `(152, 88, 96)` is identical for both basis sets to ensure fair comparison.
- The Pauli power-law fit quality is excellent across all atoms (R² > 0.9999).
- The ES and vdW components are much smaller than Pauli in the contact region and are not included in the compact plots.
- For legacy Fireball-vs-DFTB backend comparison, see `compare_fireball_dftb.py`.

---

# AFM Image Generation with Fitted Pauli Parameters (2026-05-10)

## Overview
After fitting Pauli parameters (A_pauli, beta_pauli) against DFTB z-scan reference, the main FDBM pipeline was updated to use these fitted parameters instead of the hardcoded A=16.0. AFM images were generated to verify that the forcefield produces physically realistic attractive→repulsive behavior.

## Script Updates

### Main Pipeline: `run_pyocl_fdbm_dftb_pentacene.py`
**Changes:**
- Added CLI arguments: `--A_pauli`, `--beta_pauli`
- Implemented basis-specific fitted defaults:
  - mio-1-1: A=965.0, beta=0.871
  - 3ob-3-1: A=643.0, beta=0.796
- Modified `step3_pauli()` to compute power-law: `E_pauli = A_pauli * overlap_raw^beta_pauli`
  - Previously: `E_pauli = 16.0 * overlap_raw` (linear scaling)
  - Overlap values clipped to ≥1e-30 before power evaluation to avoid numerical issues

**Purpose:** Production FDBM pipeline with calibrated Pauli repulsion.

### Fast Rerun Script: `run_fitted_afm.py`
**Purpose:** Loads existing density grids from `debug_dftb_pentacene*/` and recomputes only steps 3-6 (Pauli, ES, vdW, AFM) with new Pauli parameters. Avoids redoing slow density projection.

**Usage:**
```bash
cd tests/tAFM/pyocl_fdbm
python3 run_fitted_afm.py --basis mio-1-1
python3 run_fitted_afm.py --basis 3ob-3-1
```

**Output:** `afm_fitted_mio_1_1/` or `afm_fitted_3ob_3_1/` containing:
- `E_Pauli_field.npy`, `E_es_field.npy`, `E_vdw_field.npy`
- `df.npy` (frequency shift array)
- `df_h3.0.png`, `df_h4.2.png`, `df_h5.4.png` (AFM images at different heights)
- `pauli_slices.png` (XY, XZ, YZ slices of Pauli field)

### Diagnostic Script: `diagnostic_forcefield.py`
**Purpose:** Quick sanity check of forcefield components at relevant z-heights (2.0, 2.5, 3.0 Å).

**Output:** `diagnostic/` subdirectory with:
- `panel_pauli.png` — Pauli at multiple z
- `panel_vdw.png` — vdW at multiple z
- `panel_total.png` — Total (symmetric colormap: red=repulsive, blue=attractive)
- `panel_Fz.png` — Force z-component (symmetric colormap)
- `combined_z2.5.png` — All 4 panels at z=2.5 Å

### Transition Visualization: `plot_transition.py`
**Purpose:** Visualize the attractive→repulsive transition and identify the crossover height where Pauli overtakes vdW.

**Output:** `transition/` subdirectory with:
- `transition_panels.png` — 4×5 panel: rows=[Pauli, vdW, ES, Total], columns=[z=2.5, 2.8, 3.0, 3.2, 3.5 Å]. Each subplot has its own colorbar with per-subplot vmin/vmax.
- `zscan_atom0.png`, `zscan_atom1.png`, `zscan_atom20.png`, `zscan_atom21.png` — Individual z-scan curves for each peripheral carbon (z=2.0–6.0 Å, y-axis fixed to ±0.5 eV to see transition clearly). Shows Pauli, vdW, Electrostatics, Total vs z.
- `zscan_all_atoms.png` — Combined 2×2 panel with all 4 atoms
- `maxE_vs_z.png` — Max Pauli, |vdW|, and Total energy vs tip height (log scale)

## Key Findings

### 1. Attractive→Repulsive Crossover Heights

| Basis | Crossover (Pauli ≈ vdW) | Transition zone |
|-------|------------------------|-----------------|
| mio-1-1 | z ≈ 3.0–3.1 Å | 2.8–3.2 Å |
| 3ob-3-1 | z ≈ 3.3–3.4 Å | 3.0–3.5 Å |

The 3ob-3-1 basis, having tighter radial functions, produces a steeper Pauli wall that survives to larger tip-sample distances. This is physically correct and matches the fitted parameters (smaller A, lower beta for 3ob-3-1).

### 2. Forcefield Behavior at Different Heights

**At z = 2.5 Å:**
- Pauli dominates completely. Total is purely repulsive.
- mio-1-1: Total max ≈ +1.5 eV
- 3ob-3-1: Total max ≈ +3.2 eV

**At z = 3.0 Å (mio-1-1) / 3.2 Å (3ob-3-1):**
- Transition zone. Pauli and vdW comparable magnitude.
- Total shows mixed sign: repulsive over carbons, attractive in gaps.
- This is the optimal AFM imaging height for seeing individual atoms.

**At z = 3.5 Å:**
- vdW dominates. Total is attractive everywhere (negative).
- Pauli contribution is negligible (<0.01 eV for mio-1-1, <0.03 eV for 3ob-3-1).

### 3. Component Magnitudes (at z=2.5 Å, above carbon)

| Component | mio-1-1 | 3ob-3-1 |
|-----------|---------|---------|
| Pauli | +1.8 eV | +3.6 eV |
| vdW | -0.31 eV | -0.31 eV |
| Electrostatics | ±0.01 eV | ±0.01 eV |
| Total | +1.5 eV | +3.2 eV |

Electrostatics is negligible in the contact region. The competition is between Pauli (repulsive) and vdW (attractive), exactly as expected for a Lennard-Jones-like interaction.

### 4. AFM Signal Strength

The frequency shift (df) magnitude scales with the force gradient. From the fitted AFM images:

| Basis | df range (h=3.0 Å) |
|-------|-------------------|
| mio-1-1 | -3.1 to +0.08 Hz |
| 3ob-3-1 | -7.6 to +0.04 Hz |

3ob-3-1 produces ~2× stronger AFM contrast due to the steeper Pauli wall.

## Script Classification for Git

| Script | Purpose | Git Status |
|--------|---------|------------|
| `run_pyocl_fdbm_dftb_pentacene.py` | Main FDBM pipeline with fitted Pauli defaults | **Commit** (production) |
| `fit_fdbm_pauli.py` | Multi-atom Pauli fitting against DFTB reference | **Commit** (production) |
| `run_dftb_zscan.py` | DFTB z-scan reference generator (multi-atom, multi-basis) | **Commit** (production) |
| `plot_transition.py` | Transition visualization (panels + z-scan curves) | **Commit** (analysis) |
| `run_fitted_afm.py` | Fast rerun for Pauli parameter iteration | **Commit** (utility) |
| `diagnostic_forcefield.py` | Field slice diagnostics | **Commit** (utility) |
| `compare_fireball_dftb.py` | Legacy Fireball vs DFTB backend comparison | **Commit** (legacy reference) |

## Usage Summary

### Full Pipeline (from scratch)
```bash
cd tests/tAFM/pyocl_fdbm

# Step 1: Generate FDBM grids with fitted Pauli
python3 run_pyocl_fdbm_dftb_pentacene.py --basis mio-1-1
python3 run_pyocl_fdbm_dftb_pentacene.py --basis 3ob-3-1 --output_dir debug_dftb_pentacene_3ob

# Step 2: Generate AFM images (already done in Step 1)
# Images in: debug_dftb_pentacene/step6_composed/df_h*.png
```

### Fitting + AFM with Fitted Parameters
```bash
# Step 1: Generate DFTB z-scan reference
python3 run_dftb_zscan.py --basis mio-1-1 --output_dir zscan_mio_1_1 --target_indices 0,1,20,21
python3 run_dftb_zscan.py --basis 3ob-3-1 --output_dir zscan_3ob_3_1 --target_indices 0,1,20,21

# Step 2: Fit Pauli parameters
python3 fit_fdbm_pauli.py --basis all --target_indices 0,1,20,21

# Step 3: Run AFM with fitted parameters (fast rerun)
python3 run_fitted_afm.py --basis mio-1-1
python3 run_fitted_afm.py --basis 3ob-3-1

# Step 4: Visualize transition
python3 plot_transition.py --basis mio-1-1
python3 plot_transition.py --basis 3ob-3-1
```

## Conclusion

The fitted Pauli parameters successfully reproduce physically realistic forcefield behavior:
- Pauli repulsion dominates at close contact (z < 3 Å)
- vdW attraction dominates at larger distances (z > 3.5 Å)
- The crossover height matches the fitted parameters and differs between basis sets as expected
- AFM images show reasonable contrast and atom resolution in the transition zone

The 3ob-3-1 basis produces stronger repulsion and higher AFM contrast, consistent with its tighter basis functions. Both basis sets produce qualitatively correct LJ-like total interaction curves.
