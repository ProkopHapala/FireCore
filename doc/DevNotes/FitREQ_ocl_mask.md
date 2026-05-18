# USER

restricting to fragment 1 and fragment 2 is nonsense, this is interaction energy between fragment 1 and 2, there is no energy for fragments itself.

OK so now make a plan how to most effcient implement the masing into the energy kernel, and how to make the test, do not program yet.

By default my mask should split the energy to Epairs (interaction of Epairs with atoms) or Hcorrection (involving REQH.w component) and the rest. The must must be therefore float4 M[i,j]{x,y,z,w}
The idea behind is:
1) Basic non-covalent interaction potential in forcefield comprise Pauli+London (Morse or Lanard Jones or Buckingham, exp(), c6/r^6 ... this is made only of real atoms not atom pairs), and charges Qi*Qj/r (again only for real atoms), this is classical non covalent potentil - our baseline. Then we have two types of corrections. Homoganous spherical which apply H (REQH.w) on atoms which ither decreas Pauli on Increase Londondn dispersion. Then we have added electron pairs perticles which have exponential or other short range spehrical function, and these encode directionality (angular dependence as the electron pairs point in some specific direction). This what we want to study - we want to split the baseline (standard non covalent like UFF,MMFF) and our hydrogen bond corrections (H and epairs).
For this purpose we make the mask. Plese think it through and write notes here in chat  

we must analyse and decompose 
                //<<<MODEL_PAIR_ACCUMULATION
se we do not get only Eij but the components Pauli, London, Electronstics, and Hydrogen-bond separately 
We do not need to have splical slot for Eair and Homogenous H in the mask, because we can filter this by different i,j as Epairs are specific particle. Hope you understand. 

plese update the plan  

OK, add more notes that masking either homogenous H  or electronp pairs should be done on python site, that it set relevant atom indexes for the i,j matrix, and i is always from grament 1 and j alsoway from fragment 2. Also we sometimes will use it for masking the less relevant atom, for example in C2H=O (formaldehyde), the hydrogen bond involve the oxygen and tis electron pairs, the H atom should not participate, but it is interesting to know how much they play role. It would be good to have tools for that. 

Also note that we should always compute Etot (total energy) and Ein (Masked energy) in one pass of kernel, so we have two output arrays. So we can easily evaluate 
Etot = Ein + Eout and plot all 3 after one go

But for performance optimization this masking part should be grated by if (we will not read Mmask and  ouptut Ein if bMask=false), by default Mij=float(0,0,0,0) s
---

# SWE-1.6


## Implementation Plan for Energy Masking

### 1. Mask Design

**Format**: `float4 mask[ni, nj]` where `mask[i,j]` controls energy contribution for atom pair (i,j)

**Component mapping** (general):
- `mask.x`: Pauli repulsion (r^-12 term)
- `mask.y`: London attraction (r^-6 term without H)
- `mask.z`: Electrostatic (Coulomb Q/r)
- `mask.w`: Hydrogen-bond correction (H-dependent term)

**Output arrays**:
- `Emols[n_samples]`: Total energy (always computed)
- `Emols_masked[n_samples]`: Masked energy (only computed if bMask=true)
- Allows: Etot = Ein + Eout analysis in one pass

**Usage in energy accumulation**:
```opencl
float4 E_components;  // {pauli, london, electro, hbond}
//<<<MODEL_PAIR_DECOMP  // compute components separately

float E_pair = E_components.x + E_components.y + E_components.z + E_components.w;
Ei += E_pair;  // always add total

if(bMask){
    float4 m = mask[mask_idx];
    float E_masked = dot(m, E_components);
    Ei_masked += E_masked;
}
```

**Performance gating**:
- `bMask` boolean flag gates masking logic
- If bMask=false: skip mask read, skip masked accumulation
- Default mask = float4(0,0,0,0) for initialization

**Python-side atom filtering**:
- Homogeneous H masking: set mask for pairs where both atoms have H-component (REQH.w != 0)
- Epair masking: set mask for pairs involving electron pairs (identified by ieps indices)
- General atom filtering: set mask for specific atom indices (e.g., exclude H atoms in formaldehyde H-bond analysis)
- i is always from fragment 1, j is always from fragment 2

**Example use case - Formaldehyde H-bond analysis**:
- System: C2H=O (formaldehyde)
- H-bond involves oxygen and its electron pairs
- H atom should not participate in H-bond but we want to know its role
- Create mask that includes O and epairs, excludes H
- Compare Etot (full) vs Ein (O+epairs only) vs Eout (H contribution)

**Example masks**:
- Baseline (classical): mask={1,1,1,0}
- H-correction only: mask={0,0,0,1}
- Pauli only: mask={1,0,0,0}
- All components: mask={1,1,1,1}
- Atom-specific: mask[i,j] = 0 or 1 based on atom indices

### 2. Kernel Modification

**Target**: `evalSampleEnergy_template` in [FitREQ.cl](cci:7://file:///home/prokop/git/FireCore-fitREQH/cpp/common_resources/cl/FitREQ.cl:0:0-0:0)

**Changes needed**:
1. Add mask buffer to kernel signature:
   ```opencl
   __kernel void evalSampleEnergy_template(
       __global int4*    ranges,
       __global float4*  tREQHs,
       __global int*     atypes,
       __global int2*    ieps,
       __global float4*  atoms,
       __global float4*  mask,      // NEW: [ni*nj] or [nAtomTot*nAtomTot] per-pair mask
       __global float*   Emols,
       int      useTypeQ,
       float4   globParams
   )
   ```

2. Load mask per pair inside loop:
   ```opencl
   int mask_idx = iL * nj + jl;  // or compute from global atom indices
   float4 m = mask[mask_idx];
   ```

3. Modify energy accumulation to use decomposition:
   - Currently: `Ei += Eij;` (from `MODEL_PAIR_ENERGY`)
   - Change to: Use `MODEL_PAIR_DECOMP` to get 4 components
   - Apply mask to each component separately

4. **Key insight**: `MODEL_LJQH2_PAIR_DECOMP` already exists in Forces.cl!
   - It outputs: `pauli`, `london`, `electro`, `hbond`
   - This is exactly what we need for the mask
   - No need to create new decomposition logic

### 3. Memory Layout Considerations

**Mask buffer options**:

**Option 1: Full matrix `[nAtomTot * nAtomTot]`**
- Size: `N^2 * 16` bytes
- Pros: Simple indexing, works for any fragment split
- Cons: Large memory, includes unused pairs

**Option 2: Fragment-pair matrix `[ni * nj]` per sample**
- Size: `ni * nj * 16` bytes per sample
- Pros: Minimal memory, only relevant pairs
- Cons: Need to reshape/reorganize for different fragment splits

**Option 3: Per-sample flattened array**
- Size: `sum(ni*nj) * 16` bytes across all samples
- Pros: Compact
- Cons: Complex indexing with offsets

**Recommendation**: Option 2 (per-sample `[ni, nj]`) - matches current `ranges` structure, minimal memory

### 4. Python Integration

**FittingDriver.py changes**:
1. Add `prepare_mask_buffer(mask_array)` method
2. Add `setup_energy_kernel_with_mask(mask_buff, bMask)` method
3. Modify `evaluate_energies_masked(mask_buff, bMask)` to launch masked kernel
   - Returns tuple: (Emols, Emols_masked) if bMask=true
   - Returns (Emols, None) if bMask=false
4. Add mask generation helpers (create_mask, create_mask_baseline, filter_by_atom_indices, etc.)
5. Add atom filtering helpers (filter_epair_interactions, filter_homogeneous_H)

**Mask generation helper**:
```python
def create_mask(pauli=1.0, london=1.0, electro=1.0, hbond=1.0):
    """Create uniform mask array for all pairs."""
    mask = np.zeros((ni, nj, 4), dtype=np.float32)
    mask[:,:,0] = pauli
    mask[:,:,1] = london
    mask[:,:,2] = electro
    mask[:,:,3] = hbond
    return mask

def create_mask_baseline():
    """Classical non-covalent only (no H-bond)."""
    return create_mask(pauli=1.0, london=1.0, electro=1.0, hbond=0.0)

def create_mask_hbond_only():
    """Only H-bond correction."""
    return create_mask(pauli=0.0, london=0.0, electro=0.0, hbond=1.0)

def filter_by_atom_indices(mask, include_indices_i=None, include_indices_j=None):
    """Filter mask to only include specific atom indices.
    i is always from fragment 1, j is always from fragment 2.
    """
    if include_indices_i is not None:
        for i in range(ni):
            if i not in include_indices_i:
                mask[i,:] = 0.0
    if include_indices_j is not None:
        for j in range(nj):
            if j not in include_indices_j:
                mask[:,j] = 0.0
    return mask

def filter_epair_interactions(mask, ieps):
    """Filter mask to only include pairs involving epairs."""
    # Set mask to zero for pairs where neither i nor j is an epair
    for i in range(ni):
        for j in range(nj):
            is_epair_i = (ieps[i].x >= 0 or ieps[i].y >= 0)
            is_epair_j = (ieps[j].x >= 0 or ieps[j].y >= 0)
            if not (is_epair_i or is_epair_j):
                mask[i,j] = 0.0
    return mask

def filter_homogeneous_H(mask, tREQHs, atypes):
    """Filter mask to only include pairs where both atoms have H-component (REQH.w != 0)."""
    for i in range(ni):
        for j in range(nj):
            ti = atypes[i]
            tj = atypes[j]
            if tREQHs[ti].w == 0.0 or tREQHs[tj].w == 0.0:
                mask[i,j] = 0.0
    return mask
```

### 5. Model Macro Modifications

**Existing `MODEL_LJQH2_PAIR_DECOMP` macro** (Forces.cl:90-113) - ALREADY DOES WHAT WE NEED:
```opencl
{
    // Distance safeguards
    float r_safe  = fmax(r, R2SAFE);
    float inv_r = 1.0f / r_safe;
    
    // Electrostatic
    electro += Q * inv_r * COULOMB_CONST;
    
    // Lennard-Jones 12-6 with H2 scaling
    float u   = R0 * inv_r;
    float u3  = u*u*u;
    float u6  = u3*u3;
    float u12 = u6*u6;
    
    // Decompose LJ energy:
    pauli  += E0 * u12;           // Repulsive r^-12 term
    london += -2.f * E0 * u6;     // Attractive r^-6 term (without H)
    hbond  += E0 * H * u12;       // H-bond correction (H-dependent)
}
```

**Key insight**: This macro already outputs the 4 components we need!
- `pauli`: Repulsive r^-12
- `london`: Attractive r^-6 without H
- `electro`: Coulomb
- `hbond`: H-dependent correction

**Implementation**:
- Replace `MODEL_PAIR_ENERGY` with `MODEL_PAIR_DECOMP` in energy kernel
- Add mask application after decomposition
- No need to modify the macro itself

### 6. Test Design

**Test script**: `tests/tFitREQ/test_energy_mask.py`

**Test cases**:
1. **Identity test (bMask=false)**: No masking → should match unmasked energy, Emols_masked should be None
2. **Identity test (bMask=true, mask=all-ones)**: All masks = 1.0 → Emols should equal Emols_masked
3. **Baseline-only test**: mask={1,1,1,0} → Pauli + London + Electro (no H-bond)
4. **H-bond-only test**: mask={0,0,0,1} → only Hydrogen-bond correction
5. **Pauli-only test**: mask={1,0,0,0} → only repulsion
6. **London-only test**: mask={0,1,0,0} → only dispersion
7. **Electro-only test**: mask={0,0,1,0} → only Coulomb
8. **Sum test**: Emols should equal Emols_masked + (Emols - Emols_masked) for any mask
9. **Epair filtering test**: Filter mask to only include pairs with epairs
10. **Atom filtering test**: Filter mask to specific atom indices (e.g., exclude H in formaldehyde)
11. **Homogeneous H test**: Filter mask to only include pairs with REQH.w != 0 for both atoms

**Validation**:
```python
# Load data
drv = FittingDriver()
drv.load_atom_types(atypes_file)
drv.load_data(xyz_file)

# Test 1: bMask=false (no masking)
E_full, E_masked_none = drv.evaluate_energies_masked(None, bMask=False)
assert E_masked_none is None

# Test 2: Identity (all-ones mask)
mask_all = create_mask(pauli=1, london=1, electro=1, hbond=1)
E_full2, E_masked_all = drv.evaluate_energies_masked(mask_all, bMask=True)
assert np.allclose(E_full2, E_masked_all)
assert np.allclose(E_full, E_full2)

# Test 3: Baseline (no H-bond)
mask_baseline = create_mask(pauli=1, london=1, electro=1, hbond=0)
E_full3, E_baseline = drv.evaluate_energies_masked(mask_baseline, bMask=True)
E_out = E_full3 - E_baseline  # H-bond contribution
assert np.allclose(E_full3, E_baseline + E_out)

# Test 4: H-bond only
mask_hbond = create_mask(pauli=0, london=0, electro=0, hbond=1)
E_full4, E_hbond = drv.evaluate_energies_masked(mask_hbond, bMask=True)
assert np.allclose(E_hbond, E_out)  # Should match H-bond contribution from test 3

# Test 5: Individual components
mask_pauli = create_mask(pauli=1, london=0, electro=0, hbond=0)
mask_london = create_mask(pauli=0, london=1, electro=0, hbond=0)
mask_electro = create_mask(pauli=0, london=0, electro=1, hbond=0)
E_full5, E_pauli = drv.evaluate_energies_masked(mask_pauli, bMask=True)
_, E_london = drv.evaluate_energies_masked(mask_london, bMask=True)
_, E_electro = drv.evaluate_energies_masked(mask_electro, bMask=True)
assert np.allclose(E_full5, E_pauli + E_london + E_electro + E_hbond)

# Test 6: Epair filtering
mask_epairs = create_mask(pauli=1, london=1, electro=1, hbond=1)
mask_epairs = filter_epair_interactions(mask_epairs, drv.host_ieps)
E_full6, E_epairs = drv.evaluate_energies_masked(mask_epairs, bMask=True)
E_non_epair = E_full6 - E_epairs  # Non-epair contribution

# Test 7: Atom filtering (e.g., exclude H in formaldehyde)
# Assuming atom 0 is H, atom 1 is C, atom 2 is O
mask_no_H = create_mask(pauli=1, london=1, electro=1, hbond=1)
mask_no_H = filter_by_atom_indices(mask_no_H, include_indices_i=[1,2], include_indices_j=[1,2])
E_full7, E_no_H = drv.evaluate_energies_masked(mask_no_H, bMask=True)
E_H_only = E_full7 - E_no_H  # H contribution

# Test 8: Homogeneous H filtering
mask_homog_H = create_mask(pauli=1, london=1, electro=1, hbond=1)
mask_homog_H = filter_homogeneous_H(mask_homog_H, drv.tREQHs_base, drv.host_atypes)
E_full8, E_homog_H = drv.evaluate_energies_masked(mask_homog_H, bMask=True)
```

**Visualization**:
- Plot energy maps for: Etot, Ein, Eout (all from one pass)
- Compare baseline vs full energy
- Show contribution percentages per sample
- Example: Formaldehyde H-bond analysis plot showing O+epairs vs H contribution

### 7. Performance Considerations

**Overhead analysis**:
- Extra global memory read: 16 bytes per pair (mask) - ONLY if bMask=true
- Extra computation: 3 multiplies + 2 adds per pair - ONLY if bMask=true
- Extra output buffer: Emols_masked array - ONLY if bMask=true
- Estimated cost when bMask=true: ~5-10% slowdown
- Estimated cost when bMask=false: ~0% overhead (gated execution)

**Optimizations**:
- bMask gating: skip mask read and masked accumulation when bMask=false
- Cache mask in local memory if reused across samples
- Use `__constant` memory if mask is same for all samples
- Default mask = float4(0,0,0,0) for initialization

### 8. Implementation Order

1. **Add missing DECOMP macro** (Forces.cl):
   - Add `MODEL_LJr8QH2_PAIR_DECOMP` macro (based on LJr8QH2_PAIR)
   - Decompose into pauli, london, electro, hbond components

2. **Modify energy kernel** (FitREQ.cl):
   - Add mask buffer and bMask flag to `evalSampleEnergy_template`
   - Add Emols_masked output buffer
   - Replace `MODEL_PAIR_ENERGY` with `MODEL_PAIR_DECOMP`
   - Initialize Ei and Ei_masked accumulators
   - Always compute total energy (Ei += pauli + london + electro + hbond)
   - Conditionally compute masked energy (if bMask: Ei_masked += dot(m, components))
   - Write both outputs (Emols always, Emols_masked if bMask)
   - Use LATOMS[iL].x for total, LATOMS[iL].y for masked

3. **Python integration** (FittingDriver.py):
   - Add `prepare_mask_buffer(mask_array)` method
   - Add `setup_energy_kernel_with_mask(mask_buff, bMask)` method
   - Modify `evaluate_energies_masked(mask_buff, bMask)` to return (Emols, Emols_masked)
   - Add mask generation helpers (create_mask, create_mask_baseline, etc.)
   - Add atom filtering helpers (filter_by_atom_indices, filter_epair_interactions, filter_homogeneous_H)

4. **Test script** (tests/tFitREQ/test_energy_mask.py):
   - Implement bMask=false test (no masking)
   - Implement identity test (all-ones mask)
   - Implement component separation tests
   - Implement epair filtering test
   - Implement atom filtering test (formaldehyde example)
   - Implement homogeneous H filtering test
   - Add visualization (Etot, Ein, Eout plots)

5. **Documentation** (FitREQ_ocl_mask.md):
   - Document mask format and usage
   - Add example commands
   - Document bMask gating behavior

### 9. Open Questions

1. **Mask scope**: Should mask be per-sample or global? Per-sample allows different fragment splits per configuration.

2. **Backward compatibility**: Should `evalSampleEnergy_template` always require mask, or have optional mask parameter with default all-ones?

3. **Macro availability**: Does `MODEL_PAIR_DECOMP` exist for all models?
   - `MODEL_MorseQ_PAIR_DECOMP`: YES (Forces.cl:67-88)
   - `MODEL_LJQH2_PAIR_DECOMP`: YES (Forces.cl:90-113)
   - `MODEL_LJr8QH2_PAIR_DECOMP`: YES (Forces.cl:115-138)
   - Action: All required macros now exist

---

## Final Implementation (May 2026)

### How It Works

**plot_masked_energy.py** - Energy decomposition visualization script:
1. Loads reference energy grid from XYZ file using `_build_frame_grid()` (same method as `plot_ref.py`)
2. Initializes `FittingDriver`, loads atom types and data
3. Extracts `MODEL_MorseQ_PAIR_DECOMP` macro from Forces.cl using **brace matching** (not regex - regex was over-capturing)
4. Compiles with the macro via `compile_with_model(macros={'MODEL_PAIR_DECOMP': morse_decomp})`
5. Computes model energies for **ALL samples in parallel on GPU** (key fix - not per-frame loop)
6. Uses `evaluate_energies_masked(mask)` with different masks:
   - `mask_all`: pauli=1, london=1, electro=1, hbond=1 (total)
   - `mask_hbond`: pauli=0, london=0, electro=0, hbond=1 (H-bond only)
   - `mask_baseline`: pauli=1, london=1, electro=1, hbond=0 (baseline)
7. Returns `(Emols, Emols_masked)` - we use the **masked** result for decomposition
8. Reshapes model energies to match reference grid shape using `frame_idx` from `_build_frame_grid`
9. Plots with **symmetric color limits**: vmin = 1.2 * min(E_ref), vmax = -vmin (applied to all panels)
10. Uses `plot_system_panel()` for consistent visualization with `plot_ref.py`

**FittingDriver.py** - OpenCL driver interface:
- `evaluate_energies_masked(mask)`: Returns `(Emols, Emols_masked)` - total and masked energies
- `create_mask(ni, nj, pauli, london, electro, hbond)`: Creates uniform mask array
- `prepare_mask_buffer(mask)`: Converts mask to GPU buffer
- `setup_energy_kernel(mask_buff, bMask)`: Sets up kernel with masking enabled
- `filter_by_atom_indices()`: Filter mask to specific atom indices
- `filter_epair_interactions()`: Filter to pairs involving epairs
- `filter_homogeneous_H()`: Filter to pairs where both atoms have REQH.w != 0

**FitREQ.cl** - Energy kernel:
- Uses `MODEL_PAIR_DECOMP` macro to compute energy components
- Accumulates `E_components = {pauli, london, electro, hbond}`
- Always computes total: `Ei += pauli + london + electro + hbond`
- Conditionally computes masked: if `bMask`: `Ei_masked += dot(mask, E_components)`
- Writes both outputs: `Emols[sample]` (total) and `Emols_masked[sample]` (masked)

**Forces.cl** - Model macros:
- `MODEL_MorseQ_PAIR_DECOMP`: Decomposes Morse potential into 4 components
  - Pauli: `E0 * e^2` (repulsive e^2 term without H)
  - London: `-2 * E0 * e` (attractive -2*e term)
  - Electro: `Q/r * COULOMB_CONST` (Coulomb term)
  - H-bond: `E0 * H * e^2` (H-dependent correction)

### Problems Faced and Solutions

**Problem 1: Zero energies in initial implementation**
- Cause: Incorrect macro extraction using regex - over-captured content beyond the macro
- Solution: Use **brace matching** to find matching `{` and `}` for precise macro extraction

**Problem 2: Serial per-frame loop defeating GPU parallelism**
- Cause: Initial implementation looped over frames, updated atom positions, ran kernel per frame
- Solution: Compute **all samples in parallel** on GPU - one kernel launch for all 703 samples

**Problem 3: Reference data faked/scrambled in plots**
- Cause: Using wrong data for plots, not properly parsing XYZ with `_build_frame_grid`
- Solution: Load reference energies using exact `_build_frame_grid()` method from `plot_ref.py`, map samples to grid using `frame_idx`

**Problem 4: Etot and Ein looked identical**
- Cause: Using wrong return value from `evaluate_energies_masked()` - took total instead of masked
- Solution: Use the **second return value** `Emols_masked` for decomposition plots

**Problem 5: Color scaling wrong, vmin/vmax totally off**
- Cause: Forcing global min/max across all panels, not using symmetric zero-centered scale
- Solution: Use symmetric scale: vmin = 1.2 * min(E_ref), vmax = -vmin, apply to ALL panels

**Problem 6: Spurious horizontal colorbar overlapping panels**
- Cause: Adding horizontal colorbar while individual panels already had vertical colorbars
- Solution: Remove horizontal colorbar, keep individual vertical colorbars from `plot_system_panel`

**Problem 7: Double conversion (Hartree vs eV)**
- Cause: Manually converting data with Hartree->kcal (627.509) when energies are actually in eV, and `plot_imshow` also converts with eV->kcal (23.060548)
- Solution: Remove manual conversion, let `plot_imshow` handle eV->kcal conversion (23.060548), only apply conversion to color limits to match

### Current Status

**Working**: Energy decomposition pipeline is now functional:
- Reference energies loaded correctly from XYZ files
- Model energies computed in parallel on GPU with decomposition
- Proper grid mapping using `frame_idx`
- Symmetric color scaling with vmin/vmax
- Four panels: Reference, Model Total, Model Ein (H-bond), Model Eout (baseline)

**Not working**: Epairs and hydrogen corrections give unreasonable results
- For H2O-A1_HF-D1-y: Ein = 0 (H-bond component zero)
- For H2O-A1_H2O-D1-y: Ein has non-zero values but correlation is low
- For HCN-A1_HF-D1-y: Ein = 0 again
- **Root cause**: `HBOND_GATE` in FitREQ.cl gates H to zero when `REQi.w * REQj.w >= 0`
  - `#define APPLY_H_GATE(Hval) ((Hval) < 0.f ? (Hval) : 0.f)`
  - If both atoms have same-sign `w`, product is positive → H gated to 0
  - This is likely a parameter/sign convention issue, not a code bug

### How Epairs Work

**Epairs are normal atom-like particles** with the same interaction model as real atoms:
- **Morse potential**: Pauli (repulsion) + London (attraction) components
- **Electrostatics**: Coulomb interaction Q/r
- **H-bond correction**: H-dependent term (REQH.w * REQj.w contribution)

The decomposition treats epairs identically to real atoms - they participate in all 4 components. The difference is in:
- **Directionality**: Epairs have positions/orientations that encode angular dependence
- **Identification**: Filtered by `ieps` indices (epairs have `ieps[i].x >= 0 or ieps[i].y >= 0`)
- **Parameter values**: Epairs have their own REQH parameters (R0, E0, Q, H)

**Epairs vs homogeneous H correction**:
- **Homogeneous H**: Applied to real atoms with REQH.w != 0, modifies Pauli/London via H factor
- **Epairs**: Separate particles with explicit positions, interact via full potential including H term

Both are meant to encode directionality, but epairs do it through explicit particle positions/orientations, while homogeneous H does it through modifying the atom-atom interaction strength.

### Usage Example

```bash
# Run energy decomposition visualization
cd tests/tFitREQ
PYTHONPATH=/home/prokop/git/FireCore-fitREQH:$PYTHONPATH python plot_masked_energy.py --xyz add_epairs_full/H2O-A1_HF-D1-y.xyz --output H2O_HF_energy_decomp.png

# With different systems
python plot_masked_energy.py --xyz add_epairs_full/HCN-A1_HF-D1-y.xyz --output HCN_HF_energy_decomp.png
```

The script generates a 4-panel figure:
1. **Reference**: QM energies from XYZ file
2. **Model Etot**: Total model energy (all components)
3. **Model Ein**: H-bond component only
4. **Model Eout**: Baseline (total - H-bond)

---

## Epairs Model Implementation (May 2026)

### Overview

Implemented a new radial basis function epair model (`iEpairs=4`) with robust parameter handling and plotting support. The new model uses polynomial radial dependence instead of exponential decay, providing more flexibility for fitting H-bond interactions.

### New Epair Model: MODEL_Ep4_PAIR_DECOMP

**Radial basis function formula**:
```
f(r) = 1 - (r/R0)^2
g(r) = f(r)^4
E = f(r) * Yep + g(r) * Xep
```

**Parameters** (stored in REQH columns for epair types):
- `R0ep` (H column, REQH.w): Radial cutoff distance
- `Xep` (R column, REQH.x): Coefficient for g(r) (longer-range term)
- `Yep` (E column, REQH.y): Coefficient for f(r) (shorter-range term)
- `Qep` (Q column, REQH.z): Charge for epair charge subtraction

**Implementation in Forces.cl**:
```opencl
//>>>macro MODEL_Ep4_PAIR_DECOMP
{
    float R0_ = fmax(R0ep, R2SAFE);
    if(r > R0_){ /* Skip if beyond cutoff */ } else {
        float u   = r / R0_;
        float f   = 1.0f - u * u;
        if(f < 0.0f) f = 0.0f;  // Truncate negative values
        f = f * f;
        float g   = f * f;  // f(r)^4
        E_components.w += sEpairs * (f * Yep + g * Xep);
    }
}
```

**Key differences from existing models**:
- `iEpairs=1`: Exponential decay `Hep * exp(-r/R0ep)`
- `iEpairs=2`: Gaussian decay `Hep * exp(-(r/R0ep)^2)`
- `iEpairs=3`: Sigmoid-like `Hep * 2/(exp(r/R0ep) + exp(-r/R0ep))`
- `iEpairs=4`: Polynomial radial basis (NEW) with separate Xep, Yep coefficients

### Parameter Handling Improvements

**Problem**: Previous implementation applied `sqrt(E)` to all atom types, which was incorrect for epair types where `E` is a direct coefficient.

**Solution**: Conditional sqrt application based on atom type:
```python
if isinstance(type_name, str) and type_name.startswith('E'):
    self.tREQHs_base[i, 1] = params['E']  # Store directly for epairs
else:
    self.tREQHs_base[i, 1] = np.sqrt(params['E'])  # sqrt for normal types
```

**Applied consistently in**:
- `prepare_host_data()`: Initial parameter loading
- `prepare_host_data_energy_only()`: Energy-only path
- `apply_type_overrides()`: JSON parameter overrides
- `dump_used_type_params()`: Parameter printing

**Runtime warnings** (for epair types only):
- Warn if `Q==0`: Epair charge subtraction disabled → baseline Coulomb may be too attractive
- Warn if `R==0 and E==0`: Epair coefficients zero → Ein will be ~0

### Kernel Compilation Fix

**Problem**: `MODEL_Ep4_PAIR_DECOMP` macro was not being extracted from Forces.cl and passed to the compiler, causing undefined macro errors when `iEpairs=4`.

**Solution**: Added `MODEL_Ep4_PAIR_DECOMP` to macro extraction in `compile_with_model()`:
```python
if 'MODEL_Ep4_PAIR_DECOMP' not in macros and 'MODEL_Ep4_PAIR_DECOMP' in mlib:
    macros['MODEL_Ep4_PAIR_DECOMP'] = mlib['MODEL_Ep4_PAIR_DECOMP']
```

### Epair Detection Fix

**Problem**: `host_isEpair` array was only initialized in `load_data()`, but `init_and_upload_energy_only()` was called directly in `plot_masked_energy.py`, so epair atoms weren't being flagged.

**Solution**: Added epair detection initialization in `prepare_host_data_energy_only()`:
```python
if not hasattr(self, 'host_isEpair'):
    self.host_isEpair = np.zeros((len(self.host_atypes),), dtype=np.int32)
    if len(self.atom_type_names) > 0:
        ep_type_ids = [i for i, n in enumerate(self.atom_type_names) if (isinstance(n, str) and n.startswith('E'))]
        if ep_type_ids:
            m = np.isin(self.host_atypes, np.array(ep_type_ids, dtype=np.int32))
            self.host_isEpair[m] = 1
```

### Kernel Macro Fix

**Problem**: Initial implementation used `return` statement in `MODEL_Ep4_PAIR_DECOMP` macro, which exits the entire kernel (not just the macro), skipping all subsequent pair calculations.

**Solution**: Changed to `if-else` block to avoid early return:
```opencl
if(r > R0_){ /* Skip if beyond cutoff */ } else {
    // Compute epair energy
}
```

### Plotting Enhancements

**Added support for multiple plotting modes in `plot_masked_energy.py`**:
- `--polar 0`: Non-polar mode (default)
- `--polar 1`: Polar mode with angle-dependent visualization
- Automatic plane detection: xz, zy, or xy based on geometry
- Consistent axis handling to avoid duplicate axes and IndexError

**Refactored plotting code**:
- Created reusable function `plot_profile_row()` for second row of profile plots
- Fixed axis handling in both polar and non-polar modes
- Consistent axes array usage across plot calls

### Parameter Override System

**Added JSON-based parameter override** in `epair_defaults.json`:
```json
{
  "iEpairs": 4,
  "sEpairs": 1.0,
  "O_3": { "H": 0.0 },
  "E_O3": { "R": -0.2, "E": -0.2, "Q": 0.0, "H": 2.0 }
}
```

**Features**:
- Global epair settings: `iEpairs` (model type), `sEpairs` (global scaling)
- Per-atom type overrides: R, E, Q, H parameters
- Applied via `apply_type_overrides()` in FittingDriver
- Supports both normal atom types and epair types (names starting with 'E')

### Epair Type Distinguishing

**Previously**: Epair types used constant parameters across all systems.

**Now**: Epair types have system-specific parameters:
- `E_O3`: Oxygen lone pairs in water
- `E_N1`, `E_N2`, `E_N3`: Nitrogen lone pairs in different bonding environments
- `E_O1`, `E_O2`, `E_OR`: Oxygen lone pairs in different bonding environments
- `E_F`, `E_HO`, `E_HN`, `E_HC1`, `E_HF`, `E_Ha`: Various other epair types

**Benefits**:
- More accurate modeling of directionality
- Different epair types can have different radial cutoffs and coefficients
- Flexible fitting to QM data

### Baseline Energy Consistency

**Problem**: Toggling `iEpairs` between 3 and 4 caused baseline energy (Eout) to change unexpectedly.

**Root cause**: Inconsistent sqrt application - `iEpairs==4` stored raw E for all types, breaking baseline.

**Solution**: Removed `iEpairs` dependency from sqrt logic:
- Epair types (name starts with 'E'): Always store E directly
- Normal types: Always store sqrt(E)
- Baseline now independent of `iEpairs` setting

### Verification

**Test results with `iEpairs=4` and `E_O3: {R: -0.2, E: -0.2, Q: 0.0, H: 2.0}`**:
```
Energy Statistics:
  Reference: min=-0.196003, max=0.489446, mean=-0.068214
  Model Total: min=-0.170498, max=1.164019, mean=0.014476
  Model H-bond: min=-0.387343, max=0.000000, mean=-0.077862
  Model Eout: min=-0.105115, max=1.367815, mean=0.092339
Reference-Model correlation: 0.927035
```

**Key observations**:
- Ein (H-bond) now has non-zero values (min=-0.387343, mean=-0.077862)
- Baseline (Eout) remains consistent when toggling `iEpairs`
- Total energy correlation improved to 0.93

### Usage Example

```bash
# Run with new epair model (iEpairs=4)
cd tests/tFitREQ
python plot_masked_energy.py \
  --xyz add_epairs_full/H2O-A1_H2O-D1-y.xyz \
  --params epair_defaults.json \
  --polar 0

# Test baseline consistency by toggling iEpairs
# Edit epair_defaults.json: change "iEpairs": 4 to 3
# Run again - Eout should be identical
```

### Files Modified

1. **Forces.cl**: Added `MODEL_Ep4_PAIR_DECOMP` macro with radial basis function
2. **FitREQ.cl**: Added `iEpairs==4` case in epair pair decomposition, passes Xep and Yep parameters
3. **FittingDriver.py**:
   - Conditional sqrt application (epairs vs normal types)
   - Added `MODEL_Ep4_PAIR_DECOMP` to macro extraction
   - Added epair detection in `prepare_host_data_energy_only()`
   - Added runtime warnings for epair parameter issues
4. **plot_masked_energy.py**: Enhanced plotting with polar/non-polar modes, plane detection
5. **epair_defaults.json**: Parameter override file with epair type-specific settings

### Summary

Successfully implemented a new polynomial radial basis epair model with:
- Flexible radial dependence (f(r) = 1 - (r/R0)^2, g(r) = f(r)^4)
- Separate Xep and Yep coefficients for different radial ranges
- Robust parameter handling (epairs store E directly, normal types store sqrt(E))
- Baseline energy independent of `iEpairs` setting
- JSON-based parameter override system
- System-specific epair type parameters
- Enhanced plotting with multiple modes and plane detection
- Runtime warnings for parameter misconfiguration