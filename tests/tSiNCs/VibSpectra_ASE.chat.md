# USER

I would like to compute hessian and vibration spectra of Adamantane and sillicon analog of adamantane. Freely available codes I can use for it easily. 
DFTB+
pySCF
psi4
cp2k
anything else esaily accessible by ASE (Atomistic simulation environment)

---

# Gemini-3.5-flash

To compute the Hessian and vibrational spectra (Normal Modes and IR intensities) of **adamantane** ($C_{10}H_{16}$) and its all-silicon analog, **sila-adamantane** ($Si_{10}H_{16}$), using the Atomic Simulation Environment (ASE) is a highly efficient choice.

Since you are dealing with finite, rigid molecular cages with high symmetry ($T_d$), you can utilize ASE's native `ase.vibrations` or `ase.vibrations.Infrared` modules. These packages handle the finite-difference displacements automatically, meaning any calculator that can provide **forces** (for frequencies) or **dipole moments** (for IR intensities) will work out of the box.

## 1. Quick Structural Setup

First, define your molecules. You can easily construct them via their SMILES strings or basic coordinates.

```python
from ase.data.pubchem import pubchem_atoms
from ase.build import molecule

# Adamantane via PubChem or SMILES
adamantane = pubchem_atoms('adamantane') 

# For Sila-adamantane: simplest way is mutating the carbon backbone to silicon
sila_adamantane = adamantane.copy()
for atom in sila_adamantane:
    if atom.symbol == 'C':
        atom.symbol = 'Si'
# Note: Si-Si bonds (~2.35 Å) are longer than C-C (~1.54 Å). 
# You MUST scale or relax this structure deeply before running Hessian.
sila_adamantane.scale_positions(1.5) 

```

## 2. Choosing the Best Calculator

Depending on your precision needs and computational budget, here is how the codes you mentioned fit this specific problem:

### Method A: DFTB+ (Fastest for initial screening)

Excellent choice for quick screening. The `Skf` parameter sets like `pbc-0-3` or `mio-1-1` contain well-tested $C-H$ and $Si-H$ interactions.

* **Pros:** Extremely fast; can compute IR intensities because it yields charges/dipoles.
* **Cons:** Approximated; frequencies will have a systematic shift (scaling factor needed).

### Method B: Psi4 or PySCF (Best for Isolated Quantum Chemistry)

Because these are isolated, highly symmetric molecules, standard Gaussian-basis molecular codes are ideal.

* **Pros:** Highly accurate molecular point-group symmetry handling; exact analytic electronic gradients.
* **Cons:** PySCF's ASE calculator interface can occasionally be minimalist regarding native dipole properties, whereas **Psi4** integrates seamlessly with ASE for both energy/forces and dipole moments.

### Method C: CP2K (Best if transitioning to bulk/surfaces)

* **Pros:** Mixed Gaussian/Plane-Wave (GPW) scheme makes it incredibly fast for medium-to-large systems.
* **Cons:** It treats systems using periodic boundary conditions (PBC). You must place the adamantane in a sufficiently large vacuum box (e.g., $L = 20.0$ Å) to avoid artificial periodic interactions, and ensure wavefunctions are decoupled (Poisson solver settings). *Note: The native ASE-CP2K interface sometimes struggles to pass dipoles back to ASE's `Infrared` module; if you need IR intensities with CP2K, it is often easier to use CP2K's internal `&VIBRATIONAL_ANALYSIS` block directly.*

## 3. Implementation Blueprint (ASE + Calculator)

Here is a clean template script using **Psi4** (as an example of an ab-initio method) to fully optimize the structure, calculate the Hessian, and output the vibrational spectrum.

```python
import numpy as np
from ase.calculators.psi4 import Psi4
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.vibrations import Infrared

# 1. Setup Calculator (e.g., B3LYP or PBE0 with a split-valence basis set)
calc = Psi4(
    atoms=adamantane,
    method='b3lyp',
    basis='6-311g(d,p)',
    charge=0,
    multiplicity=1
)
adamantane.calc = calc

# 2. CRITICAL: Tight Geometry Optimization
# Frequencies are valid only at a true stationary point (minimum).
opt = BFGS(adamantane, trajectory='opt.traj')
opt.run(fmax=0.001) # Tight threshold for clean Hessian

# 3. Compute Hessian and Normal Modes via Finite Difference
# This will displace each atom by delta (default 0.01 Å) in ±x, ±y, ±z
vib = Vibrations(adamantane, name='vib_adamantane')
vib.run()

# Summary of frequencies (cm^-1)
vib.summary()

# 4. (Optional) Get Infrared Intensities
# Requires the calculator to return 'dipole' properties
ir = Infrared(adamantane, name='ir_adamantane')
ir.run()
ir.summary()

# Write a trajectory file to visualize a specific mode (e.g., mode idx 12)
vib.write_mode(12)

```

## 4. Key Physics Gotchas for $C_{10}H_{16}$ vs $Si_{10}H_{16}$

When you analyze your results, look out for these distinct physical shifts:

| Property | Adamantane ($C_{10}H_{16}$) | Sila-adamantane ($Si_{10}H_{16}$) |
| --- | --- | --- |
| **Symmetry** | $T_d$ (high degeneracy) | $T_d$ (high degeneracy) |
| **Cage Breathing Mode** | Distinctly high energy ($\sim 750\text{ cm}^{-1}$) due to rigid, strong $C-C$ bonds. | Drastically red-shifted ($\sim 300\text{--}400\text{--}500\text{ cm}^{-1}$) due to heavier $Si$ mass and weaker $Si-Si$ backbones. |
| **Exocyclic Stretch** | $C-H$ stretching ($\sim 2900\text{ cm}^{-1}$). | $Si-H$ stretching ($\sim 2100\text{ cm}^{-1}$). |
| **Imaginary Frequencies** | Should be exactly 0 (or $< 10\text{ cm}^{-1}$ numerical noise for the first 6 translational/rotational modes). | If you get large imaginary frequencies, your initial scale factor for the $Si$ cage size was insufficient—re-optimize with a looser optimizer before tightening it down. |

If you plan to scale this setup to larger diamondoids or wish to include explicit solvent/surface environments later, we can adjust the CP2K setup parameters to handle it efficiently. Which calculator framework do you have readily compiled in your current workspace?

# Implementation Review and Practical Notes

## Summary of Tested Methods

### DFTB+ (Recommended for Production)
**Status:** ✅ Working, fastest option

**Installation:**
- Requires DFTB+ binary installed (tested with `/home/prokop/git_SW/dftbplus/Build/app/dftb+/dftb+`)
- Requires Slater-Koster (SK) parameter files
- SK parameter locations used: `/home/prokop/SIMULATIONS/dftbplus/slakos/`
  - `mio-1-1/` for C, H systems
  - `matsci-0-3/` for Si, H systems

**ASE Interface:**
```python
from ase.calculators.dftb import Dftb
calc = Dftb(
    Hamiltonian_='DFTB',
    Hamiltonian_SlaterKosterFiles_Type='Directory',  # Critical: use 'Directory' not 'Set2SKFiles'
    Hamiltonian_SlaterKosterFiles_Prefix='/path/to/sk/files/',  # Must have trailing slash
    Hamiltonian_MaxAngularMomentum_C='p',
    Hamiltonian_MaxAngularMomentum_H='s',
    Hamiltonian_MaxAngularMomentum_Si='p',
)
```

**Key Issues Resolved:**
1. **SK parameter path format:** Must use `Type='Directory'` with trailing slash in path. Using `Type='Set2SKFiles'` caused file not found errors.
2. **MaxAngularMomentum:** Must only specify angular momentum for atoms actually present in the system. Specifying all possible atoms caused DFTB+ to fail.
3. **Structure source:** Used `ase.data.pubchem.pubchem_atoms_search('adamantane')` instead of hardcoded coordinates.

**Performance:**
- Adamantane (26 atoms): ~30 seconds for optimization + ~2 minutes for vibrations
- Sila-adamantane (26 atoms): ~30 seconds for optimization + ~2 minutes for vibrations
- Total: ~5 minutes for both molecules

**Accuracy:**
- Systematic frequency shift expected (DFTB+ is approximate)
- Qualitatively correct: C-H stretches at ~2900 cm⁻¹, Si-H stretches at ~2150 cm⁻¹
- IR intensities computed successfully

### PySCF (Alternative for Accuracy)
**Status:** ⚠️ Works but slow, requires additional dependencies

**Installation:**
```bash
pip install pyscf geometric
```

**Issues Encountered:**
1. **No ASE calculator interface:** PySCF doesn't have a native ASE calculator, must use PySCF's own optimization and Hessian routines
2. **Geometry optimization dependency:** Requires `geometric` package for optimization (not included by default)
3. **Slow:** HF/STO-3G test took several minutes; DFT with larger basis would be much slower
4. **No IR intensities:** PySCF's harmonic analysis doesn't easily provide IR intensities through ASE

**When to Use:**
- When you need high accuracy and can afford longer computation time
- For benchmarking against DFTB+ results
- When you need specific quantum chemistry methods not available in DFTB+

### Psi4 (Not Tested - Installation Failed)
**Status:** ❌ Not available in current environment

**Installation Issue:**
```bash
pip install psi4  # Failed: no matching distribution found
```

**Notes:**
- Would likely work well if installed (good ASE integration per documentation)
- May require conda installation instead of pip
- Would be slower than DFTB+ but more accurate

## Structure Setup

### PubChem Integration (Recommended)
```python
from ase.data.pubchem import pubchem_atoms_search
atoms = pubchem_atoms_search('adamantane')
```

**Advantages:**
- No hardcoded coordinates
- Automatic hydrogen placement
- Correct stoichiometry (C₁₀H₁₆ for adamantane)

**For Sila-adamantane:**
```python
sila_adamantane = adamantane.copy()
for atom in sila_adamantane:
    if atom.symbol == 'C':
        atom.symbol = 'Si'
sila_adamantane.positions *= 1.5  # Scale for longer Si-Si bonds
```

## Workflow Summary

### Recommended Production Workflow (DFTB+)

1. **Activate environment:**
   ```bash
   source ~/venvs/ML/bin/activate
   ```

2. **Run vibration calculation:**
   ```bash
   cd tests/tSiNCs
   python vib_spectra_dftb_working.py
   ```

3. **Plot results:**
   ```bash
   python plot_vib_spectra.py
   ```

**Output files:**
- `*_frequencies.npy` - Frequency data (numpy arrays)
- `*_optimized.xyz` - Optimized geometries
- `*_opt.traj` - Optimization trajectories
- `vib_spectra_stick.png` - Stick spectrum plot
- `vib_spectra_broadened.png` - Broadened spectrum plot

## Key Configuration Details

### DFTB+ SK Parameter Sets
- **mio-1-1:** Standard organic parameters (C, H, N, O, P, S)
- **matsci-0-3:** Materials science parameters (includes Si, Na, etc.)
- **3ob-3-1:** Third-order parameters (more accurate but slower)
- **pbc-0-3:** Periodic boundary conditions (not needed for molecules)

### Geometry Optimization
- **Critical:** Must converge to fmax < 0.001 eV/Å for valid Hessian
- **Method:** BFGS through ASE
- **Typical steps:** 50-100 steps for adamantane-sized molecules

### Hessian Calculation
- **Method:** Finite difference (ASE Vibrations module)
- **Displacement:** 0.01 Å (default)
- **Cost:** 6N force evaluations (N = number of atoms)
- **For 26 atoms:** 156 force calculations

## Troubleshooting

### DFTB+ Fails with "SK file does not exist"
**Solution:** Use `Hamiltonian_SlaterKosterFiles_Type='Directory'` with trailing slash in path.

### DFTB+ Fails with "MaxAngularMomentum error"
**Solution:** Only specify angular momentum for atoms present in your system. Build dynamically:
```python
atom_types = set(atoms.get_chemical_symbols())
max_angular = {'Hamiltonian_MaxAngularMomentum_': ''}
for atom in atom_types:
    if atom == 'H':
        max_angular['Hamiltonian_MaxAngularMomentum_H'] = 's'
    elif atom in ['C', 'Si', 'N', 'O', 'P', 'S']:
        max_angular[f'Hamiltonian_MaxAngularMomentum_{atom}'] = 'p'
```

### IR Calculation Fails
**Solution:** This is often a caching issue when reusing vibration data. Clear old cache files or skip IR if frequencies are sufficient.

### Large Imaginary Frequencies
**Solution:** Structure not fully optimized. Re-run optimization with tighter convergence or check initial geometry quality.

## File Locations

### Scripts
- `vib_spectra_dftb_working.py` - Main DFTB+ calculation script
- `plot_vib_spectra.py` - Matplotlib plotting script
- `vib_spectra_pyscf.py` - PySCF alternative (slow)
- `vib_spectra_pyscf_fast.py` - PySCF fast test (no optimization)

### SK Parameters
- `/home/prokop/SIMULATIONS/dftbplus/slakos/mio-1-1/` - C, H parameters
- `/home/prokop/SIMULATIONS/dftbplus/slakos/matsci-0-3/` - Si, H parameters

### CP2K Basis Sets
- `/home/prokop/SIMULATIONS/cp2k/basis/` - CP2K basis set files (not tested)

## Future Improvements

1. **Test CP2K:** Would be useful for comparison, especially for larger systems
2. **Install Psi4:** Would provide high-accuracy benchmark
3. **Frequency scaling:** DFTB+ frequencies may need empirical scaling factors
4. **IR intensities:** Currently computed but could be validated against experiment
5. **Normal mode visualization:** Add mode animation capabilities

---

# Unified Pipeline Implementation (2026-05)

## Consolidation to Universal Scripts

The workflow has been consolidated into a universal pipeline that handles multiple computational methods through a single interface, eliminating the need for method-specific scripts.

### Core Modules

**`vib_utils.py`** - Shared utilities for:
- Structure fetching from PubChem with local caching
- Calculator creation for DFTB+ and PySCF
- Geometry optimization with caching (relaxed geometries saved as `<mol>_<method>_relaxed.xyz`)
- Vibration calculation with caching (frequencies saved as `<mol>_<method>_freq.npy`)
- Frequency filtering (removes imaginary modes below threshold)
- Vibration mode export to XYZ with displacement vectors

**`run_vib_spectra.py`** - Universal CLI runner:
```bash
# Compute spectra for a molecule with specific methods
python run_vib_spectra.py adamantane --methods dftb_mio dftb_3ob

# Export vibration modes to XYZ
python run_vib_spectra.py adamantane --methods dftb_mio --export-modes

# Cache Hessians (PySCF only - ASE Vibrations doesn't expose Hessian)
python run_vib_spectra.py adamantane --methods pyscf_hf --cache-hessians
```

**`plot_vib_spectra.py`** - Overlay plotting:
- Loads all cached frequency files for a molecule
- Overlays spectra from different methods (not different molecules)
- Generates stick and broadened spectra in two-panel figure
- Saves as `<mol>_vib_overlay.png`

### Removed Scripts
The following method-specific scripts are no longer needed:
- `vib_spectra_dftb.py` - replaced by universal pipeline
- `vib_spectra_dftb_working.py` - replaced by universal pipeline
- `vib_spectra_pyscf.py` - replaced by universal pipeline
- `vib_spectra_pyscf_fast.py` - replaced by universal pipeline

## Vibration Mode Export to XYZ

### Multi-Frame XYZ Format

All vibration modes for a method are now exported to a single multi-frame XYZ file with displacement vectors:

**Format:**
```
26
Freq: 313.89 cm-1
C     -0.546753    -0.608974    -1.288074     0.000473     0.028413    -0.013654
C      0.588195     1.376346    -0.297679     0.050192    -0.016995     0.020453
...
26
Freq: 313.91 cm-1
C     -0.546753    -0.608974    -1.288074     0.050324     0.006136    -0.024235
...
```

**Columns:** atom_symbol x y z dx dy dz (displacement vector)

**Output file:** `<mol>_<method>_all_modes.xyz`

**Usage:**
```bash
python run_vib_spectra.py adamantane --methods dftb_mio --export-modes
```

This format allows visualization of eigenmodes in molecular viewers that support extended XYZ with vector data (e.g., VMD, ChemCraft, or custom scripts).

## DFTB+ Basis Set Comparison

### Available SK Parameter Sets

The following Slater-Koster parameter sets are available in `/home/prokop/SIMULATIONS/dftbplus/slakos/`:

| Set | Elements | Use Case |
|-----|----------|----------|
| **mio-1-1** | C, H, N, O, P, S | General organic molecules |
| **3ob-3-1** | C, H, N, O, halogens, metals | 3rd order, broader coverage |
| **matsci-0-3** | Si, C, H, N, O, F, Na, P | Materials science (includes Si) |
| **pbc-0-3** | Si, H, O | Periodic boundary systems |

### Basis Set Selection by Molecule

**Adamantane (C₁₀H₁₆):**
- Recommended: `mio-1-1` (standard organic) vs `3ob-3-1` (3rd order)
- Comparison plot shows frequency shifts between parameterizations

**Sila-adamantane (Si₁₀H₁₆):**
- Recommended: `matsci-0-3` (Si-optimized) vs `pbc-0-3` (periodic)
- Note: `3ob-3-1` lacks Si parameters, cannot be used

### Updated Max Angular Momentum

For better accuracy with Si-containing systems, Si angular momentum increased to 'd':
```python
MAX_ANG = {'H': 's', 'C': 'p', 'N': 'p', 'O': 'p', 'S': 'p', 'P': 'p', 'Si': 'd', 'Ge': 'd'}
```

### Comparison Workflow

```bash
# Compare basis sets for adamantane
python run_vib_spectra.py adamantane --methods dftb_mio dftb_3ob
python plot_vib_spectra.py adamantane

# Compare basis sets for sila-adamantane
python run_vib_spectra.py sila_adamantane --methods dftb_matsci dftb_pbc
python plot_vib_spectra.py sila_adamantane
```

## Hessian Caching - Known Limitations

### ASE Vibrations Hessian Access

**Problem:** ASE's `Vibrations` class does not expose a `get_hessian()` method in version 3.24.0, making it impossible to extract the Hessian matrix directly from DFTB+ calculations.

**Attempted Solutions:**
1. **`ase.hessian.Hessian` module:** Does not exist in ASE 3.24.0
2. **Manual reconstruction from force constants:** Not straightforward with ASE's caching structure

**Current Status:**
- DFTB+ Hessian caching: **Not implemented** (ASE limitation)
- PySCF Hessian caching: **Implemented** (PySCF has native Hessian support)

### PySCF Hessian Caching

For PySCF calculations, Hessian matrices can be cached:
```python
def compute_or_load_hessian_pyscf(mol, mf, mol_name, method_tag, workdir='.'):
    hess_obj = mf.Hessian()
    hess = hess_obj.kernel()  # Returns (n_atoms x 3 x n_atoms x 3) in Hartree/Bohr^2
    np.save(str(hess_path), hess)
```

**Usage:**
```bash
python run_vib_spectra.py adamantane --methods pyscf_hf --cache-hessians
```

### Hessian Comparison Plotting

A `plot_hessians()` function is available to compare cached Hessians side-by-side using `imshow`:
```python
plot_hessians(mol_name, method_tags, workdir='.')
```

**Output:** `<mol>_hessian_comparison.png` with RdBu_r colormap showing Hessian matrix structure.

**Note:** Only works for PySCF methods currently due to ASE limitation for DFTB+.

## Caching Strategy

All intermediate results are cached to avoid redundant computations:

1. **Structure caching:** `<mol>_pubchem.xyz` - fetched from PubChem once
2. **Relaxed geometry:** `<mol>_<method>_relaxed.xyz` - optimization result
3. **Frequencies:** `<mol>_<method>_freq.npy` - vibration frequencies
4. **Hessians:** `<mol>_<method>_hessian.npy` - Hessian matrix (PySCF only)
5. **Vibration cache:** `<mol>_<method>_vib/` - ASE Vibrations force calculations

This caching allows:
- Re-running calculations with different plotting options
- Comparing methods without recomputing
- Debugging intermediate steps
- Fast iteration on analysis