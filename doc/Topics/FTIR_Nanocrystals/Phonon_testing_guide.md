# Phonon Frequency Testing Guide

**Purpose**: This document provides a reproducible workflow for testing phonon frequencies in FireCore MMFF forcefield, enabling comparison with other computational methods (DFTB, LAMMPS, DFT/pySCF-B3LYP).

**Date**: June 2026
**Status**: Diamond phonons verified; diatomic bond stiffness tested; bond order assignment bug documented

---

## Quick Start

### Diamond Phonon Bands (PBC Mode)

```bash
cd /home/prokop/git/FireCore/tests/tMMFF
bash run.sh test_diamond_phonon_bands.py --pbc --unit THz --super-n 3
```

**Output**:
- Data: `diamond_phonon_bands_THz_PBC.npz` (cached force constants and frequencies)
- Plot: `diamond_phonon_bands_THz_PBC.png` (band structure)
- Console: Γ-point frequencies, negative frequency count

**Expected Results**:
- Γ-point optical modes: **41.25 THz** (triple-degenerate)
- No negative frequencies (< -0.01 THz)
- Acoustic modes go to zero at Γ

### Plot from Cached Data (No Re-computation)

```bash
cd /home/prokop/git/FireCore/tests/tMMFF
python3 plot_phonon_bands.py diamond_phonon_bands_THz_PBC.npz --unit THz
```

**Note**: Use standalone `plot_phonon_bands.py` instead of `--plot-only` flag to avoid ASan/matplotlib crashes.

---

## Test Scripts

### 1. Diamond Phonon Bands

**Script**: `tests/tMMFF/test_diamond_phonon_bands.py`

**Purpose**: Compute phonon band structure for diamond primitive cell using supercell finite-difference force constants.

**Key Features**:
- Supports PBC mode (`--pbc`) and cluster mode (default)
- Phonopy-compatible workflow: displace central cell atoms, read forces on all supercell atoms
- Bloch sum to k-dependent dynamical matrix
- Optional ASR (acoustic sum rule) symmetrization (`--asr`)

**Usage**:
```bash
# PBC mode (recommended for bulk crystals)
bash run.sh test_diamond_phonon_bands.py --pbc --unit THz --super-n 3

# Cluster mode (may have surface artifacts)
bash run.sh test_diamond_phonon_bands.py --unit THz --super-n 3

# With phonopy ASR symmetrization
bash run.sh test_diamond_phonon_bands.py --asr --unit THz --super-n 3
```

**Arguments**:
- `--unit`: Frequency unit (`THz` or `cm-1`, default: `cm-1`)
- `--super-n`: Supercell size (odd integer, default: 3)
- `--pbc`: Enable periodic boundary conditions
- `--asr`: Apply phonopy acoustic sum rule symmetrization
- `--dx`: Finite-difference displacement in Bohr (default: 1e-4)

**Output Files**:
- `diamond_phonon_bands_{unit}_{mode}.npz`: Cached data (k-points, frequencies, Phi blocks)
- `diamond_phonon_bands_{unit}_{mode}.png`: Band structure plot

**Data Structure (.npz)**:
- `kdist`: Cumulative k-path distance
- `freqs`: Frequencies at each k-point (shape: nkpoints × 3×natoms)
- `kpts`: k-point vectors in reciprocal space
- `k_frac`: Fractional k-point coordinates
- `prim_lvec`: Primitive cell lattice vectors
- `recip`: Reciprocal lattice vectors
- `Phi_blocks`: Real-space force constant blocks
- `super_n`: Supercell size
- `pbc`: PBC flag
- `asr`: ASR flag

### 2. Diatomic Hessian Test

**Script**: `tests/tMMFF/test_diatomic_hessian.py`

**Purpose**: Test Hessian eigenvalues for diatomic molecules to isolate bond stretch contributions.

**Key Features**:
- Uses MDL .mol format with explicit bond definitions
- Supports H₂, C₂, and C₂_sp3 test cases
- Reports raw Hessian eigenvalues, mass-weighted eigenvalues, and frequencies
- Compares maximum eigenvalue to expected bond stiffness

**Usage**:
```bash
# H2 molecule
bash run.sh test_diatomic_hessian.py --case H2 --dx 1e-4

# C2 with automatic atom type detection
bash run.sh test_diatomic_hessian.py --case C2 --dx 1e-4

# C2 with forced C_3 atom type (tests C-C single bond parameters)
bash run.sh test_diatomic_hessian.py --case C2_sp3 --dx 1e-4
```

**Arguments**:
- `--case`: Test case (`H2`, `C2`, `C2_sp3`)
- `--l0`: Bond length (Å)
- `--uff`: Use UFF forcefield instead of MMFF
- `--dx`: Finite-difference displacement (Å, default: 1e-4)

**Output**: Console printout of eigenvalues and frequencies

**Known Issue**: H₂ test fails due to MMFF sp³ forcefield requiring back-neighbors (topology limitation).

### 3. Diatomic Energy Scan

**Script**: `tests/tMMFF/test_diatomic_scan.py`

**Purpose**: Scan energy along bond direction to directly measure bond stiffness via d²E/dx².

**Key Features**:
- Displaces one atom along bond axis
- Measures energy at each displacement
- Fits quadratic curve around minimum to extract stiffness
- Compares to bond table parameters

**Usage**:
```bash
# C2 with C_3 atom types
bash run.sh test_diatomic_scan.py --case C2_sp3 --l0 1.538 --range 0.2 --nstep 41

# C2 with C_R atom types
bash run.sh test_diatomic_scan.py --case C2_R --l0 1.538 --range 0.2 --nstep 41
```

**Arguments**:
- `--case`: Test case (`H2`, `C2`, `C2_sp3`, `C2_R`)
- `--l0`: Bond length (Å)
- `--range`: Scan range (Å, default: 0.3)
- `--nstep`: Number of scan points (default: 61)

**Output**: Console printout of energy vs bond length and fitted stiffness

### 4. Ethane Gamma Hessian

**Script**: `tests/tMMFF/test_ethane_gamma.py`

**Purpose**: Test Hessian eigenvalues for ethane molecule at Γ-point.

**Key Features**:
- Compares MMFF and UFF forcefields
- `--bondsOnly` flag to disable angles and pi terms for bond-only stiffness test
- Reports stiffness eigenvalues and mass-weighted frequencies

**Usage**:
```bash
# Full forcefield
bash run.sh test_ethane_gamma.py

# Bond-only (angles and pi terms disabled)
bash run.sh test_ethane_gamma.py --bondsOnly
```

**Arguments**:
- `--bondsOnly`: Disable angles, PiSigma, and PiPiI terms

---

## Data Storage

### Directory Structure
```
/home/prokop/git/FireCore/tests/tMMFF/
├── test_diamond_phonon_bands.py
├── test_diatomic_hessian.py
├── test_diatomic_scan.py
├── test_ethane_gamma.py
├── plot_phonon_bands.py
├── run.sh
├── diamond_phonon_bands_THz_PBC.npz     # Cached diamond data
├── diamond_phonon_bands_THz_PBC.png     # Diamond band plot
└── data_UFF/                            # Forcefield parameters
    ├── BondTypes.dat
    ├── AngleTypes.dat
    └── AtomTypes.dat
```

### Cached Data Files

**Diamond PBC (current reference)**:
- File: `diamond_phonon_bands_THz_PBC.npz`
- Supercell: 3×3×3
- Mode: PBC enabled
- Γ-point optical modes: 41.25 THz
- Negative frequencies: 0

**Loading in Python**:
```python
import numpy as np
data = np.load('diamond_phonon_bands_THz_PBC.npz')
freqs = data['freqs']  # Shape: (280, 6)
kdist = data['kdist']
k_frac = data['k_frac']
```

---

## Current Status and Caveats

### Verified Results

1. **Diamond Phonons (PBC Mode)**
   - ✅ Γ-point optical modes: 41.25 THz (physically reasonable)
   - ✅ No spurious negative frequencies
   - ✅ PBC correctly propagates via `neighCell` arrays
   - ✅ Force constants extracted correctly from supercell FD

2. **Diatomic Bond Stiffness**
   - ✅ Energy scan confirms d²E/dx² matches assigned bond k
   - ✅ Hessian λ_max = 2×k (collective stretch mode of both atoms)
   - ✅ Physical interpretation: 2× factor is correct for diatomic

### Known Issues

1. **Bond Order Assignment Bug**
   - **Issue**: MMFF forcefield ignores explicit bond order from MDL .mol files
   - **Evidence**: C_3-C_3 diatomic assigned triple bond (k=36.309) instead of single bond (k=10.082)
   - **Impact**: Cannot validate single-bond stiffness using diatomic C_3 or C_R types
   - **Workaround**: Use molecules where forcefield correctly infers single bonds from coordination (e.g., ethane)

2. **H₂ Diatomic Test**
   - **Issue**: MMFF sp³ forcefield requires back-neighbors; H₂ has none
   - **Error**: `ERROR makeBackNeighs() capping atom[0] has not back-neighbor => exit()`
   - **Status**: Topology limitation, not a bug in Hessian code

3. **Matplotlib ASan Crash**
   - **Issue**: `--plot` flag crashes under AddressSanitizer due to matplotlib/ft2font incompatibility
   - **Workaround**: Use standalone `plot_phonon_bands.py` script instead of `--plot-only`

### Comparison with Other Methods

When comparing with DFTB, LAMMPS, or DFT/pySCF-B3LYP:

1. **Frequency Units**
   - FireCore reports frequencies in THz by default
   - Conversion: 1 THz = 33.356 cm⁻¹
   - Use `--unit cm-1` flag for direct comparison with DFT results

2. **Forcefield Parameters**
   - MMFF bond parameters: `cpp/common_resources/BondTypes.dat`
   - MMFF angle parameters: `cpp/common_resources/AngleTypes.dat`
   - MMFF atom types: `cpp/common_resources/AtomTypes.dat`

3. **Expected Ranges**
   - Diamond optical modes: ~40 THz (MMFF) vs ~40 THz (DFT)
   - C-C single bond stretch: ~30-40 THz (depends on molecule)
   - Note: MMFF is a classical forcefield; quantum methods may show different anharmonicity

4. **Finite-Difference Parameters**
   - Default displacement: 1e-4 Bohr (≈0.005 Å)
   - Can adjust with `--dx` flag for convergence testing
   - Smaller dx → more accurate but more numerical noise

---

## Reproduction Workflow for External Agents

### Step 1: Clone FireCore
```bash
git clone <FireCore_repo_url>
cd FireCore
```

### Step 2: Build C++ Library
```bash
cd cpp/Build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Step 3: Run Diamond Phonon Test
```bash
cd ../../tests/tMMFF
bash run.sh test_diamond_phonon_bands.py --pbc --unit THz --super-n 3
```

### Step 4: Load and Compare Results
```python
import numpy as np

# Load FireCore results
data = np.load('diamond_phonon_bands_THz_PBC.npz')
fc_freqs = data['freqs']  # THz
fc_kdist = data['kdist']

# Load your method's results (e.g., DFTB, LAMMPS, DFT)
# Assuming same k-path structure
other_freqs = np.load('other_method.npz')['freqs']

# Compare at Γ-point
print(f"FireCore Γ-point: {fc_freqs[0]} THz")
print(f"Other method Γ-point: {other_freqs[0]} THz")

# Plot comparison
import matplotlib.pyplot as plt
plt.plot(fc_kdist, fc_freqs, label='FireCore MMFF')
plt.plot(fc_kdist, other_freqs, label='Other method')
plt.legend()
plt.ylabel('Frequency (THz)')
plt.xlabel('k-path')
plt.savefig('comparison.png')
```

### Step 5: Diatomic Bond Stiffness Comparison
```bash
# Test C-C bond (note: will use triple bond due to forcefield bug)
bash run.sh test_diatomic_scan.py --case C2_sp3 --l0 1.538 --range 0.2 --nstep 41
```

Extract d²E/dx² from output and compare with:
- Bond table value: 36.309 eV/Å² (triple bond)
- Your method's bond stiffness (if available)

---

## Contact and Support

For questions about:
- Forcefield parameters: See `cpp/common_resources/*.dat` files
- Code implementation: See `cpp/common/molecular/MolWorld_sp3.h`
- Python API: See `pyBall/MMFF.py`
- Bug reports: Check `doc/Topics/FTIR_Nanocrystals/Debug_negative_phonon_freqs.md`

---

## References

1. **Forcefield Parameters**: MMFF parameter files in `cpp/common_resources/`
2. **Debug Log**: `doc/Topics/FTIR_Nanocrystals/Debug_negative_phonon_freqs.md`
3. **Phonopy**: Togo et al., "First-principles phonon calculations in materials science", Scr. Mater. 108 (2015)
