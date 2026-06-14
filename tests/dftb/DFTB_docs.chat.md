# USER

I need to put comprehesive DFT documentation into markdown

I'm especailly interersted in how
1) compute phonon spectra
2) compute reaction paths (nudged elastic band NEB)
3) how to read internal varibles like hamiltonian (eventually charges or density matrix)
4) how to integrat it with python on deeper level (as power usere, hacker, developer)

there are some sources
https://dftbplus-recipes.readthedocs.io/en/stable/ 

https://github.com/by-student-2017/DFTBplus-v.23.1-examples 

https://dftbplus-recipes.readthedocs.io/en/stable/interfaces/ase/neb/index.html

---

# Kimi 2.6

 Here is a comprehensive markdown document covering the four advanced topics you requested. This synthesizes information from the official DFTB+ recipes, recent API developments, and the ASI interface.

# Comprehensive DFTB+ Documentation for Power Users

## Table of Contents
0. [Setup and Configuration](#0-setup-and-configuration)
1. [Phonon Spectra with Phonopy](#1-phonon-spectra-with-phonopy)
2. [Reaction Paths (Nudged Elastic Band - NEB)](#2-reaction-paths-nudged-elastic-band---neb)
3. [Reading Internal Variables (Hamiltonian, Charges, Density Matrix)](#3-reading-internal-variables-hamiltonian-charges-density-matrix)
4. [Deep Python Integration (Developer/Hacker Level)](#4-deep-python-integration-developerhacker-level)

## 0. Setup and Configuration

### Environment Variables

DFTB+ integration in FireCore uses two environment variables:

- **`DFTB_EXE`**: Path to DFTB+ executable
- **`DFTB_SK_PATH`**: Path to library directory containing Slater-Koster parametrization subdirectories

Set these in `~/.bashrc`:

```bash
export DFTB_EXE=/path/to/dftb+
export DFTB_SK_PATH=/path/to/slakos/library/
```

Example:
```bash
export DFTB_EXE=/home/prokop/git_SW/dftbplus/Build/app/dftb+/dftb+
export DFTB_SK_PATH=/home/prokop/SIMULATIONS/dftbplus/slakos/library/
```

### Slater-Koster Parameter Sets

DFTB+ requires Slater-Koster (`.skf`) files for element interactions. The library-based design allows multiple parametrizations in one directory.

**Download parameter sets:**

```bash
mkdir -p ~/SIMULATIONS/dftbplus/slakos/library
cd ~/SIMULATIONS/dftbplus/slakos/library

# 3ob-3-1 (bio/organic, Br-C-Ca-Cl-F-H-I-K-Mg-N-Na-O-P-S-Zn)
wget https://github.com/dftbparams/3ob/releases/download/v3.1.0/3ob-3-1.tar.xz
tar xf 3ob-3-1.tar.xz

# mio-1-1 (simpler, H-C-N-O-S-P)
wget https://github.com/dftbparams/mio/releases/download/v1.1.0/mio-1-1.tar.xz
tar xf mio-1-1.tar.xz
```

Directory structure:
```
~/SIMULATIONS/dftbplus/slakos/library/
├── 3ob-3-1/
│   ├── C-C.skf
│   ├── C-H.skf
│   ├── H-H.skf
│   └── ...
└── mio-1-1/
    ├── C-C.skf
    ├── C-H.skf
    ├── H-H.skf
    └── ...
```

See https://dftb.org/parameters/download.html for all available sets.

### Using Parametrizations in Code

The `pyBall/dftb_utils.py` module provides a library-based interface:

```python
from pyBall import dftb_utils

# Use default parametrization (first available)
dftb_utils.run_pbc(apos, enames, lvs)

# Specify which parametrization to use
dftb_utils.run_pbc(apos, enames, lvs, sk_set='mio-1-1')
dftb_utils.run_pbc(apos, enames, lvs, sk_set='3ob-3-1')

# Check available sets
print(dftb_utils.AVAILABLE_SK_SETS)  # ['3ob-3-1', 'mio-1-1']
print(dftb_utils.DEFAULT_SK_SET)      # '3ob-3-1'
```

### Validation

The module validates at import time:
- Checks `DFTB_EXE` exists and is executable
- Checks `DFTB_SK_PATH` directory exists
- Checks for at least one parametrization subdirectory with `.skf` files

Error messages provide setup instructions and download URLs.

## 1. Phonon Spectra with Phonopy

DFTB+ interfaces with **phonopy** to calculate phonon band structures and related vibrational properties. The workflow requires DFTB+ to compute forces on displaced supercell geometries.

### Prerequisites
- `phonopy` installed (e.g., `mamba install -c conda-forge phonopy`) 
- DFTB+ compiled with standard settings

### Input Preparation

Your `dftb_in.hsd` must enable force calculation and tagged output:

```hsd
Analysis = {
  # REQUIRED: phonopy needs forces
  CalculateForces = Yes
}
Options = {
  # REQUIRED: phonopy reads results.tag
  WriteResultsTag = Yes
}
```

The starting geometry **must** be named `geo.gen` (phonopy convention for DFTB+ interface) .

### Workflow

#### Step 1: Generate Displacements
```bash
phonopy -d --dim="4 4 4" --dftb+
```
- Creates `geo.genS` (undistorted supercell) and `geo.genS-001`, `geo.genS-002`, etc. (displaced structures)
- Generates `phonopy_disp.yaml` (displacement metadata)
- **Convergence check**: Test supercell size and k-point density 

#### Step 2: Run DFTB+ on Each Displacement
For each `geo.genS-*` file:
1. Copy it to working directory as `geo.genS-001` (or appropriate name)
2. Run `dftb+` to get `results.tag`

**Batch processing example** (for multiple displacements):
```bash
for struct in geo.genS-*; do
  dir=$(echo $struct | sed 's/geo.genS-//g')
  mkdir -p $dir
  cp $struct $dir/geo.genS-001
  cp dftb_in.hsd $dir/
  cd $dir && dftb+ && cd ..
done
```

#### Step 3: Collect Forces
```bash
# Single displacement
phonopy -f results.tag --dftb+

# Multiple displacements in subdirectories
phonopy -f [0-9][0-9][0-9]/results.tag --dftb+
```
This creates `FORCE_SETS` from the `results.tag` files and `phonopy_disp.yaml` .

#### Step 4: Calculate Band Structure
Create `band.conf` specifying your Brillouin zone path (see phonopy docs), then:
```bash
phonopy -p band.conf --dim="4 4 4" --dftb+
phonopy-bandplot --gnuplot band.yaml > band.dat
```

### Advanced: Electron-Phonon Coupling
For electron-phonon coupling calculations, use **dftBephy** which builds on phonopy results. Requires `WriteHS = Yes` in DFTB+ input to export Hamiltonian/overlap matrices .

---

## 2. Reaction Paths (Nudged Elastic Band - NEB)

DFTB+ supports NEB through two main pathways: **ASE integration** (recommended for flexibility) and **socket communication** (for performance).

### Method A: ASE NEB Interface

The ASE `NEB` class handles the path optimization while DFTB+ provides energies and forces via the standard calculator interface.

#### Basic Setup
```python
from ase.io import read
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb

# Read endpoints (must have same atom ordering!)
initial = read('initial.xyz')
final = read('final.xyz')

# Create images
n_images = 7
images = [initial]
images += [initial.copy() for _ in range(n_images)]
images += [final]

# Initialize NEB
neb = NEB(images)
neb.interpolate()  # Linear interpolation between endpoints

# Set up DFTB+ calculator for each image
for image in images[1:-1]:  # Endpoints are fixed by default
    calc = Dftb(
        label='neb',
        Hamiltonian_SCC='Yes',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        ParserOptions_ParserVersion=14
    )
    image.set_calculator(calc)

# Optimize
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)  # eV/Å
```

**Critical requirements** :
- Atom ordering must be identical in initial and final states
- Endpoints should be pre-relaxed to local minima
- Start with few images, increase if convergence is poor

#### Climbing Image NEB (CI-NEB)
To find the transition state with higher precision:
```python
from ase.neb import NEB

# climb=True enables CI-NEB for the highest energy image
neb = NEB(images, climb=True)
```

### Method B: Socket Communication (High Performance)

For large systems, use Unix domain sockets to avoid file I/O overhead. This runs DFTB+ as a server and ASE sends coordinates/receives forces via memory .

```python
import os
from subprocess import Popen
from ase.io import read, write
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator

NIMAGES = 13
DFTBP_PATH = 'dftb+'

def main():
    initial = read('initial.traj')
    final = read('final.traj')
    
    images = [initial]
    images += [initial.copy() for _ in range(NIMAGES)]
    images += [final]
    
    neb = NEB(images)
    neb.interpolate()
    
    # Create socket calculators
    socketids = range(1, NIMAGES + 1)
    wdirs = [f'_calc/image_{s}' for s in socketids]
    unixsockets = [f'dftbplus_{s}' for s in socketids]
    
    # Write modified HSD with Driver = Socket for each image
    write_modhsd(socketids)
    
    calcs = [SocketIOCalculator(log='socket.log', unixsocket=us) 
             for us in unixsockets]
    
    for i, calc in enumerate(calcs):
        images[i + 1].set_calculator(calc)
    
    # Launch DFTB+ instances
    for cwd in wdirs:
        Popen(DFTBP_PATH, cwd=cwd)
    
    opt = BFGS(neb, trajectory='path.traj')
    opt.run(fmax=1.0E-02)
    
    for calc in calcs:
        calc.close()

def write_modhsd(socketids):
    for socket in socketids:
        path = f'_calc/image_{socket}'
        os.makedirs(path, exist_ok=True)
        with open(f'{path}/dftb_in.hsd', 'w') as f:
            f.write(f'''  <<+ ../../dftb_in.hsd
+Driver = +Socket {{
    !File = "dftbplus_{socket}"
}}''')
```

**HSD template** (`dftb_in.hsd`) must include:
```hsd
Driver = {}
```
which gets modified at runtime to `Driver = Socket { ... }` .

## 3. Reading Internal Variables (Hamiltonian, Charges, Density Matrix)

DFTB+ provides three levels of access to internal electronic structure data, ranging from file I/O to zero-copy memory access.

### Level 1: File-Based Export (Simplest)

Set in `dftb_in.hsd`:
```hsd
Options {
  WriteHS = Yes           # Square H and S matrices → hamsqr1.dat, oversqr.dat
  WriteRealHS = Yes       # Real-space sparse H and S → hamreal.dat, overreal.dat
}
```

**Warning**: When either option is `Yes`, DFTB+ **only** builds and writes matrices, then exits immediately. No diagonalization or SCC cycles occur .

To build H with converged charges:
```hsd
Options {
  ReadInitialCharges = Yes
  WriteHS = Yes
}
```

### Level 2: Native Python API (Intermediate)

The `dftbplus` Python package (in `tools/pythonapi/`) exposes energies, forces, and Mulliken charges .

```python
import dftbplus
import numpy as np

LIB_PATH = '/path/to/libdftbplus'

# Initialize with existing HSD input
calc = dftbplus.DftbPlus(libpath=LIB_PATH)
calc.register_input_file('dftb_in.hsd')

# Set geometry (in Bohr!)
coords = np.array([...])  # shape: (n_atoms, 3)
latvecs = np.array([...]) # shape: (3, 3) or None
calc.set_coords(coords, latvecs)

# Run calculation
calc.run()

# Extract results (atomic units!)
energy = calc.get_energy()
forces = calc.get_gradients()      # shape: (n_atoms, 3)
charges = calc.get_gross_charges() # Mulliken charges

calc.close()
```

**Key limitation**: The native Python API does **not** directly expose Hamiltonian or density matrices—only scalar/vector quantities like energy, forces, and charges .

### Level 3: ASI API — Zero-Copy Matrix Access (Advanced)

The **Atomic Simulation Interface (ASI)** provides callback-based access to Hamiltonian (H), overlap (S), and density matrices (DM) without memory copies. This is the "hacker" level for machine learning workflows and multiscale simulations .

#### Architecture
- DFTB+ compiled as shared library with ASI support
- `asi4py` Python wrapper installed (`pip install asi4py`)
- Callbacks registered before calculation; invoked by DFTB+ when matrices are ready

#### Python Example with asi4py
```python
import os
import numpy as np
from ctypes import CDLL, RTLD_GLOBAL
from ase.build import molecule
from asi4py.asecalc import ASI_ASE_calculator

# Load ASI library (must be compiled with ASI support)
ASI_LIB_PATH = os.environ['ASI_LIB_PATH']
asilib = CDLL(ASI_LIB_PATH, mode=RTLD_GLOBAL)

# ASE calculator wrapper
def init_dftb(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(
        label='cluster',
        Hamiltonian_SCC='Yes',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        Hamiltonian_ASI_='',
        Hamiltonian_ASI_AsiModifiesModel='No',
        ParserOptions_ParserVersion=14
    )
    calc.write_input(asi.atoms)
    return calc

atoms = molecule('H2O')
atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_dftb, None, atoms)

# Request matrix storage
atoms.calc.asi.keep_density_matrix = True
atoms.calc.asi.keep_hamiltonian = True
atoms.calc.asi.keep_overlap = True

# Run calculation
energy = atoms.get_potential_energy()

# Retrieve matrices (index: k-point, spin channel)
S = atoms.calc.asi.overlap_storage[(1, 1)]  # Overlap matrix
H = atoms.calc.asi.hamiltonian_storage[(1, 1)]  # Hamiltonian
DM = atoms.calc.asi.dm_storage.get((1, 1), None)  # Density matrix

# Verify electron count
n_electrons = np.sum(S * DM.T)
band_energy = np.sum(H * DM.T)
print(f'Electrons: {n_electrons:.6f}')
print(f'Band energy: {band_energy:.6f}')
```

**Matrix formats** :
- Dense matrices with BLACS descriptors (for parallel runs)
- LAPACK triangular storage
- Sparse real-space formats (boundary-condition agnostic)

#### Callback Mechanism (C-level)
For deep integration, register callbacks that DFTB+ invokes:
```c
// C pseudo-code showing callback registration
ASI_register_hamiltonian_callback(my_callback, user_data);
// DFTB+ calls my_callback(hamiltonian_ptr, desc, user_data) when H is ready
```

This enables workflows like:
- **QM/ML**: Export H/S to train neural network Hamiltonians
- **QM/MM**: Import external potentials via `ASI_register_external_potential`
- **Multiscale**: Couple DFTB+ with continuum solvers 

---

## 4. Deep Python Integration (Developer/Hacker Level)

### Compilation Requirements

For full API access, compile DFTB+ with:
```bash
cmake -DWITH_PYTHON=1 \
      -DWITH_API=1 \
      -DENABLE_DYNAMIC_LOADING=1 \
      -DBUILD_SHARED_LIBS=1 \
      -DWITH_ASI=1 \
      ..
```

Produces `libdftbplus.so` and ASI libraries .

### Installation Example (FireCore Setup)

**Library path**: `/home/prokophapala/git_SW/dftbplus/_build/src/dftbp/libdftbplus.so`

**Add to ~/.bashrc**:
```bash
export LD_LIBRARY_PATH=/home/prokophapala/git_SW/dftbplus/_build/src/dftbp:/home/prokophapala/git_SW/dftbplus/_build/external/mudpack:/home/prokophapala/git_SW/dftbplus/_build/external/libnegf/origin/src:$LD_LIBRARY_PATH
source ~/.bashrc
```

**Python API**: Install from `tools/pythonapi`:
```bash
cd /home/prokophapala/git_SW/dftbplus/tools/pythonapi
python -m pip install -e .
```

**Test import**:
```python
from dftbplus import DftbPlus
calc = DftbPlus(libpath='/home/prokophapala/git_SW/dftbplus/_build/src/dftbp/libdftbplus.so')
```

**Available methods** (Level 2 API):
- `get_energy()` - Total energy
- `get_gradients()` - Forces
- `get_gross_charges()` - Mulliken charges
- `get_cm5_charges()` - CM5 charges
- `set_geometry(coords, latvecs)` - Set atomic coordinates (coords in Bohr)
- `set_external_potential()` - Set external potential
- `register_ext_pot_generator()` - Register external potential callback
- `get_nr_atoms()` - Number of atoms

**Note**: Coordinates must be in Bohr when using `set_geometry()`.

### API Layers

| Layer | Access | Use Case |
|-------|--------|----------|
| **Native Python API** | Energy, forces, charges, external potentials | Standard ASE workflows |
| **ASI + asi4py** | H, S, DM matrices, ESP callbacks | ML training, multiscale |
| **C/Fortran API** | Full control, custom solvers | Method development |

### Advanced: Population-Dependent External Potentials

Register Python callbacks for self-consistent embedding :

```python
class PotentialCalculator:
    def __init__(self, qmcoords, mmcoords, mmcharges):
        self.qmcoords = qmcoords
        self.mmcoords = mmcoords
        self.mmcharges = mmcharges
    
    def calc_extpot(self, dqatom):
        """Return external potential at QM atoms based on current charges"""
        # dqatom: population difference from neutral atom reference
        extpot = np.zeros(len(self.qmcoords))
        # ... compute electrostatic potential from MM charges ...
        return extpot
    
    def calc_extpotgrad(self, dqatom):
        """Return gradients of external potential"""
        return np.zeros((len(self.qmcoords), 3))

# Register with DFTB+
potcalc = PotentialCalculator(qm_coords, mm_coords, mm_charges)
calc.register_ext_pot_generator(potcalc, get_extpot, get_extpotgrad)
```

### Workflow Integration Patterns

#### Pattern 1: ASE + SocketIO (Production NEB)
Use for geometry optimizations where DFTB+ runs as a persistent server:
```python
from ase.calculators.socketio import SocketIOCalculator
calc = SocketIOCalculator(unixsocket='dftb_socket')
atoms.set_calculator(calc)
# ASE sends coords, receives forces via Unix socket
```

#### Pattern 2: ASI Driver (Python Controls DFTB+)
Python drives the calculation, DFTB+ acts as a library:
```python
# Initialize DFTB+ via ASI
asi = ASIlib(ASI_LIB_PATH)
asi.init('dftb_in.hsd', 'output.log', mpi_comm)
asi.set_geometry(coords)
asi.run()

# Extract data via callbacks registered pre-run
H = asi.hamiltonian_storage[(1,1)]
```

#### Pattern 3: Custom Eigensolvers
Replace DFTB+'s diagonalization by registering matrix callbacks and solving in Python:
```python
# DFTB+ builds H and S, then calls your callback
def my_solver(h_ptr, s_ptr, n_basis, eigenvalues, eigenvectors):
    H = np.array(h_ptr[:n_basis**2]).reshape(n_basis, n_basis)
    S = np.array(s_ptr[:n_basis**2]).reshape(n_basis, n_basis)
    e, v = scipy.linalg.eigh(H, S)
    eigenvalues[:] = e
    eigenvectors[:] = v.flatten()

# Register before asi.run()
asi.register_eigenvalue_callback(my_solver)
```

## Summary Cheat Sheet

| Task | Method | Key Input/Flag |
|------|--------|----------------|
| **Phonons** | phonopy + DFTB+ | `CalculateForces = Yes`, `WriteResultsTag = Yes` |
| **NEB** | ASE NEB + DFTB+ calc | `NEB(images, climb=True)` |
| **NEB (fast)** | SocketIO + DFTB+ | `Driver = Socket { ... }` |
| **Charges** | Native Python API | `get_gross_charges()` |
| **H/S/DM matrices** | ASI + asi4py | `keep_hamiltonian = True` |
| **External potential** | Native API callback | `register_ext_pot_generator()` |
| **File export H/S** | Direct file I/O | `WriteHS = Yes` |

## References

- [DFTB+ Recipes Documentation](https://dftbplus-recipes.readthedocs.io/en/stable/)
- [DFTB+ GitHub Examples](https://github.com/by-student-2017/DFTBplus-v.23.1-examples)
- [ASI API Documentation](https://pvst.gitlab.io/asi)
- [Phonopy DFTB+ Interface](https://phonopy.github.io/phonopy/dftb%2B.html)
- Recent DFTB+ developments: J. Phys. Chem. A 2025 
- ASI paper: arXiv:2403.15625 

---

# Gemini-3-flash

This documentation provides a comprehensive guide for using DFTB+ for advanced atomistic simulations, focusing on phonon spectra, reaction paths (NEB), extracting internal variables, and deep Python integration.

---

# 1. Phonon Spectra with DFTB+ and Phonopy
DFTB+ does not compute the full phonon band structure natively in one step; it is typically used as a force calculator for **Phonopy**.

### Workflow
1.  **Relax the Structure:** Ensure your geometry is perfectly optimized (Force tolerance $< 10^{-4}$ or $10^{-5}$ au).
2.  **Generate Displacements:** Use Phonopy to create supercells with small atomic displacements.
    ```bash
    phonopy -d --dim=\"2 2 2\" --dftb+
    ```
    This generates `geo.genS` (supercell) and `geo.genS-001`, `geo.genS-002`, etc.
3.  **Run DFTB+ for Forces:** Run a single-point calculation for each displaced structure. Ensure `CalculateForces = Yes` and `WriteResultsTag = Yes` are set in `dftb_in.hsd`.
4.  **Collect Forces:**
    ```bash
    phonopy -f disp-*/results.tag --dftb+
    ```
5.  **Post-process:** Create a `band.conf` and plot:
    ```bash
    phonopy -p band.conf --dim=\"2 2 2\" --dftb+
    ```



---

# 2. Reaction Paths (Nudged Elastic Band - NEB)
For NEB, the Atomic Simulation Environment (ASE) is the preferred interface. It manages the \"images\" of the reaction path and uses DFTB+ to calculate the forces on each image.

### Python Implementation (ASE)
```python
from ase.optimize import BFGS
from ase.mep import NEB
from ase.calculators.dftb import Dftb
from ase.io import read

# Load initial and final states
initial = read('initial.gen')
final = read('final.gen')

# Create 5 intermediate images
images = [initial]
for i in range(5):
    images.append(initial.copy())
images.append(final)

# Set calculators for intermediate images
for image in images[1:-1]:
    image.calc = Dftb(Hamiltonian_SCC='Yes', 
                      Hamiltonian_MaxAngularMomentum_H='s',
                      Hamiltonian_MaxAngularMomentum_O='p')

# Interpolate positions
neb = NEB(images)
neb.interpolate()

# Optimize the path
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)
```



---

# 3. Reading Internal Variables (Hamiltonian & Density Matrix)
To access the \"guts\" of the calculation (Hamiltonian $H$, Overlap $S$, or Density Matrix $P$), you must instruct DFTB+ to write them during the SCC cycle.

### Via HSD Input
In the `Options` or `Analysis` block of `dftb_in.hsd`:
```hsd
Analysis {
  WriteHS = Yes           # Writes hamsqrN.dat and oversqr.dat
  WriteRealHS = Yes       # Writes real-valued matrices
}
```
*Note: Setting `WriteHS = Yes` often causes the code to stop after building the matrices but before the SCC cycle finishes. For the final converged density matrix, use `WriteDetailedXML = Yes`.*

### Files to Look For:
- **`oversqr.dat`**: Overlap matrix.
- **`hamsqr0.dat`**: Zero-charge Hamiltonian.
- **`charges.bin`**: Binary file containing Mulliken charges (can be read by `read_bin` scripts in the DFTB+ source).
- **`detailed.xml`**: Contains the most comprehensive set of converged data.

---

# 4. Deep Python Integration (Developer/Power User)
For hackers and developers, you should move beyond simple \"File-IO\" wrappers and use **Socket Communication** or the **C/Python API**.

### A. The Socket Interface (Performance)
Rather than restarting the DFTB+ binary for every force call (which is slow due to overhead), use the **SocketCalculator** in ASE. This keeps DFTB+ running in the background as a \"server.\"
```python
from ase.calculators.socketio import SocketIOCalculator
# This allows DFTB+ to communicate via a Unix socket or TCP
```

### B. Accessing the Library (`libdftbplus`)
If you are developing a new method (e.g., machine learning potentials), you can link against `libdftbplus.so`. 
- **The API:** DFTB+ provides a C-API (`capi.f90`) which can be wrapped using `ctypes` or `cython` in Python.
- **Internal Storage:** The `detailed.xml` file is the standard way to programmatically exchange data. Use Python's `lxml` or `xml.etree` to parse it:
```python
import xml.etree.ElementTree as ET
tree = ET.parse('detailed.xml')
# Access charges, dipole moments, or energy components
e_total = tree.find(\".//total_energy\").text
```

### C. Advanced: Data Extraction with dftbplus Python Package

The official dftbplus python tools (often found in the tools/ folder of the source code) include dataset utilities to convert binary outputs into NumPy arrays for deep learning or custom analysis.

---

# USER

do you know how to obtain charge density, molecular orbitals, and electrostaci potnetial from DFTB at some real spalce points, or projected on grid? Ideallly using python api (or ASE if necessary)

Also how to obtain hessian (deynamical matrix) that is change of energy due to samll diplacement of pairs of atoms in caresian direction

Dij = ddE/(dxi dj)

---

# Kimi 2.6

Here is a comprehensive markdown document addressing your two questions about obtaining real-space grid data (charge density, MOs, ESP) and the Hessian matrix from DFTB+, with a focus on Python/ASE workflows.

## 1. Real-Space Grid Data: Charge Density, Molecular Orbitals, and Electrostatic Potential

DFTB+ does **not** compute real-space volumetric data directly during the main calculation. Instead, it provides specialized post-processing tools and output formats that you can drive programmatically from Python.

### 1.1 Charge Density and Molecular Orbitals (via Waveplot)

The **`waveplot`** utility is the official tool for generating Gaussian cube files from DFTB+ results. It evaluates the LCAO wavefunctions on a real-space grid using the Slater-type orbital coefficients from your SK parameter set .

#### Step 1: Required DFTB+ Output

First, run DFTB+ with these flags in `dftb_in.hsd`:

```hsd
Options {
  WriteDetailedXml = Yes    # System info for waveplot
}

Analysis {
  WriteEigenvectors = Yes   # Binary eigenvectors needed for MOs
}
```

This produces `detailed.xml` and `eigenvec.bin` .

#### Step 2: Waveplot Input (`waveplot_in.hsd`)

```hsd
Options {
  TotalChargeDensity = Yes        # wp-abs2.cube
  TotalChargeDifference = Yes     # wp-abs2diff.cube
  ChargeDensity = Yes             # Per-orbital density
  RealComponent = Yes             # Real part of wavefunction
  PlottedSpins = 1 -1
  PlottedLevels = 4               # HOMO for H2O
  PlottedRegion = OptimalCuboid {}
  NrOfPoints = 50 50 50
  NrOfCachedGrids = -1
  Verbose = Yes
}

DetailedXml = "detailed.xml"
EigenvecBin = "eigenvec.bin"

Basis {
  Resolution = 0.01
  <<+ "../../slakos/wfc/wfc.mio-1-1.hsd"
}
```

**Critical**: The `wfc.*.hsd` file contains STO coefficients and must match your SK set .

#### Step 3: Python Automation Wrapper

Since waveplot is a standalone executable, wrap it in Python:

```python
import subprocess
import os
from pathlib import Path

def run_waveplot(waveplot_exe='waveplot', workdir='.', sk_wfc_path=None):
    """
    Run waveplot programmatically from Python.
    
    Parameters:
        waveplot_exe: path to waveplot binary
        workdir: working directory containing detailed.xml and eigenvec.bin
        sk_wfc_path: path to wavefunction coefficient file (e.g., wfc.mio-1-1.hsd)
    """
    cwd = Path(workdir).absolute()
    
    # Generate waveplot_in.hsd dynamically
    wfc_include = f'  <<+ "{sk_wfc_path}"' if sk_wfc_path else '  <<+ "wfc.mio-1-1.hsd"'
    
    waveplot_input = f'''
Options {{
  TotalChargeDensity = Yes
  TotalChargeDifference = Yes
  ChargeDensity = Yes
  RealComponent = Yes
  PlottedSpins = 1 -1
  PlottedLevels = 1:-1
  PlottedRegion = OptimalCuboid {{}}
  NrOfPoints = 50 50 50
  NrOfCachedGrids = -1
  Verbose = Yes
}}

DetailedXml = "detailed.xml"
EigenvecBin = "eigenvec.bin"

Basis {{
  Resolution = 0.01
{wfc_include}
}}
'''
    with open(cwd / 'waveplot_in.hsd', 'w') as f:
        f.write(waveplot_input)
    
    # Run waveplot
    result = subprocess.run(
        [waveplot_exe],
        cwd=cwd,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"Waveplot failed: {result.stderr}")
    
    # Collect output cube files
    cubes = list(cwd.glob('wp-*.cube'))
    return cubes

# Example usage
cubes = run_waveplot(
    waveplot_exe='/path/to/waveplot',
    workdir='./calc',
    sk_wfc_path='/path/to/slakos/wfc/wfc.mio-1-1.hsd'
)
```

**Output files** :
- `wp-abs2.cube` — Total charge density
- `wp-abs2diff.cube` — Density difference (molecular - atomic superposition)
- `wp-K-S-L-abs2.cube` — Charge density for orbital L, spin S, k-point K
- `wp-K-S-L-real.cube` — Real part of wavefunction

#### Reading Cube Files in Python

Use `ase.io.cube` or `scipy`:

```python
from ase.io.cube import read_cube_data
import numpy as np

# Returns (data, atoms)
data, atoms = read_cube_data('wp-abs2.cube')
print(f"Grid shape: {data.shape}")  # (50, 50, 50)
print(f"Grid origin (Bohr): {atoms.cell.origin}")

# Access value at specific grid point
value = data[25, 25, 25]
```

### 1.2 Electrostatic Potential (ESP)

DFTB+ can evaluate the electrostatic potential at arbitrary points or on a regular grid via the `Analysis` block .

#### Point-by-Point Evaluation

```hsd
Analysis {
  ElectrostaticPotential {
    OutputFile = "ESP.dat"
    Softening = 1E-6
    Points [Angstrom] = {
      0.0  0.0  0.0
      1.0  0.0  0.0
      0.0  1.0  0.0
    }
  }
}
```

#### Grid Evaluation (Recommended for Visualization)

```hsd
Analysis {
  ElectrostaticPotential {
    OutputFile = "ESP.dat"
    Softening = 1E-6
    Grid [Angstrom] = {
      GridPoints = 50 50 50
      Origin = -5.0 -5.0 -5.0
      Spacing = 0.2 0.2 0.2
      Directions = 1 0 0  0 1 0  0 0 1
    }
  }
}
```

**Output format** (`ESP.dat`) :
```
# T 1 1 1 125000
# 0.000000000000E+00 -0.200000000000E+01 -0.200000000000E+01
# 0.200000000000E+01 0.000000000000E+00 0.000000000000E+00
# 0.000000000000E+00 0.200000000000E+01 0.000000000000E+00
# 0.000000000000E+00 0.000000000000E+00 0.200000000000E+01
# Internal (V) External (V)
0.173386318927E-10 0.314737193575E+00
...
```

#### Python Parser for ESP Grid Data

```python
import numpy as np

def read_esp_grid(filename='ESP.dat'):
    """
    Read DFTB+ electrostatic potential grid output.
    Returns: (grid_data, origin, vectors, has_external_field)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    header = lines[0].split()
    has_field = header[1] == 'T'
    is_grid = [int(x) for x in header[2:5]]
    n_points = int(header[5])
    
    # Grid vectors (next 4 lines after header)
    origin = np.array([float(x) for x in lines[1].split()[1:]])
    vec1 = np.array([float(x) for x in lines[2].split()])
    vec2 = np.array([float(x) for x in lines[3].split()])
    vec3 = np.array([float(x) for x in lines[4].split()])
    
    # Read potential values
    data_lines = lines[6:]  # Skip comment line
    if is_grid == [1, 1, 1]:
        # Gridded data: no position columns
        potentials = []
        for line in data_lines:
            parts = line.split()
            if has_field:
                potentials.append([float(parts[0]), float(parts[1])])
            else:
                potentials.append([float(parts[0])])
        grid_shape = tuple(is_grid)
        potentials = np.array(potentials).reshape((*grid_shape, -1))
    else:
        # Non-gridded: has position columns
        potentials = []
        for line in data_lines:
            parts = line.split()
            pos = [float(parts[0]), float(parts[1]), float(parts[2])]
            if has_field:
                potentials.append(pos + [float(parts[3]), float(parts[4])])
            else:
                potentials.append(pos + [float(parts[3])])
        potentials = np.array(potentials)
    
    return potentials, origin, (vec1, vec2, vec3), has_field

# Usage
esp_data, origin, vectors, has_ext = read_esp_grid('ESP.dat')
internal_potential = esp_data[:,:,:,0]  # Shape: (nx, ny, nz)
```

### 1.3 Transport/Poisson Solver: Real-Space Density and Potential

For transport calculations, DFTB+ uses a real-space Poisson solver that outputs `potential.dat` and `charge_density.dat` on a grid. These can be converted to cube format using the `makecube` tool :

```hsd
Electrostatics = Poisson {
  PoissonBox [Angstrom] = 40.0 30.0 30.0
  MinimalGrid [Angstrom] = 0.5 0.5 0.5
  SavePotential = Yes
}
```

Then convert:
```bash
makecube potential.dat -r refpot.dat
```

## 2. Hessian Matrix (Dynamical Matrix)

The Hessian in DFTB+ is defined as the matrix of second derivatives of the energy with respect to atomic Cartesian coordinates:

$$D_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}$$

where indices $i, j$ run over all $3N$ Cartesian coordinates.

### 2.1 Computing the Hessian with DFTB+

Use the `SecondDerivatives` driver :

```hsd
Driver = SecondDerivatives {
    Delta = 1E-4          # Finite difference step (atomic units)
    Atoms = 1:-1          # All atoms (default)
}

Hamiltonian = DFTB {
    Scc = Yes
    SccTolerance = 1E-7   # MUST be tight for accurate frequencies!
    ...
}
```

**Critical**: For accurate vibrational frequencies, set `SccTolerance = 1E-7` or tighter. The default tolerance can yield unphysical modes .

This produces `hessian.out` in the working directory.

### 2.2 File Format of `hessian.out`

The file contains a flattened $3N \times 3N$ matrix stored in the following order :

```
∂²E/∂x₁∂x₁  ∂²E/∂y₁∂x₁  ∂²E/∂z₁∂x₁  ∂²E/∂x₂∂x₁  ∂²E/∂y₂∂x₁  ...
```

That is, for each atom $i$ and coordinate $\alpha \in \{x,y,z\}$, the row contains derivatives with respect to all atoms $j$ and coordinates $\beta \in \{x,y,z\}$. The matrix is symmetric but stored in full (not triangular) format .

For a system with $N$ atoms, the file contains $(3N)^2$ floating-point numbers, typically in scientific notation.

### 2.3 Python Reader for `hessian.out`

```python
import numpy as np

def read_hessian(filename='hessian.out', n_atoms=None):
    """
    Read DFTB+ hessian.out file.
    
    Parameters:
        filename: path to hessian.out
        n_atoms: number of atoms (if None, inferred from file size)
    
    Returns:
        hessian: (3N, 3N) numpy array in atomic units (Hartree/Bohr²)
    """
    data = np.loadtxt(filename)
    
    if n_atoms is None:
        # Infer from data size: should be (3N)^2 elements
        n_total = int(np.sqrt(len(data)))
        if n_total * n_total != len(data):
            raise ValueError(f"Cannot infer n_atoms from {len(data)} elements")
        n_atoms = n_total // 3
    
    expected = (3 * n_atoms) ** 2
    if len(data) != expected:
        raise ValueError(f"Expected {expected} elements for {n_atoms} atoms, got {len(data)}")
    
    hessian = data.reshape((3 * n_atoms, 3 * n_atoms))
    return hessian

def hessian_to_modes(hessian, masses):
    """
    Convert Hessian to mass-weighted form and diagonalize.
    
    Parameters:
        hessian: (3N, 3N) array in Hartree/Bohr²
        masses: (N,) array in atomic mass units
    
    Returns:
        frequencies: (3N,) vibrational frequencies in cm^-1
        modes: (3N, 3N) normal mode eigenvectors
    """
    n_atoms = len(masses)
    assert hessian.shape == (3*n_atoms, 3*n_atoms)
    
    # Mass-weighting: H_mw = H / sqrt(m_i * m_j)
    mass_vec = np.repeat(masses, 3)
    mass_matrix = np.sqrt(np.outer(mass_vec, mass_vec))
    hessian_mw = hessian / mass_matrix
    
    # Convert to SI for frequency calculation
    # 1 Hartree/Bohr²/amu = (Eh/a0²) / (1.6605e-27 kg)
    # ω² = eigenvalue / mass
    # Use ASE units or explicit conversion
    from ase.units import Hartree, Bohr, _amu
    hessian_si = hessian_mw * Hartree / (Bohr**2) / _amu  # J/kg = 1/s²
    
    # Diagonalize
    eigvals, eigvecs = np.linalg.eigh(hessian_si)
    
    # Convert to cm^-1
    # ω [cm^-1] = sqrt(λ) / (2πc) * 100
    from ase.units import invcm
    frequencies = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) / (2*np.pi) / invcm
    
    return frequencies, eigvecs

# Example usage
hessian = read_hessian('hessian.out', n_atoms=3)  # H2O
print(f"Hessian shape: {hessian.shape}")

# Masses in amu
masses = np.array([15.999, 1.008, 1.008])
freqs, modes = hessian_to_modes(hessian, masses)
print(f"Frequencies (cm^-1): {freqs}")
```

### 2.4 Using the `modes` Code (Post-Processing)

DFTB+ provides a separate `modes` executable for diagonalizing the Hessian and producing animated mode files :

```hsd
# modes_in.hsd
Geometry = GenFormat {
   <<< "geom.gen"
}

DisplayModes = {
  PlotModes = -20:-1
  Animate = Yes
}

SlaterKosterFiles = Type2FileNames {
  Prefix = "../../slakos/mio-ext/"
  Separator = "-"
  Suffix = ".skf"
}

Hessian = {
  <<< "hessian.out"
}

InputVersion = 3
```

Run: `modes`

Output: `modes.xyz` containing animated displacements for each mode.

### 2.5 ASE-Based Hessian Calculation (Alternative)

If you prefer pure Python without `SecondDerivatives`, ASE can compute the Hessian via finite differences:

```python
from ase.calculators.dftb import Dftb
from ase.vibrations import Vibrations
from ase.io import read

atoms = read('geo_end.gen')
atoms.calc = Dftb(
    label='hess',
    Hamiltonian_SCC='Yes',
    Hamiltonian_SccTolerance='1E-7',
    Hamiltonian_MaxAngularMomentum_='',
    Hamiltonian_MaxAngularMomentum_O='"p"',
    Hamiltonian_MaxAngularMomentum_H='"s"',
    ParserOptions_ParserVersion=14
)

vib = Vibrations(atoms)
vib.run()
vib.summary()

# Get Hessian matrix
hessian_ase = vib.get_vibrations().get_hessian_2d()
```

**Note**: ASE's `Vibrations` uses its own finite-difference scheme. The DFTB+ native `SecondDerivatives` is generally preferred for consistency with the `modes` code.

## Summary Table

| Property | Tool/Method | Output File | Python Access |
|----------|-------------|-------------|---------------|
| **Charge density** | `waveplot` | `wp-abs2.cube` | `ase.io.cube.read_cube_data` |
| **MO wavefunctions** | `waveplot` | `wp-K-S-L-real.cube` | `ase.io.cube.read_cube_data` |
| **ESP (grid)** | `Analysis/ElectrostaticPotential` | `ESP.dat` | Custom parser (see above) |
| **ESP (transport)** | `Poisson/SavePotential` | `potential.dat` | `makecube` → cube file |
| **Hessian** | `Driver = SecondDerivatives{}` | `hessian.out` | `np.loadtxt` + reshape |
| **Vibrational modes** | `modes` executable | `modes.xyz` | ASE or custom parser |

## References

- [Waveplot Tutorial](https://dftbplus-recipes.readthedocs.io/en/stable/basics/waveplot.html) 
- [Python Interface Documentation](https://dftbplus-recipes.readthedocs.io/en/stable/interfaces/pyapi/pyapi.html) 
- [DFTB+ v19.1 User Manual (ESP format, Hessian format)](http://gr.xjtu.edu.cn/c/document_library/get_file?folderId=2704407&name=DLFE-129828.pdf) 
- [Vibrational Modes Recipe](https://dftbplus-recipes.readthedocs.io/en/latest/moleculardynamics/startinggeometry.html) 
- [Recent DFTB+ Developments (ASI, transport)](https://pubs.acs.org/doi/10.1021/acs.jpca.5c01146) 

---

# Gemini-3-flash

To obtain real-space quantities and second-order derivatives in DFTB+, you typically move from \"calculation\" mode to \"analysis\" mode. While ASE is great for driving the geometry, the specific task of grid projection and Hessian construction often requires leveraging the DFTB+ native tools or specific Python post-processing.


## 1. Real-Space Projections (Density, Orbitals, & ESP)
DFTB+ uses a basis of Slater-type or Gaussian-type orbitals. To see these in real space (on a grid), you must use the `waveplot` utility or the `dpy` Python library.

### A. The `waveplot` Method
This is the most robust way to project onto a grid. You need the `tag` file (or `detailed.xml`) and the `eigenvec.bin` from your SCC calculation.
1.  **Configure `waveplot_in.hsd`**:
    ```hsd
    Options {
      TotalChargeDensity = Yes
      TotalPotential = Yes
      PlottedOrbitals = 1..10  # Indices of MOs you want
    }
    Grid {
      # Define your real-space box
      Resolution = 0.1
    }
    ```
2.  **Run**: `waveplot` (executable). This generates `.vmd` or `.cube` files.
3.  **Python Integration**: Use `ase.io.cube.read_cube` to pull these into Python as NumPy arrays for custom analysis.

### B. Pythonic Way: `dpy` (DFTB+ Python Tools)
The `dftbplus` package (available via `pip install dftbplus`) allows you to manipulate the wavefunctions directly.
```python
import dftbplus
# Load the output data
calc_data = dftbplus.data.Results(\"detailed.xml\")
# Use the Wavefunction object to evaluate at specific [x, y, z] points
# Note: This requires the Slater-Koster files to be accessible for the basis definitions.
```

## 2. Obtaining the Hessian (Dynamical Matrix)
The Hessian $D_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}$ can be obtained in two ways: analytically (if the specific DFTB+ version/model supports it) or via finite differences.

### A. The \"Native\" Way (Finite Differences)
DFTB+ has a built-in driver for this. In your `dftb_in.hsd`:
```hsd
Driver = SecondDerivatives {
  Delta = 1e-4 # Displacement size
}
```
*   **Output**: This produces a `hessian.out` file containing the $3N \times 3N$ matrix.
*   **Calculation**: It performs $6N$ single-point evaluations (displacing each atom $\pm \Delta$ in $x, y, z$).

### B. The Python/ASE Way (Flexible)
If you want the Hessian as a NumPy array directly within a script, use ASE’s `Vibrations` module. This is numerically identical to the native driver but gives you immediate programmatic access.
```python
from ase.vibrations import Vibrations
from ase.calculators.dftb import Dftb

atoms = ... # Your optimized structure
calc = Dftb(label='system', Hamiltonian_SCC='Yes', ...)
atoms.calc = calc

# This will run the displacements
vib = Vibrations(atoms)
vib.run()

# Access the Hessian (Mass-weighted or raw)
hessian = vib.get_hessian2d() # Returns 3N x 3N numpy array
```

## 3. Working with the Density Matrix Programmatically
If you need the **Density Matrix ($P$)** for a hacker-level integration (e.g., building a custom Hamiltonian or ML descriptor), you should extract it from `detailed.xml`.


```python
import xml.etree.ElementTree as ET
import numpy as np

def get_density_matrix(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Locate the density matrix node
    # Note: Structure depends on whether it's spin-polarized
    dm_node = root.find(\".//density_matrix/real_part\")
    dim = int(root.find(\".//basis_size\").text)
    
    data = np.fromstring(dm_node.text, sep=' ')
    return data.reshape((dim, dim))

P = get_density_matrix(\"detailed.xml\")
```

### Summary of Extraction 

ToolsTargetPreferred ToolFormatGrid Projectionwaveplot.cube (Gaussian Cube)Hessianase.vibrationsnp.arrayESP / Potentialswaveplot.cubeInternal Matricesdetailed.xmlXML ParsingPro Tip: If you are doing real-time hacking on the energy functional, look into the DFTB+ JSON output option (if compiled with JSON support), as it is significantly faster to parse than the heavy XML files.

---

# USER

can we find sourve code of waveplot somewhere e.g. on github so we can reimplement it ourself / itegrate to python?

---

# Gemini-3-flash

Yes, the source code for `waveplot` is part of the main **DFTB+ repository** on GitHub. Since DFTB+ is open-source (LGPL), you can dissect the logic to reimplement the grid projection or integrate it more tightly with your Python stack.

### 1. Where to find the Source Code
The `waveplot` source is typically located in the `tools/waveplot/` directory of the official repository.

*   **GitHub Link:** [DFTB+ Repository - Waveplot Tools](https://github.com/dftbplus/dftbplus/tree/main/tools/waveplot)
*   **Core Logic:** Look for `.f90` files like `waveplot.f90` and `grid.f90`.
*   **The Math:** `waveplot` works by reading the **Slater-Koster (SK) files**, which contain the radial wavefunctions, and the `eigenvec.bin` file, which contains the expansion coefficients $c_{\mu i}$. The orbital $\psi_i$ at a grid point $\mathbf{r}$ is calculated as:
$$\psi_i(\mathbf{r}) = \sum_{\mu} c_{\mu i} \phi_{\mu}(\mathbf{r} - \mathbf{R}_A)$$

where $\phi_{\mu}$ are the basis functions reconstructed from the radial tables in the SK files.

### 2. The \"Hacker\" Way: Bypassing Fortran with Python
If you want to reimplement this in Python without writing raw grid-interpolation code from scratch, there are two high-level paths:

#### A. Use the `dftbplus` Python Package (Official)
The developers provide a Python API that can handle some of these tasks. It is often used for post-processing and can be found in the `python/` folder of the repo or installed via `pip`.
```bash
pip install dftbplus
```
This package includes a `subpackage` for handling Slater-Koster files, allowing you to evaluate the radial functions directly in Python.

#### B. The `sk-gen` and Basis Reconstruction
To get the density or orbitals on a grid, you need the **basis functions**. In DFTB, these are not usually simple Gaussians; they are numerical solutions to the atomic Schrödinger equation.
1.  **Extract Radial Functions:** The radial parts are stored in the `.skf` files or separate `.wfs` files (depending on the parameter set).
2.  **Spline Interpolation:** Use `scipy.interpolate.CubicSpline` to load the radial data from the SK files so you can query the value of $\phi(r)$ at any distance $r$.

### 3. Integrating with ASE (The `Dftb` Calculator)
ASE actually has a \"recipe\" for extracting the density matrix and potentially working with wavefunctions, but for real-space grids, ASE usually expects the calculator to provide a `get_pseudo_density()` or `get_wave_function()` method.

If you are developing a deep integration, you can \"hack\" the ASE `Dftb` calculator to read the binary `eigenvec.bin` directly using `numpy.fromfile`.

### 4. Direct Hessian Access for Developers
Since you asked about the Hessian ($D_{ij}$) earlier, if you are looking at the source code, check `src/dftbplus/driver/secondderivatives.f90`. 

If you want to avoid the filesystem entirely (no `hessian.out`), you can use the **DFTB+ JSON interface**. If you compile DFTB+ with JSON support, you can pipe the output directly into a Python dictionary:
```python
import json
import subprocess

# Run DFTB+ with JSON output
result = subprocess.run(['dftb+', '--output-format', 'json'], capture_output=True)
data = json.loads(result.stdout)

# Access the Hessian directly from the JSON tree
hessian = data['results']['hessian']
```

### Recommendation for Implementation
Reimplementing the full `waveplot` logic is non-trivial because you have to handle the **angular part** (Spherical Harmonics) and the **radial part** (Splines from SK files). 

**The most efficient \"Hacker\" workflow is:**
1. Use DFTB+ to generate the `detailed.xml` and `eigenvec.bin`.
2. Use the `dftbplus` Python library to parse the SK files.
3. Use **NumPy/CuPy** to perform the basis set expansion $\sum c \phi$ over your grid for massive speedups compared to the serial Fortran `waveplot`.

---

# Kimi 2.6

Yes — **waveplot is part of the main DFTB+ source tree** on GitHub and is fully open-source under the LGPL. You can find it, study it, and reimplement or wrap it in Python.

---

## Where to Find the Source Code

### Official Repository

The DFTB+ project is hosted at:

**`https://github.com/dftbplus/dftbplus`** 

Waveplot lives inside the `app/waveplot/` directory. After cloning:

```bash
git clone https://github.com/dftbplus/dftbplus.git
cd dftbplus/app/waveplot
```

Key files you will find there:
- `waveplot.F90` — Main Fortran driver
- `input.F90` — HSD input parsing
- `gridcache.F90` — Grid caching logic for performance
- `wfc.F90` — Wavefunction coefficient handling (reads the `wfc.*.hsd` basis files)
- `slater.F90` — Evaluation of Slater-type orbitals on the grid

### License

DFTB+ is released under the **GNU Lesser General Public License (LGPL)** , so you are free to:
- Read and modify the source
- Reimplement the algorithms in Python
- Integrate it into your own tools (provided you respect LGPL terms for derivative works)

## What Waveplot Actually Does (Algorithm Overview)

To reimplement it in Python, you need to understand the core pipeline:

1. **Read `detailed.xml`** — Contains atom positions, species, lattice vectors, and orbital mapping
2. **Read `eigenvec.bin`** — Binary file with MO coefficients $c_{\mu i}$ (LCAO expansion)
3. **Read `wfc.*.hsd`** — Slater-type orbital exponents $\zeta_{nlm}$ and contraction coefficients for each element
4. **Build a real-space grid** — Typically a cuboid around the molecule
5. **Evaluate basis functions** — For each grid point $\mathbf{r}$, compute all STO basis functions $\phi_\mu(\mathbf{r})$
6. **Construct density or orbitals** —
   - Total density: $\rho(\mathbf{r}) = \sum_i n_i |\psi_i(\mathbf{r})|^2 = \sum_i n_i \left| \sum_\mu c_{\mu i} \phi_\mu(\mathbf{r}) \right|^2$
   - MO density: $|\psi_i(\mathbf{r})|^2$
   - MO wavefunction: $\psi_i(\mathbf{r}) = \sum_\mu c_{\mu i} \phi_\mu(\mathbf{r})$
7. **Write Gaussian cube format**

The STO evaluation is the trickiest part — waveplot uses numerical radial functions combined with spherical harmonics .

## Python Reimplementation Strategy

### Option 1: Direct Fortran-to-Python Wrapper (F2PY / CFFI)

Instead of rewriting the physics, compile waveplot's core routines as a Python-callable library:

```bash
# In the DFTB+ build directory, ensure shared libs are built
cmake -DBUILD_SHARED_LIBS=1 -DWITH_API=1 ..
cmake --build _build
```

Then use `ctypes` or `f2py` to call the grid evaluation routines directly from Python. This is the fastest path if you want the exact same physics.

### Option 2: Pure Python Reimplementation

If you want full control and Python-native integration, you need to implement:

#### A. Reading DFTB+ Output Files

```python
import xml.etree.ElementTree as ET
import numpy as np
import struct

def read_detailed_xml(filename='detailed.xml'):
    """Parse detailed.xml for geometry and orbital info."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    # Atom positions, species, lattice vectors
    # Orbital mapping: which atom -> which basis functions
    # This is all in the XML structure
    ...

def read_eigenvec_bin(filename='eigenvec.bin', n_orbitals=None):
    """Read binary eigenvector file."""
    # Format: typically raw Fortran binary, record-based
    # You may need to experiment with record markers
    with open(filename, 'rb') as f:
        data = f.read()
    # Parse based on n_orbitals and n_states
    ...
```

#### B. Parsing WFC Basis Files

The `wfc.*.hsd` files contain STO definitions. Example structure from `wfc.mio-1-1.hsd` :

```hsd
C {
  # s-orbital
  1 0 0 {
    1.000000  # contraction coefficient
    1.200000  # zeta (exponent)
  }
  # p-orbital
  2 1 0 {
    1.000000
    1.150000
  }
  ...
}
```

You need a parser for this HSD subset and a radial STO evaluator:

```python
import numpy as np
from scipy.special import factorial, factorial2

def slater_radial(r, n, zeta):
    """Normalized Slater-type radial function."""
    N = (2 * zeta)**n * np.sqrt(2 * zeta / factorial(2 * n))
    return N * r**(n - 1) * np.exp(-zeta * r)

def real_spherical_harmonic(l, m, theta, phi):
    """Real spherical harmonics Y_lm."""
    # Use scipy.special.sph_harm for complex, then convert to real
    from scipy.special import sph_harm
    ...
```

#### C. Grid Evaluation Engine

```python
class DFTBGridEvaluator:
    def __init__(self, detailed_xml, eigenvec_bin, wfc_file):
        self.atoms = read_detailed_xml(detailed_xml)
        self.eigenvecs = read_eigenvec_bin(eigenvec_bin)
        self.basis = parse_wfc(wfc_file)
        
    def evaluate_basis(self, r_grid):
        """Evaluate all basis functions on a grid."""
        # r_grid: (nx, ny, nz, 3) array of points
        n_points = np.prod(r_grid.shape[:3])
        n_basis = len(self.basis)
        
        phi = np.zeros((n_points, n_basis))
        
        for i_atom, atom in enumerate(self.atoms):
            r_vec = r_grid.reshape(-1, 3) - atom.position
            r = np.linalg.norm(r_vec, axis=1)
            
            for orb in atom.orbitals:
                # Evaluate STO at all points
                radial = slater_radial(r, orb.n, orb.zeta)
                angular = real_spherical_harmonic(orb.l, orb.m, r_vec)
                phi[:, orb.global_index] = radial * angular
                
        return phi
    
    def get_mo(self, state_idx, spin=0, kpoint=0):
        """Get MO coefficients for a specific state."""
        return self.eigenvecs[spin][kpoint][:, state_idx]
    
    def compute_density_grid(self, r_grid, states=None):
        """Compute total or state-resolved charge density."""
        phi = self.evaluate_basis(r_grid)  # (n_points, n_basis)
        
        if states is None:
            # Total density: sum over occupied states
            density = np.zeros(phi.shape[0])
            for state_idx, occ in enumerate(self.occupations):
                coeffs = self.get_mo(state_idx)
                psi = phi @ coeffs  # (n_points,)
                density += occ * np.abs(psi)**2
        else:
            # Single state density
            coeffs = self.get_mo(states)
            psi = phi @ coeffs
            density = np.abs(psi)**2
            
        return density.reshape(r_grid.shape[:3])
```

#### D. Writing Cube Files

```python
def write_cube(filename, data, atoms, origin, vectors):
    """Write Gaussian cube format."""
    nx, ny, nz = data.shape
    
    with open(filename, 'w') as f:
        f.write("DFTB+ density grid\n")
        f.write("Generated by custom Python waveplot\n")
        
        # Origin and grid vectors
        f.write(f"{len(atoms):5d} {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        f.write(f"{nx:5d} {vectors[0][0]:12.6f} {vectors[0][1]:12.6f} {vectors[0][2]:12.6f}\n")
        f.write(f"{ny:5d} {vectors[1][0]:12.6f} {vectors[1][1]:12.6f} {vectors[1][2]:12.6f}\n")
        f.write(f"{nz:5d} {vectors[2][0]:12.6f} {vectors[2][1]:12.6f} {vectors[2][2]:12.6f}\n")
        
        # Atom positions
        for atom in atoms:
            f.write(f"{atom.atomic_number:5d} {atom.charge:12.6f} "
                   f"{atom.position[0]:12.6f} {atom.position[1]:12.6f} {atom.position[2]:12.6f}\n")
        
        # Volumetric data
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f"{data[ix,iy,iz]:14.5e}")
                    if (iz + 1) % 6 == 0:
                        f.write("\n")
                if nz % 6 != 0:
                    f.write("\n")
```

## Option 3: Hybrid Approach (Recommended)

The most practical path for a "power user/developer" is:

1. **Use DFTB+ Python API** to run calculations and extract H, S, and density matrix via ASI 
2. **Use `pyscf` or custom code** for real-space grid evaluation, since pyscf already has efficient STO/grid machinery
3. **Wrap waveplot as a subprocess** for quick validation, but build your own grid engine for production workflows

This gives you:
- Full control over the grid (adaptive resolution, arbitrary regions)
- Integration with ML workflows (differentiable grids, GPU acceleration via JAX/PyTorch)
- No Fortran compilation dependencies for end-users

## Key References

| Resource | Location | Purpose |
|----------|----------|---------|
| Main repo | `github.com/dftbplus/dftbplus`  | Full source including waveplot |
| Waveplot source | `app/waveplot/` in repo  | Fortran implementation |
| Waveplot docs | DFTB+ Recipes  | Input/output reference |
| WFC basis files | `dftb.org/parameters/download`  | STO coefficients needed for reimplementation |
| ASI interface | Release 25.1 changelog  | For extracting H/S/DM programmatically |

## Quick Start: Clone and Inspect

```bash
# Clone DFTB+ source
git clone https://github.com/dftbplus/dftbplus.git
cd dftbplus

# Find waveplot source
ls app/waveplot/
# Output: CMakeLists.txt  gridcache.F90  input.F90  slater.F90  
#         waveplot.F90  wfc.F90  ...

# Read the main driver
less app/waveplot/waveplot.F90

# Check recent changes (MPI parallelization added in 25.1)
git log --oneline app/waveplot/
```

The Fortran source is well-structured and commented — `waveplot.F90` is only a few hundred lines and clearly shows the high-level algorithm flow.

---

# USER

What about Hamiltonian and Greens function. If I want to compute current flowing over moleculer using Caroli formula (Landaure-Buttinekr) or if I want to compute STM images (using Bardeen, Tersoff-Hamman etc), it would be good to know Green's function G(E)=(EI-H)^(-1)

can we easily get H , ideally as numpy array (without the overhead of files / disk I/O) and eventually are therer some packages like Transiesta which direcly compute electron transporrt using DFTB+ as bakcedn (hamiltonian provider)

---

# Gemini-3-flash

Accessing the Hamiltonian and Green's functions directly in Python is the \"Holy Grail\" for building custom transport (Landauer-Büttiker) or STM (Tersoff-Hamann) simulators.

### 1. Direct Hamiltonian Access (Numpy)
To get the Hamiltonian $H$ and Overlap $S$ matrices into a Numpy array without slow disk I/O, you have two primary \"hacker\" routes:

#### A. The Python API (`libdftbplus`)
The most modern way is using the **Python API** provided in the DFTB+ source. It requires compiling DFTB+ as a shared library.
```python
import dftbplus
import numpy as np

# Initialize the library with your input file
calc = dftbplus.DftbPlus(libpath='/path/to/libdftbplus.so', hsdpath='dftb_in.hsd')

# After a single-point SCC run:
# Extract matrices directly from memory
h_matrix = calc.get_hamiltonian()  # Returns a dense or sparse-like object
s_matrix = calc.get_overlap()

# Convert to standard numpy if needed
h_numpy = np.array(h_matrix)
```

#### B. The JSON Stream (Low Overhead)
If you don't want to deal with shared library compilation, use the **JSON output** flag. While it technically uses a pipe, it avoids the messy parsing of `.dat` files.
```bash
dftb+ --output-format json > results.json
```
In Python, `json.load()` gives you the full Hamiltonian and Density Matrix as nested lists which you can immediately cast to `np.array()`.

### 2. Computing the Green's Function $G(E)$
Once you have $H$ and $S$ as Numpy arrays, the retarded Green's function for a closed system is:
$$G(E) = \left[ (E + i\eta)S - H \right]^{-1}$$

In Python:
```python
def get_greens_function(E, H, S, eta=1e-5):
    identity_size = H.shape[0]
    # G = inverse of [(E + i*eta)*S - H]
    z = E + 1j * eta
    return np.linalg.inv(z * S - H)
```

### 3. Transport Packages (The \"Transiesta\" for DFTB+)
You don't need to write everything from scratch. There are specialized packages that use DFTB+ as the electronic structure \"engine\" for transport:

1.  **DFTB+ Native Transport (NEGF):**
    DFTB+ actually has a built-in **NEGF (Non-Equilibrium Green's Function)** module. It can handle semi-infinite leads and compute transmission $T(E)$ natively.
    *   **Workflow:** Lead calculation $\rightarrow$ Device calculation.
    *   **Keyword:** `Transport { ... }` block in `dftb_in.hsd`.
2.  **SMEAGOL:**
    A world-class transport code that has a long-standing interface with DFTB+. It is effectively the \"Transiesta equivalent\" for the DFTB community.
3.  **DFTBephy:**
    As found in your source interest, this is excellent for **Electron-Phonon coupling** and transport properties (mobility/conductivity), specifically using DFTB+ as the backend.
4.  **NEMO5:**
    A nanoelectronics modeling tool that supports DFTB+ Hamiltonians for massive-scale device simulations.

### 4. STM Imaging (Bardeen/Tersoff-Hamann)
For STM, you need the **Local Density of States (LDOS)** at the tip position $\mathbf{r}_0$:
$$\text{LDOS}(E, \mathbf{r}_0) = \sum_i |\psi_i(\mathbf{r}_0)|^2 \delta(E - \epsilon_i)$$

Since you can now extract $H$ and $S$, you can solve the generalized eigenvalue problem $H C = S C \epsilon$ in Python to get eigenvectors $C$. Use the `waveplot` logic mentioned previously to project these onto a 2D plane above your molecule to generate the STM map.

**Expert Tip:** If you are doing real-space STM, look into the **`pyChemia`** or **`sisl`** Python packages. `sisl` in particular is designed to handle Hamiltonians and Green's functions from DFT codes and has great support for transport and STM analysis.


---

# Kimi 2.6

Here is a comprehensive markdown document addressing your questions about extracting the Hamiltonian as numpy arrays, computing Green's functions for transport/STM, and using DFTB+ with transport packages like libNEGF.

## 1. Getting H, S, and DM as NumPy Arrays (Zero-Copy via ASI)

The **Atomic Simulation Interface (ASI)** is the modern, zero-copy pathway to extract dense Hamiltonian (H), overlap (S), and density matrix (DM) directly into Python as NumPy arrays — no file I/O overhead .

### Prerequisites

Compile DFTB+ with ASI support:
```bash
cmake -DWITH_ASI=1 -DBUILD_SHARED_LIBS=1 -DWITH_API=1 ..
```

Install the Python wrapper:
```bash
pip install asi4py
```

### Python Example: Extracting H, S, DM

```python
import os
import numpy as np
from ctypes import CDLL, RTLD_GLOBAL
from ase.build import molecule
from asi4py.asecalc import ASI_ASE_calculator

ASI_LIB_PATH = os.environ['ASI_LIB_PATH']  # e.g., /path/to/libdftbplus.so
asilib = CDLL(ASI_LIB_PATH, mode=RTLD_GLOBAL)

def init_dftb(asi):
    from ase.calculators.dftb import Dftb
    calc = Dftb(
        label='cluster',
        Hamiltonian_SCC='Yes',
        Hamiltonian_MaxAngularMomentum_='',
        Hamiltonian_MaxAngularMomentum_O='"p"',
        Hamiltonian_MaxAngularMomentum_H='"s"',
        Hamiltonian_ASI_='',
        Hamiltonian_ASI_AsiModifiesModel='No',
        ParserOptions_ParserVersion=14
    )
    calc.write_input(asi.atoms)
    return calc

atoms = molecule('H2O')
atoms.calc = ASI_ASE_calculator(ASI_LIB_PATH, init_dftb, None, atoms)

# Request matrix storage BEFORE running
atoms.calc.asi.keep_density_matrix = True
atoms.calc.asi.keep_hamiltonian = True
atoms.calc.asi.keep_overlap = True

# Run calculation
energy = atoms.get_potential_energy()

# Extract matrices as NumPy arrays
# Index (1, 1) = first k-point, first spin channel
S = atoms.calc.asi.overlap_storage[(1, 1)]      # (n_basis, n_basis)
H = atoms.calc.asi.hamiltonian_storage[(1, 1)]   # (n_basis, n_basis)
DM = atoms.calc.asi.dm_storage.get((1, 1), None) # (n_basis, n_basis)

print(f"Basis size: {atoms.calc.asi.n_basis}")
print(f"Electrons: {np.sum(S * DM.T):.6f}")
print(f"Band energy: {np.sum(H * DM.T):.6f}")
```

**Key properties** :
- `modify_hamiltonian` — You can add an external potential matrix to H before diagonalization
- `register_hamiltonian_callback()` — Register a Python function called when H is built
- `register_dm_callback()` — Register a callback when DM is updated

## 2. Green's Function and Transport in DFTB+

DFTB+ has built-in NEGF transport capabilities via **libNEGF**, an LGPL library that handles non-equilibrium Green's function calculations .

### 2.1 Built-in NEGF Transport (Landauer-Caroli Formula)

DFTB+ can compute transmission and current using the Caroli (Fisher-Lee) formula directly :

$$T_{\beta\alpha}(E) = \text{Tr}[G^r \Gamma_\beta G^a \Gamma_\alpha]$$

$$I = \frac{2e}{h} \int T(E) [f_L(E) - f_R(E)] dE$$

#### Input Setup (`dftb_in.hsd`)

**Non-SCC transmission** (fastest, for testing):
```hsd
Hamiltonian = DFTB {
  SCC = No
  SlaterKosterFiles = Type2FileNames {
    Prefix = "../../slako/"
    Separator = "-"
    Suffix = ".skf"
  }
  Solver = TransportOnly{}  # Skip diagonalization, go straight to NEGF
}

Transport {
  Device {
    AtomRange = 1 8
  }
  Contact {
    Id = "source"
    AtomRange = 9 24
    FermiLevel [eV] = -8.4123
    Potential = 0.0
  }
  Contact {
    Id = "drain"
    AtomRange = 25 40
    FermiLevel [eV] = -8.4123
    Potential = 1.0
  }
  Task = UploadContacts
}

Analysis {
  TunnelingAndDos {
    Verbosity = 101
    EnergyRange [eV] = -6.0 -3.0
    EnergyStep [eV] = 0.01
  }
}
```

**SCC NEGF** (full self-consistent):
```hsd
Hamiltonian = DFTB {
  SCC = Yes
  SccTolerance = 1.0e-7
  Solver = GreensFunction {
    FirstLayerAtoms = 1 61 92 145
    Delta [eV] = 1E-4
    ContourPoints = 20 20
    RealAxisPoints = 55
    LowestEnergy [eV] = -60.0
    FermiCutoff = 10
    EnclosedPoles = 3
    SaveSurfaceGFs = Yes
  }
  Mixer = Broyden {
    MixingParameter = 0.02  # Low mixing for NEGF stability
  }
}

Transport {
  # ... same device/contact blocks ...
}

Analysis {
  TunnelingAndDos {
    EnergyRange [eV] = -6.0 -3.0
    EnergyStep [eV] = 0.01
  }
}
```

**Output**: `transmission.dat` and `lcurrent_u.dat` (local bond currents) .

### 2.2 libNEGF: Standalone NEGF Library

**libNEGF** is a general-purpose NEGF library that can interface with any Hamiltonian provider, including DFTB+ .

**Repository**: `https://github.com/libnegf/libnegf`

**Features**:
- Ballistic and inelastic transport (electron-phonon)
- GPU acceleration (CUDA/cuBLAS)
- Python API via `asi4py` or direct C/Fortran bindings
- Handles systems with >10⁶ atoms 

**Codes interfaced**:
- DFTB+ (semi-empirical DFT)
- TiberCAD (FEM/ETB)
- bigDFT (wavelet DFT, in progress)
- Materials Studio/Biovia DMol3 (commercial) 

### 2.3 neXGf: Exascale NEGF Demonstrator

The **neXGf** project (EoCoE-II) provides a modified DFTB+ (v22.1-based) with enhanced libNEGF for GPU-accelerated transport :

```bash
git clone https://gitlab.maisondelasimulation.fr/ecoce-ii/code-demonstrators.git
cd code-demonstrators/neXGf
```

Structure:
- `dftplus/` — Modified DFTB+ source
- `external/libnegf/` — NEGF library with GPU support
- `interface/` — Example calculations (elastic/inelastic CPU, GPU)

**GPU benchmarks** show 100–400× speedup for large systems (e.g., 20k basis functions) .

## 3. STM Image Simulation (Tersoff-Hamann / Bardeen)

DFTB+ does not have a built-in STM module, but you can build one using the extracted H/S/DM and wavefunction data.

### 3.1 Tersoff-Hamann Approximation

The simplest approach: STM current is proportional to the local density of states (LDOS) at the tip position :

$$I(\mathbf{r}, V) \propto \int_{E_F}^{E_F+eV} \rho(\mathbf{r}, E) dE$$

Where $\rho(\mathbf{r}, E) = \sum_i |\psi_i(\mathbf{r})|^2 \delta(E - E_i)$.

#### Python Implementation Using Waveplot Output

```python
import numpy as np
from ase.io.cube import read_cube_data

def compute_ldos_stm(cube_files, energies, fermi_level, bias_voltage):
    """
    Compute STM image from waveplot-generated MO cube files.
    
    Parameters:
        cube_files: list of cube file paths for states in bias window
        energies: corresponding eigenvalues (eV)
        fermi_level: Fermi energy (eV)
        bias_voltage: tip bias (V)
    
    Returns:
        stm_image: 3D array of tunneling current
    """
    # Determine energy window
    e_min = fermi_level
    e_max = fermi_level + bias_voltage
    
    # Sum densities of states in window
    stm_image = None
    for cube_file, energy in zip(cube_files, energies):
        if e_min <= energy <= e_max:
            data, atoms = read_cube_data(cube_file)
            if stm_image is None:
                stm_image = np.zeros_like(data)
            stm_image += data  # |ψ|² contribution
    
    return stm_image

# Example: Generate STM from waveplot output
# Assuming wp-1-1-N-abs2.cube files for states N=1..M
import glob
cubes = sorted(glob.glob('wp-1-1-*-abs2.cube'))
# ... map to energies from detailed.xml ...
```

### 3.2 Bardeen Approximation (More Accurate)

For quantitative STM, use the Bardeen tunneling formula with explicit tip states :

$$I = \frac{2\pi e}{\hbar} \sum_{\mu,\nu} |M_{\mu\nu}|^2 [f(E_\mu) - f(E_\nu)] \delta(E_\mu - E_\nu)$$

$$M_{\mu\nu} = \frac{\hbar^2}{2m} \int_S (\psi_\mu^* \nabla \psi_\nu - \psi_\nu \nabla \psi_\mu^*) \cdot d\mathbf{S}$$

Where the integral is over a surface between sample and tip.

#### Implementation Strategy

1. **Extract sample wavefunctions** from DFTB+ via waveplot or ASI callbacks
2. **Model tip wavefunction** as an s-wave Green's function :
   $$G(\mathbf{r}, \mathbf{r}_0, E) = -\frac{m}{2\pi\hbar^2} \frac{e^{-\kappa |\mathbf{r} - \mathbf{r}_0|}}{|\mathbf{r} - \mathbf{r}_0|}$$
   where $\kappa = \sqrt{2m\phi}/\hbar$ and $\phi$ is the work function.
3. **Evaluate Bardeen integral** on a surface above the sample

```python
import numpy as np
from scipy.special import sph_harm

def tip_green_function(r, r_tip, work_function=4.5):
    """
    Spherical s-wave tip Green's function (Tersoff-Hamann).
    Work function in eV.
    """
    # κ in atomic units: sqrt(2 * φ [Hartree])
    phi_hartree = work_function / 27.2114
    kappa = np.sqrt(2 * phi_hartree)
    
    dist = np.linalg.norm(r - r_tip)
    return -0.5 / np.pi * np.exp(-kappa * dist) / dist

def bardeen_matrix_element(sample_psi, sample_grad, tip_psi, tip_grad, surface):
    """
    Compute Bardeen matrix element M_μν over a surface.
    """
    # Surface normal and area elements
    # M = (ℏ²/2m) ∫ (ψ_μ* ∇ψ_ν - ψ_ν ∇ψ_μ*) · dS
    prefactor = 0.5  # atomic units: ℏ²/2m = 0.5
    
    integrand = (np.conj(sample_psi) * tip_grad - tip_psi * np.conj(sample_grad))
    M = prefactor * np.sum(integrand * surface.normals * surface.areas)
    return M
```

### 3.3 Existing Python STM Tools

| Tool | Source | Basis | Notes |
|------|--------|-------|-------|
| **conquest_stm** | [GitHub](https://github.com/johanmcquillan/conquest_stm)  | PAO/LCAO | Bardeen + Tersoff-Hamann for CONQUEST; adaptable to DFTB+ |
| **STM-2DScan** | [GitHub](https://github.com/ShuangLeung/STM_2DScan)  | Plane waves/VASP | CHGCAR-based; could adapt to DFTB+ cube output |
| **cp2k-sumbias.py** | CP2K | GPW/GAPW | Sums bias-window orbitals; similar workflow for DFTB+ |

## 4. Complete Workflow: DFTB+ → Python → Transport/STM

### 4.1 Architecture Overview

```
┌─────────────┐     ASI API      ┌──────────────┐
│   DFTB+     │ ────────────────→│  Python      │
│  (libdftb+) │   H, S, DM       │  (asi4py)    │
└─────────────┘                  └──────┬───────┘
                                        │
                    ┌───────────────────┼───────────────────┐
                    ↓                   ↓                   ↓
              ┌─────────┐        ┌──────────┐        ┌──────────┐
              │ libNEGF │        │  Custom  │        │  STM     │
              │ (NEGF)  │        │  Green's │        │  Engine  │
              │         │        │  Function│        │          │
              └─────────┘        └──────────┘        └──────────┘
```

### 4.2 Custom Green's Function in Python

If you want to compute $G(E) = (EI - H)^{-1}$ yourself using DFTB+'s H and S:

```python
import numpy as np
from numpy.linalg import inv

def compute_green_function(H, S, energy, eta=1e-5):
    """
    Compute retarded Green's function G^r(E) = [(E+iη)S - H]^{-1}
    
    Parameters:
        H: Hamiltonian matrix (n_basis, n_basis)
        S: Overlap matrix (n_basis, n_basis)
        energy: Energy value (Hartree)
        eta: Broadening parameter (Hartree)
    
    Returns:
        G: Complex Green's function matrix
    """
    z = (energy + 1j * eta) * S - H
    G = inv(z)
    return G

# Example: LDOS from Green's function
def compute_ldos(G, S):
    """
    Local density of states: ρ(E) = -1/π Im[diag(G·S)]
    """
    return -1.0 / np.pi * np.imag(np.diag(G @ S))

# Usage with ASI-extracted matrices
E = -0.2  # Hartree (~ -5.4 eV)
G = compute_green_function(H, S, E)
ldos = compute_ldos(G, S)
```

### 4.3 Integrating with libNEGF for Full Transport

For a complete Landauer-Büttiker calculation with DFTB+ as the Hamiltonian provider:

```python
# Pseudocode for libNEGF + DFTB+ via ASI
import libnegf  # hypothetical Python bindings
from asi4py.asecalc import ASI_ASE_calculator

# 1. Run DFTB+ and extract H, S via ASI
atoms.calc = ASI_ASE_calculator(...)
atoms.calc.asi.keep_hamiltonian = True
atoms.calc.asi.keep_overlap = True
energy = atoms.get_potential_energy()
H = atoms.calc.asi.hamiltonian_storage[(1,1)]
S = atoms.calc.asi.overlap_storage[(1,1)]

# 2. Pass to libNEGF
negf = libnegf.NEGF()
negf.set_hamiltonian(H)
negf.set_overlap(S)

# 3. Define contacts (from DFTB+ Transport block)
negf.add_contact(atom_range=(0, 23), fermi_level=-8.4123, potential=0.0)   # source
negf.add_contact(atom_range=(24, 39), fermi_level=-8.4123, potential=1.0)   # drain

# 4. Compute transmission
transmission = negf.calculate_transmission(energy_grid=np.linspace(-6, -3, 300))
current = negf.calculate_current(temperature=300)
```

**Note**: libNEGF's Python API is actively developed but may require building from source. The Fortran/C API is stable and well-tested .

## 5. Hessian and Vibrational Analysis

### 5.1 Overview

DFTB+ provides built-in Hessian computation via the `SecondDerivatives` driver, which computes second derivatives of energy with respect to atomic coordinates using finite differences.

### 5.2 Method: DFTB+ Native SecondDerivatives Driver

**Example script:** `tests/dftb/example_hessian.py`

**CLI usage:**
```bash
python3 tests/dftb/example_hessian.py molecule.xyz
python3 tests/dftb/example_hessian.py molecule.xyz --modes-output modes.xyz --vector-scale 2.0
```

**Output files:**
- `hessian.out` - Raw Hessian matrix (Hartree/Bohr²)
- `modes.xyz` - Vibration modes in Jmol format
- Frequencies printed to console (cm⁻¹)

**Key functions in `pyBall/dftb_utils.py`:**
- `read_hessian()` - Read DFTB+ hessian.out file
- `hessian_hartree_bohr_to_eV_angstrom()` - Unit conversion
- `hessian_to_mass_weighted()` - Convert to dynamical matrix
- `hessian_to_frequencies()` - Compute vibrational frequencies
- `write_vibration_modes_jmol()` - Export modes for Jmol visualization

**CLI options:**
```
  molecule.xyz          XYZ file to load (default: built-in H2O)
  --use-ase             Use ASE to load XYZ (default: AtomicSystem)
  --sk-path PATH        Path to Slater-Koster files
  --delta FLOAT         Finite difference step in atomic units (default: 1e-4)
  --no-freq             Skip frequency calculation
  --workdir DIR         Working directory (default: hessian_calc)
  --modes-output FILE   Output file for vibration modes (default: vibration_modes.xyz)
  --vector-scale FLOAT  Scaling factor for vibration vectors (default: 1.0)
```

### 5.3 Hessian File Format

DFTB+ `hessian.out` contains a flattened (3N × 3N) matrix in Hartree/Bohr²:
```
∂²E/∂x₁∂x₁  ∂²E/∂y₁∂x₁  ∂²E/∂z₁∂x₁  ∂²E/∂x₂∂x₁  ...
```

The file uses formatted output with varying columns per line. The reader in `dftb_utils.py` extracts all numbers and reshapes them.

### 5.4 Unit Conversions

- **DFTB+ output:** Hartree/Bohr²
- **Converted to:** eV/Å²
- **Conversion factor:** 97.207 (27.2114 / 0.529177²)

### 5.5 Vibrational Frequencies

Computed from mass-weighted Hessian:
```
D = M^(-1/2) * H * M^(-1/2)
ω², v = eig(D)
freq = ω / (2πc)  [cm⁻¹]
```

Where:
- M is diagonal mass matrix (repeated 3× per atom)
- c = 2.998×10¹⁰ cm/s
- Frequencies in wavenumbers (cm⁻¹)

### 5.6 Jmol Visualization

Vibration modes are exported in multi-model XYZ format:
```
3
Vibration Frequency: 1600.00 cm^-1
O 0.000 0.000 0.000 0.000 0.000 0.050
H 0.957 0.000 0.000 0.000 0.100 0.000
H -0.240 0.927 0.000 0.000 0.100 0.000
```

Each line: `element x y z vx vy vz`

**Jmol commands:**
```jmol
load "modes.xyz"
vectors on
vibrate on
vectorScale 2.0
```

### 5.7 Molecular Orbitals

See the waveplot section in this document for generating cube files of molecular orbitals.

**Key functions:**
- `run_waveplot()` - Python wrapper for waveplot utility
- `read_cube()` - Read cube files using ASE
- `plot_orbital()` - Plot orbital slices

## 6. Summary and Recommendations

| Goal | Method | Python Access | Performance |
|------|--------|-------------|-------------|
| **Get H, S, DM** | ASI + `asi4py` | Direct NumPy arrays | Zero-copy, fast |
| **Built-in transport** | DFTB+ `Solver = GreensFunction{}` | File I/O (`transmission.dat`) | Native, parallel |
| **Full NEGF control** | libNEGF standalone | C/Fortran API, partial Python | GPU-accelerated |
| **Exascale/GPU transport** | neXGf (libNEGF+DFTB+) | GitLab repo | A100 GPU, 100-400× speedup |
| **STM (Tersoff-Hamann)** | Custom Python + waveplot | ASE cube I/O | Moderate |
| **STM (Bardeen)** | Custom Python + ASI callbacks | NumPy/scipy | Accurate, flexible |

### Recommended Path for Power Users

1. **Immediate**: Use ASI (`asi4py`) to extract H/S/DM as NumPy arrays and build your own Green's function / STM tools 
2. **Production transport**: Use DFTB+'s built-in `GreensFunction` solver for standard NEGF calculations 
3. **Advanced/GPU**: Integrate with libNEGF directly for inelastic scattering, large systems, or GPU acceleration 
4. **STM**: Start with Tersoff-Hamann using waveplot-generated cube files, then implement Bardeen with ASI-extracted wavefunctions for quantitative results 

---

# USER

Can we somehow integrate with these?

and these are some modules I need to install separately ?

**Advanced Features (not yet wrapped in dftb_lib.py):**
- Transport/NEGF functions
- Time-dependent TD-DFTB
- REKS multi-reference
- xTB integration
- CHIMES ML force fields
- PLUMED metadynamics
- D4 dispersion
- Custom eigensolvers (via ASI)


---

# Kimi 2.6
  
Here is a comprehensive tutorial for integrating DFTB+ with these advanced computational chemistry tools through Python interfaces.

# Advanced DFTB+ Integration Tutorial

## 1. PLUMED Metadynamics Integration

PLUMED is a plugin for free-energy calculations that works with DFTB+ for metadynamics, umbrella sampling, and enhanced sampling .

### 1.1 Building DFTB+ with PLUMED Support

**Prerequisites**: Install PLUMED2 first:

```bash
# Option A: From conda (easiest)
conda install -c conda-forge plumed

# Option B: From source (for latest features)
git clone https://github.com/plumed/plumed2.git
cd plumed2
./configure --prefix=$HOME/opt/plumed
make -j$(nproc)
make install
```

**Build DFTB+ with PLUMED**:

```bash
cd ~/git_SW/dftbplus/_build

# Clean and reconfigure with PLUMED
rm -rf *

export FC=gfortran-12
export CC=gcc-12
export CXX=g++-12
export PLUMED_PREFIX=$HOME/opt/plumed  # or conda path

cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt/dftb+ \
    -DWITH_PYTHON=1 \
    -DWITH_API=1 \
    -DWITH_ASI=1 \
    -DENABLE_DYNAMIC_LOADING=1 \
    -DBUILD_SHARED_LIBS=1 \
    -DWITH_PLUMED=1 \
    -DPLUMED_INCLUDE_DIR=$PLUMED_PREFIX/include \
    -DPLUMED_LIBRARY=$PLUMED_PREFIX/lib/libplumed.so \
    -DWITH_TRANSPORT=1 \
    -DWITH_TBLITE=1 \
    -DWITH_ARPACK=1 \
    -DWITH_OMP=1 \
    ..

cmake --build . -- -j$(nproc)
cmake --install .
```

### 1.2 Running Metadynamics with DFTB+

**Input file** (`dftb_in.hsd`):

```hsd
Driver = VelocityVerlet {
    TimeStep [fs] = 1.0
    Plumed = Yes
    Thermostat = NoseHoover {
        Temperature [Kelvin] = 400
        CouplingStrength [cm^-1] = 3050
    }
    Steps = 100000
}

Hamiltonian = DFTB {
    Scc = Yes
    SlaterKosterFiles = Type2FileNames {
        Prefix = "./slakos/mio-ext/"
        Separator = "-"
        Suffix = ".skf"
    }
    MaxAngularMomentum {
        O = "p"
        H = "s"
    }
}
```

**PLUMED input** (`plumed.dat`):

```plumed
# Define collective variables
DISTANCE ATOMS=4,9 LABEL=d1
DISTANCE ATOMS=5,9 LABEL=d2

# Metadynamics bias
METAD ...
LABEL=met
ARG=d1,d2
PACE=100
HEIGHT=3
SIGMA=0.01,0.01
FILE=HILLS
BIASFACTOR=4
TEMP=400
... METAD

PRINT ARG=d1,d2 STRIDE=100 FILE=plumed_o.dat

ENDPLUMED
```

### 1.3 Python Analysis of Metadynamics Results

```python
import numpy as np
import matplotlib.pyplot as plt

# Read HILLS file
hills = np.loadtxt('HILLS', skiprows=1)

# Extract CVs and bias
cv1, cv2, bias = hills[:, 1], hills[:, 2], hills[:, 3]

# Free energy surface (simple sum_hills approximation)
from scipy.ndimage import gaussian_filter

# Create 2D histogram
H, xedges, yedges = np.histogram2d(cv1, cv2, bins=50, 
                                    weights=bias, 
                                    range=[[0, 3], [0, 3]])

# Smooth and plot
F = -gaussian_filter(H, sigma=2)
plt.contourf(F.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
             levels=20, cmap='viridis')
plt.colorbar(label='Free Energy (kJ/mol)')
plt.xlabel('CV1: d(O1-H)')
plt.ylabel('CV2: d(O2-H)')
plt.title('Metadynamics Free Energy Surface')
plt.savefig('fes.png', dpi=150)
```

## 2. xTB via tblite (Python Native)

The `tblite` Python package provides **native** xTB calculations without needing DFTB+ at all. This is the easiest path for xTB .

### 2.1 Installation

```bash
# Using conda/mamba (recommended)
mamba create -n xtb -c conda-forge tblite-python qcelemental
mamba activate xtb

# Or pip
pip install tblite
```

### 2.2 Single-Point Calculation

```python
from tblite.interface import Calculator
import numpy as np

# Caffeine molecule (atomic numbers and positions in Bohr)
elements = np.array([6,7,6,7,6,6,6,8,7,6,8,7,6,6,1,1,1,1,1,1,1,1,1,1])
coordinates = np.array([
    [ 2.028,  0.092, -0.143],
    [ 4.750,  0.024, -0.143],
    # ... (full coordinates from tblite docs)
])

# Create calculator
calc = Calculator("GFN2-xTB", elements, coordinates)

# Run single point
res = calc.singlepoint()

print(f"Energy: {res.get('energy'):.6f} Hartree")
print(f"Dipole: {res.get('dipole')} e·Bohr")
```

### 2.3 Geometry Optimization

```python
from tblite.interface import Calculator
from berny import Berny, geomlib
import numpy as np

# Read geometry
optimizer = Berny(geomlib.readfile("caffeine.xyz"))
geom = next(optimizer)
elements = [symbol for symbol, _ in geom]
coords = np.asarray([coordinate for _, coordinate in geom])

# Setup calculator (convert Angstrom to Bohr)
calc = Calculator("GFN2-xTB", 
                  Calculator.symbols_to_numbers(elements), 
                  coords * 1.889726)

# Optimization loop
calc.set("verbosity", 0)
for geom in optimizer:
    coordinates = np.asarray([coordinate for _, coordinate in geom])
    calc.update(positions=coordinates * 1.889726)
    res = calc.singlepoint()
    energy = res.get("energy")
    gradient = res.get("gradient")
    optimizer.send((energy, gradient / 1.889726))
```

### 2.4 Extracting Hamiltonian and Density Matrix

```python
# After singlepoint calculation
res = calc.singlepoint()

# Get matrices directly as NumPy arrays
H = res.get("hamiltonian-matrix")      # Core Hamiltonian
S = res.get("overlap-matrix")           # Overlap
DM = res.get("density-matrix")          # Density matrix

print(f"H shape: {H.shape}")
print(f"Tr(S·DM) = {np.trace(S @ DM):.4f}")
```

**Note**: The Hamiltonian from `tblite` is the **core Hamiltonian**, not the full converged Hamiltonian. To get the full Hamiltonian, transform orbital energies using coefficients .

## 3. TD-DFTB Excited States

### 3.1 Via DFTB+ Input (Standard Approach)

```hsd
Hamiltonian = DFTB {
    Scc = Yes
    # ... standard DFTB settings ...
}

Analysis {
    ExcitedState {
        StateOfInterest = 1
        Symmetry = Singlet
        NrOfExcitations = 10
        WriteTransitions = Yes
    }
}
```

Output: `excitations.log` with transition energies, oscillator strengths, and dipole moments.

### 3.2 Via MLatom Python Interface

MLatom provides a high-level Python interface to TD-DFTB :

```python
import mlatom as ml

# Setup DFTB method
dftb = ml.models.methods(method="DFTB")

# Ground state
dftb.predict(molecule=mymol, 
             calculate_energy=True, 
             calculate_energy_gradients=True)

# Excited state (10th state out of 20)
dftb.predict(molecule=mymol,
             calculate_energy=True,
             calculate_energy_gradients=True,
             nstates=20,
             current_state=10)
```

### 3.3 Via DFTBaby (Pure Python)

DFTBaby is a pure-Python DFTB implementation with TD-DFTB :

```python
from DFTB.LR_TDDFTB import ExcitedStates

# Setup ground state first
# ... (DFTBaby ground state calculation) ...

# Calculate excited states
excited = ExcitedStates(ground_state, nstates=10)
excited.solve()

print(f"Transition energies: {excited.energies} eV")
print(f"Oscillator strengths: {excited.oscillator_strengths}")
```

## 4. REKS Multi-Reference Calculations

REKS (Restricted Ensemble Kohn-Sham) is built into DFTB+ for ground and low-lying excited states .

### 4.1 DFTB+ Input

```hsd
Hamiltonian = DFTB {
    Scc = Yes
    # ... standard settings ...
}

Analysis {
    REKS = SA-REKS {
        Energy = {
            Functional = {
                "PPS" "OSS"
            }
        }
        StateInteractions = Yes
        Shift = 0.3
    }
}
```

### 4.2 Python Interface

Currently, REKS is **not directly exposed** in the Python API. You must:
1. Write HSD input with REKS block
2. Run via `dftbplus` executable or Python API with `run()` method
3. Parse output files for results

```python
import dftbplus

cdftb = dftbplus.DftbPlus(libpath=LIB_PATH, hsdpath='reks_input.hsd')
cdftb.run()  # Runs calculation, no geometry changes

# Parse results from output files
with open('detailed.out', 'r') as f:
    for line in f:
        if 'REKS' in line and 'Energy' in line:
            print(line.strip())

cdftb.close()
```

## 5. CHIMES Machine Learning Force Fields

CHIMES is supported in DFTB+ for machine learning corrections .

### 5.1 Building with CHIMES

```bash
cmake \
    -DWITH_CHIMES=1 \
    -DCHIMES_INCLUDE_DIR=$HOME/opt/chimes/include \
    -DCHIMES_LIBRARY=$HOME/opt/chimes/lib/libchimes.so \
    ..
```

### 5.2 Input Usage

```hsd
Hamiltonian = DFTB {
    Scc = Yes
    # ... standard DFTB ...
    
    Chimes = Yes {
        ParamFile = "chimes_params.txt"
    }
}
```

### 5.3 Python Integration

CHIMES corrections are applied automatically during DFTB+ calculations. Access via standard Python API:

```python
cdftb = dftbplus.DftbPlus(libpath=LIB_PATH, hsdpath='chimes_input.hsd')
cdftb.set_geometry(coords)
energy = cdftb.get_energy()  # Includes CHIMES correction
forces = cdftb.get_gradients()
```

## 6. Complete Integration Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Python Workflow Layer                    │
├─────────────┬─────────────┬─────────────┬───────────────────┤
│   PLUMED    │   tblite    │  TD-DFTB    │     REKS/CHIMES   │
│  (MetaD)    │   (xTB)     │  (Excited)  │   (Multi-ref/ML)  │
├─────────────┼─────────────┼─────────────┼───────────────────┤
│  ASE/PLUMED │  tblite-py  │  MLatom/    │   DFTB+ HSD       │
│  interface  │  native     │  DFTBaby    │   + file I/O      │
├─────────────┴─────────────┴─────────────┴───────────────────┤
│              DFTB+ libdftbplus.so (your build)               │
│         WITH_PLUMED=1  WITH_TBLITE=1  WITH_CHIMES=1          │
└─────────────────────────────────────────────────────────────┘
```

## Summary Table: Integration Paths

| Feature | Python Access | Build Flag | External Dependency |
|---------|--------------|------------|---------------------|
| **PLUMED** | ASE + PLUMED.dat | `-DWITH_PLUMED=1` | PLUMED2 library |
| **xTB** | `tblite` package (native) | N/A (separate package) | `pip install tblite` |
| **TD-DFTB** | MLatom / DFTBaby | N/A | `pip install mlatom` |
| **REKS** | File I/O only | Built-in | None |
| **CHIMES** | Standard Python API | `-DWITH_CHIMES=1` | CHIMES library |
| **H/S/DM** | `tblite` or ASI | `-DWITH_ASI=1` | ASI library |

## Quick Start Recommendations

1. **For metadynamics**: Build DFTB+ with `-DWITH_PLUMED=1`, use ASE for setup, PLUMED for bias
2. **For xTB**: Skip DFTB+ entirely, use `tblite` Python package directly — it has the best Python API with native NumPy matrix access 
3. **For excited states**: Use MLatom for TD-DFTB or DFTBaby for pure-Python implementation
4. **For multi-reference**: Use DFTB+ REKS via HSD input, parse output files
5. **For ML corrections**: Build with CHIMES, use standard Python API

The `tblite` package is particularly powerful for your use case because it provides **direct NumPy access to Hamiltonian, overlap, and density matrices** without needing ASI or file I/O .

---

# USER

what can I do about this?

Waveplot is failing because it requires wavefunction coefficient files (wfc.*.hsd) that match the Slater-Koster parameter set. The test SK files we're using don't have the corresponding wavefunction coefficient files.

**Issue:** Waveplot needs wfc files that are specific to each Slater-Koster parameter set (mio, 3ob, etc.). The test directory has SK files but not the matching wfc files.

**Solution:** You need to use a complete DFTB+ parameter set that includes both SK files and wavefunction coefficient files. The script found wfc files at:

- `/home/prokophapala/SIMULATIONS/dftbplus/recipes/slakos/wfc/wfc.mio-1-1.hsd`

But these are for the mio parameter set, while we're using the test SK files which may be incompatible.

**To use waveplot properly:**

1. Install a complete DFTB+ parameter set (e.g., mio-1-1) with both SK and wfc files
2. Update the script to use the matching SK path
3. Run with the correct wfc path

The script structure is correct - it's just a file availability issue with the test parameter set. Would you like me to modify the script to use the mio parameter set instead, or document this limitation in the DFTB_docs.md?

---

# Kimi 2.6

The issue is that **waveplot requires wavefunction coefficient files (`wfc.*.hsd`) that match your Slater-Koster parameter set**. The test SK files (`testparams`) do **not** include these WFC files — they are only for testing DFTB+ core functionality, not for visualization .

The WFC files contain the Slater-type orbital exponents and contraction coefficients needed to reconstruct wavefunctions in real space. They are parameter-set-specific and must be downloaded separately .

## Solution: Download Complete Parameter Sets with WFC Files

### Option 1: Use DFTB+ Recipes Script (Recommended)

The DFTB+ recipes repository includes a script that downloads both SK files **and** WFC files:

```bash
# If you have the recipes directory
cd ~/git_SW/dftbplus/recipes
./scripts/get_slakos

# Or download manually from dftb.org
wget https://www.dftb.org/fileadmin/DFTB/public/slako-packed.tar
tar xf slako-packed.tar
```

This gives you complete parameter sets like `mio-1-1`, `3ob-3-1`, etc., each with matching `wfc.*.hsd` files .

### Option 2: Direct Download from dftb.org

```bash
mkdir -p ~/dftb_params
cd ~/dftb_params

# Download packed parameters (includes WFC files)
wget https://www.dftb.org/fileadmin/DFTB/public/slako-packed.tar
tar xf slako-packed.tar

# Or unpacked (larger download, individual files)
wget https://www.dftb.org/fileadmin/DFTB/public/slako-unpacked.tar.xz
tar xf slako-unpacked.tar.xz
```

After extraction, you'll have:
```
slako/
├── mio/
│   ├── mio-1-1/
│   │   ├── C-C.skf
│   │   ├── C-H.skf
│   │   └── ...
│   └── wfc/
│       └── wfc.mio-1-1.hsd   # ← WFC file for waveplot
├── 3ob/
│   ├── 3ob-3-1/
│   └── wfc/
│       └── wfc.3ob-3-1.hsd
└── ...
```

## Updated Python Script Using Complete Parameter Set

```python
import os
import subprocess
from pathlib import Path

# Configuration: use COMPLETE parameter set with WFC files
PARAM_SET = "mio-1-1"  # or "3ob-3-1", "ob2-1-1", etc.
PARAM_BASE = Path.home() / "dftb_params" / "slako"  # Where you extracted slako-packed.tar

SK_PATH = PARAM_BASE / PARAM_SET / PARAM_SET
WFC_PATH = PARAM_BASE / "wfc" / f"wfc.{PARAM_SET}.hsd"

def run_dftb_waveplot(workdir, geometry, sk_path, wfc_path):
    """
    Run DFTB+ calculation and waveplot with proper WFC files.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    
    # Write DFTB+ input
    dftb_input = f'''
Geometry = GenFormat {{
{geometry}
}}

Hamiltonian = DFTB {{
    Scc = Yes
    SccTolerance = 1.0e-5
    SlaterKosterFiles = Type2FileNames {{
        Prefix = "{sk_path}/"
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum = {{
        O = "p"
        H = "s"
    }}
}}

Options {{
    WriteDetailedXml = Yes
}}

Analysis {{
    WriteEigenvectors = Yes
}}

ParserOptions {{
    ParserVersion = 14
}}
'''
    with open(workdir / "dftb_in.hsd", "w") as f:
        f.write(dftb_input)
    
    # Run DFTB+
    subprocess.run(["dftb+"], cwd=workdir, check=True)
    
    # Write waveplot input with WFC file
    waveplot_input = f'''
Options {{
    TotalChargeDensity = Yes
    TotalChargeDifference = Yes
    ChargeDensity = Yes
    RealComponent = Yes
    PlottedSpins = 1 -1
    PlottedLevels = 1:-1
    PlottedRegion = OptimalCuboid {{}}
    NrOfPoints = 50 50 50
    NrOfCachedGrids = -1
    Verbose = Yes
}}

DetailedXml = "detailed.xml"
EigenvecBin = "eigenvec.bin"

Basis {{
    Resolution = 0.01
    <<+ "{wfc_path}"
}}
'''
    with open(workdir / "waveplot_in.hsd", "w") as f:
        f.write(waveplot_input)
    
    # Run waveplot
    subprocess.run(["waveplot"], cwd=workdir, check=True)
    
    # Collect output
    cubes = list(workdir.glob("wp-*.cube"))
    return cubes

# Example: H2O molecule
geometry = '''3  C
 O H
     1    1    0.00000000000E+00  -0.10000000000E+01   0.00000000000E+00
     2    2    0.00000000000E+00   0.00000000000E+00   0.78306400000E+00
     3    2    0.00000000000E+00   0.00000000000E+00  -0.78306400000E+00'''

cubes = run_dftb_waveplot(
    workdir="./h2o_waveplot",
    geometry=geometry,
    sk_path=SK_PATH,
    wfc_path=WFC_PATH
)

print(f"Generated {len(cubes)} cube files:")
for c in cubes:
    print(f"  {c.name}")
```

## Key Points

| Issue | Solution |
|-------|----------|
| **Missing WFC files** | Download complete parameter sets from dftb.org, not just testparams  |
| **WFC must match SK set** | Use `wfc.mio-1-1.hsd` with `mio-1-1` SK files, `wfc.3ob-3-1.hsd` with `3ob-3-1`, etc.  |
| **Test params incomplete** | `testparams` are for regression testing only, lack WFC and may give unphysical results  |
| **WFC location** | Usually in `slako/wfc/` directory alongside the SK parameter folders  |

---

## Quick Fix for Your Existing Setup

If you already have the recipes repository or DFTB+ source:

```bash
# The DFTB+ source includes a utility to get externals
cd ~/git_SW/dftbplus
./utils/get_opt_externals ALL  # Downloads SK + WFC + everything

# This creates external/slakos/ with complete parameter sets
ls external/slakos/
# Should show: mio-1-1/, 3ob-3-1/, wfc/, etc.
```

Then update your script to use `~/git_SW/dftbplus/external/slakos/` as the parameter base path.

---

# Gemini

This is a classic \"mismatch\" hurdle in DFTB+. Because DFTB+ is a semi-empirical method, the radial part of the wavefunction isn't stored in the `.skf` files (which only contain the pre-calculated integrals). To plot densities or orbitals, you need those supplemental `.hsd` files that describe the basis set shape.

Here is how you can resolve this and structure your documentation to save future-you from this headache.

### 1. Fix the Immediate Error: Using the `mio` set
The `mio-1-1` set is the most common \"standard\" set for organic molecules. If your test system consists of atoms supported by `mio` (H, C, N, O, S, P), you should switch to it entirely.

**Modify your `dftb_in.hsd` or script to point here:**
*   **SK Files:** `/home/prokophapala/SIMULATIONS/dftbplus/recipes/slakos/mio-1-1/`
*   **WFC Files:** `/home/prokophapala/SIMULATIONS/dftbplus/recipes/slakos/wfc/`

### 2. Understanding the Waveplot Requirement
In your `waveplot_in.hsd`, the `Basis` block must point to the directory containing those `wfc.*.hsd` files. 

```hsd
Basis {
  Resolution = 0.01
  # This folder MUST contain wfc.H.hsd, wfc.C.hsd, etc.
  ReadWfc {
    InputPrefix = \"/home/prokophapala/SIMULATIONS/dftbplus/recipes/slakos/wfc/wfc\"
  }
}
```

### 3. How to \"Hack\" or Find missing WFC files
If you are using a parameter set (like `3ob` or `pbc`) and can't find the `wfc` files:
1.  **Check [dftb.org](https://dftb.org/parameters):** When you download a parameter set, look for a \"Wavefunction coefficients\" or \"Additional files\" download link.
2.  **The `mio` fallback:** Many parameter sets use the same confinement potentials as `mio`. In a pinch, developers often use `mio` wfc files for other sets, though this technically introduces a small inconsistency in the visual representation of the density.

### 4. Updating your `DFTB_docs.md`
I recommend documenting this under a **\"Troubleshooting & Basis Sets\"** section. Here is a snippet you can copy-paste:

## Troubleshooting: Waveplot & Basis Functions
**Problem:** `Waveplot` fails with \"Missing wavefunction coefficient files.\"
**Cause:** DFTB+ requires specific radial basis function descriptions (`wfc.Element.hsd`) to project data onto a real-space grid. These are separate from the `.skf` files.

**The Fix:**
*   Ensure the `Basis` block in `waveplot_in.hsd` points to the correct prefix.
*   **Standard Path:** Always maintain a `slakos/wfc/` directory in your environment.
*   **Mapping:** 
    *   `mio-1-1` parameters $\rightarrow$ use `wfc.mio-1-1` files.
    *   `3ob` parameters $\rightarrow$ use `wfc.3ob-1-1` files.

### 5. Strategy: Modify Script vs. Documenting
**I recommend modifying the script to use `mio-1-1` for your current tests.** 
Why? Because using \"test\" SK files is only useful for verifying the code runs; for any meaningful density or orbital analysis, you need the physical consistency of a complete set like `mio`.

---

# Tutorial: Solving Waveplot Missing WFC Files

## Problem
Waveplot fails with error: "Missing wavefunction coefficient files" or parsing error when trying to read `wfc.*.hsd` files.

## Root Cause
DFTB+ requires wavefunction coefficient files (`wfc.*.hsd`) to reconstruct wavefunctions in real space. These are separate from the Slater-Koster `.skf` files and must match the parameter set being used.

## Solution Steps

### Step 1: Download Complete Parameter Set

```bash
# Download complete parameter set from dftb.org
mkdir -p ~/dftb_params
cd ~/dftb_params
wget https://www.dftb.org/fileadmin/DFTB/public/slako-packed.tar
tar xf slako-packed.tar
```

This creates a directory structure like:
```
slako/
├── mio-1-1/          # SK files
├── 3ob-3-1/          # SK files
└── wfc/              # Wavefunction coefficient files
    ├── wfc.mio-1-1.hsd
    └── wfc.3ob-3-1.hsd
```

### Step 2: Update Script Paths

Modify `example_orbitals.py` to use the complete parameter set:

```bash
python3 example_orbitals.py \
  --sk-path ~/dftb_params/slako/mio-1-1/ \
  --sk-wfc-path ~/dftb_params/slako/wfc/wfc.mio-1-1.hsd
```

### Step 3: Verify WFC File Exists

```bash
ls ~/dftb_params/slako/wfc/wfc.mio-1-1.hsd
```

### Step 4: Run the Script

```bash
cd tests/dftb
python3 example_orbitals.py \
  --sk-path ~/dftb_params/slako/mio-1-1/ \
  --sk-wfc-path ~/dftb_params/slako/wfc/wfc.mio-1-1.hsd
```

## Alternative: Use Existing DFTB+ Installation

If you have DFTB+ source installed:

```bash
cd ~/git_SW/dftbplus
./utils/get_opt_externals ALL
```

Then use:
```bash
--sk-path ~/git_SW/dftbplus/external/slakos/mio-1-1/ \
--sk-wfc-path ~/git_SW/dftbplus/external/slakos/wfc/wfc.mio-1-1.hsd
```

## Quick Reference

| Parameter Set | SK Path | WFC File |
|--------------|---------|----------|
| mio-1-1 | `mio-1-1/` | `wfc.mio-1-1.hsd` |
| 3ob-3-1 | `3ob-3-1/` | `wfc.3ob-3-1.hsd` |

---

# Tutorial: Orbital Plotting Approaches

## Overview

Two approaches are available for plotting molecular orbitals from DFTB+:

1. **Cube file approach (recommended)** - Using waveplot to generate cube files
2. **C API approach (limited)** - Direct access to DFTB+ matrices via C API

## Why C API Approach Does Not Work for Orbital Plotting

The DFTB+ C API provides access to:
- Hamiltonian matrix (H)
- Overlap matrix (S)
- Density matrix (DM)
- System information

However, **it does NOT expose basis function evaluation routines**.

### The Problem

To evaluate molecular orbitals at arbitrary points, you need:

ψ(r) = Σᵢ Cᵢ · φᵢ(r)

Where:
- Cᵢ = Molecular orbital coefficients (available via eigenvectors from H and S)
- φᵢ(r) = Basis function values at point r (**NOT available in C API**)

Without basis function evaluation, you cannot compute orbital values at arbitrary points, even if you have the coefficients.

### What the C API Can Do

The C API approach (see `example_orbitals_c_api.py`) can:
- Extract H, S, and DM matrices during SCF calculation
- Solve generalized eigenvalue problem to get eigenvectors
- Compute eigenvalues
- Calculate electron count

But it cannot:
- Evaluate orbitals at arbitrary points
- Generate orbital density maps
- Create orbital visualizations

### Why Cube File Approach Works

The cube file approach (see `example_orbitals.py`) works because:
- waveplot has **internal access to DFTB+ basis function routines**
- waveplot evaluates orbitals on a grid using these routines
- Cube files contain pre-evaluated orbital values
- Interpolation from cube files is fast and accurate

### Comparison

| Feature | C API | Cube Files |
|---------|--------|------------|
| Get Hamiltonian/Overlap | ✓ Yes | ✗ No |
| Get eigenvectors | ✓ Yes (via diagonalization) | ✗ No |
| Evaluate orbitals at points | ✗ No (basis functions not exposed) | ✓ Yes (via waveplot) |
| Generate orbital plots | ✗ No | ✓ Yes |
| Performance | Fast for matrix extraction | Fast for plotting |
| Complexity | Requires C bindings | Simple (waveplot + interpolation) |

## Recommendation

**Use the cube file approach** (`example_orbitals.py`) for orbital plotting:
- It's the only working approach for visualization
- waveplot has proper access to internal DFTB+ routines
- Interpolation from cube files is efficient
- No need to modify DFTB+ source code

**Use the C API approach** (`example_orbitals_c_api.py`) only if you need:
- Hamiltonian/Overlap matrices for analysis
- Density matrix extraction
- Custom matrix operations

## Example Usage

### Cube File Approach (Recommended)

```bash
python3 example_orbitals.py \
  --sk-path ~/dftb_params/slako/mio-1-1/ \
  --sk-wfc-path ~/dftb_params/slako/wfc/wfc.mio-1-1.hsd \
  --plot-2d \
  --plane xy \
  --points "0,0,0"
```

This generates:
- Cube files with orbital/density data
- 2D projections with atom overlays
- Proper Angstrom coordinate scaling

### C API Approach (Limited)

```bash
python3 example_orbitals_c_api.py input.dftb --lib-path /path/to/libdftbplus.so
```

This extracts matrices but cannot plot orbitals directly.

## Future Work

To enable true C API orbital evaluation, you would need to:
1. Modify DFTB+ source code to expose basis evaluation functions
2. Add C API functions for basis function evaluation
3. Recompile DFTB+ with custom API
4. Write Python bindings for the new functions

This requires DFTB+ development expertise and is not recommended for general use.

---

# Orbital Plotting Scripts

## Available Scripts

### example_orbitals.py
**Location:** `/home/prokophapala/git/FireCore/tests/dftb/example_orbitals.py`

**Purpose:** Plot molecular orbitals and electron density from DFTB+ using waveplot cube files.

**Features:**
- 2D orbital projections (xy, xz, yz planes)
- Electron density plotting with atom overlays
- Proper Angstrom coordinate scaling (Bohr → Angstrom conversion)
- Small atom markers with symbols and indices (e.g., H0, H1, O2)
- CLI support for loading molecules (H2O, PTCDA, or custom XYZ files)
- Automatic DFTB+ calculation and waveplot execution

**Usage:**
```bash
# Basic usage with built-in H2O
python3 example_orbitals.py --plot-2d --points "0,0,0"

# Load specific molecule
python3 example_orbitals.py --molecule H2O --plot-2d --points "0,0,0"
python3 example_orbitals.py --molecule PTCDA --plot-2d --points "0,0,0"

# Custom XYZ file
python3 example_orbitals.py --molecule /path/to/molecule.xyz --plot-2d --points "0,0,0"

# Select projection plane
python3 example_orbitals.py --molecule H2O --plot-2d --plane xz --points "0,0,0"

# Skip DFTB+ calculation (use existing results)
python3 example_orbitals.py --skip-dftb --skip-waveplot --plot-2d
```

**Output:**
- Cube files: `orbital_calc/wp-*.cube`
- 2D plots: `orbital_plots/*_xy.png`, `*_xz.png`, `*_yz.png`
- Density plot: `orbital_plots/density_*.png`

**Requirements:**
- DFTB+ executable
- waveplot executable
- Slater-Koster files (e.g., 3ob-3-1)
- Wavefunction coefficient file (wfc.hsd)
- ASE (for XYZ loading and cube reading)

### example_orbitals_c_api.py
**Location:** `/home/prokophapala/git/FireCore/tests/dftb/example_orbitals_c_api.py`

**Purpose:** Demonstrate DFTB+ C API for matrix extraction (limited for orbital plotting).

**Features:**
- Extract Hamiltonian, Overlap, Density matrices
- Compute molecular orbital coefficients (eigenvectors)
- Calculate eigenvalues
- **Cannot** evaluate orbitals at arbitrary points (basis functions not exposed in C API)

**Usage:**
```bash
python3 example_orbitals_c_api.py input.dftb --lib-path /path/to/libdftbplus.so
```

**Limitations:**
- Cannot generate orbital visualizations
- Cannot evaluate orbitals at arbitrary points
- Only useful for matrix analysis

**Recommendation:** Use `example_orbitals.py` for orbital plotting.

## Key Implementation Details

### Bohr to Angstrom Conversion
Cube files store coordinates in Bohr. The script automatically converts:
```python
BOHR_TO_ANG = 0.529177210903
origin_ang = origin * BOHR_TO_ANG
spacing_ang = spacing * BOHR_TO_ANG
```

### Atom Position Handling
- Original input atoms (from XYZ file) are used for plotting
- ASE's `read_cube_data` returns atoms in Angstrom
- Grid extent from cube header is converted from Bohr to Angstrom
- This ensures atoms and density/orbital data are in the same coordinate system

### Plotting Style
- Small red `+` markers for atoms (s=20)
- Labels: symbol + index (e.g., H0, H1, O2)
- Black semi-transparent box around labels for readability
- Proper extent parameter in `imshow` for Angstrom coordinates

## Molecule Library
Pre-configured molecules available via `--molecule` argument:
- **H2O**: `cpp/common_resources/xyz/H2O.xyz`
- **PTCDA**: `cpp/common_resources/xyz/PTCDA.xyz`

Custom XYZ files can be specified with full path or relative path.

## Troubleshooting

### Coordinate System Mismatch
If atoms appear at wrong positions:
- Check that Bohr → Angstrom conversion is applied to origin/spacing
- Verify original input atoms are used (not cube file atoms)
- Ensure extent parameter uses converted Angstrom values

### Waveplot WFC File Errors
If waveplot fails with "Missing wavefunction coefficient files":
- Download complete parameter set from dftb.org
- Ensure WFC file matches SK set (e.g., wfc.mio-1-1.hsd with mio-1-1 SK files)
- Use `--sk-path` and `--sk-wfc-path` arguments

### SK File Not Found
If DFTB+ fails with "SK file does not exist":
- Check `--sk-path` argument points to valid SK directory
- Verify SK files exist for all elements in molecule
- For organic molecules, use 3ob-3-1 or mio-1-1 parameter sets

| ob2-1-1 | `ob2-1-1/` | `wfc.ob2-1-1.hsd` |

## Important Notes

- **WFC files must match SK set**: Using `wfc.mio-1-1.hsd` with `3ob-3-1` SK files is not recommended
- **Test params incomplete**: The test SK files in `git_SW/asi/tests/testcases/` lack WFC files and are only for regression testing
- **Physical consistency**: For meaningful orbital/density analysis, always use complete parameter sets

---

# Code Refactoring Summary (May 2026)

## Overview

Major refactoring effort to consolidate DFTB+ Python utilities into reusable modules (`pyBall/dftb_utils.py`, `pyBall/plotUtils.py`) and eliminate code duplication across example scripts.

## Changes to pyBall/dftb_utils.py

### New Functions Added

#### Atom Selection and Path Generation
- `select_atom_index(enames, apos, symbol, axis=1, mode='abs_min', value=0.0)` - Select atom by symbol and position criteria
- `find_closest_indices(enames, apos, target_idx, symbol, n=2)` - Find N closest atoms of given symbol to target
- `make_axis_path(p0, p1, svals)` - Generate linear path between two points with spacing values
- `identify_hbond_transfer(enames, apos, ...)` - Identify H-bond transfer atoms (donor, acceptor, H atoms)

#### DFTB+ I/O and Execution
- `load_molecule(filename, use_ase=True)` - Load molecule from XYZ/GEN using ASE or AtomicSystem
- `write_dftb_input_hessian(...)` - Write DFTB+ input for Hessian calculation with SecondDerivatives driver
- `write_dftb_input_orbitals(...)` - Write DFTB+ input for orbital calculation with eigenvector output
- `parse_energy_out(fname, allow_unconverged=False)` - Parse total energy from DFTB+ OUT file (with fallback for unconverged SCC)
- `parse_forces(fname, natoms)` - Parse forces from detailed.out
- `read_relaxed_geometry(apos, do_relax)` - Read relaxed geometry from geo_end.gen or return original
- `run_pbc(...)` - Unified function for periodic DFTB+ calculations (write input, run, parse energy/forces)
- `constrained_scan(...)` - General-purpose constrained geometry scan with fixed/moved atoms
- `save_xyz_movie(results, fname, lvs=None, ...)` - Save scan trajectory as XYZ movie

#### Hessian and Vibrational Analysis
- `read_hessian(filename, n_atoms=None)` - Read Hessian matrix from hessian.out
- `hessian_hartree_bohr_to_eV_angstrom(hessian)` - Convert Hessian units (Hartree/Bohr² → eV/Å²)
- `hessian_to_mass_weighted(hessian, masses)` - Convert Hessian to mass-weighted dynamical matrix
- `hessian_to_frequencies(hessian, masses)` - Compute vibrational frequencies from Hessian
- `write_vibration_modes_jmol(...)` - Write vibration modes in Jmol XYZ format

#### Waveplot and Cube Files
- `run_waveplot(...)` - Run waveplot to generate cube files from DFTB+ results
- `read_cube(filename)` - Read Gaussian cube file using ASE
- `read_cube_with_grid(filename)` - Read cube file with full grid information (origin, spacing)
- `evaluate_orbital_at_points(...)` - Evaluate molecular orbitals at arbitrary points using DFTB+ C API
- `interpolate_orbital(cube_file, points)` - Interpolate orbital values from cube file at points

### Default Parameters Centralized
```python
DEFAULT_SK_PATH = '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/'
DFTB_EXE = '/home/prokophapala/miniconda3/bin/dftb+'
default_params = { ... }  # Standard DFTB+ parameters
```

## Changes to pyBall/plotUtils.py

### New Functions Added

#### Cube File Visualization
- `plot_cube_slice(cube_file, atoms=None, plane='xy', ...)` - Plot 2D slice of cube data with atom overlay
  - Supports xy, xz, yz projections
  - Automatic atom labeling (symbol + index)
  - Colorbar and extent handling

#### Scan Profile Plotting
- `plot_scan_profile(L_vals, E_vals, fname, title, xlabel, ylabel)` - Plot energy vs reaction coordinate

#### Geometry with Forces
- `plotGeometryWithForces(apos, enames, forces=None, ...)` - Plot geometry with force arrows
  - Supports highlighting specific atoms
  - Force vector scaling and coloring

### Helper Functions
- `_cube_slice_data(data, origin, spacing, plane)` - Extract 2D slice and metadata from 3D cube data

## Refactored Scripts

### tests/pyFireball/scan_constrained.py
**Changes:**
- Removed local geometry builders (now uses `GrapheneRibbonBuilder.build_two_ribbon_cell`)
- Removed local DFTB run wrapper (now uses `dftb_utils.constrained_scan`)
- Removed local plotting code (now uses `plotUtils.plot_scan_profile`, `plotGeometryWithForces`)
- Removed local atom selection (now uses `dftb_utils.identify_hbond_transfer`)
- Removed local path generation (now uses `dftb_utils.make_axis_path`)
- Added CLI options for flexible testing (`--rigid`, `--n_steps`, `--L_min`, `--L_max`, `--outdir`)
- Script reduced from ~260 lines to ~140 lines

**Tested:** 3-step and 10-step rigid scans work correctly.

### tests/pyFireball/scan_LHb.py
**Changes:**
- Removed local geometry builders (now uses `GrapheneRibbonBuilder.build_two_ribbon_cell`)
- Removed local DFTB run wrapper (now uses `dftb_utils.run_pbc`)
- Removed local XYZ movie saving (now uses `dftb_utils.save_xyz_movie`)
- Removed local plotting (now uses `plotUtils.plot_scan_profile`)
- Added CLI options for scan parameters (`--width`, `--Lx`, `--LHb_min`, `--LHb_max`, `--n_steps`)
- Added convergence control options (`--nk`, `--Temperature`, `--MixingParameter`, `--MaxScc`, `--SCCTolerance`)
- Added `--allow_unconverged` flag for debugging SCC failures
- Fixed SCC convergence by updating defaults: `nk=16`, `Temperature=300K`, `MixingParameter=0.1`

**Tested:** 20-point scan (1.5-2.5 Å) with 18/20 points converged successfully.

### tests/dftb/example_hessian.py
**Changes:**
- Removed redundant `load_molecule_xyz()` (now uses `dftb_utils.load_molecule`)
- Removed redundant `read_hessian()` (now uses `dftb_utils.read_hessian`)
- Removed redundant `hessian_hartree_bohr_to_eV_angstrom()` (now uses `dftb_utils.hessian_hartree_bohr_to_eV_angstrom`)
- Already used shared functions for Hessian conversion and frequency calculation

**Tested:** H2O Hessian calculation works correctly with built-in and loaded molecules.

### tests/dftb/example_orbitals.py
**Changes:**
- Already uses `dftb_utils.load_molecule()` for molecule loading
- Already uses `dftb_utils.run_waveplot()` for cube generation
- Already uses `plotUtils.plot_cube_slice()` for 2D orbital projections
- Removed specialized `plot_orbitals_2d()` and `plot_density_at_points()` (kept simpler `plot_cube_slice`)

**Tested:** H2O orbital calculation with waveplot works correctly.

## Removed Files

- `tests/pyFireball/ribbon_utils.py` - Functions moved to `doc/Topics/Kekule_Topology/GrapheneRibbonBuilder.py`
- `tests/pyFireball/scan_constrained_legacy.py` - Legacy backup of original scan_constrained.py

## Key Benefits

1. **Reduced Code Duplication** - Common DFTB+ I/O, plotting, and geometry building now in shared modules
2. **Improved Maintainability** - Bug fixes in shared functions automatically benefit all scripts
3. **Consistent API** - All scripts use same function signatures and parameter defaults
4. **Better Testing** - Shared functions can be tested independently
5. **Lightweight Scripts** - Example scripts focus on workflow, not implementation details

## Testing Summary

All refactored scripts tested with small sweeps:
- `scan_constrained.py`: 3-step and 10-step rigid scans ✓
- `scan_LHb.py`: 20-point scan (18/20 converged) ✓
- `example_hessian.py`: H2O Hessian with/without frequency calculation ✓
- `example_orbitals.py`: H2O orbitals with waveplot ✓


