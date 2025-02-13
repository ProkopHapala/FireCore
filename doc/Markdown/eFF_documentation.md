# Electron Forcefield (EFF) Documentation

## 1. Files relevant for the  Electron forcefield (EFF) project

1.  The core eFF logic is implemented in C++ (`eFF_lib.cpp`, `eFF.h`, `InteractionsGauss.h`).
2.  The Python file (`eFF.py`) provides an interface to this C++ code, allowing it to be used from Python.
3.  The `run_tests.py` file uses the Python interface to test the C++ implementation.
4.  The markdown files (`EFF_otevrena_veda.md` and `EFF_otevrena_veda_gpt.md`) provide documentation and context for the code.

The `eFF.h` and `InteractionsGauss.h` files are included in `eFF_lib.cpp` to provide the definitions and implementations of the classes and functions used in the eFF calculation. The `eFF.py` file uses `ctypes` to load the compiled `eFF_lib.cpp` file and access its functions. The `run_tests.py` file imports the `eFF.py` file and uses its functions to test the eFF implementation.

### Description of the Files:

*   `eFF.py`: This file serves as a Python wrapper for the C++ eFF library (`eFF_lib.cpp`). It uses the `ctypes` module to load the compiled C++ library and provides Python functions that can be called to access the C++ functions. This allows users to interact with the eFF library from Python scripts.
*   `eFF_lib.cpp`: This file contains the core implementation of the eFF model in C++. It includes functions for calculating the energy and forces of the system, loading and saving data, and running simulations. This file relies on `eFF.h` and `InteractionsGauss.h` for the definition of the `EFF` class and the implementation of the interaction potentials.
*   `eFF.h`: This header file defines the `EFF` class, which is the central class in the eFF implementation. It includes the data structures for representing the system (atoms, electrons, positions, sizes, etc.) and the function declarations for the methods that operate on the system.
*   `InteractionsGauss.h`: This header file implements functions for evaluating the interaction between Gaussian functions, such as the product, overlap, kinetic energy, and electrostatic interaction. These functions are used in `eFF_lib.cpp` to calculate the energy and forces of the system.
*   `run_tests.py`: This file contains a set of unit tests for the eFF library. It uses the Python interface (`eFF.py`) to call the C++ functions and check their correctness. The tests cover various aspects of the eFF implementation, such as the energy and forces, the geometry relaxation, and the core size scanning.
*   `EFF_otevrena_veda.md` and `EFF_otevrena_veda_gpt.md`: These files contain background information, theory, and description of goals for the project.

## Functions implemented in each file:

### `tests/tEFF/run_tests.py`

This file contains the tests for the eFF library. It includes functions for checking the energy and forces, relaxing molecules, and scanning the core size of atoms.

*   `check_DerivsPauli(r0=0.5,r1=1.5, s0=0.5,s1=0.5,   sj=1.0, n=10, spin=1 )`: Checks the derivatives of the Pauli potential.
*   `checkNumDerivs(name, d=0.001, bDetail=False, bPrint=True)`: Checks the numerical derivatives of the energy with respect to the degrees of freedom.
*   `check_Derivs_ie(name, ie=0, r0=0.5,r1=1.5, s0=0.5,s1=0.5, n=10 )`: Checks the derivatives of the energy with respect to the electron positions and sizes.
*   `relax_mol(name, dt=0.03,damping=0.1, bTrj=True, bResult=True, perN=1, bPrintLbonds=True, nMaxIter=10000, outE=None, outF=None, fUnits=1., bFixNuclei=False )`: Relaxes a molecule to its equilibrium geometry.
*   `scan_core_size(name, core_sizes, dt=0.03,damping=0.1, nMaxIter=10000, fUnits=1., ia=0 )`: Scans the core size of an atom and calculates the bond lengths.
*   `check_H2(bRelax=True, name="H2_eFF", bPyeff=True)`: Checks the H2 molecule properties.

### `pyBall/eFF.py`

This file provides a Python interface to the C++ eFF library. It uses `ctypes` to load the library and define the function signatures. It loads C++ libaray `eFF_lib.cpp` compiled into `libeFF.so`.

*   `setVerbosity(verbosity_=0, idebug=0)`: Sets the verbosity level for debugging.
*   `init_buffers()`: Initializes the buffers used in the C++ library.
*   `load_xyz(fname)`: Loads atomic coordinates from an XYZ file.
*   `load_fgo(fname, bVel_=False, fUnits=1.)`: Loads atomic coordinates, electron positions, and sizes from a custom FGO file format.
*   `save_fgo(filename, bVel=False, bAppend=False)`: Saves atomic coordinates, electron positions, and sizes to a custom FGO file format.
*   `save_xyz(filename, bAppend=False)`: Saves atomic coordinates to an XYZ file.
*   `setTrjName(trj_fname_="trj.xyz", savePerNsteps=1, bDel=True )`: Sets the trajectory file name for saving the simulation progress.
*   `init(na, ne, bVel_=False)`: Initializes the eFF system with `na` atoms and `ne` electrons.
*   `eval()`: Evaluates the energy and forces of the system.
*   `evalFuncDerivs(r, s, Es=None, Fs=None, Fr=None, ie=0)`: Evaluates the energy and forces for a given electron, used for testing.
*   `info()`: Prints information about the system.
*   `getEnergyTerms(sh=(7,))`: Returns the different energy terms of the system (kinetic, electron-electron, electron-atom, etc.).
*   `getDimPointer(sh=(3,))`: Returns the dimensions of the system (number of atoms, electrons, and degrees of freedom).
*   `getBuff(name, sh)`: Returns a pointer to a buffer in the C++ library.
*   `getIBuff(name,sh)`: Returns a pointer to an integer buffer in the C++ library.
*   `setPauliModel(i)`: Sets the Pauli repulsion model.
*   `setKPauli(KPauli)`: Sets the scaling factor for the Pauli potential.
*   `setSwitches(kinetic=0, coulomb=0, pauli=0, AA=0, AE=0, AECoulomb=0, AEPauli=0)`: Enables or disables different energy terms in the calculation.
*   `initOpt(dt=0.1, damping=0.1, f_limit=1000.0, bMass=False )`: Initializes the optimizer for geometry relaxation.
*   `run(nstepMax=1000, dt=None, Fconv=1e-6, ialg=0, outE=None, outF=None)`: Runs the simulation for a given number of steps or until convergence.
*   `evalNumDerivs(Fnum=None, d=0.01)`: Evaluates numerical derivatives of the energy with respect to the degrees of freedom.
*   `sample_ee(RSs, spin, FEout=None, KRSrho=[1.125,0.9,-0.2], bEvalCoulomb=True, bEvalPauli=True, iPauliModel=1 )`: Samples the electron-electron interaction energy for given distances and spins.
*   `sample_EA(RSs, FEout=None, KRSrho=[1.125,0.9,-0.2], aPar=[4.,0.1,0.1,2.0], bEvalAECoulomb=True, bCoreCoul=True, bEvalAEPauli=True)`: Samples the electron-atom interaction energy for given distances and spins.
*   `relax_mol(name, dt=0.03,damping=0.1, bTrj=True, bResult=True, perN=1, bPrintLbonds=True, nMaxIter=10000, outE=None, outF=None, fUnits=1., bFixNuclei=False )`: Relaxes a molecule to its equilibrium geometry.
*   `scan_core_size(name, core_sizes, dt=0.03,damping=0.1, nMaxIter=10000, fUnits=1., ia=0 )`: Scans the core size of an atom and calculates the bond lengths.
*   `check_H2(bRelax=True, name="H2_eFF", bPyeff=True)`: Checks the H2 molecule properties.
*   `init_eff(natom_=0, nelec_=1, s=0.5,  aQ=1.0,aQs=0.0,aP=0.0,aPs=0.1 )`: Initializes the eFF system with given parameters.
*   `scan_size(ss, ie0)`: Scans the size of an electron and calculates the energy.
*   `test_Hatom(bDerivs=False)`: Tests the Hydrogen atom properties.

### `cpp/libs/Molecular/eFF_lib.cpp`

This file is basically just C++ side of the interface exporting eFF library functions as `extern "C" {..}` so they can be called from Python or other languages ( C++ avoiding name mangling). It includes functions for evaluating the energy and forces, loading and saving data, and running simulations. This compiles into `libeFF_lib.so` which is then loaded by `pyBall/eFF.py`

*   `setVerbosity( int verbosity_, int idebug_ )`: Sets the verbosity level for debugging.
*   `eval()`: Evaluates the energy of the system.
*   `evalFuncDerivs( int ie, int n, double* r, double* s, double* Es, double* Fr, double* Fs )`: Evaluates the energy and forces for a given electron, used for testing.
*   `evalNumDerivs( double* Fnum, double d )`: Evaluates numerical derivatives of the energy with respect to the degrees of freedom.
*   `init_buffers()`: Initializes the buffers used in the C++ library.
*   `setTrjName( char* trj_fname_, int savePerNsteps_ )`: Sets the trajectory file name for saving the simulation progress.
*   `load_xyz( const char* fname )`: Loads atomic coordinates from an XYZ file.
*   `load_fgo( const char* fname, bool bVel, double fUnits )`: Loads atomic coordinates, electron positions, and sizes from a custom FGO file format.
*   `init( int na, int ne )`: Initializes the eFF system with `na` atoms and `ne` electrons.
*   `info()`: Prints information about the system.
*   `getEnergyPointer()`: Returns a pointer to the energy of the system.
*   `getDimPointer()`: Returns a pointer to the dimensions of the system (number of atoms, electrons, and degrees of freedom).
*   `initOpt( double dt, double damping, double f_limit, bool bMass )`: Initializes the optimizer for geometry relaxation.
*   `run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF )`: Runs the simulation for a given number of steps or until convergence.
*   `sample_ee( int n, double* RSs_, double* FEout_, int spin, double* KRSrho_, bool bEvalCoulomb, bool bEvalPauli, int iPauliModel )`: Samples the electron-electron interaction energy for given distances and spins.
*   `sample_EA( int n, double* RSs_, double* FEout_, double* KRSrho_,  double* aPar_,  bool bEvalAECoulomb, bool bCoreCoul, bool bEvalAEPauli )`: Samples the electron-atom interaction energy for given distances.
*   `save_fgo( char const* filename, bool bVel, bool bAppend )`: Saves atomic coordinates, electron positions, and sizes to a custom FGO file format.
*   `save_xyz( char const* filename, bool bAppend )`: Saves atomic coordinates to an XYZ file.
*   `setPauliModel(int i  )`: Sets the Pauli repulsion model.
*   `setKPauli( double KPauli )`: Sets the scaling factor for the Pauli potential.
*   `setSwitches( int bEvalKinetic, int bEvalCoulomb, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli )`: Enables or disables different energy terms in the calculation.

### `cpp/common/molecular/eFF.h`

This file defines the `EFF` class, which is the core of the eFF implementation. It includes the data structures and functions for representing the system and evaluating its energy and forces.

*   `EFF()`: Constructor for the `EFF` class.
*   `~EFF()`: Destructor for the `EFF` class.
*   `realloc(int na_, int ne_, bool bVel=false)`: Reallocates the memory for the system with `na_` atoms and `ne_` electrons.
*   `makeMasses(double*& invMasses, double m_const=-1)`: Makes the masses for the atoms and electrons.
*   `dealloc()`: Deallocates the memory for the system.
*   `fixElectron(int ie, double* vs=0)`: Fixes the position of an electron.
*   `evalKinetic()`: Evaluates the kinetic energy of all each electron in the system: $ E_k = \frac{3}{2} \frac{\hbar^2}{m_e s_i^2} $
    where $s_i$ is the size of the electron.
    Calls: `addKineticGauss_eFF()` from `InteractionsGauss.h`.
*   `evalEE()`: Evaluates the electron-electron interaction energy. The electron-electron interaction energy is calculated using the following equation:
   * $E_{ee} = E_{Coulomb} + E_{Pauli}$
      * $E_{Coulomb} = \sum_{i<j} \frac{q_i q_j}{r_{ij}} \text{erf}(\frac{r_{ij}}{s_{ij}})$ where $q_i$ and $q_j$ are the charges of the electrons, $r_{ij}$ is the distance between the electrons, $s_{ij}$ is the effective  size of the electrons, and $E_{Pauli}$ is the Pauli repulsion energy. $\text{erf}(x)$ is error function, that is integral of Gaussian function.
      * $E_{Pauli}$ is more complicated and our goal is to replace $E_{Pauli}$ with machine-learned function which better reproduce molecular geometries and energies from more accurate quentum methods.
    Calls: `addCoulombGauss()`, `addPauliGauss_New()`, `addPauliGaussVB()`, `addDensOverlapGauss_S()` from `InteractionsGauss.h`.
*   `evalAE()`: Evaluates the atom-electron interaction energy.
      * $E_{ae} = E_{Coulomb} + E_{CER}$
      * $E_{Coulomb}$ is again computed by $\frac{Q_a q_i}{r_{aj}} \text{erf}(\frac{r_{ai}}{s_i})$, where $Q_a$ is the charge of the atom, $q_i$ is the charge of the electron, $r_{ai}$ is the distance between the atom and the electron, $s_i$ is the effective size of the electron.
      * $E_{CER}$ is the core-electron repulsion energy. This energy can be calculated as $E_{ee}$ defined above. But we want to replace it by new machine-learned function which allows us to tune each atom type to reproduce molecular geometries and energies from more accurate quantum methods.
    Calls: `addCoulombGauss()`, `addPauliGauss_New()` from `InteractionsGauss.h`.
*   `evalAA()`: Evaluates the atom-atom interaction energy. This is simple Coulomb interaction between point-charges $E_{Coulomb}=\frac{Q_a Q_b}{r_{ab}}$.
*   `evalCoreCorrection()`: Evaluates the core correction energy. Calls:  `addKineticGauss_eFF()`, `addCoulombGauss()` from `InteractionsGauss.h`.
*   `eval()`: Evaluates the total energy of the system.
    Calls: `evalKinetic()`, `evalEE()`, `evalAE()`, `evalAA()`, `evalCoreCorrection()`.
*   `clearForce()`: Clears the forces on the atoms and electrons.
*   `clearForce_noAlias()`: Clears the forces on the atoms and electrons (no aliasing).
*   `move_GD(double dt)`: Moves the atoms and electrons using the gradient descent algorithm.
*   `move_GD_noAlias(double dt)`: Moves the atoms and electrons using the gradient descent algorithm (no aliasing).
*   `electronPotAtPoint( const Vec3d& pi, double si, double Q, int spini=0, bool bEvalCoulomb=True )const`: Calculates the potential at a given point due to the electrons. Calls: `addCoulombGauss()`, `addPauliGauss_New()`, `addPauliGaussVB()`, `addDensOverlapGauss_S()` from `InteractionsGauss.h`.
*   `atomsPotAtPoint( const Vec3d& pos, double s, double Q )const`: Calculates the potential at a given point due to the atoms. Calls: `addCoulombGauss()`, `addDensOverlapGauss_S()` from `InteractionsGauss.h`.
*   `evalPotAtPoints( int n, Vec3d* ps, double* out=0, double s=0.0, double Q=1.0, int spin=0, bool bAtom=True, bool bElectron=False )const`: Evaluates the potential at a given set of points.
*   `printEnergies()`: Prints the different energy terms of the system.
*   `printAtoms()`: Prints the atomic positions and charges.
*   `printElectrons()`: Prints the electron positions and sizes.
*   `info()`: Prints information about the system.
*   `Eterms2str(char* str)`: Converts the energy terms to a string.
*   `orb2str(char* str, int ie)`: Converts the orbital information to a string.
*   `orbs2str(char* str0)`: Converts the orbital information to a string.
*   `to_xyz( FILE* pFile )`: Writes the atomic and electron positions to an XYZ file.
*   `save_xyz( const char* filename, const char const* mode="w" )`: Saves the atomic and electron positions to an XYZ file.
*   `loadFromFile_xyz( const char* filename )`: Loads the atomic and electron positions from an XYZ file.
*   `loadFromFile_fgo( char const* filename, bool bVel=false, double fUnits=1. )`: Loads the atomic and electron positions from a custom FGO file.
*   `writeTo_fgo( char const* filename, bool bVel=false, const char* fmode="w" )`: Writes the atomic and electron positions to a custom FGO file.

### `cpp/common/molecular/InteractionsGauss.h`

This file implements functions for evaluating the interaction between Gaussian functions, such as the product, overlap, kinetic energy, and electrostatic interaction.

*   `addKineticGauss( double s, double& fs )`: (NOT USED)
    Calculates the kinetic energy of a Gaussian function.
*   `addKineticGauss_eFF( double s, double& fs )`:
    Calculates the kinetic energy of a Gaussian function, used in eFF.
    Called from: `EFF::evalKinetic()`, `EFF::evalCoreCorrection()`
*   `CoulombGauss( double r, double s, double& fr, double& fs, double qq )`:
    Calculates the Coulomb interaction between two Gaussian functions.
*   `CoulombGauss_FixSize( double r, double beta, double& fr, double qq )`: (NOT USED)
    Calculates the Coulomb interaction between two Gaussian functions with fixed sizes.
*   `addCoulombGauss( const Vec3d& dR, double s, Vec3d& f, double& fsi, double qq )`:
    Adds the Coulomb interaction between two Gaussian functions to the force.
    Called from: `EFF::evalAE()`, `EFF::electronPotAtPoint()`.
*   `addCoulombGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, double qq )`:
    Adds the Coulomb interaction between two Gaussian functions to the force.
    Called from: `EFF::evalEE()`, `EFF::evalCoreCorrection()`, `EFF::atomsPotAtPoint()`.
*   `DensOverlapGauss_S( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj, double si2, double sj2, double is2, double is4 )`:
    Calculates the density overlap between two Gaussian functions.
*   `DensOverlapGauss_Snorm( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj, double si2, double sj2, double is2, double is4 )`: (NOT USED)
    Calculates the normalized density overlap between two Gaussian functions.
*   `DensOverlapGauss_P( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj, double si2, double sj2, double is2, double is4 )`: (NOT USED)
    Calculates the density overlap between two Gaussian functions.
*   `addDensOverlapGauss_S( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj )`:
    Adds the density overlap between two Gaussian functions to the force.
    Called from: `EFF::evalEE()`, `EFF::atomsPotAtPoint()`, `EFF::electronPotAtPoint()`.
*   `addDensOverlapGauss_P( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj )`: (NOT USED)
    Adds the density overlap between two Gaussian functions to the force.
*   `getDeltaTGauss( double r2, double si, double sj,  double& dTr, double& dTsi, double& dTsj, double isi2, double isj2, double s2, double is2, double is4 )`:
    Calculates the kinetic energy difference between two Gaussian functions.
*   `getOverlapSGauss( double r2, double si, double sj, double& dSr, double& dSsi, double& dSsj, double si2, double sj2, double is2, double is4 )`:
    Calculates the overlap between two Gaussian functions.
*   `PauliSGauss_anti( double S, double& fS, double rho )`: Calculates the Pauli repulsion energy for anti-parallel spins.
*   `PauliSGauss_syn( double S, double& fS, double rho )`: Calculates the Pauli repulsion energy for parallel spins.
*   `addPauliGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho )`: (NOT USED)
    Adds the Pauli repulsion energy to the force.
*   `addPauliGauss_New( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, int spin, const Vec3d& KRSrho, double sc=1.0 )`:
    Adds the Pauli repulsion energy to the force (new version).
    Called from: `EFF::evalEE()`, `EFF::evalAE()`, `EFF::electronPotAtPoint()`.
*   `addPauliGaussVB( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj )`:
    Adds the Pauli repulsion energy to the force (valence bond version).
    Called from: `EFF::evalEE()`, `EFF::electronPotAtPoint()`.
*   `PauliCoreElec_Orig(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_A, double PAULI_CORE_B, double PAULI_CORE_C)`: (NOT USED)
    Calculates the Pauli repulsion energy between the core and valence electrons.
*   `PauliCoreElec(double r, double re2, double& epauli, double& frc, double& fre2, double A, double B, double C)`: (NOT USED)
    Calculates the Pauli repulsion energy between the core and valence electrons.
*   `PauliCorePElec_Orig(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_P_A, double PAULI_CORE_P_B, double PAULI_CORE_P_C, double PAULI_CORE_P_D, double PAULI_CORE_P_E)`: (NOT USED)
    Calculates the Pauli repulsion energy between the core and valence electrons.
*   `PauliCorePElec(double rc, double re2, double& epauli, double& frc, double& fre2, double A, double B, double C, double D, double E_ )`: (NOT USED)
    Calculates the Pauli repulsion energy between the core and valence electrons.


