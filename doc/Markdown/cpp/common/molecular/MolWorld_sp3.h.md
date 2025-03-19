# MolWorld_sp3.h

## File Purpose

`MolWorld_sp3.h` is a header file that defines the `MolWorld_sp3` class, which serves as a comprehensive container for managing molecular simulations. It handles various aspects such as building molecules, applying force fields (MMFF and UFF), managing constraints, and performing optimizations. The class inherits from `SolverInterface`, providing methods to initialize, run, and analyze molecular dynamics simulations.

## Includes

- `<stdlib.h>`: Standard library header for basic utilities.
- `<stdio.h>`: Standard I/O library for input/output operations.
- `<sys/stat.h>`: System-specific file status functions.
- `<string.h>`: String manipulation functions.
- `<vector>`: C++ standard library container for dynamic arrays.
- `<math.h>`: Mathematical functions and constants.
- `<omp.h>`: OpenMP API for parallel programming in C/C++.
- `IO_utils.h`: Utility functions for input/output operations.
- `fastmath.h`: Optimized mathematical functions.
- `Vec3.h`: 3D vector class for geometric calculations.
- `Mat3.h`: 3x3 matrix class for linear algebra operations.
- `Vec3Utils.h`: Utility functions for vectors.
- `MMFFparams.h`: Parameters and constants for MMFF force field.
- `constants.h`: Global constants used throughout the code.
- `Forces.h`: Force calculation utilities.
- `MMFFsp3.h`: MMFF sp3 force field implementation.
- `MMFFsp3_loc.h`: Localized version of MMFF sp3 force field.
- `MMFFf4.h`: MMFF f4 force field implementation.
- `UFF.h`: UFF (Universal Force Field) implementation.
- `NBFF.h`: Non-bonded force field interface.
- `GridFF.h`: Grid-based non-bonded force field.
- `RigidBodyFF.h`: Rigid body dynamics utilities.
- `QEq.h`: Quantum electrostatics calculations.
- `constrains.h`: Constraint handling for molecular structures.
- `molecular_utils.h`: Utility functions for molecular operations.
- `LimitedGraph.h`: Graph-based data structure for limited graph traversal.
- `Molecule.h`: Molecular structure representation.
- `MMFFBuilder.h`: Builder class for constructing molecules.
- `SMILESparser.h`: Parser for SMILES string representations of molecules.
- `DynamicOpt.h`: Dynamic optimization utilities.
- `MultiSolverInterface.h`: Interface for multiple solver implementations.
- `GlobalOptimizer.h`: Global optimization algorithms.
- `GOpt.h`: Global optimizer implementation.
- `Groups.h`: Group-based molecular structure representation.
- `datatypes_utils.h`: Utility functions for data type manipulations.
- `arrayAlgs.h`: Array algorithm utilities.
- `SVG_render.h`: SVG rendering utilities.
- `EwaldGrid.h`: Ewald summation grid for periodic boundary conditions.

---

## Types (Classes and Structs)

---

### Class `MolWorld_sp3`

#### Purpose

`MolWorld_sp3` is a comprehensive class that manages the state of molecular simulations. It handles various aspects such as building molecules, applying force fields, managing constraints, and performing optimizations. This class inherits from `SolverInterface`, providing methods for initializing, running, and analyzing molecular simulations.

#### Inheritance
- **Parent Class:** `SolverInterface`

#### Properties

##### General Settings
- `bool isInitialized`: Indicates whether the simulation has been initialized.
- `const char* data_dir`: Directory containing common resources like element types and atom types.
- `const char* xyz_name`: Name of the input XYZ file for loading molecular structures.
- `const char* surf_name`: Name of the surface file if applicable.
- `const char* smile_name`: SMILES string used to initialize molecules.
- `const char* constr_name`: Name of the constraints file.
- `int savePerNsteps`: Frequency of saving trajectory data.
- `const char* trj_fname`: File name for trajectory output.
- `int iterPerFrame`: Number of iterations per frame in visualization.
- `int iParalel`: Parallelization level (0: no parallel, 1: OpenMP).
- `int iParalelMax`: Maximum parallelization level supported.
- `int iParalelMin`: Minimum parallelization level supported.
- `int iParalel_default`: Default parallelization level.
- `bool bOcl`: Flag for using OpenCL (not used in this version).
- `bool bTricubic`: Flag for tricubic interpolation (not used).
- `bool bPlaneSurfForce`: Flag for plane surface forces (not used).
- `bool b141`: Flag for 1-4 interactions (not used).
- `bool bSimple`: Simplified rule for assigning UFF types.
- `bool bConj`: Flag for conjugation in sp3 nitrogen and oxygen atoms.
- `bool bCumulene`: Exception to avoid cumulenes.
- `bool bRigid`: Flag for rigid body dynamics (not used).
- `bool bAnimManipulation`: Flag for animation and manipulation features (not used).
- `bool bWhichAtomNotConv`: Flag to identify atoms not converging in optimization.
- `bool bCheckInit`: Flag to check initialization of the system.
- `bool bBondInitialized`: Flag indicating if bonds have been initialized.
- `float anim_speed`: Speed factor for animation (not used).
- `int itest`: Test index for debugging purposes (not used).
- `int isubs`: Index of substitution molecule, if any.
- `bool bRelax`: Flag to enable relaxation during optimization.

##### Molecular Data
- `Vec3d* apos_bak`: Backup of atomic positions before optimization.
- `Vec3i nMulPBC`: Dimensions for multiplying periodic boundary conditions (PBC).
- `Vec3d cog`: Center of geometry vector.
- `Groups groups`: Group-based molecular structure representation.
- `std::vector<int> atom2group`: Mapping of atoms to groups.
- `Mat3d bbox`: Bounding box dimensions.
- `std::vector<int> selection`: Selection of atoms for manipulation.
- `std::unordered_set<int> selection_set`: Set version of the selection list.
- `std::vector<int> constrain_list`: List of constrained bonds or angles.
- `std::vector<Vec3i> Hbonds`: List of hydrogen bonds found in the system.
- `Vec3d pivotPoint`: Pivot point for rotation and translation operations.
- `Vec3d anim_vec`: Vector for animation manipulation (not used).
- `Vec3d manipulation_p0`: Pivot point for manipulation operations.
- `Vec3d manipulation_ax`: Axis vector for rotation and translation.
- `Vec3d pick_hray`: Direction vector for picking interactions.
- `int* manipulation_sel`: Pointer to selected indices for manipulation.
- `int manipulation_nsel`: Number of selected atoms for manipulation.
- `Vec3d* picked_lvec`: Lattice vector for picking operations (not used).
- `Mat3d* dlvec`: Displacement vector for lattice changes (not used).
- `Mat3d* latscan_dlvec`: Lattice displacement vectors for scanning.
- `Vec2d bySurf_c0`: Coordinates for changing the cell.
- `Vec2d bySurf_lat[2]`: Lattice vectors for changing the cell.
- `Mat3d new_lvec`: New lattice vector after changes.
- `Mat3d debug_rot`: Debugging matrix for molecular orientation.
- `int bySurf_ia0`: Index of atom used for changing the cell.
- `bool bConstrZ`: Flag to apply constraints along z-axis.
- `double ConstrZ_xmin`: Minimum x-coordinate for constraining atoms.
- `double ConstrZ_l`: Length of the constraint region.
- `double ConstrZ_k`: Spring constant for z-axis constraints.
- `double Kfix`: Fixed constraint strength.

##### Force Fields
- `MMFFparams* params`: Pointer to MMFF parameters.
- `MMFFsp3 ff`: MMFF sp3 force field implementation.
- `MMFFsp3_loc ffl`: Localized version of the MMFF sp3 force field.
- `UFF ffu`: UFF (Universal Force Field) implementation.
- `NBFF surf`: Non-bonded force field for surface interactions.
- `NBFF nbmol`: Non-bonded force field for non-surface molecules.
- `GridFF gridFF`: Grid-based non-bonded force field.
- `EwaldGrid gewald`: Ewald summation grid for periodic boundary conditions.
- `bool bMMFF`: Flag to use MMFF force field.
- `bool bUFF`: Flag to use UFF force field.
- `bool bGridFF`: Flag for grid-based force field.
- `bool bNonBonded`: Flag to include non-bonded interactions.
- `bool bConstrains`: Flag to apply constraints during simulation.
- `double gridStep`: Grid step size for non-bonded interactions.
- `bool doBonded`: Flag to include bonded interactions.
- `bool bSurfAtoms`: Flag for surface atoms.
- `bool bPBC`: Flag for periodic boundary conditions.
- `int npbc`: Number of periodic boundary conditions.

##### Optimization & Simulation
- `DynamicOpt opt`: Dynamic optimization utilities.
- `DynamicOpt optRB`: Optimizer for rigid bodies.
- `GlobalOptimizer gopt`: Global optimizer implementation.
- `OptLog opt_log`: Logging object for optimization details.
- `bool bOptimizer`: Flag to enable optimization during simulation.
- `bool bGopt`: Flag to enable global optimization.
- `bool bCheckInvariants`: Flag to check invariants of the system.
- `bool bRelaxPi`: Flag to relax pi orbitals.
- `bool bChargeUpdated`: Flag indicating if charges have been updated.
- `bool bNonBondNeighs`: Flag for non-bonded neighbors.
- `bool bCheckStuck`: Flag to check for stuck atoms during simulation.
- `double Etot`: Total energy of the system.
- `double maxVcog`: Maximum velocity of center of geometry.
- `double maxFcog`: Maximum force on center of geometry.
- `double maxTg`: Maximum torque on center of geometry.
- `double Kmorse`: Morse potential parameter.
- `double Ftol_default`: Default force tolerance for convergence.
- `double dt_default`: Default time step size.
- `double time_per_iter`: Average time per iteration during simulation.
- `bool bEpairs`: Flag to include electron pairs in calculations.
- `double fAutoCharges`: Charge assignment parameter.
- `bool bCheckStuck`: Flag to check for stuck atoms during simulation.
- `double RStuck`: Radius threshold for detecting stuck atoms.
- `int nStuckMax`: Maximum number of iterations before considering the system as stuck.
- `int nStuckTrj`: Number of consecutive iterations without change in trajectory.
- `int nStuck`: Current count of stuck iterations.
- `GridShape MOgrid`: Grid shape for molecular orbitals.

#### Methods

##### Initialization & Setup
- `void init()`: Initializes the molecular world with specified parameters and settings.
- `void pre_loop()`: Prepares the system for the main simulation loop by setting up groups and other necessary configurations.
- `void initParams()`: Initializes global parameters required by the molecular simulation.
- `void insertSMILES()`: Inserts a SMILES string into the molecular world for further processing or simulation.
- `void buildMolecule_xyz()`: Builds a molecule from an XYZ file, assigning types and initializing force fields.
- `void makeMoleculeTopology()`: Builds the topology of the molecule, including bonds and angles, based on loaded atomic positions.
- `void assingMoleculeTopoTypes()`: Assigns types to atoms in the molecule using predefined rules and parameters.
- `void loadGeom()`: Loads a molecular structure from an XYZ file into the system with optional initialization of periodic boundary conditions.
- `void makeMMFFs()`: Generates MMFF force field representations for bonded and non-bonded interactions in the system.
- `void makeFFs()`: Sets up various force fields (MMFF, UFF) and applies necessary configurations to the molecular world.
- `void setOptimizer()`: Sets up the optimizer for dynamic optimization during simulations.
- `void initRigid()`: Initializes rigid body dynamics for molecules in the system.
- `void initWithSMILES()`: Initializes the molecular world using a SMILES string, setting up parameters and structures accordingly.
- `void clear()`: Clears the entire molecular world, including parameters and structures, optionally clearing surface-related data as well.

##### Force Field Management
- `void setNonBond()`: Sets up non-bonded interactions for the molecular system based on specified parameters.
- `void updateBuilderFromFF()`: Updates the builder object with positions and charges from the force fields.
- `void clearFFs()`: Clears all force field data structures in preparation for new simulations.

##### Energy Evaluation & Simulation Control
- `void eval()`: Evaluates the total energy of the molecular system using various force fields and interactions.
- `void run_omp_Milan()`: Runs a molecular dynamics simulation with OpenMP parallelization, optimizing for convergence based on specified criteria.
- `void relax()`: Performs relaxation of the molecular system to minimize potential energy.
- `void run()`: Executes the main molecular dynamics loop, optionally using different algorithms and parameters.
- `void pullAtom()`: Applies a spring force to pull an atom towards a target position.
- `void MDloop()`: Runs a molecular dynamics simulation with specified time steps and convergence criteria.
- `void eval_no_omp()`: Evaluates the total energy of the system without using OpenMP parallelization, useful for debugging or specific scenarios.
- `void run_no_omp()`: Executes the main molecular dynamics loop without using OpenMP parallelization, providing a non-parallelized version of the simulation.
- `void run_omp()`: Runs a molecular dynamics simulation with OpenMP parallelization, optimizing performance and convergence based on specified criteria.
- `void update_GOpt()`: Updates global optimization settings for ongoing simulations.
- `void checkStuck()`: Checks if any atoms have become stuck, indicating potential issues in the simulation.
- `void handleStuckAtom()`: Handles atoms that are detected as being stuck by logging their details and potentially restarting the simulation.
- `void checkInvariants()`: Checks if any invariants of the molecular world exceed predefined thresholds.

##### Molecular Manipulation & Analysis
- `void shift_atoms()`: Shifts selected atoms by a specified displacement vector.
- `void rotate_atoms()`: Rotates selected atoms around a specified axis and pivot point.
- `void splitAtBond()`: Splits the selection of atoms at a specified bond index.
- `void selectByType()`: Selects atoms based on their atomic types or element symbols.
- `void selectRect()`: Selects atoms within a rectangular region defined by two points in 3D space and an optional rotation matrix.
- `void selectFragment()`: Selects all atoms belonging to a specified fragment of the molecular structure.
- `void selectAllBonded()`: Selects all atoms bonded to a given atom, including those connected through multiple bonds or rings.
- `void selectionFromBuilder()`: Updates the selection list based on the current builder object's selected atoms.
- `void trySel()`: Tries to retrieve and return the current selection of atoms from the molecular world.
- `void center()`: Calculates the center of geometry for a given set of selected atoms, optionally translating them to this point.
- `void getInertiaTensor()`: Computes the inertia tensor for a specified subset of atoms in the system.
- `void deleteAtomSelection()`: Deletes all atoms from the current selection list.
- `void clearSelections()`: Clears both the global and local selection lists, resetting them to empty states.
- `void selectAll()`: Selects all atoms in the molecular world for manipulation or analysis.
- `void selectInverse()`: Inverts the current selection of atoms by selecting those not already included in the list.
- `void fragmentsByBonds()`: Groups atoms into fragments based on their connectivity, identifying separate molecular entities within the system.
- `whichAtomNotConv()`: Identifies atoms that are not converging in the current simulation or optimization process.
- `getMostDisplacedAtom()`: Finds the atom with the largest displacement during a simulation run.

##### Scanning & Trajectory Generation
- `void scanTranslation_ax()`: Performs a translation scan along a specified axis for a given number of steps, optionally saving trajectory data to an XYZ file.
- `void scanTranslation()`: Performs a translation scan between two selected atoms along the line connecting them, optionally saving trajectory data to an XYZ file.
- `void scanRotation_ax()`: Performs a rotational scan around a specified axis and pivot point for a given number of steps, optionally saving trajectory data to an XYZ file.
- `void scanRotation()`: Performs a rotational scan around the axis defined by two selected atoms for a given number of steps, optionally saving trajectory data to an XYZ file.
- `void scanAngleToAxis_ax()`: Scans angles between selected atoms and a specified axis, adjusting their positions accordingly while evaluating energy changes at each step.
- `void toXYZ()`: Converts the current molecular structure into an XYZ file format for visualization or further analysis.

##### Charge & Bond Handling
- `void autoCharges()`: Automatically calculates charges for the molecular system using predefined rules and parameters.
- `void findHbonds_PBC()`: Finds hydrogen bonds within the molecular structure, considering periodic boundary conditions.

##### System Configuration & Analysis
- `void evalPBCshifts()`: Evaluates periodic boundary condition shifts for a given set of PBC dimensions and lattice vectors.
- `void makePBCshifts()`: Creates an array of periodic boundary condition shifts based on provided dimensions and lattice vectors.
- `void printPBCshifts()`: Prints the calculated periodic boundary condition shifts to the console.
- `void change_lvec()`: Changes the lattice vector of the molecular world.
- `void add_to_lvec()`: Adds a displacement vector to the current lattice vector.
- `void change_lvec_relax()`: Relaxes the lattice vector by applying small changes iteratively.
- `void optimizeLattice_1d()`: Optimizes the lattice in one dimension using a scan approach.
- `void changeCellBySurf()`: Changes the cell dimensions based on surface parameters and optional atom shifts.
- `void findBridgeBonds()`: Finds bridge bonds within the molecular structure using limited graph traversal.
- `void PBC_multiply()`: Multiplies the periodic boundary conditions of a fragment by a given factor in each direction.
- `void setConstrains()`: Sets up and applies constraints for the molecular structure.
- `addDistConstrain()`: Adds distance constraints between atoms to the system.

##### Multi-System Management
- `void setSystemReplica()`: Sets the current system replica for parallel or distributed simulations.
- `void countSystemReplica()`: Counts the number of system replicas available in the simulation.
- `void nextSystemReplica()`: Advances to the next system replica in a sequence.
- `void prevSystemReplica()`: Moves back to the previous system replica in a sequence.
- `void getMultiSystemPointers()`: Retrieves pointers to multiple system configurations for further processing or analysis.
- `void swith_method()`: Switches between different force field methods (e.g., MMFF, UFF).

##### Global Optimization
- `void runGlobalOptimization()`: Runs a global optimization process to find optimal configurations of the molecular structure.
- `void startExploring()`: Starts exploring new configurations in global optimization.
- `void stopExploring()`: Stops exploring configurations in global optimization.
- `void getMultiConf()`: Retrieves multiple configurations from the optimizer.

##### Surface Interaction
- `void scanSurfFF()`: Scans the surface with a force field, calculating forces at each point.

##### Data Upload and Evaluation
- `void upload_pop()`: Uploads population data for further processing or analysis.
- `void evalAFMscan()`: Evaluates an atomic force microscopy (AFM) scan on a grid shape.
- `void evalAFM_FF()`: Evaluates the force field on a grid shape using AFM techniques.

##### Information and Debugging
- `void getStatusString()`: Generates a string containing status information about the current state of the molecular world.
- `void printSwitches()`: Prints out the current state of various switches used in the molecular simulation to aid debugging or configuration checks.
- `void addSnapshot()`: Adds a snapshot of the current molecular structure to a database, optionally creating new entries if necessary.
- `void printDatabase()`: Prints detailed information about the molecular database, including stored configurations and their properties.
- `void computeDistance()`: Computes the distance between two atoms in the molecular system using data from the global optimization database.
- `void info_str()`: Provides detailed information about the current settings and parameters of the molecular world.

##### Electron Pair Handling
- `void hideEPairs()`: Hides electron pairs from further calculations or visualization.
- `void unHideEPairs()`: Unhides previously hidden electron pairs, restoring them to the system.

##### Molecular Geometry
- `void getGroupPose()`: Retrieves the pose and orientation of a group within the molecule.

##### Non-Bonded Molecule Handling
- `void initNBmol()`: Initializes non-bonded molecules with specified parameters and settings.
- `void loadNBmol()`: Loads a molecular structure from an XYZ file into the system.

##### Other
- `void setFromRef()`: Sets From Ref
- `void getTitle()`: Gets title.