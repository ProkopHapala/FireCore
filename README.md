# FireCore

FireCore is an integrated simulation environment dedicated to on-surface chemistry and scanning-probe microscopy (SPM). The code aims to streamline high-throughput screening of molecular configurations and processes on the surface of inorganic crystals, including:

* exploration of non-covalent self-assembling processes
* exploration of chemical reaction pathways between organic molecules on surface
* computational prototyping of molecular designs for pre-programmed self-assembling
* simulation of high-resolution AFM / STM experiments of organic molecules

To achieve a favorable balance between speed and accuracy the code use combination of both ab-initio quantum methods and classical forcefield (QM/MM approach) and extensive GPU acceleration.

![Schematic illustration of different aspects of on-surface chemistry simulations, and methods used for their efficient description.](Software_Schematic.png?raw=true "Schematic illustration of different aspects of on-surface chemistry simulations, and methods used for their efficient description.")

# Structure of the package

FireCore consists of high-performance simulation modules written in Fortran, C/C++ & OpenCL integrated using a common Python3 interface. 

## 1. Quantum solution of electronic & geometric structure of molecules

A streamlined version of the Fireball self-consistent local-orbital ab-initio tight-binding molecular dynamics code (https://github.com/fireball-QMD/progs) is utilized for the following tasks:

* Optimization of the electronic structure of molecules and determination of molecular orbitals and density.
* Optimization of the geometric structure of molecules (the quantum component of QM/MM).
* Projection of Molecular Orbitals and Electron density to a real-space grid.

Fireball achieves rapid solutions to self-consistent electronic structure problems by utilizing a streamlined, optimized numerical basis set and employing precomputed interaction integrals. These integrals are efficiently interpolated during the construction of the Hamiltonian.

## 2. Classical Molecular Mechanics

FireCore is currently integrated with homebrew classical-forcefields acalerated using [OpenMP](https://github.com/ProkopHapala/FireCore/blob/4975228f970bc7b5aa560f86874120442f8b4d36/cpp/common/molecular/MolWorld_sp3.h#L1090) and [OpenCL](https://github.com/ProkopHapala/FireCore/blob/master/cpp/common_resources/cl/relax_multi.cl) to leverage both multi-core CPU and GPU. The main purpose is to reduce the computational cost required for the relaxation of organic molecules and the description of the interaction between organic molecules and substrate. These codes comprise:

* **[MMFF_sp3](https://github.com/ProkopHapala/FireCore/blob/master/cpp/common/molecular/MMFFsp3_loc.h)** - a simple homebrew classical molecular mechanics code with fixed boding topology (springs on bond length, angles and dihedrals)
    *  Code includes [utilities to automatically generate bonding topology](https://github.com/ProkopHapala/FireCore/blob/master/cpp/common/molecular/MMFFBuilder.h) from atomic coordinates and assign atomic charges by the charge-equilibration algorithm (depending on electronegativity and local polarization) 
* **[GridFF](https://github.com/ProkopHapala/FireCore/blob/master/cpp/common/molecular/GridFF.h)** - grid-based classical forcefield to describe the interaction of organic molecules with the rigid substrate or other rigid bodies (e.g. tip of AFM?). The code projects non-covalent forces from rigid objects (e.g. electrostatics, van der Waals attraction and Pauli repulsion) onto a real-space grid. Then during molecular-dynamics simulation, we read from the grid forces acting on each atom of the flexible (non-rigid) molecules deposited on the rigid substrate.

In future we plan integration with state-of-the-art classical forcefield packages such as LAMMPS.

## 3. High-resolution AFM simulations

The GPU-accelerated [Probe Particle Model](https://github.com/Probe-Particle/ppafm) is employed for simulating high-resolution AFM images using a flexible tip, such as a CO molecule or Xe atom. Currently, the integrated model in FireCore's GUI employs simple [original approach](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.085421) with a Lennard-Jones potential and point charges accelerated using GPU (OpenCL) to achieve interactive simulation speed.

We are working on implementation of *[full density-based model (FDBM)](https://pubs.acs.org/doi/full/10.1021/acsnano.8b08209)* using Pauli repulsion and electrostatic interaction potential from the density provided by the Fireball-DFT module through a seamless QMMM interface, all within less than one second for typical molecules like PTCDA. This capability will allow for swift screening of candidate molecular structures to align with experimental data, or for the rapid generation of extensive databases for training machine-learning models for AFM imaging.

# Instalation (compile & run)

To install all library dependencies. On Ubuntu 22.04 this can be done by running:
```
sudo apt-get install cmake g++ gfortran intel-mkl libmkl_intel_lp64 libmkl_intel_core libmkl_intel_thread libsdl2-dev libsdl2-image-dev nvidia-opencl-dev libclfft-dev python3-numpy python3-matplotlib
```
however individual parts can be installed independently. For example, it is possible to use DFT-Fireball without OpenCL or OpenGL. Therefore dependencies can be split as follows:

* C++ modules (Classical Force-fields & Molecular Manipulation tools): 
    * (without SDL or OpenCL): `cmake g++`
    * GUI 3D Graphics (SDL & OpenGL): `libsdl2-dev libsdl2-image-dev`
    * GPU acceleration (OpenCL): `nvidia-opencl-dev libclfft-dev` 
* Fortran (Fireball-DFT): `gfortran intel-mkl libmkl_intel_lp64 libmkl_intel_core libmkl_intel_thread`
* Python API: `python3-numpy python3-matplotlib`

## C/C++ Modules (classical forcefields, molecular manipulation, grids, OpenCL, GUI)

Most of estential parts of the package are written in C++, icluding MMFF, gridFF, topology generator, GUI, and OpenMP & OpenCL accelerated solvers. 

1. install C/C++ compiler: `g++ cmake`
2. install libraries:
    * for GUI 3D Graphics (SDL & OpenGL): `libsdl2-dev libsdl2-image-dev`
    * GPU acceleration (OpenCL): `nvidia-opencl-dev libclfft-dev` 
3. Create a build directory in `cpp/Build` and navigate into it
4. Configure cmake project: from inside `cpp/Build` 
5. run     `cmake .. -DWITH_SDL=ON -DWITH_OPENCL=ON`
    * options `-DWITH_SDL=OFF` and `-DWITH_OPENCL=OFF` can be eventually switched to reduce compilation time and/or limit library dependencies
6. compile by running `make` inside `cpp/Build` directory
7. Run tests:
    * Forcefield with visual GUI go to `/home/prokophapala/git/FireCore/cpp/sketches_SDL/Molecular` and run e.g. `./test_MMFFsp3` or `test_RARFFarr` or `test_SoftMolecularDynamics`
    * to run Visual GUI for QMMM dynamics (assuming both C++ and Fortran modules are compiled) got to `FireCore/tests/tQMMM_diacetylene` and run `./run.sh`

## Fortran modules (Fireball-DFT)

Only Fortran module is Fireball-DFT program. Its installation is optional. It is currently required only for operation of QMMM modules and plotting of molecular orbitals and FDBM AFM simulations. 

1. install Fortran compiler: `gfortran`
2. install intel mkl libraries and other dependencies: `intel-mkl libmkl_intel_lp64 libmkl_intel_core libmkl_intel_thread`
3. create directory `FireCore/Build`
4. run `FireCore/make.sh` which will compile source codes from `fortran` directory into `Build` directory
5. run test:
    * navigate to `FireCore/tests`
    * Make sure that `Fdata_HC_minimal` is present, otherwise download it here: https://fireball-qmd.github.io/
    * got to `tests/t01_H2` and run `./run.sh`

## Python API Binding

Python provides API to both Fireball-DFT as well as classical forcefieds implemented in C++ and OpenCL. It also provides some additional utilities for manipulating, post-processing and plotting molecular structures. All these modules are in `pyBall` package. No instalation is required. It is only recomanded to set python path evironment variable. 

* Set environement variable, e.g. `export PYTHONPATH=/home/prokop/git/FireCore/:PYTHONPATH`
* Try classical forcefield library `MMFFsp3` go tp `FireCore/tests/tMMFF` and run `./run.sh`
* Try Fireball-DFT in `FireCore/tests/FitFF` by running `./run.sh`
    
# WARNING: FireCore is Work-in-progress stage of development

This scientific software is still under active development and has not undergone comprehensive testing for bugs and stability. We have chosen to share it on GitHub in its current work-in-progress state to uphold the principles of transparency and open-source development, especially for software funded by public research organizations. However, it is not yet suitable for general use in scientific production without consulting the developers first.

If you plan to use this software for your scientific application, please reach out to us initially to discuss the readiness of specific sub-modules for your particular scientific application. We are eager to collaborate on potential applications of this code, as well as its further development. We particularly encourage bug reports and feature suggestions, encompassing both the numerical and scientific aspects of computational models, as well as practical considerations like user interface and potential integration with other software. If you have any feedback or suggestions of this nature, please submit a new **issue** in the *Issues* section on the GitHub page, or contact the main developer at ProkopHapala@gmail.com.

## Known issues

* **Computational models**
  * The classical forcefield is currently undergoing parameterization and testing. Presently, we have implemented a simplified version of the [UFF forcefield](https://pubs.acs.org/doi/10.1021/ja00051a040). However, it has not yet undergone rigorous testing against any reference. While the forcefield currently converges to reasonably optimized molecular geometries, the energetics (including bond stiffness and intramolecular angles) may still deviate significantly from physical reality.
  * Some molecules may experience drift in the evaluation of non-covalent interactions (i.e., non-zero force on the center of mass and torque) when using GPU accelerated forcefields implemented with OpenCL. To resolve this problem please switch to OpenMP parallelization (with command line argument `-iParalel 1`) or switch of paralelization completely ( `-iParalel 0`). While we though we already resolved this issue, it seems to have resurfaced recently, at least in certain situations. We are actively working to understand the root cause of this problem and rectify it.
  * OpenMP parallelization (`OMP_NUM_THREADS >1`) appears to conflict with the operation of MKL LAPACK used for matrix diagonalization in the Fireball-DFT program (similar to [this](https://community.intel.com/t5/Intel-oneAPI-Math-Kernel-Library/Incorrect-eigenvectors-from-ZHEEV-in-MKL-using-multi-threading/td-p/836860)). As a result, the SCF cycle may fail to converge, leading to nonsensical QMMM forces and Molecular Orbitals. This issue can be mitigated by setting `export OMP_NUM_THREADS 1` or compiling with OpenMP disabled (`-DWITH_OPENMP=OFF`). However, this approach significantly limits the program's performance, especially on modern machines equipped with multiple CPU cores. We are actively investigating the precise cause of this problem and seeking a resolution.
* **User interface**
  * It is relatively easy to encounter program crashes due to missing input files, incorrect input formats, or even pressing buttons in the GUI out of sequence. We are striving to handle these situations and provide appropriate error messages and warnings, but due to limited manpower and rapid changes of the software we are still far handling all situations.
  * The grids in the GUI interface may not align precisely. For instance, the visualization of molecular orbitals and AFM simulation grids may be slightly offset in relation to the position of the molecule. This discrepancy arises from these quantities being calculated in different coordinate systems, and we perform multiple transformations to align them. We need to conduct a thorough examination of the definitions and mutual relationships between these coordinate systems to achieve precise alignment.
  * While the GUI currently includes widgets (such as buttons, input boxes, and drop-down lists), these widgets may not always be in sync with the latest program functionality (which is evolving rapidly). Therefore, it is possible that interactions with these widgets may (i) have no effect, (ii) produce unexpected or nonsensical results, or (iii) cause the program to crash. The only part of the user interface currently kept relatively up-to-date with the computational core are the **keyboard shortcuts** and mouse interactions (e.g., selecting atoms, dragging the molecule). Therefore keyboard shortcuts are currently considered as main mode of interaction with the GUI. This issue is primarily due to limited resources, with our main focus on the computational core rather than the user interface. 
    * We would greatly appreciate collaboration on refining the user interface, not only from a programming perspective but more importantly from an ergonomic and beta-testing standpoint. 
  * The Python API for calling FireCore library functions may lag behind the development of the GUI. This is because, for the development of new features, we typically use the GUI where we can easily visualize debug information (such as molecular geometry, atomic types, forces, charges, etc.). Only once this is finalized do we expose the functionality to the Python interface. Therefore, it is possible that the Python API may expose calls to functions that are not up-to-date or have been deprecated, or otherwise deviate from the latest state of the computational core.      