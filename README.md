# FireCore

FireCore is an integrated simulation environment dedicated to on-surface chemistry and scanning-probe microscopy (SPM). The code aims to streamline high-throuhput screening of molecular configurations and processes on surface of inorganic crystals, including:

* exploration of non-covalent self-assembling processes
* exploration of chemical reaction pathways between organic molecules on surface
* computational prototyping of molecular desings for pre-programed self-assembling
* simulation of high-resolution of AFM / STM experiments of organic molecules

To achieve favourable compromise between speed and accuracy the code use combination of both ab-initio quantum methods and classical forcefield (QM/MM approach) and extensive GPU acceleartion.

![Schematic illustration of different aspects of on-surface chemistry simulations, and methods used for their efficient description.](Software_Schematic.png?raw=true "Schematic illustration of different aspects of on-surface chemistry simulations, and methods used for their efficient description.")

# Structure of the package

FireCore consist of high-performance simulation modules written in Fortran, C/C++ & OpenCL integrated using common Python3 interface. The modules are based on modified version of other simulation packages.

## 1. quantum solution of electronic & geometric structure of molecules

Striped down version of Fireball self-cnsistent local-orbital ab-initio tight-binding molecular dynamics code (https://github.com/fireball-QMD/progs) is used for following tasks:

* To optimize electronic structure of molecules and find molecular orbitals and density
* To optimize geometric structure of molecules (quantum part of QM/MM)
* To project Molecular Orbitals, and Electron density to real space grid

Fireball allows fast solution of self-consistent electronic structure problem thanks to small optimized numerical basiset and pre-calculated table of integrals.

## 2. Classical Molecular Mechanics

FireCore is currently integrated with homebrew classical-forcefields codes implemented in SimpleSimulationEngine (https://github.com/ProkopHapala/SimpleSimulationEngine). The main purpose is to reduce computational cost required for relaxation of organic molecules and description of interaction between organic molecule and substrate. These codes comprise of:

* **MMFF** - a simple homebrew classical molecular mechanics code with fixed boding topology (springs on bond length, angles and dihedrals)
    *  Code includes utilities to automatically generate bonding topology form atomic coordinates and assign atomic charges by charge-equalibration algorithm (depending on electronegativity and local polarization) 
* **RFF** - a simple homebrew reactive foce-field to roughly simulate chemical interactins between organic species (H,C,O,N) including formation and breaking of new bonds and changes of bond-order
* **gridFF** - grid based classical forcefield to describe interaction of organic molecules with rigid substrate or other rigid bodies (e.g. tip of AFM?). The code projects non-covalent forces from rigid objects (e.g. electrostatics, van der Waals attraction and Pauli repulsion) onto real-space grid. Then during molecular-dynamics simulation we read from the grid forces acting on each atom of the flexible (non-rigid) molecules deposited on the rigid substrate.

In future we plan integration with state-of-the art classical forcefield packages such as LAMMPS.

## 3. High-resolution AFM simulations

GPU accelerated Probe Particle Model (https://github.com/ProkopHapala/ProbeParticleModel) is used for simulation of high-resolution of AFM images using flexible tip (e.g. CO molecule or Xe atom). It was shown that variation of local electrostatic field and electron density (e.g. bond-order and free electron pairs) plays an importaint role in formation of HR-AFM contrast and its distortions. To obtain good estimate of electron density and electrostatic potential we use self-consistent electronic structure calculated using fireball as an input. With combination of classical forcefield relaxation of molecular geometry, fireball solution of self-consistent electronic structure and GPU accelerated AFM relaxation the FireCore is able do conduct series of volumetric AFM simulation (starting from molecular topology) in less then 1 second. This allows rapid screend of candidate molecular structures to match experimental data, or to rapidly generate large database for training machine-learning models for AFM imaging.

We use GPU accleration using OpenCL kernells for most performance intensive steps of the simulation, including:

* Simulation of AFM scan where probe-particle (CO apex) relax in Volumgetric potential (gridFF) comprising of electrostatics, van der Waals attraction and Pauli repulsion from rigid sample.
* Evaluation of gridFF by projection of comtributions from each atom of the sample - comprising van der Waals attraction and Pauli.
* Evaluation of electrostatic potential and forcefield on grid from electron density (i.e. solution of Poisson equation) using Fast-Fourier-Trasfrom. We use library clFFT (https://github.com/clMathLibraries/clFFT) for the purpose.
* Projection of wave-functions or electron density density form FireBall onto real space grid using FireBall numerical basiset. 
# WARRNING: FireCore is Work-in-progess stage of developement

