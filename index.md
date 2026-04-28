# FireCore: A Comprehensive Overview

## 1. Project Goal and Vision

FireCore is an integrated simulation environment for on-surface chemistry and scanning-probe microscopy (SPM). It aims to provide a powerful tool for high-throughput screening of molecular configurations and processes on surfaces. The project combines quantum mechanics (QM) and molecular mechanics (MM) methods (QM/MM) with GPU acceleration to achieve a balance between computational speed and accuracy.

The primary goals of FireCore are:

*   **Exploration of non-covalent self-assembling processes:** Simulating how molecules arrange themselves on surfaces.
*   **Exploration of chemical reaction pathways:** Studying chemical reactions between molecules on a surface.
*   **Computational prototyping of molecular designs:** Designing molecules for pre-programmed self-assembly.
*   **Simulation of high-resolution AFM/STM experiments:** Simulating images from atomic force microscopy (AFM) and scanning tunneling microscopy (STM).

## 2. Project Structure

FireCore is not a single application but a collection of interconnected subprojects, libraries, and scripts. The main components are written in C++, Fortran, and Python, each with a specific role.

*   **`cpp/`**: The C++ core of the project, containing implementations of classical force fields, molecular manipulation tools, and visualization components.
*   **`fortran/`**: The Fortran-based Fireball DFTB implementation, providing the quantum mechanics capabilities.
*   **`pyBall/`**: The Python interface that glues together the C++ and Fortran components, providing a scripting and automation layer.
*   **`tests/`**: A collection of examples and test scripts, which also serve as a practical guide to the project's features.
*   **`doc/`**: A repository of documentation, development notes, and technical descriptions.

### 2.1. `cpp/`: The C++ Core

The `cpp/` directory is the heart of FireCore's classical mechanics and visualization capabilities. It contains high-performance C++ libraries and applications for molecular mechanics, visualization, and GPU acceleration.

**Key Subdirectories:**

*   **`common/`**: Contains the core computational libraries, including molecular mechanics algorithms, mathematical utilities, and data structures.
*   **`common_SDL/`**: Home to the SDL-based GUI and visualization components.
*   **`common_resources/`**: Shared data such as force field parameters, molecular structures, and OpenCL kernels.
*   **`apps/`**: Interactive applications, including the main `MolecularEditor`.
*   **`apps_OCL/` and `apps_CUDA/`**: Applications specifically accelerated with OpenCL and CUDA.
*   **`libs/`**: Compiled libraries for molecular mechanics and other functionalities.

**Core Functionalities:**

*   **Molecular Mechanics (MM):**
    *   **`MMFFsp3`**: A molecular mechanics force field for molecules with sp3 hybridization.
    *   **`GridFF`**: A grid-based force field to describe interactions with substrates.
    *   **`MMFFBuilder`**: A tool for automatic topology generation and charge assignment.
*   **Visualization:**
    *   The `MolecularEditor` application provides real-time 3D visualization of molecular dynamics, force fields, and QM/MM calculations.
*   **GPU Acceleration:**
    *   FireCore leverages OpenCL and CUDA to accelerate force evaluations and other computationally intensive tasks.

**Building the C++ Core:**

The C++ components are built using CMake. The general process is:

1.  Create a build directory (e.g., `cpp/Build`).
2.  Run `cmake` with appropriate options (e.g., `-DWITH_SDL=ON`, `-DWITH_OPENCL=ON`).
3.  Compile the code using `make`.

### 2.2. `fortran/`: The Quantum Core (Fireball)

The `fortran/` directory contains the implementation of the Fireball Density Functional Tight Binding (DFTB) code. This is the quantum mechanics (QM) engine of the FireCore project.

**Key Features:**

*   **DFTB Calculations:** Implements the semi-empirical DFTB method for efficient quantum calculations.
*   **Self-Consistent Field (SCF):** Solves for the electronic ground state of the system.
*   **Real-Space Grids:** Can project electron density and molecular orbitals onto real-space grids for visualization and analysis.
*   **QM/MM Interface:** Provides an interface to the C++ molecular mechanics code for QM/MM simulations.

**Directory Structure:**

The `fortran/` directory is organized into several subdirectories, each responsible for a specific part of the calculation:

*   **`MAIN/`**: The main program and SCF drivers.
*   **`MODULES/`**: Fortran modules and data structures.
*   **`ASSEMBLERS/`**: Hamiltonian and matrix assembly.
*   **`INTERACTIONS/`**: Two-center and three-center interactions.
*   **`GRID/`**: Real-space grid operations.

**Building the Fortran Core:**

The Fortran part of the project is built using a shell script (`make.sh`) that compiles the source code and links it with the Intel MKL library.

### 2.3. `pyBall/`: The Python Interface

The `pyBall/` directory provides the Python interface to the C++ and Fortran components of FireCore. It also contains standalone Python implementations of various computational methods and utilities.

**Core Modules:**

*   **`FireCore.py`**: The main interface to the Fortran Fireball DFT code.
*   **`MMFF.py`**: An interface to the C++ molecular mechanics libraries.
*   **`AtomicSystem.py`**: A module for high-level molecular structure manipulation.

**Specialized Modules:**

*   **`OCL/`**: Contains pure pyOpenCL implementations of molecular mechanics and other algorithms, allowing for GPU acceleration without relying on the C++ libraries.
*   **`DFT/`**: Provides high-level utilities for performing and managing DFT calculations.
*   **`GUI/`**: Contains Python-based GUI components.

**Utilities:**

*   **`plotUtils.py`**: A collection of tools for plotting and visualization.
*   **`grid_utils.py`**: Utilities for manipulating 3D grids.

**External Integrations:**

`pyBall` also includes interfaces to other popular quantum chemistry and molecular dynamics packages, such as:

*   Psi4
*   PySCF
*   DFTB+
*   LAMMPS

**Usage:**

To use the `pyBall` package, you need to set the `PYTHONPATH` environment variable to the root of the FireCore repository. The Python interface is designed to be intuitive and automatically handles the compilation of the necessary C++ libraries.

### 2.4. `tests/`: Examples and Getting Started

The `tests/` directory is the recommended starting point for new users. It contains a wide range of examples and test scripts that demonstrate the various features of FireCore. These scripts are not only for testing but also serve as practical tutorials on how to use the different components of the project.

Each subdirectory in `tests/` is dedicated to a specific feature or module, such as:

*   **`tMMFF*/`**: Tests for the molecular mechanics force fields.
*   **`tMolGUIapp*/`**: Tests for the GUI application.
*   **`Fireball/`**: Tests for the Fireball DFT code.

By exploring the `run.sh` scripts in these directories, you can learn how to compile and run the different parts of the project, and see how they are used to perform various types of simulations.

## 3. Key Features and Capabilities

FireCore offers a wide range of features for molecular simulation, including:

*   **Hybrid QM/MM Simulations:** Combine the accuracy of quantum mechanics (with Fireball DFTB) with the speed of classical mechanics (with C++ and OpenCL accelerated force fields).
*   **GPU Acceleration:** Leverages OpenCL and CUDA to accelerate simulations, enabling the study of larger systems and longer timescales.
*   **Interactive Visualization:** An SDL-based GUI allows for real-time visualization and manipulation of molecular structures and simulation data.
*   **Python Scripting:** A comprehensive Python API (`pyBall`) provides a flexible and powerful way to automate simulations, analyze data, and integrate with other tools.
*   **High-Resolution AFM Simulation:** The project includes tools for simulating high-resolution AFM images, bridging the gap between simulation and experiment.

## 4. Current Status and Future Development

FireCore is currently in a **work-in-progress** stage of development. While many features are implemented, the software is still under active development and may contain bugs or instabilities. The developers have made it available to the public to promote transparency and open-source collaboration.

Future development will focus on:

*   **Improving the classical force fields:** Further parameterization and testing of the existing force fields.
*   **Integration with other packages:** Expanding the interfaces to other quantum chemistry and molecular dynamics codes.
*   **Enhancing the GUI:** Improving the user-friendliness and functionality of the graphical interface.
*   **Bug fixing and stabilization:** Addressing known issues and improving the overall stability of the software.

If you are interested in using or contributing to FireCore, the developers encourage you to get in touch with them.





