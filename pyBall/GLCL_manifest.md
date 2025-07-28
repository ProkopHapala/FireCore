# GLCL Simulation Framework - Detailed Plan

## Goal
Develop a flexible simulation environment that integrates OpenGL and OpenCL, configured via JSON files, with a dynamic PyQt5 GUI. The environment should abstract the management of kernels/shaders, buffers, and uniforms/parameters.

## Key Components and Files

### 1. `GLCL_manifest.md` (This File)
*   **Purpose**: Central documentation for the new GLCL simulation framework, detailing its architecture, components, and usage.

### 2. `OGLsystem.py` (New File)
*   **Purpose**: To encapsulate all OpenGL-related abstractions and utility functions. This file will handle shader compilation, program linking, buffer management (VBO, EBO, VAO), texture handling, and uniform setting. It should be independent of PyQt.
*   **Inspiration/Code to Copy**: 
    *   From `/home/prokophapala/git/FireCore/pyBall/GUI/GLGUI.py`:
        *   `class GLobject`: Adapt this structure for managing OpenGL objects (VAO, VBO, EBO, textures). Its methods like `init`, `draw`, `update_vbo`, `update_ebo` will be crucial.
        *   Utility functions like `upload_buffer`, `setup_ebo`, `setup_vbo` (or similar logic) should be extracted and generalized.
        *   Shader compilation and linking logic (e.g., `compileShader`, `linkProgram`).
    *   From `/home/prokophapala/git/FireCore/pyBall/GUI/GLSL_Simulation.py`:
        *   The concept of `GLSL_Program` and `GLSL_Simulation` classes, particularly how they manage shaders, uniforms, and textures, will be a good starting point for structuring the OpenGL program management.
*   **Key Components**: 
    *   Functions for creating/compiling shaders and linking programs.
    *   Functions for creating/updating VBOs, EBOs, VAOs.
    *   A class (e.g., `GLProgram` or `ShaderManager`) to manage a compiled GLSL program and its uniforms/textures.
    *   A class (e.g., `GLBuffer`) to abstract OpenGL buffer objects.

### 3. `OCLsystem.py` (New File)
*   **Purpose**: To encapsulate all OpenCL-related abstractions and utility functions. This file will handle OpenCL context and command queue creation, kernel compilation, buffer creation, and kernel argument setting. It should be independent of PyQt.
*   **Inspiration/Code to Copy**: 
    *   From `/home/prokophapala/git/FireCore/pyBall/OCL/OpenCLBase.py`:
        *   `init_cl()`: For OpenCL context and command queue initialization.
        *   `try_make_buff()`: For creating OpenCL buffers.
        *   `toGPU()`, `fromGPU()`: For data transfer between host and device.
        *   `callKernel()`: For executing OpenCL kernels.
    *   From `/home/prokophapala/git/FireCore/pyBall/OCL/HubbardSolver.py`:
        *   Observe how kernels are loaded and arguments are set (e.g., `init_cl_kernels`, `run_cl_kernels`).
*   **Key Components**: 
    *   Functions for OpenCL context and queue setup.
    *   Functions for compiling OpenCL kernels from source files.
    *   A class (e.g., `CLKernel` or `KernelManager`) to manage OpenCL kernels and their arguments.
    *   A class (e.g., `CLBuffer`) to abstract OpenCL memory objects.

### 4. `GLCLGUI.py` (New File)
*   **Purpose**: To integrate `OGLsystem.py` and `OCLsystem.py` within a PyQt `QOpenGLWidget`. This file will define the core `GLCLWidget` that handles the rendering loop, simulation updates, and manages the OpenGL and OpenCL resources for a specific simulation. It should primarily focus on the OpenGL rendering and OpenCL computation, with minimal direct GUI elements.
*   **Inspiration/Code to Copy**: 
    *   From `/home/prokophapala/git/FireCore/pyBall/GUI/GLGUI.py`:
        *   The structure of `class GLView(QOpenGLWidget)` and its methods like `initializeGL`, `paintGL`, `resizeGL`. This will be the base for our `GLCLWidget`.
        *   Event handling methods (e.g., `mousePressEvent`, `mouseMoveEvent`) can be adapted if needed for interaction within the GL widget.
*   **Key Components**: 
    *   `class GLCLWidget(QOpenGLWidget)`:
        *   Initializes OpenGL and OpenCL systems.
        *   Manages the simulation state.
        *   Calls rendering functions from `OGLsystem.py` and computation functions from `OCLsystem.py`.
        *   Handles the main loop for simulation and rendering.

### 5. `GLCLBrowser.py` (New File)
*   **Purpose**: This will be the main application entry point, inheriting from `BaseGUI`. It will dynamically load simulation configurations from JSON files, instantiate the `GLCLWidget`, and create dynamic GUI controls for simulation parameters.
*   **Inspiration/Code to Copy**: 
    *   From `/home/prokophapala/git/FireCore/pyBall/OGL/BaseGUI.py`:
        *   `class BaseGUI(QMainWindow)`: This will be the base class.
        *   `populate_params_from_json` (or `populate_params_from_dict`): This method is crucial for dynamically creating GUI widgets from JSON-defined parameters. We will reuse this directly.
    *   From `/home/prokophapala/git/FireCore/pyBall/GUI/GLSL_GUI.py`:
        *   Observe how `GLSL_GUI` integrates `GLSL_Simulation` and uses `populate_params_from_json` to create controls. This will guide the structure of `GLCLBrowser`.
*   **Key Components**: 
    *   `class GLCLBrowser(BaseGUI)`:
        *   Loads JSON configuration files.
        *   Instantiates `GLCLWidget` and sets it as the central widget.
        *   Calls `populate_params_from_json` to create GUI controls for simulation parameters.
        *   Provides mechanisms to load different simulations.

### 6. JSON Configuration Structure (Example: NBody Simulation)
*   **Purpose**: To define the simulation's OpenCL kernels, OpenGL shaders, buffers, and parameters in a declarative way.
*   **Proposed Structure**:
    ```json
    {
        "name": "NBody Simulation",
        "description": "A simple N-body simulation using OpenCL for physics and OpenGL for rendering.",
        "parameters": {
            "dt": ["float", 0.01, 0.001, 0.1], 
            "damping": ["float", 0.99, 0.9, 1.0],
            "num_particles": ["int", 1024]
        },
        "buffers": {
            "positions": {"size": "num_particles * 4 * 4", "type": "float", "usage": "read_write"},
            "velocities": {"size": "num_particles * 4 * 4", "type": "float", "usage": "read_write"}
        },
        "opencl_kernels": [
            {"name": "update_positions", "file": "nbody_update.cl", "entry_point": "update_positions_kernel", "args": ["positions", "velocities", "dt"]},
            {"name": "calculate_forces", "file": "nbody_forces.cl", "entry_point": "calculate_forces_kernel", "args": ["positions", "velocities", "damping"]}
        ],
        "opengl_shaders": [
            {"name": "particle_render", "vertex_file": "particle.vert", "fragment_file": "particle.frag", "uniforms": ["positions", "view_matrix", "projection_matrix"]},
            {"name": "debug_lines", "vertex_file": "line.vert", "fragment_file": "line.frag", "uniforms": ["positions", "view_matrix", "projection_matrix"]}
        ],
        "pipeline": [
            {"type": "opencl", "kernel": "calculate_forces", "global_size": "num_particles"},
            {"type": "opencl", "kernel": "update_positions", "global_size": "num_particles"},
            {"type": "opengl", "shader": "particle_render", "draw_mode": "points", "count": "num_particles"},
            {"type": "opengl", "shader": "debug_lines", "draw_mode": "lines", "count": "num_particles"}
        ]
    }
    ```
*   **Note**: The `size` in buffers and `global_size` in kernels can reference `parameters` for dynamic sizing.

### 7. NBody Test Case (JSON, Shaders, Kernels)
*   **Purpose**: To demonstrate the new framework with a concrete example.
*   **Files**:
    *   `simulations/nbody.json`: The configuration file for the NBody simulation.
    *   `shaders/nbody/particle.vert`, `shaders/nbody/particle.frag`: OpenGL shaders for rendering particles.
    *   `cl/nbody/nbody_update.cl`, `cl/nbody/nbody_forces.cl`: OpenCL kernels for NBody physics.
*   **Inspiration**: `/home/prokophapala/git/FireCore/pyBall/GUI/NBody_glcl.py` will be the source for the actual physics and rendering logic to be translated into OpenCL kernels and OpenGL shaders.

## Development Steps:
1.  Create `GLCL_manifest.md` (this file) and write this plan into it.
2.  Implement `OGLsystem.py` with basic shader/program/buffer management.
3.  Implement `OCLsystem.py` with basic context/queue/kernel/buffer management.
4.  Implement `GLCLGUI.py` defining `GLCLWidget` that integrates `OGLsystem` and `OCLsystem`.
5.  Implement `GLCLBrowser.py` inheriting from `BaseGUI`, loading JSON, and managing the `GLCLWidget` and dynamic controls.
6.  Create the `nbody.json` configuration file.
7.  Translate the NBody logic from `/home/prokophapala/git/FireCore/pyBall/GUI/NBody_glcl.py` into OpenCL kernels and OpenGL shaders, placing them in the `cl/nbody` and `shaders/nbody` directories.
