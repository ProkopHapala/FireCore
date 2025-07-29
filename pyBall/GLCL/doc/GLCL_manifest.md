# Mission: A "Shadertoy" for Scientific Simulation

The primary goal of the GLCL project is to create a powerful, yet easy-to-use browser for developing and visualizing scientific simulations. Inspired by Shadertoy, this framework aims to eliminate boilerplate code related to GUI, OpenGL/OpenCL context management, and data handling.

This allows users—researchers, students, and hobbyists—to focus purely on their core algorithms, expressed through:
1.  **OpenCL kernels** for computation.
2.  **GLSL shaders** for visualization.
3.  A simple **Python configuration script** to wire everything together.

The browser dynamically loads and executes these user-provided components, providing an interactive environment for rapid prototyping and exploration of complex physical systems.

---

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
*   **Current Status/Progress**:
    *   **OpenGL Uniform and VBO Handling**: Fixed `GLError: invalid operation` by ensuring `glUseProgram` is called before `glGetUniformLocation` and `glUniform4fv` in `initializeGL`. Corrected `TypeError` in `glBindBuffer` by ensuring `glGenBuffers` and `glGenVertexArrays` are called only once within `initializeGL` (by initializing `self.gl_vbo` and `self.vao` to `0` and checking against `0`). Added a check in `update_particle_vbo` to only update if `self.gl_vbo` is valid.
    *   **Code Clean-up**: Removed redundant local OpenGL imports within methods to ensure consistent global imports.
    *   **Deferred OpenGL Initialization**: Moved `setup_particle_vbo()` calls from `set_particle_data()` to `initializeGL()`. This ensures OpenGL functions are called only after the OpenGL context is properly set up, resolving `OpenGL.error.Error: Attempt to retrieve context when no valid context`.
    *   **Particle Data Handling**: `set_particle_data()` now primarily stores data and triggers a repaint, with the actual VBO setup occurring in `initializeGL`.
*   **Inspiration/Code to Copy**: 
    *   From `/home/prokophapala/git/FireCore/pyBall/GUI/GLGUI.py`:

### 5. `GLCLBrowser.py` and `GLCLBrowser_bak.py` (New Files)
*   **Purpose**: Main application entry point and GUI for the GLCL simulation framework. Manages the overall application flow, including loading simulation configurations (JSON), setting up the `GLCLWidget`, and handling user interactions. `GLCLBrowser_bak.py` is the currently active development version.
*   **Current Status/Progress**:
    *   **Shader Source Loading Timing**: Resolved the "Warning: Shader sources not provided to GLCLWidget." by ensuring `apply_simulation_config` is called early in the `__init__` method, before `GLCLWidget.initializeGL` is triggered, to provide shader sources in time.
    *   **GUI Initialization Order**: Resolved `AttributeError` by ensuring GUI layout components (e.g., `self.params_layout`) are initialized before being accessed by `apply_simulation_config`.
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

### 6. Python Script Configuration Structure (Example: NBody Simulation)
*   **Purpose**: To define the simulation's OpenCL kernels, OpenGL shaders, buffers, and parameters, as well as initial data, using a Python script.
*   **Proposed Structure**:
    A Python script (e.g., `nbody.py`) will define a `config` dictionary and an `init()` function.
    The `config` dictionary will contain:
    *   `simulation_name`: Name of the simulation.
    *   `description`: Description of the simulation.
    *   `parameters`: Dictionary of simulation parameters (e.g., `dt`, `particle_count`). Each parameter can specify its value, type, and step for GUI controls.
    *   `buffers`: Dictionary defining the buffers, including their size, stride, and type.
    *   `opencl_source`: List of OpenCL kernel source files.
    *   `kernels`: Dictionary mapping kernel names to their local size, global size, associated buffers, and parameters.
    *   `kernel_pipeline`: List defining the sequence of OpenCL kernels to be executed per frame.
    *   `opengl_shaders`: Dictionary mapping shader program names to their vertex and fragment shader files, and a list of uniforms.
    *   `render_pipeline`: List defining the sequence of OpenGL render operations (draw calls), specifying the shader, element count, vertex buffer, and optional index buffer.

    The `init()` function will be a callback that returns a dictionary of initial data for the buffers. This allows for arbitrary initialization logic.

*   **Example (`nbody.py` simplified structure)**:
    ```python
    import numpy as np

    config = {
        "simulation_name":  "NBody Simulation",
        "description":      "NBody simulation using OpenCL for computation and OpenGL for rendering.",
        "parameters": {
            "particle_count": (2048,  "int" , 1      ),
            "dt":             (0.001, "float", 0.0001 )
        },
        "buffers":{
            "positions":  (2048, 4, "f4"),
            "velocities": (2048, 4, "f4")
        },
        "opencl_source": ["nbody_sim.cl"],
        "kernels": {
            "nbody_sim" : ( (32,), ("particle_count"), ["positions", "velocities"], ["dt"] )
        },
        "kernel_pipeline": ["nbody_sim"],
        "opengl_shaders": {
            "nbody_render" : ("points.glslv", "monocolor.glslf", ["positions"])
        },
        "render_pipeline":   [
            ( "nbody_render", "particle_count", "positions",     None ),
        ]
    }

    def init():
        particle_count = config["particle_count"]
        positions      = np.random.rand(particle_count, 4) * 2 - 1
        velocities     = np.random.rand(particle_count, 4) * 0.1
        return {
            "positions":  positions.astype(np.float32),
            "velocities": velocities.astype(np.float32)
        }
    ```
*   **Note**: This approach provides maximum generality, allowing for complex initialization logic and dynamic definition of simulation components without hardcoding names or structures.

### 7. NBody Test Case (JSON, Shaders, Kernels)
*   **Purpose**: To demonstrate the new framework with a concrete example.
*   **Files**:
    *   `simulations/nbody.json`: The configuration file for the NBody simulation.
    *   `shaders/nbody/particle.vert`, `shaders/nbody/particle.frag`: OpenGL shaders for rendering particles.
    *   `cl/nbody/nbody_update.cl`, `cl/nbody/nbody_forces.cl`: OpenCL kernels for NBody physics.
*   **Inspiration**: `/home/prokophapala/git/FireCore/pyBall/GUI/NBody_glcl.py` will be the source for the actual physics and rendering logic to be translated into OpenCL kernels and OpenGL shaders.

## Development Steps:
1.  Create `GLCL_manifest.md` (this file) and write this plan into it. **(Completed)**
2.  Implement `OGLsystem.py` with basic shader/program/buffer management. **(Completed - Basic implementation for NBody)**
3.  Implement `OCLsystem.py` with basic context/queue/kernel/buffer management. **(Completed - Significant bug fixes in argument parsing)**
4.  Implement `GLCLGUI.py` defining `GLCLWidget` that integrates `OGLsystem` and `OCLsystem`. **(Completed - Fixed OpenGL context initialization)**
5.  Implement `GLCLBrowser.py` inheriting from `BaseGUI`, loading JSON, and managing the `GLCLWidget` and dynamic controls. **(Completed - Automatic JSON loading, GUI parameter binding)**
6.  Create the `nbody.json` configuration file. **(Completed - `nbody_simulation.json` created)**
7.  Translate the NBody logic from `/home/prokophapala/git/FireCore/pyBall/GUI/NBody_glcl.py` into OpenCL kernels and OpenGL shaders, placing them in the `cl/nbody` and `shaders/nbody` directories. **(Completed - `nbody_sim.cl`, `nbody_vertex.glsl`, `nbody_fragment.glsl` created and integrated)**

## Current Status and Progress Summary

We have made significant progress in integrating the NBody simulation into the PyOpenCL-OpenGL framework with a PyQt5 GUI. The application now launches without immediate console errors, and the core components are in place.

### Key Changes and Rationale:

1.  **`GLCLBrowser.py` Refinements:**
    *   **Automatic JSON Loading:** Modified to accept a JSON filepath argument, enabling automatic loading of `nbody_simulation.json` on startup. This streamlines testing and configuration.
    *   **Improved GUI Layout:** Redesigned the main window to feature a large OpenGL viewport and a narrow control panel, enhancing usability.
    *   **Dynamic Parameter Binding:** Implemented `update_sim_uniforms` to dynamically link GUI spin box values (e.g., `dt`, `particle_count`) to OpenCL kernel parameters, allowing real-time simulation control.
    *   **Robust Error Handling:** Removed silent `try-except` blocks during JSON loading and initialization to ensure errors propagate, aiding debugging.

2.  **`BaseGUI.py` Adjustments:**
    *   **`populate_params_from_json` Fixes:** Corrected a critical bug where `populate_params_from_json` passed a list instead of a scalar float to `QDoubleSpinBox.setValue`, resolving `setValue` errors.
    *   **Layout Correction:** Fixed improper `QFormLayout` usage by ensuring widgets are added with `addRow` instead of `addLayout`.

3.  **`GLCLGUI.py` OpenGL Context Fixes:**
    *   **Deferred OpenGL Initialization:** Moved `setup_particle_vbo()` calls from `set_particle_data()` to `initializeGL()`. This ensures OpenGL functions are called only after the OpenGL context is properly set up, resolving `OpenGL.error.Error: Attempt to retrieve context when no valid context`.
    *   **Particle Data Handling:** `set_particle_data()` now primarily stores data and triggers a repaint, with the actual VBO setup occurring in `initializeGL`.

4.  **`OCLsystem.py` Robustness Improvements:**
    *   **`create_buffer` Enhancement:** Modified `create_buffer` to accept an optional `hostbuf` argument, allowing direct copying of host data to the OpenCL buffer during creation, resolving `INVALID_HOST_PTR` errors.
    *   **Kernel Argument Parsing Fixes:**
        *   Refined `_extract_kernel_headers` to accurately extract only the kernel signature (up to the first `{`), preventing the parser from misinterpreting kernel body code as arguments.
        *   Improved `_parse_kernel_header` to robustly distinguish between `__global` pointers (buffers) and scalar arguments, resolving the `ValueError: Required buffer 'particle_count' for kernel 'nbody_sim' is not allocated or is zero-sized.`

### Known Issues / Next Steps:

*   **Visual Verification:** The application launches without console errors, but visual confirmation of the NBody simulation rendering and dynamic parameter updates is still required.
*   **OpenCL-OpenGL Interoperability:** Implement proper OpenCL-OpenGL buffer sharing for performance optimization (currently, data is copied).
*   **Further Refinement:** Continue modularizing and documenting the framework for easier extension and maintenance.

## Progress in Current Session

This session focused on debugging the persistent shader compilation and linking issues:

1.  **Shader Compilation Error Resolution:**
    *   **Missing OpenGL Imports**: Resolved `NameError: name 'glGetShaderiv' is not defined` by adding `glGetShaderiv`, `GL_COMPILE_STATUS`, and `glDetachShader` to imports in `OGLsystem.py`.
    *   **Error Log Decoding**: Corrected `AttributeError: 'str' object has no attribute 'decode'` by removing erroneous `.decode('utf-8')` calls, as PyOpenGL's info logs are already strings.
    *   **Improved Error Reporting**: Enhanced `OGLsystem.py` to explicitly print shader and program info logs before raising `RuntimeError` to ensure any available error messages are visible.

2.  **Diagnosis of Empty Shader Logs:**
    *   Despite valid shader source files (`nbody_vertex.glsl`, `nbody_fragment.glsl` confirmed via `cat`), shader linking continued to fail silently with empty error logs.
    *   This strongly indicates that the OpenGL context was not active or properly initialized when shader compilation/linking was attempted.
    *   The `load_shader_program` call in `GLCLBrowser.py`'s `__init__` method likely executes *before* `GLCLWidget`'s OpenGL context is ready (i.e., before `initializeGL()` is called).

3.  **Locating `GLCLWidget`:**
    *   Initial attempts to locate `GLCLWidget.py` or `GLCLGUI.py` were unsuccessful due to incorrect file paths.
    *   `list_dir` confirmed that `GLGUI.py` in `/home/prokophapala/git/FireCore/pyBall/GUI` is the most probable file containing the `GLCLWidget` class.

### Next Steps (Identified during this session):

*   **Defer Shader Compilation**: The primary task is to move the `ogl_system.load_shader_program` call from `GLCLBrowser.py`'s `__init__` to `GLCLWidget`'s `initializeGL()` method (or a method called from it) to ensure a valid OpenGL context is active during shader compilation.
