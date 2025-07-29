# Mission: A "Shadertoy" for Scientific Simulation

The primary goal of the GLCL project is to create a powerful, yet easy-to-use browser for developing and visualizing scientific simulations. Inspired by Shadertoy, this framework aims to eliminate boilerplate code related to GUI, OpenGL/OpenCL context management, and data handling.

This allows users—researchers, students, and hobbyists—to focus purely on their core algorithms, expressed through:
1.  **OpenCL kernels** for computation.
2.  **GLSL shaders** for visualization.
3.  A simple **Python configuration script** to wire everything together.

The browser dynamically loads and executes these user-provided components, providing an interactive environment for rapid prototyping and exploration of complex physical systems.

## Example file

There is example of user in `pyBall/GLCL/scripts/nbody.py` which is using 
* `pyBall/GLCL/scripts/nbody.cl` for OpenCL kernel
* `pyBall/GLCL/scripts/points.glslv` for vertex shader
* `pyBall/GLCL/scripts/points.glslf` for fragment shader
notice the directories, we need to search shaders and .cl source code in those directories.

this script is user script, which we load dynamically in the `GLCLBrowser.py`
We should not modify this script, nor the shaders or kernels it uses, nor any other user scripts!!!
this script should not have any dependencies and the format of config should stay as it is.

## Goal specification

make design of `GLCLBrowser.py` so that it can process (dynamically load) any config scripts of similar format. 
you are not supposed to create any new opencl kenrels or shares, they are responsibility of the user. 

You need to understand that our ultimate goal is not to implement nbody simulation per se. That is just an example of simulations which user may want to create. 
It is just one of may examples of many user scripts which user may eventually crate and want to run using our Browser.

We are making only the browser. The purpose of browser is to execute any script, shader and kernel which user provides. The purpose is to eliminate the boilerplate, allowing user focus on OpenCL and GLSL programing of physics simulations and not thinking about buffer management, GUI and this stuff. Basically we are making something like shaderoy in python augumented by OpenCL.

when you want to run test (which you should always do after code changes) do it like this:

`prokophapala@carbsisYoga:~/git/FireCore$ python -u -m pyBall.GLCL.GLCLBrowser`

---

# GLCL Simulation Framework Manifest (Updated)

## 1. Project Context and Goal

**Objective**: Refactor the GLCL molecular simulation framework to eliminate hardcoded buffer and uniform names. Instead of using JSON configuration files, we will use flexible Python scripts. We cannot assume any particular buffer or uniform names in the code, we can know what names of buffers and kernels and other object we get only after loading the script. This enables dynamic simulation parameters, buffer specifications, kernel pipelines, and rendering instructions, allowing for fully general loading and execution of arbitrary OpenCL/OpenGL simulations.

## Relevant files

* We can edit these files:
    * `GLCLBrowser.py` - main browser, our **primary file to edit**
    * `GLCLWidget.py` - main widget, integrates with OpenGL and OpenCL and PyQt5
    * `OGLsystem.py` - OpenGL system, handles OpenGL context and rendering
    * `OCLsystem.py` - OpenCL system, handles OpenCL context and kernels
* We should not much edit these dependecies:
    * `/home/prokophapala/git/FireCore/pyBall/GLCL/GLCLBrowser_old.py` - old browser, which works with `nbody.json` and `nbody.cl` and `points.glslv` and `monocolor.glslf`, but it is too hardcoded that is why we need to create new browser, but it works so use it as reference
    * `BaseGUI.py` - base GUI class, handles GUI creation and management
* Inspiration and examples :
    * Something similar for GLSL we have here (it also uses baking of render pipelines described using .json )
       * `GLSL_GUI.py` - GLSL GUI class, handles GLSL rendering and management
       * `GLSL_Simulation.py` - GLSL simulation class, handles GLSL simulation and rendering


## 2. Current State and Progress

### 2.1 Shift from JSON to Python Script Configuration

*   **Previous System**: Used JSON files (e.g., `nbody.json`). This was problematic as it did not allowed to execute any user-specific code - such as initialization of buffers.
*   **New Approach**: Adopted Python script configurations (e.g., `nbody.py`) that define a `config` dictionary (replacing the previous JSON file) and an `init()` function for buffer initialization. This provides maximum flexibility for dynamic simulation definitions and arbitrary initialization code.

### 2.2 Key Modifications in `GLCLBrowser.py`

*   **Dynamic Module Import**: Replaced JSON loading with dynamic Python module import using `importlib.util` in `load_and_apply_script` and `load_simulation_script`.
*   **Generalized `apply_simulation_config`**: This method now dynamically parses the `config` dictionary from the Python script to:
    *   Populate GUI parameters (uniforms) dynamically from the `parameters` section.
    *   Load OpenCL kernel sources and OpenGL shaders based on paths in the `config`.
    *   Create OpenCL buffers dynamically with names and dimensions specified in the `config` and support for initial data from `init()`.
    *   Set kernel parameters dynamically from GUI-controlled parameters.
    *   Pass the full configuration to the `GLCLWidget` for rendering setup.
*   **Dynamic `update_simulation`**: Implemented logic to iterate through the `kernel_pipeline` defined in `self.config` and execute OpenCL kernels. It also includes a preliminary mechanism to update OpenGL buffers (currently focused on 'positions').
*   **Dynamic `update_sim_uniforms`**: Generalized this method to dynamically update OpenCL kernel parameters based on the `self.config` dictionary and current GUI widget values. It includes special handling for `particle_count` to re-initialize buffers if the count changes.
*   **Main Execution Block**: Updated to load `nbody.py` as an example of simulation script. But other scripts can be loaded as well defining completely different simulations.

### 2.3 Dependencies and APIs

*   **PyQt5**: For GUI and OpenGL widget integration.
*   **PyOpenCL**: For OpenCL context, kernel, and buffer management.
*   **NumPy**: For numerical data and buffer initialization.
*   **`importlib.util`**: For dynamic Python module loading.

### 3.2 Generalization and Refinement

1.  **Generalize OpenGL Buffer Updates**: The `update_simulation` method currently only updates the 'positions' buffer for OpenGL rendering. This needs to be generalized to handle any relevant buffers specified in the `config`'s `render_pipeline`. We cannot rely on any hard-coded names in our program, everything needs to be dynamic responsive to the config from the loaded script.
2.  **Dynamic OpenGL Uniforms**: Implement dynamic uniform updates for OpenGL shaders in `update_sim_uniforms` based on GUI parameter changes. Currently, it primarily focuses on OpenCL kernel arguments.
3.  **Complete `GLCLWidget` Integration**: Ensure that the `GLCLWidget` correctly interprets and utilizes the full `config` passed to it for rendering, including draw calls, buffer bindings, and uniform settings.
4.  **Error Handling and Validation**: Add more robust error handling and input validation for configuration files and runtime operations.
5.  **Documentation**: Document the expected structure and requirements for Python simulation configuration scripts to guide future development and extension.

## 4. Proposed New Architecture: The Baking Approach

To address the issues of inefficiency, hardcoded references, and lack of generality, a new "baking" architecture will be implemented. This approach will process the simulation configuration once at the start to create highly optimized, executable pipelines, eliminating runtime checks and string comparisons.

### 4.1 Core Concept: Baking

Instead of interpreting the `config` dictionary in real-time every frame, we will process it once in `apply_simulation_config` (and within `GLCLWidget`) to build two "baked" pipelines:

*   **Baked Kernel Pipeline**: An ordered list of callable Python objects, where each object represents a fully prepared OpenCL kernel launch.
*   **Baked Render Pipeline**: An ordered list of callable objects for the `GLCLWidget`, where each object represents a prepared OpenGL draw call (a render pass).

### 4.2 Baking the Kernel Pipeline (in `GLCLBrowser.py`)

A new mechanism will be introduced to create `BakedKernelCall` objects. Each `BakedKernelCall` will encapsulate:

*   The actual `pyopencl.Kernel` object (fetched once from `OCLSystem`).
*   Pre-resolved global and local work sizes.
*   A prepared list of arguments, with kernel arguments (uniforms, buffers) pre-resolved (i.e. we do not store names, but handles to OpenCL objects).

The `update_simulation` function will be simplified to iterate through this `self.baked_kernel_pipeline` and execute each `BakedKernelCall` directly.

### 4.3 Baking the Render Pipeline (in `GLCLGUI.py`)

The `GLCLWidget` will receive the full `config` and, during its initialization (e.g., in `set_config` or `initializeGL`), will bake a render pipeline. This involves creating `BakedRenderPass` objects, each containing:

*   The compiled OpenGL shader program.
*   The Vertex Array Object (VAO) to use.
*   The drawing mode (e.g., `GL_POINTS`).
*   uniforms for shaders does not need to be baked, because they can be set to each shader at the start, no need to do it every frame, only when uniforms are updated.

The `paintGL` method will then simply iterate through `self.baked_render_pipeline` and execute each `BakedRenderPass`.

### 4.4 Efficient Parameter Updates

The problematic `update_sim_uniforms` method will be removed. Instead, during the initial baking process, a `param_update_map` will be constructed. This map will link GUI parameter names directly to the specific attributes or closures within the `BakedKernelCall` and `BakedRenderPass` objects that need to be updated.

**Current Issue**: Parameter values from config are being reset to 0.0 during initialization instead of using the correct values from the simulation script.

**Root Cause**: The parameter dictionary is being modified somewhere between config loading and kernel baking. The expected flow should be:
1. Config dictionary contains original values (e.g., dt=0.001)
2. Controls are initialized FROM the config dictionary
3. Only when user changes controls, the config dictionary should update
4. Kernel baking should always use the original config values directly.

**Fix Required**: Ensure create_parameter_controls() only reads from config, never modifies it, and that bake_kernels() uses the original config values directly.

### 4.5 Benefits

*   **Generality**: No hardcoded names of buffers or uniforms, kernels or shaders assumed in the code; all names are derived from the dynamic `config`.
*   **Efficiency**: All expensive lookups and processing occur once at startup, not every frame.
*   **Modularity**: Clear separation between configuration parsing/baking and runtime execution.
*   **Readability**: The main simulation loop becomes much cleaner and easier to understand.

## Reference Hard-coded implementation

As a reference please look into `GLCL/GLCLBrowser_old.py` 

## Debugging Strategy

Efficient debbuging is our priority of the code desing. We are not creating app for layman users, but for scientists and engineers who need insight into the code and the simulations. We do not want to hide anything. Therefoore we avoid any try-excpet hading which can hide errors and make debugging more difficult. We should never skip silently unexpected events or errors. We should always prefer to crash loudly and provide as much information as possible to the user, especially the stack trace. This also makes our debugind process eaiser. Instead of if (file-not-exist): return, we should do raise FileNotFoundError. We should never catch exceptions and print it without crashing, or even handle it silently.

### 4.6 Runtime Debug Flags (added 2025-07-28)

To streamline troubleshooting of user scripts the browser now supports runtime debug switches:

* **bDebugCL** (default False)
  * Kernels run as usual.
  * After each kernel in the baked pipeline completes, the first few elements of every OpenCL buffer are downloaded and printed.
  * The main loop automatically stops after **nDebugFrames** iterations to keep logs manageable.

* **bDebugGL** (default False)
  * All baked kernels are skipped (buffers remain in their initial state).
  * OpenGL rendering continues, letting you verify that initial data are displayed correctly (e.g. check that `positions` appear as points for `nbody.py`).

These flags are optional parameters of `GLCLBrowser(python_script_path, bDebugCL=False, bDebugGL=False, nDebugFrames=5)`.

At the end of the setup phase the browser prints a concise summary of:

* Allocated buffer names and shapes.
* Baked kernels with their resolved `global_size` and `local_size`.

This confirms that configuration parsing and baking worked before the simulation loop starts.

## 5. Dynamic Buffer Management System

### 5.1 Problem Statement
Current implementation has hardcoded buffer names like 'positions' in GLCLBrowser.py and GLCLGUI.py. This violates the principle of a generic browser that should work with any user-defined buffer names.

### 5.2 Solution Architecture
Instead of hardcoding buffer names, we need:

1. **Dynamic Buffer Registry**: A system that creates GLobject instances for every buffer mentioned in user scripts
2. **Name-based Lookup**: Use dictionary mapping buffer names to GLobject instances
3. **Render Pipeline Binding**: Connect shader inputs to buffers by name as specified in render_pipeline configuration
4. **Generic Data Passing**: Remove all hardcoded references to specific buffer names

### 5.3 Implementation Details - Baking Strategy

#### Baking Phase (Setup)
- **Parse render pipeline** from user script to identify all render passes
- **Create GLobject instances** for each render pass, pre-configured with:
  - Shader program
  - Vertex buffer data
  - Element count
  - Render mode (GL_POINTS, GL_TRIANGLES, etc.)
- **Store in list**: `render_objects = [GLobject, GLobject, ...]`

#### Runtime Phase
- **Simple iteration**: For each GLobject in render_objects, call `obj.draw()`
- **No runtime lookups**: All buffer binding and configuration is pre-baked
- **Render pipeline becomes**: A simple list of GLobjects ready for drawing

#### GLobject Pre-configuration
Each GLobject is created with all necessary state:
- Shader program ID
- VAO/VBO configuration
- Vertex data uploaded
- Element count set
- Render mode specified

### 5.4 Baking Process
1. **Parse user script** to extract render pipeline configuration
2. **Create GLobject** for each render pass
3. **Pre-configure** each GLobject with all necessary OpenGL state
4. **Store in list** for runtime iteration
5. **Runtime**: Simply call `draw()` on each pre-baked GLobject

### 5.4 Required Changes

#### GLCLBrowser.py
- Remove hardcoded 'positions' buffer handling
- Implement dynamic buffer registry creation
- Update render pipeline to use name-based buffer binding

#### GLCLGUI.py
- Remove hardcoded particle_count and positions references
- Implement generic buffer data passing system
- Use GLobject instances dynamically created based on user script

### 5.5 Benefits
- Complete abstraction from specific buffer names
- Support for any user-defined buffer configuration
- True generic browser capability
- Eliminates need for browser updates when user scripts use different buffer names

### 5.6 Kernel Baking System

**Problem**: Current kernel execution has runtime overhead - parameters are looked up and type-converted every frame.

**Solution**: Implement true kernel baking that pre-computes everything at setup time.

**Design Pattern**:
```python
# Kernel description format
"kernel_name" : ( (local_size,), "global_size_expr", [buffer_args], [scalar_args] )

# Baked result: (kernel_func, global_size, local_size, args_tuple)
# Execution: kernel_func(queue, global_size, local_size, *args_tuple)
```

**Implementation**:
1. `bake_kernel()` - creates pre-computed (kernel, global_size, local_size, args) tuple
2. `execute_baked()` - simple unpacking and execution
3. `setup_kernel_params()` - initializes all parameters immediately after baking

**Benefits**:
- Zero runtime overhead
- Simple execution: `kernel(*args)`
- Pre-computed global size resolution
- Pre-converted parameter types
- Immediate parameter initialization

### 5.7 Critical: OpenGL Initialization Timing

**IMPORTANT**: When working with OpenGL, certain operations must be deferred until after the OpenGL context is fully initialized:

- **Shader compilation** (glCreateShader, glCompileShader)
- **Buffer creation** (glGenBuffers, glBufferData)
- **Texture creation** (glGenTextures)
- **Framebuffer operations** (glGenFramebuffers)

**Symptoms of premature initialization**:
- glCreateShader returns 0
- OpenGL calls silently fail
- Random crashes or black screens

**Solution**:
- Move all OpenGL resource creation to GLCLWidget.initializeGL()
- Use explicit callbacks between GLCLBrowser and GLCLWidget
- Add debug prints to verify OpenGL context availability

**Debugging Tips**:
1. Check OpenGL version info at startup
2. Verify all OpenGL calls succeed
3. Never assume OpenGL is ready before initializeGL()
4. Add explicit error checking for all OpenGL calls

This issue has caused multiple hours of debugging - always verify OpenGL context is ready before making OpenGL calls!