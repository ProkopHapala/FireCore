This is good start but I ralized I already have several classes which I would like to reuse (use a base for our new viewer) so we do not reinvent the wheel and don't start from scratch

Please look on folloding classes/ files

We have:
- `GLGUI.py` - Genergal Base class for PyQt GUI in combination with OpenGL
- `xyz_view_new.py` - this is example using `GLGUI.py`, we applied it to visualization of electron-forcefield  (`EFF`), which has just balls (no bonds), but what we want to do will be probably rather similar.
- `AtomicSystem.py` - this is class for opreation with molecular geometry and topology
   - it uses `atomicUtils.py` as backend for operation for atomic geometry, topology, file I/O, etc. but we do not need to know details, the functions and variables exposed by AtomicSystem should be enough for us


# OpenGL Viewer Documentation

This document provides an overview of the OpenGL viewer components within the `FireCore` project, focusing on `pyBall/GUI/GLGUI.py`, `tests/tEFF/xyz_view_new.py`, and their associated GLSL shaders. This documentation aims to facilitate both the use and further development of existing applications like `xyz_view_new.py`, as well as the creation of new OpenGL applications based on the `GLGUI.py` platform.

## 1. Core OpenGL GUI Framework (`pyBall/GUI/GLGUI.py`)

`GLGUI.py` serves as the foundational layer for OpenGL-based graphical user interfaces. It provides a `QOpenGLWidget` subclass (`GLGUI`) that encapsulates essential OpenGL setup, rendering loop management, and user interaction for camera control.

### 1.1. `GLGUI` Class

- **Purpose:** A `QOpenGLWidget` subclass providing the core OpenGL context, rendering loop, and camera interaction.
- **Key Methods & Attributes:**
    - `initializeGL()`: Initializes OpenGL functions, compiles shaders (in derived classes), and sets up default rendering states.
    - `paintGL()`: The main rendering function, called repeatedly to draw the scene. It sets up camera uniforms and calls `draw_scene()`.
    - `resizeGL(width, height)`: Handles viewport resizing and updates the projection matrix.
    - `mousePressEvent(event)`, `mouseMoveEvent(event)`, `wheelEvent(event)`: Implement trackball-like camera rotation and zoom functionality.
    - `_project_to_sphere(x, y)`: Helper for trackball rotation, projecting screen coordinates onto a virtual sphere.
    - `_rotation_to_gl_matrix(rotation_obj)`: Converts `scipy.spatial.transform.Rotation` objects to OpenGL 4x4 matrices.
    - `draw_scene()`: An abstract method that *must* be overridden by derived classes to implement specific drawing logic. The shader program is already bound when this is called.
    - `projection_matrix`, `view_matrix`, `model_matrix`: `QMatrix4x4` objects managing the camera and global scene transformations.
    - `orientation`: `scipy.spatial.transform.Rotation` object storing the current camera orientation.

### 1.2. `Mesh` Class

- **Purpose:** Manages vertex data (positions, normals) and optional index data for 3D models. It handles the creation and binding of OpenGL Vertex Array Objects (VAOs), Vertex Buffer Objects (VBOs), and Element Buffer Objects (EBOs).
- **Key Methods:**
    - `__init__(vertices, normals=None, indices=None)`: Constructor, takes numpy arrays for vertex data.
    - `setup_buffers()`: Generates and populates VAO, VBOs, and EBO (if indices are provided).
    - `draw(mode=GL_TRIANGLES)`: Renders the mesh using `glDrawElements` (if indexed) or `glDrawArrays`.
    - `cleanup()`: Releases OpenGL resources (buffers, VAO).
    - `bind()`, `unbind()`: Binds/unbinds the mesh's VAO.
    - `draw_instanced(num_instances, mode=GL_TRIANGLES)`: Renders multiple instances of the mesh.

### 1.3. `InstancedData` Class

- **Purpose:** Facilitates efficient instanced rendering by managing per-instance data (e.g., positions, colors, scales) and associating it with a base `Mesh`.
- **Key Methods:**
    - `__init__(base_attrib_location)`: Initializes with a starting attribute location for instance data.
    - `associate_mesh(mesh_object)`: Links the instance data manager to a `Mesh` object, setting up its VAO for instanced drawing.
    - `setup_instance_vbos(attribs)`: Creates VBOs for instance attributes. `attribs` is a list of tuples `(name, components, type)`.
    - `update(buffs, n=-1)`: Updates the data in the instance VBOs. `buffs` is a dictionary mapping attribute names to numpy arrays.
    - `draw(mode=GL_TRIANGLES)`: Renders the associated mesh using `glDrawElementsInstanced` or `glDrawArraysInstanced`.
    - `cleanup()`: Releases OpenGL resources.

### 1.4. Utility Functions

- `set_ogl_blend_mode(mode, depth_test=True)`: Configures OpenGL blending for transparency effects based on predefined modes.
- `disable_blend(depth_test=True)`: Disables blending and re-enables depth testing.
- `setup_vbo(vertices, array_indx, components=3, usage=GL_STATIC_DRAW)`: Helper to create and setup a single VBO.
- `setup_ebo(indices)`: Helper to create and setup an EBO.
- `create_sphere_mesh(radius=1.0, segments=16, rings=16)`: Generates vertex and normal data for a sphere mesh.
- `alpha_blend_modes`: A dictionary defining common OpenGL blending configurations.

### 1.5. `AppWindow` Class

- **Purpose:** A `QMainWindow` subclass providing a standard application window for OpenGL applications.
- **Key Method:**
    - `launch(*args, **kwargs)` (classmethod): A convenient static method to create and run a Qt application, instantiating the `AppWindow` (or derived class) and showing it.

## 2. Molecular Viewer Application (`tests/tEFF/xyz_view_new.py`)

`xyz_view_new.py` is a concrete example of an OpenGL application built upon the `GLGUI.py` framework. It specializes in visualizing molecular trajectories and electron force field data.

### 2.1. `MolViewerWidget` Class

- **Purpose:** Extends `BaseGLWidget` (from `GLGUI.py`) to provide specific rendering logic for molecular data.
- **Inheritance:** Inherits from `BaseGLWidget` (likely `GLGUI.GLGUI` or a similar base class).
- **Data Management:** Stores trajectory data (`self.trj`), current frame, and pre-processed atom/electron data (`self.frames_atoms`, `self.frames_elecs`).
- **Shader Loading & Management:**
    - Loads and compiles two primary fragment shaders: `sphere.glslf` (for ray-traced spheres) and `sphere_max.glslf` (for volumetric spheres).
    - Defines `self.render_modes` dictionary, mapping user-friendly names (e.g., "Shiny", "Volumetric") to specific shader programs, blending modes, and depth test settings.
- **Instanced Rendering Setup:**
    - Initializes two `InstancedData` objects: `self.atom_instances` and `self.elec_instances`.
    - Associates a shared sphere mesh (`self.default_sphere_mesh`) with both instance managers.
    - Configures instance VBOs for `positions`, `radii`, and `colors` for both atoms and electrons.
- **Key Methods:**
    - `initializeGL()`: Overrides base class to load specific shaders, set up `InstancedData`, and define `render_modes`.
    - `cleanupGL()`: Releases OpenGL resources specific to the viewer.
    - `precalculate_frames_data()`: Processes raw trajectory data into optimized numpy arrays for atoms and electrons, including calculating colors and radii.
    - `update_instance_data()`: Updates the VBOs of `self.atom_instances` and `self.elec_instances` with data for the current frame, including applying opacity to electron colors.
    - `draw_scene()`: Overrides base class to render atoms (opaque, using `sphere.glslf`) and then electrons (potentially transparent, using selected shader and blend mode from `render_modes`). It manages shader switching and blending for these two passes.
    - `use_shader(shader_prog_id)`: Sets the currently active shader program.
    - `set_frame(frame_idx)`: Updates the current trajectory frame and triggers a re-render.
    - `set_opacity(opacity_percent)`: Adjusts the opacity of electron instances and triggers a re-render.
    - `set_render_mode(mode_key)`: Changes the rendering mode for electrons based on the `render_modes` dictionary.

### 2.2. `MolViewer` Class

- **Purpose:** Extends `AppWindow` (from `GLGUI.py`) to create the main application window for the molecular viewer, integrating the `MolViewerWidget` with UI controls.
- **Key Methods:**
    - `__init__(trj)`: Sets up the main window, creates the `MolViewerWidget`, and adds UI elements like frame and opacity sliders, and a render mode combo box.
    - `on_render_mode_changed(text)`: Callback for the render mode combo box, updating the `MolViewerWidget`'s render mode.

## 3. GLSL Shaders

The viewer utilizes several GLSL shaders for rendering, demonstrating different rendering techniques for spherical objects.

### 3.1. `pyBall/GUI/shaders/instances.glslv` (Vertex Shader)

- **Purpose:** Transforms the vertices of a base sphere mesh (acting as a bounding box) into clip space, preparing per-instance data for the fragment shader.
- **Inputs:**
    - `aPos` (location 0): Base mesh vertex position.
    - `aNormal` (location 1): Base mesh vertex normal.
    - `instancePosition_model` (location 2): Center of the sphere instance in model space.
    - `instanceActualSphereRadius` (location 3): Actual radius of the sphere for ray-tracing/volumetric rendering.
    - `instanceColor` (location 4): RGBA color of the instance.
    - `uniform mat4 projection`, `view`, `model`: Standard transformation matrices.
- **Outputs (to Fragment Shader):**
    - `fpos_world`: World-space position of the fragment on the bounding mesh.
    - `sphere_obj_world`: `vec4(sphere_center_world.xyz, sphere_radius_world.w)` - the actual sphere data for intersection tests.
    - `fColor`: Pass-through instance color.
- **Key Operations:** Calculates world-space sphere center and radius, scales base mesh vertices, and computes `gl_Position` for the bounding box.

### 3.2. `pyBall/GUI/shaders/sphere_max.glslf` (Fragment Shader - Volumetric)

- **Purpose:** Renders spheres with a soft, volumetric appearance, suitable for visualizing electron density or other continuous distributions. It uses a ray-marching like approach to determine fragment opacity.
- **Inputs:**
    - `sphere_obj_world`, `fpos_world`, `fColor`: From vertex shader.
    - `uniform vec3 viewPos`: Camera position.
    - `uniform mat4 projection`, `view`: For `gl_FragDepth` calculation.
- **Output:** `FragColor` (RGBA).
- **Key Operations:**
    - Defines a ray from `viewPos` to `fpos_world`.
    - Uses `rayPointDist` to find the closest approach of the ray to the sphere's center.
    - If the ray intersects the sphere, calculates a `density` based on the normalized distance from the sphere's center (exponential fall-off).
    - Sets fragment `alpha` as `density * fColor.a`.
    - `discard`s fragments that miss the sphere.
    - Explicitly sets `gl_FragDepth` for correct transparency rendering.

### 3.3. `pyBall/GUI/shaders/sphere.glslf` (Fragment Shader - Ray-traced Sphere)

- **Purpose:** Renders spheres with a solid, shiny appearance using ray-sphere intersection and Phong lighting model (ambient, diffuse, specular).
- **Inputs:**
    - `sphere_obj_world`, `fpos_world`, `fColor`: From vertex shader.
    - `uniform vec3 viewPos`, `lightPos`, `lightColor`.
    - `uniform mat4 projection`, `view`: For `gl_FragDepth` calculation.
- **Output:** `FragColor` (RGBA).
- **Key Operations:**
    - Defines a ray from `viewPos` to `fpos_world`.
    - Uses `raySphere` to find the intersection point on the sphere.
    - `discard`s fragments that miss the sphere or are behind the camera.
    - Calculates surface normal (`sphereNormal`) at the intersection point.
    - Applies Phong lighting: ambient, diffuse (based on `dot(N, L)`), and specular (based on `pow(dot(V, R), shininess)`).
    - Sets `FragColor` with the calculated lighting and `fColor.a`.
    - Explicitly sets `gl_FragDepth` for correct depth buffering.

## 4. Building New OpenGL Applications

To build new OpenGL applications or editors based on the `GLGUI.py` platform, follow these steps:

1.  **Inherit from `GLGUI.GLGUI` (or `BaseGLWidget`):** Create a new class that inherits from `GLGUI.GLGUI` (or `BaseGLWidget` if `xyz_view_new.py` is used as a direct template). This class will be your custom OpenGL rendering widget.
2.  **Override `initializeGL()`:** In your derived class, load and compile your specific GLSL shaders. Set up any `Mesh` or `InstancedData` objects you need for your scene.
3.  **Override `draw_scene()`:** Implement your custom rendering logic here. This is where you will bind your shaders, set uniforms, and call `draw()` or `draw_instanced()` on your `Mesh` or `InstancedData` objects.
4.  **Create an `AppWindow` (or inherit from it):** For a complete application, create a `QMainWindow` that hosts your custom `GLGUI` widget. You can inherit from `GLGUI.AppWindow` for convenience and add your UI controls (sliders, buttons, etc.) to interact with your OpenGL scene.
5.  **Use `AppWindow.launch()`:** Use the `AppWindow.launch()` class method to easily start your application.

### Example Structure for a New Application:

```python
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from pyBall.GUI.GLGUI import GLGUI, AppWindow, Mesh, InstancedData
from OpenGL.GL.shaders import compileProgram, compileShader
from OpenGL.GL import *

# Define your custom GLSL shaders (or load from files)
VERT_SHADER_SOURCE = """#version 330 core\n..."""
FRAG_SHADER_SOURCE = """#version 330 core\n..."""

class MyCustomGLWidget(GLGUI):
    def initializeGL(self):
        super().initializeGL() # Call base class init
        self.shader_program = self.compile_shader_program(VERT_SHADER_SOURCE, FRAG_SHADER_SOURCE)

        # Example: Setup a simple mesh
        vertices = np.array([...], dtype=np.float32)
        self.my_mesh = Mesh(vertices)
        self.my_mesh.setup_buffers()

    def draw_scene(self):
        glUseProgram(self.shader_program)
        # Set uniforms specific to your shader
        # glUniformMatrix4fv(glGetUniformLocation(self.shader_program, "model"), 1, GL_FALSE, self.model_matrix.data())
        # ... other uniforms ...

        self.my_mesh.draw() # Draw your mesh

    def cleanupGL(self):
        super().cleanupGL() # Call base class cleanup
        self.my_mesh.cleanup()

class MyCustomAppWindow(AppWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("My Custom OpenGL App")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.gl_widget = MyCustomGLWidget()
        layout.addWidget(self.gl_widget)

        # Add custom UI controls here if needed
        # slider = QSlider(Qt.Horizontal)
        # slider.valueChanged.connect(self.gl_widget.some_method_to_control_rendering)
        # layout.addWidget(slider)

        self.show()

if __name__ == '__main__':
    MyCustomAppWindow.launch()
```

This documentation provides a solid foundation for understanding and extending the `FireCore` OpenGL viewer components. Remember to consult the actual source code for the most up-to-date and detailed implementation specifics.



## GLGUI.py


```python
import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse

from PyQt5.QtWidgets import ( QOpenGLWidget)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from PyQt5.QtWidgets import (QMainWindow)
#from PyQt5.QtOpenGLWidgets import QOpenGLWidget
from PyQt5.QtGui import QSurfaceFormat
from PyQt5.QtWidgets import QApplication

from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

# It's good practice to keep shaders in separate files or as multi-line strings
# For simplicity here, they are embedded as strings.

alpha_blend_modes={
    "standard"    :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    "standard2"    :(GL_FUNC_ADD, GL_ONE, GL_ONE_MINUS_SRC_ALPHA),
    #"standard2"   :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE),
    #"standard3"   :(GL_FUNC_ADD, GL_ONE,       GL_ONE),
    #"minimum1"    :(GL_MIN,      GL_ONE,       GL_ONE),
    #"minimum2"    :(GL_MIN,      GL_SRC_ALPHA, GL_ONE),
    #"minimum3"    :(GL_MIN,      GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    #"minimum4"    :(GL_MIN,GL_ONE, GL_ONE,   GL_FUNC_ADD,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA ),
    #"additive"   :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE),
    "subtractive" :(GL_FUNC_REVERSE_SUBTRACT, GL_SRC_ALPHA, GL_ONE),
    #"subtractive2" :(GL_FUNC_REVERSE_SUBTRACT, GL_ONE, GL_ONE),
    #"subtractive2" :(GL_FUNC_REVERSE_SUBTRACT, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    #"subtractive":(GL_FUNC_SUBTRACT, GL_SRC_ALPHA, GL_ONE),
    #"minimum"    :(GL_FUNC_MIN, GL_ONE, GL_ONE),
    #"maximum"    :(GL_MAX, GL_ONE, GL_ONE)
}

def set_ogl_blend_mode(mode, depth_test=True):
    #print("set_ogl_blend_mode", mode)
    glEnable(GL_BLEND)
    if len(mode) == 3:
        glBlendEquation(mode[0])
        glBlendFunc(mode[1], mode[2])
    else:
        glBlendEquationSeparate(mode[0], mode[3])
        glBlendFuncSeparate    (mode[1], mode[2], mode[4], mode[5])
    #glDepthMask(GL_FALSE)
    #glDepthTest(GL_FALSE)
    if depth_test:
        glDisable(GL_DEPTH_TEST)

def disable_blend( depth_test=True):
    glDisable(GL_BLEND)
    #glDepthMask(GL_TRUE)
    #glDepthTest(GL_TRUE)
    if depth_test:
        glEnable(GL_DEPTH_TEST)

def setup_vbo(vertices, array_indx, components=3, usage=GL_STATIC_DRAW):
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, usage)
    glVertexAttribPointer(array_indx, components, GL_FLOAT, GL_FALSE, 0, None) # Location 0 for positions
    glEnableVertexAttribArray(array_indx)
    return vbo

def setup_ebo(indices):
    ebo = glGenBuffers(1)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_STATIC_DRAW)
    return ebo

class Mesh:
    def __init__(self, vertices, normals=None, indices=None):
        self.vertices = vertices
        self.normals = normals
        self.indices = indices

        self.vao = None
        self.vbo_vertices = None
        self.vbo_normals = None
        self.ebo = None

        self.has_indices = indices is not None
        self.vertex_count = len(vertices) // 3
        self.index_count = len(indices) if self.has_indices else 0

    def setup_buffers(self):
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vbo_vertices = setup_vbo(self.vertices, 0)
        if self.normals is not None:
            self.vbo_normals = setup_vbo(self.normals, 1)
        if self.has_indices:
            self.ebo = setup_ebo(self.indices)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0) # Unbind VAO

    def draw(self, mode=GL_TRIANGLES):
        glBindVertexArray(self.vao)
        if self.has_indices:
            glDrawElements(mode, self.index_count, GL_UNSIGNED_INT, None)
        else:
            glDrawArrays(mode, 0, self.vertex_count)
        glBindVertexArray(0)

    def cleanup(self):
        if self.vbo_vertices: glDeleteBuffers(1, [self.vbo_vertices])
        if self.vbo_normals:  glDeleteBuffers(1, [self.vbo_normals])
        if self.has_indices and self.ebo: glDeleteBuffers(1, [self.ebo])
        if self.vao:          glDeleteVertexArrays(1, [self.vao])

    def bind(self):
        glBindVertexArray(self.vao)

    def unbind(self):
        glBindVertexArray(0)

    def draw_instanced(self, num_instances, mode=GL_TRIANGLES):
        # Assumes the correct VAO (e.g., from InstancedData) is already bound by the caller.
        if self.has_indices:
            glDrawElementsInstanced(mode, self.index_count, GL_UNSIGNED_INT, None, num_instances)
        else:
            glDrawArraysInstanced(mode, 0, self.vertex_count, num_instances)


class InstanceData:
    def __init__(self, base_attrib_location):
        self.vbos = {}
        self.base_attrib_location = base_attrib_location # Starting location for instance attributes
        self.associated_mesh = None
        self.num_instances = 0 # Now stored here
        self.vao = None # Each InstancedData will have its own VAO
        self._instance_attrib_configs = [] # To store VBO configs

def make_vbo( num_components, base_attrib_location, attrib_idx_offset, usage=GL_DYNAMIC_DRAW):
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, 0, None, usage) # Initial empty allocation
    attrib_loc = base_attrib_location + attrib_idx_offset
    glVertexAttribPointer(attrib_loc, num_components, GL_FLOAT, GL_FALSE, 0, None)
    glEnableVertexAttribArray(attrib_loc)
    glVertexAttribDivisor(attrib_loc, 1) # Per-instance
    return vbo

class InstancedData:
    def __init__(self, base_attrib_location):
        self.vbos = []
        self.vbos_inds = {}
        self.base_attrib_location = base_attrib_location # Starting location for instance attributes
        self.associated_mesh = None
        self.num_instances = 0 # Now stored here
        self.vao = None # Each InstancedData will have its own VAO
        self._instance_attrib_configs = [] # To store VBO configs

    def cleanup(self):
        for vbo in self.vbos.values():
            glDeleteBuffers(1, [vbo])
        self.vbos.clear()
        if self.vao:
            glDeleteVertexArrays(1, [self.vao])
            self.vao = None

    def associate_mesh(self, mesh_object):
        self.associated_mesh = mesh_object
        if self.associated_mesh is None:
            if self.vao: # Cleanup old VAO if mesh is disassociated
                glDeleteVertexArrays(1, [self.vao])
                self.vao = None
            return

        if self.vao is None:
            self.vao = glGenVertexArrays(1)

        glBindVertexArray(self.vao)
        # Bind mesh's vertex buffer and set attribute pointer for location 0 (aPos)
        glBindBuffer(GL_ARRAY_BUFFER, self.associated_mesh.vbo_vertices)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)

        # Bind mesh's normal buffer and set attribute pointer for location 1 (aNormal)
        if self.associated_mesh.vbo_normals is not None:
            glBindBuffer(GL_ARRAY_BUFFER, self.associated_mesh.vbo_normals)
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, None)
            glEnableVertexAttribArray(1)
        
        # If mesh has EBO, bind it to this VAO as well
        if self.associated_mesh.ebo is not None:
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.associated_mesh.ebo)

        glBindVertexArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0) # Unbind GL_ARRAY_BUFFER
        if self.associated_mesh.ebo is not None:
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) # Unbind GL_ELEMENT_ARRAY_BUFFER

    def setup_instance_vbos(self, attribs):
        if self.vao is None:
            print("Error: InstancedData VAO not initialized. Call associate_mesh first.")
            return
        glBindVertexArray(self.vao)
        self._instance_attrib_configs = [] # Clear previous configs
        for i, (name, attrib_idx_offset, num_components) in enumerate(attribs):
            vbo = make_vbo(num_components, self.base_attrib_location, attrib_idx_offset)
            self.vbos.append(vbo)
            self.vbos_inds[name]=i  # this should be faster than using name as key
        glBindVertexArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        return self.get_vbo_inds( [name for name, _, _ in attribs] )

    def get_vbo_inds(self, names):
        return [self.vbos_inds[name] for name in names]

    def update(self, buffs, n=-1):
        if n<0: n=len(buffs["positions"])
        self.num_instances = n
        for name, buff_data in buffs.items():
            i = self.vbos_inds[name]
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[i])
            glBufferData(GL_ARRAY_BUFFER, buff_data.nbytes, buff_data, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def update_list(self, buffs, n=-1):  # this should be faster than using name as key
        if n<0: n=len(buffs[0])
        self.num_instances = n
        for i, buff_data in enumerate(buffs):
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[i])
            glBufferData(GL_ARRAY_BUFFER, buff_data.nbytes, buff_data, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def draw(self, mode=GL_TRIANGLES):
        if self.vao is None or self.associated_mesh is None or self.num_instances == 0:
            return

        glBindVertexArray(self.vao)
        # The VAO is already configured with base mesh EBO if it exists (done in associate_mesh)
        self.associated_mesh.draw_instanced(self.num_instances, mode)
        glBindVertexArray(0)

# ================= Free Utility Functions =================

def create_sphere_mesh(radius=1.0, segments=16, rings=16):
    """Generates vertices and normals for a sphere."""
    vertices = []
    normals = []
    indices = []

    for ring_num in range(rings + 1):
        theta = ring_num * np.pi / rings
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        for seg_num in range(segments + 1):
            phi = seg_num * 2 * np.pi / segments
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)

            x = radius * cos_phi * sin_theta
            y = radius * sin_phi * sin_theta
            z = radius * cos_theta
            
            nx, ny, nz = x/radius, y/radius, z/radius # Normals are just normalized positions for a sphere centered at origin

            vertices.extend([x, y, z])
            normals.extend([nx,ny,nz])

    for ring_num in range(rings):
        for seg_num in range(segments):
            first  = (ring_num * (segments + 1)) + seg_num
            second = first + segments + 1
            
            indices.extend([first, second, first + 1])
            indices.extend([second, second + 1, first + 1])
            
    return np.array(vertices, dtype=np.float32), np.array(normals, dtype=np.float32), np.array(indices, dtype=np.uint32)

def append_normalized( v_list, radius, vertices_list, normals_list ):
    for v in v_list:
        v_norm = v / np.linalg.norm(v)
        vertices_list.extend(v_norm * radius)
        normals_list .extend(v_norm)

def octahedron_sphere_face(c1, c2, c3, vertices_list, normals_list, nsub=2,  radius=1.0):
    """
    Generates vertices and normals for one face of an octahedron, subdivided,
    for non-indexed drawing. Appends directly to vertices_list and normals_list.
    """
    # Interpolation vectors along the edges of the base triangle face
    da = (c2 - c1) / nsub
    db = (c3 - c1) / nsub
    for i in range(nsub):      # Iterate along the c1-c2 edge direction
        for j in range(nsub - i):  # Iterate along the c1-c3 edge direction
            # Define the 3 vertices of the "bottom-left" triangle in the current cell
            p00 = c1 + da*i     + db*j
            p10 = c1 + da*(i+1) + db*j
            p01 = c1 + da*i     + db*(j+1)
            append_normalized([p00, p10, p01], radius, vertices_list, normals_list)
            # Define the 3 vertices of the "top-right" triangle in the current cell (if it exists)
            if j < nsub - 1 - i: # Check if there's a second triangle in this parallelogram
                p11 = c1 + da*(i+1) + db*(j+1)
                append_normalized([p10, p11, p01], radius, vertices_list, normals_list)

def octahedron_sphere_mesh(radius=1.0, nsub=2):
    """
    Generates vertices and normals for a sphere by subdividing an octahedron,
    vectorized with NumPy for non-indexed GL_TRIANGLES drawing.
    'subdivisions' is recursion depth: 0 for base octahedron (8 tris), 1 for 32 tris, etc.
    Returns flat vertices, flat normals, and None for indices.
    """
    pv_np = np.array([
        [ 1.0,  0.0,  0.0], [-1.0,  0.0,  0.0], [ 0.0,  1.0,  0.0], 
        [ 0.0, -1.0,  0.0], [ 0.0,  0.0,  1.0], [ 0.0,  0.0, -1.0]
    ], dtype=np.float32)
    octa_face_indices = np.array([ # CCW from outside
        [4, 1, 3], [4, 3, 0], [4, 0, 2], [4, 2, 1], # Top pyramid (pz, mx, my), ...
        [5, 3, 1], [5, 0, 3], [5, 2, 0], [5, 1, 2]   # Bottom pyramid (mz, my, mx), ...
    ], dtype=np.int32)
    vertices_list = []
    normals_list  = []
    for face_def_indices in octa_face_indices:
        c1, c2, c3 = pv_np[face_def_indices,:]
        octahedron_sphere_face(c1, c2, c3, vertices_list, normals_list, nsub=nsub, radius=radius)
    final_vertices = np.array(vertices_list, dtype=np.float32)
    final_normals  = np.array(normals_list,  dtype=np.float32)
    # print("final_vertices.shape", final_vertices.shape ) #, final_vertices ) # Debug: too verbose
    # print("final_normals.shape",  final_normals.shape ) #, final_normals )  # Debug: too verbose
    return final_vertices, final_normals

class BaseGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.zoom_factor    = 10.0  # Initial zoom
        self.orientation    = R.identity()
        self.last_mouse_pos = None
        self.shader_program = None
        self.default_sphere_mesh = None # For common sphere rendering

        # Camera and light
        self.view_matrix       = QMatrix4x4()
        self.projection_matrix = QMatrix4x4()
        self.light_pos         = QVector3D(15.0, 15.0, 30.0)
        self.light_color       = QVector3D( 1.0,  1.0,  1.0)
        self.camera_pos        = QVector3D( 0  ,  0  , self.zoom_factor) # Will be updated by zoom
        self.current_shader_program_id = None # To be set by derived class if it manages shaders

    def initializeGL_base(self, vertex_shader_src, fragment_shader_src, bPrint=False):
        # Modern OpenGL context should be requested via QSurfaceFormat in main
        if bPrint:
            print(f"OpenGL Version: {glGetString(GL_VERSION).decode()}")
            print(f"GLSL Version: {glGetString(GL_SHADING_LANGUAGE_VERSION).decode()}")

        glClearColor(1.0, 1.0, 1.0, 1.0) # White background
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_CULL_FACE) # Cull back faces


        # Compile shaders if sources are provided
        if vertex_shader_src and fragment_shader_src:
            self.shader_program = self.compile_shader_program(vertex_shader_src, fragment_shader_src)
        elif not hasattr(self, 'shader_program') or not self.shader_program:
             print("BaseGLWidget: No shader sources provided and no shader_program pre-set.")
             # self.shader_program should be set by derived class if sources are None

        # Initialize a default sphere mesh
        #sphere_v, sphere_n, sphere_idx = create_sphere_mesh(radius=1.0)
        #self.default_sphere_mesh = Mesh(vertices=sphere_v, normals=sphere_n, indices=sphere_idx)

        sphere_v, sphere_n = octahedron_sphere_mesh(radius=1.3, nsub=2)
        self.default_sphere_mesh = Mesh(vertices=sphere_v, normals=sphere_n)
        self.default_sphere_mesh.setup_buffers()

    def compile_shader_program(self, vertex_shader_src, fragment_shader_src):
        try:
            return compileProgram(
                compileShader(vertex_shader_src,   GL_VERTEX_SHADER),
                compileShader(fragment_shader_src, GL_FRAGMENT_SHADER)
            )
        except RuntimeError as e:
            print(f"Shader Compilation Error: {e}")
            sys.exit(1)

    def cleanupGL_base(self):
        if self.shader_program:
            glDeleteProgram(self.shader_program)
        self.shader_program = None
        if self.default_sphere_mesh:
            self.default_sphere_mesh.cleanup()

    def paintGL_base(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        

        # Camera / View transformation
        self.camera_pos = QVector3D(0, 0, self.zoom_factor)
        self.view_matrix.setToIdentity()
        self.view_matrix.lookAt(self.camera_pos, QVector3D(0, 0, 0), QVector3D(0, 1, 0))

        # Model transformation (from trackball)
        model_matrix_np = self._rotation_to_gl_matrix(self.orientation)
        model_qmatrix = QMatrix4x4(model_matrix_np.reshape(4,4).T.flatten().tolist())

        programs_to_update = self.all_shader_programs
        
        for prog_id in programs_to_update:
            #print(f"BaseGLWidget.paintGL_base: Updating uniforms for program: {prog_id}")
            #if prog_id is None or prog_id == 0: # Check for 0 as it's an invalid program ID
            #    print(f"BaseGLWidget.paintGL_base: Skipping invalid prog_id: {prog_id}")
            #    continue
            
            # print(f"BaseGLWidget.paintGL_base: Using program {prog_id} for common uniforms")
            glUseProgram(prog_id)
            
            # Check and set common uniforms
            # It's good to check locations, though if a shader doesn't use a uniform, glGetUniformLocation returns -1, and glUniform* with -1 is a no-op.
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "projection"), 1, GL_FALSE, self.projection_matrix.data()) # pyArgs: (location, count, transpose, value_ptr)
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "view"),       1, GL_FALSE, self.view_matrix.data())
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "model"),      1, GL_FALSE, model_qmatrix.data())
            
            glUniform3fv(glGetUniformLocation(prog_id, "lightPos"),   1, [self.light_pos.x(),   self.light_pos.y(),   self.light_pos.z()])
            glUniform3fv(glGetUniformLocation(prog_id, "viewPos"),    1, [self.camera_pos.x(),  self.camera_pos.y(),  self.camera_pos.z()])
            glUniform3fv(glGetUniformLocation(prog_id, "lightColor"), 1, [self.light_color.x(), self.light_color.y(), self.light_color.z()])
        
        self.draw_scene() # Call specific drawing method of derived class

        glUseProgram(0) # Unbind shader after this widget's drawing pass is complete

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        self.projection_matrix.setToIdentity()
        aspect_ratio = width / height if height > 0 else 1
        self.projection_matrix.perspective(45.0, aspect_ratio, 0.1, 200.0)

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.last_mouse_pos = event.pos()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.LeftButton and self.last_mouse_pos:
            p1 = self._project_to_sphere(self.last_mouse_pos.x(), self.last_mouse_pos.y())
            p2 = self._project_to_sphere(event.x(), event.y())
            
            rotation_axis = np.cross(p1, p2)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-6:
                rotation_axis /= rotation_axis_norm
                dot_product = np.dot(p1, p2)
                rotation_angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
                delta_rotation = R.from_rotvec(rotation_angle * rotation_axis)
                self.orientation = delta_rotation * self.orientation
            
            self.last_mouse_pos = event.pos()
            self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom_factor -= delta * 0.01 
        self.zoom_factor = max(1.0, min(self.zoom_factor, 100.0)) # Clamp zoom
        self.update()

    def _project_to_sphere(self, x, y):
        width = self.width()
        height = self.height()
        radius = min(width, height) * 0.4 

        vx = (x - width / 2.0) / radius
        vy = (height / 2.0 - y) / radius 

        vz_squared = 1.0 - vx * vx - vy * vy
        if vz_squared > 0:
            vz = np.sqrt(vz_squared)
        else:
            norm_xy = np.sqrt(vx*vx + vy*vy)
            if norm_xy > 1e-6: 
                vx /= norm_xy
                vy /= norm_xy
            vz = 0.0
        
        p = np.array([vx, vy, vz])
        norm_p = np.linalg.norm(p)
        return p / norm_p if norm_p > 1e-6 else np.array([0.0, 0.0, 1.0])

    def _rotation_to_gl_matrix(self, rotation_obj):
        mat3x3 = rotation_obj.as_matrix()
        mat4x4 = np.eye(4, dtype=np.float32)
        mat4x4[:3, :3] = mat3x3
        return mat4x4.flatten(order='F')

    def draw_scene(self):
        # This method MUST be overridden by derived classes to perform actual drawing.
        # It is called after common uniforms and transformations are set.
        # The shader program is already in use.
        pass

class AppWindow(QMainWindow):
    def __init__(self, parent=None, **kwargs ): # Accept kwargs to pass to QMainWindow if needed
        super().__init__(parent, **kwargs)
        # Common window initialization can go here if any, e.g., default title
        # self.setWindowTitle("Application Window") 

    @classmethod
    def launch(cls, *args, **kwargs):
        """
        Creates and runs the Qt application with an instance of this window class (or a derived class).
        This function will block until the GUI is closed.
        Any *args and **kwargs will be passed to the constructor of `cls`.
        """
        app = QApplication.instance()
        if not app: # Create QApplication if it doesn't exist
            gl_format = QSurfaceFormat()
            gl_format.setVersion(3, 3)
            gl_format.setProfile(QSurfaceFormat.CoreProfile)
            QSurfaceFormat.setDefaultFormat(gl_format)
            # Pass sys.argv if available, otherwise an empty list for robustness
            app = QApplication(sys.argv if hasattr(sys, 'argv') and sys.argv is not None else [])
        
        main_window = cls(*args, **kwargs) # Instantiate the class `cls`
        main_window.show()
        return app.exec_()
```


## xyz_view_new.py

```python
from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import ( QVBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
import argparse

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes

from OpenGL.GL import glUseProgram

class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.trj = []
        self.current_frame_index = 0
        self.opacity      = 0.5
        self.frames_atoms = []
        self.frames_elecs = []
        # self.blend_mode = alpha_blend_modes["standard"] # Replaced by render_modes
        self.shader_program_sphere_raytrace = None
        self.shader_program_sphere_max_vol = None
        self.instance_data_dirty = True
        self.electron_names = set(["E", "e", "e2", "e+", "e-"])
        self.render_modes = {}
        self.current_render_mode_key = None

    @property
    def all_shader_programs(self):
        programs = []
        if self.shader_program_sphere_raytrace: programs.append(self.shader_program_sphere_raytrace)
        if self.shader_program_sphere_max_vol:  programs.append(self.shader_program_sphere_max_vol)
        return programs

    def initializeGL(self):
        #print("MolViewerWidget.initializeGL()")
        # Call base class initialization for shaders, basic GL setup
        shader_folder="../../pyBall/GUI/shaders"
        vert_instances   = open(shader_folder + "/instances.glslv").read()
        frag_ray_sphere  = open(shader_folder + "/sphere.glslf").read()
        frag_max_sphere  = open(shader_folder + "/sphere_max.glslf").read() 
        #FRAGMENT_SHADER_MAX_VOL_SOURCE = open(shader_folder + "/sphere_simple.glslf").read() 
        self.shader_program_sphere_raytrace = self.compile_shader_program( vert_instances, frag_ray_sphere)
        self.shader_program_sphere_max_vol  = self.compile_shader_program( vert_instances, frag_max_sphere)

        # Set the default shader program for BaseGLWidget's uniform setup
        self.shader_program = self.shader_program_sphere_raytrace
        if not self.shader_program: # Fallback if primary fails
            self.shader_program = self.shader_program_sphere_max_vol
        
        # Define render modes after shaders are compiled
        self.render_modes = {
            "Shiny"      : (self.shader_program_sphere_raytrace, alpha_blend_modes["standard"],    True  ),
            "Shiny2"     : (self.shader_program_sphere_raytrace, alpha_blend_modes["standard"],    False ),
            "Volumetric" : (self.shader_program_sphere_max_vol,  alpha_blend_modes["standard"],    True  ), # Example
            "Volumetric2": (self.shader_program_sphere_max_vol,  alpha_blend_modes["subtractive"], True  ),
        }
        # Set a default render mode
        self.current_render_mode_key = "Shiny" 

        self.main_window.render_mode_combo.addItems(self.render_modes.keys())
        self.main_window.render_mode_combo.setCurrentText(self.current_render_mode_key)
        #print("!!!!!! self.render_modes", self.render_modes)
        #print("!!!!!! self.render_mode_combo.currentText()", self.main_window.render_mode_combo.currentText())

        super().initializeGL_base(None, None) # Pass None for shaders, as we compile them here
        sphere_mesh=self.default_sphere_mesh
        self.atom_instances = InstancedData(base_attrib_location=2)
        self.atom_instances.associate_mesh(sphere_mesh)
        atribs = [
            ("positions", 0, 3),
            ("radii",     1, 1),
            ("colors",    2, 4),
        ]
        self.sphere_vbo_inds = self.atom_instances.setup_instance_vbos(atribs)
        self.elec_instances = InstancedData(base_attrib_location=2) # Same attrib locations
        self.elec_instances.associate_mesh(sphere_mesh)
        self.elec_instances.setup_instance_vbos(atribs)
    
    def cleanupGL(self):
        # super().cleanupGL_base() # BaseGLWidget.shader_program will be one of these
        if self.shader_program_sphere_raytrace:
            glDeleteProgram(self.shader_program_sphere_raytrace)
            self.shader_program_sphere_raytrace = None
        if self.shader_program_sphere_max_vol:
            glDeleteProgram(self.shader_program_sphere_max_vol)
            self.shader_program_sphere_max_vol = None
        if self.atom_instances:
            self.atom_instances.cleanup()
        if self.elec_instances:
            self.elec_instances.cleanup()

    def precalculate_frames_data(self, dtype=np.float32):
        for frame_idx in range(len(self.trj)):
            es,ps,qs, rs, comment = self.trj[frame_idx]
            na = len(es)
            # atoms = SphereInstancesData()
            # elecs = SphereInstancesData()
            apos,arad,acol=[],[],[]
            epos,erad,ecol=[],[],[]
            if na > 0:
                for i in range(na):
                    atom_symbol = es[i]
                    current_radius = rs[i]
                    element_data = elements.ELEMENT_DICT[atom_symbol]
                    hex_color_str = element_data[elements.index_color]
                    r_col, g_col, b_col = elements.hex_to_float_rgb(hex_color_str)
                    if atom_symbol in self.electron_names: # Electron or transparent particle
                        epos.append(ps[i])
                        erad.append(current_radius)
                        ecol.append([r_col, g_col, b_col, self.opacity])
                    else: 
                        apos.append(ps[i])
                        arad.append(current_radius)
                        acol.append([r_col, g_col, b_col, 1.0])
            
            self.frames_atoms.append( [np.array(apos,dtype=dtype), np.array(arad,dtype=dtype), np.array(acol,dtype=dtype)] )
            self.frames_elecs.append( [np.array(epos,dtype=dtype), np.array(erad,dtype=dtype), np.array(ecol,dtype=dtype)] )
            # if frame_idx == 0:
            #     print("arad",   self.frames_atoms[-1][1])
            #     print("apos\n", self.frames_atoms[-1][0])
            #     print("acol\n", self.frames_atoms[-1][2])
            #     print("erad",   self.frames_elecs[-1][1])
            #     print("epos\n", self.frames_elecs[-1][0])
            #     print("ecol\n", self.frames_elecs[-1][2])

    def update_instance_data(self):
        elecs = self.frames_elecs[self.current_frame_index]
        elecs[2][:,3] = self.opacity
        self.atom_instances.update_list( self.frames_atoms[self.current_frame_index] )
        self.elec_instances.update_list( elecs )
        # self.atom_instances     .update( {
        #     "positions": self.frames_atoms[self.current_frame_index][0],
        #     "radii":     self.frames_atoms[self.current_frame_index][1],
        #     "colors":    self.frames_atoms[self.current_frame_index][2]
        #     })
        # self.elec_instances.update( {
        #     "positions": self.frames_elecs[self.current_frame_index][0],
        #     "radii":     self.frames_elecs[self.current_frame_index][1],
        #     "colors":    self.frames_elecs[self.current_frame_index][2]
        #     })
        self.instance_data_dirty = False

    def paintGL(self):
        # Calls paintGL_base, which in turn calls draw_scene
        super().paintGL_base() 

    def draw_scene(self):

        shader_id_electrons, blend_params_electrons, depth_test = self.render_modes[self.current_render_mode_key]

        if self.instance_data_dirty:
            self.update_instance_data()
        self.use_shader(self.shader_program_sphere_raytrace)
        disable_blend() # Opaque objects don't need blending
        self.atom_instances.draw()
        
        # Transparent electrons (shader and blend mode depend on selection)
        
        self.use_shader(shader_id_electrons)
        set_ogl_blend_mode(blend_params_electrons, depth_test)
        self.elec_instances.draw()
        
    def use_shader(self, shader_prog_id):
        glUseProgram(shader_prog_id)
        self.current_shader_program_id = shader_prog_id # For BaseGLWidget to set uniforms

    def set_frame(self, frame_idx):
        if 0 <= frame_idx < len(self.trj):
            self.current_frame_index = frame_idx
            self.instance_data_dirty = True
            self.update()

    def set_opacity(self, opacity_percent):
        self.opacity = opacity_percent / 100.0
        self.instance_data_dirty = True # Need to update colors VBO
        self.update()

    def set_render_mode(self, mode_key):
        if mode_key in self.render_modes:
            self.current_render_mode_key = mode_key
            self.update() # Trigger repaint
        else:
            print(f"Warning: Render mode '{mode_key}' not found.")


class MolViewer(AppWindow):
    def __init__(self, trj):
        super().__init__()
        self.setWindowTitle("Modern OpenGL Molecular Viewer")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.gl_widget = MolViewerWidget()
        layout.addWidget(self.gl_widget, 1) # <-- Add stretch factor here

        self.frame_slider = QSlider(Qt.Horizontal)
        self.frame_slider.valueChanged.connect(self.gl_widget.set_frame)
        layout.addWidget(QLabel("Frame:"))
        layout.addWidget(self.frame_slider)

        self.opacity_slider = QSlider(Qt.Horizontal)
        self.opacity_slider.setMinimum(0)
        self.opacity_slider.setMaximum(100)
        self.opacity_slider.setValue(int(self.gl_widget.opacity * 100))
        self.opacity_slider.valueChanged.connect(self.gl_widget.set_opacity)
        layout.addWidget(QLabel("Opacity:"))
        layout.addWidget(self.opacity_slider)

        self.gl_widget.trj = trj
        self.gl_widget.update()
        if self.gl_widget.trj:
            self.frame_slider.setMinimum(0)
            self.frame_slider.setMaximum(len(self.gl_widget.trj) - 1)
        self.gl_widget.precalculate_frames_data()
        self.gl_widget.update()

        self.render_mode_combo = QComboBox()

        self.gl_widget.main_window = self

        self.render_mode_combo.currentTextChanged.connect(self.on_render_mode_changed)
        layout.addWidget(QLabel("Electron Render Mode:"))
        layout.addWidget(self.render_mode_combo)

        self.show()

    def on_render_mode_changed(self, text):
        self.gl_widget.set_render_mode(text)


if __name__ == '__main__':
    from pyBall import atomicUtils as au
    #app = QApplication(sys.argv)
    parser = argparse.ArgumentParser(description="Modern OpenGL Molecular Viewer")
    parser.add_argument("-f", "--file", type=str, help="Path to the XYZ trajectory file", default=None) # Default to None
    args = parser.parse_args()

    trj = au.load_xyz_movie(args.file)
    trj = au.trj_to_ename(trj)
    trj = au.trj_fill_radius(trj, bVdw=True, rFactor=0.001, rmin=0.05) # Adjusted rFactor for visibility
    #trj = au.trj_fill_radius(trj, bVdw=False, rFactor=1.0)
    print( "trj.enames", trj[0])
    MolViewer.launch(trj=trj)

    #main_window = MolViewer( trj=trj)
    #sys.exit(app.exec_())

```

## /AtomicSystem.md

# `AtomicSystem.py` Documentation

This document provides a comprehensive overview of the `AtomicSystem.py` module, which defines the `AtomicSystem` class. This class serves as a fundamental data structure for representing and manipulating atomic and molecular systems within the `FireCore` project. It leverages many utility functions from `atomicUtils.py` (documented in `atomicUtils.md`).

## 1. `AtomicSystem` Class

### 1.1. Purpose

The `AtomicSystem` class encapsulates all relevant data for an atomic system, including atomic coordinates, types, element names, charges, radii, bond connectivity, neighbor lists, and simulation cell lattice vectors. It provides a rich set of methods for loading, saving, querying, analyzing, and manipulating these systems.

### 1.2. Initialization (`__init__`)

```python
AtomicSystem(fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True, bPreinit=True)
```

- **`fname` (str, optional):** Path to a file (`.mol`, `.mol2`, `.xyz`) to load the atomic system from.
- **`apos` (np.ndarray, optional):** (N, 3) array of atomic positions.
- **`atypes` (np.ndarray, optional):** (N,) array of atomic types (integer indices).
- **`enames` (list of str, optional):** List of element names (e.g., `['C', 'H', 'O']`).
- **`lvec` (np.ndarray, optional):** (3, 3) array representing lattice vectors for periodic boundary conditions.
- **`qs` (np.ndarray, optional):** (N,) array of atomic charges.
- **`Rs` (np.ndarray, optional):** (N,) array of atomic radii.
- **`bonds` (np.ndarray, optional):** (M, 2) array of bond indices `[[i, j], ...]`. If not provided, bonds can be found later using `findBonds()`.
- **`ngs` (list of dict, optional):** List of dictionaries representing neighbor lists for each atom.
- **`bReadN` (bool, optional):** If `True`, reads the number of atoms from the file header (for some formats). Default is `True`.
- **`bPreinit` (bool, optional):** If `True` (default), calls `preinitialize_atomic_properties()` after loading/initialization to set up `iZs` (atomic numbers), `masses`, and `Rs` (radii) based on `enames`.

**Note:** When initialized from a file, `AtomicSystem` internally uses functions from `atomicUtils.py` (e.g., `au.loadMol`, `au.load_xyz`).

### 1.3. Attributes

- **`apos` (np.ndarray):** Atomic positions (N, 3).
- **`atypes` (np.ndarray):** Atomic types (N,).
- **`enames` (list of str):** Element names (N,).
- **`qs` (np.ndarray):** Atomic charges (N,).
- **`Rs` (np.ndarray):** Atomic radii (N,).
- **`bonds` (np.ndarray):** Bond connectivity (M, 2).
- **`ngs` (list of dict):** Neighbor lists for each atom.
- **`lvec` (np.ndarray):** Lattice vectors (3, 3).
- **`iZs` (np.ndarray):** Atomic numbers (N,).
- **`masses` (np.ndarray):** Atomic masses (N,).
- **`aux_labels` (list, optional):** Auxiliary labels for atoms.

### 1.4. File I/O Methods

These methods facilitate saving the current atomic system to various file formats.

- **`saveXYZ(fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None)`:**
    Saves the system to an XYZ file. Can include lattice vectors in the comment and optionally save charges and radii. Internally uses `atomicUtils.saveXYZ`.
- **`save_mol(fname, title="Avogadro")`:**
    Saves the system to a MOL file. Internally uses `atomicUtils.save_mol`.
- **`save_mol2(fname, comment="")`:**
    Saves the system to a MOL2 file. Internally uses `atomicUtils.save_mol2`.
- **`toLines()`:**
    Returns a list of strings, each representing an atom in a format suitable for geometry files (e.g., `"C   0.000   0.000   0.000"`). Internally uses `atomicUtils.geomLines`.
- **`toXYZ(fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False)`:**
    Writes the system data to an already open file object `fout` in XYZ format. Internally uses `atomicUtils.writeToXYZ`.

### 1.5. Information & Query Methods

These methods provide ways to inspect and retrieve information about the atomic system.

- **`print()`:**
    Prints a summary of the system, including the number of atoms and details for each atom (index, type, element name, position, and auxiliary labels if present).
- **`getValenceElectrons()`:**
    Returns a NumPy array containing the number of valence electrons for each atom in the system.
- **`subtractValenceE(f0=-1.0, f=+1.0)`:**
    Adjusts the `qs` (charges) array by subtracting or adding valence electron counts, scaled by `f0` and `f`.
- **`printBonds()`:**
    Prints a list of all defined bonds in the system.
- **`printNeighs()`:**
    Prints the neighbor list for each atom.
- **`findBonds(Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True)`:**
    Identifies and stores bonds within the system based on interatomic distances and optional Van der Waals radii. Internally uses `atomicUtils.findBondsNP`.
- **`findHBonds(Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=au.neg_types_set, bPrint=False, bHbase=False)`:**
    Identifies hydrogen bonds within the system based on distance and angular criteria. Internally uses `atomicUtils.findHBondsNP`.
- **`findBondsOfAtom(ia, bAtom=False)`:**
    Returns a list of bonds connected to atom `ia`. If `bAtom` is `True`, returns atom indices, otherwise bond indices.
- **`neighs(bBond=True)`:**
    Generates and returns the neighbor list (`self.ngs`). If `bBond` is `True`, neighbors are determined by existing bonds; otherwise, by proximity (using `atomicUtils.neigh_atoms`).
- **`find_groups()`:**
    Identifies and returns a list of atom indices for each disconnected molecular fragment (bonded cluster) in the system. Internally uses `atomicUtils.selectBondedCluster`.
- **`select_by_ename(elist)`:**
    Returns a boolean mask or list of indices for atoms whose element names are present in `elist`.
- **`getNeighsOfType(selection, typ='N')`:**
    Given a `selection` of atom indices, returns the indices of their neighbors that match the specified `typ`. Internally uses `atomicUtils.findNeighsOfType`.
- **`select_by_neighType(neighs, typ='N', neighTyps={'H':(1,2)})`:**
    Selects atoms based on the types of their neighbors. Internally uses `atomicUtils.findTypeNeigh_`.
- **`findAngles(select=None, ngs=None)`:**
    Finds and returns bond angles within the system. `select` can specify a subset of atoms.
- **`findDihedral(select=None, ngs=None, neighTyp={'H'})`:**
    Finds and returns dihedral angles within the system.
- **`findCOG(apos, byBox=False)`:**
    Calculates the center of geometry (COG) for a given set of positions `apos`. Internally uses `atomicUtils.findCOG`.
- **`projectAlongBondDir(i0, i1)`:**
    Projects atomic positions along the direction defined by the bond between atoms `i0` and `i1`. Internally uses `atomicUtils.projectAlongBondDir`.
- **`store_bond_lengths()`:**
    Calculates and stores the current lengths of all bonds in `self.bond_lengths`.
- **`restore_bond_length(ij, L=None)`:**
    Restores the length of a specific bond `ij` to its stored value or a provided `L`.

### 1.6. System Manipulation & Transformation Methods

These methods allow for modification, transformation, and assembly of atomic systems.

- **`clonePBC(nPBC=(1,1,1))`:**
    Clones the current system across periodic boundary conditions specified by `nPBC` (e.g., `(2,2,1)` for a 2x2x1 supercell).
- **`symmetrized(d=0.1)`:**
    Applies symmetry operations to the system (details depend on internal implementation).
- **`selectSubset(inds)`:**
    Returns a new `AtomicSystem` object containing only the atoms specified by `inds`.
- **`selectBondedCluster(s)`:**
    Selects a bonded cluster of atoms starting from a seed set `s` and returns a new `AtomicSystem` object for that cluster. Internally uses `atomicUtils.selectBondedCluster`.
- **`makeRotMat(ip1, ip2, _0=1)`:**
    Creates a rotation matrix based on two atom indices `ip1` and `ip2`. Internally uses `atomicUtils.makeRotMat`.
- **`orient_mat(rot, p0=None, bCopy=False)`:**
    Orients the system by applying a given rotation matrix `rot`. `p0` specifies the origin of rotation. If `bCopy` is `True`, returns a new `AtomicSystem`.
- **`orient_vs(fw, up, p0=None, trans=None, bCopy=False)`:**
    Orients the system such that the `fw` (forward) and `up` vectors align with the system's axes. Internally uses `atomicUtils.orient_vs`.
- **`orient(i0, b1, b2, _0=1, trans=None, bCopy=False)`:**
    Orients the molecule by aligning three specified atom indices (`i0`, `b1`, `b2`) to a target orientation. Internally uses `atomicUtils.orient`.
- **`orientPCA(perm=None)`:**
    Orients the system by aligning its principal components of inertia with the coordinate axes. Internally uses `atomicUtils.orientPCA`.
- **`shift(vec, sel=None)`:**
    Shifts all or a `sel`ected subset of atoms by a given `vec`tor.
- **`rotate_ax(ang, ax=(0,1), p0=None)`:**
    Rotates all or a selected subset of atoms by `ang` around a specified `ax`is passing through `p0`.
- **`delete_atoms(lst)`:**
    Deletes atoms specified by the list of indices `lst` from the system. It reindexes bonds and neighbor lists to maintain consistency. Internally uses `atomicUtils.reindex_bonds` and `atomicUtils.make_reindex`.
- **`add_atom(pos, ename, atype=-1, q=0.0, R=1.0)`:**
    Adds a new atom to the system with specified `pos`ition, `ename` (element name), `atype` (type), `q` (charge), and `R` (radius).
- **`add_bond(i, j)`:**
    Adds a new bond between atoms `i` and `j`.
- **`merge(other_system, rot=None, trans=None)`:**
    Merges another `AtomicSystem` object (`other_system`) into the current system. The `other_system` can optionally be rotated and translated before merging. It adjusts bond indices accordingly.
- **`attach_group(G, i0, i1, iup, bond, up=(0., 0., 1.), _0=1, pre="A")`:**
    Attaches an end-group (`G`, which is another `AtomicSystem` object) to the current system's backbone at a specified bond. Internally uses `atomicUtils.attach_group`.
- **`attach_group_by_marker(G, markerX="Xe", markerY="He", _0=1, pre="X")`:**
    Attaches an end-group (`G`) to the current system using marker atoms (`markerX`, `markerY`) to define the attachment point and orientation. Internally uses `atomicUtils.attach_group_by_marker`.

### 1.7. Electron Pair (`Epair`) Geometry Methods

These methods are used for calculating and placing electron pairs (e.g., lone pairs, pi-electron pairs) based on atomic bonding configurations.

- **`get_atomi_pi_direction(i)`:**
    Determines the direction of a pi-orbital for atom `i`, typically used for atoms involved in double or triple bonds.
- **`make_epair_geom(i, npi, nb)`:**
    Calculates and places electron pairs around atom `i` based on its number of pi-bonds (`npi`) and number of neighbors (`nb`). It handles different geometries (e.g., like NH3, H2O, =N-, =O).
- **`place_electron_pair(i, direction, distance=0.5, ename='E', atype=200, qs=0.0, Rs=1.0)`:**
    Adds a new "electron pair" to the system. This is represented as a new atom (typically with `ename='E'`) at a calculated position relative to atom `i` and updates `apos`, `atypes`, `enames`, `qs`, and `Rs` arrays. It also adds a bond between atom `i` and the new electron pair.

## 2. Usage Example (Conceptual)

```python
import numpy as np
from pyBall.AtomicSystem import AtomicSystem

# 1. Initialize from a file
sys1 = AtomicSystem(fname="my_molecule.xyz")
sys1.findBonds() # Find bonds if not loaded from file
sys1.print()

# 2. Initialize from arrays
positions = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0]
], dtype=np.float32)
e_names = ['C', 'H', 'H']
sys2 = AtomicSystem(apos=positions, enames=e_names)
sys2.add_bond(0, 1)
sys2.add_bond(0, 2)
sys2.printBonds()

# 3. Manipulate the system
sys1.shift(np.array([10.0, 0.0, 0.0]))
sys1.rotate_ax(np.pi/2, ax=(0,1)) # Rotate around Z-axis

# 4. Merge systems
merged_sys = AtomicSystem(apos=np.array([]), enames=[]) # Create an empty system
merged_sys.merge(sys1)
merged_sys.merge(sys2, trans=np.array([5.0, 5.0, 0.0])) # Merge sys2 with translation
merged_sys.saveXYZ("merged_system.xyz")

# 5. Find and place electron pairs (example for a specific atom)
# Assuming sys1 has an oxygen atom at index 0 and bonds are defined
# sys1.make_epair_geom(0, npi=0, nb=2) # Example for H2O-like oxygen
# sys1.saveXYZ("molecule_with_epairs.xyz")
```

This documentation provides a solid foundation for understanding and utilizing the `AtomicSystem` class for various atomic and molecular system manipulations. For detailed parameter descriptions and specific behaviors, always refer to the source code of `AtomicSystem.py` and `atomicUtils.py`.


## AtomicSystem.py

```python
#!/usr/bin/python

from random import random
import numpy as np
from . import elements
#import elements
#import numpy as np
import copy
from . import atomicUtils as au

VALENCE_DICT={
#        nBond  nEpair
'O':   ( 2,     2  ),
'N':   ( 3,     1  ),
}

class AtomicSystem( ):

    def __init__(self,fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True, bPreinit=True ) -> None:
        self.apos    = apos
        self.atypes  = atypes
        self.enames  = enames
        self.qs      = qs
        self.Rs      = Rs
        self.bonds   = bonds
        self.ngs     = ngs 
        self.lvec    = lvec
        self.aux_labels = None
        if fname is not None:
            ext = fname.split('.')[-1]
            #print( f"AtomicSystem.__init__({fname}) ext=", ext  )
            if( 'mol' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds = au.loadMol(fname=fname, bReadN=bReadN )
            elif ( 'mol2' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds, self.lvec = au.loadMol2(fname=fname, bReadN=bReadN )
            elif ( 'xyz' == ext ):
                self.apos,self.atypes,self.enames,self.qs, comment = au.load_xyz(fname=fname, bReadN=bReadN )
                if comment is not None:
                    if comment[:3] == 'lvs':      
                        self.lvec = au.string_to_matrix( comment, nx=3,ny=3, bExactSize=False )
                        #print( f"AtomicSystem.__init__({fname}) lvec=\n", self.lvec   )
                #print( f"AtomicSystem.__init__({fname}) comment=", comment  )
            else:
                self.apos,self.atypes,self.enames,self.qs = au.loadAtomsNP(fname=fname , bReadN=bReadN )
            if bPreinit:
                self.preinitialize_atomic_properties()

    def saveXYZ(self, fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None ):
        if blvec and (self.lvec is not None):
            #print( self.lvec )
            comment= ( "lvs %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f" %(self.lvec[0,0],self.lvec[0,1],self.lvec[0,2],  self.lvec[1,0],self.lvec[1,1],self.lvec[1,2],  self.lvec[2,0],self.lvec[2,1],self.lvec[2,2]   ) ) + comment
        qs = self.qs
        if(not bQs): qs=None
        au.saveXYZ( self.enames, self.apos, fname, qs=qs, Rs=self.Rs, mode=mode, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def save_mol(self, fname, title="Avogadro"):
        au.save_mol(fname, self.enames, self.apos, self.bonds, title="Avogadro")

    def save_mol2(self, fname, comment=""):
        au.save_mol2(fname, self.enames, self.apos, self.bonds, comment="")

    def toLines(self):
        #lines = []
        #for i,pos in enumerate(self.apos):
        #    lines.append(  "%s %3.5f %3.5f %3.5f\n" %(self.enames[i], pos[0],pos[1],pos[2]) )
        return au.geomLines( self.apos, self.enames )

    def toXYZ(self, fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False ):
        au.writeToXYZ( fout, self.enames, self.apos, qs=self.qs, Rs=self.Rs, bHeader=bHeader, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def print(self):
        print( len(self.atypes), len(self.enames), len(self.apos) )
        for i in range(len(self.apos)):
            print( "[%i] %i=%s p(%10.5f,%10.5f,%10.5f)" %( i, self.atypes[i],self.enames[i], self.apos[i,0], self.apos[i,1], self.apos[i,2] ), end =" " )
            if(self.aux_labels is not None): print(self.aux_labels[i], end =" ")
            print("")

    def getValenceElectrons( self ):
        return  np.array( [ elements.ELEMENT_DICT[e][9] for e in self.enames ] )

    def subtractValenceE(self, f0=-1.0, f=+1.0 ):
        self.qs[:] = self.qs[:]*f0 + self.getValenceElectrons()*f       

    def printBonds(self):
        print("AtomicSystem.printBonds():")
        if self.bonds is None:
            print("No bonds defined")
            return
        for i, (a, b) in enumerate(self.bonds):
            print(f"[{i}] ({a},{b}) ({self.enames[a]},{self.enames[b]})")

    def printNeighs(self):
        print("AtomicSystem.printNeighs():")
        if self.neighs is None:
            print("No neighs defined")
            return
        for i, ngi in enumerate(self.ngs):
            print(f"ngs[{i}]: ", end="")
            for j,ia in enumerate(ngi):
                print(ia, end=" ")
            print("")

    def findBonds(self, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True ):
        if self.atypes is None:
            self.atypes = [ elements.ELEMENT_DICT[e][0] for e in self.enames ]
        self.bonds, rs = au.findBondsNP( self.apos, self.atypes, Rcut=Rcut, RvdwCut=RvdwCut, RvdWs=RvdWs, byRvdW=byRvdW )
        return self.bonds, rs

    def findHBonds(self, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=au.neg_types_set, bPrint=False, bHbase=False ):
        return au.findHBondsNP( self.apos, atypes=self.enames, Rb=Rb, Rh=Rh, angMax=angMax, typs1=typs1, typs2=typs2, bPrint=bPrint,  bHbase=bHbase )

    def findBondsOfAtom(self, ia, bAtom=False ):
        if bAtom: 
            return [ b[1] for b in self.bonds if(b[0]==ia) ] + [ b[0] for b in self.bonds if(b[1]==ia) ] 
        else:
            return [i for i,b in enumerate(self.bonds) if (b[0]==ia) or (b[1]==ia) ]

    def neighs( self, bBond=True ):
        if(self.bonds is None):
            self.findBonds()
        self.ngs = au.neigh_bonds( len(self.apos), self.bonds )
        return self.ngs

    def find_groups(self):
        if self.ngs is None: self.neighs()
        ngs = self.ngs
        #print( ngs )
        groups = { }
        for inod in range(len(self.apos)):
            if len(ngs[inod]) > 1: groups[inod] = [inod]
        for inod,g in groups.items():
            inod = g[0] 
            g += [ ia for ia in ngs[inod].keys() if ia not in groups ] 
        return groups

    def select_by_ename( self, elist ):
        return [ i for i,e in enumerate(self.enames) if e in elist ]

    def getNeighsOfType( self, selection, typ='N'):
        if self.ngs is None: self.neighs()
        return au.findNeighsOfType( selection, self.enames, self.ngs, typ=typ ) 

    def select_by_neighType( self, neighs, typ='N', neighTyps={'H':(1,2)} ):
        return au.findTypeNeigh_( self.enames, neighs, typ=typ, neighTyps=neighTyps )

    def findAngles(self, select=None, ngs=None, ):
        if ngs is None:
            ngs = self.neighs()
        return au.findAngles( self.apos, select=select, neighs=ngs )

    def findDihedral( self, select=None, ngs=None, neighTyp={'H'} ):
        if ngs is None:
            ngs = self.neighs()
        return au.findDihedral( self.apos, self.enames, ngs, select=select, neighTyp=neighTyp ) 

    def findCOG(self, apos, byBox=False ):
        return au.findCOG( apos, byBox=byBox )
    
    def projectAlongBondDir( self, i0, i1 ):
        return au.projectAlongBondDir( self.apos, i0, i1 )

    def store_bond_lengths(self):
        bond_lengths = {}
        bonds = self.findBonds()  # Get all bonds in the system
        for bond in bonds:
            ia,ja = bond
            if ia>ja: ia,ja = ja,ia
            length = np.linalg.norm(self.apos[ia]-self.apos[ja])
            bond_lengths[(ia,ja)] = length
        self.bond_legths = bond_lengths
        return bond_lengths

    def restore_bond_length(self, ij, L=None ):
        ia,ja= ij
        d = self.apos[ja] - self.apos[ia]
        Lnow = np.linalg.norm(d)
        if L is None:
            if ia>ja: i,j = ja,ia
            else:     i,j = ia,ja
            L = self.bond_lengths[(i,j)]
        f = L / Lnow
        self.apos[ia] = self.apos[ja] + d * f


    def clonePBC(self,nPBC=(1,1,1) ):
        nx,ny,nz= nPBC
        nxyz=nx*ny*nz
        na = len(self.apos)
        apos   = np.zeros((na*nxyz,3))
        #print( "clonePBC ", na, len(self.atypes) )
        if self.atypes is not None: 
            atypes = np.zeros(na*nxyz,np.int32)
        else:
            atypes = None

        if self.enames is not None: 
            enames = []
        else:
            enames = None

        if self.qs is not None: 
            qs = np.zeros(na*nxyz) 
        else:
            qs = None

        #print( nxyz, na, apos.shape, atypes.shape )
        if( nxyz > 1 ):
            lvec   = np.array([ self.lvec[0,:]*nx,self.lvec[1,:]*ny,self.lvec[2,:]*nz ]) 
            i0=0
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        shift = self.lvec[0,:]*ix  + self.lvec[1,:]*iy + self.lvec[2,:]*iz
                        apos  [i0:i0+na,:] = self.apos[:,:] + shift[None,:]
                        if atypes is not None: atypes[i0:i0+na  ] = self.atypes
                        if qs     is not None: qs    [i0:i0+na  ] = self.qs    
                        if enames is not None: enames[i0:i0+na  ] = self.enames
                        #if enames is not None: enames += self.enames
                        i0+=na
        else:
            lvec=self.lvec
            apos  [:,:] = self.apos[:,:]
            if atypes is not None: atypes[:] = self.atypes[:]
            if qs     is not None: qs    [:] = self.qs    [:]  
            if enames is not None: enames[:] = self.enames[:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def symmetrized(self, d=0.1 ):
        # def atoms_symmetrized( atypes, apos, lvec, qs=None, REQs=None, d=0.1):
        atypes, apos, qs, REQs, ws = au.atoms_symmetrized( self.atypes, self.apos, self.lvec, qs=self.qs, d=d );
        enames = au.iz2enames( atypes )
        return AtomicSystem( apos=apos, atypes=atypes, enames=enames, lvec=self.lvec.copy(), qs=qs ), ws 

    def selectSubset(self, inds ):
        if self.atypes is not None: 
                atypes = self.atypes[inds]
        else:
            atypes = None

        if self.enames is not None: 
            enames = [ self.enames[i] for i in inds ]
        else:
            enames = None

        if self.qs is not None: 
            qs = self.qs[inds]
        else:
            qs = None

        lvec=self.lvec
        apos  = self.apos[inds,:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def selectBondedCluster( self, s ):
        na = len(self.apos)
        if self.bonds is None: self.findBonds()
        s     = au.selectBondedCluster( s, self.bonds )
        ins  = [ i for i in range(na) if (i in s) ]
        outs = [ i for i in range(na) if (i not in s) ] 
        return ins,outs

    def makeRotMat( self, ip1, ip2, _0=1 ):
        fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        return au.makeRotMat( fw, up )

    def orient_mat(self, rot, p0=None, bCopy=False ):
        apos=self.apos  
        if(bCopy): apos=apos.copy()
        if p0  is not None: apos[:,:]-=p0[None,:]
        if rot is not None: au.mulpos( apos, rot )
        return apos

    def orient_vs(self, fw, up, p0=None, trans=None, bCopy=False ):
        if fw is None:
            rot = None
        else:
            rot = au.makeRotMat( fw, up )
            if trans is not None: rot=rot[trans,:]
        return self.orient_mat( rot, p0, bCopy )

    def orient( self, i0, b1, b2, _0=1, trans=None, bCopy=False ):
        #print( "orient i0 ", i0, " ip1 ", ip1, " ip2 ",ip2 )
        # p0  = self.apos[i0-_0]
        # fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        # up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        p0, fw, up = au.makeVectros( self.apos, i0, b1, b2, _0=_0 )
        return self.orient_vs( fw, up, p0, trans=trans, bCopy=bCopy )
    
    def orientPCA(self, perm=None):
        au.orientPCA(self.apos, perm=perm )

    def shift(self, vec, sel=None ):
        if sel is None: 
            self.apos[:,0] += vec[0]
            self.apos[:,1] += vec[1]
            self.apos[:,2] += vec[2]
        else:
            self.apos[sel,0] += vec[0]
            self.apos[sel,1] += vec[1]
            self.apos[sel,2] += vec[2]

    def rotate_ax(self, ang, ax=(0,1), p0=None ):
        rot = au.makeRotMatAng( ang, ax=ax )
        if p0  is not None: self.apos[:,:]-=p0[None,:]
        au.mulpos( self.apos, rot )
        if p0  is not None: self.apos[:,:]+=p0[None,:]

    def delete_atoms(self, lst ):
        st = set(lst)
        if( self.apos   is not None ): self.apos   =  np.delete( self.apos,   lst, axis=0 )
        if( self.atypes is not None ): self.atypes =  np.delete( self.atypes, lst )
        if( self.qs     is not None ): self.qs     =  np.delete( self.qs,     lst )
        if( self.Rs     is not None ): self.Rs     =  np.delete( self.Rs,     lst )
        if( self.enames is not None ): self.enames =  np.delete( self.enames, lst )
        if( self.aux_labels is not None ): self.aux_labels = [ v for i,v in enumerate(self.aux_labels) if i not in st ] 


    def preinitialize_atomic_properties(self):
        """
        Preinitialize per-atom arrays for an AtomicSystem.
        
        This function assumes that the system’s atypes (or enames) have been set.
        It uses the global 'elements.ELEMENTS' (a list of lists) to set default values:
        - qs: set to the element’s default valence electron count (column index 9)
        - Rs: set to the element’s van der Waals radius (column index 7)
        - aux_labels: set to a default label (simply the atom’s index as a string)
        
        Parameters:
        atomicSystem (AtomicSystem): an instance of AtomicSystem.
        
        Raises:
        ValueError: if atomicSystem.atypes is not defined.
        """
        natoms = len(self.apos)
        if self.atypes is None:   raise ValueError("The system does not have atypes defined. Please initialize the system’s atypes (or enames) first.")
        if self.qs is None:   # Assume atypes is an array of atomic numbers (e.g. 6 for carbon, etc.)
            qs = []
            for z in self.atypes:   # our ELEMENTS list is zero-based: for atomic number z, use ELEMENTS[z-1]
                qs.append(elements.ELEMENTS[z-1][9])
            self.qs = np.array(qs)
        if self.Rs is None:  # For each atom, use the vdW radius (column index 7)
            Rs = []
            for z in self.atypes: Rs.append(elements.ELEMENTS[z-1][7])
            self.Rs = np.array(Rs)
        # Initialize aux_labels if not defined.
        if self.aux_labels is None: self.aux_labels = [str(i) for i in range(natoms)]
        self.neighs()
        # (If you have other arrays you want to preinitialize, do it here.)
        #print(f"Pre-initialized atomic properties for {natoms} atoms.")

        
    def check_atomic_properties(atomicSystem):
        """
        Check that the per-atom arrays (qs, Rs, aux_labels) are defined and
        have the correct length. If not, raise an error telling the user
        to run preinitialize_atomic_properties().
        """
        natoms = len(atomicSystem.apos)
        if (atomicSystem.qs is None or len(atomicSystem.qs) != natoms or
            atomicSystem.Rs is None or len(atomicSystem.Rs) != natoms or
            atomicSystem.aux_labels is None or len(atomicSystem.aux_labels) != natoms):
            raise ValueError("Not all per-atom arrays are initialized correctly. Please call preinitialize_atomic_properties() on your system.")
                            
                            
    # Example: modify append_atoms() to check rather than auto-initialize
    def append_atoms(self, B, pre="A"):
        # Ensure both systems have been pre-initialized:
        self.check_atomic_properties()
        B.check_atomic_properties()
        
        # Number of atoms in self and in B
        nA = len(self.apos)
        nB = len(B.apos)
        
        self.apos   = np.append(self.apos,   B.apos, axis=0)
        self.atypes = np.append(self.atypes, B.atypes)
        self.qs     = np.append(self.qs,     B.qs)
        self.Rs     = np.append(self.Rs,     B.Rs)
        self.enames = np.append(self.enames, B.enames)
        
        # Extend the aux_labels list:
        self.aux_labels.extend(B.aux_labels)


    def remap( self, lst ):
        dct = {   key:value for (value,key) in enumerate(self.aux_labels) }
        return [ dct.get(key,-1) for key in lst ]

    def addBonds( self, added_bonds ):
        """
        Add bonds to the current system.
        added_bonds : list List of bonds to be added. Each bond is a tuple (i,j), where i and j are the indices of the atoms in the current system that are connected by the bond.
        """
        self.bonds = np.append(self.bonds, np.array(added_bonds), axis=0)


    def addSystems(self, other, pos=None, rot=None, added_bonds=None, _0=1 ):
        """
        Add a new system to the current system with optional position and orientation.
        
        Parameters:
        -----------
        other : AtomicSystem
            The system to be added
        pos : np.ndarray, optional
            Position vector (3,) where to place the other system
        rot : np.ndarray, optional
            Rotation matrix (3,3) defining orientation of the other system
            
        This is a simplified merge operation that:
        1. Optionally transforms the coordinates of the other system (rotation and translation)
        2. Adjusts bond indices of the other system by adding offset
        3. Merges atomic data using merge_arrays
        No atoms are removed in this process.
        """

        offset = len(self.apos) # Get current number of atoms as offset for bond indices

        # Make a copy of the other system to avoid modifying it
        other_copy = copy.deepcopy(other)

        # Transform coordinates if needed
        if rot is not None:
            rot=np.array(rot)
            au.mulpos(other_copy.apos, rot)
        if pos is not None:
            pos=np.array(pos)
            other_copy.apos += pos[None,:]  # Broadcasting the translation to all atoms
            
        # Let merge_arrays handle the bond index adjustment
        self.merge_arrays(other_copy, offset)

        # Add the bonds that connect the two systems ( we assume i is from self, j is from other )
        if added_bonds is not None:
            added_bonds = [ (i-_0,j-_0 + offset) for i,j in added_bonds ] 
            self.addBonds( added_bonds )
        
        # Reset neighbor list since we modified atomic indices
        self.ngs = None

    def attach_group( self, G,  i0, i1, iup,   bond,  up=(0.,0.,1.),  _0=1, pre="A"  ): 
        """
        Attach an end–group (G) to the backbone (self) at a specified bond.
        
        The attachment is done in two steps:
        1. **Internal Orientation of the Group:**  
            The group is reoriented in its own frame by calling:
                G.orient(i0, (i0, i1), iup, _0=_0)
            - *i0*: the index (or indices) for the pivot atom in the group. This atom
                    is moved to the attachment position.
            - *(i0, i1)*: a tuple defining a bond in the group that determines the
                        forward (direction) vector. The forward vector is computed as
                        the difference between the positions of the atom at i1 and i0.
                        The atom corresponding to i1 is then deleted (replaced) in the group.
            - *iup*: a tuple (or list) of two indices that defines the up vector in the group.
                    The up vector is computed (typically as the difference between the
                    positions of the atoms provided) and is used to fix the rotation about
                    the forward axis.
        
        2. **Alignment to the Backbone:**  
            The backbone provides the attachment bond (given by `bond`) and a backbone
            up vector (given by `up`). A rotation matrix is computed with:
                rot = rotmat_from_points(self.apos, ifw=bond, up=up, _0=_0)
            This matrix aligns the backbone’s forward vector (computed from the bond) with
            the group’s forward vector. The group is then rotated by this matrix (via
                G.orient_mat(rot)
            ) and translated so that the pivot atom of the group coincides with the backbone’s
            attachment position.
        
        Parameters:
        G      : AtomicSystem
                The end–group to attach. It must have its atoms pre‐oriented as per the
                expected coordinate system.
        i0     : int or iterable
                The index (or indices) of the pivot atom in G (1-based indexing if _0=1).
        i1     : int
                The index (1-based) of the atom in G used to define the forward vector.
                This atom will be removed after orientation.
        iup    : tuple (i_up0, i_up1)
                A pair of indices (1-based) in G whose difference defines the up vector.
        bond   : tuple (i_backbone1, i_backbone2)
                A pair of atom indices (1-based) in the backbone that define the bond where
                the end–group is attached. The forward vector on the backbone is computed as
                the vector from i_backbone1 to i_backbone2.
        up     : 3-tuple or array, optional (default=(0.,0.,1.))
                The up vector for the backbone. This is used to fix the rotation about the
                forward axis.
        _0     : int, optional (default=1)
                An offset to account for whether the provided indices are 0-based or 1-based.
        pre    : str, optional (default="A")
                A prefix for labeling the atoms that come from the group.
        
        After executing, the group G is reoriented, rotated, and translated so that its
        pivot atom is placed at the backbone’s attachment site. The atom used for forward
        definition (i1) is deleted, and the group’s atoms are appended to the backbone.
        """
        up  = np.array( up )
        rot = au.rotmat_from_points( self.apos, ifw=bond, up=up, _0=1 );   
        rot = rot.transpose()
        p0  = self.apos[bond[0]-_0]
        
        if( G.aux_labels is None ): G.aux_labels = [ pre+str(i) for  i in range(len(G.apos)) ]

        G.orient( i0,(i0,i1),iup, _0=_0 )
        G.orient_mat( rot ); 
        G.apos[:,:]+=p0[None,:]
        G.delete_atoms( [i1-_0] )

        self.append_atoms( G )

    def reindex_bonds(self, old_to_new_map, to_remove=None ):
        #print ("self.reindex_bonds: old_to_new_map \n", old_to_new_map)
        #print ("self.reindex_bonds: to_remove      \n", to_remove)
        self.bonds = au.reindex_bonds( self.bonds, old_to_new_map, to_remove )
        self.ngs   = None

    def extract_marker_pairs(self, markerX, markerY, remove=True):
        """Legacy method that combines finding marker pairs without removal."""
        pairs = self.find_marker_pairs(markerX, markerY)
        return pairs

    def find_marker_pairs(self, markerX, markerY):
        """
        Find marker pairs in this system based on element types and bonding information.
        For each atom with element equal to markerX, look for a bonded neighbor with element equal to markerY.
        Returns a list of tuples (iX, iY) where iX is the index of a markerX atom and iY is the index of its bonded markerY neighbor.
        """
        if self.ngs is None:
            self.neighs(bBond=True)
        mks = []
        for i, ename in enumerate(self.enames):
            if ename == markerX:    # pivot atom marker-X
                ngi = self.ngs[i]    
                for j in ngi:
                    if self.enames[j] == markerY:
                        i2=j   # index of a markerY-typed bonded to the pivot atom
                    else:
                        i3=j   # atom to which pivot atom is bonded, but is not a marker (i.e. Anchor-atom)  
                mks.append( (i, i2, i3))                
        return mks

    def ensure_numpy_arrays(self):
        """Ensure position arrays are numpy arrays."""
        if not isinstance(self.apos, np.ndarray):
            self.apos = np.array(self.apos)
        
    def filter_system(self, mask, bInverted=False):
        """Create a new filtered system without marker atoms.
        Parameters:
            mask : set Set of atom indices to keep ( or remove if inverted)
            bInverted : bool Invert the mask
        Returns:
            tuple : (filtered arrays, old_to_new index_map)
        """
        n = len(self.apos)
        #if bInverted: mask = set(range(len(self.apos))).difference(mask)
        if bInverted: mask = [i for i in range(n) if i not in mask]
        old_to_new = au.make_reindex( n, mask, bInverted=False)    
        if self.bonds is not  None:
            bonds = [ b for b in self.bonds if (b[0] in mask) and (b[1] in mask) ]
            bonds = [ (old_to_new[b[0]], old_to_new[b[1]]) for b in bonds  ]
        else:
            bonds = None
        filtered = AtomicSystem(
            apos  =self.apos  [mask], 
            atypes=self.atypes[mask],
            enames=self.enames[mask],
            lvec  =self.lvec,
            qs    =self.qs[mask] if self.qs is not None else None,
            Rs    =self.Rs[mask] if self.Rs is not None else None,
            #ngs   =self.ngs[mask] if self.ngs is not None else None,
            bonds = bonds
        )
          
        return filtered, old_to_new
        
    #def merge_arrays(self, other, other_bonds, offset):
    def merge_arrays(self, other, offset=None):
        """Merge arrays from other system into self.
        
        Parameters:
            other_arrays : dict
                Filtered arrays from other system
            other_bonds : ndarray
                Reindexed bonds from other system
            offset : int
                Offset for bond indices
        """
        if offset is None: offset = len(self.apos)

        self.apos   = np.concatenate([self.apos,   other.apos], axis=0)
        self.atypes = np.concatenate([self.atypes, other.atypes])
        self.enames = np.concatenate([self.enames, other.enames])
        if self.qs         is not None and other.qs         is not None: self.qs         = np.concatenate([self.qs, other.qs])
        if self.Rs         is not None and other.Rs         is not None: self.Rs         = np.concatenate([self.Rs, other.Rs])
        if self.aux_labels is not None and other.aux_labels is not None: 
            self.aux_labels = np.concatenate([self.aux_labels, other.aux_labels])
        else:
            self.aux_labels = None
            
        # Merge bonds
        if other.bonds is not None:
            adjusted_bonds = np.array(other.bonds) + offset
            if self.bonds is not None:
                self.bonds = np.array(self.bonds)
                self.bonds = np.concatenate([self.bonds, adjusted_bonds], axis=0)
            else:
                self.bonds = adjusted_bonds

    def add_bond(self, b):
        #print ("add_bond", b)
        #print( "bonds.shape", self.bonds.shape )
        self.bonds = np.concatenate((self.bonds, np.array([b])), axis=0)

    def merge_geometries(self, other, group_mk, backbone_mk ):
        """Merge another AtomicSystem into this one using the provided group marker pair for alignment.
        
        This implementation appends the atoms and bonds from 'other' into self, adjusting indices appropriately.
        The process follows these steps:
        1. Find neighbors of marker atoms in both systems
        2. Remove marker atoms from the group system
        3. Merge the remaining atoms and bonds
        4. Create new bonds between fragments based on marker neighbors
        
        Parameters:
            other : AtomicSystem
                   The system to merge into this one
            group_marker_pair : tuple
                   (iX, iY) marker pair from the group system used for alignment
        """
        # Ensure numpy arrays
        self.ensure_numpy_arrays()
        other.ensure_numpy_arrays()
        
        removed = set(group_mk[:2])
        other_filtered, old_to_new = other.filter_system( removed, bInverted=True            )     # Filter group system without markers 
        
        # Merge arrays with offset
        offset = len(self.apos)
        #self.merge_arrays(other_filtered, other_bonds, offset)
        self.merge_arrays(other_filtered, offset)

        #print( "old_to_new", old_to_new )
        i2 = old_to_new[group_mk[2]]

        #print( "-----BEFORE self.bonds", self.bonds )
        self.add_bond( (backbone_mk[2], i2 + offset) )
        #print( "-----AFTER self.bonds", self.bonds )

        # Clear neighbor list since it needs to be rebuilt
        self.ngs = None

    def compute_group_orientation(self, G, backbone_pair, group_pair, _0=1):
        """Compute the orientation transformation for attaching a group to the backbone.
        
        Parameters:
            G : AtomicSystem
                The group to orient
            backbone_pair : tuple
                (iX, iY) indices of marker atoms in backbone
            group_pair : tuple
                (iX, iY) indices of marker atoms in group
            _0 : int
                Offset for index conversion
                
        Returns:
            tuple: (R, X_b, A2)
                R - rotation matrix
                X_b - translation point (backbone marker position)
                A2 - attachment point in group (for translation)
        """
        iX_b, iY_b,_  = backbone_pair
        iX_g, iY_g,_ = group_pair
        
        # Compute frames for backbone and group
        X_b, A1, M_target = au.compute_attachment_frame_from_indices(self.apos, iX_b, iY_b, self, bFlipFw=False, _0=_0)
        X_g, A2, M_group = au.compute_attachment_frame_from_indices(G.apos, iX_g, iY_g, G, bFlipFw=True, _0=_0)
        
        # Compute rotation matrix R = M_target @ (M_group)ᵀ
        R = M_target @ M_group.T
        
        return R, X_b, A2

    def delete_atoms(self, to_remove):
            rem = sorted(to_remove, reverse=True)
            for idx in rem:
                self.apos   = np.delete(self.apos,   idx, axis=0)
                self.atypes = np.delete(self.atypes, idx)
                self.enames = np.delete(self.enames, idx)
                if self.qs is not None:
                    self.qs = np.delete(self.qs, idx)
                if self.Rs is not None:
                    self.Rs = np.delete(self.Rs, idx)
                if self.aux_labels is not None:
                    self.aux_labels = np.delete(self.aux_labels, idx)

    def attach_group_by_marker(self, G, markerX="Xe", markerY="He", _0=1):
        """Attach an end–group G to this backbone using marker atoms and connectivity.
        Steps:
          1. Find marker pairs in both the backbone and the group.
          2. Ensure the group has exactly one marker pair.
          3. Orient and transform the group.
          4. Merge the transformed group.
          5. Remove marker atoms from the backbone.
          6. Update the neighbor list.
        """
        # 1. Find marker pairs
        backbone_inds = self.find_marker_pairs(markerX, markerY)
        group_inds    = G.find_marker_pairs(markerX, markerY)
        #print( "backbone_inds ", backbone_inds )
        #print( "group_inds    ", group_inds )
        if not backbone_inds:    raise ValueError(f"No marker pair ({markerX}, {markerY}) found in backbone")
        if len(group_inds) != 1: raise ValueError(f"Group must have exactly one marker pair, found {len(group_inds)}")
            
        # 2. Get orientation transformation
        R, X_b, A2 = self.compute_group_orientation(G, backbone_inds[0], group_inds[0], _0)
        
        # 3. Make a deep copy of G and transform it
        G_copy = copy.deepcopy(G)
        G_copy.apos = (R @ (G_copy.apos - A2).T).T + X_b
        
        gind = group_inds[0]
        bind = backbone_inds[0]

        for ii, bind in enumerate(backbone_inds):

            bind = self.find_marker_pairs(markerX, markerY)[0]

            # 2. Get orientation transformation
            R, X_b, A2 = self.compute_group_orientation(G, bind, gind, _0)
            
            # 3. Make a deep copy of G and transform it
            G_copy = copy.deepcopy(G)
            G_copy.apos = (R @ (G_copy.apos - A2).T).T + X_b

            self.merge_geometries(G_copy, gind, bind )
            to_remove = set( [bind[0], bind[1]] )

            old_to_new = {}
            new_idx = 0
            for old_idx in range(len(self.apos)):
                if old_idx not in to_remove:
                    old_to_new[old_idx] = new_idx
                    new_idx += 1
            self.reindex_bonds(old_to_new, to_remove)

            self.delete_atoms( to_remove )

        # 6. Update neighbor list
        self.neighs(bBond=True)

# ========= Adding Electron Pairs

    def add_electron_pairs(self):
        """
        Add electron pairs to atoms (N, O) based on their chemical neighborhood.
        """
        for i, ename in enumerate(self.enames):
            if ename not in VALENCE_DICT:
                continue
            nb     = VALENCE_DICT[ename][0]
            nep    = VALENCE_DICT[ename][1]
            nsigma = len(self.ngs[i])
            npi    = nb - nsigma
            #print( "Atom %i: %s, npi = %i, nsigma = %i, nep = %i" % (i, ename, npi, nsigma, nep) )
            if nep > 0:
                self.make_epair_geom(i, npi, nsigma)

    def get_atomi_pi_direction(self, i):
        """
        Get the pi-direction for atom i.
        """
        # we should go over 2-3 neighbors (depending how many there is) and compute cross product, take average and normalize
        #if self.bonds is None or i not in self.ngs or len(self.ngs[i]) < 2:
        #    return np.array([0.0, 0.0, 1.0])  # Default direction if not enough neighbors
        neighbors = self.ngs[i]  # Take up to 3 neighbors
        vectors = [ au.normalize(self.apos[j] - self.apos[i]) for j in neighbors if j != -1]
        dir = np.zeros(3)
        for a, b in zip(vectors, vectors[1:] + [vectors[0]]):
            dir += au.normalize(np.cross(a, b))
        return au.normalize(dir)

    def make_epair_geom(self, i, npi, nb ):
        """
        Add electron pairs to atom i based on configuration using vector operations.
        """
        pos = self.apos[i]
        # Get neighbor positions as list of numpy arrays
        neighbors = [self.apos[j] for j in self.ngs[i]]
        if nb > 0: v1 = au.normalize( neighbors[0] - pos )
        if nb > 1: v2 = au.normalize( neighbors[1] - pos )
        if nb > 2: v3 = au.normalize( neighbors[2] - pos )
        #print( f"make_epair_geom() ia: {i} npi: {npi} nb: {nb}" )
        if npi == 0:
            if nb == 3:   # like NH3
                #print( f"make_epair_geom() like NH3 {self.enames[i]}    ia: {i} npi: {npi} nb: {nb}" )
                base = np.cross(v2 - v1, v3 - v1)
                base = au.normalize(base)
                if np.dot(base, (v1 + v2 + v3)) > 0:  base = -base
                self.place_electron_pair(i, base)
            elif nb == 2: # like H2O
                #print( f"make_epair_geom() like H2O {self.enames[i]}    ia: {i} npi: {npi} nb: {nb}" )
                m_c = au.normalize(v1 + v2)  # Average (bisector) direction
                m_b = np.cross( v1, v2 )
                m_b = au.normalize(m_b)
                cc = 0.57735026919  # sqrt(1/3)
                cb = 0.81649658092  # sqrt(2/3)
                ep1 = au.normalize(m_c * -cc + m_b * cb)
                ep2 = au.normalize(m_c * -cc - m_b * cb)
                self.place_electron_pair(i, ep1)
                self.place_electron_pair(i, ep2)
        elif npi == 1:
            #print( "make_epair_geom() PI=1 ia: %i npi: %i nb: %i" % (i, npi, nb ) )
            if nb == 2: # like =N-
                #print( f"make_epair_geom() like =N- {self.enames[i]}    ia: {i} npi: {npi} nb: {nb}" )
                m_c = au.normalize(v1 + v2)  # Bisector
                self.place_electron_pair(i, m_c*-1.)
            elif nb == 1:  # like =O 
                #print( f"make_epair_geom() like =O {self.enames[i]}    ia: {i} npi: {npi} nb: {nb}" )
                m_b = self.get_atomi_pi_direction( self.ngs[i][0] )
                m_c = au.normalize(np.cross( v1, m_b ))
                self.place_electron_pair(i, v1*-0.5 + m_c*0.86602540378)
                self.place_electron_pair(i, v1*-0.5 - m_c*0.86602540378)
        # elif npi == 2:
        #     if nb == 1:
        #         m_c = normalize(neighbors[0] - pos)
        #         self._place_electron_pair(i, -m_c)

    def place_electron_pair(self, i, direction, distance=0.5, ename='E', atype=200, qs=0.0, Rs=1.0):
        """
        Place an electron pair in the specified direction from atom i.
        Adds new atom with ename='E', atype=200 to apos, atypes, enames arrays.
        """
        # Calculate electron pair position
        ep_pos = self.apos[i] + direction * distance

        # Append to arrays
        self.apos = np.append(self.apos, [ep_pos], axis=0)
        self.atypes = np.append(self.atypes, atype)
        self.enames.append(ename)

        # Initialize charge to zero
        if self.qs is not None:
            self.qs = np.append(self.qs, qs)

        # Initialize Rs with default value
        if self.Rs is not None:
            self.Rs = np.append(self.Rs, Rs)  # Default radius for electron pairs

        # Update bonds and neighbors if needed
        if self.bonds is not None:
            self.bonds = np.append(self.bonds, [[i, len(self.apos) - 1]], axis=0)
        if self.ngs is not None:
            self.ngs[i][len(self.apos) - 1] = 1  # Add to neighbors
            self.ngs.append({i: 1})  # Add new neighbor list for electron pair
```