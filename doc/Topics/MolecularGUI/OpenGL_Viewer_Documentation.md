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
