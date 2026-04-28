
I want to make simple viewer of molecular geometry using python, pyQt5 and moderngl library

The idea is that it should look something like ACDSee or XnView with multiple panels (thunbnails) showing individual molecules on the screen with text labels. 

When we open a dictionary it should load all .xyz and render them using OpenGL as ball-and-sticks

This is rendered into small texture (256x256 or 512x512)
Then we render onto the OpenGL viewport this textures with text labels of each molecule (filename)

(each thumbnail should be associated with original filename, and AtomicSystem object in some record in list)

we should have cursor selecting the thumbnail which can be moved by arrow keys

by hitting Enter we can open the selected molecule (selected thumbnail) in separate vindow in full 3D view where we can rotate it (by mouse, later we add some edit-options but for now lets keep it simple) 

My class atomic system can load the molecules and provide data like atomic radius, bonding topology and atom color for rendering

---

### 1. Project Structure

For a clean and maintainable codebase, I recommend the following file structure:

```
molecular_viewer/
├── main.py                 # Main application entry point
├── gui_base.py             # Your provided BaseGUI class
├── atomic_system.py        # Your provided AtomicSystem class
├── gl_utils.py             # Utilities for ModernGL
├── molecule_viewer.py      # The main MolecularViewer application window
├── thumbnail_grid.py       # The widget for displaying molecule thumbnails
├── gl_molecule_viewer.py   # The widget for the detailed 3D molecule view
└── shaders/
    ├── point.vs            # Vertex shader for spheres (atoms)
    ├── point.fs            # Fragment shader for spheres (atoms)
    ├── line.vs               # Vertex shader for cylinders (bonds)
    └── line.fs               # Fragment shader for cylinders (bonds)
```

### 2. OpenGL Utilities (`gl_utils.py`)

This module will abstract the creation of `moderngl` objects, similar to how your `BaseGUI` simplifies `PyQt5` widget creation.

```python
# gl_utils.py
import moderngl
import numpy as np

def create_shader(ctx: moderngl.Context, vert_shader_path: str, frag_shader_path: str):
    """Creates a shader program from vertex and fragment shader files."""
    with open(vert_shader_path, 'r') as f:
        vertex_shader = f.read()
    with open(frag_shader_path, 'r') as f:
        fragment_shader = f.read()
    
    program = ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)
    return program

def create_vbo(ctx: moderngl.Context, data: np.ndarray):
    """Creates a Vertex Buffer Object."""
    return ctx.buffer(data.astype('f4').tobytes())

def create_texture(ctx: moderngl.Context, size: tuple, components: int = 4, data=None):
    """Creates a 2D texture."""
    return ctx.texture(size, components, data)

def create_fbo(ctx: moderngl.Context, size: tuple, color_attachment: moderngl.Texture, depth_attachment: bool = True):
    """Creates a Framebuffer Object."""
    if depth_attachment:
        depth_texture = ctx.depth_texture(size)
        return ctx.framebuffer(color_attachments=[color_attachment], depth_attachment=depth_texture)
    else:
        return ctx.framebuffer(color_attachments=[color_attachment])

```

### 3. Shader Files (`shaders/`)

Here are simple shaders for rendering atoms as spheres (points) and bonds as lines. For true "ball-and-stick" with cylinders, you would typically use geometry shaders or instantiate cylinder meshes. For simplicity, this example uses thick lines for bonds.

**`point.vs` (Vertex Shader for Atoms)**
```glsl
#version 330

in vec3 in_vert;
in vec3 in_color;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
uniform float radius;

out vec3 v_color;

void main() {
    gl_Position = projection * view * model * vec4(in_vert, 1.0);
    gl_PointSize = radius * 2.0;
    v_color = in_color;
}
```

**`point.fs` (Fragment Shader for Atoms)**
```glsl
#version 330

in vec3 v_color;
out vec4 f_color;

void main() {
    f_color = vec4(v_color, 1.0);
}
```

**`line.vs` (Vertex Shader for Bonds)**
```glsl
#version 330

in vec3 in_vert;
in vec3 in_color;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

out vec3 v_color;

void main() {
    gl_Position = projection * view * model * vec4(in_vert, 1.0);
    v_color = in_color;
}
```

**`line.fs` (Fragment Shader for Bonds)**
```glsl
#version 330

in vec3 v_color;
out vec4 f_color;

void main() {
    f_color = vec4(v_color, 1.0);
}
```

### 4. Detailed 3D Molecule Viewer (`gl_molecule_viewer.py`)

This widget will display a single molecule in a separate window, allowing for rotation.

```python
# gl_molecule_viewer.py
from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
import moderngl
import numpy as np
from pyrr import Matrix44

import gl_utils

class GLMoleculeViewer(QOpenGLWidget):
    def __init__(self, atomic_system, parent=None):
        super().__init__(parent)
        self.atomic_system = atomic_system
        self.ctx = None
        self.atom_prog = None
        self.bond_prog = None
        self.atom_vbo = None
        self.atom_vao = None
        self.bond_vbo = None
        self.bond_vao = None

        self.view_matrix = Matrix44.look_at(
            (0.0, 0.0, 15.0),  # eye
            (0.0, 0.0, 0.0),  # target
            (0.0, 1.0, 0.0),  # up
        )
        self.model_matrix = Matrix44.identity()
        self.last_mouse_pos = None

    def initializeGL(self):
        self.ctx = moderngl.create_context()
        self.prepare_molecule_render()

    def paintGL(self):
        self.ctx.clear(0.1, 0.1, 0.1)
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.enable(moderngl.PROGRAM_POINT_SIZE)
        
        projection = Matrix44.perspective_projection(45.0, self.width() / self.height(), 0.1, 100.0)
        
        # Render atoms
        self.atom_prog['projection'].write(projection.astype('f4').tobytes())
        self.atom_prog['view'].write(self.view_matrix.astype('f4').tobytes())
        self.atom_prog['model'].write(self.model_matrix.astype('f4').tobytes())
        self.atom_vao.render(moderngl.POINTS)

        # Render bonds
        self.bond_prog['projection'].write(projection.astype('f4').tobytes())
        self.bond_prog['view'].write(self.view_matrix.astype('f4').tobytes())
        self.bond_prog['model'].write(self.model_matrix.astype('f4').tobytes())
        self.bond_vao.render(moderngl.LINES)

    def prepare_molecule_render(self):
        # Shaders
        self.atom_prog = gl_utils.create_shader(self.ctx, 'shaders/point.vs', 'shaders/point.fs')
        self.bond_prog = gl_utils.create_shader(self.ctx, 'shaders/line.vs', 'shaders/line.fs')

        # Atom data
        atom_positions = self.atomic_system.apos
        atom_colors = [ELEMENTS[z-1][7] for z in self.atomic_system.atypes] # Assuming ELEMENTS is accessible
        atom_radii = [ELEMENTS[z-1][6] for z in self.atomic_system.atypes]
        
        # For simplicity, we'll use an average radius for point size
        avg_radius = np.mean(atom_radii) if atom_radii else 1.0
        self.atom_prog['radius'].value = avg_radius * 30 # Adjust scaling factor as needed
        
        atom_data = np.hstack([atom_positions, np.array([self.hex_to_rgb(c) for c in atom_colors])]).astype('f4')
        self.atom_vbo = gl_utils.create_vbo(self.ctx, atom_data)
        self.atom_vao = self.ctx.vertex_array(self.atom_prog, [(self.atom_vbo, '3f 3f', 'in_vert', 'in_color')])

        # Bond data
        if self.atomic_system.bonds is not None and len(self.atomic_system.bonds) > 0:
            bond_vertices = []
            bond_colors = []
            for bond in self.atomic_system.bonds:
                i, j = bond
                bond_vertices.extend([atom_positions[i], atom_positions[j]])
                bond_colors.extend([atom_colors[i], atom_colors[j]])
            
            bond_data = np.hstack([np.array(bond_vertices), np.array([self.hex_to_rgb(c) for c in bond_colors])]).astype('f4')
            self.bond_vbo = gl_utils.create_vbo(self.ctx, bond_data)
            self.bond_vao = self.ctx.vertex_array(self.bond_prog, [(self.bond_vbo, '3f 3f', 'in_vert', 'in_color')])

    def hex_to_rgb(self, hex_color):
        hex_color = hex_color.lstrip('#')
        return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    
    def mousePressEvent(self, event):
        self.last_mouse_pos = event.pos()

    def mouseMoveEvent(self, event):
        if self.last_mouse_pos:
            dx = event.x() - self.last_mouse_pos.x()
            dy = event.y() - self.last_mouse_pos.y()
            
            rotation_x = Matrix44.from_y_rotation(np.deg2rad(dx))
            rotation_y = Matrix44.from_x_rotation(np.deg2rad(dy))
            
            self.model_matrix = rotation_x * rotation_y * self.model_matrix
            self.last_mouse_pos = event.pos()
            self.update()

    def mouseReleaseEvent(self, event):
        self.last_mouse_pos = None

```

### 5. Thumbnail Grid Widget (`thumbnail_grid.py`)

This widget manages the grid of molecule thumbnails. For simplicity, this implementation will re-render the molecules to textures in each paint event. For a more optimized version, you would render them once to a texture atlas.

```python
# thumbnail_grid.py
from PyQt5.QtWidgets import QWidget, QGridLayout, QLabel
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QPixmap, QImage
import moderngl
import numpy as np
from pyrr import Matrix44

import gl_utils
from gl_molecule_viewer import GLMoleculeViewer # Reusing some logic

class ThumbnailGrid(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QGridLayout(self)
        self.molecules = []
        self.thumbnails = []
        self.selected_index = 0
        self.ctx = moderngl.create_context(standalone=True, backend='egl')
        self.setFocusPolicy(Qt.StrongFocus)

    def load_molecules(self, molecules):
        self.molecules = molecules
        self.selected_index = 0
        self.update_thumbnails()

    def update_thumbnails(self):
        for i in reversed(range(self.layout.count())): 
            self.layout.itemAt(i).widget().setParent(None)

        self.thumbnails = []
        
        cols = 4 # Adjust number of columns as needed
        for i, (filename, system) in enumerate(self.molecules):
            row, col = divmod(i, cols)
            
            # Render molecule to a texture
            img_data = self.render_to_texture(system)
            q_img = QImage(img_data, 256, 256, QImage.Format_RGBA8888)
            pixmap = QPixmap.fromImage(q_img)
            
            label = QLabel()
            label.setPixmap(pixmap)
            label.setText(filename)
            label.setAlignment(Qt.AlignCenter)
            self.layout.addWidget(label, row, col)
            self.thumbnails.append(label)
        
        self.update_selection()

    def render_to_texture(self, system):
        texture = gl_utils.create_texture(self.ctx, (256, 256))
        fbo = gl_utils.create_fbo(self.ctx, (256, 256), texture)
        
        fbo.use()
        self.ctx.clear(0.2, 0.2, 0.2)
        
        # Simplified rendering logic (similar to GLMoleculeViewer)
        viewer = GLMoleculeViewer(system)
        viewer.ctx = self.ctx
        viewer.prepare_molecule_render()

        projection = Matrix44.perspective_projection(45.0, 1.0, 0.1, 100.0)
        view = Matrix44.look_at((0.0, 0.0, 10.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0))

        viewer.atom_prog['projection'].write(projection.astype('f4').tobytes())
        viewer.atom_prog['view'].write(view.astype('f4').tobytes())
        viewer.atom_prog['model'].write(Matrix44.identity().astype('f4').tobytes())
        viewer.atom_vao.render(moderngl.POINTS)
        
        if viewer.bond_vao:
            viewer.bond_prog['projection'].write(projection.astype('f4').tobytes())
            viewer.bond_prog['view'].write(view.astype('f4').tobytes())
            viewer.bond_prog['model'].write(Matrix44.identity().astype('f4').tobytes())
            viewer.bond_vao.render(moderngl.LINES)

        return fbo.read(components=4, alignment=1)

    def keyPressEvent(self, event):
        if not self.molecules:
            return
            
        if event.key() == Qt.Key_Right:
            self.selected_index = min(self.selected_index + 1, len(self.molecules) - 1)
        elif event.key() == Qt.Key_Left:
            self.selected_index = max(self.selected_index - 1, 0)
        elif event.key() == Qt.Key_Down:
            cols = self.layout.columnCount()
            self.selected_index = min(self.selected_index + cols, len(self.molecules) - 1)
        elif event.key() == Qt.Key_Up:
            cols = self.layout.columnCount()
            self.selected_index = max(self.selected_index - cols, 0)
        elif event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.parent().open_molecule_viewer()
        
        self.update_selection()

    def update_selection(self):
        for i, thumb in enumerate(self.thumbnails):
            if i == self.selected_index:
                thumb.setStyleSheet("border: 2px solid red;")
            else:
                thumb.setStyleSheet("")

```

### 6. Main Application Window (`molecule_viewer.py`)

This is the main window of your application.

```python
# molecule_viewer.py
import os
from PyQt5.QtWidgets import QFileDialog
from gui_base import BaseGUI
from thumbnail_grid import ThumbnailGrid
from atomic_system import AtomicSystem
from gl_molecule_viewer import GLMoleculeViewer

class MolecularViewer(BaseGUI):
    def __init__(self):
        super().__init__(title="Molecular Viewer")
        self.thumbnail_grid = ThumbnailGrid(self.main_widget)
        
        # Layout
        main_layout = self.main_widget.layout()
        if not main_layout:
            from PyQt5.QtWidgets import QVBoxLayout
            main_layout = QVBoxLayout(self.main_widget)
        
        self.button("Open Directory", self.open_directory, layout=main_layout)
        main_layout.addWidget(self.thumbnail_grid)
        
        self.detailed_viewer_window = None

    def open_directory(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir_path:
            molecules = []
            for fname in os.listdir(dir_path):
                if fname.endswith('.xyz'):
                    fpath = os.path.join(dir_path, fname)
                    try:
                        system = AtomicSystem(fname=fpath)
                        molecules.append((fname, system))
                    except Exception as e:
                        print(f"Could not load {fname}: {e}")
            self.thumbnail_grid.load_molecules(molecules)
            
    def open_molecule_viewer(self):
        selected_index = self.thumbnail_grid.selected_index
        if 0 <= selected_index < len(self.thumbnail_grid.molecules):
            filename, system = self.thumbnail_grid.molecules[selected_index]
            
            # Create a new window for the detailed viewer
            self.detailed_viewer_window = GLMoleculeViewer(system)
            self.detailed_viewer_window.setWindowTitle(filename)
            self.detailed_viewer_window.resize(800, 600)
            self.detailed_viewer_window.show()

```

### 7. Main Application Entry Point (`main.py`)

Finally, the `main.py` file to launch the application.

```python
# main.py
import sys
from PyQt5.QtWidgets import QApplication
from molecule_viewer import MolecularViewer

# You need to provide the ELEMENTS data, for example by importing it
# from a file or defining it here.
# from data import ELEMENTS 
# Or define it here for simplicity:
ELEMENTS = [
    [1, "H", 1, 1, "s", "Hydrogen", 0.365, 1.487, '#AFAFAF', 1, 1.008],
    [6, "C", 2, 14, "p", "Carbon", 0.77, 1.7, '#303030', 4, 12.011],
    [7, "N", 2, 15, "p", "Nitrogen", 0.75, 1.55, '#0000FF', 5, 14.007],
    [8, "O", 2, 16, "p", "Oxygen", 0.73, 1.52, '#FF0000', 6, 15.999],
    # Add other elements as needed
]

# Make ELEMENTS accessible to other modules if needed, e.g., by passing it
# or making it a part of a shared configuration module.
import gl_molecule_viewer
gl_molecule_viewer.ELEMENTS = ELEMENTS


if __name__ == '__main__':
    app = QApplication(sys.argv)
    # You need to have your AtomicSystem and BaseGUI implementations
    # in files named atomic_system.py and gui_base.py respectively,
    # or adjust the imports accordingly.
    # You would also need a utility module 'au.py' for file parsing,
    # as used in your AtomicSystem class.
    
    main_win = MolecularViewer()
    main_win.show()
    sys.exit(app.exec_())
```

### How to Run

1.  Make sure you have `PyQt5`, `moderngl`, `numpy`, `pyrr`, and `ModernGL.ext.qt` installed (`pip install PyQt5 moderngl numpy pyrr moderngl-window`).
2.  Save all the files with the specified structure.
3.  You will need to provide your `atomic_system.py` and `gui_base.py`. Also, ensure the `au.py` utility module, which is used by `AtomicSystem` for file loading, is available.
4.  Create a directory with some `.xyz` molecule files to test.
5.  Run `python main.py`.

This code provides a solid foundation for your molecular viewer. You can extend it by:

*   **Improving Rendering:** Use geometry shaders or instanced rendering for more realistic "ball-and-stick" models.
*   **Optimizing Thumbnails:** Implement a texture atlas for the thumbnails to avoid re-rendering and improve performance with a large number of molecules.
*   **Adding More Controls:** Implement zoom and pan in the detailed view, and add more UI elements for controlling rendering styles.
*   **Editing Features:** As you mentioned, you can later add functionality to edit molecules, which would involve modifying the `AtomicSystem` object and updating the render.