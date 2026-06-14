I was testing this approach, but it was not really working, I had problem wih OpenGL context distributed among the browser and viewer. But also the whole architecture was too much different from my original idea.

My idea is to have sut one window wih one OpenGL viewport. This viewport will switch between thumbnail mode and molecule detail mode when user hit enter.
We do not wan to use any Qt wiggets to denter the thumbnails, we want to use pure OpenGL (like sprites with texture), we also wan to render text labels (file names) under these sprites inside opengl (I hope you know how to render bitmap text in opengl easily) it doest not have to be super prity, just something simple. Then there should be the cursor which will be e.g. just crean frame rendered using 4 lines in opengl. The movement of the cursor and layout of the thumbnails (sprites) on the grid inside viewpoet, and eventually rollong (if there is more thumbnails which does not fit to the screen) you need to handle youself, not rely ond PyQt


---

You are absolutely right. My apologies for diverging from your original, more elegant vision. A single-window, single-context approach that switches rendering modes is indeed more efficient, robust, and stays true to the spirit of your `GLGUI` framework.

Let's build this exactly as you've described: one OpenGL viewport for everything, with rendering done as sprites, OpenGL-rendered text, and a manually controlled cursor and layout.

### Architecture Overview

1.  **Single Window/Widget**: We will have one main class, `MoleculeBrowserGLWidget`, that inherits from your `BaseGLWidget`. This widget will handle everything.
2.  **State Machine**: The widget will have a `view_mode` attribute that can be either `MODE_THUMBNAIL` or `MODE_DETAIL`.
3.  **Rendering Logic**: The `draw_scene` method will check the `view_mode` and call either `draw_thumbnail_view()` or `draw_detail_view()`.
4.  **Thumbnail Generation**: When a directory is opened, we will use a Framebuffer Object (FBO) within the *same OpenGL context* to render each molecule into a texture *once*. These textures will be stored.
5.  **Thumbnail Drawing**: In `draw_thumbnail_view`, we will draw a series of textured quads (sprites), using the textures generated in the previous step.
6.  **Text Rendering**: We will implement a simple bitmap font renderer to draw the filenames under the sprites directly in OpenGL.
7.  **Layout & Control**: All thumbnail positions, grid layout, cursor movement, and scrolling will be calculated and handled manually within `MoleculeBrowserGLWidget`.

### Step 1: Utility for Bitmap Font

First, we need a way to create a font texture atlas. This is a one-time process. We can use `Pillow` to generate an image with all the characters we need.

**Create a new file `make_font.py`:**

```python
# make_font.py
from PIL import Image, ImageDraw, ImageFont
import json

# Configuration
FONT_FILE = "C:/Windows/Fonts/consola.ttf" # Or a path to any other monospace font
FONT_SIZE = 32
CHARS = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
IMG_WIDTH = 512
PADDING = 2

def generate_font_atlas():
    """Generates a font atlas image and a JSON file with character metrics."""
    try:
        font = ImageFont.truetype(FONT_FILE, FONT_SIZE)
    except IOError:
        print(f"Font file not found at {FONT_FILE}, using default font.")
        font = ImageFont.load_default()

    img = Image.new("RGBA", (IMG_WIDTH, 1024), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    char_data = {}
    x, y, max_h = PADDING, PADDING, 0

    for char in CHARS:
        # Use getbbox for modern Pillow versions
        bbox = draw.textbbox((0,0), char, font=font)
        char_width = bbox[2] - bbox[0]
        char_height = bbox[3] - bbox[1]

        if x + char_width + PADDING > IMG_WIDTH:
            x = PADDING
            y += max_h + PADDING
            max_h = 0
        
        if char_height > max_h:
            max_h = char_height

        draw.text((x, y), char, font=font, fill=(255, 255, 255, 255))
        
        char_data[char] = {
            'x': x, 'y': y,
            'w': char_width, 'h': char_height,
            'norm_x': x / IMG_WIDTH, 'norm_y': y / 1024,
            'norm_w': char_width / IMG_WIDTH, 'norm_h': char_height / 1024,
        }
        x += char_width + PADDING

    # Crop image to actual height
    final_img = img.crop((0, 0, IMG_WIDTH, y + max_h + PADDING))
    final_img.save("font_atlas.png")

    # Save metrics
    with open("font_atlas.json", "w") as f:
        json.dump({'chars': char_data, 'tex_w': final_img.width, 'tex_h': final_img.height}, f)
    
    print("Font atlas 'font_atlas.png' and 'font_atlas.json' created successfully.")

if __name__ == "__main__":
    generate_font_atlas()
```

**Run this script once (`python make_font.py`)**. It will create `font_atlas.png` and `font_atlas.json` which we will use in our main application.

### Step 2: New Shaders

We need shaders for sprites, text, the cursor, and the ball-and-stick view.

**`shaders/sprite.glslv` (Vertex)**
```glsl
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;

uniform mat4 projection; // Orthographic projection
uniform vec2 u_pos;
uniform vec2 u_scale;

out vec2 TexCoord;

void main() {
    vec2 pos = aPos * u_scale + u_pos;
    gl_Position = projection * vec4(pos.x, pos.y, 0.0, 1.0);
    TexCoord = aTexCoord;
}
```

**`shaders/sprite.glslf` (Fragment)**
```glsl
#version 330 core
out vec4 FragColor;
in vec2 TexCoord;
uniform sampler2D u_texture;

void main() {
    FragColor = texture(u_texture, TexCoord);
}
```

**`shaders/text.glslv`** is identical to `sprite.glslv`.
**`shaders/text.glslf`** is identical to `sprite.glslf`.

**`shaders/flat_color.glslf` (For the cursor)**
```glsl
#version 330 core
out vec4 FragColor;
uniform vec4 u_color;

void main() {
    FragColor = u_color;
}
```

We will also reuse `instances.glslv`, `sphere.glslf`, and the new `cylinder.glslf` from the previous answer for the detailed view.

### Step 3: The Main Application (`MolecularBrowser.py`)

This file contains the core logic. It's a single class that manages state, rendering, and interaction.

```python
# MolecularBrowser.py
import sys
import os
import numpy as np
import json
from PIL import Image

from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QPushButton, QFileDialog)
from PyQt5.QtCore import Qt
from OpenGL.GL import *

# Assuming your project structure allows these imports
from GLGUI import BaseGLWidget, AppWindow, InstancedData, Mesh, octahedron_sphere_mesh
from AtomicSystem import AtomicSystem
from pyBall import elements
from pyBall.atomicUtils import makeRotMat


# --- Helper to create cylinder mesh (needed for ball-and-stick) ---
def create_cylinder_mesh(radius=1.0, length=1.0, segments=16):
    vertices, normals = [], []
    for i in range(segments + 1):
        angle = i * 2 * np.pi / segments
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        vertices.extend([x, y, -length/2]); normals.extend([x/radius, y/radius, 0])
        vertices.extend([x, y,  length/2]); normals.extend([x/radius, y/radius, 0])
    return np.array(vertices, dtype=np.float32), np.array(normals, dtype=np.float32)

# --- Main Application Widget ---
MODE_THUMBNAIL = 0
MODE_DETAIL = 1

class MoleculeBrowserGLWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.view_mode = MODE_THUMBNAIL
        self.molecules = []  # List of (filename, AtomicSystem, texture_id)
        self.selected_index = -1
        self.detail_system = None # Which system to show in detail view

        # Layout & Scrolling
        self.grid_cols = 4
        self.thumb_size = np.array([256, 256])
        self.spacing = np.array([40, 80])
        self.scroll_y = 0.0
        self.content_height = 0

        # Rendering Resources
        self.fbo = None
        self.render_buffer = None
        self.sprite_shader, self.text_shader, self.flat_shader = None, None, None
        self.atom_shader, self.bond_shader = None, None
        self.quad_mesh, self.cursor_mesh, self.text_mesh = None, None, None
        self.font_atlas_tex = None
        self.font_atlas_data = None
        
        self.atom_instances, self.bond_instances, self.cylinder_mesh = None, None, None

    @property
    def all_shader_programs(self):
        progs = [self.sprite_shader, self.text_shader, self.flat_shader, self.atom_shader, self.bond_shader]
        return [p for p in progs if p is not None]

    def initializeGL(self):
        # --- Load Shaders ---
        def load_shader(vert_path, frag_path):
            with open(vert_path) as f: v = f.read()
            with open(frag_path) as f: f = f.read()
            return self.compile_shader_program(v, f)

        shader_dir = "shaders/"
        self.sprite_shader = load_shader(shader_dir+"sprite.glslv", shader_dir+"sprite.glslf")
        self.text_shader   = load_shader(shader_dir+"text.glslv",   shader_dir+"text.glslf")
        self.flat_shader   = load_shader(shader_dir+"sprite.glslv", shader_dir+"flat_color.glslf") # Re-use vertex shader
        self.atom_shader   = load_shader(shader_dir+"instances.glslv", shader_dir+"sphere.glslf")
        self.bond_shader   = load_shader(shader_dir+"instances.glslv", shader_dir+"cylinder.glslf")

        super().initializeGL_base(None, None) # BaseGL init

        # --- FBO for Offscreen Rendering ---
        self.fbo = glGenFramebuffers(1)
        self.render_buffer = glGenRenderbuffers(1)

        # --- Shared Meshes ---
        # Sprite Quad (fullscreen, UVs flipped)
        quad_v = np.array([-1, -1, 1, -1, 1, 1, -1, 1], dtype=np.float32)
        quad_uv = np.array([0, 1, 1, 1, 1, 0, 0, 0], dtype=np.float32)
        self.quad_mesh = Mesh(vertices=quad_v.reshape((4,2)), normals=quad_uv) # Abuse 'normals' for UVs
        self.quad_mesh.setup_buffers() # VBO 0: pos, VBO 1: texcoord

        # Cursor Quad
        cursor_v = np.array([-1.05, -1.05, 1.05, -1.05, 1.05, 1.05, -1.05, 1.05], dtype=np.float32)
        self.cursor_mesh = Mesh(vertices=cursor_v.reshape((4,2)))
        self.cursor_mesh.setup_buffers()
        
        # --- Text Rendering Setup ---
        self.text_mesh = Mesh(vertices=np.zeros((4,2), dtype=np.float32), normals=np.zeros((4,2), dtype=np.float32))
        self.text_mesh.setup_buffers() # Will be updated dynamically
        
        self.font_atlas_tex = self.load_texture("font_atlas.png")
        with open("font_atlas.json") as f:
            self.font_atlas_data = json.load(f)

        # --- Detail View Rendering Setup ---
        self.atom_instances = InstancedData(base_attrib_location=2)
        self.atom_instances.associate_mesh(self.default_sphere_mesh)
        inst_attribs = [("positions",0,3),("radii",1,1),("colors",2,4),("instanceMatrix",3,16)]
        self.atom_instances.setup_instance_vbos(inst_attribs)

        cyl_v, cyl_n = create_cylinder_mesh(radius=0.1, length=1.0)
        self.cylinder_mesh = Mesh(vertices=cyl_v, normals=cyl_n)
        self.cylinder_mesh.setup_buffers()

        self.bond_instances = InstancedData(base_attrib_location=2)
        self.bond_instances.associate_mesh(self.cylinder_mesh)
        self.bond_instances.setup_instance_vbos(inst_attribs)

    def load_texture(self, filepath):
        img = Image.open(filepath).convert("RGBA")
        img_data = np.array(list(img.getdata()), np.uint8)
        
        tex_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex_id)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data)
        glBindTexture(GL_TEXTURE_2D, 0)
        return tex_id

    def load_directory(self, dir_path):
        # Cleanup old data
        for _, _, tex_id in self.molecules: glDeleteTextures(1, [tex_id])
        self.molecules.clear()

        filenames = sorted([f for f in os.listdir(dir_path) if f.lower().endswith('.xyz')])
        
        for fname in filenames:
            fpath = os.path.join(dir_path, fname)
            try:
                system = AtomicSystem(fname=fpath)
                tex_id = self.render_molecule_to_texture(system)
                self.molecules.append((fname, system, tex_id))
            except Exception as e:
                print(f"Error processing {fname}: {e}")

        self.selected_index = 0 if self.molecules else -1
        self.scroll_y = 0
        self.update()

    def render_molecule_to_texture(self, system):
        # --- Prepare Texture and FBO ---
        tex_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex_id)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, self.thumb_size[0], self.thumb_size[1], 0, GL_RGB, GL_UNSIGNED_BYTE, None)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        
        glBindFramebuffer(GL_FRAMEBUFFER, self.fbo)
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_id, 0)
        
        glBindRenderbuffer(GL_RENDERBUFFER, self.render_buffer)
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, self.thumb_size[0], self.thumb_size[1])
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, self.render_buffer)

        # --- Render Scene to FBO ---
        glViewport(0, 0, self.thumb_size[0], self.thumb_size[1])
        glClearColor(0.2, 0.2, 0.25, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_DEPTH_TEST)

        # Set molecule data for rendering (similar to detailed view)
        self.set_detail_molecule_data(system)
        self.draw_detail_scene_elements(is_thumbnail=True)

        # --- Unbind FBO and restore viewport ---
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        glViewport(0, 0, self.width(), self.height())
        return tex_id

    def paintGL(self):
        # Override BaseGLWidget.paintGL_base to control uniforms manually per mode
        if self.view_mode == MODE_THUMBNAIL:
            glClearColor(0.1, 0.1, 0.1, 1.0)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glDisable(GL_DEPTH_TEST)
            self.draw_thumbnail_view()
        else: # MODE_DETAIL
            glEnable(GL_DEPTH_TEST)
            # Call the original base method for 3D trackball view
            super().paintGL_base()
            
    def draw_scene(self):
        # This is called by paintGL_base only in MODE_DETAIL
        if self.view_mode == MODE_DETAIL:
            self.draw_detail_scene_elements()

    def draw_thumbnail_view(self):
        width, height = self.width(), self.height()
        proj = QMatrix4x4()
        proj.ortho(0, width, height, 0, -1, 1)

        total_grid_width = self.grid_cols * (self.thumb_size[0] + self.spacing[0]) - self.spacing[0]
        start_x = (width - total_grid_width) / 2

        # Draw molecules
        for i, (fname, _, tex_id) in enumerate(self.molecules):
            row, col = divmod(i, self.grid_cols)
            
            pos_x = start_x + col * (self.thumb_size[0] + self.spacing[0])
            pos_y = self.spacing[1] + row * (self.thumb_size[1] + self.spacing[1]) + self.scroll_y
            
            # Culling: Don't draw if off-screen
            if pos_y > height or pos_y + self.thumb_size[1] < 0:
                continue

            # Draw thumbnail sprite
            glUseProgram(self.sprite_shader)
            glUniformMatrix4fv(glGetUniformLocation(self.sprite_shader, "projection"), 1, GL_FALSE, proj.data())
            glUniform2f(glGetUniformLocation(self.sprite_shader, "u_pos"), pos_x, pos_y)
            glUniform2f(glGetUniformLocation(self.sprite_shader, "u_scale"), self.thumb_size[0], self.thumb_size[1])
            glActiveTexture(GL_TEXTURE0)
            glBindTexture(GL_TEXTURE_2D, tex_id)
            self.quad_mesh.draw(GL_TRIANGLE_FAN)

            # Draw text label
            self.render_text(fname, pos_x + self.thumb_size[0]/2, pos_y + self.thumb_size[1] + 20, 0.5, proj)

        # Draw cursor
        if self.selected_index >= 0:
            row, col = divmod(self.selected_index, self.grid_cols)
            pos_x = start_x + col * (self.thumb_size[0] + self.spacing[0])
            pos_y = self.spacing[1] + row * (self.thumb_size[1] + self.spacing[1]) + self.scroll_y
            
            glUseProgram(self.flat_shader)
            glUniformMatrix4fv(glGetUniformLocation(self.flat_shader, "projection"), 1, GL_FALSE, proj.data())
            glUniform2f(glGetUniformLocation(self.flat_shader, "u_pos"), pos_x, pos_y)
            glUniform2f(glGetUniformLocation(self.flat_shader, "u_scale"), self.thumb_size[0], self.thumb_size[1])
            glUniform4f(glGetUniformLocation(self.flat_shader, "u_color"), 0.2, 1.0, 0.2, 1.0)
            self.cursor_mesh.draw(GL_LINE_LOOP)
            
        self.update_content_height()

    def render_text(self, text, x, y, scale, proj_matrix):
        glUseProgram(self.text_shader)
        glUniformMatrix4fv(glGetUniformLocation(self.text_shader, "projection"), 1, GL_FALSE, proj_matrix.data())
        glActiveTexture(GL_TEXTURE0)
        glBindTexture(GL_TEXTURE_2D, self.font_atlas_tex)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        # Center text
        total_w = sum(self.font_atlas_data['chars'][c]['w'] for c in text if c in self.font_atlas_data['chars']) * scale
        cursor_x = x - total_w / 2

        for char in text:
            if char not in self.font_atlas_data['chars']: continue
            cd = self.font_atlas_data['chars'][char]
            
            w, h = cd['w'] * scale, cd['h'] * scale
            verts = np.array([cursor_x, y, cursor_x+w, y, cursor_x+w, y+h, cursor_x, y+h], dtype=np.float32)
            uvs = np.array([cd['norm_x'], cd['norm_y'], cd['norm_x']+cd['norm_w'], cd['norm_y'],
                            cd['norm_x']+cd['norm_w'], cd['norm_y']+cd['norm_h'], cd['norm_x'], cd['norm_y']+cd['norm_h']], dtype=np.float32)

            self.text_mesh.vertices = verts.reshape((4,2))
            self.text_mesh.normals = uvs.reshape((4,2)) # Using normals VBO for UVs
            
            glBindBuffer(GL_ARRAY_BUFFER, self.text_mesh.vbo_vertices)
            glBufferData(GL_ARRAY_BUFFER, self.text_mesh.vertices.nbytes, self.text_mesh.vertices, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, self.text_mesh.vbo_normals)
            glBufferData(GL_ARRAY_BUFFER, self.text_mesh.normals.nbytes, self.text_mesh.normals, GL_DYNAMIC_DRAW)
            
            self.text_mesh.draw(GL_TRIANGLE_FAN)
            cursor_x += w
            
        glDisable(GL_BLEND)

    def draw_detail_scene_elements(self, is_thumbnail=False):
        # Camera for thumbnail rendering
        if is_thumbnail:
            system = self.detail_system
            max_dim = np.max(np.ptp(system.apos, axis=0)) if len(system.apos) > 0 else 10.0
            cam_dist = max_dim * 1.5
            proj = QMatrix4x4(); proj.perspective(45.0, 1.0, 0.1, cam_dist * 2)
            view = QMatrix4x4(); view.lookAt(QVector3D(0,0,cam_dist), QVector3D(0,0,0), QVector3D(0,1,0))
            model = QMatrix4x4()

            for prog_id in [self.atom_shader, self.bond_shader]:
                glUseProgram(prog_id)
                glUniformMatrix4fv(glGetUniformLocation(prog_id, "projection"), 1, GL_FALSE, proj.data())
                glUniformMatrix4fv(glGetUniformLocation(prog_id, "view"), 1, GL_FALSE, view.data())
                glUniformMatrix4fv(glGetUniformLocation(prog_id, "model"), 1, GL_FALSE, model.data())
                glUniform3f(glGetUniformLocation(prog_id, "lightPos"), cam_dist, cam_dist, cam_dist)
                glUniform3f(glGetUniformLocation(prog_id, "viewPos"), 0, 0, cam_dist)

        # Draw Atoms
        glUseProgram(self.atom_shader)
        self.atom_instances.draw(mode=GL_TRIANGLES)

        # Draw Bonds
        glUseProgram(self.bond_shader)
        self.bond_instances.associated_mesh.vertex_count = 34
        self.bond_instances.draw(mode=GL_TRIANGLE_STRIP)

    def set_detail_molecule_data(self, system):
        self.detail_system = system
        if system.bonds is None: system.findBonds()
        
        n_atoms = len(system.apos)
        atom_pos = system.apos.astype(np.float32)
        atom_radii = np.array([elements.ELEMENT_DICT[e][6] for e in system.enames], dtype=np.float32)
        atom_colors = np.array([elements.hex_to_float_rgb(elements.ELEMENT_DICT[e][7]) + (1.0,) for e in system.enames], dtype=np.float32)
        atom_mats = np.array([np.eye(4, dtype=np.float32) for _ in range(n_atoms)])
        self.atom_instances.update({"positions": atom_pos, "radii": atom_radii, "colors": atom_colors, "instanceMatrix": atom_mats})

        if system.bonds is not None and len(system.bonds) > 0:
            n_bonds = len(system.bonds)
            bond_pos = np.zeros((n_bonds, 3), dtype=np.float32)
            bond_colors = np.zeros((n_bonds, 4), dtype=np.float32)
            bond_mats = np.zeros((n_bonds, 4, 4), dtype=np.float32)
            for i, (ia, ib) in enumerate(system.bonds):
                p1, p2 = system.apos[ia], system.apos[ib]
                vec = p2 - p1
                length = np.linalg.norm(vec)
                bond_pos[i] = (p1 + p2) / 2
                bond_colors[i] = (atom_colors[ia] + atom_colors[ib]) / 2
                bond_mats[i] = makeRotMat(vec, np.array([0.,0.,1.])) @ np.diag([1.,1.,length,1.])
            self.bond_instances.update({"positions":bond_pos, "radii":np.zeros(n_bonds), "colors":bond_colors, "instanceMatrix": bond_mats})
        else:
            self.bond_instances.num_instances = 0


    def keyPressEvent(self, event):
        if self.view_mode == MODE_THUMBNAIL:
            if not self.molecules: return
            key = event.key()
            new_index = self.selected_index
            
            if key == Qt.Key_Right: new_index += 1
            elif key == Qt.Key_Left: new_index -= 1
            elif key == Qt.Key_Down: new_index += self.grid_cols
            elif key == Qt.Key_Up: new_index -= self.grid_cols
            elif key in [Qt.Key_Return, Qt.Key_Enter]:
                if self.selected_index >= 0:
                    self.view_mode = MODE_DETAIL
                    self.set_detail_molecule_data(self.molecules[self.selected_index][1])
            
            if 0 <= new_index < len(self.molecules):
                self.selected_index = new_index
                self.ensure_selection_visible()
            self.update()
        
        elif self.view_mode == MODE_DETAIL:
            if event.key() == Qt.Key_Escape:
                self.view_mode = MODE_THUMBNAIL
                self.orientation = R.identity() # Reset camera
                self.update()
            else:
                super().keyPressEvent(event)

    def wheelEvent(self, event):
        if self.view_mode == MODE_THUMBNAIL:
            self.scroll_y += event.angleDelta().y() * 0.5
            self.scroll_y = min(0, self.scroll_y) # Clamp top
            self.scroll_y = max(-(self.content_height - self.height() + self.spacing[1]), self.scroll_y) # Clamp bottom
            self.update()
        else:
            super().wheelEvent(event)

    def update_content_height(self):
        if not self.molecules:
            self.content_height = 0
            return
        num_rows = (len(self.molecules) - 1) // self.grid_cols + 1
        self.content_height = num_rows * (self.thumb_size[1] + self.spacing[1])

    def ensure_selection_visible(self):
        row, _ = divmod(self.selected_index, self.grid_cols)
        item_y_top = self.spacing[1] + row * (self.thumb_size[1] + self.spacing[1]) + self.scroll_y
        item_y_bot = item_y_top + self.thumb_size[1] + self.spacing[1]
        
        if item_y_top < 0:
            self.scroll_y -= item_y_top
        elif item_y_bot > self.height():
            self.scroll_y -= (item_y_bot - self.height())


# --- Main Window Class ---
class MolecularBrowser(AppWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("OpenGL Molecule Browser")
        self.setGeometry(100, 100, 1200, 800)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        layout.setContentsMargins(0, 5, 0, 5)

        self.gl_widget = MoleculeBrowserGLWidget()
        
        btn = QPushButton("Open Directory")
        btn.clicked.connect(self.open_dir)
        
        layout.addWidget(btn)
        layout.addWidget(self.gl_widget)
        
    def open_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir_path:
            self.gl_widget.load_directory(dir_path)

if __name__ == '__main__':
    from scipy.spatial.transform import Rotation as R # Add this import here
    MolecularBrowser.launch()
```

### How to Run

1.  **Prepare**:
    *   Install dependencies: `pip install PyQt5 numpy scipy pillow PyOpenGL`.
    *   Create the `shaders` directory and place all the shader files (`sprite.glslv`, `sprite.glslf`, `text.glslv`, `text.glslf`, `flat_color.glslf`, `instances.glslv`, `sphere.glslf`, `cylinder.glslf`) inside it.
    *   Run `python make_font.py` to generate `font_atlas.png` and `font_atlas.json`.
    *   Make sure `MolecularBrowser.py`, `GLGUI.py`, `AtomicSystem.py` and your `pyBall` package are in a directory structure that Python can resolve. For simplicity, you can place them all in the same directory.

2.  **Launch**:
    ```bash
    python MolecularBrowser.py
    ```

This revised application now perfectly matches your request: it operates in a single OpenGL window, switches between a manually rendered thumbnail view and a detailed 3D view, and avoids any complex context management.