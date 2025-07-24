# MolecularBrowser.py
import sys
import os
import numpy as np
import json
from PIL import Image
from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QPushButton, QFileDialog, QApplication)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, Mesh, octahedron_sphere_mesh
from pyBall.AtomicSystem import AtomicSystem
from pyBall import elements
from pyBall.atomicUtils import makeRotMat

# --- Helper to create cylinder mesh (needed for ball-and-stick) ---

def make_bond_transform_matrix(p1, p2):
    """Create a 4x4 transformation matrix for a bond between p1 and p2.
    This positions and orients a cylinder to connect the two points.
    """
    # Get vector from p1 to p2
    direction = p2 - p1
    
    # Bond length
    length = np.linalg.norm(direction)
    if length < 1e-6:  # Avoid division by zero
        return np.identity(4, dtype=np.float32)
    
    # Normalize direction
    direction = direction / length
    
    # Default cylinder is along Z-axis, so we need to rotate from Z to our direction
    # Find rotation axis and angle
    z_axis = np.array([0, 0, 1], dtype=np.float32)
    
    # Check if direction is parallel to z_axis
    if np.allclose(direction, z_axis) or np.allclose(direction, -z_axis):
        # Special case: direction is parallel to z-axis
        if np.allclose(direction, -z_axis):
            # If pointing down, rotate 180 degrees around X-axis
            rot_matrix = np.array([
                [1, 0, 0],
                [0, -1, 0],
                [0, 0, -1]
            ], dtype=np.float32)
        else:
            # Already aligned with z-axis
            rot_matrix = np.identity(3, dtype=np.float32)
    else:
        # Normal case: compute rotation from z-axis to direction
        rotation_axis = np.cross(z_axis, direction)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        
        # Compute angle
        cos_angle = np.dot(z_axis, direction)
        sin_angle = np.linalg.norm(np.cross(z_axis, direction))
        
        # Build rotation matrix using Rodrigues' formula
        K = np.array([
            [0, -rotation_axis[2], rotation_axis[1]],
            [rotation_axis[2], 0, -rotation_axis[0]],
            [-rotation_axis[1], rotation_axis[0], 0]
        ], dtype=np.float32)
        
        rot_matrix = np.identity(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
    
    # Scale along z-axis for bond length
    scale_matrix = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, length]
    ], dtype=np.float32)
    
    # Combine rotation and scaling
    transform_3x3 = np.dot(rot_matrix, scale_matrix)
    
    # Convert to 4x4 matrix with translation to midpoint
    transform_4x4 = np.identity(4, dtype=np.float32)
    transform_4x4[:3, :3] = transform_3x3
    
    # Position at midpoint between atoms
    midpoint = (p1 + p2) / 2
    transform_4x4[:3, 3] = midpoint
    
    return transform_4x4
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
    def __init__(self, parent=None, dir_path=None):
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
        self.shaders = {}
        self.sphere_mesh = None
        self.cylinder_mesh = None
        self.atom_instances = None
        self.bond_instances = None
        self.fbo = None
        self.proj_matrix = None
        self.view_matrix = None
        
        # Create a simple Carbon atom system for debugging
        apos = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        atypes = np.array([elements.ELEMENT_DICT['C'][0]], dtype=np.int32)
        enames = ['C']
        qs = np.array([0.0], dtype=np.float32)
        self._debug_system = AtomicSystem(apos=apos, atypes=atypes, enames=enames, qs=qs, bPreinit=True)
        self.selected_index = 0
        self.scroll_y = 0

        # Store the directory path for later loading after OpenGL is initialized
        self.initial_dir_path = dir_path

    @property
    def all_shader_programs(self):
        progs = [self.sprite_shader, self.text_shader, self.flat_shader, self.atom_shader, self.bond_shader]
        return [p for p in progs if p is not None]

    def initializeGL(self):
        print("DEBUG: initializeGL called - OpenGL initialization starting")

        # --- Load Shaders ---
        shader_dir = os.path.join(os.path.dirname(__file__), "shaders") + os.sep
        def load_shader(vert_name, frag_name):
            try:
                with open(shader_dir + vert_name, 'r') as f: v = f.read()
                with open(shader_dir + frag_name, 'r') as f: f = f.read()
                return self.compile_shader_program(v, f)
            except Exception as e:
                print(f"Error loading shader {vert_name}/{frag_name}: {e}")
                return None

        self.sprite_shader = load_shader("sprite.glslv", "sprite.glslf")
        self.text_shader   = load_shader("text.glslv",   "text.glslf")
        self.flat_shader   = load_shader("sprite.glslv", "flat_color.glslf")
        self.atom_shader   = load_shader("instances.glslv", "sphere.glslf")
        self.bond_shader   = load_shader("instances.glslv", "cylinder.glslf")

        super().initializeGL_base(None, None)

        # --- Initialize Meshes ---

        # Sphere for atoms
        verts, norms = octahedron_sphere_mesh(nsub=2)
        sphere_mesh = Mesh(verts, norms)
        sphere_mesh.setup_buffers()

        # Quad for sprites
        quad_verts = np.array([-1.0,-1.0,0.0, 1.0,-1.0,0.0, 1.0,1.0,0.0, -1.0,1.0,0.0], dtype=np.float32).reshape(4,3)
        self.quad_mesh = Mesh(quad_verts)
        self.quad_mesh.setup_buffers()

        # Outline for cursor
        cursor_verts = np.array([0.0,0.0,0.0, 0.0,1.0,0.0, 1.0,1.0,0.0, 1.0,0.0,0.0], dtype=np.float32).reshape(4,3)
        self.cursor_mesh = Mesh(cursor_verts)
        self.cursor_mesh.setup_buffers()

        # Cylinder for bonds
        verts, norms = create_cylinder_mesh()
        cylinder_mesh = Mesh(verts, norms)
        cylinder_mesh.setup_buffers()

        # --- Initialize Instanced Data ---
        # For atoms (spheres)
        self.atom_instances = InstancedData(base_attrib_location=2)
        self.atom_instances.associate_mesh(sphere_mesh)
        atom_attribs = [
            ("positions", 0, 3),
            ("radii",     1, 1),
            ("colors",    2, 4),
        ]
        self.atom_instances.setup_instance_vbos(atom_attribs)

        # For bonds (cylinders)
        self.bond_instances = InstancedData(base_attrib_location=2)
        self.bond_instances.associate_mesh(cylinder_mesh)
        bond_attribs = [
            ("instanceMatrix_row0", 0, 4),
            ("instanceMatrix_row1", 1, 4),
            ("instanceMatrix_row2", 2, 4),
            ("instanceMatrix_row3", 3, 4),
            ("colors",              4, 4),
        ]
        self.bond_instances.setup_instance_vbos(bond_attribs)

        # --- FBO for Offscreen Rendering ---
        self.fbo = glGenFramebuffers(1)
        self.render_buffer = glGenRenderbuffers(1)

        # DEBUG: render hardcoded system after GL initialization
        tex_id = self.render_molecule_to_texture(self._debug_system)
        self.molecules.append(("Carbon.xyz", self._debug_system, tex_id))

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
        
        font_atlas_dir = os.path.join(os.path.dirname(__file__), "shaders") + os.sep # Ensure trailing slash
        self.font_atlas_tex = self.load_texture(os.path.join(font_atlas_dir, "font_atlas.png"))
        with open(os.path.join(font_atlas_dir, "font_atlas.json")) as f:
            self.font_atlas_data = json.load(f)

        # --- Detail View Rendering Setup ---
        self.atom_instances = InstancedData(base_attrib_location=2)
        self.atom_instances.associate_mesh(self.default_sphere_mesh)
        inst_attribs = [
        ("positions",0,3),
        ("radii",1,1),
        ("colors",2,4),
        ("instanceMatrix_row0",3,4),
        ("instanceMatrix_row1",4,4),
        ("instanceMatrix_row2",5,4),
        ("instanceMatrix_row3",6,4)
    ]
        self.atom_instances.setup_instance_vbos(inst_attribs)

        cyl_v, cyl_n = create_cylinder_mesh(radius=0.1, length=1.0)
        self.cylinder_mesh = Mesh(vertices=cyl_v, normals=cyl_n)
        self.cylinder_mesh.setup_buffers()

        self.bond_instances = InstancedData(base_attrib_location=2)
        self.bond_instances.associate_mesh(self.cylinder_mesh)
        self.bond_instances.setup_instance_vbos(inst_attribs)

        print("OpenGL Initialized")
        
        # Now that OpenGL is initialized, load the directory if one was specified
        if hasattr(self, 'initial_dir_path') and self.initial_dir_path is not None:
            print(f"DEBUG: Loading directory {self.initial_dir_path} now that OpenGL is initialized")
            try:
                self.load_directory(self.initial_dir_path)
                print(f"DEBUG: Directory loaded successfully")
            except Exception as e:
                print(f"ERROR loading directory: {type(e)}: {e}")

    def load_texture(self, filepath):
        try:
            img = Image.open(filepath).convert("RGBA")
        except FileNotFoundError:
            print(f"Error: Texture file not found at {filepath}")
            return 0

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
            
            system = AtomicSystem(fname=fpath)
            tex_id = self.render_molecule_to_texture(system)
            self.molecules.append((fname, system, tex_id))

        self.selected_index = 0 if self.molecules else -1
        self.scroll_y = 0
        self.update()

    def render_molecule_to_texture(self, system):
        print(f"DEBUG: render_molecule_to_texture called for {system}")
        if self.fbo is None:
            print("ERROR: FBO is None. OpenGL not initialized yet.")
            # Create a temporary texture ID as fallback
            return glGenTextures(1)
            
        # --- Prepare Texture and FBO ---
        tex_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex_id)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, self.thumb_size[0], self.thumb_size[1], 0, GL_RGB, GL_UNSIGNED_BYTE, None)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        
        print(f"DEBUG: Using FBO: {self.fbo}")
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

        # --- DEBUG: Inspect the rendered texture ---
        debug_texture_to_numpy(tex_id, self.thumb_size[0], self.thumb_size[1])

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

        print(f"DEBUG: cursor_mesh type: {type(self.cursor_mesh)}")
        print(f"DEBUG: cursor_mesh is None: {self.cursor_mesh is None}")
        print(f"DEBUG: quad_mesh type: {type(self.quad_mesh)}")
        print(f"DEBUG: quad_mesh is None: {self.quad_mesh is None}")
        

        
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
        # --- Camera Setup ---
        if is_thumbnail:
            # Special camera for thumbnail rendering to fit the molecule
            system = self.detail_system
            if system and system.apos is not None and len(system.apos) > 0:
                max_dim = np.max(np.ptp(system.apos, axis=0))
            else:
                max_dim = 10.0 # Default zoom if no atoms
            cam_dist = max_dim * 2.0 # Zoom out to see the whole molecule
            proj_matrix = QMatrix4x4()
            proj_matrix.perspective(45.0, 1.0, 0.1, cam_dist * 2)
            view_matrix = QMatrix4x4()
            view_matrix.lookAt(QVector3D(0,0,cam_dist), QVector3D(0,0,0), QVector3D(0,1,0))
            light_pos = QVector3D(cam_dist, cam_dist, cam_dist)
            view_pos  = QVector3D(0, 0, cam_dist)
        else:
            # Use the main trackball camera for detailed view
            proj_matrix = self.cam.pre_proj * self.cam.proj
            view_matrix = self.cam.pre_view * self.cam.view
            # These would need to be derived from the camera's state
            light_pos = QVector3D(10.0, 10.0, 10.0) 
            view_pos = QVector3D(0.0, 0.0, 10.0) # Simplified

        # --- Draw Atoms ---
        print(f"DEBUG: atom_instances: {self.atom_instances}")
        if self.atom_instances: print(f"DEBUG: atom_instances.num_instances = {self.atom_instances.num_instances}")
        if self.atom_instances and self.atom_instances.num_instances > 0:
            glUseProgram(self.atom_shader)
            loc_view = glGetUniformLocation(self.atom_shader, "view")
            loc_proj = glGetUniformLocation(self.atom_shader, "proj")
            loc_light = glGetUniformLocation(self.atom_shader, "lightPos")
            loc_viewer = glGetUniformLocation(self.atom_shader, "viewPos")
            
            print(f"DEBUG: atom_shader uniforms - view:{loc_view}, proj:{loc_proj}, light:{loc_light}, viewer:{loc_viewer}")

            glUniformMatrix4fv(loc_view, 1, GL_FALSE, view_matrix.data())
            glUniformMatrix4fv(loc_proj, 1, GL_FALSE, proj_matrix.data())
            glUniform3f(loc_light, light_pos.x(), light_pos.y(), light_pos.z())
            glUniform3f(loc_viewer, view_pos.x(), view_pos.y(), view_pos.z())
            self.atom_instances.draw()

        # --- Draw Bonds ---
        print(f"DEBUG: bond_instances: {self.bond_instances}")
        if self.bond_instances: print(f"DEBUG: bond_instances.num_instances = {self.bond_instances.num_instances}")
        if self.bond_instances and self.bond_instances.num_instances > 0:
            glUseProgram(self.bond_shader)
            loc_view = glGetUniformLocation(self.bond_shader, "view")
            loc_proj = glGetUniformLocation(self.bond_shader, "proj")
            loc_light = glGetUniformLocation(self.bond_shader, "lightPos")
            loc_viewer = glGetUniformLocation(self.bond_shader, "viewPos")

            print(f"DEBUG: bond_shader uniforms - view:{loc_view}, proj:{loc_proj}, light:{loc_light}, viewer:{loc_viewer}")

            glUniformMatrix4fv(loc_view, 1, GL_FALSE, view_matrix.data())
            glUniformMatrix4fv(loc_proj, 1, GL_FALSE, proj_matrix.data())
            glUniform3f(loc_light, light_pos.x(), light_pos.y(), light_pos.z())
            glUniform3f(loc_viewer, view_pos.x(), view_pos.y(), view_pos.z())
            self.bond_instances.draw()

    def update_content_height(self):
        # Placeholder for now. Implement actual logic later if needed.
        pass

    def set_detail_molecule_data(self, system):
        print(f"DEBUG: set_detail_molecule_data called with system: {system}")
        print(f"DEBUG: system.apos type: {type(system.apos)}, shape: {system.apos.shape if hasattr(system.apos, 'shape') else 'no shape'}, dtype: {system.apos.dtype if hasattr(system.apos, 'dtype') else 'no dtype'}")
        print(f"DEBUG: system.enames: {system.enames}")
        
        self.detail_system = system
        if system.bonds is None: 
            print(f"DEBUG: Finding bonds for system")
            system.findBonds()
            print(f"DEBUG: Bonds found: {system.bonds}")

        # --- Atoms ---
        n_atoms = len(system.apos)
        print(f"DEBUG: Number of atoms: {n_atoms}")
        try:
            atom_pos = system.apos.astype(np.float32)
            print(f"DEBUG: atom_pos shape: {atom_pos.shape}, dtype: {atom_pos.dtype}")
            
            # Check element names
            print(f"DEBUG: Element names: {system.enames}")
            try:
                atom_radii = np.array([elements.ELEMENT_DICT[e][6] for e in system.enames], dtype=np.float32)
                print(f"DEBUG: atom_radii: {atom_radii}")
            except Exception as e:
                print(f"ERROR in atom_radii: {e}")
                # Fallback with default radius
                atom_radii = np.full(n_atoms, 0.5, dtype=np.float32)
                
            try:
                atom_colors = np.array([elements.hex_to_float_rgb(elements.ELEMENT_DICT[e][8]) + (1.0,) for e in system.enames], dtype=np.float32)
                print(f"DEBUG: atom_colors shape: {atom_colors.shape}, dtype: {atom_colors.dtype}")
            except Exception as e:
                print(f"ERROR in atom_colors: {e}")
                # Fallback with default color (gray)
                atom_colors = np.tile([0.5, 0.5, 0.5, 1.0], (n_atoms, 1)).astype(np.float32)
            
            print(f"DEBUG: About to update atom_instances")
            print(f"DEBUG: atom_instances: {self.atom_instances}")
            self.atom_instances.update({"positions": atom_pos, "radii": atom_radii, "colors": atom_colors})
            print(f"DEBUG: Successfully updated atom_instances")
        except Exception as e:
            print(f"ERROR updating atom data: {e}")

        # --- Bonds ---
        try:
            bonds = system.bonds
            print(f"DEBUG: Bonds: {bonds}")
            if bonds is not None and len(bonds) > 0:
                try:
                    print(f"DEBUG: Processing {len(bonds)} bonds")
                    bond_indices1 = np.array([b[0] for b in bonds])
                    bond_indices2 = np.array([b[1] for b in bonds])
                    print(f"DEBUG: Bond indices1: {bond_indices1}, indices2: {bond_indices2}")
                    
                    bond_pos1 = atom_pos[bond_indices1]
                    bond_pos2 = atom_pos[bond_indices2]
                    print(f"DEBUG: bond_pos1 shape: {bond_pos1.shape}, bond_pos2 shape: {bond_pos2.shape}")
                    
                    try:
                        # Use our custom function that creates 4x4 matrices directly
                        bond_mats = np.array([make_bond_transform_matrix(p1, p2) for p1, p2 in zip(bond_pos1, bond_pos2)], dtype=np.float32)
                        print(f"DEBUG: bond_mats shape after transformation: {bond_mats.shape}")
                        print(f"DEBUG: bond_mats shape: {bond_mats.shape}")
                    except Exception as e:
                        print(f"ERROR creating bond matrices: {e}")
                        return
                        
                    bond_colors = np.array([[0.5,0.5,0.5,1.0]]*len(bonds), dtype=np.float32)
                    print(f"DEBUG: bond_colors shape: {bond_colors.shape}")

                    # Prepare buffer data for the mat4 - no need to reshape since we already have 4x4 matrices
                    mats_r = bond_mats.reshape(len(bonds) * 4, 4)  # This should work now with 4x4 matrices
                    print(f"DEBUG: mats_r shape: {mats_r.shape}")
                    
                    # Calculate bond positions as the midpoints between atoms (required by InstancedData.update)
                    bond_positions = np.array([(p1 + p2) / 2 for p1, p2 in zip(bond_pos1, bond_pos2)], dtype=np.float32)
                    print(f"DEBUG: bond_positions shape: {bond_positions.shape}")
                    
                    # Include positions in bond_buffs (required by update method)
                    bond_buffs = {
                        "positions":           bond_positions,      # Required by InstancedData.update
                        "instanceMatrix_row0": mats_r[0::4],
                        "instanceMatrix_row1": mats_r[1::4],
                        "instanceMatrix_row2": mats_r[2::4],
                        "instanceMatrix_row3": mats_r[3::4],
                        "colors":              bond_colors
                    }
                    print(f"DEBUG: About to update bond_instances")
                    print(f"DEBUG: bond_instances: {self.bond_instances}")
                    self.bond_instances.update(bond_buffs)
                    print(f"DEBUG: Successfully updated bond_instances")
                except Exception as e:
                    print(f"ERROR processing bonds: {e}")
            else:
                print("DEBUG: No bonds to process")
        except Exception as e:
            print(f"ERROR in bond section: {e}")


def debug_texture_to_numpy(tex_id, width, height):
    """Downloads a texture from GPU to a NumPy array and prints min/max values."""
    glBindTexture(GL_TEXTURE_2D, tex_id)
    # To use glGetTexImage, the texture must be bound to the active texture unit
    glActiveTexture(GL_TEXTURE0)
    # Allocate numpy array to hold the texture data
    img_data = np.zeros((height, width, 4), dtype=np.uint8)
    # Get the texture image data
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data)
    glBindTexture(GL_TEXTURE_2D, 0)

    if img_data.size > 0:
        min_vals = img_data.min(axis=(0, 1))
        max_vals = img_data.max(axis=(0, 1))
        print(f"--- Texture Debug (ID: {tex_id}) ---")
        print(f"Shape: {img_data.shape}")
        print(f"Min (R,G,B,A): {min_vals}")
        print(f"Max (R,G,B,A): {max_vals}")
        print("---------------------------")
    else:
        print("--- Texture Debug: Error, no data downloaded. ---")
    return img_data

if __name__ == '__main__':


    # Run Like this:
    #   python -m pyBall.GUI.MolecularBrowser

    # relative path to the directory with the molecules
    # /home/prokop/git/FireCore/cpp/common_resources/xyz/
    dir_path = "../../cpp/common_resources/xyz/"
    dir_path = "../../cpp/common_resources/xyz_mini/"
    this_path = os.path.dirname(os.path.abspath(__file__))
    # resolve the full absolute path
    dir_path = os.path.abspath(os.path.join(this_path, dir_path))
    
    app = QApplication(sys.argv)
    window = MoleculeBrowserGLWidget(dir_path=dir_path)
    window.show()
    sys.exit(app.exec_())