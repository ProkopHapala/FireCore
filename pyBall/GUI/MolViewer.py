from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QVBoxLayout, QHBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
import argparse
import os

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes, GLobject, rotation_to_gl_matrix
from pyBall.OGL.BaseGUI import BaseGUI
from pyBall.atomicUtils import findAllBonds, findBondsNP

import OpenGL.GL as GL
import json
import ctypes
from PIL import Image

np.set_printoptions(linewidth=np.inf)

class FrameData:
    """Unified container for all rendering data for a single frame"""
    def __init__(self):
        self.atoms = None      # [positions, radii, colors, names]
        self.electrons = None  # [positions, radii, colors, names]
        self.bonds = None      # list of (i,j) pairs
        self.bond_vecs = None  # list of (length, direction) tuples
        
    def compute_bonds(self, Rcut=3.0, RvdwCut=0.7):
        """Compute bonds using findAllBonds"""
        if self.atoms is None:
            return
        #print("self.atoms[0] ", self.atoms[0])
        #print("self.atoms[1] ", self.atoms[1])
        RvdWs = np.array( [ elements.ELEMENT_DICT[e][7] for e in self.atoms[3] ])
        #print("RvdWs ", RvdWs)
        #bonds, bond_vecs = findAllBonds(self.atoms[0], self.atoms[1] )
        bonds, bond_vecs = findBondsNP(self.atoms[0], RvdWs=RvdWs,  )
        #print("bonds ", bonds)
        self.bonds     = bonds
        self.bond_vecs = bond_vecs

class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.trj = []
        self.current_frame_index = 0
        self.opacity      = 0.5
        self.frames_data = []  # List of FrameData objects
        self.atom_instances = InstancedData(base_attrib_location=1)
        self.elec_instances = InstancedData(base_attrib_location=1)
        # self.blend_mode = alpha_blend_modes["standard"] # Replaced by render_modes
        self.shader_program_sphere_raytrace = None
        self.shader_program_sphere_max_vol = None
        self.instance_data_dirty = True
        self.electron_names = set(["E", "e", "e2", "e+", "e-"])
        self.render_modes = {}
        self.current_render_mode_key = None

        self.atom_labels = None
        #self.atom_labels = GLobject()

        self.num_label_verts = 0
        self.label_data_dirty = True

        # Bond rendering setup
        self.bonds_vao = None
        self.bonds_vbo = None
        self.bonds_ebo = None
        self.bond_shader = None
        self.frames_bonds = []  # Store bonds for each frame
        self.bond_data_dirty = True

        # FBO related attributes
        self.fbo_id = None
        self.fbo_texture_id = None
        self.fbo_depth_buffer_id = None
        self.fbo_width = 0
        self.fbo_height = 0

    def initializeGL(self):
        #print("MolViewerWidget.initializeGL()")
        # Call base class initialization for shaders, basic GL setup
        shader_folder=self.shader_folder
        print("shader_folder", shader_folder)
        vert_instances   = open(shader_folder + "/instances.glslv").read()
        frag_ray_sphere  = open(shader_folder + "/sphere.glslf").read()
        frag_max_sphere  = open(shader_folder + "/sphere_max.glslf").read() 
        #FRAGMENT_SHADER_MAX_VOL_SOURCE = open(shader_folder + "/sphere_simple.glslf").read() 
        self.shader_program_sphere_raytrace = self.compile_shader_program( vert_instances, frag_ray_sphere)
        self.shader_program_sphere_max_vol  = self.compile_shader_program( vert_instances, frag_max_sphere)

        # Set the default shader program for BaseGLWidget's uniform setup
        self.shader_program = self.shader_program_sphere_raytrace
        if not self.shader_program: self.shader_program = self.shader_program_sphere_max_vol
        
        self.all_shader_programs.append(self.shader_program_sphere_raytrace)
        self.all_shader_programs.append(self.shader_program_sphere_max_vol)
            
        # Initialize bond rendering
        vert_bond = open(shader_folder + "/cylinder.glslv").read()
        frag_bond = open(shader_folder + "/cylinder.glslf").read()
        self.bond_shader = self.compile_shader_program(vert_bond, frag_bond)

        self.all_shader_programs.append(self.bond_shader)

        self.bonds = GLobject(components=[3], mode=GL.GL_LINES)
        
        # Define render modes after shaders are compiled
        self.render_modes = {
            "raytrace": (self.shader_program_sphere_raytrace, alpha_blend_modes["standard"], True),
            "max_vol":  (self.shader_program_sphere_max_vol,  alpha_blend_modes["standard"], False)
        }
        self.current_render_mode_key = "raytrace"  # Set default render mode

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

        self.font_atlas_tex = self.load_texture(shader_folder + "/font_atlas.png")
        print(f"DEBUG: font_atlas_tex loaded: {self.font_atlas_tex}")
        with open(shader_folder + "/font_atlas.json") as f:
            self.font_atlas_data = json.load(f)
        print(f"DEBUG: font_atlas_data loaded: {self.font_atlas_data}")

    def cleanupGL(self):
        super().cleanupGL()
        if self.text_shader: GL.glDeleteProgram(self.text_shader)
        if self.atom_labels: self.atom_labels.cleanup()
        if self.shader_program_sphere_raytrace: GL.glDeleteProgram(self.shader_program_sphere_raytrace)
        if self.shader_program_sphere_max_vol: GL.glDeleteProgram(self.shader_program_sphere_max_vol)
        self.shader_program_sphere_raytrace = None
        self.shader_program_sphere_max_vol = None
        if self.atom_instances:  self.atom_instances.cleanup()
        if self.elec_instances:  self.elec_instances.cleanup()
        self._deinit_fbo()
        # if self.bonds_vao:       GL.glDeleteVertexArrays(1, [self.bonds_vao])
        # if self.bonds_vbo:       GL.glDeleteBuffers(1, [self.bonds_vbo])
        # if self.bonds_ebo:       GL.glDeleteBuffers(1, [self.bonds_ebo])
        #if self.bond_shader:     GL.glDeleteProgram(self.bond_shader)

    def _init_fbo(self, width, height):
        # De-initialize existing FBO first to avoid resource leaks
        self._deinit_fbo()

        self.fbo_width = width
        self.fbo_height = height

        # Create FBO
        self.fbo_id = GL.glGenFramebuffers(1)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo_id)

        # Create texture to render to
        self.fbo_texture_id = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.fbo_texture_id)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, self.fbo_width, self.fbo_height, 0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, None)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D, self.fbo_texture_id, 0)

        # Create depth buffer
        self.fbo_depth_buffer_id = GL.glGenRenderbuffers(1)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, self.fbo_depth_buffer_id)
        GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT24, self.fbo_width, self.fbo_height)
        GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT, GL.GL_RENDERBUFFER, self.fbo_depth_buffer_id)

        # Check if FBO is complete
        if GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER) != GL.GL_FRAMEBUFFER_COMPLETE:
            print("Error: Framebuffer is not complete!")

        # Unbind FBO to return to default framebuffer
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)

    def _deinit_fbo(self):
        if self.fbo_id:
            GL.glDeleteFramebuffers(1, [self.fbo_id])
            self.fbo_id = None
        if self.fbo_texture_id:
            GL.glDeleteTextures(1, [self.fbo_texture_id])
            self.fbo_texture_id = None
        if self.fbo_depth_buffer_id:
            GL.glDeleteRenderbuffers(1, [self.fbo_depth_buffer_id])
            self.fbo_depth_buffer_id = None

    def hex2rgb(self, hex_color_str):
        return [int(hex_color_str[i:i+2], 16) / 255.0 for i in (0, 2, 4)]

    def precalculate_frames_data(self, dtype=np.float32):
        self.frames_data = []
        #index_color = elements.index_color
        for frame_idx in range(len(self.trj)):
            frame_data = FrameData()
            es = self.trj[frame_idx][0]
            na = len(es)
            ps = self.trj[frame_idx][1]
            apos,arad,acol,aname=[],[],[],[]
            epos,erad,ecol,ename=[],[],[],[]
            if na > 0:
                for i in range(na):
                    atom_symbol = es[i]
                    print(f"i: {i}, atom_symbol: {atom_symbol}")
                    current_radius = self.trj[frame_idx][3][i]
                    rec = elements.ELEMENT_DICT[atom_symbol]
                    #r_col, g_col, b_col = elements.ELEMENTS[ elements.ELEMENT_DICT[atom_symbol] ][1][:3]
                    r_col, g_col, b_col = self.hex2rgb(rec[elements.index_color][1:])
                    if atom_symbol == "e":
                        epos.append(ps[i])
                        erad.append(current_radius)
                        ecol.append([r_col, g_col, b_col, self.opacity])
                        ename.append(atom_symbol)
                    else: 
                        apos.append(ps[i])
                        arad.append(current_radius)
                        acol.append([r_col, g_col, b_col, self.opacity])
                        aname.append(atom_symbol)
            frame_data.atoms = [np.array(apos,dtype=dtype), np.array(arad,dtype=dtype), np.array(acol,dtype=dtype), aname]
            frame_data.electrons = [np.array(epos,dtype=dtype), np.array(erad,dtype=dtype), np.array(ecol,dtype=dtype), ename]
            frame_data.compute_bonds()
            self.frames_data.append(frame_data)

    def update_instance_data(self):
        atoms = self.frames_data[self.current_frame_index].atoms
        elecs = self.frames_data[self.current_frame_index].electrons
        self.atom_instances.update_list( atoms[:3] )
        self.elec_instances.update_list( elecs[:3] )
        self.instance_data_dirty = False

    def paintGL(self):
        # Calls paintGL_base, which in turn calls draw_scene
        super().paintGL_base() 

    def render_bonds_shader(self):
        print("render_bonds_shader()")

        frame_data = self.frames_data[self.current_frame_index]
        bonds = frame_data.bonds
        self.use_shader(self.bond_shader)
        self.set_default_uniforms()
                
        # Upload bond indices (point pairs)
        bond_indices = np.array(bonds, dtype=np.uint32).flatten()
        
        # Draw bonds
        GL.glBindVertexArray(self.bonds_vao)
        GL.glDrawElements(GL.GL_LINES, len(bond_indices), GL.GL_UNSIGNED_INT, None)
        GL.glBindVertexArray(0)

    def draw_scene(self):
        #print("!!!!!!! -self.current_render_mode_key", self.current_render_mode_key)
        if self.current_render_mode_key is None: return
        shader, blend_mode, use_depth_mask = self.render_modes[self.current_render_mode_key]
        if self.instance_data_dirty:  self.update_instance_data()
        GL.glUseProgram(self.shader_program_sphere_raytrace)
        GL.glDisable(GL.GL_BLEND) # Opaque objects don't need blending
        self.atom_instances.draw()
        
        # Transparent electrons (shader and blend mode depend on selection)
        GL.glUseProgram(shader)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendEquation(blend_mode[0])
        GL.glBlendFunc(blend_mode[1], blend_mode[2])
        self.elec_instances.draw()

        #self.render_bonds_shader()
        #self.render_bonds_simple()

        if (self.bonds is None) or (self.bonds.dirty):
            frame = self.frames_data[self.current_frame_index]
            bonds      = np.array(frame.bonds, dtype=np.uint32)
            self.bonds.upload_vbo_ebo( bonds, frame.atoms[0] )
        self.use_simple_shader()
        self.bonds.draw()

        # --- Render Atom Labels ---
        if (self.atom_labels is None) or (self.atom_labels.dirty):
            atom_positions, _, _, atom_enames = self.frames_data[self.current_frame_index].atoms
            self.atom_labels = self.make_labels( atom_positions, atom_enames )

        if self.atom_labels.nelements > 0:

            # --- DEBUG: Force GL state to a known-good configuration ---
            GL.glDisable(GL.GL_CULL_FACE)
            GL.glDisable(GL.GL_SCISSOR_TEST)
            GL.glColorMask(GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE)
            GL.glDisable(GL.GL_DEPTH_TEST)

            self.use_text_shader()
            self.atom_labels.draw()

            GL.glEnable(GL.GL_DEPTH_TEST) # <<< DEBUG: Re-enable depth test
            GL.glDisable(GL.GL_BLEND)

    def use_shader(self, shader_prog_id):
        GL.glUseProgram(shader_prog_id)
        self.current_shader_program_id = shader_prog_id # For BaseGLWidget to set uniforms

    def set_frame(self, frame_idx):
        if 0 <= frame_idx < len(self.trj):
            self.current_frame_index = frame_idx
            self.instance_data_dirty = True
            self.label_data_dirty    = True
            self.update()

    def set_opacity(self, opacity_percent):
        self.opacity = opacity_percent / 100.0
        self.instance_data_dirty = True # Need to update colors VBO
        self.update()

    def save_gl_frame_to_image(self, filename="molviewer_output.png"):
        """Renders the current scene to an offscreen FBO and saves it as a PNG image."""
        
        self.makeCurrent()

        # Ensure FBO is initialized and sized correctly for the current viewport
        if self.fbo_id is None or self.fbo_width != self.width() or self.fbo_height != self.height():
            self._init_fbo(self.width(), self.height())

        # Bind the FBO for offscreen rendering
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo_id)
        GL.glViewport(0, 0, self.fbo_width, self.fbo_height)

        # Clear the FBO's color and depth buffers
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        # Render the scene to the FBO
        self.draw_scene()

        # Read pixels from the FBO's color attachment
        pixels = GL.glReadPixels(0, 0, self.fbo_width, self.fbo_height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)

        # Unbind the FBO and restore the default viewport
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
        GL.glViewport(0, 0, self.width(), self.height())

        self.doneCurrent()

        # Create a PIL Image from the pixel data
        image = Image.frombytes("RGBA", (self.fbo_width, self.fbo_height), pixels)

        # OpenGL's glReadPixels reads from bottom-left, so flip the image vertically
        image = image.transpose(Image.FLIP_TOP_BOTTOM)

        # Save the image
        image.save(filename)
        print(f"Saved OpenGL frame to {filename}")

        # Trigger a repaint to ensure the screen is updated after FBO operations,
        # as some GL state might have been altered.
        self.update()

    def set_render_mode(self, mode_key):
        if mode_key in self.render_modes:
            self.current_render_mode_key = mode_key
            self.update() # Trigger repaint
        else:
            print(f"Warning: Render mode '{mode_key}' not found.")


class MolViewer(BaseGUI):
    def __init__(self, trj=None):
        super().__init__("Modern OpenGL Molecular Viewer")
        self.setGeometry(100, 100, 800, 600)

        main_layout = QHBoxLayout(self.main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self.gl_widget = MolViewerWidget()
        main_layout.addWidget(self.gl_widget, 1)

        self.controls_panel = QWidget()
        self.controls_layout = QVBoxLayout(self.controls_panel)
        self.controls_layout.setContentsMargins(0, 0, 0, 0)
        self.controls_layout.setSpacing(0)
        main_layout.addWidget(self.controls_panel, 0)

        self.gl_widget.trj = trj
        self.gl_widget.update()
        if self.gl_widget.trj:
            # The slider range needs to be set after trj is loaded and its length is known
            self.frame_slider = self.slider(0, 1, 100, 0, len(self.gl_widget.trj) - 1, layout=self.controls_layout, label="Frame", callback=self.gl_widget.set_frame)
        else:
            self.frame_slider = self.slider(0, 1, 100, 0, 100, layout=self.controls_layout, label="Frame", callback=self.gl_widget.set_frame) # Default range
        self.opacity_slider = self.slider(50, 1, 100, 0, 100, layout=self.controls_layout, label="Opacity", callback=self.gl_widget.set_opacity)
        self.gl_widget.precalculate_frames_data()
        self.gl_widget.update()
        self.gl_widget.main_window = self
        self.render_mode_combo = self.comboBox(items=self.gl_widget.render_modes.keys(), callback=self.on_render_mode_changed, layout=self.controls_layout)

        # Add Save Image and Quit buttons
        self.button("Save Image", callback=lambda: self.gl_widget.save_gl_frame_to_image(), layout=self.controls_layout)
        self.button("Quit", callback=self.close, layout=self.controls_layout)
        self.controls_layout.addStretch(1)

        # self.gl_widget.initializeGL()
        # self.render_mode_combo.addItems(self.gl_widget.render_modes.keys())
        # self.render_mode_combo.setCurrentText(self.gl_widget.current_render_mode_key)

        self.show()

    def on_render_mode_changed(self, text):
        self.gl_widget.set_render_mode(text)


if __name__ == '__main__':
    import os
    import sys
    import argparse
    from .. import atomicUtils as au

    # Run Like this:
    #   python -u -m pyBall.GUI.MolViewer -f ./tests/tEFF/H2O_relaxation.xyz         | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/uracil.xyz   | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/thymine.xyz  | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/adenine.xyz  | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/citosine.xyz | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/guanine.xyz  | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/CG.xyz       | tee OUT
    #   python -u -m pyBall.GUI.MolViewer -f ./cpp/common_resources/xyz/PTCDA.xyz    | tee OUT

    # relative path to the directory with the molecules
    # /home/prokop/git/FireCore/cpp/common_resources/xyz/
    dir_path = "../../cpp/common_resources/xyz/"
    dir_path = "../../cpp/common_resources/xyz_mini/"
    this_path = os.path.dirname(os.path.abspath(__file__))
    # resolve the full absolute path
    dir_path = os.path.abspath(os.path.join(this_path, dir_path))
    
    parser = argparse.ArgumentParser(description="Modern OpenGL Molecular Viewer")
    parser.add_argument("-f", "--file", type=str, help="Path to the XYZ trajectory file", default=None) # Default to None
    args = parser.parse_args()

    trj = au.load_xyz_movie(args.file)
    trj = au.trj_to_ename(trj)
    trj = au.trj_fill_radius(trj, bVdw=True, rFactor=0.001, rmin=0.05) # Adjusted rFactor for visibility
    #trj = au.trj_fill_radius(trj, bVdw=False, rFactor=1.0)
    print( "trj[0]: ", trj[0])

    app = QApplication(sys.argv)
    viewer = MolViewer(trj=trj)
    viewer.show()
    sys.exit(app.exec_())
