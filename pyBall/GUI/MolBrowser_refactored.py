from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QVBoxLayout, QHBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
import argparse
import os

sys.path.append("../../") # To find pyBall
from pyBall import elements
from .GLGUI import make_labels
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
        RvdWs = np.array( [ elements.ELEMENT_DICT[e][7] for e in self.atoms[3] ])
        bonds, bond_vecs = findBondsNP(self.atoms[0], RvdWs=RvdWs,  )
        self.bonds     = bonds
        self.bond_vecs = bond_vecs

class MoleculeRender:
    def __init__(self, parent_widget, frame_data, font_atlas_data):
        self.parent = parent_widget
        self.frame_data = frame_data

        # GL Objects
        self.atom_instances = InstancedData(base_attrib_location=1)
        self.elec_instances = InstancedData(base_attrib_location=1)
        self.atom_labels    = make_labels( frame_data.atoms[0], frame_data.atoms[3], font_atlas_data )
        self.bonds          = GLobject(components=[3], mode=GL.GL_LINES)
        
        bonds_array = np.array(frame_data.bonds, dtype=np.uint32)
        self.bonds.upload_vbo_ebo( bonds_array, frame_data.atoms[0] )

    def render(self, shader, blend_mode, use_depth_mask):
        self.parent.use_shader(shader)
        self.parent.set_default_uniforms()
        #self.atom_instances.draw()
        #self.elec_instances.draw()
        self.bonds.draw()
        self.atom_labels.draw()
    
    def cleanup(self):
        self.atom_instances.cleanup()
        self.elec_instances.cleanup()
        self.atom_labels.cleanup()
        self.bonds.cleanup()

class MolBrowserWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.molecule_renders = []
        self.render_modes = {}
        self.current_render_mode_key = None
        self.fbo_id = None
        self.fbo_texture_id = None
        self.fbo_depth_buffer_id = None
        self.fbo_width = 0
        self.fbo_height = 0
        self.font_atlas_data = None

    def initializeGL(self):
        super().initializeGL()
        font_atlas_path = os.path.join(self.shader_folder, 'font_atlas.json')
        if os.path.exists(font_atlas_path):
            with open(font_atlas_path, 'r') as f:
                self.font_atlas_data = json.load(f)
        else:
            print(f"Warning: font_atlas.json not found at {font_atlas_path}")
        shader_folder = self.shader_folder
        vert_instances = open(shader_folder + "/instances.glslv").read()
        frag_ray_sphere = open(shader_folder + "/sphere.glslf").read()
        frag_max_sphere = open(shader_folder + "/sphere_max.glslf").read()
        self.shader_program_sphere_raytrace = self.compile_shader_program(vert_instances, frag_ray_sphere)
        self.shader_program_sphere_max_vol = self.compile_shader_program(vert_instances, frag_max_sphere)
        self.bond_shader = self.compile_shader_program(open(shader_folder + "/cylinder.glslv").read(), open(shader_folder + "/cylinder.glslf").read())
        self.render_modes = {
            'pretty': (self.shader_program_sphere_raytrace, alpha_blend_modes["standard"], True),
            'fast': (self.shader_program_sphere_max_vol, alpha_blend_modes["additive"], False),
        }
        self.current_render_mode_key = 'pretty'

    def cleanupGL(self):
        for mr in self.molecule_renders:
            mr.cleanup()
        GL.glDeleteProgram(self.shader_program_sphere_raytrace)
        GL.glDeleteProgram(self.shader_program_sphere_max_vol)
        GL.glDeleteProgram(self.bond_shader)
        self._deinit_fbo()
        super().cleanupGL()

    def set_molecule(self, frame_data):
        if self.font_atlas_data is None:
            print("Error: font_atlas_data is not loaded. Cannot create molecule render.")
            return
        mr = MoleculeRender(self, frame_data, self.font_atlas_data)
        self.molecule_renders.append(mr)
        self.update()

    def use_shader(self, shader_prog_id):
        GL.glUseProgram(shader_prog_id)
        self.current_shader_program_id = shader_prog_id

    def paintGL(self):
        GL.glClearColor(0.1, 0.1, 0.1, 1.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for mr in self.molecule_renders:
            mr.render(self.shader_program_sphere_raytrace, None, True)

    def _init_fbo(self, width, height):
        self._deinit_fbo()
        self.fbo_width = width
        self.fbo_height = height
        self.fbo_id = GL.glGenFramebuffers(1)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fbo_id)
        self.fbo_texture_id = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.fbo_texture_id)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, width, height, 0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, None)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D, self.fbo_texture_id, 0)
        self.fbo_depth_buffer_id = GL.glGenRenderbuffers(1)
        GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, self.fbo_depth_buffer_id)
        GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT24, width, height)
        GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT, GL.GL_RENDERBUFFER, self.fbo_depth_buffer_id)
        if GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER) != GL.GL_FRAMEBUFFER_COMPLETE:
            print("Error: Framebuffer is not complete!")
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)

    def _deinit_fbo(self):
        if self.fbo_id: GL.glDeleteFramebuffers(1, [self.fbo_id])
        if self.fbo_texture_id: GL.glDeleteTextures(1, [self.fbo_texture_id])
        if self.fbo_depth_buffer_id: GL.glDeleteRenderbuffers(1, [self.fbo_depth_buffer_id])
        self.fbo_id = self.fbo_texture_id = self.fbo_depth_buffer_id = None

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


class MolBrowser(BaseGUI):
    def __init__(self):
        super().__init__("Molecule Browser")
        main_layout = QHBoxLayout(self.main_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self.gl_widget = MolBrowserWidget()
        main_layout.addWidget(self.gl_widget, 1)

        self.controls_panel = QWidget()
        self.controls_layout = QVBoxLayout(self.controls_panel)
        self.controls_layout.setContentsMargins(0, 0, 0, 0)
        self.controls_layout.setSpacing(0)
        main_layout.addWidget(self.controls_panel, 0)

        self.button("Save Image", callback=lambda: self.gl_widget.save_gl_frame_to_image(), layout=self.controls_layout)
        self.button("Quit", callback=self.close, layout=self.controls_layout)

        self.render_mode_combo = self.comboBox(items=self.gl_widget.render_modes.keys(), callback=self.on_render_mode_changed, layout=self.controls_layout)
        if self.gl_widget.current_render_mode_key:
            self.render_mode_combo.setCurrentText(self.gl_widget.current_render_mode_key)

        self.show()

    def on_render_mode_changed(self, text):
        self.gl_widget.set_render_mode(text)

    def cleanupGL(self):
        self.gl_widget.cleanupGL()

    def render_directory_to_images(self, directory_path, output_dir=None):
        from .. import atomicUtils as au
        if not os.path.isdir(directory_path):
            print(f"Error: Provided path '{directory_path}' is not a directory.")
            return

        if output_dir is None:
            output_dir = directory_path
        os.makedirs(output_dir, exist_ok=True)

        xyz_files = [f for f in os.listdir(directory_path) if f.endswith('.xyz')]
        print(f"Rendering {len(xyz_files)} molecules from {directory_path}")

        for filename in xyz_files:
            filepath = os.path.join(directory_path, filename)
            print(f"Processing {filepath}...")

            trj = au.load_xyz_movie(filepath)
            if not trj:
                print(f"Warning: Could not load {filepath}")
                continue

            # Create FrameData for the first frame
            frame_data = FrameData()
            atoms = trj[0]
            apos = np.array([a[1] for a in atoms], dtype=np.float32)
            enames = [a[0] for a in atoms]
            colors = [elements.ELEMENT_DICT.get(e, elements.ELEMENT_DICT['H'])[8] for e in enames]
            radii = np.array([elements.ELEMENT_DICT.get(e, elements.ELEMENT_DICT['H'])[6] for e in enames], dtype=np.float32)

            frame_data.atoms = (apos, radii, colors, enames)
            frame_data.compute_bonds()

            # Set molecule and render
            self.gl_widget.set_molecule(frame_data)
            
            output_filename = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}.png")
            self.gl_widget.save_gl_frame_to_image(output_filename)

        print("Batch rendering finished.")


def process_trj_to_frame_data(trj):
    if not trj:
        return None
    frame_data = FrameData()
    frame = trj[0]  # Get the first frame
    enames, apos, _, _, _ = frame
    colors = [elements.ELEMENT_DICT.get(e, elements.ELEMENT_DICT['H'])[8] for e in enames]
    radii = np.array([elements.ELEMENT_DICT.get(e, elements.ELEMENT_DICT['H'])[6] for e in enames], dtype=np.float32)
    frame_data.atoms = (apos, radii, colors, enames)
    frame_data.compute_bonds()
    return frame_data

if __name__ == '__main__':
    from .. import atomicUtils as au
    parser = argparse.ArgumentParser(description="Molecular Browser")
    parser.add_argument("-i", "--input", type=str, help="Path to an XYZ file or a directory of XYZ files", required=True)
    args = parser.parse_args()

    app = QApplication(sys.argv)
    browser = MolBrowser()

    all_frames_data = []
    if os.path.isdir(args.input):
        xyz_files = [f for f in os.listdir(args.input) if f.endswith('.xyz')]
        for filename in xyz_files:
            filepath = os.path.join(args.input, filename)
            trj = au.load_xyz_movie(filepath)
            frame_data = process_trj_to_frame_data(trj)
            if frame_data:
                all_frames_data.append(frame_data)
    elif os.path.isfile(args.input) and args.input.endswith('.xyz'):
        trj = au.load_xyz_movie(args.input)
        frame_data = process_trj_to_frame_data(trj)
        if frame_data:
            all_frames_data.append(frame_data)

    for frame_data in all_frames_data:
        browser.gl_widget.set_molecule(frame_data)

    browser.show()
    sys.exit(app.exec_())
