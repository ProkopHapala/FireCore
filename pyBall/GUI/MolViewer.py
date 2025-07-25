from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import (QVBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
import argparse
import os

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes

import OpenGL.GL as GL
import json
import ctypes

np.set_printoptions(linewidth=np.inf)


class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.trj = []
        self.current_frame_index = 0
        self.opacity      = 0.5
        self.frames_atoms = []
        self.frames_elecs = []
        self.atom_instances = InstancedData(base_attrib_location=1)
        self.elec_instances = InstancedData(base_attrib_location=1)
        # self.blend_mode = alpha_blend_modes["standard"] # Replaced by render_modes
        self.shader_program_sphere_raytrace = None
        self.shader_program_sphere_max_vol = None
        self.instance_data_dirty = True
        self.electron_names = set(["E", "e", "e2", "e+", "e-"])
        self.render_modes = {}
        self.current_render_mode_key = None

        # --- Attributes for Text Rendering
        self.text_shader = None
        self.text_vao = None
        self.text_vbo = None
        self.text_ebo = None
        self.font_atlas_tex = None
        self.font_atlas_data = None
        self.num_label_verts = 0
        self.label_data_dirty = True

    @property
    def all_shader_programs(self):
        programs = []
        if self.shader_program_sphere_raytrace: programs.append(self.shader_program_sphere_raytrace)
        if self.shader_program_sphere_max_vol:  programs.append(self.shader_program_sphere_max_vol)
        if self.text_shader: programs.append(self.text_shader)
        return programs

    def initializeGL(self):
        #print("MolViewerWidget.initializeGL()")
        # Call base class initialization for shaders, basic GL setup
        shader_folder=os.path.dirname(os.path.abspath(__file__)) + "/shaders"
        print("shader_folder", shader_folder)
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

        # --- Text Rendering Setup ---
        vert_text = open(shader_folder + "/text_billboard.glslv").read()
        frag_text = open(shader_folder + "/text_billboard.glslf").read()
        self.text_shader = self.compile_shader_program(vert_text, frag_text)

        self.font_atlas_tex = self.load_texture(shader_folder + "/font_atlas.png")
        with open(shader_folder + "/font_atlas.json") as f:
            self.font_atlas_data = json.load(f)

        self.text_vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(self.text_vao)

        self.text_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.text_vbo)

        # All attributes are now in a single interleaved VBO
        # Stride = (3_pos + 2_offset + 2_uv) * sizeof(float)
        stride = (3 + 2 + 2) * 4

        # aPos3D
        GL.glEnableVertexAttribArray(0)
        GL.glVertexAttribPointer(0, 3, GL.GL_FLOAT, GL.GL_FALSE, stride, ctypes.c_void_p(0))

        # aLocalOffset
        GL.glEnableVertexAttribArray(1)
        GL.glVertexAttribPointer(1, 2, GL.GL_FLOAT, GL.GL_FALSE, stride, ctypes.c_void_p(3 * 4))

        # aTexCoord
        GL.glEnableVertexAttribArray(2)
        GL.glVertexAttribPointer(2, 2, GL.GL_FLOAT, GL.GL_FALSE, stride, ctypes.c_void_p((3 + 2) * 4))

        self.text_ebo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.text_ebo)

        GL.glBindVertexArray(0)
    
    def update_atom_labels_data(self):
        # if not self.trj or self.font_atlas_data is None:
        #     self.num_label_verts = 0
        #     return

        atom_positions, _, _, atom_enames = self.frames_atoms[self.current_frame_index]

        print("atom_enames", atom_enames)

        tile_w = self.font_atlas_data['tile_w']  # width of one character tile in atlas
        tile_h = self.font_atlas_data['tile_h']  # height of one character tile in atlas
        tex_w  = self.font_atlas_data['tex_w']   # full atlas texture width

        all_vertex_data = []
        szx  = 1.0  # width of one character quad in world units
        szy  = 2.0*0.5  # height of one character quad in world units
        yoff = 0.0

        margin_lx = 0.1
        margin_rx = 0.7
        margin_y  = 0.0
        space     = 0.0

        conde_min = 32
        conde_max = 126

        for pos_3d, ename in zip(atom_positions, atom_enames):
            #ename += "_Hey"  # debug longer text
            symbol = ename.decode('utf-8') if isinstance(ename, bytes) else str(ename) 
            #xoff     = -len(symbol)* szx / 2.0 # centered text
            xoff     = 0 # left aligned text

            for i_char, char in enumerate(symbol):
                code = ord(char)
                if code < conde_min or code > conde_max: continue

                # 1D atlas UVs  - using normalized coordinates
                u0 = (code - conde_min + margin_lx ) * tile_w / tex_w
                u1 = (code - conde_min + margin_rx ) * tile_w / tex_w
                #u1 = u0               + tile_w*(1-margin_x) / tex_w
                v1 = 0.0 + margin_y
                v0 = 1.0 - margin_y

                # Same 3D position for every character in this label
                base_3d = pos_3d + np.array([0.0, 0.0, 0.0])

                # Local screen-space offset for this character
                x = xoff + i_char * (szx*(1.0+space))
                local_offsets = np.array([
                    [ x       , yoff-szy ],
                    [ x + szx , yoff-szy ],
                    [ x + szx , yoff+szy ],
                    [ x       , yoff+szy ]
                ], dtype=np.float32)
                uvs = np.array([[u0, v0], [u1, v0], [u1, v1], [u0, v1]], dtype=np.float32)

                for j in range(4):
                    vertex_data = np.concatenate((base_3d, local_offsets[j], uvs[j]))
                    all_vertex_data.append(vertex_data)

        self.num_label_verts = len(all_vertex_data)
        #print(f"DEBUG: update_atom_labels_data created {self.num_label_verts} vertices.")
        
        if self.num_label_verts > 0:
            # Generate and upload index data for the EBO
            num_quads = self.num_label_verts // 4
            indices = np.zeros(num_quads * 6, dtype=np.uint32)
            for i in range(num_quads):
                base = i * 4
                indices[i*6:i*6+6] = [base, base + 1, base + 2, base, base + 2, base + 3]
            
            vertex_data_np = np.array(all_vertex_data, dtype=np.float32)

            #print("vertex_data_np #legend:  pos(x,y,z)   local_offset(x,y)   uv(x,y)\n", vertex_data_np)
            #print("indices: ", indices)
            
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.text_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, vertex_data_np.nbytes, vertex_data_np, GL.GL_DYNAMIC_DRAW)

            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.text_ebo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)

        self.label_data_dirty = False

    def cleanupGL(self):
        super().cleanupGL()
        if self.text_shader: GL.glDeleteProgram(self.text_shader)
        if self.text_vao: GL.glDeleteVertexArrays(1, [self.text_vao])
        if self.text_vbo: GL.glDeleteBuffers(1, [self.text_vbo])
        if self.text_ebo: GL.glDeleteBuffers(1, [self.text_ebo])
        if self.font_atlas_tex: GL.glDeleteTextures(1, [self.font_atlas_tex])
        if self.shader_program_sphere_raytrace:
            GL.glDeleteProgram(self.shader_program_sphere_raytrace)
            self.shader_program_sphere_raytrace = None
        if self.shader_program_sphere_max_vol:
            GL.glDeleteProgram(self.shader_program_sphere_max_vol)
            self.shader_program_sphere_max_vol = None
        if self.atom_instances:
            self.atom_instances.cleanup()
        if self.elec_instances:
            self.elec_instances.cleanup()

    def precalculate_frames_data(self, dtype=np.float32):
        for frame_idx in range(len(self.trj)):
            es,ps,qs, rs, comment = self.trj[frame_idx]
            na = len(es)
            GL.glDisable(GL.GL_BLEND)
            # atoms = SphereInstancesData()
            # elecs = SphereInstancesData()
            apos,arad,acol,aname=[],[],[],[]
            epos,erad,ecol,ename=[],[],[],[]
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
                        ename.append(atom_symbol)
                    else: 
                        apos.append(ps[i])
                        arad.append(current_radius)
                        acol.append([r_col, g_col, b_col, 1.0])
                        aname.append(atom_symbol)
            
            self.frames_atoms.append( [np.array(apos,dtype=dtype), np.array(arad,dtype=dtype), np.array(acol,dtype=dtype), aname] )
            self.frames_elecs.append( [np.array(epos,dtype=dtype), np.array(erad,dtype=dtype), np.array(ecol,dtype=dtype), ename] )

    def update_instance_data(self):
        atoms = self.frames_atoms[self.current_frame_index][:3]
        elecs = self.frames_elecs[self.current_frame_index][:3]
        self.atom_instances.update_list( atoms )
        self.elec_instances.update_list( elecs )
        self.instance_data_dirty = False

    def paintGL(self):
        # Calls paintGL_base, which in turn calls draw_scene
        super().paintGL_base() 

    def draw_scene(self):
        shader, blend_mode, use_depth_mask = self.render_modes[self.current_render_mode_key]

        if self.instance_data_dirty:
            self.update_instance_data()
        GL.glUseProgram(self.shader_program_sphere_raytrace)
        GL.glDisable(GL.GL_BLEND) # Opaque objects don't need blending
        self.atom_instances.draw()
        
        # Transparent electrons (shader and blend mode depend on selection)
        GL.glUseProgram(shader)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendEquation(blend_mode[0])
        GL.glBlendFunc(blend_mode[1], blend_mode[2])
        self.elec_instances.draw()

        # --- Render Atom Labels ---
        if self.label_data_dirty:
            self.update_atom_labels_data()

        if self.num_label_verts > 0:
            #print(f"DEBUG: draw_scene attempting to draw {(self.num_label_verts // 4) * 4} vertices for {(self.num_label_verts // 4)} quads.")
            self.use_shader(self.text_shader)
            self.set_default_uniforms() # Set camera matrices

            self.bind_texture('fontAtlas', self.font_atlas_tex, 0)
            GL.glDisable(GL.GL_DEPTH_TEST) # <<< DEBUG: Disable depth test
            #set_ogl_blend_mode(alpha_blend_modes["standard"], True)
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
            GL.glUniform1f(GL.glGetUniformLocation(self.text_shader, "labelScale"), 0.25*self.zoom_factor ) # Adjust scale as needed
            GL.glUniform4f(GL.glGetUniformLocation(self.text_shader, "textColor"), 1.0, 0.0, 0.0, 1.0) # purple
            GL.glUniform2f(GL.glGetUniformLocation(self.text_shader, "offset"), 0.5, 0.0) # purple

            # --- DEBUG: Force GL state to a known-good configuration ---
            GL.glDisable(GL.GL_CULL_FACE)
            GL.glDisable(GL.GL_SCISSOR_TEST)
            GL.glColorMask(GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE)
            # -----------------------------------------------------------

            GL.glDisable(GL.GL_DEPTH_TEST) # --- DEBUG: Disable depth test
            GL.glBindVertexArray(self.text_vao)
            num_indices = (self.num_label_verts // 4) * 6
            GL.glDrawElements(GL.GL_TRIANGLES, num_indices, GL.GL_UNSIGNED_INT, None)
            GL.glBindVertexArray(0)
            GL.glEnable(GL.GL_DEPTH_TEST) # <<< DEBUG: Re-enable depth test
            #print(f"DEBUG: text_shader={self.text_shader}, fontAtlas uniform loc={GL.glGetUniformLocation(self.text_shader, 'fontAtlas')}, tex_id={self.font_atlas_tex}")
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

    def set_render_mode(self, mode_key):
        if mode_key in self.render_modes:
            self.current_render_mode_key = mode_key
            self.update() # Trigger repaint
        else:
            print(f"Warning: Render mode '{mode_key}' not found.")


class MolViewer(AppWindow):
    def __init__(self, trj=None):
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

        #self.gl_widget.initializeGL()
        #self.render_mode_combo.addItems(self.gl_widget.render_modes.keys())
        #self.render_mode_combo.setCurrentText(self.gl_widget.current_render_mode_key)

        self.show()

    def on_render_mode_changed(self, text):
        self.gl_widget.set_render_mode(text)


if __name__ == '__main__':
    import os
    import sys
    import argparse
    from .. import atomicUtils as au

    # Run Like this:
    #   python -m pyBall.GUI.MolViewer -f /tests/tEFF/H2O_relaxation.xyz | tee OUT

    # relative path to the directory with the molecules
    # /home/prokop/git/FireCore/cpp/common_resources/xyz/
    dir_path = "../../cpp/common_resources/xyz/"
    dir_path = "../../cpp/common_resources/xyz_mini/"
    this_path = os.path.dirname(os.path.abspath(__file__))
    # resolve the full absolute path
    dir_path = os.path.abspath(os.path.join(this_path, dir_path))
    
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
