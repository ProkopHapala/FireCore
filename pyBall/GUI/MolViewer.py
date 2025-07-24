from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import ( QVBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
import argparse
import os

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes

from OpenGL.GL import (
    glUseProgram, glGenVertexArrays, glBindVertexArray, glGenBuffers, glBindBuffer,
    glBufferData, glVertexAttribPointer, glEnableVertexAttribArray, glDrawArrays, glDrawElements,
    GL_FLOAT, GL_FALSE, GL_DYNAMIC_DRAW, GL_STATIC_DRAW, GL_TRIANGLES, GL_UNSIGNED_INT,
    glDeleteVertexArrays, glDeleteBuffers, glUniformMatrix4fv, glGetUniformLocation,
    glUniform1f, glUniform4f, GL_ARRAY_BUFFER, GL_ELEMENT_ARRAY_BUFFER
)
import json
import ctypes

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
        self.text_pos_vbo = None
        self.text_quad_data_vbo = None
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

        self.text_vao = glGenVertexArrays(1)
        glBindVertexArray(self.text_vao)

        self.text_pos_vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.text_pos_vbo)
        glEnableVertexAttribArray(0) # aPos3D
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)

        self.text_quad_data_vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.text_quad_data_vbo)
        stride = 4 * 4 # 4 floats (pos.xy, uv.xy)
        glEnableVertexAttribArray(1) # aLocalOffset
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(0))
        glEnableVertexAttribArray(2) # aTexCoord
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(2 * 4))

        self.text_ebo = glGenBuffers(1)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.text_ebo)

        glBindVertexArray(0)
    
    def update_atom_labels_data(self):
        if not self.trj or self.font_atlas_data is None:
            self.num_label_verts = 0
            return

        frame_atoms_data = self.frames_atoms[self.current_frame_index]
        atom_positions, _, _ = frame_atoms_data # Unpack [positions, radii, colors]

        if atom_positions.shape[0] == 0:
            self.num_label_verts = 0
            return
        atom_enames = self.trj[self.current_frame_index][2]

        all_char_base_positions_3d = []
        all_char_quad_data = []
        label_scale = 0.3  # World-space size of the label quad
        label_y_offset = 0.4 # Offset above the atom

        num_atoms = len(atom_positions)
        for i in range(num_atoms):
            pos_3d = atom_positions[i]
            element_symbol = atom_enames[i].decode('utf-8') if isinstance(atom_enames[i], bytes) else str(atom_enames[i])
            
            # Simple centering based on number of characters
            label_world_width = len(element_symbol) * label_scale * 0.5
            cursor_x_start = -label_world_width / 2.0

            for i_char, char in enumerate(element_symbol):
                if char not in self.font_atlas_data['chars']: continue
                cd = self.font_atlas_data['chars'][char]

                # Calculate offset for this specific character
                cursor_x_offset = cursor_x_start + i_char * (label_scale * 0.5)

                # Base 3D position for this character's quad
                char_base_pos = pos_3d + np.array([cursor_x_offset, label_y_offset, 0.0], dtype=np.float32)
                all_char_base_positions_3d.extend([char_base_pos] * 4) # 4 vertices for the quad

                # Local quad offsets (from character's base position)
                local_offsets = np.array([
                    [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]
                ], dtype=np.float32)

                # UVs for the character (Y is often flipped in atlases)
                uvs = np.array([
                    [cd['norm_x'], cd['norm_y'] + cd['norm_h']],
                    [cd['norm_x'] + cd['norm_w'], cd['norm_y'] + cd['norm_h']],
                    [cd['norm_x'] + cd['norm_w'], cd['norm_y']],
                    [cd['norm_x'], cd['norm_y']]
                ], dtype=np.float32)

                # Interleave local_offset and uv for each vertex
                for j in range(4):
                    all_char_quad_data.append(np.concatenate((local_offsets[j], uvs[j])))

        self.num_label_verts = len(all_char_quad_data)
        if self.num_label_verts > 0:
            # Generate and upload index data for the EBO
            num_quads = self.num_label_verts // 4
            indices = np.zeros(num_quads * 6, dtype=np.uint32)
            for i in range(num_quads):
                base = i * 4
                indices[i*6:i*6+6] = [base, base + 1, base + 2, base, base + 2, base + 3]
            
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.text_ebo)
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_DYNAMIC_DRAW)
            positions_data = np.array(all_char_base_positions_3d, dtype=np.float32)
            quad_data = np.array(all_char_quad_data, dtype=np.float32)

            glBindBuffer(GL_ARRAY_BUFFER, self.text_pos_vbo)
            glBufferData(GL_ARRAY_BUFFER, positions_data.nbytes, positions_data, GL_DYNAMIC_DRAW)

            glBindBuffer(GL_ARRAY_BUFFER, self.text_quad_data_vbo)
            glBufferData(GL_ARRAY_BUFFER, quad_data.nbytes, quad_data, GL_STATIC_DRAW)

        self.label_data_dirty = False

    def cleanupGL(self):
        super().cleanupGL()
        if self.text_shader: glDeleteProgram(self.text_shader)
        if self.text_vao: glDeleteVertexArrays(1, [self.text_vao])
        if self.text_pos_vbo: glDeleteBuffers(1, [self.text_pos_vbo])
        if self.text_quad_data_vbo: glDeleteBuffers(1, [self.text_quad_data_vbo])
        if self.text_ebo: glDeleteBuffers(1, [self.text_ebo])
        if self.font_atlas_tex: glDeleteTextures(1, [self.font_atlas_tex])
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
        #print( "update_instance_data() elecs.shape", elecs.shape )
        print( "update_instance_data() elecs:\n", elecs )
        #elecs[2][:,3] = self.opacity
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

        # --- Render Atom Labels ---
        if self.label_data_dirty:
            self.update_atom_labels_data()

        if self.num_label_verts > 0:
            self.use_shader(self.text_shader)
            self.set_default_uniforms() # Set camera matrices
            self.bind_texture('fontAtlas', self.font_atlas_tex, 0)
            set_ogl_blend_mode(alpha_blend_modes["standard"], True)
            glUniform1f(glGetUniformLocation(self.text_shader, "labelScale"), 0.3) # Adjust scale as needed
            glUniform4f(glGetUniformLocation(self.text_shader, "textColor"), 1.0, 1.0, 1.0, 1.0) # White

            glBindVertexArray(self.text_vao)
            num_indices = (self.num_label_verts // 4) * 6
            glDrawElements(GL_TRIANGLES, num_indices, GL_UNSIGNED_INT, None)
            glBindVertexArray(0)
            disable_blend()
        
    def use_shader(self, shader_prog_id):
        glUseProgram(shader_prog_id)
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
