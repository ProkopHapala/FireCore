from re import S
import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
# from PyQt5.QtGui     import QSurfaceFormat
import argparse

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes

from OpenGL.GL import glUseProgram, glGetUniformLocation
#     glDisable, glEnable, GL_BLEND,
#     glBlendEquation, glBlendFunc, 
#     GL_FUNC_ADD, GL_FUNC_REVERSE_SUBTRACT, GL_MIN, GL_MAX,
#     GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE
# )

# class SphereInstancesData:
#     def __init__(self, positions=[], radii=[], colors=[]):
#         self.positions  = positions
#         self.radii      = radii
#         self.colors_rgb = colors
#         #self.nu         = len(self.positions)

#     def add_instance(self, position, radius, color):
#         self.positions .append(position)
#         self.radii     .append(radius)
#         self.colors_rgb.append(color)

#     def to_array(self, dtype=np.float32):
#         self.positions  = np.array(self.positions,  dtype=dtype)
#         self.radii      = np.array(self.radii,      dtype=dtype)
#         self.colors_rgb = np.array(self.colors_rgb, dtype=dtype)

class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.trj = []
        self.current_frame_index = 0
        self.opacity = 0.5
        self.frames_atoms = []
        self.frames_elecs = []
        self.blend_mode = alpha_blend_modes["standard"] # Default blend mode for transparent spheres
        self.shader_program_sphere_raytrace = None
        self.shader_program_sphere_max_vol = None
        self.instance_data_dirty = True
        self.electron_names = set(["E", "e", "e2", "e+", "e-"])

    @property
    def all_shader_programs(self):
        programs = []
        if self.shader_program_sphere_raytrace: programs.append(self.shader_program_sphere_raytrace)
        if self.shader_program_sphere_max_vol:  programs.append(self.shader_program_sphere_max_vol)
        return programs

    def initializeGL(self):
        # Call base class initialization for shaders, basic GL setup
        shader_folder="../../pyBall/GUI/shaders"
        VERTEX_SHADER_SOURCE   = open(shader_folder + "/instances.glslv").read()
        FRAGMENT_SHADER_SOURCE = open(shader_folder + "/sphere.glslf").read()
        # Use sphere_simple.glslf for FRAGMENT_SHADER_MAX_VOL_SOURCE for this diagnostic
        FRAGMENT_SHADER_MAX_VOL_SOURCE = open(shader_folder + "/sphere_max.glslf").read() 
        #FRAGMENT_SHADER_MAX_VOL_SOURCE = open(shader_folder + "/sphere_simple.glslf").read() 
        self.shader_program_sphere_raytrace = self.compile_shader_program( VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE)
        self.shader_program_sphere_max_vol  = self.compile_shader_program( VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_MAX_VOL_SOURCE)

        # --- Diagnostic: Check uniform locations ---
        print(f"MolViewerWidget.initializeGL: shader_program_sphere_raytrace ID = {self.shader_program_sphere_raytrace}")
        if self.shader_program_sphere_raytrace:
            glUseProgram(self.shader_program_sphere_raytrace) # Must use program to query active uniforms
            print(f"  Raytrace - model loc: {glGetUniformLocation(self.shader_program_sphere_raytrace, 'model')}")
            print(f"  Raytrace - view loc: {glGetUniformLocation(self.shader_program_sphere_raytrace, 'view')}")
            print(f"  Raytrace - projection loc: {glGetUniformLocation(self.shader_program_sphere_raytrace, 'projection')}")

        print(f"MolViewerWidget.initializeGL: shader_program_sphere_max_vol ID =  {self.shader_program_sphere_max_vol}")
        if self.shader_program_sphere_max_vol:
            glUseProgram(self.shader_program_sphere_max_vol) # Must use program
            print(f"  MaxVol/Simple - model loc: {glGetUniformLocation(self.shader_program_sphere_max_vol, 'model')}")
            print(f"  MaxVol/Simple - view loc: {glGetUniformLocation(self.shader_program_sphere_max_vol, 'view')}")
            print(f"  MaxVol/Simple - projection loc: {glGetUniformLocation(self.shader_program_sphere_max_vol, 'projection')}")
        glUseProgram(0) # Unbind
        # --- End Diagnostic ---

        print(f"MolViewerWidget.initializeGL: shader_program_sphere_raytrace ID = {self.shader_program_sphere_raytrace}")
        print(f"MolViewerWidget.initializeGL: shader_program_sphere_max_vol ID =  {self.shader_program_sphere_max_vol}")

        # Set the default shader program for BaseGLWidget's uniform setup
        self.shader_program = self.shader_program_sphere_raytrace
        if not self.shader_program: # Fallback if primary fails
            self.shader_program = self.shader_program_sphere_max_vol

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
            if frame_idx == 0:
                print("arad",   self.frames_atoms[-1][1])
                print("apos\n", self.frames_atoms[-1][0])
                print("acol\n", self.frames_atoms[-1][2])

                print("erad",   self.frames_elecs[-1][1])
                print("epos\n", self.frames_elecs[-1][0])
                print("ecol\n", self.frames_elecs[-1][2])

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
        if self.instance_data_dirty:
            self.update_instance_data()
        self.use_shader(self.shader_program_sphere_raytrace)
        disable_blend() # Opaque objects don't need blending
        self.atom_instances.draw()
        
        self.use_shader(self.shader_program_sphere_max_vol)
        set_ogl_blend_mode(self.blend_mode)
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

        self.blend_mode_combo = QComboBox()
        self.blend_mode_combo.addItems(alpha_blend_modes.keys())
        self.blend_mode_combo.setCurrentText("subtractive")
        self.gl_widget.blend_mode = alpha_blend_modes[self.blend_mode_combo.currentText()]
        
        self.blend_mode_combo.currentTextChanged.connect(self.on_blend_mode_changed)
        layout.addWidget(QLabel("Transparent Blend Mode:"))
        layout.addWidget(self.blend_mode_combo)
        
        self.gl_widget.trj = trj
        self.gl_widget.update()
        if self.gl_widget.trj:
            self.frame_slider.setMinimum(0)
            self.frame_slider.setMaximum(len(self.gl_widget.trj) - 1)

        self.gl_widget.precalculate_frames_data()
        self.gl_widget.update()
        self.show()

    def on_blend_mode_changed(self, text):
        self.gl_widget.blend_mode = alpha_blend_modes[text]
        self.gl_widget.update()


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
