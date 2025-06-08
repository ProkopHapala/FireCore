import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, QSlider, QLabel, QComboBox)
from PyQt5.QtCore    import Qt
# from PyQt5.QtGui     import QSurfaceFormat
import argparse

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import BaseGLWidget, AppWindow, InstancedData, set_ogl_blend_mode, disable_blend, alpha_blend_modes


#from OpenGL.GL import (
#     glDisable, glEnable, GL_BLEND,
#     glBlendEquation, glBlendFunc, 
#     GL_FUNC_ADD, GL_FUNC_REVERSE_SUBTRACT, GL_MIN, GL_MAX,
#     GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE
# )

# def set_ogl_blend_mode(mode_str):
#     """Sets the OpenGL blend function and equation based on a mode string."""
#     if mode_str == "standard":
#         glBlendEquation(GL_FUNC_ADD)
#         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
#     elif mode_str == "additive":
#         glBlendEquation(GL_FUNC_ADD)
#         glBlendFunc(GL_SRC_ALPHA, GL_ONE)
#     elif mode_str == "subtractive_alpha_one":
#         glBlendEquation(GL_FUNC_REVERSE_SUBTRACT)
#         glBlendFunc(GL_SRC_ALPHA, GL_ONE)
#     elif mode_str == "minimum" or mode_str == "maximum":
#         glBlendEquation(GL_MIN if mode_str == "minimum" else GL_MAX)
#         glBlendFunc(GL_ONE, GL_ONE) # Common for min/max, takes component-wise min/max

class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.trj = []
        self.current_frame_index = 0
        self.opacity = 0.5
        self.opaque_sphere_instances = None
        self.transparent_sphere_instances = None
        self.blend_mode = alpha_blend_modes["standard"] # Default blend mode for transparent spheres
        self.instance_data_dirty = True

    def initializeGL(self):
        # Call base class initialization for shaders, basic GL setup
        shader_folder="../../pyBall/GUI/shaders"
        VERTEX_SHADER_SOURCE   = open(shader_folder + "/instances.glslv").read()
        FRAGMENT_SHADER_SOURCE = open(shader_folder + "/sphere.glslf").read()
        super().initializeGL_base(VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE)
        sphere_mesh=self.default_sphere_mesh

        self.opaque_sphere_instances = InstancedData(base_attrib_location=2)
        self.opaque_sphere_instances.associate_mesh(sphere_mesh)
        self.opaque_sphere_instances.setup_instance_vbos([
            ("positions", 0, 3),
            ("radii",     1, 1),
            ("colors",    2, 4),
        ])

        self.transparent_sphere_instances = InstancedData(base_attrib_location=2) # Same attrib locations
        self.transparent_sphere_instances.associate_mesh(sphere_mesh)
        self.transparent_sphere_instances.setup_instance_vbos([
            ("positions", 0, 3),
            ("radii",     1, 1),
            ("colors",    2, 4),
        ])
    def cleanupGL(self):
        super().cleanupGL_base() # Clean up shader program
        if self.opaque_sphere_instances:
            self.opaque_sphere_instances.cleanup()
        if self.transparent_sphere_instances:
            self.transparent_sphere_instances.cleanup()

    def update_instance_data(self):
        if not self.trj or not (0 <= self.current_frame_index < len(self.trj)):
            print("MolViewerWidget.update_instance_data: Invalid trajectory or frame index")
            self.opaque_sphere_instances.update({}, 0)
            self.transparent_sphere_instances.update({}, 0)
            self.instance_data_dirty = False
            return

        es, apos, qs, rs, comment = self.trj[self.current_frame_index]
        na = len(es)
        print(f"MolViewerWidget.update_instance_data: Frame {self.current_frame_index}, Total particles: {na}")
        if na == 0: # Handle empty frames
            self.opaque_sphere_instances.update({}, 0)
            self.transparent_sphere_instances.update({}, 0)
            self.instance_data_dirty = False
            return

        opaque_positions_list = []
        opaque_radii_list = []
        opaque_colors_list = []

        transparent_positions_list = []
        transparent_radii_list = []
        transparent_colors_list = []

        for i in range(na):
            atom_symbol    = es[i]
            current_radius = rs[i]
            element_data  = elements.ELEMENT_DICT[atom_symbol]
            hex_color_str = element_data[elements.index_color]
            r_col, g_col, b_col = elements.hex_to_float_rgb(hex_color_str)

            if atom_symbol == "E": # Assuming "E" denotes an electron or transparent particle
                #print(f"  Transparent ({atom_symbol}): idx={i}, pos={apos[i]}, radius={current_radius}, alpha={self.opacity}")
                transparent_positions_list.append(apos[i])
                transparent_radii_list.append(current_radius)
                transparent_colors_list.append([r_col, g_col, b_col, self.opacity])
            else:
                #print(f"  Opaque      ({atom_symbol}): idx={i}, pos={apos[i]}, radius={current_radius}, alpha=1.0 color=({r_col},{g_col},{b_col})")
                opaque_positions_list.append(apos[i])
                opaque_radii_list.append(current_radius)
                opaque_colors_list.append([r_col, g_col, b_col, 1.0]) # Opaque atoms have alpha = 1.0

        num_opaque      = len(opaque_positions_list)
        num_transparent = len(transparent_positions_list)
        #print(f"MolViewerWidget.update_instance_data: num_opaque_processed={num_opaque}, num_transparent_processed={num_transparent}")

        if num_opaque > 0:
            opaque_positions = np.array(opaque_positions_list, dtype=np.float32)
            opaque_radii     = np.array(opaque_radii_list, dtype=np.float32)
            opaque_colors    = np.array(opaque_colors_list, dtype=np.float32)
            # print(f"  Opaque radii being sent to GPU: {opaque_radii.flatten()}")
            self.opaque_sphere_instances.update({"positions": opaque_positions, "radii": opaque_radii, "colors": opaque_colors}, num_opaque)
        else:
            self.opaque_sphere_instances.update({}, 0) # Clear if no opaque instances
        if num_transparent > 0:
            transparent_positions = np.array(transparent_positions_list, dtype=np.float32)
            transparent_radii     = np.array(transparent_radii_list, dtype=np.float32)
            transparent_colors    = np.array(transparent_colors_list, dtype=np.float32)
            self.transparent_sphere_instances.update({"positions": transparent_positions, "radii": transparent_radii, "colors": transparent_colors}, num_transparent)
        else:
            self.transparent_sphere_instances.update({}, 0) # Clear if no transparent instances

        self.instance_data_dirty = False

    def paintGL(self):
        # Calls paintGL_base, which in turn calls draw_scene
        super().paintGL_base() 

    def draw_scene(self):
        if self.instance_data_dirty:
            self.update_instance_data()

        # Set blend mode for opaque spheres (standard)
        # This ensures if previous draw call changed blend equation, it's reset.
        disable_blend() # Opaque objects don't need blending
        self.opaque_sphere_instances.draw()

        # # Set blend mode for transparent spheres based on selection
        set_ogl_blend_mode(self.blend_mode)
        self.transparent_sphere_instances.draw()
        # set_ogl_blend_mode(self.blend_mode)
        # self.transparent_sphere_instances.draw()
    # def load_trajectory(self, filename):
    #     self.trj = au.load_xyz_movie(filename)
    #     self.trj = au.trj_to_ename(self.trj) # Ensure element names
    #     self.trj = au.trj_fill_radius(self.trj, bVdw=False, rFactor=1.0, rmin=0.1) # Use covalent radii
    #     self.current_frame_index = 0
    #     self.instance_data_dirty = True
    #     self.update()
    #     return len(self.trj)

    def set_frame(self, frame_idx):
        if 0 <= frame_idx < len(self.trj):
            self.current_frame_index = frame_idx
            self.instance_data_dirty = True
            self.update()

    def set_opacity(self, opacity_percent):
        self.opacity = opacity_percent / 100.0
        self.instance_data_dirty = True # Need to update colors VBO
        self.update()

    # def set_blend_mode(self, mode_str):
    #     if mode_str in ["standard", "additive", "subtractive_alpha_one", "minimum", "maximum"]:
    #         self.blend_mode = mode_str
    #         self.update() # Trigger repaint


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
        self.blend_mode_combo.currentTextChanged.connect(self.on_blend_mode_changed)
        layout.addWidget(QLabel("Transparent Blend Mode:"))
        layout.addWidget(self.blend_mode_combo)
        
        self.gl_widget.trj = trj
        self.gl_widget.update()
        # if filepath is not None:
        #     self.filepath = filepath
        #     self.gl_widget.load_trajectory(self.filepath)
        if self.gl_widget.trj:
            self.frame_slider.setMinimum(0)
            self.frame_slider.setMaximum(len(self.gl_widget.trj) - 1)
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
