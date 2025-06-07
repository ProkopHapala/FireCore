import sys
import numpy as np
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, QSlider, QLabel)
from PyQt5.QtCore    import Qt
# from PyQt5.QtGui     import QSurfaceFormat
import argparse

sys.path.append("../../") # To find pyBall
from pyBall import elements
from pyBall.GUI.GLGUI import VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE, InstancedData, BaseGLWidget, AppWindow


class MolViewerWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.trj = []
        self.current_frame_index = 0
        self.opacity = 0.5
        self.sphere_instances = None # Manages VBOs for instance attributes
        self.instance_data_dirty = True

    def initializeGL(self):
        # Call base class initialization for shaders, basic GL setup
        super().initializeGL_base(VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE)
        self.sphere_instances = InstancedData(base_attrib_location=2)
        sphere_mesh=self.default_sphere_mesh
        self.sphere_instances.associate_mesh(sphere_mesh)
        self.sphere_instances.setup_instance_vbos([
            ("positions", 0, 3),
            ("radii",     1, 1),
            ("colors",    2, 4),
        ])

    def cleanupGL(self):
        super().cleanupGL_base() # Clean up shader program
        if self.sphere_instances:
            self.sphere_instances.cleanup()

    def update_instance_data(self):
        if not self.trj or not (0 <= self.current_frame_index < len(self.trj)):
            print("Invalid trajectory or frame index")
            return

        es, apos, qs, rs, comment = self.trj[self.current_frame_index]
        na = len(es)
        positions  = np.zeros((na, 3), dtype=np.float32)
        radii      = np.zeros( na,     dtype=np.float32)
        colors     = np.zeros((na, 4), dtype=np.float32)
        for i in range(na):
            positions[i] = apos[i]
            radii[i]     = rs[i] 
            atom_symbol  = es[i]
            opacity  = self.opacity if(atom_symbol == "E") else 1.0
            try:
                element_data  = elements.ELEMENT_DICT[atom_symbol]
                hex_color_str = element_data[elements.index_color]
                r_col, g_col, b_col = elements.hex_to_float_rgb(hex_color_str)
            except KeyError:
                r_col, g_col, b_col = 0.5, 0.5, 0.5 # Gray for unknowns
            colors[i,0] = r_col
            colors[i,1] = g_col
            colors[i,2] = b_col
            colors[i,3] = opacity

        self.sphere_instances.update({"positions": positions, "radii": radii, "colors": colors}, na)
        self.instance_data_dirty = False

    def paintGL(self):
        # Calls paintGL_base, which in turn calls draw_scene
        super().paintGL_base() 

    def draw_scene(self):
        if self.instance_data_dirty:
            self.update_instance_data()
        self.sphere_instances.draw() # Pass the number of instances to draw

    # --- Trajectory specific methods ---
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
        
        self.gl_widget.trj = trj
        self.gl_widget.update()
        # if filepath is not None:
        #     self.filepath = filepath
        #     self.gl_widget.load_trajectory(self.filepath)
        if self.gl_widget.trj:
            self.frame_slider.setMinimum(0)
            self.frame_slider.setMaximum(len(self.gl_widget.trj) - 1)
        self.show()


if __name__ == '__main__':
    from pyBall import atomicUtils as au
    #app = QApplication(sys.argv)
    parser = argparse.ArgumentParser(description="Modern OpenGL Molecular Viewer")
    parser.add_argument("-f", "--file", type=str, help="Path to the XYZ trajectory file", default=None) # Default to None
    args = parser.parse_args()

    trj = au.load_xyz_movie(args.file)
    trj = au.trj_to_ename(trj)
    trj = au.trj_fill_radius(trj, bVdw=True, rFactor=0.005, rmin=0.1)
    #trj = au.trj_fill_radius(trj, bVdw=False, rFactor=1.0)
    print( "trj.enames", trj[0])
    MolViewer.launch(trj=trj)

    #main_window = MolViewer( trj=trj)
    #sys.exit(app.exec_())
