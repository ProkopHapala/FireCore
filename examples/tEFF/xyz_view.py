import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget, QSlider, QLabel, QOpenGLWidget)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QSurfaceFormat
# from PyQt5.QtOpenGL import QOpenGLWidget # Incorrect import
from OpenGL.GL import *
from OpenGL.GLU import *

sys.path.append("../../")
from pyBall import atomicUtils as au  # Assuming atomicUtils.py is in the same dir

class GLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.traj = []
        self.current_frame = 0
        self.opacity = 0.5
        self.zoom = -10.0  # Adjusted default zoom for better initial view
        self.orientation = R.identity() # Store orientation as a Scipy Rotation object
        self.last_pos = None

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)
        glClearColor(1.0, 1.0, 1.0, 1.0)  # White background
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        # Lighting setup
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_NORMALIZE) # Important for correct lighting with scaling

        # Light properties
        light_position = [15.0, 15.0, 15.0, 1.0] # Positional light
        light_ambient  = [0.3, 0.3, 0.3, 1.0]
        light_diffuse  = [0.8, 0.8, 0.8, 1.0]
        light_specular = [0.9, 0.9, 0.9, 1.0] # Bright specular light

        glLightfv(GL_LIGHT0, GL_POSITION, light_position)
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)

        # Enable color material, so glColor affects material properties (ambient & diffuse)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, self.zoom)
        glMultMatrixf(rotation_to_gl_matrix(self.orientation))
        self._draw_scene()
    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(width) / height, 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.last_pos = event.pos()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.LeftButton and self.last_pos:
            p1 = self._project_to_sphere(self.last_pos.x(), self.last_pos.y())
            p2 = self._project_to_sphere(event.x(), event.y())

            rotation_axis = np.cross(p1, p2)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-6: # Avoid division by zero if no significant movement
                rotation_axis /= rotation_axis_norm
                dot_product = np.dot(p1, p2)
                rotation_angle = np.arccos(np.clip(dot_product, -1.0, 1.0))

                # Create a delta rotation from axis-angle
                delta_rotation = R.from_rotvec(rotation_angle * rotation_axis)
                self.orientation = delta_rotation * self.orientation
            self.last_pos = event.pos()
            self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom += delta * 0.01
        self.update()

    def load_trajectory(self, filename):
        self.traj = au.load_xyz_movie(filename)
        self.update()

    def set_frame(self, frame):
        self.current_frame = frame
        self.update()

    def set_opacity(self, opacity):
        self.opacity = opacity
        self.update()

    def _draw_scene(self):
        if self.traj and 0 <= self.current_frame < len(self.traj):
            es, apos, qs, rs, comment = self.traj[self.current_frame]
            for i, (atom_symbol, pos) in enumerate(zip(es, apos)):
                radius = rs[i]
                if np.isnan(radius): # Default radius if not specified or NaN
                    radius = 0.5
                # Get color from elements.py
                try:
                    element_data = au.elements.ELEMENT_DICT[atom_symbol]
                    hex_color_str = element_data[au.elements.index_color]
                    r, g, b = au.elements.hex_to_float_rgb(hex_color_str)
                except KeyError: # Fallback for unknown elements
                    r, g, b = 0.5, 0.5, 0.5 # Gray
                
                color_with_opacity = (r, g, b, self.opacity)
                self.draw_sphere(pos, radius, color_with_opacity)

    def draw_sphere(self, pos, radius, color, slices=16, stacks=16):
        # Material properties for specular reflection
        mat_specular   = [0.7, 0.7, 0.7, 1.0] # Specular color
        mat_shininess  = [60.0]                # Shininess factor (0-128)
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular)
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess)

        glColor4f(*color)
        glPushMatrix()
        glTranslatef(*pos)
        quad = gluNewQuadric()
        gluQuadricNormals(quad, GLU_SMOOTH)
        gluSphere(quad, radius, slices, stacks)
        gluDeleteQuadric(quad)
        glPopMatrix()

    def _project_to_sphere(self, x, y):
        width = self.width()
        height = self.height()
        radius = min(width, height) * 0.4 # Effective radius of the trackball

        vx = (x - width / 2.0) / radius
        vy = (height / 2.0 - y) / radius # Y is inverted

        vz_squared = 1.0 - vx * vx - vy * vy
        if vz_squared > 0:
            vz = np.sqrt(vz_squared)
        else:
            norm_xy = np.sqrt(vx*vx + vy*vy)
            if norm_xy > 1e-6: # Avoid division by zero
                # Project to the edge of the sphere if outside
                vx /= norm_xy
                vy /= norm_xy
            vz = 0.0
        
        p = np.array([vx, vy, vz])
        norm_p = np.linalg.norm(p)
        if norm_p > 1e-6:
            return p / norm_p
        return np.array([0.0, 0.0, 1.0]) # Default if at center

    def _rotation_to_gl_matrix(self, rotation_obj):
        # Scipy's as_matrix() returns a 3x3 rotation matrix
        mat3x3 = rotation_obj.as_matrix()
        # OpenGL needs a 4x4 matrix in column-major order
        mat4x4 = np.eye(4, dtype=np.float32)
        mat4x4[:3, :3] = mat3x3
        return mat4x4.flatten(order='F')

class MainWindow(QMainWindow):
    def __init__(self, filepath=None, trj=None ):
        super().__init__()
        self.setWindowTitle("Simple Molecular Viewer")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.gl_widget = GLWidget()
        layout.addWidget(self.gl_widget)

        self.frame_slider = QSlider(Qt.Horizontal)
        self.frame_slider.valueChanged.connect(self.gl_widget.set_frame)
        layout.addWidget(self.frame_slider)

        self.opacity_slider = QSlider(Qt.Horizontal)
        self.opacity_slider.setMinimum(0)
        self.opacity_slider.setMaximum(100)
        self.opacity_slider.setValue(int(self.gl_widget.opacity * 100))
        self.opacity_slider.valueChanged.connect(self.set_opacity)
        layout.addWidget(self.opacity_slider)
        
        if trj is not None:
            self.gl_widget.traj = trj
            self.gl_widget.update()
        if filepath is not None:
            self.filepath = filepath
            self.gl_widget.load_trajectory(self.filepath)
        if self.gl_widget.traj:
            self.frame_slider.setMinimum(0)
            self.frame_slider.setMaximum(len(self.gl_widget.traj) - 1)

    def set_opacity(self, value):
        opacity = value / 100.0
        self.gl_widget.set_opacity(opacity)

if __name__ == '__main__':
    app = QApplication(sys.argv)

    parser = argparse.ArgumentParser(description="Simple Molecular Viewer")
    parser.add_argument("-f", "--file", type=str, help="Path to the XYZ trajectory file", default="/home/prokophapala/git/FireCore/tests/tEFF/relaxation.xyz")
    args = parser.parse_args()

    traj = au.load_xyz_movie(args.file)
    traj = au.trj_to_ename(traj)
    #traj = au.trj_fill_radius(traj, bVdw=True, rFactor=0.5)
    traj = au.trj_fill_radius(traj, bVdw=False, rFactor=1.0)
    print( "trj.enames", traj[0])

    # Configure OpenGL
    format = QSurfaceFormat()
    format.setDepthBufferSize(24)
    format.setStencilBufferSize(8)
    format.setVersion(2, 1)
    format.setProfile(QSurfaceFormat.CoreProfile)
    QSurfaceFormat.setDefaultFormat(format)

    main_window = MainWindow(trj=traj)
    main_window.show()
    sys.exit(app.exec_())