import sys
import numpy as np
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
        self.zoom = -5.0
        self.x_rot = 0
        self.y_rot = 0
        self.last_pos = None

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.0, 0.0, 0.0, 1.0)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, self.zoom)
        glRotatef(self.x_rot, 1, 0, 0)
        glRotatef(self.y_rot, 0, 1, 0)

        if self.traj and 0 <= self.current_frame < len(self.traj):
            es, apos, _, _, _ = self.traj[self.current_frame]
            for i, (atom, pos) in enumerate(zip(es, apos)):
                radius = 0.5  # Example radius, you might want to vary by atom type
                color = self.get_color(atom) + (self.opacity,)
                self.draw_sphere(pos, radius, color)

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
            dx = event.x() - self.last_pos.x()
            dy = event.y() - self.last_pos.y()
            self.x_rot += dy * 0.5
            self.y_rot += dx * 0.5
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

    def draw_sphere(self, pos, radius, color, slices=16, stacks=16):
        glColor4f(*color)
        glPushMatrix()
        glTranslatef(*pos)
        quad = gluNewQuadric()
        gluQuadricNormals(quad, GLU_SMOOTH)
        gluSphere(quad, radius, slices, stacks)
        gluDeleteQuadric(quad)
        glPopMatrix()

    def get_color(self, atom):
        # Simple color mapping based on element (you can expand this)
        if atom == 'O' or atom == '8':
            return (1.0, 0.0, 0.0)  # Red for Oxygen
        elif atom == 'H' or atom == '1':
            return (1.0, 1.0, 1.0)  # White for Hydrogen
        else:
            return (0.5, 0.5, 0.5)  # Gray for others

class MainWindow(QMainWindow):
    def __init__(self, filepath):
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
        
        self.filepath = filepath
        if self.filepath:
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


    # Configure OpenGL
    format = QSurfaceFormat()
    format.setDepthBufferSize(24)
    format.setStencilBufferSize(8)
    format.setVersion(2, 1)
    format.setProfile(QSurfaceFormat.CoreProfile)
    QSurfaceFormat.setDefaultFormat(format)

    main_window = MainWindow(filepath=args.file)
    main_window.show()
    sys.exit(app.exec_())