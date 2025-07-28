import sys
import numpy as np

from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D

from OpenGL.GL import glClear, glClearColor, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, glEnable, GL_DEPTH_TEST, glViewport, glUseProgram, glGetUniformLocation, glUniformMatrix4fv, GL_FALSE

from .OGLsystem import GLobject, Mesh, InstancedData, compile_shader_program, create_sphere_mesh # Import necessary OpenGL utilities
from .OCLsystem import OCLSystem # Import OpenCL system

class GLCLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.ogl_system = None # Will be initialized later, or passed in
        self.ocl_system = None # Will be initialized later, or passed in
        self.shader_program = None
        self.view_matrix = QMatrix4x4()
        self.projection_matrix = QMatrix4x4()
        self.model_matrix = QMatrix4x4()
        self.camera_pos = QVector3D(0, 0, 10)
        self.light_pos = QVector3D(15.0, 15.0, 30.0)
        self.light_color = QVector3D(1.0, 1.0, 1.0)

        self.setFocusPolicy(Qt.StrongFocus)
        self.gl_vbo = None
        self.vao = None
        self.particle_count = 0
        self.positions = None

    def set_shader_sources(self, name, vertex_src, fragment_src):
        self.shader_name = name
        self.vertex_shader_src = vertex_src
        self.fragment_shader_src = fragment_src

    def set_shader_program(self, program):
        self.shader_program = program

    def set_particle_data(self, positions, particle_count):
        self.positions = positions
        self.particle_count = particle_count
        self.update() # Trigger a repaint, setup_particle_vbo will be called in initializeGL or paintGL for updates

    def setup_particle_vbo(self):
        from OpenGL.GL import glGenBuffers, glBindBuffer, glBufferData, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW, glGenVertexArrays, glBindVertexArray, glVertexAttribPointer, glEnableVertexAttribArray, GL_FLOAT, GL_FALSE

        if self.gl_vbo is None:
            self.gl_vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferData(GL_ARRAY_BUFFER, self.positions.nbytes, self.positions, GL_DYNAMIC_DRAW)

        if self.vao is None:
            self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)
        glBindVertexArray(0)

    def update_particle_vbo(self, new_positions):
        from OpenGL.GL import glBindBuffer, glBufferSubData, GL_ARRAY_BUFFER
        self.positions = new_positions
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferSubData(GL_ARRAY_BUFFER, 0, self.positions.nbytes, self.positions)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

    def initializeGL(self):
        from OpenGL.GL import glClearColor, glEnable, GL_DEPTH_TEST
        glClearColor(0.1, 0.1, 0.1, 1.0)
        glEnable(GL_DEPTH_TEST)

        # Initialize OGLSystem and OCLSystem here or pass them in
        if self.ogl_system is None:
            print("Warning: OGLSystem not provided. Some functionality may be limited.")
        
        if self.ocl_system is None:
            print("Warning: OCLSystem not provided. Some functionality may be limited.")

        # Load shaders after OpenGL context is ready
        if hasattr(self, 'vertex_shader_src') and hasattr(self, 'fragment_shader_src'):
            print(f"Loading shader: {self.shader_name}")
            self.ogl_system.load_shader_program(self.shader_name, self.vertex_shader_src, self.fragment_shader_src)
            self.shader_program = self.ogl_system.get_shader_program(self.shader_name)
        else:
            print("Warning: Shader sources not provided to GLCLWidget.")

        # Initial camera and projection setup
        self.update_matrices()

        # Setup particle VBO here, after OpenGL context is ready
        if self.positions is not None and self.particle_count > 0:
            self.setup_particle_vbo()

    def paintGL(self):
        from OpenGL.GL import glClear, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, glUseProgram, glBindVertexArray, glEnable, GL_PROGRAM_POINT_SIZE, glDrawArrays, GL_POINTS

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        if self.shader_program and self.vao and self.particle_count > 0:
            glUseProgram(self.shader_program)
            self.update_matrices()
            glBindVertexArray(self.vao)
            glEnable(GL_PROGRAM_POINT_SIZE)
            glDrawArrays(GL_POINTS, 0, self.particle_count)
            glBindVertexArray(0)

    def resizeGL(self, width, height):
        from OpenGL.GL import glViewport
        glViewport(0, 0, width, height)
        self.update_matrices()

    def update_matrices(self):
        from PyQt5.QtGui import QMatrix4x4, QVector3D
        from OpenGL.GL import glGetUniformLocation, glUniformMatrix4fv, GL_FALSE

        # Model matrix (identity for now, or set by simulation)
        self.model_matrix.setToIdentity()

        # View matrix
        self.view_matrix.setToIdentity()
        self.view_matrix.lookAt(self.camera_pos, QVector3D(0, 0, 0), QVector3D(0, 1, 0))

        # Projection matrix
        self.projection_matrix.setToIdentity()
        aspect_ratio = self.width() / self.height()
        self.projection_matrix.perspective(45.0, aspect_ratio, 0.1, 100.0)

        if self.shader_program:
            model_loc = glGetUniformLocation(self.shader_program, "model")
            view_loc = glGetUniformLocation(self.shader_program, "view")
            proj_loc = glGetUniformLocation(self.shader_program, "projection")

            if model_loc != -1: glUniformMatrix4fv(model_loc, 1, GL_FALSE, self.model_matrix.data())
            if view_loc != -1: glUniformMatrix4fv(view_loc, 1, GL_FALSE, self.view_matrix.data())
            if proj_loc != -1: glUniformMatrix4fv(proj_loc, 1, GL_FALSE, self.projection_matrix.data())

    def set_systems(self, ogl_system, ocl_system):
        self.ogl_system = ogl_system
        self.ocl_system = ocl_system
