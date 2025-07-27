import os
import sys
import numpy as np
import ctypes
from OpenGL.GL import *
from OpenGL.GL.shaders import compileShader, compileProgram
from OpenGL.GLU import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QOpenGLWidget
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QSurfaceFormat, QMatrix4x4, QVector3D

#from .GLGUI import upload_buffer

def upload_buffer( index, buffer_id, data, mode=GL_DYNAMIC_DRAW):
    glBindBuffer    (GL_SHADER_STORAGE_BUFFER, buffer_id)
    glBufferData    (GL_SHADER_STORAGE_BUFFER, data.nbytes, data, mode)
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, index, buffer_id)
    glBindBuffer    (GL_SHADER_STORAGE_BUFFER, 0)


# --- Compute Shader GLSL Code ---
COMPUTE_SHADER_SOURCE = """
#version 430 core

layout (local_size_x = 32) in;  // Smaller workgroup size for better performance

layout(std430, binding = 0) buffer Positions  { vec4 positions[];  };
layout(std430, binding = 1) buffer Velocities { vec4 velocities[]; };
layout(std430, binding = 2) buffer Masses     { float masses[];    };

uniform float dt;
uniform float G; // gravitational constant
uniform int particleCount;

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= particleCount) return;
    
    vec3  pos     = positions[i].xyz;
    vec3  vel     = velocities[i].xyz;
    float mi      = masses[i];
    vec3  force   = vec3(0.0);
    
    // Calculate forces from all other particles
    for (int j = 0; j < particleCount; j++) {
        if (i == j) continue;
        vec3  d    = positions[j].xyz - pos;
        float r2   = dot(d, d) + 0.01;
        float ir   = inversesqrt(r2);
        float ir3  = ir * ir * ir;
        float mj   = masses[j];
        force     += G * mi * mj * d * ir3;
    }
    
    // Update velocity and position
    vel += (force / mi) * dt;
    pos += vel * dt;
    
    positions[i].xyz = pos;
    velocities[i].xyz = vel;
}
"""

# Vertex shader for rendering particles
VERTEX_SHADER_SOURCE = """
#version 430 core

layout(location = 0) in vec4 position;
uniform mat4 mvp;

void main() {
    gl_Position = mvp * vec4(position.xyz, 1.0);
    gl_PointSize = 2.0;
}
"""

# Fragment shader for rendering particles
FRAGMENT_SHADER_SOURCE = """
#version 430 core

out vec4 FragColor;

void main() {
    FragColor = vec4(1.0, 1.0, 1.0, 1.0); // White color
}
"""

class MyOpenGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.compute_program = None
        self.render_program  = None
        self.ssbo_ids        = None  # Will be initialized in _init_compute_shader
        self.vao             = 0
        
        # Camera parameters
        self.mvp_matrix      = QMatrix4x4()
        self.camera_pos      = [0, 0, 3]
        self.camera_target   = [0, 0, 0]
        self.camera_up       = [0, 1, 0]
        self.fov             = 45.0
        
        # N-body simulation parameters
        self.particle_count  = 1024
        self.dt              = 0.0001
        self.G               = 1.0
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._update_simulation)
        self.timer.start(16)

    def initializeGL(self):
        print(f"OpenGL Context: {glGetString(GL_VERSION).decode()}")
        
        # Initialize separate arrays for positions, velocities, and masses
        self.positions = np.zeros((self.particle_count, 4), dtype=np.float32)
        self.velocities = np.zeros((self.particle_count, 4), dtype=np.float32)
        self.masses = np.ones(self.particle_count, dtype=np.float32)
        
        # Initialize particles in a 2D spiral pattern
        center_x, center_y = 0.0, 0.0
        for i in range(self.particle_count):
            angle  = 2.0 * np.pi * i / self.particle_count
            radius = 0.5 + 0.5 * (i / self.particle_count)
            x      = center_x + radius * np.cos(angle)
            y      = center_y + radius * np.sin(angle)
            z      = 0.0
            
            # Give particles initial orbital velocities
            vx     = -radius * np.sin(angle) * 0.1
            vy     =  radius * np.cos(angle) * 0.1
            vz     =  0.0
            
            self.positions[i]  = [x, y, z, 1.0]
            self.velocities[i] = [vx, vy, vz, 0.0]
            self.masses[i]     = 1.0 + 0.5 * np.random.random()
        
        # Increase gravitational constant for more visible movement
        self.init_compute_shader()  # Initialize compute shader first
        self.init_render_pipeline()  # Initialize render pipeline
        self.perform_compute_operation()  # Run initial compute operation
        
    def init_compute_shader(self):
        print("init_compute_shader()")
        
        # Compile compute shader
        compute_shader       = compileShader(COMPUTE_SHADER_SOURCE, GL_COMPUTE_SHADER)
        self.compute_program = compileProgram(compute_shader)
        
        # Create separate SSBOs for positions, velocities, and masses
        self.position_ssbo = glGenBuffers(1)
        self.velocity_ssbo = glGenBuffers(1)
        self.mass_ssbo     = glGenBuffers(1)
        
        upload_buffer(0, self.position_ssbo, self.positions)
        upload_buffer(1, self.velocity_ssbo, self.velocities)
        upload_buffer(2, self.mass_ssbo, self.masses)
        
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0)
        
        self.ssbo_ids = [self.position_ssbo, self.velocity_ssbo, self.mass_ssbo]

        self.setup_compute_shader()
        
        return self.ssbo_ids
        
    def init_render_pipeline(self):
        print("init_render_pipeline() ")
        
        # Compile shaders
        vertex_shader       = compileShader(VERTEX_SHADER_SOURCE, GL_VERTEX_SHADER)
        fragment_shader     = compileShader(FRAGMENT_SHADER_SOURCE, GL_FRAGMENT_SHADER)
        self.render_program = compileProgram(vertex_shader, fragment_shader)
        
        # Create VAO
        self.vao = glGenVertexArrays(1)
        
        # Setup VAO to use the position buffer for rendering
        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.position_ssbo)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(0))
        glEnableVertexAttribArray(0)
        glBindVertexArray(0)
        
        # Setup MVP matrix
        self.setup_projection()

    def setup_compute_shader(self):
        self.makeCurrent()
        glUseProgram(self.compute_program)
        
        # --- Set uniforms
        dt_loc = glGetUniformLocation(self.compute_program, "dt")
        if dt_loc != -1: glUniform1f(dt_loc, self.dt)
        
        G_loc = glGetUniformLocation(self.compute_program, "G")
        if G_loc != -1: glUniform1f(G_loc, self.G)
        
        count_loc = glGetUniformLocation(self.compute_program, "particleCount")
        if count_loc != -1: glUniform1i(count_loc, self.particle_count)

        glUseProgram(0)

    def perform_compute_operation(self):
        #print(f"perform_compute_operation() particle_count={self.particle_count}")
        self.makeCurrent()
        glUseProgram(self.compute_program)        
        glDispatchCompute(int(np.ceil(self.particle_count/32)), 1, 1)  # run compute shader
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)                 # Ensure compute shader has finished
        glUseProgram(0) # Deactivate the compute program

    def update_mvp_matrix(self):
        self.mvp_matrix.setToIdentity()
        self.mvp_matrix.translate(0, 0, -5)
        self.mvp_matrix.rotate(45, QVector3D(1, 1, 0))

    def setup_projection(self):
        """Setup orthographic projection for 2D rendering"""
        self.mvp_matrix = QMatrix4x4()
        self.mvp_matrix.setToIdentity()
        self.mvp_matrix.ortho(-2.0, 2.0, -2.0, 2.0, -1.0, 100.0)
        
    def paintGL(self):
        # Clear screen
        glClearColor(0.1, 0.1, 0.1, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        if not self.render_program or not self.vao:
            raise Exception("ERROR: Missing render program or VAO")
            
        glUseProgram(self.render_program)
        mvp_loc = glGetUniformLocation(self.render_program, "mvp")
        if mvp_loc != -1: glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, self.mvp_matrix.data())
        
        # Draw particles
        glBindVertexArray(self.vao)
        glPointSize(2.0)
        glDrawArrays(GL_POINTS, 0, self.particle_count)
        glBindVertexArray(0)
        glUseProgram(0)

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)

    def closeEvent(self, event):
        print("Cleaning up OpenGL resources...")
        if self.compute_program: glDeleteProgram(self.compute_program)
        if self.render_program:  glDeleteProgram(self.render_program)
        if self.ssbo_ids:        glDeleteBuffers(len(self.ssbo_ids), self.ssbo_ids) # Check if buffers were generated
        if self.vao:             glDeleteVertexArrays(1, [self.vao]) # Check if VAO was generated 
        super().closeEvent(event)

    def _update_simulation(self):
        self.perform_compute_operation()
        self.update()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("N-Body Simulation with Compute Shaders")
        self.setGeometry(100, 100, 800, 600)

        # Request an OpenGL 4.3 Core Profile context
        fmt = QSurfaceFormat()
        fmt.setVersion(4, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        QSurfaceFormat.setDefaultFormat(fmt) # Apply the format as default for all new contexts

        self.gl_widget = MyOpenGLWidget(self)
        self.setCentralWidget(self.gl_widget)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())