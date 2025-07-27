import os
import sys
import numpy as np
import pyopencl as cl
from OpenGL.GL import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QOpenGLWidget
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QSurfaceFormat

# OpenCL Kernel for N-Body Simulation
OPENCL_SOURCE = """
__kernel void nbody_sim(
    __global float4* positions,
    __global float4* velocities,
    const float dt,
    const int particle_count)
{
    int i = get_global_id(0);
    if (i >= particle_count) return;

    float3 my_pos = positions[i].xyz;
    float3 my_vel = velocities[i].xyz;
    float3 force = (float3)(0.0f, 0.0f, 0.0f);
    float G = 1.0f;

    for (int j = 0; j < particle_count; j++) {
        if (i == j) continue;

        float3 other_pos = positions[j].xyz;
        float3 diff = other_pos - my_pos;
        float dist_sq = dot(diff, diff) + 0.1f; // Softening factor to avoid singularity
        float dist = sqrt(dist_sq);
        float3 direction = diff / dist;
        
        // Simplified gravity calculation for demonstration
        force += direction / dist_sq;
    }

    my_vel += force * G * dt;
    my_pos += my_vel * dt;

    positions[i].xyz = my_pos;
    velocities[i].xyz = my_vel;
}
"""

# Vertex Shader for rendering particles
VERTEX_SHADER_SOURCE = """
#version 330 core
layout (location = 0) in vec4 position;

void main()
{
    gl_Position = vec4(position.xyz, 1.0);
    gl_PointSize = 2.0;
}
"""

# Fragment shader for rendering particles
FRAGMENT_SHADER_SOURCE = """
#version 330 core
out vec4 FragColor;

void main()
{
    FragColor = vec4(1.0, 1.0, 1.0, 1.0); // White color
}
"""

class MyOpenGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.particle_count = 1024
        self.dt = 0.001

        # OpenCL and OpenGL objects
        self.cl_context     = None
        self.cl_queue       = None
        self.cl_program     = None
        self.gl_vbo         = None
        self.cl_gl_buffer   = None
        self.cl_velocities  = None
        self.vao            = None
        self.render_program = None
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(16)

    def initializeGL(self):
        glClearColor(0.1, 0.1, 0.1, 1.0)
        self.positions  = np.zeros((self.particle_count, 4), dtype=np.float32)
        self.velocities = np.zeros((self.particle_count, 4), dtype=np.float32)
        for i in range(self.particle_count):
            angle = 2.0 * np.pi * i / self.particle_count
            radius = 0.5 * np.random.rand() + 0.2
            self.positions[i] = [radius * np.cos(angle), radius * np.sin(angle), 0, 1.0]
        self.setup_render_pipeline()
        self.setup_opencl()

    def setup_render_pipeline(self):
        # Create VBO
        self.gl_vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferData(GL_ARRAY_BUFFER, self.positions.nbytes, self.positions, GL_DYNAMIC_DRAW)

        # Create VAO
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)
        glBindVertexArray(0)

        # Compile shaders
        vertex_shader   = self.compile_shader(VERTEX_SHADER_SOURCE, GL_VERTEX_SHADER)
        fragment_shader = self.compile_shader(FRAGMENT_SHADER_SOURCE, GL_FRAGMENT_SHADER)
        self.render_program = glCreateProgram()
        glAttachShader(self.render_program, vertex_shader)
        glAttachShader(self.render_program, fragment_shader)
        glLinkProgram (self.render_program)

    def setup_opencl(self):
        # Create regular OpenCL context (fallback from GL sharing)
        self.cl_context = cl.Context()
        self.cl_queue   = cl.CommandQueue(self.cl_context)

        # Create OpenCL buffers for positions and velocities ... Fallback to regular buffers since GL sharing isn't working
        self.cl_positions   = cl.Buffer (self.cl_context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.positions)
        self.cl_velocities  = cl.Buffer (self.cl_context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.velocities)
        self.cl_program     = cl.Program(self.cl_context, OPENCL_SOURCE).build()   # Compile OpenCL kernel
        self.use_gl_sharing = False   # Flag to indicate we're not using GL sharing

    def update_simulation(self):
        self.makeCurrent()
        if self.cl_queue is None: raise Exception("ERROR: OpenCL queue is None")
            
        # Run OpenCL kernel on position and velocity buffers
        self.cl_program.nbody_sim(
            self.cl_queue, 
            (self.particle_count,), 
            None,
            self.cl_positions, 
            self.cl_velocities,
            np.float32(self.dt), 
            np.int32(self.particle_count)
        )
        
        # Copy updated positions back to OpenGL buffer
        cl.enqueue_copy(self.cl_queue, self.positions, self.cl_positions)
        
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferSubData(GL_ARRAY_BUFFER, 0, self.positions.nbytes, self.positions) # Update VBO with new positions
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        
        self.cl_queue.finish()
        self.update()           # Trigger paintGL

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT)
        glUseProgram(self.render_program)
        glBindVertexArray(self.vao)
        glDrawArrays(GL_POINTS, 0, self.particle_count)
        glBindVertexArray(0)

    def compile_shader(self, source, shader_type):
        shader = glCreateShader(shader_type)
        glShaderSource(shader, source)
        glCompileShader(shader)
        if not glGetShaderiv(shader, GL_COMPILE_STATUS):
            raise Exception(glGetShaderInfoLog(shader).decode())
        return shader

    def closeEvent(self, event):
        self.timer.stop()
        if self.gl_vbo:         glDeleteBuffers(1, [self.gl_vbo])
        if self.vao:            glDeleteVertexArrays(1, [self.vao])
        if self.render_program: glDeleteProgram(self.render_program)
        super().closeEvent(event)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyOpenCL N-Body Simulation")
        self.setGeometry(100, 100, 800, 800)

        fmt = QSurfaceFormat()
        fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        QSurfaceFormat.setDefaultFormat(fmt)

        self.gl_widget = MyOpenGLWidget(self)
        self.setCentralWidget(self.gl_widget)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())