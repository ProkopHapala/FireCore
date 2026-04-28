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

# --- How to Run This Script Correctly on a Hybrid Graphics Laptop (Linux) ---
# To ensure the high-performance NVIDIA GPU is used, you MUST launch the
# script with these special environment variables:
#
# $ __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python your_script_name.py
#
# This instructs the operating system to use the NVIDIA card for OpenGL.
# -----------------------------------------------------------------------------


def upload_buffer(index, buffer_id, data, mode=GL_DYNAMIC_DRAW):
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buffer_id)
    glBufferData(GL_SHADER_STORAGE_BUFFER, data.nbytes, data, mode)
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, index, buffer_id)
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0)


# --- Compute Shader GLSL Code ---
COMPUTE_SHADER_SOURCE = """
#version 430 core

layout (local_size_x = 256) in; // A common workgroup size

layout(std430, binding = 0) buffer Positions  { vec4 positions[];  };
layout(std430, binding = 1) buffer Velocities { vec4 velocities[]; };
layout(std430, binding = 2) buffer Masses     { float masses[];    };

uniform float dt;
uniform float G;
uniform int particleCount;

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= particleCount) return;
    
    vec3 my_pos   = positions[i].xyz;
    vec3 my_vel   = velocities[i].xyz;
    float my_mass = masses[i];
    vec3 force    = vec3(0.0);
    
    // N-Squared force calculation
    for (int j = 0; j < particleCount; j++) {
        if (i == j) continue;
        
        vec3 other_pos  = positions[j].xyz;
        vec3 d_vec      = other_pos - my_pos;
        
        // Use a softening factor to prevent extreme forces at close range
        float dist_sq   = dot(d_vec, d_vec) + 0.01; 
        float inv_dist  = inversesqrt(dist_sq);
        float inv_dist3 = inv_dist * inv_dist * inv_dist;
        
        float other_mass = masses[j];
        force += d_vec * other_mass * inv_dist3;
    }
    
    force *= G * my_mass;
    
    // Update velocity and position using Forward Euler integration
    my_vel += (force / my_mass) * dt;
    my_pos += my_vel * dt;
    
    positions[i].xyz  = my_pos;
    velocities[i].xyz = my_vel;
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
        self.position_ssbo, self.velocity_ssbo, self.mass_ssbo = 0, 0, 0
        self.vao = 0
        
        self.mvp_matrix = QMatrix4x4()
        
        self.particle_count = 2048
        self.dt = 0.0001
        self.G = 0.1 # Adjusted for more stable simulation
        
        self.timer = QTimer(self)
        self.timer.timeout.connect(self._update_simulation)
        self.timer.start(16) # ~60 FPS

    def initializeGL(self):
        # --- GPU DIAGNOSTICS ---
        print("--- Initializing OpenGL Context ---")
        try:
            vendor = glGetString(GL_VENDOR).decode('utf-8')
            renderer = glGetString(GL_RENDERER).decode('utf-8')
            version = glGetString(GL_VERSION).decode('utf-8')
            print(f"OpenGL Vendor:   {vendor}")
            print(f"OpenGL Renderer: {renderer}")
            print(f"OpenGL Version:  {version}")

            # Check if we are running on the desired GPU
            if "nvidia" not in renderer.lower():
                print("\n\033[91mCRITICAL WARNING: OpenGL is NOT running on the NVIDIA GPU.\033[0m")
                print("\033[93mThe application was likely launched on the integrated Intel GPU by default.\033[0m")
                print("\033[92mTO FIX THIS, you MUST relaunch the script using the following command:\033[0m")
                print(f"\033[92m$ __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python {sys.argv[0]}\033[0m\n")
            else:
                print("\n\033[92mSUCCESS: OpenGL is running on the NVIDIA GPU.\033[0m\n")
        except Exception as e:
            print(f"Could not retrieve OpenGL device info: {e}")

        # Initialize particle data
        self.positions = np.zeros((self.particle_count, 4), dtype=np.float32)
        self.velocities = np.zeros((self.particle_count, 4), dtype=np.float32)
        
        # Initialize particles in a disk for a galaxy-like simulation
        radii = np.random.rand(self.particle_count) * 1.5
        angles = np.random.rand(self.particle_count) * 2.0 * np.pi
        self.positions[:, 0] = radii * np.cos(angles)
        self.positions[:, 1] = radii * np.sin(angles)
        self.positions[:, 3] = 1.0 # W-component
        
        # Give initial orbital velocity
        self.velocities[:, 0] = -self.positions[:, 1] * 0.1
        self.velocities[:, 1] =  self.positions[:, 0] * 0.1

        self.masses = np.ones(self.particle_count, dtype=np.float32) * 10.0 # Make particles heavier
        
        self.init_compute_shader()
        self.init_render_pipeline()
        
    def init_compute_shader(self):
        print("Initializing Compute Shader and SSBOs...")
        self.compute_program = compileProgram(compileShader(COMPUTE_SHADER_SOURCE, GL_COMPUTE_SHADER))
        
        self.position_ssbo = glGenBuffers(1)
        self.velocity_ssbo = glGenBuffers(1)
        self.mass_ssbo = glGenBuffers(1)
        
        upload_buffer(0, self.position_ssbo, self.positions)
        upload_buffer(1, self.velocity_ssbo, self.velocities)
        upload_buffer(2, self.mass_ssbo, self.masses)
        
        self.setup_compute_uniforms()
        print("Compute setup complete.")
        
    def init_render_pipeline(self):
        print("Initializing Render Pipeline...")
        self.render_program = compileProgram(
            compileShader(VERTEX_SHADER_SOURCE, GL_VERTEX_SHADER),
            compileShader(FRAGMENT_SHADER_SOURCE, GL_FRAGMENT_SHADER)
        )
        
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        
        # The VAO's vertex attribute should point to the position SSBO
        glBindBuffer(GL_ARRAY_BUFFER, self.position_ssbo)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(0))
        glEnableVertexAttribArray(0)
        
        glBindVertexArray(0)
        self.setup_projection_matrix()
        print("Render setup complete.")

    def setup_compute_uniforms(self):
        glUseProgram(self.compute_program)
        glUniform1f(glGetUniformLocation(self.compute_program, "dt"), self.dt)
        glUniform1f(glGetUniformLocation(self.compute_program, "G"), self.G)
        glUniform1i(glGetUniformLocation(self.compute_program, "particleCount"), self.particle_count)
        glUseProgram(0)

    def perform_compute_operation(self):
        self.makeCurrent()
        glUseProgram(self.compute_program)
        
        # Calculate the number of workgroups needed
        workgroup_size = 256
        num_workgroups = (self.particle_count + workgroup_size - 1) // workgroup_size
        glDispatchCompute(num_workgroups, 1, 1)
        
        # Barrier to ensure compute shader finishes writing before rendering reads
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)
        glUseProgram(0)

    def setup_projection_matrix(self):
        self.mvp_matrix.setToIdentity()
        self.mvp_matrix.perspective(45.0, self.width() / self.height(), 0.1, 100.0)
        self.mvp_matrix.lookAt(QVector3D(0, 0, 5), QVector3D(0, 0, 0), QVector3D(0, 1, 0))
        
    def paintGL(self):
        glClearColor(0.1, 0.1, 0.1, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        if not self.render_program or not self.vao:
            return # Don't draw if setup is incomplete
            
        glUseProgram(self.render_program)
        mvp_loc = glGetUniformLocation(self.render_program, "mvp")
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, self.mvp_matrix.data())
        
        glBindVertexArray(self.vao)
        glDrawArrays(GL_POINTS, 0, self.particle_count)
        glBindVertexArray(0)
        glUseProgram(0)

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        self.setup_projection_matrix() # Update projection matrix on resize

    def closeEvent(self, event):
        print("Cleaning up OpenGL resources...")
        self.makeCurrent()
        glDeleteProgram(self.compute_program)
        glDeleteProgram(self.render_program)
        glDeleteBuffers(3, [self.position_ssbo, self.velocity_ssbo, self.mass_ssbo])
        glDeleteVertexArrays(1, [self.vao])
        super().closeEvent(event)

    def _update_simulation(self):
        self.perform_compute_operation()
        self.update() # Triggers paintGL

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("N-Body Simulation with OpenGL Compute Shaders")
        self.setGeometry(100, 100, 800, 800)

        fmt = QSurfaceFormat()
        fmt.setVersion(4, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        QSurfaceFormat.setDefaultFormat(fmt)

        self.gl_widget = MyOpenGLWidget(self)
        self.setCentralWidget(self.gl_widget)

if __name__ == "__main__":
    # Run Like this:
    #   python -u -m pyBall.GUI.NBody_glsl
    #   __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u -m pyBall.GUI.NBody_glsl
    
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())