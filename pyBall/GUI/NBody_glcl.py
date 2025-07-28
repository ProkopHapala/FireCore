import os
import sys
import numpy as np
import pyopencl as cl
from OpenGL.GL import *
from PyQt5.QtWidgets import QApplication, QMainWindow, QOpenGLWidget
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QSurfaceFormat
import pyopencl.tools

# --- How to Run This Script Correctly on a Hybrid Graphics Laptop (Linux) ---
# To force BOTH OpenGL and OpenCL onto the NVIDIA GPU for high-performance
# data sharing, you MUST launch the script with these environment variables:
#
# $ __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u your_script_name.py
#
# This is the required OS-level instruction.
# -----------------------------------------------------------------------------

# Platform-specific imports for GL/CL sharing
if sys.platform == "linux":
    from OpenGL.GLX import glXGetCurrentContext, glXGetCurrentDisplay
elif sys.platform == "win32":
    from OpenGL.WGL import wglGetCurrentContext, wglGetCurrentDC

# OpenCL Kernel, Shaders, etc. remain the same...
OPENCL_SOURCE = """
__kernel void nbody_sim(__global float4* positions, __global float4* velocities, const float dt, const int particle_count) {
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
        float dist_sq = dot(diff, diff) + 0.1f;
        float dist = sqrt(dist_sq);
        float3 direction = diff / dist;
        force += direction / dist_sq;
    }
    my_vel += force * G * dt;
    my_pos += my_vel * dt;
    positions[i].xyz = my_pos;
    velocities[i].xyz = my_vel;
}
"""
VERTEX_SHADER_SOURCE = """
#version 330 core
layout (location = 0) in vec4 position;
void main() { gl_Position = vec4(position.xyz, 1.0); gl_PointSize = 2.0; }
"""
FRAGMENT_SHADER_SOURCE = """
#version 330 core
out vec4 FragColor;
void main() { FragColor = vec4(1.0, 1.0, 1.0, 1.0); }
"""

def find_nvidia_opencl_device():
    """ Scans and automatically selects an NVIDIA OpenCL device. """
    print("--- Scanning for OpenCL Devices ---")
    platforms = cl.get_platforms()
    if not platforms: raise RuntimeError("FATAL: No OpenCL platforms found!")
    all_devices = [dev for p in platforms for dev in p.get_devices()]
    if not all_devices: raise RuntimeError("FATAL: No OpenCL devices found!")

    print("Available devices:")
    for i, dev in enumerate(all_devices):
        print(f"  [{i}] {dev.name} (Platform: {dev.platform.name})")

    nvidia_device = next((dev for dev in all_devices if "nvidia" in dev.vendor.lower()), None)
    
    if nvidia_device is None:
        raise RuntimeError("FATAL: Could not automatically find an NVIDIA OpenCL device.")
    
    print(f"\n---> Automatically selected NVIDIA device: {nvidia_device.name}\n")
    return nvidia_device


class MyOpenGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.particle_count = 2048
        self.dt = 0.001
        self.cl_context, self.cl_queue, self.nbody_kernel = None, None, None
        self.gl_vbo, self.vao, self.render_program = None, None, None
        self.cl_gl_buffer, self.cl_velocities, self.cl_positions = None, None, None
        self.use_gl_sharing = False

        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(16)

    def initializeGL(self):
        print("--- Initializing OpenGL Context ---")
        self.gl_vendor = glGetString(GL_VENDOR).decode('utf-8')
        self.gl_renderer = glGetString(GL_RENDERER).decode('utf-8')
        print(f"OpenGL Vendor: {self.gl_vendor}")
        print(f"OpenGL Renderer: {self.gl_renderer}")

        if "nvidia" not in self.gl_renderer.lower():
            print("\n\033[91mCRITICAL WARNING: OpenGL is NOT running on the NVIDIA GPU.\033[0m")
            print("\033[92mTO FIX THIS, you MUST relaunch the script using:\033[0m")
            print("\033[92m$ __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python your_script.py\033[0m\n")

        glClearColor(0.1, 0.1, 0.1, 1.0)
        self.positions = (np.random.rand(self.particle_count, 4) * 2 - 1).astype(np.float32)
        self.positions[:, 2] = 0.0; self.positions[:, 3] = 1.0
        self.velocities = (np.random.rand(self.particle_count, 4) * 0.1).astype(np.float32)
        
        try:
            self.setup_render_pipeline()
            self.setup_opencl()
        except Exception as e:
            print(f"\n\033[91mFATAL ERROR during setup: {e}\033[0m")
            self.parent().close()

    def setup_render_pipeline(self):
        print("--- Setting up OpenGL Render Pipeline ---")
        self.gl_vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glBufferData(GL_ARRAY_BUFFER, self.positions.nbytes, self.positions, GL_DYNAMIC_DRAW)
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)
        glBindVertexArray(0)
        self.render_program = self.compile_program(VERTEX_SHADER_SOURCE, FRAGMENT_SHADER_SOURCE)
        print("OpenGL pipeline setup complete.")

    def setup_opencl(self):
        print("\n--- Setting up OpenCL ---")
        self.cl_device = find_nvidia_opencl_device()
        self.cl_platform = self.cl_device.platform

        # *** THE FIX IS HERE ***
        try:
            print("Attempting to create OpenGL-OpenCL shared context...")
            properties = [(cl.context_properties.PLATFORM, self.cl_platform)]
            
            if sys.platform == "linux":
                print("Querying for GLX context and display handles...")
                properties.extend([
                    (cl.context_properties.GL_CONTEXT_KHR, glXGetCurrentContext()),
                    (cl.context_properties.GLX_DISPLAY_KHR, glXGetCurrentDisplay())
                ])
            elif sys.platform == "win32":
                print("Querying for WGL context and device context handles...")
                properties.extend([
                    (cl.context_properties.GL_CONTEXT_KHR, wglGetCurrentContext()),
                    (cl.context_properties.WGL_HDC_KHR, wglGetCurrentDC())
                ])
            else:
                raise RuntimeError(f"Unsupported platform for GL-CL sharing: {sys.platform}")

            self.cl_context = cl.Context(properties=properties, devices=[self.cl_device])
            self.cl_queue = cl.CommandQueue(self.cl_context)
            self.cl_gl_buffer = cl.GLBuffer(self.cl_context, cl.mem_flags.READ_WRITE, int(self.gl_vbo))
            self.use_gl_sharing = True
            print("\033[92mSUCCESS: OpenGL-OpenCL shared context created. Using high-performance path.\033[0m")

        except Exception as e:
            print(f"\033[91mERROR: Could not create shared GL-CL context.\033[0m")
            print(f"Details: {e}")
            print("\033[93mFalling back to slow, manual memory copy.\033[0m")
            self.use_gl_sharing = False
            self.cl_context = cl.Context(devices=[self.cl_device])
            self.cl_queue = cl.CommandQueue(self.cl_context)
            self.cl_positions = cl.Buffer(self.cl_context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.positions)

        self.cl_velocities = cl.Buffer(self.cl_context, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR, hostbuf=self.velocities)
        print("Compiling OpenCL kernel...")
        program = cl.Program(self.cl_context, OPENCL_SOURCE).build()
        self.nbody_kernel = program.nbody_sim
        print("OpenCL setup complete.")

    def update_simulation(self):
        if not self.nbody_kernel: return
        self.makeCurrent()
        pos_buffer = self.cl_gl_buffer if self.use_gl_sharing else self.cl_positions
        if self.use_gl_sharing:
            cl.enqueue_acquire_gl_objects(self.cl_queue, [self.cl_gl_buffer])
        self.nbody_kernel(self.cl_queue, (self.particle_count,), None, pos_buffer, self.cl_velocities, np.float32(self.dt), np.int32(self.particle_count))
        if self.use_gl_sharing:
            cl.enqueue_release_gl_objects(self.cl_queue, [self.cl_gl_buffer])
        else:
            cl.enqueue_copy(self.cl_queue, self.positions, pos_buffer)
            glBindBuffer(GL_ARRAY_BUFFER, self.gl_vbo)
            glBufferSubData(GL_ARRAY_BUFFER, 0, self.positions.nbytes, self.positions)
            glBindBuffer(GL_ARRAY_BUFFER, 0)
        self.cl_queue.finish()
        self.update()

    def paintGL(self):
        if not self.render_program: return
        glClear(GL_COLOR_BUFFER_BIT)
        glUseProgram(self.render_program)
        glBindVertexArray(self.vao)
        glEnable(GL_PROGRAM_POINT_SIZE)
        glDrawArrays(GL_POINTS, 0, self.particle_count)
        glBindVertexArray(0)

    def compile_program(self, vs_src, fs_src):
        def compile_shader(source, shader_type):
            shader = glCreateShader(shader_type)
            glShaderSource(shader, source)
            glCompileShader(shader)
            if not glGetShaderiv(shader, GL_COMPILE_STATUS):
                raise RuntimeError(glGetShaderInfoLog(shader).decode())
            return shader
        vs = compile_shader(vs_src, GL_VERTEX_SHADER)
        fs = compile_shader(fs_src, GL_FRAGMENT_SHADER)
        prog = glCreateProgram()
        glAttachShader(prog, vs); glAttachShader(prog, fs)
        glLinkProgram(prog)
        if not glGetProgramiv(prog, GL_LINK_STATUS):
            raise RuntimeError(glGetProgramInfoLog(prog).decode())
        return prog

    def closeEvent(self, event):
        print("Closing application and cleaning up resources...")
        self.timer.stop()
        super().closeEvent(event)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("N-Body Simulation - High Performance")
        self.setGeometry(100, 100, 800, 800)
        fmt = QSurfaceFormat()
        fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        QSurfaceFormat.setDefaultFormat(fmt)
        self.setCentralWidget(MyOpenGLWidget(self))

if __name__ == "__main__":
    # run like this:
    #  python -u -m pyBall.GUI.NBody_glcl
    #  __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia python -u -m pyBall.GUI.NBody_glcl

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())