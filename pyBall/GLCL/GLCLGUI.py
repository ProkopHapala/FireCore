import sys
import numpy as np

from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from OpenGL.GL import (
    glClear, glClearColor, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, glEnable, GL_DEPTH_TEST, glViewport,
    glUseProgram, glGetUniformLocation, glUniform4fv, glUniformMatrix4fv, GL_FALSE,
    glGenBuffers, glBindBuffer, glBufferData, GL_ARRAY_BUFFER, GL_DYNAMIC_DRAW,
    glGenVertexArrays, glBindVertexArray, glVertexAttribPointer, glEnableVertexAttribArray, GL_FLOAT,
    glBufferSubData, glDrawArrays, GL_POINTS, GL_PROGRAM_POINT_SIZE
)
import OpenGL.GL as GL

from .OGLsystem import GLobject, Mesh, InstancedData, compile_shader_program, create_sphere_mesh # Import necessary OpenGL utilities
from .OCLsystem import OCLSystem # Import OpenCL system

# These calls will fail at module load time - moved to initializeGL
# print("GL_VERSION =",   GL.glGetString(GL.GL_VERSION).decode())
# print("GLSL_VERSION =", GL.glGetString(GL.GL_SHADING_LANGUAGE_VERSION).decode())

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
        self.gl_vbo = 0
        self.vao = 0
        # Dynamic buffer registry
        self.gl_objects = {}  # Dictionary mapping buffer names to GLobject instances
        self.render_pipeline_info = []  # List of render pass configurations
        self.buffer_data = {}  # Buffer data for deferred GLobject creation
        self.render_objects = []  # List of pre-baked GLobjects for rendering
        self.current_element_count = 0

        self.color_loc = None
        self.model_loc = None
        self.view_loc  = None
        self.proj_loc  = None

    def set_shader_sources(self, name, vertex_src, fragment_src):
        self.shader_name = name
        self.vertex_shader_src = vertex_src
        self.fragment_shader_src = fragment_src

    def set_shader_program(self, program):
        self.shader_program = program

    def set_buffer_data(self, buffer_data, render_pipeline_info):
        """Set buffer data and render pipeline for baking GLobjects."""
        self.buffer_data = buffer_data
        self.render_pipeline_info = render_pipeline_info
        # GLobjects will be baked in initializeGL

    def bake_render_objects(self):
        """Bake render objects once OpenGL context is ready."""
        print("GLCLWidget::bake_render_objects() Starting baking process...")
        print(f"GLCLWidget::bake_render_objects() render_pipeline_info: {self.render_pipeline_info}")
        print(f"GLCLWidget::bake_render_objects() buffer_data keys: {list(self.buffer_data.keys())}")
        
        if not self.render_pipeline_info:
            print("GLCLWidget::bake_render_objects() No render pipeline defined")
            return
            
        print(f"GLCLWidget::bake_render_objects() Number of render passes: {len(self.render_pipeline_info)}")
        
        # Force OpenGL context to be current
        self.makeCurrent()
        
        # Check OpenGL version
        try:
            from OpenGL.GL import glGetString, GL_VERSION
            version = glGetString(GL_VERSION)
            print(f"OpenGL version: {version.decode('utf-8') if version else 'Unknown'}")
        except Exception as e:
            print(f"Could not get OpenGL version: {e}")
        
        for i, render_pass in enumerate(self.render_pipeline_info):
            shader_name, count_expr, vertex_buffer_name, uniforms = render_pass
            print(f"GLCLWidget::bake_render_objects() Processing render pass {i}: {render_pass}")
            print(f"GLCLWidget::bake_render_objects() shader_name: {shader_name}")
            print(f"GLCLWidget::bake_render_objects() vertex_buffer_name: {vertex_buffer_name}")
            
            # First ensure we have the shader source
            if not hasattr(self, 'vertex_shader_src') or not hasattr(self, 'fragment_shader_src'):
                print("GLCLWidget::bake_render_objects() ERROR: No shader sources available")
                continue
                
            # Ensure shader program is compiled and linked
            print(f"GLCLWidget::bake_render_objects() Loading shader '{shader_name}' into OGLSystem...")
            self.ogl_system.load_shader_program(shader_name, self.vertex_shader_src, self.fragment_shader_src)
            shader_program = self.ogl_system.get_shader_program(shader_name)
            
            if not shader_program:
                print(f"GLCLWidget::bake_render_objects() ERROR: Failed to create shader program '{shader_name}'")
                continue
                
            print(f"GLCLWidget::bake_render_objects() Successfully created shader program: {shader_program}")
            
            # Get vertex buffer data
            if vertex_buffer_name not in self.buffer_data:
                print(f"GLCLWidget::bake_render_objects() ERROR: Vertex buffer '{vertex_buffer_name}' not found")
                print(f"GLCLWidget::bake_render_objects() Available buffers: {list(self.buffer_data.keys())}")
                continue
                
            buffer_info = self.buffer_data[vertex_buffer_name]
            if buffer_info['data'] is None:
                print(f"GLCLWidget::bake_render_objects() ERROR: Vertex buffer '{vertex_buffer_name}' has no data")
                continue
            
            print(f"GLCLWidget::bake_render_objects() Creating GLobject for buffer '{vertex_buffer_name}' with {buffer_info['nelements']} elements")
            
            # Create simple GLobject with no initial setup
            gl_obj = GLobject(
                nelements=buffer_info['nelements'],
                mode=GL_POINTS
            )
            
            # Create VAO
            gl_obj.vao = GL.glGenVertexArrays(1)
            print(f"Created VAO: {gl_obj.vao}")
            
            # Create VBO
            gl_obj.vbo = GL.glGenBuffers(1)
            print(f"Created VBO: {gl_obj.vbo}")
            
            # Binding sequence is important:
            # 1. Bind VAO first
            GL.glBindVertexArray(gl_obj.vao)
            print("Bound VAO")
            
            # 2. Bind VBO and upload data
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, gl_obj.vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, buffer_info['data'].nbytes, buffer_info['data'], GL.GL_DYNAMIC_DRAW)
            print(f"Uploaded {buffer_info['data'].nbytes} bytes to VBO")
            
            # 3. Configure vertex attributes WHILE VAO is bound
            GL.glEnableVertexAttribArray(0)  # Position attribute at layout location 0
            GL.glVertexAttribPointer(0, buffer_info['components'], GL.GL_FLOAT, GL.GL_FALSE, 0, None)
            print(f"Configured vertex attribute: components={buffer_info['components']}")
            
            # 4. Unbind VAO and VBO when done
            GL.glBindVertexArray(0)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
            
            # Store shader program for drawing
            gl_obj.shader_program = shader_program
            print(f"Associated shader program {shader_program} with GLobject")
            
            # Test the shader by temporarily using it
            GL.glUseProgram(shader_program)
            color_loc = GL.glGetUniformLocation(shader_program, "color")
            print(f"Color uniform location: {color_loc}")
            if color_loc >= 0:
                # Set a default color if the uniform exists
                GL.glUniform4f(color_loc, 1.0, 0.0, 0.0, 1.0)  # Red color
            GL.glUseProgram(0)
            
            self.render_objects.append(gl_obj)
            print(f"GLCLWidget::bake_render_objects() SUCCESS: Baked GLobject #{len(self.render_objects)} for render pass '{shader_name}' with {buffer_info['nelements']} elements")

    def update_buffer_data(self, buffer_name, new_data):
        """Update buffer data for a specific buffer."""
        if buffer_name in self.gl_objects and new_data is not None:
            gl_obj = self.gl_objects[buffer_name]
            gl_obj.upload_vbo(new_data.astype(np.float32))
            gl_obj.nelements = len(new_data)
            if new_data.ndim > 1:
                gl_obj.components = new_data.shape[1]
            else:
                gl_obj.components = 1

    def initializeGL(self):
        glClearColor(0.1, 0.1, 0.1, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_PROGRAM_POINT_SIZE)
        
        # Initial camera and projection setup
        self.update_matrices()

        # Load shaders after OpenGL context is ready
        if hasattr(self, 'vertex_shader_src') and hasattr(self, 'fragment_shader_src'):
            print(f"Loading shader: {self.shader_name}")
            self.ogl_system.load_shader_program(self.shader_name, self.vertex_shader_src, self.fragment_shader_src)
            self.shader_program = self.ogl_system.get_shader_program(self.shader_name)
            # Ensure shader program is used before setting uniforms and getting locations
            if self.shader_program:
                glUseProgram(self.shader_program)
                self.get_default_uniform_locations(self.shader_program)
        else:
            print("Warning: Shader sources not provided to GLCLWidget.")

        # Initial camera and projection setup
        self.update_matrices()

    def paintGL(self):
        """Render the scene using pre-baked GLobjects."""
        if not self.render_objects:
            print("GLCLWidget::paintGL() No pre-baked render objects available")
            return
        print(f"GLCLWidget::paintGL() Drawing {len(self.render_objects)} pre-baked GLobjects")
        
        # Simple iteration over pre-baked GLobjects
        for gl_obj in self.render_objects:
            if gl_obj.nelements > 0:
                # Use pre-baked configuration
                glUseProgram(gl_obj.shader_program)
                glBindVertexArray(gl_obj.vao)
                glDrawArrays(GL_POINTS, 0, gl_obj.nelements)
                glBindVertexArray(0)
                glUseProgram(0)

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        self.update_matrices()

    def get_default_uniform_locations(self, shader_program):
        print("GLCLWidget::get_default_uniform_locations() Getting uniform locations")
        self.color_loc = glGetUniformLocation(shader_program, "color")
        print(f"GLCLWidget::get_default_uniform_locations() color_loc: {self.color_loc}")
        self.model_loc = glGetUniformLocation(shader_program, "model")
        self.view_loc  = glGetUniformLocation(shader_program, "view")
        self.proj_loc  = glGetUniformLocation(shader_program, "projection")

    def update_matrices(self):
        # Model matrix (identity for now, or set by simulation)
        self.model_matrix.setToIdentity()

        # View matrix
        self.view_matrix.setToIdentity()
        self.view_matrix.lookAt(self.camera_pos, QVector3D(0, 0, 0), QVector3D(0, 1, 0))

        # Projection matrix
        self.projection_matrix.setToIdentity()
        aspect_ratio = self.width() / self.height()
        self.projection_matrix.perspective(45.0, aspect_ratio, 0.1, 100.0)

        self.color = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float32)

        if self.shader_program:

            if self.color_loc is None:
                self.get_default_uniform_locations(self.shader_program)

            #if self.color_loc != -1: glUniform4f(self.color_loc, 1.0, 1.0, 1.0, 1.0)
            if self.color_loc != -1: glUniform4fv(self.color_loc, 1, self.color)
            if self.model_loc != -1: glUniformMatrix4fv(self.model_loc, 1, GL_FALSE, self.model_matrix.data())
            if self.view_loc  != -1: glUniformMatrix4fv(self.view_loc, 1, GL_FALSE,  self.view_matrix.data())
            if self.proj_loc  != -1: glUniformMatrix4fv(self.proj_loc, 1, GL_FALSE,  self.projection_matrix.data())

    def set_systems(self, ogl_system, ocl_system):
        self.ogl_system = ogl_system
        self.ocl_system = ocl_system
