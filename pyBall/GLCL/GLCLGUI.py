import sys
import numpy as np

from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.arrays import vbo
import OpenGL.GL as GL

from .OGLsystem import GLobject, Mesh, InstancedData, compile_shader_program, create_sphere_mesh # Import necessary OpenGL utilities
from .OCLsystem import OCLSystem # Import OpenCL system

# These calls will fail at module load time - moved to initializeGL
# print("GL_VERSION =",   GL.glGetString(GL.GL_VERSION).decode())
# print("GLSL_VERSION =", GL.glGetString(GL.GL_SHADING_LANGUAGE_VERSION).decode())

class GLCLWidget(QOpenGLWidget):
    def __init__(self, parent=None, enable_opengl_debug=False):
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
        self.enable_opengl_debug = enable_opengl_debug

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

    def on_exception(self, e):
        print(f"Simulation error: {e}")
        import traceback
        traceback.print_exc()
        self.simulation_running = False
        self.sim_timer.stop()
        os._exit(1)  # Forcefully terminate to avoid PyQt swallowing exceptions
        
    def set_shader_sources(self, name, vertex_src, fragment_src):
        """Store shader sources for backward compatibility, but shaders are now compiled by GLCLBrowser."""
        self.shader_name = name
        self.vertex_shader_src = vertex_src
        self.fragment_shader_src = fragment_src
        print(f"GLCLWidget::set_shader_sources() Stored shader sources for '{name}' (for backward compatibility)")

    def set_shader_program(self, program):
        self.shader_program = program
        print(f"GLCLWidget::set_shader_program() Set shader program: {program}")

    def set_shader_programs(self, shader_programs):
        """Receive compiled shader programs from GLCLBrowser."""
        self.shader_programs = shader_programs
        print(f"GLCLWidget::set_shader_programs() Received {len(shader_programs)} shader programs: {list(shader_programs.keys())}")

    def set_buffer_data(self, buffer_data, render_pipeline_info):
        """Set buffer data and render pipeline for baking GLobjects."""
        self.buffer_data = buffer_data
        self.render_pipeline_info = render_pipeline_info
        # GLobjects will be baked in initializeGL

    def bake_render_objects(self):
        """Create GLobjects for each render pass using precompiled shaders from GLCLBrowser."""
        self.render_objects.clear()
        
        if not self.render_pipeline_info:
            print("GLCLWidget::bake_render_objects() No render pipeline defined")
            return
            
        print(f"GLCLWidget::bake_render_objects() Starting baking process using pipeline approach...")
        print(f"GLCLWidget::bake_render_objects() render_pipeline_info: {self.render_pipeline_info}")
        print(f"GLCLWidget::bake_render_objects() buffer_data keys: {list(self.buffer_data.keys())}")
        print(f"GLCLWidget::bake_render_objects() Number of render passes: {len(self.render_pipeline_info)}")
        
        # Check if we have shader programs
        if not hasattr(self, 'shader_programs'):
            print("GLCLWidget::bake_render_objects() No shader programs available")
            return
        
        # 2. For each render pass, get the compiled shader program handle from GLCLBrowser
        for i, render_pass in enumerate(self.render_pipeline_info):
            shader_name, element_count_name, vertex_buffer_name = render_pass[0], render_pass[1], render_pass[2]
            index_buffer_name = render_pass[3] if len(render_pass) > 3 else None
            
            print(f"GLCLWidget::bake_render_objects() Processing render pass {i}: shader='{shader_name}', vertex_buffer='{vertex_buffer_name}', index_buffer='{index_buffer_name}'")
            
            # Get precompiled shader program from GLCLBrowser
            shader_program = self.shader_programs.get(shader_name)
            if not shader_program:
                shader_program = self.ogl_system.get_shader_program(shader_name)
                if not shader_program:
                    #print(f"GLCLWidget::bake_render_objects() ERROR: Shader program '{shader_name}' not found")
                    raise Exception(f"GLCLWidget::bake_render_objects() ERROR: Shader program '{shader_name}' not found")
            
            # 3. Create GLobject with VAO/VBO for the buffer data
            buffer_info = self.buffer_data.get(vertex_buffer_name)
            if buffer_info is None or buffer_info['data'] is None:
                #print(f"GLCLWidget::bake_render_objects() ERROR: Buffer data for '{vertex_buffer_name}' not found")
                raise Exception(f"GLCLWidget::bake_render_objects() ERROR: Buffer data for '{vertex_buffer_name}' not found")
            
            gl_obj = GLobject(nelements=buffer_info['nelements'], mode=GL_POINTS)
            gl_obj.vao = GL.glGenVertexArrays(1)
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
            gl_obj.shader_name = shader_name  # Store shader name for debugging
            print(f"Associated shader program {shader_program} with GLobject")
            
            # Store the baked render object as (GLobject, shader_program) tuple
            # Note: these are handles (values) not names (keys)
            self.render_objects.append((gl_obj, shader_program))
            print(f"GLCLWidget::bake_render_objects() SUCCESS: Baked render object #{len(self.render_objects)} for render pass '{shader_name}' with {buffer_info['nelements']} elements")

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
        
        # Setup OpenGL debug output if requested
        if hasattr(self, 'enable_opengl_debug') and self.enable_opengl_debug:
            from .OGLsystem import setup_opengl_debug
            setup_opengl_debug(enable=True, synchronous=True)
        
        # Initial camera and projection setup
        self.update_matrices()

        # Now that OpenGL context is ready, compile shaders and bake render objects
        if hasattr(self, 'buffer_data') and hasattr(self, 'render_pipeline_info'):
            print("GLCLWidget::initializeGL() Compiling shaders and baking render objects after OpenGL context is ready")
            # Compile shaders now that OpenGL context is available
            if hasattr(self, 'browser') and self.browser:
                print("GLCLWidget::initializeGL() Calling browser.compile_shaders()")
                self.browser.compile_shaders()
                print("GLCLWidget::initializeGL() Shader compilation complete")
            else:
                print("GLCLWidget::initializeGL() No browser reference available")
            self.bake_render_objects()
        else:
            print("Warning: No buffer data or render pipeline provided to GLCLWidget.")

        # Initial camera and projection setup
        self.update_matrices()

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glClearColor(0.2, 0.2, 0.2, 1.0)
        
        # Debug: Print matrix details
        view_translation = self.view_matrix.column(3)
        #print(f"View matrix translation: ({view_translation.x()}, {view_translation.y()}, {view_translation.z()})")
        #print(f"Projection matrix near: {self.projection_matrix[3,2]}")
        
        # Check OpenGL errors before drawing
        err = glGetError()
        if err != GL_NO_ERROR:
            print(f"OpenGL error before draw: {err}")
        
        if hasattr(self, 'render_objects') and self.render_objects:
            # Debug print OpenGL buffer data before rendering
            if hasattr(self.parent(), 'bDebugGL') and self.parent().bDebugGL:
                print(f"=== OpenGL Rendering Debug ===")
            
            for i, (gl_obj, shader_program) in enumerate(self.render_objects):
                #print(f"RenderObject {i} has {gl_obj.nelements} elements in VBO {gl_obj.vbo}")
                
                glUseProgram(shader_program)
                # Set large point size and red color
                point_size_loc = glGetUniformLocation(shader_program, "point_size")
                color_loc = glGetUniformLocation(shader_program, "color")
                #print(f"Shader uniforms - point_size: {point_size_loc}, color: {color_loc}")
                
                glUniform1f(point_size_loc, 4.0)  # Very large points
                glUniform4f(color_loc, 1.0, 0.0, 0.0, 1.0)
                glDisable(GL_DEPTH_TEST)
                
                # Debug print buffer data being rendered
                if hasattr(self.parent(), 'bDebugGL') and self.parent().bDebugGL:
                    print(f"Rendering object with {gl_obj.nelements} vertices")
                    if gl_obj.vbo is not None:
                        # Read back buffer data for debugging
                        glBindBuffer(GL_ARRAY_BUFFER, gl_obj.vbo)
                        buffer_size = glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE)
                        if buffer_size > 0:
                            data = glGetBufferSubData(GL_ARRAY_BUFFER, 0, min(buffer_size, 64))
                            import numpy as np
                            float_data = np.frombuffer(data, dtype=np.float32)
                            print(f"  First 4 vertex values: {float_data[:4]}")
                
                glBindVertexArray(gl_obj.vao)
                glDrawArrays(GL_POINTS, 0, gl_obj.nelements)
                glBindVertexArray(0)
                
                # Check OpenGL errors after drawing
                err = glGetError()
                if err != GL_NO_ERROR:
                    print(f"OpenGL error after draw: {err}")
                
                glEnable(GL_DEPTH_TEST)
                glUseProgram(0)
                #print(f"Drew {gl_obj.nelements} elements with size 25px")
        else:
            print("WARNING: No render objects to draw")

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
            # Ensure shader program is active before setting uniforms
            current_program = glGetIntegerv(GL_CURRENT_PROGRAM)
            need_to_bind = current_program != self.shader_program
            if need_to_bind:
                print(f"GLCLWidget::update_matrices() Binding shader program {self.shader_program} before setting uniforms")
                glUseProgram(self.shader_program)

            if self.color_loc is None:
                self.get_default_uniform_locations(self.shader_program)

            #if self.color_loc != -1: glUniform4f(self.color_loc, 1.0, 1.0, 1.0, 1.0)
            if self.color_loc != -1: glUniform4fv(self.color_loc, 1, self.color)
            if self.model_loc != -1: glUniformMatrix4fv(self.model_loc, 1, GL_FALSE, self.model_matrix.data())
            if self.view_loc  != -1: glUniformMatrix4fv(self.view_loc, 1, GL_FALSE,  self.view_matrix.data())
            if self.proj_loc  != -1: glUniformMatrix4fv(self.proj_loc, 1, GL_FALSE,  self.projection_matrix.data())

            # Restore previous program if we changed it
            if need_to_bind and current_program != 0:
                glUseProgram(current_program)

    def set_systems(self, ogl_system, ocl_system, browser=None):
        if self.color_loc is None:
            self.get_default_uniform_locations(self.shader_program)

        #if self.color_loc != -1: glUniform4f(self.color_loc, 1.0, 1.0, 1.0, 1.0)
        if self.color_loc != -1: glUniform4fv(self.color_loc, 1, self.color)
        if self.model_loc != -1: glUniformMatrix4fv(self.model_loc, 1, GL_FALSE, self.model_matrix.data())
        if self.view_loc  != -1: glUniformMatrix4fv(self.view_loc, 1, GL_FALSE,  self.view_matrix.data())
        if self.proj_loc  != -1: glUniformMatrix4fv(self.proj_loc, 1, GL_FALSE,  self.projection_matrix.data())
        
        # Restore previous program if we changed it
        if need_to_bind and current_program != 0:
            glUseProgram(current_program)

    def set_systems(self, ogl_system: object, ocl_system: object, browser: object = None) -> None:
        """Set references to the OpenGL and OpenCL systems."""
        self.ogl_system = ogl_system
        self.ocl_system = ocl_system
        self.browser = browser
        print(f"GLCLWidget::set_systems() Set systems: ogl={ogl_system}, ocl={ocl_system}")
