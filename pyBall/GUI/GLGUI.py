import os
import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse
from PIL import Image

from PyQt5.QtWidgets import ( QOpenGLWidget)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from PyQt5.QtWidgets import (QMainWindow)
#from PyQt5.QtOpenGLWidgets import QOpenGLWidget
from PyQt5.QtGui import QSurfaceFormat
from PyQt5.QtWidgets import QApplication

from OpenGL.GL import (
    glUseProgram, glClear, glClearColor, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT, 
    glEnable, glDepthFunc, GL_LESS, GL_DEPTH_TEST, glViewport, glCreateProgram, 
    glShaderSource, glCompileShader, glAttachShader, glLinkProgram, glGetProgramiv, 
    GL_LINK_STATUS, glGetShaderiv, GL_COMPILE_STATUS, glGetProgramInfoLog, glGetShaderInfoLog, 
    glDeleteShader, glGetUniformLocation, glUniformMatrix4fv, GL_TRUE, glUniform3fv, GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, 
    glUniform1f, glUniform4f, glUniform1i, glGenBuffers, glBindBuffer, glBufferData, 
    GL_ARRAY_BUFFER, GL_STATIC_DRAW, GL_DYNAMIC_DRAW, glVertexAttribPointer, glEnableVertexAttribArray, 
    glDrawArrays, GL_TRIANGLES, glDeleteBuffers, glGenVertexArrays, glBindVertexArray, glVertexAttribDivisor, 
    glDeleteVertexArrays, GL_FLOAT, GL_FALSE, glDrawArraysInstanced, glBlendFunc, 
    glBlendEquation, GL_BLEND, GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, 
    glBlendEquationSeparate, glBlendFuncSeparate, GL_MIN, GL_MAX, GL_FUNC_SUBTRACT, 
    GL_FUNC_REVERSE_SUBTRACT, GL_ONE, glActiveTexture, GL_TEXTURE0, glBindTexture, GL_TEXTURE_2D,
    glGenTextures, glTexParameteri, GL_TEXTURE_WRAP_S, GL_TEXTURE_WRAP_T, GL_REPEAT, 
    GL_TEXTURE_MIN_FILTER, GL_TEXTURE_MAG_FILTER, GL_LINEAR, GL_RGBA, GL_UNSIGNED_BYTE, 
    glPixelStorei, GL_UNPACK_ALIGNMENT, GL_CULL_FACE, glTexImage2D, glDisable
)
from OpenGL.GL.shaders import compileProgram, compileShader
import OpenGL.GL as GL
import ctypes

# It's good practice to keep shaders in separate files or as multi-line strings
# For simplicity here, they are embedded as strings.

alpha_blend_modes={
    "standard"    :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    "standard2"    :(GL_FUNC_ADD, GL_ONE, GL_ONE_MINUS_SRC_ALPHA),
    #"standard2"   :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE),
    #"standard3"   :(GL_FUNC_ADD, GL_ONE,       GL_ONE),
    #"minimum1"    :(GL_MIN,      GL_ONE,       GL_ONE),
    #"minimum2"    :(GL_MIN,      GL_SRC_ALPHA, GL_ONE),
    #"minimum3"    :(GL_MIN,      GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    #"minimum4"    :(GL_MIN,GL_ONE, GL_ONE,   GL_FUNC_ADD,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA ),
    "additive"   :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE),
    "subtractive" :(GL_FUNC_REVERSE_SUBTRACT, GL_SRC_ALPHA, GL_ONE),
    #"subtractive2" :(GL_FUNC_REVERSE_SUBTRACT, GL_ONE, GL_ONE),
    #"subtractive2" :(GL_FUNC_REVERSE_SUBTRACT, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    #"subtractive":(GL_FUNC_SUBTRACT, GL_SRC_ALPHA, GL_ONE),
    #"minimum"    :(GL_FUNC_MIN, GL_ONE, GL_ONE),
    #"maximum"    :(GL_MAX, GL_ONE, GL_ONE)
}

def rotation_to_gl_matrix( rotation_obj):
    mat3x3 = rotation_obj.as_matrix()
    mat4x4 = np.eye(4, dtype=np.float32)
    mat4x4[:3, :3] = mat3x3
    #return mat4x4.flatten(order='F').copy()
    return mat4x4.copy()

def set_ogl_blend_mode(mode, depth_test=True):
    #print("set_ogl_blend_mode", mode)
    glEnable(GL_BLEND)
    if len(mode) == 3:
        glBlendEquation(mode[0])
        glBlendFunc(mode[1], mode[2])
    else:
        glBlendEquationSeparate(mode[0], mode[3])
        glBlendFuncSeparate    (mode[1], mode[2], mode[4], mode[5])
    #glDepthMask(GL_FALSE)
    #glDepthTest(GL_FALSE)
    if depth_test:
        glDisable(GL_DEPTH_TEST)

def disable_blend( depth_test=True):
    glDisable(GL_BLEND)
    #glDepthMask(GL_TRUE)
    #glDepthTest(GL_TRUE)
    if depth_test:
        glEnable(GL_DEPTH_TEST)

def setup_vbo(vertices, array_indx, components=3, usage=GL_STATIC_DRAW):
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, usage)
    glVertexAttribPointer(array_indx, components, GL_FLOAT, GL_FALSE, 0, None) # Location 0 for positions
    glEnableVertexAttribArray(array_indx)
    return vbo

def setup_ebo(indices):
    ebo = glGenBuffers(1)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL_STATIC_DRAW)
    return ebo

class GLobject:

    def alloc_vao_vbo_ebo(self, components):
        self.vao = GL.glGenVertexArrays(1);  GL.glBindVertexArray(self.vao)
        self.vbo = GL.glGenBuffers(1);       GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vbo)
        self.ebo = GL.glGenBuffers(1);       GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ebo)
        if len(components) == 1:
            GL.glEnableVertexAttribArray(0)
            GL.glVertexAttribPointer(0, components[0], GL.GL_FLOAT, GL.GL_FALSE, 0, None)
        else:
            total_stride = sum(components) * ctypes.sizeof(GL.GLfloat) # Calculate total stride in bytes
            offset = 0
            for i, n_comp in enumerate(components):
                GL.glEnableVertexAttribArray(i)
                GL.glVertexAttribPointer(i, n_comp, GL.GL_FLOAT, GL.GL_FALSE, total_stride, ctypes.c_void_p(offset))
                offset += n_comp * ctypes.sizeof(GL.GLfloat) # Increment offset by size of current attribute
        GL.glBindVertexArray(0)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)
        return self.vao, self.vbo, self.ebo

    def __init__(self, vao=None, vbo=None, ebo=None, dirty=True, nelements=0, components=None, mode=GL_TRIANGLES ):
        if vao is not None: self.vao = vao
        if vbo is not None: self.vbo = vbo
        if ebo is not None: self.ebo = ebo
        self.dirty = True
        if components is not None: self.alloc_vao_vbo_ebo(components)
        self.nelements = nelements
        self.mode = mode

    def draw(self, n=-1 ):
        GL.glBindVertexArray(self.vao)
        GL.glEnableVertexAttribArray(0) # Ensure attribute 0 is enabled for drawing
        if n == -1: n= self.nelements
        GL.glDrawElements(self.mode, n, GL.GL_UNSIGNED_INT, None)
        GL.glBindVertexArray(0)

    def upload_ebo(self, indices):
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ebo)
        #print("indices", indices.shape, indices.nbytes)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)

    def upload_vbo(self, vertices):
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vbo)
        #print("vertices", vertices.shape, vertices.nbytes)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

    def upload_vbo_ebo(self, indices, vertices ):
        self.upload_ebo(indices)
        self.upload_vbo(vertices)
        self.nelements = indices.size
        #print("self.nelements", self.nelements)
class Mesh:
    def __init__(self, vertices, normals=None, indices=None):
        self.vertices = vertices
        self.normals = normals
        self.indices = indices

        self.vao = None
        self.vbo_vertices = None
        self.vbo_normals = None
        self.ebo = None

        self.has_indices = indices is not None
        self.vertex_count = len(vertices) // 3
        self.index_count = len(indices) if self.has_indices else 0

    def setup_buffers(self):
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vbo_vertices = setup_vbo(self.vertices, 0)
        if self.normals is not None:
            self.vbo_normals = setup_vbo(self.normals, 1)
        if self.has_indices:
            self.ebo = setup_ebo(self.indices)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0) # Unbind VAO

    def draw(self, mode=GL_TRIANGLES):
        glBindVertexArray(self.vao)
        if self.has_indices:
            glDrawElements(mode, self.index_count, GL_UNSIGNED_INT, None)
        else:
            glDrawArrays(mode, 0, self.vertex_count)
        glBindVertexArray(0)

    def cleanup(self):
        if self.vbo_vertices: glDeleteBuffers(1, [self.vbo_vertices])
        if self.vbo_normals:  glDeleteBuffers(1, [self.vbo_normals])
        if self.has_indices and self.ebo: glDeleteBuffers(1, [self.ebo])
        if self.vao:          glDeleteVertexArrays(1, [self.vao])

    def bind(self):
        glBindVertexArray(self.vao)

    def unbind(self):
        glBindVertexArray(0)

    def draw_instanced(self, num_instances, mode=GL_TRIANGLES):
        # Assumes the correct VAO (e.g., from InstancedData) is already bound by the caller.
        if self.has_indices:
            glDrawElementsInstanced(mode, self.index_count, GL_UNSIGNED_INT, None, num_instances)
        else:
            glDrawArraysInstanced(mode, 0, self.vertex_count, num_instances)


class InstanceData:
    def __init__(self, base_attrib_location):
        self.vbos = {}
        self.base_attrib_location = base_attrib_location # Starting location for instance attributes
        self.associated_mesh = None
        self.num_instances = 0 # Now stored here
        self.vao = None # Each InstancedData will have its own VAO
        self._instance_attrib_configs = [] # To store VBO configs

def make_vbo( num_components, base_attrib_location, attrib_idx_offset, usage=GL_DYNAMIC_DRAW):
    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, 0, None, usage) # Initial empty allocation
    attrib_loc = base_attrib_location + attrib_idx_offset
    glVertexAttribPointer(attrib_loc, num_components, GL_FLOAT, GL_FALSE, 0, None)
    glEnableVertexAttribArray(attrib_loc)
    glVertexAttribDivisor(attrib_loc, 1) # Per-instance
    return vbo

class InstancedData:
    def __init__(self, base_attrib_location):
        self.vbos = []
        self.vbos_inds = {}
        self.base_attrib_location = base_attrib_location # Starting location for instance attributes
        self.associated_mesh = None
        self.num_instances = 0 # Now stored here
        self.vao = None # Each InstancedData will have its own VAO
        self._instance_attrib_configs = [] # To store VBO configs

    def cleanup(self):
        for vbo in self.vbos.values():
            glDeleteBuffers(1, [vbo])
        self.vbos.clear()
        if self.vao:
            glDeleteVertexArrays(1, [self.vao])
            self.vao = None

    def associate_mesh(self, mesh_object):
        self.associated_mesh = mesh_object
        if self.associated_mesh is None:
            if self.vao: # Cleanup old VAO if mesh is disassociated
                glDeleteVertexArrays(1, [self.vao])
                self.vao = None
            return

        if self.vao is None:
            self.vao = glGenVertexArrays(1)

        glBindVertexArray(self.vao)
        # Bind mesh's vertex buffer and set attribute pointer for location 0 (aPos)
        glBindBuffer(GL_ARRAY_BUFFER, self.associated_mesh.vbo_vertices)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)

        # Bind mesh's normal buffer and set attribute pointer for location 1 (aNormal)
        if self.associated_mesh.vbo_normals is not None:
            glBindBuffer(GL_ARRAY_BUFFER, self.associated_mesh.vbo_normals)
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, None)
            glEnableVertexAttribArray(1)
        
        # If mesh has EBO, bind it to this VAO as well
        if self.associated_mesh.ebo is not None:
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.associated_mesh.ebo)

        glBindVertexArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0) # Unbind GL_ARRAY_BUFFER
        if self.associated_mesh.ebo is not None:
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0) # Unbind GL_ELEMENT_ARRAY_BUFFER

    def setup_instance_vbos(self, attribs):
        if self.vao is None:
            print("Error: InstancedData VAO not initialized. Call associate_mesh first.")
            return
        glBindVertexArray(self.vao)
        self._instance_attrib_configs = [] # Clear previous configs
        for i, (name, attrib_idx_offset, num_components) in enumerate(attribs):
            vbo = make_vbo(num_components, self.base_attrib_location, attrib_idx_offset)
            self.vbos.append(vbo)
            self.vbos_inds[name]=i  # this should be faster than using name as key
        glBindVertexArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        return self.get_vbo_inds( [name for name, _, _ in attribs] )

    def get_vbo_inds(self, names):
        return [self.vbos_inds[name] for name in names]

    def update(self, buffs, n=-1):
        if n<0: n=len(buffs["positions"])
        self.num_instances = n
        for name, buff_data in buffs.items():
            i = self.vbos_inds[name]
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[i])
            glBufferData(GL_ARRAY_BUFFER, buff_data.nbytes, buff_data, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def update_list(self, buffs, n=-1):  # this should be faster than using name as key
        if n<0: n=len(buffs[0])
        self.num_instances = n
        for i, buff_data in enumerate(buffs):
            glBindBuffer(GL_ARRAY_BUFFER, self.vbos[i])
            glBufferData(GL_ARRAY_BUFFER, buff_data.nbytes, buff_data, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

    def draw(self, mode=GL_TRIANGLES):
        if self.vao is None or self.associated_mesh is None or self.num_instances == 0:
            return

        glBindVertexArray(self.vao)
        # The VAO is already configured with base mesh EBO if it exists (done in associate_mesh)
        self.associated_mesh.draw_instanced(self.num_instances, mode)
        glBindVertexArray(0)

# ================= Free Utility Functions =================

def create_sphere_mesh(radius=1.0, segments=16, rings=16):
    """Generates vertices and normals for a sphere."""
    vertices = []
    normals = []
    indices = []

    for ring_num in range(rings + 1):
        theta = ring_num * np.pi / rings
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        for seg_num in range(segments + 1):
            phi = seg_num * 2 * np.pi / segments
            sin_phi = np.sin(phi)
            cos_phi = np.cos(phi)

            x = radius * cos_phi * sin_theta
            y = radius * sin_phi * sin_theta
            z = radius * cos_theta
            
            nx, ny, nz = x/radius, y/radius, z/radius # Normals are just normalized positions for a sphere centered at origin

            vertices.extend([x, y, z])
            normals.extend([nx,ny,nz])

    for ring_num in range(rings):
        for seg_num in range(segments):
            first  = (ring_num * (segments + 1)) + seg_num
            second = first + segments + 1
            
            indices.extend([first, second, first + 1])
            indices.extend([second, second + 1, first + 1])
            
    return np.array(vertices, dtype=np.float32), np.array(normals, dtype=np.float32), np.array(indices, dtype=np.uint32)

def append_normalized( v_list, radius, vertices_list, normals_list ):
    for v in v_list:
        v_norm = v / np.linalg.norm(v)
        vertices_list.extend(v_norm * radius)
        normals_list .extend(v_norm)

def octahedron_sphere_face(c1, c2, c3, vertices_list, normals_list, nsub=2,  radius=1.0):
    """
    Generates vertices and normals for one face of an octahedron, subdivided,
    for non-indexed drawing. Appends directly to vertices_list and normals_list.
    """
    # Interpolation vectors along the edges of the base triangle face
    da = (c2 - c1) / nsub
    db = (c3 - c1) / nsub
    for i in range(nsub):      # Iterate along the c1-c2 edge direction
        for j in range(nsub - i):  # Iterate along the c1-c3 edge direction
            # Define the 3 vertices of the "bottom-left" triangle in the current cell
            p00 = c1 + da*i     + db*j
            p10 = c1 + da*(i+1) + db*j
            p01 = c1 + da*i     + db*(j+1)
            append_normalized([p00, p10, p01], radius, vertices_list, normals_list)
            # Define the 3 vertices of the "top-right" triangle in the current cell (if it exists)
            if j < nsub - 1 - i: # Check if there's a second triangle in this parallelogram
                p11 = c1 + da*(i+1) + db*(j+1)
                append_normalized([p10, p11, p01], radius, vertices_list, normals_list)

def octahedron_sphere_mesh(radius=1.0, nsub=2):
    """
    Generates vertices and normals for a sphere by subdividing an octahedron,
    vectorized with NumPy for non-indexed GL_TRIANGLES drawing.
    'subdivisions' is recursion depth: 0 for base octahedron (8 tris), 1 for 32 tris, etc.
    Returns flat vertices, flat normals, and None for indices.
    """
    pv_np = np.array([
        [ 1.0,  0.0,  0.0], [-1.0,  0.0,  0.0], [ 0.0,  1.0,  0.0], 
        [ 0.0, -1.0,  0.0], [ 0.0,  0.0,  1.0], [ 0.0,  0.0, -1.0]
    ], dtype=np.float32)
    octa_face_indices = np.array([ # CCW from outside
        [4, 1, 3], [4, 3, 0], [4, 0, 2], [4, 2, 1], # Top pyramid (pz, mx, my), ...
        [5, 3, 1], [5, 0, 3], [5, 2, 0], [5, 1, 2]   # Bottom pyramid (mz, my, mx), ...
    ], dtype=np.int32)
    vertices_list = []
    normals_list  = []
    for face_def_indices in octa_face_indices:
        c1, c2, c3 = pv_np[face_def_indices,:]
        octahedron_sphere_face(c1, c2, c3, vertices_list, normals_list, nsub=nsub, radius=radius)
    final_vertices = np.array(vertices_list, dtype=np.float32)
    final_normals  = np.array(normals_list,  dtype=np.float32)
    # print("final_vertices.shape", final_vertices.shape ) #, final_vertices ) # Debug: too verbose
    # print("final_normals.shape",  final_normals.shape ) #, final_normals )  # Debug: too verbose
    return final_vertices, final_normals

def make_labels(points, labels, font_atlas_data, text_object=None ):
    if text_object is None:
        text_object = GLobject( components=[3,2,2] )
    #if not self.trj or self.font_atlas_data is None:
    #    self.num_label_verts = 0
    #    return
    #atom_positions, _, _, atom_enames = self.frames_data[self.current_frame_index].atoms
    #print("atom_enames", atom_enames)
    tile_w = font_atlas_data['tile_w']  # width of one character tile in atlas
    tile_h = font_atlas_data['tile_h']  # height of one character tile in atlas
    tex_w  = font_atlas_data['tex_w']   # full atlas texture width
    all_vertex_data = []
    szx  = 1.0  # width of one character quad in world units
    szy  = 2.0*0.5  # height of one character quad in world units
    yoff = 0.0
    margin_lx = 0.1
    margin_rx = 0.7
    margin_y  = 0.0
    space     = 0.0
    conde_min = 32
    conde_max = 126
    for pos_3d, label in zip(points, labels):
        #ename += "_Hey"  # debug longer text
        symbol = label.decode('utf-8') if isinstance(label, bytes) else str(label) 
        #xoff     = -len(symbol)* szx / 2.0 # centered text
        xoff     = 0 # left aligned text
        for i_char, char in enumerate(symbol):
            code = ord(char)
            if code < conde_min or code > conde_max: continue
            # 1D atlas UVs  - using normalized coordinates
            u0 = (code - conde_min + margin_lx ) * tile_w / tex_w
            u1 = (code - conde_min + margin_rx ) * tile_w / tex_w
            #u1 = u0               + tile_w*(1-margin_x) / tex_w
            v1 = 0.0 + margin_y
            v0 = 1.0 - margin_y
            # Same 3D position for every character in this label
            base_3d = pos_3d + np.array([0.0, 0.0, 0.0])
            # Local screen-space offset for this character
            x = xoff + i_char * (szx*(1.0+space))
            local_offsets = np.array([
                [ x       , yoff-szy ],
                [ x + szx , yoff-szy ],
                [ x + szx , yoff+szy ],
                [ x       , yoff+szy ]
            ], dtype=np.float32)
            uvs = np.array([[u0, v0], [u1, v0], [u1, v1], [u0, v1]], dtype=np.float32)
            for j in range(4):
                vertex_data = np.concatenate((base_3d, local_offsets[j], uvs[j]))
                all_vertex_data.append(vertex_data)
    num_label_verts = len(all_vertex_data)
    #print(f"DEBUG: update_atom_labels_data created {self.num_label_verts} vertices.")
    if num_label_verts > 0:
        # Generate and upload index data for the EBO
        num_quads = num_label_verts // 4
        indices = np.zeros(num_quads * 6, dtype=np.uint32)
        for i in range(num_quads):
            base = i * 4
            indices[i*6:i*6+6] = [base, base + 1, base + 2, base, base + 2, base + 3]
        vertex_data_np = np.array(all_vertex_data, dtype=np.float32)
        #print("vertex_data_np #legend:  pos(x,y,z)   local_offset(x,y)   uv(x,y)\n", vertex_data_np)
        #print("indices: ", indices)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, text_object.vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, vertex_data_np.nbytes, vertex_data_np, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, text_object.ebo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        text_object.nelements = num_quads*6
    text_object.dirty = False
    return text_object

class BaseGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.shader_folder=os.path.dirname(os.path.abspath(__file__)) + "/shaders"

        self.text_shader     = None
        self.font_atlas_tex  = None
        self.font_atlas_data = None

        self.zoom_factor    = 10.0  # Initial zoom
        self.orientation    = R.identity()
        self.last_mouse_pos = None
        self.shader_program = None
        self.default_sphere_mesh = None # For common sphere rendering

        # Camera and light
        self.view_matrix       = QMatrix4x4()
        self.projection_matrix = QMatrix4x4()
        self.model_matrix      = QMatrix4x4()
        self.light_pos         = QVector3D(15.0, 15.0, 30.0)
        self.light_color       = QVector3D( 1.0,  1.0,  1.0)
        self.camera_pos        = QVector3D( 0  ,  0  , self.zoom_factor) # Will be updated by zoom
        self.current_shader_program_id = None # To be set by derived class if it manages shaders

        self.all_shader_programs = []


    # Qt will call paintGL(); delegate to our internal method
    def paintGL(self):
        from PyQt5.QtGui import QOpenGLContext
        ctx=QOpenGLContext.currentContext()
        print(f"BaseGLWidget.paintGL called, currentContext valid: {ctx is not None}")
        # Note: paintGL is called frequently; keep prints minimal
        self.paintGL_base()

    def initializeGL_base(self, vertex_shader_src, fragment_shader_src, bPrint=False):
        # Modern OpenGL context should be requested via QSurfaceFormat in main
        if bPrint:
            print(f"OpenGL Version: {glGetString(GL_VERSION).decode()}")
            print(f"GLSL Version: {glGetString(GL_SHADING_LANGUAGE_VERSION).decode()}")

        glClearColor(1.0, 1.0, 1.0, 1.0) # White background
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glEnable(GL_CULL_FACE) # Cull back faces

        self.simple_shader = self.compile_shader_program(
            """
            #version 330 core
            layout(location = 0) in vec3 position;
            uniform mat4 projection;
            uniform mat4 view;
            uniform mat4 model;
            
            void main() {
                gl_Position = projection * view * model * vec4(position, 1.0);
            }
            """,
            """
            #version 330 core
            out vec4 FragColor;
            uniform vec4 color;
            
            void main() {
                FragColor = color;
            }
            """
        )
        self.all_shader_programs.append(self.simple_shader)
        print(f"Added simple_shader {self.simple_shader} to all_shader_programs")
        
        # # --- Text Rendering Setup ---
        vert_text = open(self.shader_folder + "/text_billboard.glslv").read()
        frag_text = open(self.shader_folder + "/text_billboard.glslf").read()
        self.text_shader = self.compile_shader_program(vert_text, frag_text)
        self.all_shader_programs.append(self.text_shader)
        print(f"Added text_shader {self.text_shader} to all_shader_programs")


        # Compile shaders if sources are provided
        # if vertex_shader_src and fragment_shader_src:
        #     self.shader_program = self.compile_shader_program(vertex_shader_src, fragment_shader_src)
        # elif not hasattr(self, 'shader_program') or not self.shader_program:
        #      print("BaseGLWidget: No shader sources provided and no shader_program pre-set.")
             # self.shader_program should be set by derived class if sources are None

        # Initialize a default sphere mesh
        #sphere_v, sphere_n, sphere_idx = create_sphere_mesh(radius=1.0)
        #self.default_sphere_mesh = Mesh(vertices=sphere_v, normals=sphere_n, indices=sphere_idx)

        sphere_v, sphere_n = octahedron_sphere_mesh(radius=1.3, nsub=2)
        self.default_sphere_mesh = Mesh(vertices=sphere_v, normals=sphere_n)
        self.default_sphere_mesh.setup_buffers()

    def compile_shader_program(self, vertex_shader_src, fragment_shader_src):
        return compileProgram(
            compileShader(vertex_shader_src,   GL_VERTEX_SHADER),
            compileShader(fragment_shader_src, GL_FRAGMENT_SHADER)
        )

    def cleanupGL_base(self):
        if self.shader_program: glDeleteProgram(self.shader_program)
        self.shader_program = None
        if self.default_sphere_mesh: self.default_sphere_mesh.cleanup()
        if self.text_shader:    glDeleteProgram(self.text_shader)
        if self.font_atlas_tex: glDeleteTextures(1, [self.font_atlas_tex])

    def set_default_uniforms(self):
        if self.current_shader_program_id is None: print("BaseGLWidget.set_default_uniforms: current_shader_program_id is None") 
        glUniformMatrix4fv(glGetUniformLocation(self.current_shader_program_id, "projection" ), 1, GL_FALSE, self.projection_matrix.data()) # projection_matrix is QMatrix4x4
        glUniformMatrix4fv(glGetUniformLocation(self.current_shader_program_id, "view"       ), 1, GL_FALSE, self.view_matrix.data())       # view_matrix is QMatrix4x4
        # Ensure model_matrix is a numpy array before passing to OpenGL
        if isinstance(self.model_matrix, np.ndarray):
            glUniformMatrix4fv(glGetUniformLocation(self.current_shader_program_id, "model"), 1, GL_FALSE, self.model_matrix)
        else:
            # Handle QMatrix4x4 or other types by converting to numpy
            model_array = np.array(self.model_matrix.data(), dtype=np.float32)
            glUniformMatrix4fv(glGetUniformLocation(self.current_shader_program_id, "model"), 1, GL_FALSE, model_array)

    def paintGL_base(self):
        #print("BaseGLWidget.paintGL_base()")
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        # Camera / View transformation
        self.camera_pos = QVector3D(0, 0, self.zoom_factor)
        self.view_matrix.setToIdentity()
        self.view_matrix.lookAt(self.camera_pos, QVector3D(0, 0, 0), QVector3D(0, 1, 0))
        self.model_matrix = rotation_to_gl_matrix(self.orientation)

        programs_to_update = self.all_shader_programs
        
        for prog_id in programs_to_update:
            #print(f"BaseGLWidget.paintGL_base: Updating uniforms for program: {prog_id}")
            #if prog_id is None or prog_id == 0: # Check for 0 as it's an invalid program ID
            #    print(f"BaseGLWidget.paintGL_base: Skipping invalid prog_id: {prog_id}")
            #    continue
            
            # print(f"BaseGLWidget.paintGL_base: Using program {prog_id} for common uniforms")
            glUseProgram(prog_id)
            
            # Check and set common uniforms
            # It's good to check locations, though if a shader doesn't use a uniform, glGetUniformLocation returns -1, and glUniform* with -1 is a no-op.
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "projection"), 1, GL_FALSE, self.projection_matrix.data()) # pyArgs: (location, count, transpose, value_ptr)
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "view"),       1, GL_FALSE, self.view_matrix.data())
            glUniformMatrix4fv(glGetUniformLocation(prog_id, "model"),      1, GL_FALSE, self.model_matrix)
            
            glUniform3fv(glGetUniformLocation(prog_id, "lightPos"),   1, [self.light_pos.x(),   self.light_pos.y(),   self.light_pos.z()])
            glUniform3fv(glGetUniformLocation(prog_id, "viewPos"),    1, [self.camera_pos.x(),  self.camera_pos.y(),  self.camera_pos.z()])
            glUniform3fv(glGetUniformLocation(prog_id, "lightColor"), 1, [self.light_color.x(), self.light_color.y(), self.light_color.z()])
        
        self.draw_scene() # Call specific drawing method of derived class

        glUseProgram(0) # Unbind shader after this widget's drawing pass is complete

    def make_labels(self, points, labels ):
        text_object = GLobject( components=[3,2,2] )
        #if not self.trj or self.font_atlas_data is None:
        #    self.num_label_verts = 0
        #    return
        #atom_positions, _, _, atom_enames = self.frames_data[self.current_frame_index].atoms
        #print("atom_enames", atom_enames)
        tile_w = self.font_atlas_data['tile_w']  # width of one character tile in atlas
        tile_h = self.font_atlas_data['tile_h']  # height of one character tile in atlas
        tex_w  = self.font_atlas_data['tex_w']   # full atlas texture width
        all_vertex_data = []
        szx  = 1.0  # width of one character quad in world units
        szy  = 2.0*0.5  # height of one character quad in world units
        yoff = 0.0
        margin_lx = 0.1
        margin_rx = 0.7
        margin_y  = 0.0
        space     = 0.0
        conde_min = 32
        conde_max = 126
        for pos_3d, label in zip(points, labels):
            #ename += "_Hey"  # debug longer text
            symbol = label.decode('utf-8') if isinstance(label, bytes) else str(label) 
            #xoff     = -len(symbol)* szx / 2.0 # centered text
            xoff     = 0 # left aligned text
            for i_char, char in enumerate(symbol):
                code = ord(char)
                if code < conde_min or code > conde_max: continue
                # 1D atlas UVs  - using normalized coordinates
                u0 = (code - conde_min + margin_lx ) * tile_w / tex_w
                u1 = (code - conde_min + margin_rx ) * tile_w / tex_w
                #u1 = u0               + tile_w*(1-margin_x) / tex_w
                v1 = 0.0 + margin_y
                v0 = 1.0 - margin_y
                # Same 3D position for every character in this label
                base_3d = pos_3d + np.array([0.0, 0.0, 0.0])
                # Local screen-space offset for this character
                x = xoff + i_char * (szx*(1.0+space))
                local_offsets = np.array([
                    [ x       , yoff-szy ],
                    [ x + szx , yoff-szy ],
                    [ x + szx , yoff+szy ],
                    [ x       , yoff+szy ]
                ], dtype=np.float32)
                uvs = np.array([[u0, v0], [u1, v0], [u1, v1], [u0, v1]], dtype=np.float32)
                for j in range(4):
                    vertex_data = np.concatenate((base_3d, local_offsets[j], uvs[j]))
                    all_vertex_data.append(vertex_data)
        num_label_verts = len(all_vertex_data)
        #print(f"DEBUG: update_atom_labels_data created {self.num_label_verts} vertices.")
        if num_label_verts > 0:
            # Generate and upload index data for the EBO
            num_quads = num_label_verts // 4
            indices = np.zeros(num_quads * 6, dtype=np.uint32)
            for i in range(num_quads):
                base = i * 4
                indices[i*6:i*6+6] = [base, base + 1, base + 2, base, base + 2, base + 3]
            vertex_data_np = np.array(all_vertex_data, dtype=np.float32)
            #print("vertex_data_np #legend:  pos(x,y,z)   local_offset(x,y)   uv(x,y)\n", vertex_data_np)
            #print("indices: ", indices)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, text_object.vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, vertex_data_np.nbytes, vertex_data_np, GL.GL_DYNAMIC_DRAW)
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, text_object.ebo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
            text_object.nelements = num_quads*6
        text_object.dirty = False
        return text_object

    def use_text_shader(self, labelScale=0.25, textColor=(1.0, 0.0, 0.0, 1.0), offset=(0.5, 0.0)):
        self.use_shader(self.text_shader)
        self.set_default_uniforms() # Set camera matrices
        self.bind_texture('fontAtlas', self.font_atlas_tex, 0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUniform1f(GL.glGetUniformLocation(self.text_shader, "labelScale"), labelScale * self.zoom_factor )
        GL.glUniform4f(GL.glGetUniformLocation(self.text_shader, "textColor"),  textColor[0], textColor[1], textColor[2], textColor[3] )
        GL.glUniform2f(GL.glGetUniformLocation(self.text_shader, "offset"),     offset[0], offset[1] ) 

    def use_simple_shader(self, color=(1.0, 0.0, 1.0, 1.0)):
        self.use_shader(self.simple_shader)
        self.set_default_uniforms()
        color_loc = GL.glGetUniformLocation(self.simple_shader, "color")
        GL.glUniform4f(color_loc, color[0], color[1], color[2], color[3])

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        self.projection_matrix.setToIdentity()
        aspect_ratio = width / height if height > 0 else 1
        self.projection_matrix.perspective(45.0, aspect_ratio, 0.1, 200.0)

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.last_mouse_pos = event.pos()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.LeftButton and self.last_mouse_pos:
            p1 = self._project_to_sphere(self.last_mouse_pos.x(), self.last_mouse_pos.y())
            p2 = self._project_to_sphere(event.x(), event.y())
            
            rotation_axis = np.cross(p1, p2)
            rotation_axis_norm = np.linalg.norm(rotation_axis)

            if rotation_axis_norm > 1e-6:
                rotation_axis /= rotation_axis_norm
                dot_product = np.dot(p1, p2)
                rotation_angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
                delta_rotation = R.from_rotvec(rotation_angle * rotation_axis)
                self.orientation = delta_rotation * self.orientation
            
            self.last_mouse_pos = event.pos()
            self.update()

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom_factor -= delta * 0.01 
        self.zoom_factor = max(1.0, min(self.zoom_factor, 100.0)) # Clamp zoom
        self.update()

    def _project_to_sphere(self, x, y):
        width = self.width()
        height = self.height()
        radius = min(width, height) * 0.4 

        vx = (x - width / 2.0) / radius
        vy = (height / 2.0 - y) / radius 

        vz_squared = 1.0 - vx * vx - vy * vy
        if vz_squared > 0:
            vz = np.sqrt(vz_squared)
        else:
            norm_xy = np.sqrt(vx*vx + vy*vy)
            if norm_xy > 1e-6: 
                vx /= norm_xy
                vy /= norm_xy
            vz = 0.0
        
        p = np.array([vx, vy, vz])
        norm_p = np.linalg.norm(p)
        return p / norm_p if norm_p > 1e-6 else np.array([0.0, 0.0, 1.0])

    def draw_scene(self):
        # This method MUST be overridden by derived classes to perform actual drawing.
        # It is called after common uniforms and transformations are set.
        # The shader program is already in use.
        pass

    def load_texture(self, filepath):
        try:
            img = Image.open(filepath).convert("RGBA")
        except FileNotFoundError:
            print(f"Error: Texture file not found at {filepath}")
            return 0

        img_data = np.array(list(img.getdata()), np.uint8)

        tex_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex_id)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img_data)
        glBindTexture(GL_TEXTURE_2D, 0)
        return tex_id

    def bind_texture(self, uniform_name, tex_id, texture_unit):
        loc = glGetUniformLocation(self.current_shader_program_id, uniform_name)
        if loc != -1:
            glActiveTexture(GL_TEXTURE0 + texture_unit)
            glBindTexture(GL_TEXTURE_2D, tex_id)
            glUniform1i(loc, texture_unit)

# class GLGUI(BaseGLWidget):
#     def __init__(self, parent=None):
#         super().__init__(parent)
#         self.simple_shader = None  # Will hold our default shader program

#     def initializeGL(self):
#         super().initializeGL_base(None, None)
#         # Initialize simple shader


class AppWindow(QMainWindow):
    def __init__(self, parent=None, **kwargs ): # Accept kwargs to pass to QMainWindow if needed
        super().__init__(parent, **kwargs)
        # Common window initialization can go here if any, e.g., default title
        # self.setWindowTitle("Application Window") 

    @classmethod
    def launch(cls, *args, **kwargs):
        """
        Creates and runs the Qt application with an instance of this window class (or a derived class).
        This function will block until the GUI is closed.
        Any *args and **kwargs will be passed to the constructor of `cls`.
        """
        app = QApplication.instance()
        if not app: # Create QApplication if it doesn't exist
            gl_format = QSurfaceFormat()
            gl_format.setVersion(3, 3)
            gl_format.setProfile(QSurfaceFormat.CoreProfile)
            QSurfaceFormat.setDefaultFormat(gl_format)
            # Pass sys.argv if available, otherwise an empty list for robustness
            app = QApplication(sys.argv if hasattr(sys, 'argv') and sys.argv is not None else [])
        
        main_window = cls(*args, **kwargs) # Instantiate the class `cls`
        main_window.show()
        return app.exec_()