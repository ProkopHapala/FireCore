import numpy as np
from scipy.spatial.transform import Rotation as R

from PyQt5.QtGui import QMatrix4x4, QVector3D # These are used for convenience, can be replaced by numpy if strict PyQt independence is needed for math types

from OpenGL.GL import (
    glClear, glClearColor, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT,
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
    glPixelStorei, GL_UNPACK_ALIGNMENT, GL_CULL_FACE, glTexImage2D, glDisable,
    GL_SHADER_STORAGE_BUFFER, glBindBufferBase
)
from OpenGL.GL.shaders import compileProgram, compileShader
import OpenGL.GL as GL
import ctypes

alpha_blend_modes={
    "standard"    :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA),
    "standard2"    :(GL_FUNC_ADD, GL_ONE, GL_ONE_MINUS_SRC_ALPHA),
    "additive"   :(GL_FUNC_ADD, GL_SRC_ALPHA, GL_ONE),
    "subtractive" :(GL_FUNC_REVERSE_SUBTRACT, GL_SRC_ALPHA, GL_ONE),
}

def rotation_to_gl_matrix( rotation_obj):
    mat3x3 = rotation_obj.as_matrix()
    mat4x4 = np.eye(4, dtype=np.float32)
    mat4x4[:3, :3] = mat3x3
    return mat4x4.copy()

def set_ogl_blend_mode(mode, depth_test=True):
    glEnable(GL_BLEND)
    if len(mode) == 3:
        glBlendEquation(mode[0])
        glBlendFunc(mode[1], mode[2])
    else:
        glBlendEquationSeparate(mode[0], mode[3])
        glBlendFuncSeparate    (mode[1], mode[2], mode[4], mode[5])
    if depth_test:
        glDisable(GL_DEPTH_TEST)

def disable_blend( depth_test=True):
    glDisable(GL_BLEND)
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

def upload_buffer( index, buffer_id, data, mode=GL_DYNAMIC_DRAW):
    glBindBuffer    (GL_SHADER_STORAGE_BUFFER, buffer_id)
    glBufferData    (GL_SHADER_STORAGE_BUFFER, data.nbytes, data, mode)
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, index, buffer_id)
    glBindBuffer    (GL_SHADER_STORAGE_BUFFER, 0)

def compile_shader_program(vertex_shader_src, fragment_shader_src):
    try:
        program = compileProgram(
            compileShader(vertex_shader_src, GL_VERTEX_SHADER),
            compileShader(fragment_shader_src, GL_FRAGMENT_SHADER)
        )
        return program
    except Exception as e:
        print(f"Shader compilation/linking error: {e}")
        return None

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
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indices.nbytes, indices, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)

    def upload_vbo(self, vertices):
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

    def upload_vbo_ebo(self, indices, vertices ):
        self.upload_ebo(indices)
        self.upload_vbo(vertices)
        self.nelements = indices.size

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
        if self.has_indices:
            glDrawElementsInstanced(mode, self.index_count, GL_UNSIGNED_INT, None, num_instances)
        else:
            glDrawArraysInstanced(mode, 0, self.vertex_count, num_instances)

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
        for vbo in self.vbos:
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

    def update_list(self, buffs, n=-1):
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
        self.associated_mesh.draw_instanced(self.num_instances, mode)
        glBindVertexArray(0)


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
    da = (c2 - c1) / nsub
    db = (c3 - c1) / nsub
    for i in range(nsub):      # Iterate along the c1-c2 edge direction
        for j in range(nsub - i):  # Iterate along the c1-c3 edge direction
            p00 = c1 + da*i     + db*j
            p10 = c1 + da*(i+1) + db*j
            p01 = c1 + da*i     + db*(j+1)
            append_normalized([p00, p10, p01], radius, vertices_list, normals_list)
            if j < nsub - 1 - i:
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
        [0, 2, 4], [0, 4, 3], [0, 3, 5], [0, 5, 2],
        [1, 4, 2], [1, 3, 4], [1, 5, 3], [1, 2, 5]
    ], dtype=np.int32)

    vertices_list = []
    normals_list = []
    for i in range(8):
        c1 = pv_np[octa_face_indices[i,0]]
        c2 = pv_np[octa_face_indices[i,1]]
        c3 = pv_np[octa_face_indices[i,2]]
        octahedron_sphere_face(c1, c2, c3, vertices_list, normals_list, nsub=nsub, radius=radius)

    return np.array(vertices_list, dtype=np.float32), np.array(normals_list, dtype=np.float32), None

class OGLSystem:
    """
    A system for managing OpenGL shader programs and other global GL states.
    """
    def __init__(self):
        self.shader_programs = {}

    def load_shader_program(self, name, vertex_shader_src, fragment_shader_src):
        """
        Compile and store an OpenGL shader program.
        Args:
            name (str): A unique name for this shader program.
            vertex_shader_src (str): Vertex shader source code.
            fragment_shader_src (str): Fragment shader source code.
        """
        program = compile_shader_program(vertex_shader_src, fragment_shader_src)
        if program:
            self.shader_programs[name] = program
            print(f"OGLSystem::load_shader_program() Successfully loaded shader '{name}'")
            return True
        return False

    def get_shader_program(self, name):
        """
        Get a compiled shader program by name.
        """
        return self.shader_programs.get(name)
