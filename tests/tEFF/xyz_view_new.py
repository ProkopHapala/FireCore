import sys
import numpy as np
from scipy.spatial.transform import Rotation as R
import argparse

from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget,
                             QSlider, QLabel, QOpenGLWidget)
from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QSurfaceFormat, QMatrix4x4, QVector3D

from OpenGL.GL import *
from OpenGL.GL.shaders import compileProgram, compileShader

# It's good practice to keep shaders in separate files or as multi-line strings
# For simplicity here, they are embedded as strings.

VERTEX_SHADER_SOURCE = """
#version 330 core

layout (location = 0) in vec3 aPos;      // Vertex position of the base sphere
layout (location = 1) in vec3 aNormal;   // Vertex normal of the base sphere

// Per-instance attributes
layout (location = 2) in vec3 instancePosition; // Center of the sphere instance
layout (location = 3) in float instanceRadius;  // Radius of the sphere instance
layout (location = 4) in vec4 instanceColor;    // Color (RGBA) of the sphere instance

out vec3 FragPos;
out vec3 Normal;
out vec4 AtomColor;

uniform mat4 projection;
uniform mat4 view;
uniform mat4 model; // Overall model orientation (trackball)

void main()
{
    // Create a model matrix for this specific instance (sphere)
    // It translates to instancePosition and scales by instanceRadius
    mat4 instanceModelMatrix = mat4(1.0);
    instanceModelMatrix[3][0] = instancePosition.x;
    instanceModelMatrix[3][1] = instancePosition.y;
    instanceModelMatrix[3][2] = instancePosition.z;
    
    instanceModelMatrix[0][0] = instanceRadius;
    instanceModelMatrix[1][1] = instanceRadius;
    instanceModelMatrix[2][2] = instanceRadius;

    // Apply the overall model rotation first, then the instance-specific transform
    mat4 finalModelMatrix = model * instanceModelMatrix;

    FragPos = vec3(finalModelMatrix * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(finalModelMatrix))) * aNormal; // Transform normals correctly
    AtomColor = instanceColor;

    gl_Position = projection * view * finalModelMatrix * vec4(aPos, 1.0);
}
"""

FRAGMENT_SHADER_SOURCE = """
#version 330 core

out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec4 AtomColor; // Received from vertex shader (includes opacity)

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;

void main()
{
    // Ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Specular
    float specularStrength = 0.9;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 64); // Shininess
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * AtomColor.rgb;
    FragColor = vec4(result, AtomColor.a); // Use alpha from AtomColor
}
"""

sys.path.append("../../") # To find pyBall
from pyBall import atomicUtils as au

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


class ModernGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.traj = []
        self.current_frame_index = 0
        self.opacity = 0.5
        self.zoom_factor = 20.0  # Initial zoom

        self.orientation = R.identity()
        self.last_mouse_pos = None

        self.shader_program = None
        self.sphere_vao = None
        self.sphere_vbo_vertices = None
        self.sphere_vbo_normals = None
        self.sphere_ebo = None
        self.instance_vbo_positions = None
        self.instance_vbo_radii = None
        self.instance_vbo_colors = None
        
        self.sphere_vertices = None
        self.sphere_normals  = None
        self.sphere_indices  = None
        self.num_sphere_indices = 0

        self.num_instances = 0
        self.instance_data_dirty = True # Flag to update instance VBOs

        # Camera and light
        self.view_matrix = QMatrix4x4()
        self.projection_matrix = QMatrix4x4()
        self.light_pos = QVector3D(15.0, 15.0, 30.0)
        self.light_color = QVector3D(1.0, 1.0, 1.0)
        self.camera_pos = QVector3D(0,0,self.zoom_factor) # Will be updated by zoom

    def initializeGL(self):
        # Modern OpenGL context should be requested via QSurfaceFormat in main
        print(f"OpenGL Version: {glGetString(GL_VERSION).decode()}")
        print(f"GLSL Version: {glGetString(GL_SHADING_LANGUAGE_VERSION).decode()}")

        glClearColor(1.0, 1.0, 1.0, 1.0)  # White background
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        #glEnable(GL_CULL_FACE) # Optional: for slight performance gain if backfaces are not needed

        # Compile shaders
        try:
            self.shader_program = compileProgram(
                compileShader(VERTEX_SHADER_SOURCE, GL_VERTEX_SHADER),
                compileShader(FRAGMENT_SHADER_SOURCE, GL_FRAGMENT_SHADER)
            )
        except RuntimeError as e:
            print(f"Shader Compilation Error: {e}")
            sys.exit(1)

        # Generate base sphere mesh
        self.sphere_vertices, self.sphere_normals, self.sphere_indices = create_sphere_mesh(radius=1.0) # Base radius is 1, will be scaled per instance
        self.num_sphere_indices = len(self.sphere_indices)

        # Create VAO
        self.sphere_vao = glGenVertexArrays(1)
        glBindVertexArray(self.sphere_vao)

        # VBO for base sphere vertices
        self.sphere_vbo_vertices = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.sphere_vbo_vertices)
        glBufferData(GL_ARRAY_BUFFER, self.sphere_vertices.nbytes, self.sphere_vertices, GL_STATIC_DRAW)
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)

        # VBO for base sphere normals
        self.sphere_vbo_normals = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.sphere_vbo_normals)
        glBufferData(GL_ARRAY_BUFFER, self.sphere_normals.nbytes, self.sphere_normals, GL_STATIC_DRAW)
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(1)
        
        # EBO for sphere indices
        self.sphere_ebo = glGenBuffers(1)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.sphere_ebo)
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, self.sphere_indices.nbytes, self.sphere_indices, GL_STATIC_DRAW)

        # VBOs for instance data (will be filled later)
        self.instance_vbo_positions = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_positions)
        # glBufferData will be called in update_instance_data
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(2)
        glVertexAttribDivisor(2, 1) # Tell OpenGL this is an instanced vertex attribute.

        self.instance_vbo_radii = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_radii)
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(3)
        glVertexAttribDivisor(3, 1)

        self.instance_vbo_colors = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_colors)
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 0, None) # RGBA
        glEnableVertexAttribArray(4)
        glVertexAttribDivisor(4, 1)

        glBindBuffer(GL_ARRAY_BUFFER, 0) # Unbind VBO
        glBindVertexArray(0)             # Unbind VAO

    def update_instance_data(self):
        if not self.traj or not (0 <= self.current_frame_index < len(self.traj)):
            self.num_instances = 0
            return

        es, apos, qs, rs, comment = self.traj[self.current_frame_index]
        self.num_instances = len(es)

        if self.num_instances == 0:
            return

        instance_positions = np.array(apos, dtype=np.float32).flatten()
        
        instance_radii_list = []
        for r_val in rs:
            instance_radii_list.append(r_val if not np.isnan(r_val) else 0.5)
        instance_radii = np.array(instance_radii_list, dtype=np.float32)

        instance_colors_list = []
        for atom_symbol in es:
            try:
                element_data = au.elements.ELEMENT_DICT[atom_symbol]
                hex_color_str = element_data[au.elements.index_color]
                r, g, b = au.elements.hex_to_float_rgb(hex_color_str)
            except KeyError:
                r, g, b = 0.5, 0.5, 0.5 # Gray for unknowns
            instance_colors_list.extend([r, g, b, self.opacity])
        instance_colors = np.array(instance_colors_list, dtype=np.float32)

        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_positions)
        glBufferData(GL_ARRAY_BUFFER, instance_positions.nbytes, instance_positions, GL_DYNAMIC_DRAW)
        
        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_radii)
        glBufferData(GL_ARRAY_BUFFER, instance_radii.nbytes, instance_radii, GL_DYNAMIC_DRAW)

        glBindBuffer(GL_ARRAY_BUFFER, self.instance_vbo_colors)
        glBufferData(GL_ARRAY_BUFFER, instance_colors.nbytes, instance_colors, GL_DYNAMIC_DRAW)
        
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        self.instance_data_dirty = False


    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        if self.shader_program is None or self.num_sphere_indices == 0:
            return

        glUseProgram(self.shader_program)

        # Update instance data if needed
        if self.instance_data_dirty:
            self.update_instance_data()
        
        if self.num_instances == 0:
            glUseProgram(0)
            return

        # Camera / View transformation
        self.camera_pos = QVector3D(0, 0, self.zoom_factor)
        self.view_matrix.setToIdentity()
        self.view_matrix.lookAt(self.camera_pos, QVector3D(0, 0, 0), QVector3D(0, 1, 0))

        # Model transformation (from trackball)
        model_matrix_np = self._rotation_to_gl_matrix(self.orientation) # This is column-major
        # QMatrix4x4 constructor expects row-major, or provide data pointer for column-major
        model_qmatrix = QMatrix4x4(model_matrix_np.reshape(4,4).T.flatten().tolist())


        # Set uniforms
        glUniformMatrix4fv(glGetUniformLocation(self.shader_program, "projection"), 1, GL_FALSE, self.projection_matrix.data())
        glUniformMatrix4fv(glGetUniformLocation(self.shader_program, "view"), 1, GL_FALSE, self.view_matrix.data())
        glUniformMatrix4fv(glGetUniformLocation(self.shader_program, "model"), 1, GL_FALSE, model_qmatrix.data())
        
        glUniform3fv(glGetUniformLocation(self.shader_program, "lightPos"), 1, [self.light_pos.x(), self.light_pos.y(), self.light_pos.z()])
        glUniform3fv(glGetUniformLocation(self.shader_program, "viewPos"), 1, [self.camera_pos.x(), self.camera_pos.y(), self.camera_pos.z()])
        glUniform3fv(glGetUniformLocation(self.shader_program, "lightColor"), 1, [self.light_color.x(), self.light_color.y(), self.light_color.z()])

        # Draw
        glBindVertexArray(self.sphere_vao)
        #glDrawArraysInstanced(GL_TRIANGLES, 0, self.num_sphere_vertices, self.num_instances)
        glDrawElementsInstanced(GL_TRIANGLES, self.num_sphere_indices, GL_UNSIGNED_INT, None, self.num_instances)
        glBindVertexArray(0)
        glUseProgram(0)

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
            self.update() # Trigger repaint

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

    def _rotation_to_gl_matrix(self, rotation_obj):
        mat3x3 = rotation_obj.as_matrix()
        mat4x4 = np.eye(4, dtype=np.float32)
        mat4x4[:3, :3] = mat3x3
        return mat4x4.flatten(order='F') # OpenGL expects column-major

    def load_trajectory(self, filename):
        self.traj = au.load_xyz_movie(filename)
        self.traj = au.trj_to_ename(self.traj)
        self.traj = au.trj_fill_radius(self.traj, bVdw=False, rFactor=1.0) # Use covalent radii
        self.current_frame_index = 0
        self.instance_data_dirty = True
        self.update()
        return len(self.traj)

    def set_frame(self, frame_idx):
        if 0 <= frame_idx < len(self.traj):
            self.current_frame_index = frame_idx
            self.instance_data_dirty = True
            self.update()

    def set_opacity(self, opacity_percent):
        self.opacity = opacity_percent / 100.0
        self.instance_data_dirty = True # Need to update colors VBO
        self.update()


class ModernMainWindow(QMainWindow):
    def __init__(self, filepath=None):
        super().__init__()
        self.setWindowTitle("Modern OpenGL Molecular Viewer")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.gl_widget = ModernGLWidget()
        layout.addWidget(self.gl_widget, 1) # <-- Add stretch factor here

        self.frame_slider = QSlider(Qt.Horizontal)
        self.frame_slider.valueChanged.connect(self.gl_widget.set_frame)
        layout.addWidget(QLabel("Frame:"))
        layout.addWidget(self.frame_slider)

        self.opacity_slider = QSlider(Qt.Horizontal)
        self.opacity_slider.setMinimum(0)
        self.opacity_slider.setMaximum(100)
        self.opacity_slider.setValue(int(self.gl_widget.opacity * 100))
        self.opacity_slider.valueChanged.connect(self.gl_widget.set_opacity)
        layout.addWidget(QLabel("Opacity:"))
        layout.addWidget(self.opacity_slider)
        
        if filepath:
            num_frames = self.gl_widget.load_trajectory(filepath)
            if num_frames > 0:
                self.frame_slider.setMinimum(0)
                self.frame_slider.setMaximum(num_frames - 1)
            else:
                print(f"Warning: No frames loaded from {filepath}")
        else:
            print("No trajectory file provided on startup.")


if __name__ == '__main__':
    app = QApplication(sys.argv)

    # --- Crucial for Modern OpenGL with PyQt ---
    # Request a specific OpenGL version and profile.
    # 3.3 Core is a good baseline for modern features.
    gl_format = QSurfaceFormat()
    gl_format.setVersion(3, 3)
    gl_format.setProfile(QSurfaceFormat.CoreProfile)
    # gl_format.setSamples(4) # Optional: for multisampling/antialiasing
    QSurfaceFormat.setDefaultFormat(gl_format)
    # ---

    parser = argparse.ArgumentParser(description="Modern OpenGL Molecular Viewer")
    parser.add_argument("-f", "--file", type=str, help="Path to the XYZ trajectory file", default="/home/prokophapala/git/FireCore/tests/tEFF/relaxation.xyz")
    args = parser.parse_args()

    main_window = ModernMainWindow(filepath=args.file)
    main_window.show()
    sys.exit(app.exec_())
