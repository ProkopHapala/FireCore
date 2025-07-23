# detailed_viewer.py
import sys
import numpy as np
from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget
from OpenGL.GL import *

#sys.path.append(".") # Ensure parent directory is in path
from .GLGUI import BaseGLWidget, AppWindow, InstancedData, Mesh, octahedron_sphere_mesh
from ..AtomicSystem import AtomicSystem, elements
from ..atomicUtils import makeRotMat

def create_cylinder_mesh(radius=1.0, length=1.0, segments=16):
    """Generates vertices and normals for a cylinder mesh oriented along the Z-axis."""
    vertices = []
    normals = []
    
    # Cylinder wall
    for i in range(segments + 1):
        angle = i * 2 * np.pi / segments
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        
        # Bottom vertex
        vertices.extend([x, y, -length / 2])
        normals.extend([x/radius, y/radius, 0])
        
        # Top vertex
        vertices.extend([x, y, length / 2])
        normals.extend([x/radius, y/radius, 0])

    # Convert to flat numpy array for GL_TRIANGLE_STRIP
    final_vertices = np.array(vertices, dtype=np.float32)
    final_normals = np.array(normals, dtype=np.float32)

    return final_vertices, final_normals

class DetailedGLWidget(BaseGLWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.atom_shader = None
        self.bond_shader = None
        self.atom_instances = None
        self.bond_instances = None
        self.cylinder_mesh = None


    @property
    def all_shader_programs(self):
        return [prog for prog in [self.atom_shader, self.bond_shader] if prog]

    def initializeGL(self):
        print( "DetailedGLWidget.initializeGL()" )
        import os
        script_dir = os.path.dirname(__file__)
        shader_folder = os.path.join(script_dir, "shaders")
        with open(f"{shader_folder}/instances2.glslv") as f: vert_shader = f.read()
        with open(f"{shader_folder}/sphere.glslf") as f: frag_atom_shader = f.read()
        with open(f"{shader_folder}/cylinder.glslf") as f: frag_bond_shader = f.read()

        self.atom_shader = self.compile_shader_program(vert_shader, frag_atom_shader)
        self.bond_shader = self.compile_shader_program(vert_shader, frag_bond_shader)
        
        super().initializeGL_base(None, None) # Base init after shaders are compiled

        # Setup atom instances
        self.atom_instances = InstancedData(base_attrib_location=2)
        self.atom_instances.associate_mesh(self.default_sphere_mesh)
        inst_attribs = [
            ("positions",  0, 3), ("radii", 1, 1), ("colors", 2, 4), ("instanceMatrix", 3, 16) # mat4
        ]
        self.atom_instances.setup_instance_vbos(inst_attribs)

        # Setup bond instances (Cylinders)
        cyl_v, cyl_n = create_cylinder_mesh(radius=0.15, length=1.0) # Base cylinder
        self.cylinder_mesh = Mesh(vertices=cyl_v, normals=cyl_n)
        self.cylinder_mesh.setup_buffers()

        self.bond_instances = InstancedData(base_attrib_location=2)
        self.bond_instances.associate_mesh(self.cylinder_mesh)
        self.bond_instances.setup_instance_vbos(inst_attribs)
        if hasattr(self, '_pending_system'):
            self.set_molecule(self._pending_system)
            del self._pending_system

    def set_molecule(self, system: AtomicSystem):
        print( "DetailedGLWidget.set_molecule()" )
        if system.bonds is None:
            system.findBonds()

        # 1. Prepare Atom Data
        n_atoms = len(system.apos)
        atom_pos = system.apos.astype(np.float32)
        atom_radii = np.array([elements.ELEMENT_DICT[e][6] for e in system.enames], dtype=np.float32)
        atom_colors = np.array([elements.hex_to_float_rgb(elements.ELEMENT_DICT[e][8]) + (1.0,) for e in system.enames], dtype=np.float32)
        atom_matrices = np.array([np.eye(4, dtype=np.float32) for _ in range(n_atoms)])

        self.atom_instances.update({
            "positions": atom_pos, "radii": atom_radii, "colors": atom_colors, "instanceMatrix": atom_matrices
        })

        # 2. Prepare Bond Data
        if system.bonds is not None and len(system.bonds) > 0:
            n_bonds = len(system.bonds)
            bond_pos = np.zeros((n_bonds, 3), dtype=np.float32)
            bond_colors = np.zeros((n_bonds, 4), dtype=np.float32)
            bond_matrices = np.zeros((n_bonds, 4, 4), dtype=np.float32)

            for i, (i_a, i_b) in enumerate(system.bonds):
                p1, p2 = system.apos[i_a], system.apos[i_b]
                midpoint = (p1 + p2) / 2
                bond_vector = p2 - p1
                length = np.linalg.norm(bond_vector)

                bond_pos[i] = midpoint
                # Average color of the two atoms
                avg_color = (atom_colors[i_a] + atom_colors[i_b]) / 2
                bond_colors[i] = avg_color

                # Create rotation and scale matrix
                scale_mat = np.diag([1.0, 1.0, length, 1.0])
                rot_mat = makeRotMat(bond_vector, np.array([0., 0., 1.])) # Align Z-axis with bond
                bond_matrices[i] = rot_mat @ scale_mat
            
            self.bond_instances.update({
                "positions": bond_pos, "radii": np.zeros(n_bonds, dtype=np.float32), "colors": bond_colors, "instanceMatrix": bond_matrices
            })
        self.update()

    def draw_scene(self):
        # Draw atoms
        glUseProgram(self.atom_shader)
        self.atom_instances.draw(mode=GL_TRIANGLES)

        # Draw bonds
        glUseProgram(self.bond_shader)
        self.bond_instances.associated_mesh.vertex_count = 34 # Set the correct count for GL_TRIANGLE_STRIP
        self.bond_instances.draw(mode=GL_TRIANGLE_STRIP)

    def cleanupGL(self):
        super().cleanupGL_base()
        if self.atom_shader: glDeleteProgram(self.atom_shader)
        if self.bond_shader: glDeleteProgram(self.bond_shader)
        if self.atom_instances: self.atom_instances.cleanup()
        if self.bond_instances: self.bond_instances.cleanup()
        if self.cylinder_mesh: self.cylinder_mesh.cleanup()

class DetailedViewerWindow(QMainWindow):
    def __init__(self, system: AtomicSystem):
        super().__init__()
        self.setWindowTitle(f"Molecule Viewer - {system.fname.split('/')[-1]}")
        self.setGeometry(150, 150, 800, 600)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        self.gl_widget = DetailedGLWidget()
        self.gl_widget._pending_system = system
        layout.addWidget(self.gl_widget)
        
        