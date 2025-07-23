# molecular_browser_app.py
import os
import sys
import numpy as np
import moderngl
from PyQt5.QtWidgets import (QWidget, QGridLayout, QLabel, QPushButton, QVBoxLayout, QFileDialog, QScrollArea,QApplication )
from PyQt5.QtCore import Qt, QSize, pyqtSignal
from PyQt5.QtGui import QPixmap, QImage

from .GLGUI import AppWindow, octahedron_sphere_mesh
from .detailed_viewer import DetailedViewerWindow, create_cylinder_mesh
from ..AtomicSystem import AtomicSystem
from .. import elements
from ..atomicUtils import makeRotMat

class ThumbnailRenderer:
    """Handles off-screen rendering of molecules to textures."""
    def __init__(self, size=(256, 256)):
        self.size = size
        self.ctx = moderngl.create_context(standalone=True, backend='egl')
        print(f"ThumbnailRenderer: Created moderngl context EGL: {self.ctx}, Context info: {self.ctx.info}")
        
        script_dir = os.path.dirname(__file__)
        shader_folder = os.path.join(script_dir, "shaders")
        with open(f"{shader_folder}/instances2.glslv") as f: vert_shader = f.read()
        with open(f"{shader_folder}/sphere.glslf") as f: frag_atom_shader = f.read()
        with open(f"{shader_folder}/cylinder.glslf") as f: frag_bond_shader = f.read()

        self.atom_prog = self.ctx.program(vertex_shader=vert_shader, fragment_shader=frag_atom_shader)
        self.bond_prog = self.ctx.program(vertex_shader=vert_shader, fragment_shader=frag_bond_shader)

        # Debug: Print all exposed attributes
        print("=== ATOM PROGRAM ATTRIBUTES ===")
        for attr_name in ['aPos', 'aNormal', 'instancePosition_model', 'instanceActualSphereRadius', 'instanceColor', 'instanceMatrix']:
            try:
                loc = self.atom_prog[attr_name].location
                print(f"{attr_name}: location {loc}")
            except KeyError:
                print(f"{attr_name}: NOT FOUND")
        
        print("=== BOND PROGRAM ATTRIBUTES ===")
        for attr_name in ['aPos', 'aNormal', 'instancePosition_model', 'instanceActualSphereRadius', 'instanceColor', 'instanceMatrix']:
            try:
                loc = self.bond_prog[attr_name].location
                print(f"{attr_name}: location {loc}")
            except KeyError:
                print(f"{attr_name}: NOT FOUND")



        # Geometry
        sphere_v, sphere_n = octahedron_sphere_mesh(radius=1.3, nsub=2)
        print(f"sphere_v.shape: {sphere_v.shape}, sphere_n.shape: {sphere_n.shape}")
        cyl_v, cyl_n = create_cylinder_mesh(radius=0.15, length=1.0)
        print(f"cyl_v.shape: {cyl_v.shape}, cyl_n.shape: {cyl_n.shape}")

        self.sphere_vbo_vert = self.ctx.buffer(sphere_v.tobytes())
        self.sphere_vbo_norm = self.ctx.buffer(sphere_n.tobytes())
        self.cyl_vbo_vert = self.ctx.buffer(cyl_v.tobytes())
        self.cyl_vbo_norm = self.ctx.buffer(cyl_n.tobytes())

        # Framebuffer
        self.fbo = self.ctx.framebuffer(
            color_attachments=[self.ctx.texture(self.size, 4)],
            depth_attachment=self.ctx.depth_texture(self.size)
        )

    def render_system(self, system: AtomicSystem):
        print(f"ThumbnailRenderer: render_system start for {system.fname if hasattr(system,'fname') else system}.")
        #print( "render_system() na=", len(system.enames) )
        self.fbo.use()
        self.ctx.clear(0.95, 0.95, 0.95)
        self.ctx.enable(moderngl.DEPTH_TEST)

        if system.bonds is None: system.findBonds()

        # Atom data
        n_atoms = len(system.apos)
        atom_pos = system.apos.astype(np.float32)
        atom_radii = np.array([elements.ELEMENT_DICT[e][6] for e in system.enames], dtype=np.float32)
        #print( "elements.ELEMENT_DICT['H'][7] ", elements.ELEMENT_DICT["H"][8] )
        atom_colors = np.array([elements.hex_to_float_rgb(elements.ELEMENT_DICT[e][8]) + (1.0,) for e in system.enames], dtype=np.float32)
        atom_matrices = np.array([np.eye(4, dtype=np.float32) for _ in range(n_atoms)])

        vbo_apos = self.ctx.buffer(atom_pos.tobytes())
        vbo_arad = self.ctx.buffer(atom_radii.tobytes())
        vbo_acol = self.ctx.buffer(atom_colors.tobytes())
        vbo_amat = self.ctx.buffer(atom_matrices.tobytes())
        
        # Sphere shader doesn't use aNormal (uses ray-tracing) but needs instanceColor for fColor
        vao_content_atoms = [
            (self.sphere_vbo_vert, '3f', 'aPos'),
            (vbo_apos, '3f /i', 'instancePosition_model'),
            (vbo_arad, '1f /i', 'instanceActualSphereRadius'),
            (vbo_acol, '4f /i', 'instanceColor'),
            (vbo_amat, '16f /i', 'instanceMatrix')
        ]
        vao_atoms = self.ctx.vertex_array(self.atom_prog, vao_content_atoms)

        # Bond data
        if system.bonds is not None and len(system.bonds) > 0:
            n_bonds = len(system.bonds)
            bond_pos = np.zeros((n_bonds, 3), dtype=np.float32)
            bond_colors = np.zeros((n_bonds, 4), dtype=np.float32)
            bond_matrices = np.zeros((n_bonds, 4, 4), dtype=np.float32)

            for i, (i_a, i_b) in enumerate(system.bonds):
                p1, p2 = system.apos[i_a], system.apos[i_b]
                midpoint, vec, length = (p1 + p2) / 2, p2 - p1, np.linalg.norm(p2 - p1)
                bond_pos[i], bond_colors[i] = midpoint, (atom_colors[i_a] + atom_colors[i_b]) / 2
                # Create 4x4 rotation matrix from 3x3
                #print(f"i_a: {i_a}, i_b: {i_b}, vec: {vec}")
                rot_mat_3x3 = makeRotMat(vec, np.array([0.,0.,1.]))
                rot_mat_4x4 = np.eye(4, dtype=np.float32)
                rot_mat_4x4[:3, :3] = rot_mat_3x3
                # Apply scaling
                scale_mat = np.diag([1., 1., length, 1.]).astype(np.float32)
                bond_matrices[i] = rot_mat_4x4 @ scale_mat

            vbo_bpos = self.ctx.buffer(bond_pos.tobytes())
            vbo_bcol = self.ctx.buffer(bond_colors.tobytes())
            vbo_bmat = self.ctx.buffer(bond_matrices.tobytes())
            vbo_brad = self.ctx.buffer(np.zeros(n_bonds, dtype=np.float32).tobytes()) # Dummy data

            # Cylinder shader uses aNormal and instanceColor but not instanceActualSphereRadius
            vao_content_bonds = [
                (self.cyl_vbo_vert, '3f', 'aPos'),
                (self.cyl_vbo_norm, '3f', 'aNormal'),
                (vbo_bpos, '3f /i', 'instancePosition_model'),
                (vbo_bcol, '4f /i', 'instanceColor'),
                (vbo_bmat, '16f /i', 'instanceMatrix')
            ]
            vao_bonds = self.ctx.vertex_array(self.bond_prog, vao_content_bonds)

        # Camera
        max_dim = np.max(np.ptp(system.apos, axis=0)) if n_atoms > 0 else 10.0
        cam_dist = max_dim * 2.0
        
        # Create perspective projection matrix
        fov, aspect, near, far = 45.0, 1.0, 0.1, cam_dist * 2
        f = 1.0 / np.tan(np.radians(fov) / 2.0)
        projection = np.array([
            [f/aspect, 0, 0, 0],
            [0, f, 0, 0],
            [0, 0, (far+near)/(near-far), (2*far*near)/(near-far)],
            [0, 0, -1, 0]
        ], dtype=np.float32)
        
        # Create look-at matrix
        eye, target, up = np.array([0, 0, cam_dist]), np.array([0, 0, 0]), np.array([0, 1, 0])
        f = (target - eye) / np.linalg.norm(target - eye)
        s = np.cross(f, up) / np.linalg.norm(np.cross(f, up))
        u = np.cross(s, f)
        look_at = np.array([
            [s[0], u[0], -f[0], -np.dot(s, eye)],
            [s[1], u[1], -f[1], -np.dot(u, eye)],
            [s[2], u[2], -f[2], np.dot(f, eye)],
            [0, 0, 0, 1]
        ], dtype=np.float32)
        
        # Set uniforms and render
        for prog in [self.atom_prog, self.bond_prog]:
            prog['projection'].write(projection.tobytes())
            prog['view'].write(look_at.tobytes())
            prog['model'].write(np.eye(4, dtype=np.float32).tobytes())
            prog['lightPos'].value = (cam_dist, cam_dist, cam_dist)
            prog['viewPos'].value = (0, 0, cam_dist)

        vao_atoms.render(instances=n_atoms)
        if system.bonds is not None and len(system.bonds) > 0:
            vao_bonds.render(moderngl.TRIANGLE_STRIP, instances=n_bonds)

        # Read pixels and create QImage
        image_data = self.fbo.read(components=4, alignment=1)
        image = QImage(image_data, self.size[0], self.size[1], QImage.Format_RGBA8888)
        return image.mirrored() # Flip vertically

class ThumbnailGridWidget(QWidget):
    molecule_selected = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_layout = QVBoxLayout(self)
        
        self.scroll_area = QScrollArea(self)
        self.scroll_area.setWidgetResizable(True)
        self.grid_container = QWidget()
        self.grid_layout = QGridLayout(self.grid_container)
        self.scroll_area.setWidget(self.grid_container)

        self.main_layout.addWidget(self.scroll_area)
        
        self.molecules = []
        self.thumb_widgets = []
        self.selected_index = -1
        self.setFocusPolicy(Qt.StrongFocus)

    def populate(self, molecules_data):
        self.clear_grid()
        self.molecules = molecules_data
        self.renderer = ThumbnailRenderer()
        
        cols = 4
        for i, (fname, system) in enumerate(self.molecules):
            row, col = divmod(i, cols)
            
            pixmap = QPixmap.fromImage(self.renderer.render_system(system))
            
            thumb_widget = QWidget()
            thumb_layout = QVBoxLayout(thumb_widget)
            
            img_label = QLabel()
            img_label.setPixmap(pixmap)
            img_label.setFixedSize(QSize(256, 256))
            
            name_label = QLabel(os.path.basename(fname))
            name_label.setAlignment(Qt.AlignCenter)
            
            thumb_layout.addWidget(img_label)
            thumb_layout.addWidget(name_label)
            
            self.grid_layout.addWidget(thumb_widget, row, col)
            self.thumb_widgets.append(thumb_widget)
        
        if self.molecules:
            self.selected_index = 0
            self.update_selection_visuals()

    def clear_grid(self):
        for widget in self.thumb_widgets:
            widget.deleteLater()
        self.thumb_widgets.clear()
        self.molecules.clear()
        self.selected_index = -1

    def update_selection_visuals(self):
        for i, widget in enumerate(self.thumb_widgets):
            if i == self.selected_index:
                widget.setStyleSheet("QWidget { border: 2px solid #0078d7; }")
            else:
                widget.setStyleSheet("QWidget { border: none; }")

    def keyPressEvent(self, event):
        if not self.molecules: return
        
        key = event.key()
        cols = self.grid_layout.columnCount()
        new_index = self.selected_index

        if key == Qt.Key_Right: new_index += 1
        elif key == Qt.Key_Left: new_index -= 1
        elif key == Qt.Key_Down: new_index += cols
        elif key == Qt.Key_Up: new_index -= cols
        elif key in [Qt.Key_Return, Qt.Key_Enter]:
            self.molecule_selected.emit(self.selected_index)
            return
            
        if 0 <= new_index < len(self.molecules):
            self.selected_index = new_index
            self.update_selection_visuals()
            self.scroll_area.ensureWidgetVisible(self.thumb_widgets[self.selected_index])

class MoleculeBrowserApp(AppWindow):
    def __init__(self, dir_path=None):
        super().__init__()
        self.setWindowTitle("Molecule Browser")
        self.setGeometry(100, 100, 1200, 800)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        self.layout = QVBoxLayout(central_widget)

        self.btn_open = QPushButton("Open Directory")
        self.btn_open.clicked.connect(self.open_directory)
        self.layout.addWidget(self.btn_open)
            
        self.thumbnail_grid = ThumbnailGridWidget()
        self.thumbnail_grid.molecule_selected.connect(self.open_detailed_viewer)
        self.layout.addWidget(self.thumbnail_grid)
        
        self.open_viewers = []

        if dir_path:
            print( f"MoleculeBrowserApp.open_directory({dir_path})" )
            self.open_directory(dir_path)

    def open_directory(self, dir_path=None):
        if not dir_path:
            dir_path = QFileDialog.getExistingDirectory(self, "Select Directory")
        if dir_path:
            molecules = []
            for fname in sorted(os.listdir(dir_path)):
                if fname.lower().endswith('.xyz'):
                    fpath = os.path.join(dir_path, fname)
                    try:
                        system = AtomicSystem(fname=fpath)
                        system.fname = fpath # Store the full path
                        molecules.append((fname, system))
                    except Exception as e:
                        print(f"Error loading {fname}: {e}")
            self.thumbnail_grid.populate(molecules)

    def open_detailed_viewer(self, index):
        if 0 <= index < len(self.thumbnail_grid.molecules):
            _, system = self.thumbnail_grid.molecules[index]
            viewer = DetailedViewerWindow(system)
            viewer.show()
            self.open_viewers.append(viewer)

if __name__ == '__main__':
    from PyQt5.QtGui import QSurfaceFormat
    fmt = QSurfaceFormat()
    fmt.setVersion(3, 3)
    fmt.setProfile(QSurfaceFormat.CoreProfile)
    QSurfaceFormat.setDefaultFormat(fmt)
    print('Set default QSurfaceFormat to OpenGL 3.3 Core')

    # Run Like this:
    #   python -m pyBall.GUI.molecular_browser

    # relative path to the directory with the molecules
    # /home/prokop/git/FireCore/cpp/common_resources/xyz/
    dir_path = "../../cpp/common_resources/xyz/"
    dir_path = "../../cpp/common_resources/xyz_mini/"
    this_path = os.path.dirname(os.path.abspath(__file__))
    # resolve the full absolute path
    dir_path = os.path.abspath(os.path.join(this_path, dir_path))
    
    app = QApplication(sys.argv)
    window = MoleculeBrowserApp(dir_path=dir_path)
    window.show()
    sys.exit(app.exec_())