#!/usr/bin/env python3
"""
Interactive Matplotlib-based Molecular Viewer
Supports loading, rotating, and duplicating molecules for symmetry arrangements

# Claude Sonet 4.5

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyBall.AtomicSystem import AtomicSystem
from pyBall import atomicUtils as au

class MoleculeViewer:
    def __init__(self, fname=None, lattice_a=32.7):
        self.lattice_a = lattice_a  # Lattice constant in Angstroms
        self.molecules = []  # List of (apos, enames) tuples
        self.colors = []     # Colors for each molecule
        self.pivot = np.array([0.0, 0.0, 0.0])  # Rotation pivot
        
        # Rotation angles (in degrees)
        self.rot_x = 0.0
        self.rot_y = 0.0
        self.rot_z = 0.0
        
        # Current molecule index
        self.current_mol_idx = 0
        
        # Load initial molecule if provided
        if fname is not None:
            self.load_molecule(fname)
        
        self.setup_gui()
    
    def load_molecule(self, fname):
        """Load a molecule from file"""
        print(f"Loading molecule from {fname}")
        sys = AtomicSystem(fname=fname)
        
        if sys.apos is None or len(sys.apos) == 0:
            raise ValueError(f"Failed to load molecule from {fname}")
        
        # Store original molecule
        self.original_apos = sys.apos.copy()
        self.original_enames = sys.enames.copy()
        
        # Center molecule at origin
        self.center_of_mass = np.mean(self.original_apos, axis=0)
        self.original_apos -= self.center_of_mass[None, :]
        
        # Initialize with single molecule
        self.molecules = [(self.original_apos.copy(), self.original_enames.copy())]
        self.colors = ['blue']
        self.current_mol_idx = 0
        
        # Reset rotations
        self.rot_x = 0.0
        self.rot_y = 0.0
        self.rot_z = 0.0
        self.pivot = np.array([0.0, 0.0, 0.0])
        
        print(f"Loaded {len(self.original_apos)} atoms")
        return sys
    
    def get_rotation_matrix(self, angles_deg):
        """Compute combined rotation matrix from Euler angles (degrees)"""
        ax, ay, az = np.radians(angles_deg)
        
        # Rotation around X
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(ax), -np.sin(ax)],
            [0, np.sin(ax), np.cos(ax)]
        ])
        
        # Rotation around Y
        Ry = np.array([
            [np.cos(ay), 0, np.sin(ay)],
            [0, 1, 0],
            [-np.sin(ay), 0, np.cos(ay)]
        ])
        
        # Rotation around Z
        Rz = np.array([
            [np.cos(az), -np.sin(az), 0],
            [np.sin(az), np.cos(az), 0],
            [0, 0, 1]
        ])
        
        # Combined rotation: Rz * Ry * Rx
        return Rz @ Ry @ Rx
    
    def apply_rotation(self):
        """Apply current rotation to the current molecule"""
        if len(self.molecules) == 0:
            return
        
        rot_mat = self.get_rotation_matrix([self.rot_x, self.rot_y, self.rot_z])
        
        # Rotate around pivot
        apos = self.original_apos.copy()
        apos -= self.pivot[None, :]
        apos = (rot_mat @ apos.T).T
        apos += self.pivot[None, :]
        
        # Update current molecule
        self.molecules[self.current_mol_idx] = (apos, self.original_enames.copy())
    
    def duplicate_rotational(self, n_copies, axis='z'):
        """Duplicate molecule rotationally around pivot"""
        if len(self.molecules) == 0:
            return
        
        # Get current molecule
        apos_base, enames_base = self.molecules[self.current_mol_idx]
        
        # Clear all molecules except the base
        self.molecules = [(apos_base.copy(), enames_base.copy())]
        self.colors = ['blue']
        
        # Generate n_copies around the specified axis
        angle_step = 360.0 / n_copies
        
        color_cycle = ['blue', 'red', 'green', 'orange', 'purple', 'cyan', 'magenta', 'yellow']
        
        for i in range(1, n_copies):
            angle = i * angle_step
            
            # Create rotation matrix for this copy
            if axis == 'x':
                rot_mat = self.get_rotation_matrix([angle, 0, 0])
            elif axis == 'y':
                rot_mat = self.get_rotation_matrix([0, angle, 0])
            else:  # z
                rot_mat = self.get_rotation_matrix([0, 0, angle])
            
            # Apply rotation around pivot
            apos_new = apos_base.copy()
            apos_new -= self.pivot[None, :]
            apos_new = (rot_mat @ apos_new.T).T
            apos_new += self.pivot[None, :]
            
            self.molecules.append((apos_new, enames_base.copy()))
            self.colors.append(color_cycle[i % len(color_cycle)])
        
        print(f"Created {n_copies} rotational copies around {axis}-axis")
    
    def setup_gui(self):
        """Setup matplotlib GUI"""
        self.fig = plt.figure(figsize=(14, 10))
        
        # 3D plot
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.fig.subplots_adjust(left=0.05, bottom=0.35, right=0.95, top=0.95)
        
        # Sliders for rotation
        ax_rot_x = self.fig.add_axes([0.15, 0.25, 0.7, 0.02])
        ax_rot_y = self.fig.add_axes([0.15, 0.20, 0.7, 0.02])
        ax_rot_z = self.fig.add_axes([0.15, 0.15, 0.7, 0.02])
        
        self.slider_x = Slider(ax_rot_x, 'Rot X (°)', -180, 180, valinit=0, valstep=1)
        self.slider_y = Slider(ax_rot_y, 'Rot Y (°)', -180, 180, valinit=0, valstep=1)
        self.slider_z = Slider(ax_rot_z, 'Rot Z (°)', -180, 180, valinit=0, valstep=1)
        
        # Pivot position
        ax_pivot_x = self.fig.add_axes([0.15, 0.10, 0.7, 0.02])
        ax_pivot_y = self.fig.add_axes([0.15, 0.05, 0.7, 0.02])
        ax_pivot_z = self.fig.add_axes([0.15, 0.00, 0.7, 0.02])
        
        self.slider_pivot_x = Slider(ax_pivot_x, 'Pivot X', -50, 50, valinit=0, valstep=0.5)
        self.slider_pivot_y = Slider(ax_pivot_y, 'Pivot Y', -50, 50, valinit=0, valstep=0.5)
        self.slider_pivot_z = Slider(ax_pivot_z, 'Pivot Z', -50, 50, valinit=0, valstep=0.5)
        
        # Buttons
        ax_btn_load = self.fig.add_axes([0.05, 0.90, 0.10, 0.04])
        ax_btn_reset = self.fig.add_axes([0.05, 0.85, 0.10, 0.04])
        ax_btn_dup3 = self.fig.add_axes([0.20, 0.90, 0.10, 0.04])
        ax_btn_dup6 = self.fig.add_axes([0.20, 0.85, 0.10, 0.04])
        ax_btn_save = self.fig.add_axes([0.35, 0.90, 0.10, 0.04])
        ax_btn_clear = self.fig.add_axes([0.35, 0.85, 0.10, 0.04])
        
        self.btn_load = Button(ax_btn_load, 'Load XYZ')
        self.btn_reset = Button(ax_btn_reset, 'Reset')
        self.btn_dup3 = Button(ax_btn_dup3, '3-fold Dup')
        self.btn_dup6 = Button(ax_btn_dup6, '6-fold Dup')
        self.btn_save = Button(ax_btn_save, 'Save XYZ')
        self.btn_clear = Button(ax_btn_clear, 'Clear Copies')
        
        # Text box for file path
        ax_text = self.fig.add_axes([0.50, 0.90, 0.45, 0.04])
        self.text_file = TextBox(ax_text, 'File:', initial='')
        
        # Connect callbacks
        self.slider_x.on_changed(self.on_rotation_change)
        self.slider_y.on_changed(self.on_rotation_change)
        self.slider_z.on_changed(self.on_rotation_change)
        self.slider_pivot_x.on_changed(self.on_pivot_change)
        self.slider_pivot_y.on_changed(self.on_pivot_change)
        self.slider_pivot_z.on_changed(self.on_pivot_change)
        
        self.btn_load.on_clicked(self.on_load)
        self.btn_reset.on_clicked(self.on_reset)
        self.btn_dup3.on_clicked(lambda event: self.on_duplicate(event, 3))
        self.btn_dup6.on_clicked(lambda event: self.on_duplicate(event, 6))
        self.btn_save.on_clicked(self.on_save)
        self.btn_clear.on_clicked(self.on_clear)
        
        self.update_plot()
    
    def on_rotation_change(self, val):
        """Callback for rotation sliders"""
        self.rot_x = self.slider_x.val
        self.rot_y = self.slider_y.val
        self.rot_z = self.slider_z.val
        self.apply_rotation()
        self.update_plot()
    
    def on_pivot_change(self, val):
        """Callback for pivot sliders"""
        self.pivot[0] = self.slider_pivot_x.val
        self.pivot[1] = self.slider_pivot_y.val
        self.pivot[2] = self.slider_pivot_z.val
        self.apply_rotation()
        self.update_plot()
    
    def on_load(self, event):
        """Load molecule from file"""
        fname = self.text_file.text.strip()
        if not fname:
            print("ERROR: No file path provided")
            return
        
        if not os.path.exists(fname):
            print(f"ERROR: File not found: {fname}")
            return
        
        try:
            self.load_molecule(fname)
            self.update_plot()
            print(f"Successfully loaded {fname}")
        except Exception as e:
            print(f"ERROR loading file: {e}")
    
    def on_reset(self, event):
        """Reset rotations to zero"""
        self.slider_x.set_val(0)
        self.slider_y.set_val(0)
        self.slider_z.set_val(0)
        self.slider_pivot_x.set_val(0)
        self.slider_pivot_y.set_val(0)
        self.slider_pivot_z.set_val(0)
    
    def on_duplicate(self, event, n_copies):
        """Duplicate molecule n times rotationally"""
        self.duplicate_rotational(n_copies, axis='z')
        self.update_plot()
    
    def on_save(self, event):
        """Save all molecules to XYZ file"""
        fname = self.text_file.text.strip()
        if not fname:
            fname = "output.xyz"
        
        # Change extension to .xyz if needed
        if not fname.endswith('.xyz'):
            fname = fname.rsplit('.', 1)[0] + '_out.xyz'
        
        # Combine all molecules
        all_apos = []
        all_enames = []
        
        for apos, enames in self.molecules:
            all_apos.append(apos)
            all_enames.extend(enames)
        
        if len(all_apos) > 0:
            all_apos = np.vstack(all_apos)
            
            # Save using atomicUtils
            au.saveXYZ(all_enames, all_apos, fname, comment=f"Generated by MoleculeViewer, lattice_a={self.lattice_a}")
            print(f"Saved {len(all_apos)} atoms to {fname}")
        else:
            print("ERROR: No molecules to save")
    
    def on_clear(self, event):
        """Clear all copies, keep only original"""
        if len(self.molecules) > 0:
            self.molecules = [self.molecules[0]]
            self.colors = ['blue']
            self.current_mol_idx = 0
            self.update_plot()
    
    def update_plot(self):
        """Update 3D visualization"""
        self.ax.clear()
        
        if len(self.molecules) == 0:
            self.ax.set_xlabel('X (Å)')
            self.ax.set_ylabel('Y (Å)')
            self.ax.set_zlabel('Z (Å)')
            self.ax.set_title('No molecule loaded')
            plt.draw()
            return
        
        # Plot all molecules
        for i, (apos, enames) in enumerate(self.molecules):
            color = self.colors[i] if i < len(self.colors) else 'gray'
            
            # Plot atoms
            self.ax.scatter(apos[:, 0], apos[:, 1], apos[:, 2], 
                          c=color, s=20, alpha=0.6, label=f'Mol {i+1}')
        
        # Plot pivot point
        self.ax.scatter([self.pivot[0]], [self.pivot[1]], [self.pivot[2]], 
                       c='black', s=100, marker='x', linewidths=3, label='Pivot')
        
        # Set labels and title
        self.ax.set_xlabel('X (Å)')
        self.ax.set_ylabel('Y (Å)')
        self.ax.set_zlabel('Z (Å)')
        self.ax.set_title(f'Molecule Viewer (Lattice a={self.lattice_a} Å)')
        
        # Equal aspect ratio
        all_apos = np.vstack([apos for apos, _ in self.molecules])
        max_range = np.array([all_apos[:, 0].max() - all_apos[:, 0].min(),
                             all_apos[:, 1].max() - all_apos[:, 1].min(),
                             all_apos[:, 2].max() - all_apos[:, 2].min()]).max() / 2.0
        
        mid_x = (all_apos[:, 0].max() + all_apos[:, 0].min()) * 0.5
        mid_y = (all_apos[:, 1].max() + all_apos[:, 1].min()) * 0.5
        mid_z = (all_apos[:, 2].max() + all_apos[:, 2].min()) * 0.5
        
        self.ax.set_xlim(mid_x - max_range, mid_x + max_range)
        self.ax.set_ylim(mid_y - max_range, mid_y + max_range)
        self.ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        if len(self.molecules) > 1:
            self.ax.legend(loc='upper right', fontsize=8)
        
        plt.draw()
    
    def show(self):
        """Show the GUI"""
        plt.show()


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Interactive Molecular Viewer')
    parser.add_argument('file', nargs='?', help='XYZ file to load')
    parser.add_argument('--lattice', '-a', type=float, default=32.7, 
                       help='Lattice constant in Angstroms (default: 32.7)')
    
    args = parser.parse_args()
    
    viewer = MoleculeViewer(fname=args.file, lattice_a=args.lattice)
    viewer.show()


if __name__ == '__main__':
    main()
