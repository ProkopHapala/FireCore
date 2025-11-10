#!/usr/bin/env python3
"""
Interactive GUI for UFF scanning.
Allows real-time adjustment of scan parameters and visualization of results.
"""
import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from mpl_toolkits.mplot3d import Axes3D

# Add FireCore to the Python path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(base_path)

from pyBall import MMFF_multi as uff
from functions import *

# Define the path to common resource files
data_dir = os.path.join(base_path, "cpp/common_resources")

class ScanGUI:
    def __init__(self, args):
        self.args = args
        self.data_dir = data_dir
        
        # Initialize UFF library
        self.init_uff()
        
        # Get atom information
        self.natoms = uff.natoms
        self.apos = uff.apos.copy()
        self.atom_types = [uff.getTypeName(i, fromFF=True) for i in range(self.natoms)]
        # REQs columns: [R=radius, E=epsilon, Q=charge, H=?]
        self.atom_charges = uff.REQs[:, 2].copy()  # Q is the 3rd column (index 2)
        
        # Scan parameters
        self.scan_atom = args.scan_atom
        self.z_start, self.z_end, self.z_step = 0.0, 10.0, 0.02  # Very fine Z-scan (500 points)
        self.xy_start, self.xy_end, self.xy_step = 0.0, 10.0, 0.02
        self.scan_x, self.scan_y = 0.0, 0.0  # for z-scan
        self.scan_z = 2.0  # for xy-scan
        
        # Component flags
        self.component_flags = build_component_flags(
            base_components=None if args.comps is None else list(args.comps),
            enable=args.enable,
            disable=args.disable,
            preset=args.preset,
        )
        self.bNonBonded = args.non_bonded
        self.bGridFF = args.grid_ff
        
        # Create GUI
        self.create_gui()
        
        # Initial scan
        self.update_scans()
    
    def init_uff(self):
        """Initialize the UFF library."""
        print("--- Initializing MMFF_multi library ---")
        surf_name = os.path.join(self.data_dir, self.args.surf) if self.args.surf else None
        
        uff.init(
            nSys_=self.args.nsys,
            xyz_name=self.args.molecule,
            surf_name=surf_name,
            sElementTypes=os.path.join(self.data_dir, "ElementTypes.dat"),
            sAtomTypes=os.path.join(self.data_dir, "AtomTypes.dat"),
            sBondTypes=os.path.join(self.data_dir, "BondTypes.dat"),
            sAngleTypes=os.path.join(self.data_dir, "AngleTypes.dat"),
            sDihedralTypes=os.path.join(self.data_dir, "DihedralTypes.dat"),
            bMMFF=True,
            bUFF=True,
            nExplore=1000,
            nRelax=500,
            T=300,
            gamma=0.1
        )
        uff.getBuffs_UFF()
    
    def create_gui(self):
        """Create the GUI layout."""
        self.fig = plt.figure(figsize=(16, 10))
        
        # Create subplots
        gs = self.fig.add_gridspec(3, 3, left=0.05, right=0.95, top=0.95, bottom=0.35, hspace=0.3, wspace=0.3)
        self.ax_z = self.fig.add_subplot(gs[0, 0])
        self.ax_xy = self.fig.add_subplot(gs[0, 1])
        self.ax_3d = self.fig.add_subplot(gs[0, 2], projection='3d')
        
        # Info panel
        self.ax_info = self.fig.add_subplot(gs[1, :])
        self.ax_info.axis('off')
        
        # Sliders area (4 columns for all sliders)
        # Create two separate areas for Z-scan and XY-scan controls
        # Z-scan controls (left box)
        z_control_area = self.fig.add_gridspec(4, 10, left=0.05, right=0.48, top=0.2, bottom=0.05, hspace=0.3, wspace=0.02)
        # XY-scan controls (right box)
        xy_control_area = self.fig.add_gridspec(4, 10, left=0.52, right=0.95, top=0.2, bottom=0.05, hspace=0.3, wspace=0.02)

        
        # === Z-SCAN CONTROLS ===
        # Z start
        ax_z_start = self.fig.add_subplot(z_control_area[0, 0:8])
        ax_z_start_txt = self.fig.add_subplot(z_control_area[0, 8:10])
        # Z end
        ax_z_end = self.fig.add_subplot(z_control_area[1, 0:8])
        ax_z_end_txt = self.fig.add_subplot(z_control_area[1, 8:10])
        # Z-scan X position
        ax_scan_x = self.fig.add_subplot(z_control_area[2, 0:8])
        ax_scan_x_txt = self.fig.add_subplot(z_control_area[2, 8:10])
        # Z-scan Y position
        ax_scan_y = self.fig.add_subplot(z_control_area[3, 0:8])
        ax_scan_y_txt = self.fig.add_subplot(z_control_area[3, 8:10])

        # === XY-SCAN CONTROLS ===
        # XY start
        ax_xy_start = self.fig.add_subplot(xy_control_area[0, 0:8])
        ax_xy_start_txt = self.fig.add_subplot(xy_control_area[0, 8:10])
        # XY end
        ax_xy_end = self.fig.add_subplot(xy_control_area[1, 0:8])
        ax_xy_end_txt = self.fig.add_subplot(xy_control_area[1, 8:10])
        # XY-scan Z height
        ax_scan_z = self.fig.add_subplot(xy_control_area[2, 0:8])
        ax_scan_z_txt = self.fig.add_subplot(xy_control_area[2, 8:10])

        # Create sliders (without labels - we'll add them separately inside the boxes)
        self.slider_z_start = Slider(ax_z_start, '', -5.0, 15.0, valinit=self.z_start, valstep=0.1)
        self.slider_z_end = Slider(ax_z_end, '', -5.0, 15.0, valinit=self.z_end, valstep=0.1)
        self.slider_scan_x = Slider(ax_scan_x, '', -10.0, 20.0, valinit=self.scan_x, valstep=0.5)
        self.slider_scan_y = Slider(ax_scan_y, '', -10.0, 20.0, valinit=self.scan_y, valstep=0.5)
        self.slider_xy_start = Slider(ax_xy_start, '', -10.0, 20.0, valinit=self.xy_start, valstep=0.5)
        self.slider_xy_end = Slider(ax_xy_end, '', -10.0, 20.0, valinit=self.xy_end, valstep=0.5)
        self.slider_scan_z = Slider(ax_scan_z, '', -5.0, 15.0, valinit=self.scan_z, valstep=0.1)

        # Add slider labels inside the boxes
        ax_z_start.text(0.02, 1.3, 'Z start [Å]', transform=ax_z_start.transAxes, fontsize=9, va='top')
        ax_z_end.text(0.02, 1.3, 'Z end [Å]', transform=ax_z_end.transAxes, fontsize=9, va='top')
        ax_scan_x.text(0.02, 1.3, 'Z-scan X [Å]', transform=ax_scan_x.transAxes, fontsize=9, va='top')
        ax_scan_y.text(0.02, 1.3, 'Z-scan Y [Å]', transform=ax_scan_y.transAxes, fontsize=9, va='top')
        ax_xy_start.text(0.02, 1.3, 'XY start [Å]', transform=ax_xy_start.transAxes, fontsize=9, va='top')
        ax_xy_end.text(0.02, 1.3, 'XY end [Å]', transform=ax_xy_end.transAxes, fontsize=9, va='top')
        ax_scan_z.text(0.02, 1.3, 'XY-scan Z [Å]', transform=ax_scan_z.transAxes, fontsize=9, va='top')

        # Create textboxes for direct value input
        self.textbox_z_start = TextBox(ax_z_start_txt, '', initial=str(self.z_start))
        self.textbox_z_end = TextBox(ax_z_end_txt, '', initial=str(self.z_end))
        self.textbox_scan_x = TextBox(ax_scan_x_txt, '', initial=str(self.scan_x))
        self.textbox_scan_y = TextBox(ax_scan_y_txt, '', initial=str(self.scan_y))
        self.textbox_xy_start = TextBox(ax_xy_start_txt, '', initial=str(self.xy_start))
        self.textbox_xy_end = TextBox(ax_xy_end_txt, '', initial=str(self.xy_end))
        self.textbox_scan_z = TextBox(ax_scan_z_txt, '', initial=str(self.scan_z))

        # Add visual boxes around control areas
        import matplotlib.patches as mpatches
        # Z-scan box
        z_box = mpatches.FancyBboxPatch((0.04, 0.04), 0.45, 0.22, 
                                         boxstyle="round,pad=0.01", 
                                         edgecolor="blue", facecolor="none", 
                                         linewidth=2, transform=self.fig.transFigure, zorder=0)
        self.fig.patches.append(z_box)
        # XY-scan box
        xy_box = mpatches.FancyBboxPatch((0.51, 0.04), 0.45, 0.22, 
                                          boxstyle="round,pad=0.01", 
                                          edgecolor="green", facecolor="none", 
                                          linewidth=2, transform=self.fig.transFigure, zorder=0)
        self.fig.patches.append(xy_box)
        # Add labels
        self.fig.text(0.06, 0.24, "Z-SCAN CONTROLS", fontsize=10, weight="bold", color="blue")
        self.fig.text(0.53, 0.24, "XY-SCAN CONTROLS", fontsize=10, weight="bold", color="green")
        
        # Connect sliders to update function
        self.slider_z_start.on_changed(self.on_slider_change)
        self.slider_z_end.on_changed(self.on_slider_change)
        self.slider_scan_x.on_changed(self.on_slider_change)
        self.slider_scan_y.on_changed(self.on_slider_change)
        self.slider_xy_start.on_changed(self.on_slider_change)
        self.slider_xy_end.on_changed(self.on_slider_change)
        self.slider_scan_z.on_changed(self.on_slider_change)

        # Connect textboxes to update function
        self.textbox_z_start.on_submit(lambda text: self.on_textbox_change('z_start', text))
        self.textbox_z_end.on_submit(lambda text: self.on_textbox_change('z_end', text))
        self.textbox_scan_x.on_submit(lambda text: self.on_textbox_change('scan_x', text))
        self.textbox_scan_y.on_submit(lambda text: self.on_textbox_change('scan_y', text))
        self.textbox_xy_start.on_submit(lambda text: self.on_textbox_change('xy_start', text))
        self.textbox_xy_end.on_submit(lambda text: self.on_textbox_change('xy_end', text))
        self.textbox_scan_z.on_submit(lambda text: self.on_textbox_change('scan_z', text))
        
        # Atom selector buttons (simple prev/next for now)
        ax_prev = plt.axes([0.35, 0.28, 0.08, 0.03])
        ax_next = plt.axes([0.45, 0.28, 0.08, 0.03])
        self.btn_prev = Button(ax_prev, 'Prev Atom')
        self.btn_next = Button(ax_next, 'Next Atom')
        self.btn_prev.on_clicked(self.prev_atom)
        self.btn_next.on_clicked(self.next_atom)
        
        plt.show(block=False)
    
    def on_slider_change(self, val):
        """Called when any slider changes."""
        """Called when any slider changes."""
        self.z_start = self.slider_z_start.val
        self.z_end = self.slider_z_end.val
        self.scan_x = self.slider_scan_x.val
        self.scan_y = self.slider_scan_y.val
        self.xy_start = self.slider_xy_start.val
        self.xy_end = self.slider_xy_end.val
        self.scan_z = self.slider_scan_z.val
        self.update_scans()
    

    def on_textbox_change(self, param_name, text):
        """Handle textbox value changes."""
        try:
            value = float(text)
            
            # Update parameter and corresponding slider
            if param_name == 'z_start':
                self.z_start = value
                self.slider_z_start.set_val(value)
            elif param_name == 'z_end':
                self.z_end = value
                self.slider_z_end.set_val(value)
            elif param_name == 'scan_x':
                self.scan_x = value
                self.slider_scan_x.set_val(value)
            elif param_name == 'scan_y':
                self.scan_y = value
                self.slider_scan_y.set_val(value)
            elif param_name == 'xy_start':
                self.xy_start = value
                self.slider_xy_start.set_val(value)
            elif param_name == 'xy_end':
                self.xy_end = value
                self.slider_xy_end.set_val(value)
            elif param_name == 'scan_z':
                self.scan_z = value
                self.slider_scan_z.set_val(value)
            
            # Update plots (slider will trigger on_slider_change)
        except ValueError:
            print(f"Invalid value for {param_name}: {text}")

    def prev_atom(self, event):
        """Select previous atom."""
        self.scan_atom = (self.scan_atom - 1) % self.natoms
        self.update_scans()
    
    def next_atom(self, event):
        """Select next atom."""
        self.scan_atom = (self.scan_atom + 1) % self.natoms
        self.update_scans()
    
    def update_scans(self):
        """Perform scans and update all plots."""
        print(f"Updating scans for atom {self.scan_atom} ({self.atom_types[self.scan_atom]})...")
        
        # Perform Z-scan
        z_range_str = f"{self.z_start},{self.z_end},{self.z_step}"
        scan_pos_str = f"{self.scan_x},{self.scan_y}"
        confs_z, _, _, nconf_z = generate_confs(
            self.apos, 1000, 123, self.scan_atom, 'z',
            z_range_str, "", scan_pos_str
        )
        F_z = self.scan_uff(confs_z)
        
        # Perform XY-scan
        xy_range_str = f"{self.xy_start},{self.xy_end},{self.xy_step},{self.xy_start},{self.xy_end},{self.xy_step}"
        scan_pos_str = f"{self.scan_z}"
        confs_xy, _, _, nconf_xy = generate_confs(
            self.apos, 1000, 123, self.scan_atom, 'xy',
            "", xy_range_str, scan_pos_str
        )
        F_xy = self.scan_uff(confs_xy)
        
        # Update plots
        self.plot_z_scan(F_z)
        self.plot_xy_scan(F_xy)
        self.plot_3d_view()
        self.update_info()
        
        self.fig.canvas.draw_idle()
    
    def scan_uff(self, confs):
        """Run UFF scan on configurations."""
        components_to_switches(self.component_flags, bNonBonded=self.bNonBonded, bGridFF=self.bGridFF)
        F_gpu = uff.scan(confs, iParalel=2)
        return F_gpu
    
    def plot_z_scan(self, F):
        """Plot Z-scan results."""
        self.ax_z.clear()
        Z = np.arange(self.z_start, self.z_end, self.z_step)
        Fz = F[:len(Z), self.scan_atom, 2]
        self.ax_z.plot(Z, Fz, '-', linewidth=1.5, markersize=2)  # Smaller markers
        self.ax_z.set_xlabel('Z [Å]', fontsize=10)
        self.ax_z.set_ylabel('Fz [eV/Å]', fontsize=10)
        self.ax_z.set_title(f'Z-scan at X={self.scan_x:.1f}, Y={self.scan_y:.1f}', fontsize=11)
        self.ax_z.grid(True, alpha=0.3)
    
    def plot_xy_scan(self, F):
        """Plot XY-scan results."""
        self.ax_xy.clear()
        x_values = np.arange(self.xy_start, self.xy_end, self.xy_step)
        y_values = np.arange(self.xy_start, self.xy_end, self.xy_step)
        if len(x_values) == 0 or len(y_values) == 0:
            return
        X, Y = np.meshgrid(x_values, y_values)
        Fz = F[:len(x_values)*len(y_values), self.scan_atom, 2]
        Fz_grid = Fz.reshape(len(y_values), len(x_values))
        im = self.ax_xy.contourf(X, Y, Fz_grid, levels=20, cmap='RdBu_r')
        self.ax_xy.set_xlabel('X [Å]', fontsize=10)
        self.ax_xy.set_ylabel('Y [Å]', fontsize=10)
        self.ax_xy.set_title(f'XY-scan at Z={self.scan_z:.1f}', fontsize=11)
        self.ax_xy.set_aspect('equal', adjustable='box')  # Force square aspect ratio
        if not hasattr(self, 'cbar_xy'):
            self.cbar_xy = self.fig.colorbar(im, ax=self.ax_xy, label='Fz [eV/Å]')
        else:
            self.cbar_xy.update_normal(im)
    
    def plot_3d_view(self):
        """Plot 3D view of molecule and scan atom."""
        self.ax_3d.clear()
        # Plot all atoms with size proportional to charge magnitude
        colors = ['red' if i == self.scan_atom else 'blue' for i in range(self.natoms)]
        # Size proportional to |charge|, with minimum size
        charge_sizes = np.abs(self.atom_charges) * 100 + 20
        sizes = [200 if i == self.scan_atom else charge_sizes[i] for i in range(self.natoms)]
        self.ax_3d.scatter(self.apos[:, 0], self.apos[:, 1], self.apos[:, 2], c=colors, s=sizes, alpha=0.6)

        # Highlight scan atom
        scan_pos = self.apos[self.scan_atom]
        self.ax_3d.scatter([scan_pos[0]], [scan_pos[1]], [scan_pos[2]], c='red', s=250, marker='*', edgecolors='black', linewidths=2)

        # Draw scan lines
        # Z-scan line
        z_line = np.linspace(self.z_start, self.z_end, 20)
        x_line = np.full_like(z_line, self.scan_x)
        y_line = np.full_like(z_line, self.scan_y)
        self.ax_3d.plot(x_line, y_line, z_line, 'g--', alpha=0.5, linewidth=2, label='Z-scan line')

        # Draw XY-scan plane at scan_z
        xy_range = max(self.xy_end - self.xy_start, 5.0)
        x_plane = np.linspace(self.xy_start, self.xy_end, 10)
        y_plane = np.linspace(self.xy_start, self.xy_end, 10)
        X_plane, Y_plane = np.meshgrid(x_plane, y_plane)
        Z_plane = np.full_like(X_plane, self.scan_z)
        self.ax_3d.plot_surface(X_plane, Y_plane, Z_plane, alpha=0.2, color='orange', label='XY-scan plane')

        self.ax_3d.set_xlabel('X [Å]', fontsize=9)
        self.ax_3d.set_ylabel('Y [Å]', fontsize=9)
        self.ax_3d.set_zlabel('Z [Å]', fontsize=9)
        self.ax_3d.set_title('Molecule & Scan Geometry\n(atom size ∝ |charge|)', fontsize=10)
        self.ax_3d.legend(fontsize=8)
    
    def update_info(self):
        """Update info panel."""
        self.ax_info.clear()
        self.ax_info.axis('off')
        info_text = f"Selected Atom: {self.scan_atom} | Type: {self.atom_types[self.scan_atom]} | Charge: {self.atom_charges[self.scan_atom]:.3f} e\n"
        info_text += f"Position: ({self.apos[self.scan_atom, 0]:.2f}, {self.apos[self.scan_atom, 1]:.2f}, {self.apos[self.scan_atom, 2]:.2f}) Å\n"
        info_text += f"Total atoms: {self.natoms} | Surface: {os.path.basename(self.args.surf) if self.args.surf else 'None'}"
        self.ax_info.text(0.5, 0.5, info_text, ha='center', va='center', fontsize=12, family='monospace',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))


def components_to_switches(components, *, bNonBonded=False, bGridFF=False):
    """Set UFF component switches."""
    def resolve(name: str) -> int:
        val = components.get(name, -1)
        if val is None:
            return -1
        return int(val)
    
    DoBond = resolve('bonds')
    DoAngle = resolve('angles')
    DoDihedral = resolve('dihedrals')
    DoInversion = resolve('inversions')
    DoAssemble = 1 if any(v > 0 for v in (DoBond, DoAngle, DoDihedral, DoInversion)) else -1
    SubtractBondNonBond = 1
    ClampNonBonded = 1
    uff.setSwitches2(
        NonBonded=1 if bNonBonded else -1,
        SurfAtoms=1 if bGridFF else -1,
        GridFF=1 if bGridFF else -1,
    )
    uff.setSwitchesUFF(
        DoBond=1 if DoBond > 0 else -1,
        DoAngle=DoAngle,
        DoDihedral=DoDihedral,
        DoInversion=DoInversion,
        DoAssemble=DoAssemble,
        SubtractBondNonBond=SubtractBondNonBond,
        ClampNonBonded=ClampNonBonded
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Interactive UFF Scanning GUI')
    default_mol = os.path.join(data_dir, 'xyz', 'xylitol_WO_gridFF.xyz')
    
    parser.add_argument('-m', '--molecule', type=str, default=default_mol, help='Molecule file (.mol2, .xyz)')
    parser.add_argument('--nsys', type=int, default=100, help='Number of GPU replicas (systems)')
    parser.add_argument('--scan_atom', type=int, default=5, help='Initial atom to scan')
    parser.add_argument('--preset', choices=['all', 'bonded', 'bonded-only', 'none', 'grid-only', 'nonbonded-only'], default='grid-only', help='Component preset')
    parser.add_argument('--comps', nargs='*', choices=UFF_COMPONENTS, default=None, help='Explicit list of bonded UFF components')
    parser.add_argument('--enable', nargs='*', choices=UFF_COMPONENTS, default=None, help='Additional UFF components to enable')
    parser.add_argument('--disable', nargs='*', choices=UFF_COMPONENTS, default=None, help='UFF components to disable')
    parser.add_argument('--non-bonded', action='store_true', help='Enable non-bonded interactions')
    parser.add_argument('--grid-ff', action='store_true', default=True, help='Enable GridFF interactions')
    parser.add_argument('--surf', type=str, required=True, help='Surface file (.xyz)')
    
    args = parser.parse_args()
    
    # Cleanup old files
    cleanup_xyz_files(["scan_relaxed_cpu_*.xyz", "scan_relaxed_gpu_*.xyz", "trj_multi.xyz", "relaxed_gpu_*.xyz"])
    
    # Create and run GUI
    gui = ScanGUI(args)
    plt.show()
