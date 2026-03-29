import numpy as np
import os
import json
import time
import matplotlib.pyplot as plt
import sys

# Add path to import plotUtils
sys.path.append('../../')
from pyBall import plotUtils
from pyBall import atomicUtils

# ============================================================================
# SHARED PLOTTING UTILITIES
# ============================================================================

def load_element_types(molecule_xyz="../../cpp/common_resources/xyz/PTCDA.xyz"):
    """Load element types from XYZ file and convert to atomic numbers"""
    element_types = []
    
    print(f"🔍 LOADING ELEMENT TYPES FROM: {molecule_xyz}")
    
    if not os.path.exists(molecule_xyz):
        raise FileNotFoundError(f"XYZ file not found: {molecule_xyz}")
    
    with open(molecule_xyz, 'r') as f:
        lines = f.readlines()
    
    natoms = int(lines[0].strip())
    print(f"📊 EXPECTING {natoms} ATOMS")
    
    # Parse element symbols
    for i, line in enumerate(lines[2:2+natoms]):  # Skip header and comment line
        parts = line.strip().split()
        if len(parts) < 1:
            raise ValueError(f"Line {i+2} in XYZ file is empty or malformed: '{line}'")
        
        element_symbol = parts[0]
        print(f"  ATOM {i}: Symbol='{element_symbol}'")
        
        # Convert element symbol to atomic number - NO FALLBACKS!
        if element_symbol == 'H':
            element_num = 1
        elif element_symbol == 'C':
            element_num = 6
        elif element_symbol == 'N':
            element_num = 7
        elif element_symbol == 'O':
            element_num = 8
        elif element_symbol == 'S':
            element_num = 16
        else:
            raise ValueError(f"UNKNOWN ELEMENT SYMBOL: '{element_symbol}' in atom {i}")
        
        element_types.append(element_num)
    
    print(f"✅ LOADED {len(element_types)} ELEMENT TYPES: {element_types}")
    return element_types

def load_substrate(substrate_shift=None):
    """Load substrate atoms and apply shift if provided"""
    
    if substrate_shift is None:
        print("⚠️ NO SUBSTRATE SHIFT PROVIDED - SKIPPING SUBSTRATE LOADING")
        return None, None, None
    
    print(f"🔍 LOADING SUBSTRATE WITH SHIFT: {substrate_shift}")
    
    # Load substrate from common path
    substrate_path = "../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz"
    
    if not os.path.exists(substrate_path):
        raise FileNotFoundError(f"Substrate file not found: {substrate_path}")
    
    print(f"📂 LOADING SUBSTRATE FROM: {substrate_path}")
    
    with open(substrate_path, 'r') as f:
        lines = f.readlines()
    
    # Parse XYZ format
    natoms = int(lines[0].strip())
    print(f"📊 EXPECTING {natoms} SUBSTRATE ATOMS")
    
    substrate_pos = []
    substrate_types = []
    
    for i, line in enumerate(lines[2:2+natoms]):  # Skip header and comment line
        parts = line.strip().split()
        if len(parts) < 4:
            raise ValueError(f"Substrate line {i+2} is malformed: '{line}'")
        
        atom_type = parts[0]  # Ca or F
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        substrate_pos.append([x, y, z])
        substrate_types.append(atom_type)
        
        if i < 5:  # Print first few atoms for debugging
            print(f"  SUBSTRATE ATOM {i}: {atom_type} ({x:.3f}, {y:.3f}, {z:.3f})")
    
    substrate_pos = np.array(substrate_pos)
    print(f"✅ LOADED {len(substrate_pos)} SUBSTRATE ATOMS")
    
    # Count atom types
    type_counts = {}
    for atom_type in substrate_types:
        type_counts[atom_type] = type_counts.get(atom_type, 0) + 1
    print(f"📊 SUBSTRATE COMPOSITION: {type_counts}")
    
    # Apply shift
    substrate_shifted = substrate_pos.copy()
    substrate_shifted[:, 0] += substrate_shift[0]  # x shift
    substrate_shifted[:, 1] += substrate_shift[1]  # y shift
    substrate_shifted[:, 2] += substrate_shift[2]  # z shift
    
    print(f"🔄 APPLIED SHIFT: {substrate_shift}")
    
    return substrate_pos, substrate_types, substrate_shifted

def create_molecule_bonds(positions, element_types):
    """Create molecule bonds using findAllBonds with proper element types"""
    natoms = len(positions)
    
    print(f"🔗 CREATING MOLECULE BONDS FOR {natoms} ATOMS")
    
    if element_types is None:
        raise ValueError("ELEMENT TYPES REQUIRED - NO FALLBACKS ALLOWED")
    
    if len(element_types) != natoms:
        raise ValueError(f"ELEMENT TYPES LENGTH MISMATCH: {len(element_types)} vs {natoms} positions")
    
    print(f"📊 ELEMENT TYPES FOR BONDS: {element_types}")
    
    atoms_for_bonds = np.zeros((natoms, 4))
    atoms_for_bonds[:, 0] = element_types  # Use proper element types
    atoms_for_bonds[:, 1:4] = positions  # Use positions for bond detection
    
    # Create bonds using atomicUtils with original working parameters
    # Use original Rcut=3.0, RvdwCut=0.6 but filter out H-H bonds only
    bonds_result = atomicUtils.findAllBonds(atoms_for_bonds, Rcut=3.0, RvdwCut=0.6)
    
    print(f"✅ findAllBonds RETURNED: {type(bonds_result)}")
    print(f"🔍 RETURNED {len(bonds_result)} ITEMS")
    
    # findAllBonds returns (bonds, bond_vectors) - we only want the bonds
    if isinstance(bonds_result, (list, tuple)) and len(bonds_result) >= 1:
        bonds = bonds_result[0]  # Take only the first item (atom pairs)
        print(f"✅ EXTRACTED BONDS: {len(bonds)} atom pairs")
    else:
        bonds = bonds_result  # Fallback if format is different
    
    print(f"🔍 BONDS TYPE: {type(bonds)}")
    print(f"🔍 FIRST BOND TYPE: {type(bonds[0]) if len(bonds) > 0 else 'None'}")
    
    # Print first few bonds for debugging
    for i, bond in enumerate(bonds[:5]):
        print(f"  BOND {i}: type={type(bond)}, value={bond}")
        if hasattr(bond, '__len__'):
            print(f"    Length: {len(bond)}")
            if len(bond) >= 2:
                print(f"    Atoms: {bond[0]} - {bond[1]}")
    
    # Convert to simple list of atom pairs and filter H-H bonds
    clean_bonds = []
    for bond in bonds:
        if isinstance(bond, (list, tuple, np.ndarray)) and len(bond) >= 2:
            atom1 = int(bond[0])
            atom2 = int(bond[1]) 
            # Validate atom indices
            if atom1 < natoms and atom2 < natoms and atom1 >= 0 and atom2 >= 0:
                # Filter out H-H bonds but keep all other bonds (including O-O, O-C, etc.)
                elem1 = element_types[atom1] if atom1 < len(element_types) else "UNKNOWN"
                elem2 = element_types[atom2] if atom2 < len(element_types) else "UNKNOWN"
                
                # Skip H-H bonds only
                if elem1 == 1 and elem2 == 1:
                    print(f"    ❌ FILTERING H-H BOND: {atom1} - {atom2}")
                    continue
                
                clean_bonds.append([atom1, atom2])
                elem1_name = {1: 'H', 6: 'C', 8: 'O'}.get(elem1, f'TYPE{elem1}')
                elem2_name = {1: 'H', 6: 'C', 8: 'O'}.get(elem2, f'TYPE{elem2}')
                print(f"    ✅ VALID BOND: {atom1}({elem1_name}) - {atom2}({elem2_name})")
            else:
                print(f"⚠️  INVALID BOND: {atom1} - {atom2} (natoms={natoms})")
        else:
            print(f"⚠️  MALFORMED BOND: {bond}")
    
    bonds = clean_bonds
    print(f"✅ FINAL BONDS: {len(bonds)} proper atom pairs")
    
    # Print element types for first few bonded atoms to verify
    print(f"🔍 CHECKING ALL {len(bonds)} BONDS FOR H-H BONDS:")
    h_bonds = []
    for i, (a1, a2) in enumerate(bonds):
        elem1 = element_types[a1] if a1 < len(element_types) else "UNKNOWN"
        elem2 = element_types[a2] if a2 < len(element_types) else "UNKNOWN"
        elem1_name = {1: 'H', 6: 'C', 8: 'O'}.get(elem1, f'TYPE{elem1}')
        elem2_name = {1: 'H', 6: 'C', 8: 'O'}.get(elem2, f'TYPE{elem2}')
        
        if elem1 == 1 or elem2 == 1:  # If either atom is hydrogen
            h_bonds.append((i, a1, elem1_name, a2, elem2_name))
        
        if i < 10:  # Show first 10 bonds
            print(f"  BOND {i}: {a1}({elem1_name}) - {a2}({elem2_name})")
    
    print(f"🚨 FOUND {len(h_bonds)} HYDROGEN BONDS:")
    for bond_info in h_bonds:
        i, a1, e1, a2, e2 = bond_info
        print(f"  H-BOND {i}: {a1}({e1}) - {a2}({e2})")
    
    return bonds

def plot_substrate_atoms(ax, substrate_shifted, substrate_types, substrate_shift=None):
    """Plot substrate atoms with proper filtering and coloring"""
    if substrate_shifted is None:
        return
    
    # Find top surface atoms
    top_z = np.max(substrate_shifted[:, 2])
    height_cutoff = top_z - 3.0
    visible_mask = substrate_shifted[:, 2] >= height_cutoff
    
    if np.any(visible_mask):
        visible_pos = substrate_shifted[visible_mask]
        visible_types = [substrate_types[i] for i in range(len(substrate_types)) if visible_mask[i]]
        
        # Separate Ca and F atoms
        ca_mask = np.array([t == 'Ca' for t in visible_types])
        f_mask = np.array([t == 'F' for t in visible_types])
        
        # Determine axes from subplot
        axes = getattr(ax, '_plot_axes', (0, 1))  # Default to XY
        
        # Plot with smaller sizes for background
        if np.any(ca_mask):
            if len(axes) == 2:
                ax.scatter(visible_pos[ca_mask, axes[0]], visible_pos[ca_mask, axes[1]],  c='magenta', s=10, alpha=0.3, label='Ca atoms')
        if np.any(f_mask):
            if len(axes) == 2:
                ax.scatter(visible_pos[f_mask, axes[0]], visible_pos[f_mask, axes[1]], c='green', s=5, alpha=0.3, label='F atoms')

def plot_molecule_bonds(ax, bonds, positions, axes=(0, 1), initial_color='gray', final_color='black', lw=2.0):
    """Plot molecule bonds for initial and final positions"""
    if bonds is not None:
        # Set axes for plotUtils
        ax._plot_axes = axes
        plotUtils.plotBonds(links=bonds, ps=positions[0], colors=initial_color, lws=lw, axes=axes)
        plotUtils.plotBonds(links=bonds, ps=positions[-1], colors=final_color, lws=lw, axes=axes)

def plot_trajectory(ax, APs, ia_anchor, ia_opposite, axes=(0, 1), color_anchor='red', color_opposite='blue', 
                   lw=2.0, alpha=0.9, markersize=4, label_anchor=None, label_opposite=None):
    """Plot trajectory for anchor and opposite atoms"""
    if len(axes) == 2:
        ax.plot(APs[:, ia_anchor, axes[0]], APs[:, ia_anchor, axes[1]], 'o-', color=color_anchor, linewidth=lw, markersize=markersize,   alpha=alpha, label=label_anchor or f'Anchor {ia_anchor}')
        ax.plot(APs[:, ia_opposite, axes[0]], APs[:, ia_opposite, axes[1]],  'o-', color=color_opposite, linewidth=lw, markersize=markersize,  alpha=alpha, label=label_opposite or f'Opposite {ia_opposite}')

def plot_all_improvements(ax, improvements, ia_anchor, ia_opposite, axes=(0, 1),  bold_best=True, color_anchor='red', color_opposite='blue'):
    """Plot all improvement trajectories with best one bolded"""
    for b in improvements:
        AP = b['APs']
        lw = 1.0
        alpha = 0.2
        if bold_best and (b is improvements[-1]):
            lw = 2.5
            alpha = 0.9
        
        if len(axes) == 2:
            ax.plot(AP[:, ia_anchor, axes[0]], AP[:, ia_anchor, axes[1]], 
                   color=color_anchor, lw=lw, alpha=alpha)
            ax.plot(AP[:, ia_opposite, axes[0]], AP[:, ia_opposite, axes[1]], 
                   color=color_opposite, lw=lw, alpha=alpha)

def plot_control_points(ax, control_pts, axes=(0, 1), color='red', markersize=10, 
                       markeredgewidth=2, add_numbers=True):
    """Plot spline control points with optional numbering"""
    if control_pts is not None and len(axes) == 2:
        control_array = np.array(control_pts)
        ax.plot(control_array[:, axes[0]], control_array[:, axes[1]], 
               'x', color=color, markersize=markersize, markeredgewidth=markeredgewidth, 
               label='Spline control points')
        
        if add_numbers:
            for i in range(len(control_array)):
                x, y = control_array[i, axes[0]], control_array[i, axes[1]]
                ax.text(x, y + 0.3, str(i), fontsize=8, ha='center', color=color, fontweight='bold')

def plot_target(ax, target_pos, axes=(0, 1), color='black', markersize=20, markeredgewidth=3):
    """Plot target position"""
    if len(axes) == 2:
        ax.plot(target_pos[axes[0]], target_pos[axes[1]], 
               'o', color=color, markersize=markersize, fillstyle='none', 
               markeredgewidth=markeredgewidth, label='Target')

def add_fitness_text(ax, fitness_text, transform=None):
    """Add fitness information text box"""
    if transform is None:
        transform = ax.transAxes
    ax.text(0.02, 0.98, fitness_text, transform=transform, verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

def add_penalty_decomposition_text(ax, E_total, l_pos, l_force, w_pos, w_force, fmax, f_safe, transform=None):
    """Add detailed penalty decomposition text box"""
    if transform is None:
        transform = ax.transAxes
    
    penalty_text = (
        f"PENALTY DECOMPOSITION:\n"
        f"Total: {E_total:.6f}\n"
        f"├─ Position: {w_pos:.3f} × {l_pos:.6f} = {w_pos * l_pos:.6f}\n"
        f"└─ Force:    {w_force:.3f} × {l_force:.6f} = {w_force * l_force:.6f}\n"
        f"Force details: fmax={fmax:.3f}, f_safe={f_safe:.3f}"
    )
    
    ax.text(0.02, 0.98, penalty_text, transform=transform, verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9), 
           fontsize=9, family='monospace')

# ============================================================================
# MAIN PLOTTING FUNCTIONS
# ============================================================================

def plot_single_improvement_comprehensive(APs, Fcon, target_pos, ia_anchor, ia_opposite, 
                                       attempt, out_dir=".", substrate_shift=None, 
                                       surface_name=None, plot_substrate=True, 
                                       plot_mol_atoms=False, plot_mol_bonds=True, 
                                       opt_target=None, current_control_pts=None, 
                                       suffix="", title_extra="", Es=None, 
                                       mmff_instance=None, bold_best=True,
                                       penalty_params=None):
    """Comprehensive 4-panel plot for single improvement"""
    
    # Fixed axis ranges for stable movie
    xlim = (-5, 20)
    ylim = (-5, 10)
    zlim = (5, 20)
    Flim = 3.0
    Elim = 3.0
    
    # Set fixed layout policy to prevent automatic adjustments
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.constrained_layout.use'] = False
    
    # Load shared resources
    element_types = load_element_types()
    substrate_pos, substrate_types, substrate_shifted = load_substrate(substrate_shift)
    bonds = create_molecule_bonds(APs[0], element_types)
    
    # Create figure
    fig = plt.figure(figsize=(16, 12))
    
    # Panel 1: Top view (XY)
    ax1 = fig.add_subplot(2, 2, 1)
    ax1._plot_axes = (0, 1)
    
    if plot_substrate:
        plot_substrate_atoms(ax1, substrate_shifted, substrate_types, substrate_shift)
    
    if plot_mol_bonds:
        plot_molecule_bonds(ax1, bonds, APs, axes=(0, 1), initial_color='gray', final_color='black')
    
    plot_trajectory(ax1, APs, ia_anchor, ia_opposite, axes=(0, 1))
    plot_control_points(ax1, current_control_pts, axes=(0, 1))
    plot_target(ax1, target_pos, axes=(0, 1))
    
    ax1.set_xlabel('x (Å)')
    ax1.set_ylabel('y (Å)')
    ax1.set_title(f'Improvement {attempt} - Top View{title_extra}')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Side view (XZ)
    ax2 = fig.add_subplot(2, 2, 2)
    ax2._plot_axes = (0, 2)
    
    if plot_substrate:
        plot_substrate_atoms(ax2, substrate_shifted, substrate_types, substrate_shift)
    
    if plot_mol_bonds:
        plot_molecule_bonds(ax2, bonds, APs, axes=(0, 2), initial_color='gray', final_color='black')
    
    plot_trajectory(ax2, APs, ia_anchor, ia_opposite, axes=(0, 2))
    plot_control_points(ax2, current_control_pts, axes=(0, 2))
    plot_target(ax2, target_pos, axes=(0, 2))
    
    ax2.set_xlabel('x (Å)')
    ax2.set_ylabel('z (Å)')
    ax2.set_title(f'Improvement {attempt} - Side View{title_extra}')
    ax2.set_xlim(xlim)
    ax2.set_ylim(zlim)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Energy along trajectory
    ax3 = fig.add_subplot(2, 2, 3)
    
    if Es is not None:
        ax3.plot(range(len(Es)), Es, '.-', linewidth=1.5, label='System Energy')
        ax3.set_ylabel('Energy (eV)')
        ax3.set_ylim(-Elim, Elim)  # Fixed energy range for stable frames
    else:
        raise ValueError("No energy data provided to plotting function!")
    
    ax3.set_title(f'Improvement {attempt} - Energy Along Trajectory')
    ax3.set_xlabel('Configuration')
    ax3.set_xlim(0, len(Es)-1)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Panel 4: Force components on anchor atom
    ax4 = fig.add_subplot(2, 2, 4)
    
    # Extract force components for anchor atom
    if Fcon.shape[1] == 3:  # Anchor forces only
        F_anchor = Fcon
    else:  # All atom forces - extract anchor
        n_atoms = Fcon.shape[1] // 3
        F_anchor = Fcon[:, ia_anchor*3:(ia_anchor+1)*3]
    
    F_mag = np.sqrt(np.sum(F_anchor**2, axis=1))
    
    ax4.plot(range(len(F_anchor)), F_anchor[:, 0], 'r-', label='Fx', linewidth=1.5)
    ax4.plot(range(len(F_anchor)), F_anchor[:, 1], 'g-', label='Fy', linewidth=1.5)
    ax4.plot(range(len(F_anchor)), F_anchor[:, 2], 'b-', label='Fz', linewidth=1.5)
    ax4.plot(range(len(F_anchor)), F_mag, '.-k', label='|F|', linewidth=2)
    
    ax4.set_xlabel('Configuration index')
    ax4.set_ylabel('Force (eV/Å)')
    ax4.set_title(f'Improvement {attempt} - Forces on Anchor {ia_anchor}')
    ax4.set_xlim(0, len(F_anchor[:, 0])-1)
    ax4.set_ylim(-Flim, Flim)  # Fixed y-axis for consistent frames
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    # Add penalty decomposition text to Panel 4 if penalty params provided
    if penalty_params is not None:
        add_penalty_decomposition_text(ax4, **penalty_params)
    
    # NO tight_layout - keep fixed layout for stable movies
    filename = os.path.join(out_dir, f'improvement_{attempt:05d}{suffix}.png')
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"Saved comprehensive improvement plot: {filename}")

def plot_improvements_summary(improvements, ia_anchor, ia_opposite, target_pos, 
                            substrate_shift=None, bold_best=True, out_dir=".", 
                            suffix=""):
    """Summary plot showing all improvements with substrate and molecule bonds"""
    
    # Fixed axis ranges for stable movie
    xlim = (-5, 20)
    ylim = (-5, 10)
    zlim = (5, 25)
    
    # Load shared resources
    element_types = load_element_types()
    substrate_pos, substrate_types, substrate_shifted = load_substrate(substrate_shift)
    
    if len(improvements) == 0:
        return
    
    # Create bonds from best improvement
    final_AP = improvements[-1]['APs']
    bonds = create_molecule_bonds(final_AP[0], element_types)
    
    # Create figure
    fig = plt.figure(figsize=(12, 5))
    
    # XY panel
    ax1 = fig.add_subplot(1, 2, 1)
    ax1._plot_axes = (0, 1)
    
    if substrate_shifted is not None:
        plot_substrate_atoms(ax1, substrate_shifted, substrate_types, substrate_shift)
    
    plot_all_improvements(ax1, improvements, ia_anchor, ia_opposite, axes=(0, 1), bold_best=bold_best)
    
    if bonds is not None:
        plot_molecule_bonds(ax1, bonds, final_AP, axes=(0, 1))
    
    plot_target(ax1, target_pos, axes=(0, 1))
    
    ax1.set_title('Optimization improvements (xy)')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.legend(loc='upper right', fontsize=8)
    
    # XZ panel
    ax2 = fig.add_subplot(1, 2, 2)
    ax2._plot_axes = (0, 2)
    
    if substrate_shifted is not None:
        plot_substrate_atoms(ax2, substrate_shifted, substrate_types, substrate_shift)
    
    plot_all_improvements(ax2, improvements, ia_anchor, ia_opposite, axes=(0, 2), bold_best=bold_best)
    
    if bonds is not None:
        plot_molecule_bonds(ax2, bonds, final_AP, axes=(0, 2))
    
    plot_target(ax2, target_pos, axes=(0, 2))
    
    ax2.set_title('Optimization improvements (xz)')
    ax2.set_xlim(xlim)
    ax2.set_ylim(zlim)
    
    # NO tight_layout - keep fixed layout for stable movies
    filename = os.path.join(out_dir, f'improvements{suffix}.png')
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"Saved improvements summary plot: {filename}")

# ============================================================================
# LEGACY FUNCTIONS (for backward compatibility)
# ============================================================================

def _plot_single_improvement(APs, Fcon, target_pos, ia_anchor, ia_opposite, attempt, out_dir, 
                           mmff_instance, substrate_shift, surface_name=None, 
                           plot_substrate=True, plot_mol_atoms=True, 
                           plot_mol_bonds=True, opt_target=None, 
                           current_control_pts=None, suffix="", title_extra="", Es=None,
                           penalty_params=None):
    """Legacy wrapper for backward compatibility"""
    plot_single_improvement_comprehensive(
        APs=APs, Fcon=Fcon, target_pos=target_pos, ia_anchor=ia_anchor, 
        ia_opposite=ia_opposite, attempt=attempt, out_dir=out_dir, 
        substrate_shift=substrate_shift, surface_name=surface_name, 
        plot_substrate=plot_substrate, plot_mol_atoms=plot_mol_atoms, 
        plot_mol_bonds=plot_mol_bonds, opt_target=opt_target, 
        current_control_pts=current_control_pts, suffix=suffix, 
        title_extra=title_extra, Es=Es, mmff_instance=mmff_instance,
        penalty_params=penalty_params)

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def load_spline_dat(fname):
    pts = []
    with open(fname, 'r') as f:
        for line in f:
            s = line.strip()
            if (not s) or s.startswith('#'):
                continue
            ws = s.split()
            if len(ws) < 3:
                continue
            pts.append([float(ws[0]), float(ws[1]), float(ws[2])])
    pts = np.array(pts, dtype=np.float64)
    print(f"load_spline_dat({fname}) npts={len(pts)}")
    if pts.shape[0] < 5:
        raise ValueError(f"Spline must have at least 5 points to freeze first2+last2; got {pts.shape[0]}")
    return pts


def save_spline_dat(fname, pts, header=None):
    with open(fname, 'w') as f:
        if header is not None:
            f.write(str(header).rstrip() + "\n")
        for p in pts:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
    print(f"save_spline_dat({fname}) npts={len(pts)}")


def _json_default(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.float32, np.float64)):
        return float(x)
    if isinstance(x, (np.int32, np.int64)):
        return int(x)
    return str(x)

# ============================================================================
# OPTIMIZER CLASS
# ============================================================================

class TipSplineSAOptimizer:
    def __init__(
        self,
        target_pos,
        ia_anchor,
        ia_opposite,
        w_pos=1.0,
        w_force=0.0,
        f_safe=5.0,
        temp0=0.0,
        cooling=1.0,
        step0=10.0,
        step_cooling=0.99,
        freeze_n_ends=2,
        seed=42,
        bbox_x=(-50, 50),
        bbox_y=(-50, 50),
        bbox_z=(10, 20)
    ):
        """
        Initialize optimizer with bounding box constraints
        
        Parameters:
        -----------
        bbox_x, bbox_y, bbox_z: tuples defining min/max bounds for control points
        """
        self.target_pos = np.array(target_pos, dtype=np.float64)
        self.ia_anchor = ia_anchor
        self.ia_opposite = ia_opposite
        self.w_pos = float(w_pos)
        self.w_force = float(w_force)
        self.f_safe = float(f_safe)
        self.temp0 = float(temp0)
        self.cooling = float(cooling)
        self.step0 = float(step0)
        self.step_cooling = float(step_cooling)
        self.freeze_n_ends = int(freeze_n_ends)
        self.rng = np.random.default_rng(int(seed))
        
        # Bounding box constraints
        self.bbox_x = bbox_x
        self.bbox_y = bbox_y  
        self.bbox_z = bbox_z

        if self.target_pos.shape != (3,):
            raise ValueError(f"target_pos must be (3,); got {self.target_pos.shape}")

    def _loss_components(self, APs, Fcon):
        rB = np.array(APs[-1, self.ia_opposite, :3], dtype=np.float64)
        d = rB - self.target_pos
        l_pos = float(np.dot(d, d))

        fc = np.array(Fcon[:, :3], dtype=np.float64)
        fmag = np.linalg.norm(fc, axis=1)
        fmax = float(np.max(fmag))
        over = max(0.0, fmax - self.f_safe)
        l_force = float(over * over)

        # Debug: Show force calculation details
        if self.w_force > 0.01:  # Only show if force weighting is significant
            print(f"🔍 FORCE PENALTY: fmax={fmax:.3f}, f_safe={self.f_safe:.3f}, over={over:.3f}")
            print(f"🔍 FORCE CONTRIBUTION: w_force={self.w_force:.3f} × l_force={l_force:.6f} = {self.w_force * l_force:.6f}")
            if over == 0.0:
                print(f"⚠️  FORCE PENALTY IS ZERO: fmax ({fmax:.3f}) < f_safe ({self.f_safe:.3f})")
                print(f"💡 Try reducing --opt-fsafe to {fmax:.1f} or lower to activate force penalty")

        l_tot = self.w_pos * l_pos + self.w_force * l_force
        return l_tot, l_pos, l_force, fmax, rB

    def _clamp_to_bbox(self, pts):
        """Clamp control points to bounding box"""
        pts_clamped = pts.copy()
        original = pts.copy()
        
        # Clamp x coordinates
        pts_clamped[:, 0] = np.clip(pts_clamped[:, 0], self.bbox_x[0], self.bbox_x[1])
        # Clamp y coordinates  
        pts_clamped[:, 1] = np.clip(pts_clamped[:, 1], self.bbox_y[0], self.bbox_y[1])
        # Clamp z coordinates
        pts_clamped[:, 2] = np.clip(pts_clamped[:, 2], self.bbox_z[0], self.bbox_z[1])
        
        # Debug: Check if any points were clamped
        if not np.array_equal(original, pts_clamped):
            print(f"🔳 BBOX CLAMPING: Points constrained to bounds:")
            print(f"   x: [{self.bbox_x[0]}, {self.bbox_x[1]}]")
            print(f"   y: [{self.bbox_y[0]}, {self.bbox_y[1]}]")  
            print(f"   z: [{self.bbox_z[0]}, {self.bbox_z[1]}]")
        
        return pts_clamped

    def _propose(self, pts, temp, current_step):
        pts2 = pts.copy()
        n = pts2.shape[0]
        i0 = self.freeze_n_ends
        i1 = n - self.freeze_n_ends
        if i1 <= i0:
            raise ValueError(f"Not enough internal points to mutate: n={n} freeze_n_ends={self.freeze_n_ends}")

        ip = int(self.rng.integers(i0, i1))
        # Step size independent of temperature - temperature only affects acceptance
        # Use current step size that decays over attempts
        step = current_step * 0.25  # Reduce by 4x from current step size
        print(f"🔍 DEBUG MUTATION: temp={temp}, current_step={current_step}, step={step}")
        # Use uniform distribution for more controlled mutations
        d = self.rng.uniform(-1.0, 1.0, size=3) * step
        pts2[ip, :] += d
        
        # Apply bounding box constraints
        pts2 = self._clamp_to_bbox(pts2)
        return pts2, ip, d

    def _accept(self, Eold, Enew, temp):
        # Metropolis acceptance criterion:
        # 1. Always accept better solutions
        # 2. Accept worse solutions with probability exp(-ΔE/temp)
        # 3. Temperature controls acceptance of worse solutions ONLY
        
        if Enew <= Eold:
            return True  # Better or equal to current - always accept
        
        # Worse solution: accept with probability based on temperature
        if temp <= 0:
            return False  # Zero temperature - never accept worse
        
        delta_E = Enew - Eold
        prob_accept = np.exp(-delta_E / temp)
        return self.rng.random() < prob_accept

    def optimize(
        self,
        mmff_instance,
        ts,
        spline_init_fname=None,
        out_dir=None,
        n_attempts=100,
        scan_kwargs=None,
        save_every_improvement=True,
        plot_improvements=True,
        plot_all_trials=False,
        bold_best=True,
        substrate_shift=None,
        surface_name=None,
        plot_substrate=True,
        plot_molecule_atoms=False,
        plot_molecule_bonds=True,
        opt_target=None,
        nconf=40,
        **kwargs
    ):
        """Run simulated annealing optimization with comprehensive plotting"""
        
        print(f"=== Starting TipSplineSAOptimizer.optimize() ===")
        print(f"n_attempts={n_attempts}, out_dir={out_dir}")
        print(f"plot_improvements={plot_improvements}, plot_all_trials={plot_all_trials}")
        print(f"substrate_shift={substrate_shift}")
        
        # Initialize scan_kwargs if not provided
        if scan_kwargs is None:
            scan_kwargs = {}
        
        # Load initial spline
        if spline_init_fname is not None and os.path.exists(spline_init_fname):
            pts = load_spline_dat(spline_init_fname)
            print(f"Loaded initial spline from {spline_init_fname}")
        else:
            # Create default initial spline
            pts = np.array([
                [0.0, 0.0, 14.0],
                [5.0, 0.0, 14.0],
                [10.0, 0.0, 14.0],
                [15.0, 0.0, 14.0],
                [20.0, 0.0, 14.0]
            ], dtype=np.float64)
            print(f"Created default initial spline")
        
        # Ensure output directory exists
        os.makedirs(out_dir, exist_ok=True)
        
        # CRITICAL: Save initial atom positions for consistent reset
        initial_atom_pos = mmff_instance.apos[:, :].copy()
        
        # Save initial spline
        spline_current = os.path.join(out_dir, "spline_current.dat")
        save_spline_dat(spline_current, pts)
        print(f"Saved initial spline to {spline_current}")
        
        # Evaluate initial configuration
        print("TipSplineSAOptimizer.optimize(): evaluating initial spline")
        
        # CRITICAL: Reset molecule to initial position for consistent initialization
        # This ensures attempt -1 uses the same starting state as all other attempts
        mmff_instance.apos[:, :] = initial_atom_pos
        
        Es, FAs, APs, Fcon = mmff_instance.scan_manipulation(ts, spline_fname=spline_current, **scan_kwargs)
        
        # Calculate initial fitness
        Ebest, lpos_best, lforce_best, fmax_best, rB_best = self._loss_components(APs, Fcon)
        
        # Initialize best configuration
        best = {
            'E': Ebest,
            'l_pos': lpos_best,
            'l_force': lforce_best,
            'fmax': fmax_best,
            'r_opposite_final': rB_best,
            'pts': pts.copy(),
            'APs': np.array(APs, copy=True),
            'Fcon': np.array(Fcon, copy=True),
            'attempt': -1,
        }
        
        # Always plot initial trajectory as baseline reference
        if plot_improvements:
            print("Plotting initial trajectory as baseline...")
            _plot_single_improvement(
                best['APs'], 
                best['Fcon'],
                self.target_pos, 
                self.ia_anchor, 
                self.ia_opposite, 
                -1,  # Use -1 for initial/baseline
                out_dir,
                mmff_instance,
                substrate_shift,
                surface_name,
                plot_substrate,
                plot_molecule_atoms,
                plot_molecule_bonds,
                opt_target,
                current_control_pts=pts,  # Add the actual control points!
                suffix="",
                title_extra="(Initial)",
                Es=Es,  # Add actual energy data!
                penalty_params={
                    'E_total': best['E'],
                    'l_pos': best.get('l_pos', 0.0),
                    'l_force': best.get('l_force', 0.0),
                    'w_pos': self.w_pos,
                    'w_force': self.w_force,
                    'fmax': best.get('fmax', 0.0),
                    'f_safe': self.f_safe
                }
            )

        improvements = []
        if save_every_improvement:
            improvements.append(best)

        # Optimization loop
        print(f"=== Starting optimization loop: {n_attempts} attempts ===")
        Ecur = Ebest
        temp = self.temp0
        current_step = self.step0
        
        for it in range(int(n_attempts)):
            pts_new, ip, d = self._propose(pts, temp, current_step)
            
            # CRITICAL: Use unique spline filename for each attempt
            # This prevents file caching/overwrite issues that could cause non-deterministic behavior
            spline_attempt = f"{out_dir}/spline_attempt_{it:05d}.dat"
            save_spline_dat(spline_attempt, pts_new)

            print(f"\n=== OPT attempt {it}/{n_attempts} temp={temp:.6g} mut_i={ip} d={d} ===")

            # CRITICAL: Reset molecule to initial position before each evaluation
            # This prevents accumulated drift and ensures deterministic evaluations
            mmff_instance.apos[:, :] = initial_atom_pos

            Es, FAs, APs, Fcon = mmff_instance.scan_manipulation(ts, spline_fname=spline_attempt, **scan_kwargs)
            
            # DEBUG: Check t=0 positions for both anchor and opposite atoms
            anchor_t0 = APs[0, self.ia_anchor, :]  # t=0, anchor atom
            opposite_t0 = APs[0, self.ia_opposite, :]  # t=0, opposite atom
            print(f"=== DEBUG ATTEMPT {it}: t=0 positions ===")
            print(f"    Anchor:   ({anchor_t0[0]:.6f}, {anchor_t0[1]:.6f}, {anchor_t0[2]:.6f})")
            print(f"    Opposite: ({opposite_t0[0]:.6f}, {opposite_t0[1]:.6f}, {opposite_t0[2]:.6f})")
            
            Enew, lpos, lforce, fmax, rB = self._loss_components(APs, Fcon)
            
            # Print fitness details every step
            print(f"FITNESS: E_total={Enew:.6f} (pos={lpos:.6f} force={lforce:.6f}) fmax={fmax:.6f}")
            print(f"CURRENT: E_cur={Ecur:.6f} BEST: E_best={best['E']:.6f}")
            
            # Show improvement/worsening clearly
            if Enew < best['E'] - 1e-6:
                print(f"✅ NEW BEST! ({best['E']:.6f} → {Enew:.6f})")
            elif Enew < Ecur - 1e-6:
                print(f"⬇️  Better than current ({Ecur:.6f} → {Enew:.6f}) - NOT BEST EVER")
            elif Enew > Ecur * 1.1:
                print(f"❌ MUCH WORSE than current ({Ecur:.6f} → {Enew:.6f})")
            else:
                print(f"➡️  Slightly worse ({Ecur:.6f} → {Enew:.6f})")

            accepted = self._accept(Ecur, Enew, temp)
            print(f"ACCEPTED: {accepted}")
            improved = False
            
            # Check for current improvement BEFORE updating Ecur
            is_current_improvement = (Enew < Ecur)
            
            if accepted:
                pts = pts_new
                Ecur = Enew

            if Enew < best['E'] - 1e-6:  
                best = {
                    'E': Enew,
                    'l_pos': lpos,
                    'l_force': lforce,
                    'fmax': fmax,
                    'r_opposite_final': rB,
                    'pts': pts_new.copy(),
                    'APs': np.array(APs, copy=True),
                    'Fcon': np.array(Fcon, copy=True),
                    'attempt': it,
                }
                # Update best energy ever
                self._best_E_ever = Enew
                improved = True
                print(f"NEW BEST: E={Enew:.6g} (pos={lpos:.6g} force={lforce:.6g} fmax={fmax:.6g})")

                if save_every_improvement:
                    improvements.append(best)
                    fn = os.path.join(out_dir, f"best_{it:05d}.dat")
                    save_spline_dat(fn, pts_new)
                    np.save(os.path.join(out_dir, f"best_{it:05d}_APs.npy"), best['APs'])
                    np.save(os.path.join(out_dir, f"best_{it:05d}_Fcon.npy"), best['Fcon'])
                    
                    # Plot individual improvement if requested
                    if plot_improvements:
                        _plot_single_improvement(
                            best['APs'], 
                            best['Fcon'],
                            self.target_pos, 
                            self.ia_anchor, 
                            self.ia_opposite, 
                            it, 
                            out_dir,
                            mmff_instance,
                            substrate_shift,
                            surface_name,
                            plot_substrate,
                            plot_molecule_atoms,
                            plot_molecule_bonds,
                            opt_target,
                            current_control_pts=pts_new,  # Add the actual control points!
                            suffix="",
                            title_extra=f"Improvement {it}: E={Enew:.6g}",
                            Es=Es,  # Add actual energy data!
                            penalty_params={
                                'E_total': Enew,
                                'l_pos': lpos,
                                'l_force': lforce,
                                'w_pos': self.w_pos,
                                'w_force': self.w_force,
                                'fmax': fmax,
                                'f_safe': self.f_safe
                            }
                        )

            # Plot current improvements (not best ever) using the saved flag
            # REMOVED: Only plot best-ever improvements for clarity
            
            if plot_all_trials:
                # Use the actual control points from this trial (pts_new)
                # NOT from file - use the optimizer variable directly!
                actual_control_pts = pts_new.copy()
                
                # DEBUG: Print what pts_new contains
                if len(pts_new) > 3:
                    pt3 = pts_new[3]
                    print(f"🔍 DEBUG TRIAL {it}: pts_new[3] = ({pt3[0]:.6f}, {pt3[1]:.6f}, {pt3[2]:.6f})")
                else:
                    print(f"🔍 DEBUG TRIAL {it}: pts_new too short: {len(pts_new)}")
                
                trial_data = {
                    'E': Enew,
                    'l_pos': lpos,
                    'l_force': lforce,
                    'fmax': fmax,
                    'r_opposite_final': rB,
                    'pts': actual_control_pts.copy(),
                    'APs': np.array(APs, copy=True),
                    'Fcon': np.array(Fcon, copy=True),
                    'attempt': it,
                    'accepted': accepted,
                    'improved': False  # Will be updated below if needed
                }
                
                _plot_single_improvement(
                    trial_data['APs'], 
                    trial_data['Fcon'],
                    self.target_pos, 
                    self.ia_anchor, 
                    self.ia_opposite, 
                    it, 
                    out_dir, 
                    mmff_instance, 
                    substrate_shift=substrate_shift,
                    surface_name=surface_name,
                    plot_substrate=plot_substrate,
                    plot_molecule_atoms=plot_molecule_atoms,
                    plot_molecule_bonds=plot_molecule_bonds,
                    opt_target=opt_target,
                    current_control_pts=actual_control_pts,  # Use the actual control points from file!
                    suffix=f"trial_{it:05d}",
                    title_extra=f"Trial {it}: E={Enew:.3f} {'ACC' if accepted else 'REJ'}",
                    Es=Es,  # Add actual energy data!
                    penalty_params={
                        'E_total': Enew,
                        'l_pos': lpos,
                        'l_force': lforce,
                        'w_pos': self.w_pos,
                        'w_force': self.w_force,
                        'fmax': fmax,
                        'f_safe': self.f_safe
                    }
                )

            temp *= self.cooling
            current_step *= self.step_cooling

        # Save final best spline
        best_spline = os.path.join(out_dir, "best_final.dat")
        save_spline_dat(best_spline, best['pts'])
        np.save(os.path.join(out_dir, "best_final_APs.npy"), best['APs'])
        np.save(os.path.join(out_dir, "best_final_Fcon.npy"), best['Fcon'])

        # Plot improvements summary using clean modular function
        if plot_improvements and (len(improvements) > 0):
            print(f"Plotting improvements: {len(improvements)}")
            plot_improvements_summary(
                improvements=improvements,
                ia_anchor=self.ia_anchor,
                ia_opposite=self.ia_opposite,
                target_pos=self.target_pos,
                substrate_shift=substrate_shift,
                bold_best=bold_best,
                out_dir=out_dir,
                suffix=""
            )
        
        return best_spline, best
