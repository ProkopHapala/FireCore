#!/usr/bin/env python3

import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import argparse

sys.path.append("../../")
from pyBall import MMFF as mmff

def parse_args():
    """Parse command line arguments for spline scan simulation."""
    parser = argparse.ArgumentParser(
        description='Spline-anchored molecular dynamics scan of PTCDA molecule',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Molecule and surface setup
    parser.add_argument('--xyz',        type=str, default='common_resources/xyz/PTCDA', help='XYZ file for molecule (without extension)')
    parser.add_argument('--surface',    type=str, default='common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top',help='Surface definition file (set to None for no substrate)')
    parser.add_argument('--no-surface', type=int, default=0, help='Disable surface interaction (equivalent to --surface None)')
    
    # Spline and constraint parameters
    parser.add_argument('--spline',  type=str, default='tipSpline_x.dat',help='Spline path definition file')
    parser.add_argument('--anchor',  type=int, default=28, help='Anchor atom index (0-based)')
    parser.add_argument('--kanchor', type=float, default=20.0,help='Spring constant for anchor constraint (eV/A²)')

    # Trajectory optimization parameters
    parser.add_argument('--optimize',       type=int,   default=0, help='Optimize spline control points (1=yes, 0=no)')
    parser.add_argument('--opt-outdir',     type=str,   default='opt_tipSpline', help='Output directory for optimizer logs and best splines')
    parser.add_argument('--opt-attempts',   type=int,   default=100, help='Number of spline mutations (attempts)')
    parser.add_argument('--target-xyz',     type=float, nargs=3, default=[5.0, 5.0, 10.0],  help='Target position for opposite atom (x y z) in Angstrom')
    parser.add_argument('--opt-wpos',       type=float, default=1.0, help='Weight of position loss term')
    parser.add_argument('--opt-wforce',     type=float, default=0.0, help='Weight of force loss term')
    parser.add_argument('--opt-fsafe',      type=float, default=5.0, help='Safe force threshold on anchor (eV/A)')
    parser.add_argument('--opt-temp0',      type=float, default=0.0, help='Initial temperature for simulated annealing')
    parser.add_argument('--opt-cooling',         type=float, default=1.0, help='Cooling factor per attempt')
    parser.add_argument('--opt-step-cooling',    type=float, default=1.0, help='Cooling factor per attempt')
    parser.add_argument('--opt-step0',      type=float, default=5.0, help='Mutation step size scale (Angstrom)')
    parser.add_argument('--opt-seed',       type=int,   default=45454, help='Random seed for optimization')
    parser.add_argument('--opt-plot-improvements', type=int, default=1, help='Plot trajectory for each improvement (1=yes, 0=no)')
    parser.add_argument('--opt-plot-all-trials', type=int, default=0, help='Plot trajectory for EVERY trial (1=yes, 0=no) - WARNING: many plots!')

    # Substrate and plotting parameters
    parser.add_argument('--substrate-shift', type=float, nargs=3, default=[-5.0, -10.0, 2.5], help='Shift of substrate atoms (x y z) in Angstroms')
    parser.add_argument('--plot-substrate', type=int, default=1, help='Plot substrate atoms (1=yes, 0=no)')
    parser.add_argument('--plot-molecule-atoms', type=int, default=0, help='Plot molecule atoms (1=yes, 0=no)')
    parser.add_argument('--plot-molecule-bonds', type=int, default=1, help='Plot molecule bonds (1=yes, 0=no)')
    parser.add_argument('--plot-connectors', type=int, default=0, help='Plot connector lines between trajectories (1=yes, 0=no)')
    
    # Molecular dynamics parameters
    parser.add_argument('--dt',    type=float, default=0.05,help='Time step for MD integration (fs)')
    parser.add_argument('--niter', type=int,   default=1000,help='Maximum number of MD iterations per spline point')
    parser.add_argument('--fconv', type=float, default=1e-4,  help='Force convergence criterion (eV/A)')
    parser.add_argument('--flim',  type=float, default=1000.0, help='Force limit to prevent instabilities (eV/A)')
    
    # Molecule positioning
    parser.add_argument('--shift', type=float, nargs=3, default=[-1.0, -1.0, 10.0],  help='Initial shift of molecule (x y z) in Angstroms')
    
    # Trajectory and output
    parser.add_argument('--trj',        type=str, default='trj_debug_relax.xyz',  help='Base name for trajectory output files')
    parser.add_argument('--trj-steps',  type=int, default=0, help='Save every MD step to trajectory (debugging, huge files)')
    parser.add_argument('--nconf',      type=int, default=40,  help='Number of spline configurations to run')
    parser.add_argument('--plot',       type=int, default=1, help='Generate plots of atom trajectories and forces')
    
    # Force field parameters
    parser.add_argument('--invert-coulomb', type=int, default=1, help='Invert Coulomb interaction sign')
    parser.add_argument('--mmff',   type=int, default=1, help='Use MMFF force field')
    parser.add_argument('--epairs', type=int, default=0, help='Include electron pairs in force field')
    
    # Debug options
    parser.add_argument('--verbose', type=int, default=1,  help='Verbosity level (0=quiet, 1=normal, 2=debug)')
    parser.add_argument('--debug',   type=int, default=0,  help='Debug level')
    
    return parser.parse_args()

# Parse arguments
args = parse_args()

# Convert arguments to script variables
xyz_name = args.xyz
surf_name = None if args.no_surface else args.surface
spline_file = args.spline
iAnchor = args.anchor
Kanchor = args.kanchor
shift = np.array(args.shift)
niter_max = args.niter
dt = args.dt
Fconv = args.fconv
Flim = args.flim
trj_name = args.trj
trj_steps = args.trj_steps
nconf = args.nconf
bMMFF = bool(args.mmff)
bEpairs = bool(args.epairs)
verbosity = args.verbose
idebug = args.debug

# New plotting and substrate parameters
substrate_shift = np.array(args.substrate_shift)  # Now a vector (x,y,z)
plot_substrate = bool(args.plot_substrate)
plot_molecule_atoms = bool(args.plot_molecule_atoms)
plot_molecule_bonds = bool(args.plot_molecule_bonds)
plot_connectors = bool(args.plot_connectors)

do_optimize = bool(args.optimize)
opt_outdir = args.opt_outdir
opt_attempts = args.opt_attempts
opt_target = args.target_xyz
opt_wpos = args.opt_wpos
opt_wforce = args.opt_wforce
opt_fsafe = args.opt_fsafe
opt_temp0 = args.opt_temp0
opt_cooling = args.opt_cooling
opt_step_cooling = args.opt_step_cooling
opt_step0 = args.opt_step0
opt_seed = args.opt_seed
opt_plot_improvements = bool(args.opt_plot_improvements)
opt_plot_all_trials = bool(args.opt_plot_all_trials)

# Clean only specific trajectory files we're writing to
for f in [trj_name]:
    if os.path.exists(f):
        print(f"Removing old trajectory file: {f}")
        os.remove(f)

# Always clean step-by-step trajectory files (C++ creates them regardless)
for f in [f"{trj_name}_steps.xyz"]:
    if os.path.exists(f):
        print(f"Removing old step trajectory file: {f}")
        os.remove(f)

print(f"=== Spline Scan Configuration ===")
print(f"Molecule: {xyz_name}")
print(f"Surface: {surf_name}")
print(f"Anchor atom: {iAnchor}, Kanchor: {Kanchor} eV/A²")
print(f"MD: dt={dt} fs, niter={niter_max}, Fconv={Fconv} eV/A")
print(f"Shift: {shift} Å")
print(f"Trajectory: {trj_name}")
print("="*40)

# Headless test reimplementing MolGUIapp command line:
# ./$name -x common_resources/xyz/PTCDA -g common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top -nPBC 0,0,0 -shift 0.0,0.0,5.0 -plq_factors 1.0,1.0,-1.0,1.0 -tipSpline tipSpline.dat -tipAnchor 28,20.0

mmff.setVerbosity(verbosity=verbosity, idebug=idebug)

# init (MMFF.py init expects filenames without extension in many places; the C++ init() adds .xyz internally)
mmff.init(
    xyz_name = xyz_name,
    surf_name= surf_name,
    bMMFF=bMMFF,
    bEpairs=bEpairs,
    bUFF=False,
    nPBC=(0,0,0),
    gridStep=0.1,
)
mmff.getBuffs()

# DEBUG: Check if molecule is already scrambled in ffl.apos
print("=== DEBUG: Checking ffl.apos right after initialization ===")
print("First 5 atoms positions in ffl.apos:")
for i in range(5):
    print(f"  Atom {i}: ({mmff.apos[i,0]:.3f}, {mmff.apos[i,1]:.3f}, {mmff.apos[i,2]:.3f})")

# shift molecule to specified position
print(f"Before shift: anchor atom {iAnchor} position: ({mmff.apos[iAnchor,0]:.3f}, {mmff.apos[iAnchor,1]:.3f}, {mmff.apos[iAnchor,2]:.3f})")
mmff.shift_atoms_ax(shift, sel=None)
print(f"After shift: anchor atom {iAnchor} position: ({mmff.apos[iAnchor,0]:.3f}, {mmff.apos[iAnchor,1]:.3f}, {mmff.apos[iAnchor,2]:.3f})")

# Also check if the molecule structure looks like PTCDA at all
print("First 5 atoms positions:")
for i in range(5):
    print(f"  Atom {i}: ({mmff.apos[i,0]:.3f}, {mmff.apos[i,1]:.3f}, {mmff.apos[i,2]:.3f})")

# -plq_factors 1.0,1.0,-1.0,1.0  => invert Coulomb sign
if args.invert_coulomb:
    mmff.PLQs[:,2] *= -3.0

# DEBUG: Run spline scan with specified parameters
ts = np.linspace(0, 1, nconf)  # Generate time points for spline

if do_optimize:
    if opt_target is None:
        raise RuntimeError("--optimize=1 requires --opt-target x y z")
    print("=== DEBUG: Starting spline optimization ===")
    print(f"opt_outdir={opt_outdir} opt_attempts={opt_attempts} opt_target={opt_target}")
    print(f"opt_wpos={opt_wpos} opt_wforce={opt_wforce} opt_fsafe={opt_fsafe}")
    print(f"opt_temp0={opt_temp0} opt_cooling={opt_cooling} opt_step0={opt_step0} opt_seed={opt_seed}")
    from TipSplineOptimizer import TipSplineSAOptimizer
    opt = TipSplineSAOptimizer(
        ia_anchor=iAnchor,
        ia_opposite=27-1,
        target_pos=opt_target,
        w_pos=opt_wpos,
        w_force=opt_wforce,
        f_safe=opt_fsafe,
        temp0=opt_temp0,
        cooling=opt_cooling,
        step0=opt_step0,
        step_cooling=opt_step_cooling,
        freeze_n_ends=3,
        seed=opt_seed,
    )

    scan_kwargs = {
        'iAnchor': iAnchor,
        'Kanchor': Kanchor,
        'trjName': None,
        'nPBC': (0,0,0),
        'niter_max': niter_max,
        'dt': dt,
        'Fconv': Fconv,
        'Flim': Flim,
    }
    best_spline, best = opt.optimize(
        mmff_instance=mmff,
        ts=ts,
        spline_init_fname=spline_file,
        out_dir=opt_outdir,
        n_attempts=opt_attempts,
        scan_kwargs=scan_kwargs,
        save_every_improvement=True,
        plot_improvements=opt_plot_improvements,
        plot_all_trials=opt_plot_all_trials,
        bold_best=True,
        substrate_shift=substrate_shift,
        surface_name=surf_name,
        plot_substrate=plot_substrate,
        plot_molecule_atoms=plot_molecule_atoms,
        plot_molecule_bonds=plot_molecule_bonds,
        opt_target=opt_target,
    )
    print(f"=== DEBUG: Optimization DONE best_spline={best_spline} best_E={best['E']:.6g} best_attempt={best['attempt']} ===")
    spline_file = best_spline

print("=== DEBUG: Running spline scan ===")
print(f"Initial anchor atom {iAnchor} position: ({mmff.apos[iAnchor,0]:.3f}, {mmff.apos[iAnchor,1]:.3f}, {mmff.apos[iAnchor,2]:.3f})")

# Always pass trjName to save final configurations
# Step-by-step saving is controlled separately in C++ by trj_fname
Es, FAs, APs, Fcon = mmff.scan_manipulation(
    ts,
    spline_fname=spline_file,
    iAnchor=iAnchor,
    Kanchor=Kanchor,
    trjName=trj_name,  # Always save final configurations
    nPBC=(0,0,0),
    niter_max=niter_max,
    dt=dt,
    Fconv=Fconv,
    Flim=Flim,
)

print(f"Final anchor atom {iAnchor} position: ({mmff.apos[iAnchor,0]:.3f}, {mmff.apos[iAnchor,1]:.3f}, {mmff.apos[iAnchor,2]:.3f})")
print(f"Energy: {Es[0]:.3f} eV")
print(f"Constraint force: ({Fcon[0,0]:.3f}, {Fcon[0,1]:.3f}, {Fcon[0,2]:.3f}) eV/A")
if trj_steps:
    print(f"Check step-by-step trajectory file: {trj_name}_steps.xyz")
    print(f"It contains {niter_max} MD steps per configuration")
else:
    print(f"Check final configurations trajectory file: {trj_name}")
    print(f"It contains {nconf} final relaxed configurations")

# Optional plotting
if args.plot and nconf > 1:
    print("Generating plots...")
    # NOTE: user asked anchored atom (29) and opposite (27) by their 1-based indices => convert to 0-based
    iaA = 29-1
    iaB = 27-1
    
    # Import plotUtils for molecule visualization
    import sys
    sys.path.append('../../')
    from pyBall import plotUtils
    
    # Get substrate atom positions if surface is loaded
    substrate_pos = None
    substrate_types = None
    
    if not args.no_surface:
        # Load actual substrate atoms from XYZ file - NO FALLBACKS!
        if surf_name and os.path.exists(f"{surf_name}.xyz"):
            print(f"Loading substrate atoms from {surf_name}.xyz")
            with open(f"{surf_name}.xyz", 'r') as f:
                lines = f.readlines()
            
            # Parse XYZ format
            natoms = int(lines[0].strip())
            substrate_pos = []
            substrate_types = []
            for line in lines[2:2+natoms]:  # Skip header and comment line
                parts = line.strip().split()
                if len(parts) >= 4:
                    atom_type = parts[0]  # Ca or F
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    substrate_pos.append([x, y, z])
                    substrate_types.append(atom_type)
            
            substrate_pos = np.array(substrate_pos)
            print(f"Loaded {len(substrate_pos)} substrate atoms from {surf_name}.xyz")
            print(f"Atom types: {dict(zip(*np.unique(substrate_types, return_counts=True)))}")
        else:
            print(f"ERROR: Surface file {surf_name}.xyz not found!")
            substrate_pos = None
            substrate_types = None
    
    # Get molecule bonds for visualization
    # Use the first configuration for initial position, last for final position
    initial_pos = APs[0]  # First configuration
    final_pos = APs[-1]   # Last configuration
    
    # Create simple bond list based on distance cutoff (you may need to adjust this)
    def create_bonds(positions, cutoff=1.6):
        """Create simple bond list based on distance cutoff"""
        natoms = len(positions)
        bonds = []
        for i in range(natoms):
            for j in range(i+1, natoms):
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist < cutoff:
                    bonds.append([i, j])
        return bonds
    
    bonds = create_bonds(initial_pos)
    
    # Plot trajectories and forces
    plt.figure(figsize=(12, 5))
    
    # Top view (xy)
    plt.subplot(1, 2, 1)
    plt.plot(APs[:, iaA, 0], APs[:, iaA, 1], 'r-', label=f'Anchor {iAnchor}', linewidth=2)
    plt.plot(APs[:, iaB, 0], APs[:, iaB, 1], 'b-', label=f'Opposite {iAnchor-2}', linewidth=2)
    
    # Add substrate atoms if available
    if substrate_pos is not None and plot_substrate:
        # Apply substrate shift (vector x,y,z)
        substrate_pos_shifted = substrate_pos.copy()
        substrate_pos_shifted[:, 0] += substrate_shift[0]  # x shift
        substrate_pos_shifted[:, 1] += substrate_shift[1]  # y shift  
        substrate_pos_shifted[:, 2] += substrate_shift[2]  # z shift
        
        # Find top surface atoms (highest z positions)
        top_z = np.max(substrate_pos_shifted[:, 2])
        print(f"Top surface z = {top_z:.3f} Å (after shift {substrate_shift})")
        
        # Filter atoms: only show those within 3.0Å of top surface
        height_cutoff = 3.0  # Å
        visible_mask = substrate_pos_shifted[:, 2] >= (top_z - height_cutoff)
        
        if not np.any(visible_mask):
            print("WARNING: No substrate atoms visible within 3.0Å of top surface!")
        else:
            visible_pos = substrate_pos_shifted[visible_mask]
            visible_types = np.array(substrate_types)[visible_mask]
            heights = visible_pos[:, 2]
            
            # Normalize heights for sizing (0.3 to 1.0 range)
            min_height = np.min(heights)
            max_height = np.max(heights)
            if max_height > min_height:
                normalized_sizes = 0.3 + 0.7 * (heights - min_height) / (max_height - min_height)
            else:
                normalized_sizes = np.ones_like(heights) * 0.5
            
            # Base sizes for different atom types
            ca_base_size = 40
            f_base_size = 20
            
            # Separate Ca and F atoms for different colors (magenta/green)
            ca_mask = visible_types == 'Ca'
            f_mask = visible_types == 'F'
            
            if np.any(ca_mask):
                ca_sizes = ca_base_size * normalized_sizes[ca_mask]
                plt.scatter(visible_pos[ca_mask, 0], visible_pos[ca_mask, 1], 
                           c='magenta', s=ca_sizes, alpha=0.7, label='Ca atoms')
            if np.any(f_mask):
                f_sizes = f_base_size * normalized_sizes[f_mask]
                plt.scatter(visible_pos[f_mask, 0], visible_pos[f_mask, 1], 
                           c='green', s=f_sizes, alpha=0.5, label='F atoms')
            
            print(f"Showing {np.sum(visible_mask)} substrate atoms (top {height_cutoff:.1f}Å)")
    
    # Add molecule geometry based on plotting parameters
    if plot_molecule_atoms:
        # Initial position (light gray)
        plotUtils.plotAtoms(initial_pos, colors='lightgray', sizes=30, axes=(0,1))
        # Final position (black)
        plotUtils.plotAtoms(final_pos, colors='black', sizes=40, axes=(0,1))
    
    if plot_molecule_bonds:
        # Initial position bonds (thin gray)
        plotUtils.plotBonds(links=bonds, ps=initial_pos, colors='lightgray', lws=0.5, axes=(0,1))
        # Final position bonds (thick black)
        plotUtils.plotBonds(links=bonds, ps=final_pos, colors='black', lws=1.5, axes=(0,1))
    
    # Add connecting lines between corresponding points on red and blue trajectories
    if plot_connectors:
        for i in range(nconf):
            # Connect anchor (red) to opposite (blue) atom at each configuration
            anchor_pos = APs[i, iaA, :2]  # x,y only
            opposite_pos = APs[i, iaB, :2]  # x,y only
            
            # Draw connecting line
            plt.plot([anchor_pos[0], opposite_pos[0]], 
                    [anchor_pos[1], opposite_pos[1]], 
                    'k--', alpha=0.5, linewidth=1)
    
    plt.xlabel('x (Å)')
    plt.ylabel('y (Å)')
    plt.title('Top View (xy)')
    plt.legend()
    plt.axis('equal')
    
    # Side view (xz)
    plt.subplot(1, 2, 2)
    plt.plot(APs[:, iaA, 0], APs[:, iaA, 2], 'r-', label=f'Anchor {iAnchor}', linewidth=2)
    plt.plot(APs[:, iaB, 0], APs[:, iaB, 2], 'b-', label=f'Opposite {iAnchor-2}', linewidth=2)
    
    # Add substrate atoms if available
    if substrate_pos is not None and plot_substrate:
        # Apply substrate shift (same vector as top view)
        substrate_pos_shifted = substrate_pos.copy()
        substrate_pos_shifted[:, 0] += substrate_shift[0]  # x shift
        substrate_pos_shifted[:, 1] += substrate_shift[1]  # y shift  
        substrate_pos_shifted[:, 2] += substrate_shift[2]  # z shift
        
        # Filter atoms: only show those within 3.0Å of top surface
        top_z = np.max(substrate_pos_shifted[:, 2])
        height_cutoff = 3.0  # Å
        visible_mask = substrate_pos_shifted[:, 2] >= (top_z - height_cutoff)
        
        if np.any(visible_mask):
            visible_pos = substrate_pos_shifted[visible_mask]
            visible_types = np.array(substrate_types)[visible_mask]
            heights = visible_pos[:, 2]
            
            # Normalize heights for sizing (0.3 to 1.0 range)
            min_height = np.min(heights)
            max_height = np.max(heights)
            if max_height > min_height:
                normalized_sizes = 0.3 + 0.7 * (heights - min_height) / (max_height - min_height)
            else:
                normalized_sizes = np.ones_like(heights) * 0.5
            
            # Base sizes for different atom types
            ca_base_size = 40
            f_base_size = 20
            
            # Separate Ca and F atoms for different colors (magenta/green)
            ca_mask = visible_types == 'Ca'
            f_mask = visible_types == 'F'
            
            if np.any(ca_mask):
                ca_sizes = ca_base_size * normalized_sizes[ca_mask]
                plt.scatter(visible_pos[ca_mask, 0], visible_pos[ca_mask, 2], 
                           c='magenta', s=ca_sizes, alpha=0.7, label='Ca atoms')
            if np.any(f_mask):
                f_sizes = f_base_size * normalized_sizes[f_mask]
                plt.scatter(visible_pos[f_mask, 0], visible_pos[f_mask, 2], 
                           c='green', s=f_sizes, alpha=0.5, label='F atoms')
    
    # Add molecule geometry based on plotting parameters
    if plot_molecule_atoms:
        # Initial position (light gray)
        plotUtils.plotAtoms(initial_pos, colors='lightgray', sizes=30, axes=(0,2))
        # Final position (black)
        plotUtils.plotAtoms(final_pos, colors='black', sizes=40, axes=(0,2))
    
    if plot_molecule_bonds:
        # Initial position bonds (thin gray)
        plotUtils.plotBonds(links=bonds, ps=initial_pos, colors='lightgray', lws=0.5, axes=(0,2))
        # Final position bonds (thick black)
        plotUtils.plotBonds(links=bonds, ps=final_pos, colors='black', lws=1.5, axes=(0,2))
    
    # Add connecting lines between corresponding points on red and blue trajectories
    if plot_connectors:
        for i in range(nconf):
            # Connect anchor (red) to opposite (blue) atom at each configuration
            anchor_pos = APs[i, iaA, [0, 2]]  # x,z only
            opposite_pos = APs[i, iaB, [0, 2]]  # x,z only
            
            # Draw connecting line
            plt.plot([anchor_pos[0], opposite_pos[0]], 
                    [anchor_pos[1], opposite_pos[1]], 
                    'k--', alpha=0.5, linewidth=1)
    
    plt.xlabel('x (Å)')
    plt.ylabel('z (Å)')
    plt.title('Side View (xz)')
    plt.legend()
    plt.axis('equal')
    
    plt.tight_layout()
    plt.savefig(f'{trj_name}_trajectories.png', dpi=150)
    print(f"Saved trajectories plot to {trj_name}_trajectories.png")
    
    # Plot energies
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(Es, 'k-')
    plt.xlabel('Configuration')
    plt.ylabel('Energy (eV)')
    plt.title('Energy vs Configuration')
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(Fcon[:, 0], 'r-', label='Fx')
    plt.plot(Fcon[:, 1], 'g-', label='Fy')
    plt.plot(Fcon[:, 2], 'b-', label='Fz')
    plt.xlabel('Configuration')
    plt.ylabel('Constraint Force (eV/A)')
    plt.title('Constraint Forces vs Configuration')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(f'{trj_name}_energies.png', dpi=150)
    print(f"Saved energy plot to {trj_name}_energies.png")

print("=== Simulation Complete ===")
