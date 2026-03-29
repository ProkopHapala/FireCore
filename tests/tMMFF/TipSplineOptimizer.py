import numpy as np
import os
import json
import time
import matplotlib.pyplot as plt
import sys

# Add path to import plotUtils
sys.path.append('../../')
from pyBall import plotUtils


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


def _plot_single_improvement(APs, Fcon, target_pos, ia_anchor, ia_opposite, attempt, out_dir, mmff_instance, substrate_shift, surface_name=None, plot_substrate=True, plot_molecule_atoms=True, plot_molecule_bonds=True, opt_target=None, current_control_pts=None, suffix="", title_extra="", Es=None):
    """Plot comprehensive 4-panel trajectory for a single improvement"""
    # Fixed axis ranges for stable movie
    xlim = (-5, 20)
    ylim = (-5, 10)
    zlim = (5, 20)
    Flim = 3.0
    
    fig = plt.figure(figsize=(16, 12))
    
    # Load substrate atoms from XYZ file (same as main script)
    substrate_pos = None
    substrate_types = None
    if plot_substrate and surface_name:
        xyz_file = f"{surface_name}.xyz"
        if os.path.exists(xyz_file):
            print(f"Loading substrate atoms from {xyz_file}")
            with open(xyz_file, 'r') as f:
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
            print(f"Loaded {len(substrate_pos)} substrate atoms from {xyz_file}")
        else:
            print(f"WARNING: Substrate file {xyz_file} not found!")
    
    # Create molecule bonds (same as main script)
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
    
    bonds = create_bonds(APs[0])  # Use initial position for bond detection
    
    # Load initial spline for reference
    initial_spline_file = os.path.join(os.path.dirname(out_dir), '..', 'tipSpline_x.dat')
    initial_pts = None
    initial_trajectory = None
    if os.path.exists(initial_spline_file):
        initial_pts = load_spline_dat(initial_spline_file)
        print(f"Loaded initial spline with {len(initial_pts)} control points")
        
        # Load initial trajectory data to plot baseline
        initial_spline_tmp = os.path.join(out_dir, "initial_spline.dat")
        save_spline_dat(initial_spline_tmp, initial_pts, header="# Initial spline from tipSpline_x.dat")
        try:
            # Run scan with initial spline to get baseline trajectory
            print("Computing initial trajectory for baseline...")
            initial_Es, initial_FAs, initial_APs, initial_Fcon = mmff_instance.scan_manipulation(
                np.linspace(0, 1, 40), spline_fname=initial_spline_tmp, **{}
            )
            initial_trajectory = np.array(initial_APs)
            print(f"Initial trajectory computed: {len(initial_trajectory)} configurations")
        except Exception as e:
            print(f"WARNING: Could not compute initial trajectory: {e}")
            initial_trajectory = None
    else:
        print(f"WARNING: Initial spline file {initial_spline_file} not found!")
    
    # Calculate fitness values for display
    n_conf = len(APs)
    final_opposite_pos = APs[-1, ia_opposite]
    pos_error = np.linalg.norm(final_opposite_pos - target_pos)  # Distance to target
    
    # Get force magnitude (assuming Fcon contains constraint forces)
    if Fcon.shape[1] == 3:  # Anchor forces only
        max_force = np.max(np.sqrt(np.sum(Fcon**2, axis=1)))
    else:  # All atom forces - extract anchor
        F_anchor = Fcon[:, ia_anchor*3:(ia_anchor+1)*3]
        max_force = np.max(np.sqrt(np.sum(F_anchor**2, axis=1)))
    
    total_fitness = pos_error + max_force
    
    # Panel 1: Top view (xy) with substrate and molecule bonds
    ax1 = fig.add_subplot(2, 2, 1)
    
    # Plot substrate if available
    if substrate_pos is not None:
        # Apply substrate shift (vector x,y,z)
        substrate_pos_shifted = substrate_pos.copy()
        substrate_pos_shifted[:, 0] += substrate_shift[0]  # x shift
        substrate_pos_shifted[:, 1] += substrate_shift[1]  # y shift  
        substrate_pos_shifted[:, 2] += substrate_shift[2]  # z shift
        
        # Find top surface atoms (highest z positions)
        top_z = np.max(substrate_pos_shifted[:, 2])
        
        # Filter atoms: only show those within 3.0Å of top surface
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
                ax1.scatter(visible_pos[ca_mask, 0], visible_pos[ca_mask, 1],  c='magenta', s=ca_sizes, alpha=0.7, label='Ca atoms')
            if np.any(f_mask):
                f_sizes = f_base_size * normalized_sizes[f_mask]
                ax1.scatter(visible_pos[f_mask, 0], visible_pos[f_mask, 1],  c='green', s=f_sizes, alpha=0.5, label='F atoms')
    
    # Plot molecule bonds (not atoms) as requested
    if plot_molecule_bonds:
        # Initial position bonds (thin gray)
        plotUtils.plotBonds(links=bonds, ps=APs[0], colors='gray',   lws=2.0, axes=(0,1))
        # Final position bonds (thick black)
        plotUtils.plotBonds(links=bonds, ps=APs[-1], colors='black', lws=2.0, axes=(0,1))
    
    # Plot initial spline reference if available
    if initial_pts is not None:
        initial_pts_array = np.array(initial_pts)
        ax1.plot(initial_pts_array[:, 0], initial_pts_array[:, 1], 'k.--',  linewidth=1, markersize=4, alpha=0.5, label='Initial spline')
    
    # Use the passed control points from optimizer - DO NOT LOAD FROM FILE!
    if current_control_pts is not None:
        print(f"🔍 DEBUG PLOTTING: Using provided control points: {len(current_control_pts)} points")
        # DEBUG: Print control point 3 to verify it's changing
        if len(current_control_pts) > 3:
            pt3 = current_control_pts[3]
            print(f"🔍 DEBUG PLOTTING: Control point 3 = ({pt3[0]:.6f}, {pt3[1]:.6f}, {pt3[2]:.6f})")
        else:
            print(f"🔍 DEBUG PLOTTING: Control points array too short: {len(current_control_pts)}")
    else:
        print("🔍 DEBUG PLOTTING: WARNING: No control points provided to plot!")
        current_control_pts = None
    
    # Plot initial trajectory as baseline reference
    if initial_trajectory is not None:
        ax1.plot(initial_trajectory[:, ia_anchor, 0], initial_trajectory[:, ia_anchor, 1],  'g--', linewidth=1, alpha=0.5, label='Initial trajectory (anchor)')
        ax1.plot(initial_trajectory[:, ia_opposite, 0], initial_trajectory[:, ia_opposite, 1],  'c--', linewidth=1, alpha=0.5, label='Initial trajectory (opposite)')
    
    # Plot current trajectory (interpolated from spline)
    ax1.plot(APs[:, ia_anchor, 0], APs[:, ia_anchor, 1], 'o-',  color='red', linewidth=2, markersize=4, label=f'Anchor {ia_anchor}')
    ax1.plot(APs[:, ia_opposite, 0], APs[:, ia_opposite, 1], 'o-',  color='blue', linewidth=2, markersize=4, label=f'Opposite {ia_opposite}')
    
    # Plot actual spline control points (THIS IS WHAT YOU WANTED!)
    if current_control_pts is not None:
        current_control_array = np.array(current_control_pts)
        ax1.plot(current_control_array[:, 0], current_control_array[:, 1], 'rx', markersize=10, markeredgewidth=2, label='Spline control points')
        
        # Add ordering numbers to control points
        for i in range(len(current_control_array)):
            x, y = current_control_array[i, 0], current_control_array[i, 1]
            ax1.text(x, y + 0.3, str(i), fontsize=8, ha='center', color='red', fontweight='bold')
    
    # Plot big empty blue circle for target
    if opt_target is not None:
        ax1.plot(opt_target[0], opt_target[1], 'o', color='blue', markersize=20,   fillstyle='none', markeredgewidth=3, label='Target')
    
    # Add fitness text
    fitness_text = f'Fitness: {total_fitness:.3f}\nPos: {pos_error:.3f}\nForce: {max_force:.3f}'
    ax1.text(0.02, 0.98, fitness_text, transform=ax1.transAxes,  verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax1.set_title(f'Improvement {attempt} - Top View (xy)')
    ax1.set_xlabel('x (Å)')
    ax1.set_ylabel('y (Å)')
    ax1.set_xlim(xlim)
    ax1.set_ylim(ylim)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Side view (xz) with substrate and molecule bonds
    ax2 = fig.add_subplot(2, 2, 2)
    
    # Plot substrate if available
    if substrate_pos is not None:
        if np.any(visible_mask):
            visible_pos_xz = substrate_pos_shifted[visible_mask]
            ca_mask_xz = visible_types == 'Ca'
            f_mask_xz = visible_types == 'F'
            
            if np.any(ca_mask_xz):
                ca_sizes_xz = ca_base_size * normalized_sizes[ca_mask_xz]
                ax2.scatter(visible_pos_xz[ca_mask_xz, 0], visible_pos_xz[ca_mask_xz, 2], c='magenta', s=ca_sizes_xz, alpha=0.7)
            if np.any(f_mask_xz):
                f_sizes_xz = f_base_size * normalized_sizes[f_mask_xz]
                ax2.scatter(visible_pos_xz[f_mask_xz, 0], visible_pos_xz[f_mask_xz, 2],  c='green', s=f_sizes_xz, alpha=0.5)
    
    # Plot molecule bonds (not atoms) as requested
    if plot_molecule_bonds:
        # Initial position bonds (thin gray)
        plotUtils.plotBonds(links=bonds, ps=APs[0], colors='lightgray', lws=0.5, axes=(0,2))
        # Final position bonds (thick black)
        plotUtils.plotBonds(links=bonds, ps=APs[-1], colors='black', lws=1.5, axes=(0,2))
    
    # Plot initial spline reference if available
    if initial_pts is not None:
        ax2.plot(initial_pts_array[:, 0], initial_pts_array[:, 2], 'k.--',   linewidth=1, markersize=4, alpha=0.5, label='Initial spline')
    
    # Plot initial trajectory as baseline reference
    if initial_trajectory is not None:
        ax2.plot(initial_trajectory[:, ia_anchor, 0], initial_trajectory[:, ia_anchor, 2],      'g--', linewidth=1, alpha=0.5, label='Initial trajectory (anchor)')
        ax2.plot(initial_trajectory[:, ia_opposite, 0], initial_trajectory[:, ia_opposite, 2],  'c--', linewidth=1, alpha=0.5, label='Initial trajectory (opposite)')
    
    # Plot current trajectory (interpolated from spline)
    ax2.plot(APs[:, ia_anchor, 0], APs[:, ia_anchor, 2],     'o-',  color='red', linewidth=2, markersize=4, label=f'Anchor {ia_anchor}')
    ax2.plot(APs[:, ia_opposite, 0], APs[:, ia_opposite, 2], 'o-', color='blue', linewidth=2, markersize=4, label=f'Opposite {ia_opposite}')
    
    # Plot actual spline control points (THIS IS WHAT YOU WANTED!)
    if current_control_pts is not None:
        current_control_array = np.array(current_control_pts)
        ax2.plot(current_control_array[:, 0], current_control_array[:, 2], 'rx', markersize=10, markeredgewidth=2, label='Spline control points')
        
        # Add ordering numbers to control points
        for i in range(len(current_control_array)):
            x, z = current_control_array[i, 0], current_control_array[i, 2]
            ax2.text(x, z + 0.3, str(i), fontsize=8, ha='center', color='red', fontweight='bold')
    
    # Plot big empty blue circle for target
    if opt_target is not None:
        # Use actual z position from opt_target
        ax2.plot(opt_target[0], opt_target[2], 'o', color='blue', markersize=20,  fillstyle='none', markeredgewidth=3, label='Target')
    
    # Add fitness text
    ax2.text(0.02, 0.98, fitness_text, transform=ax2.transAxes, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax2.set_title(f'Improvement {attempt} - Side View (xz)')
    ax2.set_xlabel('x (Å)')
    ax2.set_ylabel('z (Å)')
    ax2.set_xlim(xlim)
    ax2.set_ylim(zlim)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Energy along trajectory - Use actual energy data from scan_manipulation
    ax3 = fig.add_subplot(2, 2, 3)
    
    # Use actual energy data from scan_manipulation
    if Es is not None:
        ax3.plot(range(len(Es)), Es, '.-', linewidth=1.5, label='System Energy')
        ax3.set_ylabel('Energy (eV)')
        # Auto-scale y-axis based on actual energy range
        e_min, e_max = np.min(Es), np.max(Es)
        e_range = e_max - e_min
        ax3.set_ylim(e_min - 0.1*e_range, e_max + 0.1*e_range)
    else:
        raise ValueError("No energy data provided to plotting function!")
    
    ax3.set_title(f'Improvement {attempt} - Energy Along Trajectory')
    ax3.set_xlabel('Configuration')
    ax3.set_xlim(0, len(Es)-1)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # Panel 4: Force components on anchor atom
    ax4 = fig.add_subplot(2, 2, 4)
    
    # Extract force components for anchor atom from constraint forces
    # Fcon has shape (n_conf, n_atoms*3) in some formats, or (n_conf, 3) for anchor only
    if Fcon.shape[1] == 3:  # Anchor forces only
        F_anchor = Fcon
        F_mag = np.sqrt(np.sum(F_anchor**2, axis=1))
        
        ax4.plot(range(n_conf), F_anchor[:, 0], 'r-', label='Fx', linewidth=1.5)
        ax4.plot(range(n_conf), F_anchor[:, 1], 'g-', label='Fy', linewidth=1.5)
        ax4.plot(range(n_conf), F_anchor[:, 2], 'b-', label='Fz', linewidth=1.5)
        ax4.plot(range(n_conf), F_mag, 'k-', label='|F|', linewidth=2)
    else:  # All atom forces - extract anchor
        n_atoms = Fcon.shape[1] // 3
        F_anchor = Fcon[:, ia_anchor*3:(ia_anchor+1)*3]
        F_mag = np.sqrt(np.sum(F_anchor**2, axis=1))
        
        ax4.plot(range(n_conf), F_anchor[:, 0], 'r-', label='Fx', linewidth=1.5)
        ax4.plot(range(n_conf), F_anchor[:, 1], 'g-', label='Fy', linewidth=1.5)
        ax4.plot(range(n_conf), F_anchor[:, 2], 'b-', label='Fz', linewidth=1.5)
        ax4.plot(range(n_conf), F_mag, '.-k', label='|F|', linewidth=2)
    
    ax4.set_xlabel('Configuration index')
    ax4.set_ylabel('Force (eV/Å)')
    ax4.set_title(f'Improvement {attempt} - Forces on Anchor {ia_anchor}')
    ax4.set_xlim(0, len(F_anchor[:, 0])-1)  # Fixed x-axis
    ax4.set_ylim(-Flim,Flim)  # Fixed y-axis for consistent movie frames
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    fig.tight_layout()
    filename = os.path.join(out_dir, f'improvement_{attempt:05d}{suffix}.png')
    fig.savefig(filename, dpi=150)
    plt.close(fig)
    print(f"Saved comprehensive improvement plot: {filename}")


class TipSplineSAOptimizer:
    def __init__(
        self,
        ia_anchor,
        ia_opposite,
        target_pos,
        w_pos=1.0,
        w_force=1.0,
        f_safe=5.0,
        temp0=1.0,
        cooling=0.995,
        step0=0.2,
        step_cooling=1.0,
        freeze_n_ends=3,
        seed=0,
    ):
        self.ia_anchor = int(ia_anchor)
        self.ia_opposite = int(ia_opposite)
        self.target_pos = np.array(target_pos, dtype=np.float64)
        self.w_pos = float(w_pos)
        self.w_force = float(w_force)
        self.f_safe = float(f_safe)
        self.temp0 = float(temp0)
        self.cooling = float(cooling)
        self.step0 = float(step0)
        self.step_cooling = float(step_cooling)
        self.freeze_n_ends = int(freeze_n_ends)
        self.rng = np.random.default_rng(int(seed))

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

        l_tot = self.w_pos * l_pos + self.w_force * l_force
        return l_tot, l_pos, l_force, fmax, rB

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
        return pts2, ip, d

    def _accept(self, Eold, Enew, temp):
        # Metropolis acceptance criterion:
        # 1. Always accept better solutions
        # 2. Accept worse solutions with probability exp(-ΔE/temp)
        # 3. Temperature controls acceptance of worse solutions ONLY
        
        if Enew <= Eold:
            return True  # Better or equal to current - always accept
        
        # Worse solution: accept with probability based on temperature
        if temp > 1e-6:  # Only if temperature is not zero
            dE = Enew - Eold
            prob = np.exp(-dE / temp)
            return self.rng.random() < prob
        
        return False  # Never accept worse when temp=0

    def optimize(
        self,
        mmff_instance,
        ts,
        spline_init_fname,
        out_dir,
        n_attempts,
        scan_kwargs,
        save_every_improvement=False,
        plot_improvements=False,
        plot_all_trials=False,
        bold_best=False,
        substrate_shift=None,
        surface_name=None,
        plot_substrate=True,
        plot_molecule_atoms=True,
        plot_molecule_bonds=True,
        opt_target=None,
        spline_header="# tipSpline.dat (x y z) in Angstrom",
    ):
        scan_kwargs = {} if scan_kwargs is None else dict(scan_kwargs)
        os.makedirs(out_dir, exist_ok=True)

        log_jsonl = os.path.join(out_dir, "attempts.jsonl")
        meta_json = os.path.join(out_dir, "meta.json")

        pts0 = load_spline_dat(spline_init_fname)
        pts = pts0.copy()

        # CRITICAL: Save initial atom positions for consistent reset
        initial_atom_pos = mmff_instance.apos[:, :].copy()
        
        # Evaluate initial
        spline_tmp = os.path.join(out_dir, "spline_current.dat")
        save_spline_dat(spline_tmp, pts, header=spline_header)

        print("TipSplineSAOptimizer.optimize(): evaluating initial spline")
        
        # CRITICAL: Reset molecule to initial position for consistent initialization
        # This ensures attempt -1 uses the same starting state as all other attempts
        mmff_instance.apos[:, :] = initial_atom_pos
        
        Es, FAs, APs, Fcon = mmff_instance.scan_manipulation(ts, spline_fname=spline_tmp, **scan_kwargs)
        Ebest, lpos_best, lforce_best, fmax_best, rB_best = self._loss_components(APs, Fcon)
        
        # Initialize best energy ever tracking
        self._best_E_ever = Ebest

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
                Es=Es  # Add actual energy data!
            )

        improvements = []
        if save_every_improvement:
            improvements.append(best)

        meta = {
            'time_start': time.time(),
            'spline_init': spline_init_fname,
            'n_attempts': int(n_attempts),
            'ia_anchor': self.ia_anchor,
            'ia_opposite': self.ia_opposite,
            'target_pos': self.target_pos,
            'weights': {'w_pos': self.w_pos, 'w_force': self.w_force},
            'force': {'f_safe': self.f_safe},
            'sa': {'temp0': self.temp0, 'cooling': self.cooling, 'step0': self.step0},
            'freeze_n_ends': self.freeze_n_ends,
        }
        with open(meta_json, 'w') as f:
            json.dump(meta, f, indent=2, default=_json_default)

        with open(log_jsonl, 'w') as flog:
            rec0 = {
                'attempt': -1,
                'accepted': True,
                'improved': True,
                'temp': self.temp0,
                'mut_i': None,
                'mut_d': [0.0, 0.0, 0.0],
                'E_total': Ebest,
                'E_pos': lpos_best,
                'E_force': lforce_best,
                'fmax': fmax_best,
                'r_opposite_final': rB_best,
            }
            flog.write(json.dumps(rec0, default=_json_default) + "\n")

            Ecur = Ebest
            temp = self.temp0
            current_step = self.step0

            for it in range(int(n_attempts)):
                pts_new, ip, d = self._propose(pts, temp, current_step)
                
                # CRITICAL: Use unique spline filename for each attempt
                # This prevents file caching/overwrite issues that could cause non-deterministic behavior
                spline_attempt = f"{out_dir}/spline_attempt_{it:05d}.dat"
                save_spline_dat(spline_attempt, pts_new, header=spline_header)

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
                        save_spline_dat(fn, pts_new, header=spline_header)
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
                                Es=Es  # Add actual energy data!
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
                        Es=Es  # Add actual energy data!
                    )

                rec = {
                    'attempt': it,
                    'accepted': bool(accepted),
                    'improved': bool(improved),
                    'temp': float(temp),
                    'mut_i': int(ip),
                    'mut_d': d,
                    'E_total': float(Enew),
                    'E_pos': float(lpos),
                    'E_force': float(lforce),
                    'fmax': float(fmax),
                    'r_opposite_final': rB,
                }
                flog.write(json.dumps(rec, default=_json_default) + "\n")
                flog.flush()

                temp *= self.cooling
                current_step *= self.step_cooling

        # Save final best spline
        best_spline = os.path.join(out_dir, "best_final.dat")
        save_spline_dat(best_spline, best['pts'], header=spline_header)
        np.save(os.path.join(out_dir, "best_final_APs.npy"), best['APs'])
        np.save(os.path.join(out_dir, "best_final_Fcon.npy"), best['Fcon'])
        with open(os.path.join(out_dir, "best_final.json"), 'w') as f:
            json.dump({k: best[k] for k in ['E', 'l_pos', 'l_force', 'fmax', 'r_opposite_final', 'attempt']}, f, indent=2, default=_json_default)

        if plot_improvements and (len(improvements) > 0):
            print(f"Plotting improvements: {len(improvements)}")
            fig = plt.figure(figsize=(12, 5))

            # XY
            ax1 = fig.add_subplot(1, 2, 1)
            for b in improvements:
                AP = b['APs']
                lw = 1.0
                alpha = 0.2
                if bold_best and (b is improvements[-1]):
                    lw = 2.5
                    alpha = 0.9
                ax1.plot(AP[:, self.ia_anchor, 0], AP[:, self.ia_anchor, 1], 'r-', lw=lw, alpha=alpha)
                ax1.plot(AP[:, self.ia_opposite, 0], AP[:, self.ia_opposite, 1], 'b-', lw=lw, alpha=alpha)
            ax1.plot(self.target_pos[0], self.target_pos[1], 'kx')
            ax1.set_title('Optimization improvements (xy)')
            ax1.axis('equal')

            # XZ
            ax2 = fig.add_subplot(1, 2, 2)
            for b in improvements:
                AP = b['APs']
                lw = 1.0
                alpha = 0.2
                if bold_best and (b is improvements[-1]):
                    lw = 2.5
                    alpha = 0.9
                ax2.plot(AP[:, self.ia_anchor, 0], AP[:, self.ia_anchor, 2], 'r-', lw=lw, alpha=alpha)
                ax2.plot(AP[:, self.ia_opposite, 0], AP[:, self.ia_opposite, 2], 'b-', lw=lw, alpha=alpha)
            ax2.plot(self.target_pos[0], self.target_pos[2], 'kx')
            ax2.set_title('Optimization improvements (xz)')
            ax2.axis('equal')

            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, 'improvements.png'), dpi=150)
            plt.close(fig)

        return best_spline, best
