import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_truss( points, bonds, ax=None, edge_color='k', edge_alpha=1.0, point_color='b', point_size=20, color_by_stiffness=False, cmap='viridis',  show_colorbar=True, ks=None, label=None):
    """
    Plot the truss efficiently using LineCollection.
    
    Args:
        points: Nx3 array of point coordinates
        bonds: List of (i,j) tuples defining connections
        ax: matplotlib axis to plot on (creates new if None)
        edge_color: color for edges (ignored if color_by_stiffness=True)
        edge_alpha: transparency for edges
        point_color: color for points
        point_size: size of points
        color_by_stiffness: if True, color edges by their stiffness
        cmap: colormap to use when color_by_stiffness=True
        show_colorbar: whether to show colorbar when color_by_stiffness=True
        ks: array of stiffness values for each bond (required if color_by_stiffness=True)
    """
    if ax is None: _, ax = plt.subplots(figsize=(10, 10))
    
    # Prepare line segments for LineCollection
    segments = []
    for i, j in bonds:
        segments.append([points[i, :2], points[j, :2]])
    
    # Create LineCollection
    lc = LineCollection(segments, linewidths=1)
    
    if color_by_stiffness:
        if ks is None:
            raise ValueError("ks (stiffness values) must be provided when color_by_stiffness=True")
        # Normalize colors by stiffness
        lc.set_array(np.array(ks))
        lc.set_cmap(cmap)
        if show_colorbar:
            plt.colorbar(lc, ax=ax, label='Stiffness')
    else:
        lc.set_color(edge_color)
    
    lc.set_alpha(edge_alpha)
    ax.add_collection(lc)
    
    if point_size is not None:
        ax.scatter(points[:, 0], points[:, 1], c=point_color, s=point_size, zorder=2, label=label)
    
    # Set equal aspect ratio and grid
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    # Update limits
    ax.autoscale()
    
    return ax


class MoleculeTrajectoryVisualizer:
    """
    A class for visualizing molecular trajectories and structures in 2D projections.
    
    This class provides functionality to:
    1. Read and parse XYZ trajectory files
    2. Plot trajectories of selected atoms in xy, xz, and yz projections
    3. Visualize molecular structures at specified intervals
    4. Combine trajectory and structure visualization
    """
    
    def __init__(self):
        """Initialize the MoleculeTrajectoryVisualizer."""
        self.trajectory_data = None
        self.scan_parameters = None
        self.atom_positions = None
        self.atom_types = None
        self.num_frames = 0
        self.num_atoms = 0
        
    def read_xyz_trajectory(self, filepath):
        """
        Read an XYZ trajectory file.
        
        Parameters:
        - filepath: path to the XYZ trajectory file
        
        Returns:
        - self (for method chaining)
        """
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse the first frame to get number of atoms
        num_atoms = int(lines[0].strip())
        self.num_atoms = num_atoms
        
        # Calculate number of frames
        frame_length = num_atoms + 2  # num_atoms + header + comment line
        self.num_frames = len(lines) // frame_length
        
        # Initialize arrays
        self.atom_positions = np.zeros((self.num_frames, num_atoms, 3))
        self.atom_types = []
        self.scan_parameters = np.zeros(self.num_frames)
        
        # Parse each frame
        for frame in range(self.num_frames):
            start_line = frame * frame_length
            
            # Parse comment line for scan parameter (t value)
            comment_line = lines[start_line + 1].strip()
            t_match = re.search(r't\s*=\s*([0-9.]+)', comment_line)
            if t_match:
                self.scan_parameters[frame] = float(t_match.group(1))
            
            # Parse atom positions
            if frame == 0:
                # Only need to parse atom types once
                self.atom_types = []
                for i in range(num_atoms):
                    line = lines[start_line + 2 + i].strip().split()
                    self.atom_types.append(line[0])
                    self.atom_positions[frame, i, :] = [float(x) for x in line[1:4]]
            else:
                for i in range(num_atoms):
                    line = lines[start_line + 2 + i].strip().split()
                    self.atom_positions[frame, i, :] = [float(x) for x in line[1:4]]
        
        return self
    
    def _get_axis_limits(self, positions=None):
        """
        Calculate common axis limits across all dimensions.
        
        Parameters:
        - positions: Optional numpy array of positions. If None, uses all frames.
        
        Returns:
        - min_val: minimum value for all axes
        - max_val: maximum value for all axes
        """
        if positions is None:
            positions = self.atom_positions
        
        # Reshape to 2D array if working with all frames
        if positions.ndim == 3:
            positions = positions.reshape(-1, 3)
        
        # Find global min and max across all dimensions
        min_val = positions.min()
        max_val = positions.max()
        
        # Add small padding (5%)
        padding = 0.05 * (max_val - min_val)
        return min_val - padding, max_val + padding

    def _set_common_limits(self, axes_array, min_val, max_val):
        """
        Set common limits for all axes in the subplot.
        """
        for ax in axes_array:
            ax.set_xlim(min_val, max_val)
            ax.set_ylim(min_val, max_val)
    
    def plot_atom_trajectories(self, atom_indices, projections=None, figsize=(15, 5),
                              colors=None, markers='o', linestyles='-', markersize=4,
                              labels=None, title="Selected Atom Trajectories"):
        """
        Plot trajectories of selected atoms in 2D projections.
        
        Parameters:
        - atom_indices: list of atom indices to plot (0-based)
        - projections: list of projections to show, options: 'xy', 'xz', 'yz' (default: all)
        - figsize: figure size
        - colors: colors for each atom trajectory
        - markers: markers for each atom trajectory
        - linestyles: line styles for each atom trajectory
        - markersize: size of markers
        - labels: labels for each atom trajectory
        - title: global title for the figure
        
        Returns:
        - fig: matplotlib figure
        - axes: dictionary of matplotlib axes for each projection
        """
        if self.atom_positions is None:
            raise ValueError("No trajectory data loaded. Call read_xyz_trajectory first.")
        
        if projections is None:
            projections = ['xy', 'xz', 'yz']
        
        # Set up default colors if not provided
        if colors is None:
            cmap = plt.cm.get_cmap('tab10')
            colors = [cmap(i % 10) for i in range(len(atom_indices))]
        elif isinstance(colors, str):
            colors = [colors] * len(atom_indices)
            
        # Set up default labels if not provided
        if labels is None:
            if self.atom_types:
                labels = [f"{self.atom_types[i]} (Atom {i})" for i in atom_indices]
            else:
                labels = [f"Atom {i}" for i in atom_indices]
        
        # Create figure and axes - 1 row, 3 columns
        n_plots = len(projections)
        fig, axes_array = plt.subplots(1, n_plots, figsize=figsize)
        if n_plots == 1:
            axes_array = [axes_array]
            
        axes = {}
        
        # Plot each projection
        for i, proj in enumerate(projections):
            ax = axes_array[i]
            axes[proj] = ax
            
            # Determine which coordinates to plot
            if proj == 'xy':
                x_idx, y_idx = 0, 1
                xlabel, ylabel = 'X (Å)', 'Y (Å)'
                subtitle = "XY Plane"
            elif proj == 'xz':
                x_idx, y_idx = 0, 2
                xlabel, ylabel = 'X (Å)', 'Z (Å)'
                subtitle = "XZ Plane"
            elif proj == 'yz':
                x_idx, y_idx = 1, 2
                xlabel, ylabel = 'Y (Å)', 'Z (Å)'
                subtitle = "YZ Plane"
            
            # Plot each atom trajectory
            for j, atom_idx in enumerate(atom_indices):
                color = colors[j] if j < len(colors) else 'blue'
                label = labels[j] if j < len(labels) else f"Atom {atom_idx}"
                
                ax.plot(
                    self.atom_positions[:, atom_idx, x_idx],
                    self.atom_positions[:, atom_idx, y_idx],
                    color=color,
                    marker=markers if isinstance(markers, str) else markers[j % len(markers)],
                    linestyle=linestyles if isinstance(linestyles, str) else linestyles[j % len(linestyles)],
                    markersize=markersize,
                    label=label
                )
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(subtitle)
            ax.grid(True)
            ax.legend()
            ax.set_aspect('equal')
        
        # Set common axis limits
        min_val, max_val = self._get_axis_limits()
        self._set_common_limits(axes_array, min_val, max_val)
        
        # Add a global title
        fig.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Make room for the global title
        return fig, axes
    
    def visualize_structures(self, frame_indices=None, projections=None, figsize=(15, 5),
                           bond_length_threshold=2.0, atom_size=10, substrate_indices=None,
                           substrate_atom_size=30, colormap='jet', title="Molecular Structures"):
        """
        Visualize molecular structures at specified frames in 2D projections.
        
        Parameters:
        - frame_indices: indices of frames to visualize (default: evenly spaced frames)
        - projections: list of projections to show, options: 'xy', 'xz', 'yz' (default: all)
        - figsize: figure size
        - bond_length_threshold: maximum distance to consider atoms as bonded (in Angstroms)
        - atom_size: size of atoms in the plot (smaller size for molecule atoms)
        - substrate_indices: list of atom indices that belong to the substrate (will be displayed larger)
        - substrate_atom_size: size of substrate atoms in the plot (larger size for substrate atoms)
        - colormap: colormap for structures at different frames
        - title: global title for the figure
        
        Returns:
        - fig: matplotlib figure
        - axes: dictionary of matplotlib axes for each projection
        """
        if self.atom_positions is None:
            raise ValueError("No trajectory data loaded. Call read_xyz_trajectory first.")
        
        if projections is None:
            projections = ['xy', 'xz', 'yz']
        
        # Determine which frames to visualize
        if frame_indices is None:
            n_structures = min(5, self.num_frames)  # Default to 5 structures or fewer
            frame_indices = np.linspace(0, self.num_frames-1, n_structures, dtype=int)
        
        # Create figure and axes - 1 row, 3 columns
        n_plots = len(projections)
        fig, axes_array = plt.subplots(1, n_plots, figsize=figsize)
        if n_plots == 1:
            axes_array = [axes_array]
            
        axes = {}
        
        # Set up colormap
        cmap = plt.cm.get_cmap(colormap)
        colors = [cmap(i) for i in np.linspace(0, 1, len(frame_indices))]
        
        # Plot each projection
        for i, proj in enumerate(projections):
            ax = axes_array[i]
            axes[proj] = ax
            
            # Determine which coordinates to plot
            if proj == 'xy':
                x_idx, y_idx = 0, 1
                xlabel, ylabel = 'X (Å)', 'Y (Å)'
                subtitle = "XY Plane"
            elif proj == 'xz':
                x_idx, y_idx = 0, 2
                xlabel, ylabel = 'X (Å)', 'Z (Å)'
                subtitle = "XZ Plane"
            elif proj == 'yz':
                x_idx, y_idx = 1, 2
                xlabel, ylabel = 'Y (Å)', 'Z (Å)'
                subtitle = "YZ Plane"
            
            # Plot each structure
            for j, frame_idx in enumerate(frame_indices):
                color = colors[j]
                positions = self.atom_positions[frame_idx]
                
                # Plot atoms with smaller size to reduce clutter
                ax.scatter(
                    positions[:, x_idx],
                    positions[:, y_idx],
                    color=color,
                    s=atom_size,
                    alpha=0.2,
                    label=f"z = {self.scan_parameters[frame_idx]:.2f} Å"
                )
                
                # Connect atoms with bonds
                for a in range(positions.shape[0]):
                    # Calculate distances to all other atoms
                    distances = np.sqrt(np.sum((positions - positions[a])**2, axis=1))
                    distances[a] = np.inf  # Exclude self
                    
                    # Find all neighbors within the threshold distance
                    neighbors = np.where(distances < bond_length_threshold)[0]
                    
                    # Draw bonds to all nearest neighbors
                    for neighbor in neighbors:
                        if a < neighbor:  # Avoid double drawing
                            ax.plot(
                                [positions[a, x_idx], positions[neighbor, x_idx]],
                                [positions[a, y_idx], positions[neighbor, y_idx]],
                                color=color,
                                alpha=0.5
                            )
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(subtitle)
            ax.grid(True)
            ax.legend()
            ax.set_aspect('equal')
        
        # Set common axis limits
        min_val, max_val = self._get_axis_limits()
        self._set_common_limits(axes_array, min_val, max_val)
        
        # Add a global title
        fig.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Make room for the global title
        return fig, axes
    
    def visualize_combined(self, atom_indices, frame_indices=None, projections=None,
                         figsize=(15, 5), bond_length_threshold=2.0, atom_size=10,
                         trajectory_colors=None, structure_colormap='jet',
                         trajectory_markers='o', trajectory_linestyles='-',
                         trajectory_markersize=4, labels=None,
                         title="Combined Atom Trajectories and Molecular Structures"):
        """
        Combine atom trajectory and structure visualization in 2D projections.
        
        Parameters:
        - atom_indices: list of atom indices to plot trajectories for (0-based)
        - frame_indices: indices of frames to visualize structures for (default: evenly spaced)
        - projections: list of projections to show, options: 'xy', 'xz', 'yz' (default: all)
        - figsize: figure size
        - bond_length_threshold: maximum distance to consider atoms as bonded (in Angstroms)
        - atom_size: size of atoms in structure plots (smaller size to reduce clutter)
        - trajectory_colors: colors for atom trajectories
        - structure_colormap: colormap for structures at different frames
        - trajectory_markers: markers for atom trajectories
        - trajectory_linestyles: line styles for atom trajectories
        - trajectory_markersize: size of markers in trajectories
        - labels: labels for atom trajectories
        - title: global title for the figure
        
        Returns:
        - fig: matplotlib figure
        - axes: dictionary of matplotlib axes for each projection
        """
        if self.atom_positions is None:
            raise ValueError("No trajectory data loaded. Call read_xyz_trajectory first.")
        
        if projections is None:
            projections = ['xy', 'xz', 'yz']
        
        # Determine which frames to visualize
        if frame_indices is None:
            n_structures = min(5, self.num_frames)  # Default to 5 structures or fewer
            frame_indices = np.linspace(0, self.num_frames-1, n_structures, dtype=int)
        
        # Create figure and axes - 1 row, 3 columns
        n_plots = len(projections)
        fig, axes_array = plt.subplots(1, n_plots, figsize=figsize)
        if n_plots == 1:
            axes_array = [axes_array]
            
        axes = {}
        
        # Set up colormap for structures
        struct_cmap = plt.cm.get_cmap(structure_colormap)
        struct_colors = [struct_cmap(i) for i in np.linspace(0, 1, len(frame_indices))]
        
        # Set up default colors for trajectories if not provided
        if trajectory_colors is None:
            traj_cmap = plt.cm.get_cmap('tab10')
            trajectory_colors = [traj_cmap(i % 10) for i in range(len(atom_indices))]
        elif isinstance(trajectory_colors, str):
            trajectory_colors = [trajectory_colors] * len(atom_indices)
            
        # Set up default labels if not provided
        if labels is None:
            if self.atom_types:
                labels = [f"{self.atom_types[i]} (Atom {i})" for i in atom_indices]
            else:
                labels = [f"Atom {i}" for i in atom_indices]
        
        # Plot each projection
        for i, proj in enumerate(projections):
            ax = axes_array[i]
            axes[proj] = ax
            
            # Determine which coordinates to plot
            if proj == 'xy':
                x_idx, y_idx = 0, 1
                xlabel, ylabel = 'X (Å)', 'Y (Å)'
                subtitle = "XY Plane"
            elif proj == 'xz':
                x_idx, y_idx = 0, 2
                xlabel, ylabel = 'X (Å)', 'Z (Å)'
                subtitle = "XZ Plane"
            elif proj == 'yz':
                x_idx, y_idx = 1, 2
                xlabel, ylabel = 'Y (Å)', 'Z (Å)'
                subtitle = "YZ Plane"
            
            # First plot structures with smaller atom size to reduce clutter
            for j, frame_idx in enumerate(frame_indices):
                color = struct_colors[j]
                positions = self.atom_positions[frame_idx]
                
                # Plot atoms with smaller size
                ax.scatter(
                    positions[:, x_idx],
                    positions[:, y_idx],
                    color=color,
                    s=atom_size,
                    alpha=0.2,
                    label=f"Structure: z = {self.scan_parameters[frame_idx]:.2f} Å"
                )
                
                # Connect atoms with bonds
                for a in range(positions.shape[0]):
                    # Calculate distances to all other atoms
                    distances = np.sqrt(np.sum((positions - positions[a])**2, axis=1))
                    distances[a] = np.inf  # Exclude self
                    
                    # Find all neighbors within the threshold distance
                    neighbors = np.where(distances < bond_length_threshold)[0]
                    
                    # Draw bonds to all nearest neighbors
                    for neighbor in neighbors:
                        if a < neighbor:  # Avoid double drawing
                            ax.plot(
                                [positions[a, x_idx], positions[neighbor, x_idx]],
                                [positions[a, y_idx], positions[neighbor, y_idx]],
                                color=color,
                                alpha=0.5
                            )
            
            # Then plot atom trajectories
            for j, atom_idx in enumerate(atom_indices):
                color = trajectory_colors[j] if j < len(trajectory_colors) else 'blue'
                label = labels[j] if j < len(labels) else f"Atom {atom_idx}"
                
                ax.plot(
                    self.atom_positions[:, atom_idx, x_idx],
                    self.atom_positions[:, atom_idx, y_idx],
                    color=color,
                    marker=trajectory_markers if isinstance(trajectory_markers, str) else trajectory_markers[j % len(trajectory_markers)],
                    linestyle=" " if isinstance(trajectory_linestyles, str) else trajectory_linestyles[j % len(trajectory_linestyles)],
                    markersize=trajectory_markersize-3,
                    label=f"Trajectory: {label}"
                )
            
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(subtitle)
            ax.grid(True)
            ax.legend()
            ax.set_aspect('equal')
        
        # Set common axis limits
        min_val, max_val = self._get_axis_limits()
        self._set_common_limits(axes_array, min_val, max_val)
        
        # Add a global title
        fig.suptitle(title, fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.95])  # Make room for the global title
        return fig, axes
    
    def save_visualization(self, fig, filename, dpi=300):
        """
        Save the visualization to a file.
        
        Parameters:
        - fig: matplotlib figure to save
        - filename: output filename
        - dpi: resolution in dots per inch
        """
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Visualization saved to {filename}")


def compare_qgrids(qgrid_cpu, qgrid_gpu, grid_step_cpu=None, grid_step_gpu=None, title_prefix=""):
    """
    Compares CPU and GPU Qgrid arrays (charge density).

    Args:
        qgrid_cpu (np.ndarray): Charge density grid from CPU (expected ZYX order: nz, ny, nx).
        qgrid_gpu (np.ndarray): Charge density grid from GPU (expected ZYX order: nz, ny, nx).
        grid_step_cpu (float, optional): Grid spacing for CPU grid for plot labels.
        grid_step_gpu (float, optional): Grid spacing for GPU grid for plot labels.
        title_prefix (str, optional): Prefix for plot titles.
    """
    print("--- Qgrid Comparison ---")

    # 1. Shape Check
    print(f"CPU Qgrid shape: {qgrid_cpu.shape}")
    print(f"GPU Qgrid shape: {qgrid_gpu.shape}")

    if qgrid_cpu.shape != qgrid_gpu.shape:
        print("ERROR: Qgrid shapes do not match! Cannot perform detailed comparison.")
        # Attempt to give some basic stats if dimensions are different but plausible
        print(f"CPU Qgrid stats: min={qgrid_cpu.min():.4e}, max={qgrid_cpu.max():.4e}, sum_abs={np.sum(np.abs(qgrid_cpu)):.4e}")
        print(f"GPU Qgrid stats: min={qgrid_gpu.min():.4e}, max={qgrid_gpu.max():.4e}, sum_abs={np.sum(np.abs(qgrid_gpu)):.4e}")
        return

    nz, ny, nx = qgrid_cpu.shape
    print(f"Grid dimensions (nz, ny, nx): ({nz}, {ny}, {nx})")

    # 2. Numerical Comparison
    diff = qgrid_gpu - qgrid_cpu
    abs_diff_sum = np.sum(np.abs(diff))
    max_abs_diff = np.max(np.abs(diff))
    mean_abs_diff = np.mean(np.abs(diff))

    print("\nNumerical Stats:")
    print(f"  CPU Qgrid: min={qgrid_cpu.min():.4e}, max={qgrid_cpu.max():.4e}, mean={qgrid_cpu.mean():.4e}, std={qgrid_cpu.std():.4e}, sum_abs={np.sum(np.abs(qgrid_cpu)):.4e}")
    print(f"  GPU Qgrid: min={qgrid_gpu.min():.4e}, max={qgrid_gpu.max():.4e}, mean={qgrid_gpu.mean():.4e}, std={qgrid_gpu.std():.4e}, sum_abs={np.sum(np.abs(qgrid_gpu)):.4e}")
    print(f"  Difference (GPU - CPU):")
    print(f"    Min Diff: {diff.min():.4e}")
    print(f"    Max Diff: {diff.max():.4e}")
    print(f"    Mean Abs Diff: {mean_abs_diff:.4e}")
    print(f"    Max Abs Diff: {max_abs_diff:.4e}")
    print(f"    Sum of Abs Diff: {abs_diff_sum:.4e}")

    if abs_diff_sum < 1e-5: # Adjust tolerance as needed
        print("Qgrids are numerically very similar.")
    else:
        print("Qgrids show notable numerical differences.")

    # 3. Visual Comparison (Slices)
    mid_z, mid_y, mid_x = nz // 2, ny // 2, nx // 2
    
    # Determine common color limits for CPU and GPU plots for better comparison
    vmin = min(qgrid_cpu.min(), qgrid_gpu.min())
    vmax = max(qgrid_cpu.max(), qgrid_gpu.max())

    # Determine color limits for difference plot (centered at 0)
    diff_abs_max = np.max(np.abs(diff))
    if diff_abs_max < 1e-9: # Avoid issues if diff is all zeros
        diff_abs_max = 1e-9
    diff_vmin = -diff_abs_max
    diff_vmax = diff_abs_max

    plot_slices = [
        ('XY', mid_z, (slice(None), slice(None), mid_z), (0, 1), ("Y", "X") if grid_step_cpu else ("Y index", "X index")), # GPU Qgrid is (nz,ny,nx), so Q[mid_z,:,:]
        ('XZ', mid_y, (slice(None), mid_y, slice(None)), (0, 2), ("Z", "X") if grid_step_cpu else ("Z index", "X index")), # Q[:,mid_y,:]
        ('YZ', mid_x, (mid_x, slice(None), slice(None)), (1, 2), ("Z", "Y") if grid_step_cpu else ("Z index", "Y index"))  # Q[:,:,mid_x] - careful, this is ZY if imshow is used directly
    ]
    
    # Adjusting for imshow behavior with ZYX data
    # For XY: data[mid_z, :, :] -> imshow(data[mid_z, :, :]) -> Y is rows, X is columns
    # For XZ: data[:, mid_y, :] -> imshow(data[:, mid_y, :]) -> Z is rows, X is columns
    # For YZ: data[:, :, mid_x] -> imshow(data[:, :, mid_x]) -> Z is rows, Y is columns

    for plane_name, slice_idx, slice_obj_zyx, (dim1_idx, dim2_idx), (label1, label2) in plot_slices:
        
        if plane_name == 'XY':
            cpu_slice = qgrid_cpu[slice_idx, :, :]
            gpu_slice = qgrid_gpu[slice_idx, :, :]
            diff_slice = diff[slice_idx, :, :]
            aspect_ratio = (ny * (grid_step_cpu if grid_step_cpu else 1)) / (nx * (grid_step_cpu if grid_step_cpu else 1)) if grid_step_cpu else 'auto'
            extent = [0, nx*(grid_step_cpu if grid_step_cpu else 1), 0, ny*(grid_step_cpu if grid_step_cpu else 1)] if grid_step_cpu else None
            xlabel, ylabel = label2, label1 # X, Y
        elif plane_name == 'XZ':
            cpu_slice = qgrid_cpu[:, slice_idx, :]
            gpu_slice = qgrid_gpu[:, slice_idx, :]
            diff_slice = diff[:, slice_idx, :]
            aspect_ratio = (nz * (grid_step_cpu if grid_step_cpu else 1)) / (nx * (grid_step_cpu if grid_step_cpu else 1)) if grid_step_cpu else 'auto'
            extent = [0, nx*(grid_step_cpu if grid_step_cpu else 1), 0, nz*(grid_step_cpu if grid_step_cpu else 1)] if grid_step_cpu else None
            xlabel, ylabel = label2, label1 # X, Z
        elif plane_name == 'YZ': # Z vs Y
            cpu_slice = qgrid_cpu[:, :, slice_idx]
            gpu_slice = qgrid_gpu[:, :, slice_idx]
            diff_slice = diff[:, :, slice_idx]
            aspect_ratio = (nz * (grid_step_cpu if grid_step_cpu else 1)) / (ny * (grid_step_cpu if grid_step_cpu else 1)) if grid_step_cpu else 'auto'
            extent = [0, ny*(grid_step_cpu if grid_step_cpu else 1), 0, nz*(grid_step_cpu if grid_step_cpu else 1)] if grid_step_cpu else None
            xlabel, ylabel = label2, label1 # Y, Z


        fig, axs = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f'{title_prefix}Qgrid Comparison - {plane_name} Plane (Slice at {"zyx"[dim1_idx if plane_name != "YZ" else dim2_idx]}={slice_idx})', fontsize=16)

        im_cpu = axs[0].imshow(cpu_slice, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax, extent=extent, aspect=aspect_ratio)
        axs[0].set_title('CPU Qgrid')
        axs[0].set_xlabel(xlabel)
        axs[0].set_ylabel(ylabel)
        divider = make_axes_locatable(axs[0])
        cax_cpu = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im_cpu, cax=cax_cpu, label='Charge Density')

        im_gpu = axs[1].imshow(gpu_slice, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax, extent=extent, aspect=aspect_ratio)
        axs[1].set_title('GPU Qgrid')
        axs[1].set_xlabel(xlabel)
        axs[1].set_ylabel(ylabel)
        divider = make_axes_locatable(axs[1])
        cax_gpu = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im_gpu, cax=cax_gpu, label='Charge Density')

        im_diff = axs[2].imshow(diff_slice, origin='lower', cmap='coolwarm', vmin=diff_vmin, vmax=diff_vmax, extent=extent, aspect=aspect_ratio)
        axs[2].set_title('Difference (GPU - CPU)')
        axs[2].set_xlabel(xlabel)
        axs[2].set_ylabel(ylabel)
        divider = make_axes_locatable(axs[2])
        cax_diff = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im_diff, cax=cax_diff, label='Difference')

        plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust for suptitle
        plt.show()

# --- Example Usage (assuming you have qgrid_cpu and qgrid_gpu loaded) ---
# # compare_qgrids(qgrid_cpu, qgrid_gpu, grid_step_cpu=dg_val, grid_step_gpu=dg_val, title_prefix="Initial ")


def print_qgrid_boundary_slices(qgrid_cpu, qgrid_gpu, num_edge_slices=5, num_line_elements=5, title_prefix=""):
    """
    Prints numerical values from the boundary slices of CPU and GPU Qgrids
    for detailed comparison, focusing on potential 1-pixel shifts.

    Args:
        qgrid_cpu (np.ndarray): Charge density grid from CPU (expected ZYX order: nz, ny, nx).
        qgrid_gpu (np.ndarray): Charge density grid from GPU (expected ZYX order: nz, ny, nx).
        num_edge_slices (int): Number of slices to print from each edge (e.g., Z=0,1,2,... and Z=nz-1,nz-2,...).
        num_line_elements (int): Number of elements to print along a line for XZ and YZ planes.
        title_prefix (str): Prefix for print section titles.
    """
    if qgrid_cpu.shape != qgrid_gpu.shape:
        print(f"{title_prefix}ERROR: Qgrid shapes do not match! CPU: {qgrid_cpu.shape}, GPU: {qgrid_gpu.shape}")
        return

    nz, ny, nx = qgrid_cpu.shape
    mid_y, mid_x = ny // 2, nx // 2
    
    # Ensure num_line_elements is odd for centering, and calculate half-width
    if num_line_elements % 2 == 0:
        num_line_elements +=1 
    half_line = num_line_elements // 2

    print(f"\n--- {title_prefix}Numerical Boundary Slice Comparison (CPU vs GPU) ---")
    print(f"Grid dimensions (nz, ny, nx): ({nz}, {ny}, {nx})")
    print(f"Comparing {num_edge_slices} slices from Z-edges.")
    print(f"For XZ/YZ planes, showing {num_line_elements} elements around center line.")
    print(f"Values formatted to 6 decimal places (scientific notation if small).\n")

    # --- XZ Plane (fixed mid_y, line along X) ---
    print(f"--- {title_prefix}XZ Plane (fixed Y={mid_y}, X indices [{mid_x-half_line}:{mid_x+half_line+1}]) ---")
    x_slice = slice(mid_x - half_line, mid_x + half_line + 1)
    for i in range(num_edge_slices):
        z_indices = [i, nz - 1 - i]
        for z_idx in z_indices:
            if z_idx < 0 or z_idx >= nz: continue # Should not happen with typical num_edge_slices
            print(f" Z = {z_idx}:")
            cpu_vals = qgrid_cpu[z_idx, mid_y, x_slice]
            gpu_vals = qgrid_gpu[z_idx, mid_y, x_slice]
            print(f"   CPU: " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in cpu_vals]))
            print(f"   GPU: " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in gpu_vals]))
            if i == num_edge_slices -1 and z_idx == nz -1 -i : print("-" * 20) # Separator after last bottom slice
        if i == num_edge_slices -1 : print("\n")


    # --- YZ Plane (fixed mid_x, line along Y) ---
    print(f"--- {title_prefix}YZ Plane (fixed X={mid_x}, Y indices [{mid_y-half_line}:{mid_y+half_line+1}]) ---")
    y_slice = slice(mid_y - half_line, mid_y + half_line + 1)
    for i in range(num_edge_slices):
        z_indices = [i, nz - 1 - i]
        for z_idx in z_indices:
            if z_idx < 0 or z_idx >= nz: continue
            print(f" Z = {z_idx}:")
            cpu_vals = qgrid_cpu[z_idx, y_slice, mid_x]
            gpu_vals = qgrid_gpu[z_idx, y_slice, mid_x]
            print(f"   CPU: " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in cpu_vals]))
            print(f"   GPU: " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in gpu_vals]))
            if i == num_edge_slices -1 and z_idx == nz -1 -i : print("-" * 20)
        if i == num_edge_slices -1 : print("\n")

    # --- XY Plane (fixed Z, small patch) ---
    print(f"--- {title_prefix}XY Plane (Z-slices, patch around Y={mid_y}, X={mid_x}) ---")
    y_patch_slice = slice(mid_y - half_line, mid_y + half_line + 1)
    x_patch_slice = slice(mid_x - half_line, mid_x + half_line + 1)
    
    for i in range(num_edge_slices):
        z_indices = [i, nz - 1 - i]
        for z_idx in z_indices:
            if z_idx < 0 or z_idx >= nz: continue
            print(f" Z = {z_idx} (Y indices [{mid_y-half_line}:{mid_y+half_line+1}], X indices [{mid_x-half_line}:{mid_x+half_line+1}]):")
            cpu_patch = qgrid_cpu[z_idx, y_patch_slice, x_patch_slice]
            gpu_patch = qgrid_gpu[z_idx, y_patch_slice, x_patch_slice]
            
            print("  CPU Patch:")
            for row_idx in range(cpu_patch.shape[0]):
                print("    " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in cpu_patch[row_idx, :]]))
            print("  GPU Patch:")
            for row_idx in range(gpu_patch.shape[0]):
                print("    " + " ".join([f"{v:10.6f}" if abs(v) > 1e-4 or v==0 else f"{v:10.3e}" for v in gpu_patch[row_idx, :]]))
            print("-" * 10)
            if i == num_edge_slices -1 and z_idx == nz -1 -i : print("-" * 20)
        if i == num_edge_slices -1 : print("\n")

    print(f"--- {title_prefix}End of Numerical Boundary Slice Comparison ---\n")

# --- Example of how to use it in your test script ---
# print_qgrid_boundary_slices(qgrid_cpu_test, qgrid_gpu_test, num_edge_slices=3, num_line_elements=5, title_prefix="Test ")

