# Molecule Trajectory Visualization Guide

This guide explains how to use the `MoleculeTrajectoryVisualizer` class to visualize molecular trajectories and structures in 2D projections.

## Overview

The `MoleculeTrajectoryVisualizer` class provides functionality to:

1. Read and parse XYZ trajectory files
2. Plot trajectories of selected atoms in xy, xz, and yz projections
3. Visualize molecular structures at specified intervals
4. Combine trajectory and structure visualization

## Installation

The `MoleculeTrajectoryVisualizer` class is part of the `plot_utils.py` module. Make sure you have the following dependencies installed:

- NumPy
- Matplotlib

## Usage

### Basic Usage

```python
from plot_utils import MoleculeTrajectoryVisualizer
import matplotlib.pyplot as plt

# Create a visualizer instance
visualizer = MoleculeTrajectoryVisualizer()

# Read a trajectory file
visualizer.read_xyz_trajectory("path/to/trajectory.xyz")

# Plot trajectories of specific atoms
atom_indices = [26, 29]  # 0-based indexing
fig, axes = visualizer.plot_atom_trajectories(atom_indices)
plt.show()
```

### Reading Trajectory Files

The class can read XYZ trajectory files with the following format:

```
number_of_atoms
t = value, E = value eV
atom_type x y z
atom_type x y z
...
number_of_atoms
t = value, E = value eV
atom_type x y z
atom_type x y z
...
```

Example:

```
422
t = 1.600, E = 46794960.618817 eV
C     2.619100     2.814374     8.869142
C     8.222985     2.468584     9.127351
...
```

To read a trajectory file:

```python
visualizer.read_xyz_trajectory("path/to/trajectory.xyz")
```

### Plotting Atom Trajectories

To plot the trajectories of specific atoms in 2D projections:

```python
# Plot trajectories of atoms 26 and 29
atom_indices = [26, 29]  # 0-based indexing
fig, axes = visualizer.plot_atom_trajectories(
    atom_indices=atom_indices,
    projections=['xy', 'xz', 'yz'],  # Show all three 2D projections in one row
    figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
    colors=['red', 'blue'],
    markers='o',
    linestyles='-',
    markersize=4,
    title="Selected Atom Trajectories"  # Global title for the figure
)
plt.show()
```

#### Parameters:

- `atom_indices`: List of atom indices to plot (0-based)
- `projections`: List of projections to show, options: 'xy', 'xz', 'yz' (default: all)
- `figsize`: Figure size (default: (15, 5) for 1 row, 3 columns layout)
- `colors`: Colors for each atom trajectory
- `markers`: Markers for each atom trajectory
- `linestyles`: Line styles for each atom trajectory
- `markersize`: Size of markers
- `labels`: Labels for each atom trajectory
- `title`: Global title for the figure

### Visualizing Molecular Structures

To visualize molecular structures at specified frames:

```python
import numpy as np

# Select 5 evenly spaced frames from the trajectory
n_structures = 5
frame_indices = np.linspace(0, visualizer.num_frames-1, n_structures, dtype=int)

fig, axes = visualizer.visualize_structures(
    frame_indices=frame_indices,
    projections=['xy', 'xz', 'yz'],
    figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
    bond_length_threshold=2.0,  # Maximum distance to consider atoms as bonded (in Angstroms)
    atom_size=10,  # Smaller atom size to reduce clutter
    colormap='jet',
    title="Molecular Structures"  # Global title for the figure
)
plt.show()
```

#### Parameters:

- `frame_indices`: Indices of frames to visualize (default: evenly spaced frames)
- `projections`: List of projections to show, options: 'xy', 'xz', 'yz' (default: all)
- `figsize`: Figure size (default: (15, 5) for 1 row, 3 columns layout)
- `bond_length_threshold`: Maximum distance to consider atoms as bonded (in Angstroms)
- `atom_size`: Size of atoms in the plot (smaller size to reduce clutter)
- `colormap`: Colormap for structures at different frames
- `title`: Global title for the figure

### Combined Visualization

To combine atom trajectory and structure visualization:

```python
fig, axes = visualizer.visualize_combined(
    atom_indices=[26, 29],
    frame_indices=frame_indices,
    projections=['xy', 'xz', 'yz'],
    figsize=(15, 5),  # Wider figure for 1 row, 3 columns layout
    bond_length_threshold=2.0,
    atom_size=10,  # Smaller atom size to reduce clutter
    trajectory_colors=['red', 'blue'],
    structure_colormap='jet',
    trajectory_markers='o',
    trajectory_linestyles='-',
    trajectory_markersize=4,
    title="Combined Atom Trajectories and Molecular Structures"  # Global title for the figure
)
plt.show()
```

#### Parameters:

- `atom_indices`: List of atom indices to plot trajectories for (0-based)
- `frame_indices`: Indices of frames to visualize structures for (default: evenly spaced)
- `projections`: List of projections to show, options: 'xy', 'xz', 'yz' (default: all)
- `figsize`: Figure size (default: (15, 5) for 1 row, 3 columns layout)
- `bond_length_threshold`: Maximum distance to consider atoms as bonded (in Angstroms)
- `atom_size`: Size of atoms in structure plots (smaller size to reduce clutter)
- `trajectory_colors`: Colors for atom trajectories
- `structure_colormap`: Colormap for structures at different frames
- `trajectory_markers`: Markers for atom trajectories
- `trajectory_linestyles`: Line styles for atom trajectories
- `trajectory_markersize`: Size of markers in trajectories
- `labels`: Labels for atom trajectories
- `title`: Global title for the figure

### Saving Visualizations

To save visualizations to files:

```python
# Save the figure
fig.savefig("visualization.png", dpi=300, bbox_inches='tight')

# Or use the convenience method
visualizer.save_visualization(fig, "visualization.png", dpi=300)
```

## Example Script

See the `molecule_visualization_example.py` file for a complete example of how to use the `MoleculeTrajectoryVisualizer` class.

## Customization

The visualizations can be customized using standard Matplotlib functions after obtaining the figure and axes:

```python
fig, axes = visualizer.plot_atom_trajectories(atom_indices)

# Customize the xy projection
axes['xy'].set_title("Custom Title")
axes['xy'].set_xlim(-10, 10)
axes['xy'].set_ylim(-10, 10)

plt.show()