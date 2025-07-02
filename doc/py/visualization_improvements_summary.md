# Visualization Script Improvements Summary

## Overview
I've implemented several improvements to the `visualize_top_layer_xy.py` script based on your requirements. Here's a detailed summary of all the changes made:

## 1. Fixed Polygon Display for --samples=0

**Problem**: Previously, polygons were always shown regardless of the `--samples` parameter value.

**Solution**: 
- Modified the convex hull drawing logic to only display polygons when `n_sample_structures > 0`
- When `--samples 0` is specified, no polygons are drawn, only the atom trajectories and molecular snapshots are shown

**Code Changes**:
```python
# Only draw polygons if n_sample_structures > 0
if n_sample_structures > 0:
    for fr in frame_sel:
        # ... polygon drawing code ...
```

## 2. Consistent Subplot Dimensions

**Problem**: Different projections (xy, xz, yz) had different axis limits, making them look inconsistent.

**Solution**:
- Added `_calculate_common_limits()` function that calculates consistent axis limits across all dimensions
- Applied the same limits to all subplots to ensure uniform appearance
- Added padding around the data for better visualization

**Code Changes**:
```python
def _calculate_common_limits(vis, substrate_indices, molecule_indices):
    # Calculate common axis limits for consistent subplot dimensions
    all_indices = np.concatenate([substrate_indices, molecule_indices])
    all_positions = vis.atom_positions[:, all_indices, :]
    
    padding = 2.0  # Angstroms
    limits = {}
    for i, coord in enumerate(['x', 'y', 'z']):
        coord_min = np.min(all_positions[:, :, i]) - padding
        coord_max = np.max(all_positions[:, :, i]) + padding
        limits[coord] = (coord_min, coord_max)
    
    return limits
```

## 3. Multiple Substrate Layers for XZ and YZ Projections

**Problem**: XZ and YZ projections only showed the top layer, missing the 3D structure of the substrate.

**Solution**:
- Added `_multi_layer_indices()` function to select multiple substrate layers
- For XY projection: Uses only the top layer (as before)
- For XZ and YZ projections: Uses 3 layers to show the 3D stacking structure

**Code Changes**:
```python
def _multi_layer_indices(positions, indices, num_layers=3, tol=0.35):
    """Return subset of indices that belong to the top num_layers z-layers."""
    zs = positions[indices, 2]
    z_max = zs.max()
    layer_threshold = z_max - (num_layers * tol)
    return indices[zs >= layer_threshold]

# In the plotting loop:
if proj == "xy":
    # For XY projection, use only top layer
    substrate_layer = _top_layer_indices(...)
else:
    # For XZ and YZ projections, use multiple layers
    substrate_layer = _multi_layer_indices(..., num_layers=3)
```

## 4. Improved Bond Detection for 2D Projections

**Problem**: In XZ and YZ projections, atoms that appeared close in 2D but were far apart in 3D were incorrectly bonded.

**Solution**:
- Modified bond detection to use 3D distances for determining bonds
- Bonds are still drawn in 2D, but only if atoms are actually close in 3D space
- This prevents false bonds that appear due to projection effects

**Code Changes**:
```python
# bonds in 2-D projection - use 3D distances but plot in 2D
# This prevents false bonds that appear close in 2D but are far in 3D
for i in range(pos_3d.shape[0]):
    # Calculate 3D distances to avoid false bonds in projections
    d2_3d = np.sum((pos_3d - pos_3d[i]) ** 2, axis=1)
    neigh = np.where((d2_3d < bond_length_thresh**2) & (d2_3d > 0))[0]
    for j in neigh:
        if i < j:
            ax.plot(
                (pos_2d[i, 0], pos_2d[j, 0]),
                (pos_2d[i, 1], pos_2d[j, 1]),
                # ... plotting parameters ...
            )
```

## Usage Examples

### No Polygons (samples=0)
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --projections xy,xz,yz \
    --out output.png
```

### With Polygons (samples>0)
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 4 \
    --projections xy,xz,yz \
    --out output.png
```

## Benefits

1. **Better Control**: You can now choose whether to show polygons or not
2. **Consistent Appearance**: All subplots have the same dimensions and scale
3. **Better 3D Representation**: XZ and YZ projections now show the layered structure of the substrate
4. **Accurate Bonding**: No more false bonds in 2D projections
5. **Cleaner Code**: Removed unused imports and variables

## Testing

The improvements have been tested with your original command:
```bash
python visualize_top_layer_xy.py \
    --traj /home/indranil/git/FireCore/tests/tMMFF/PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --num-mol-snapshots 2 \
    --out 1d_z_trajectory_improved.png \
    --projections xy,xz,yz
```

Both with `--samples 0` (no polygons) and `--samples 4` (with polygons) work correctly.
