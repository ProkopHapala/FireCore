# Additional Visualization Script Improvements

## Overview
I've implemented three additional improvements to the `visualize_top_layer_xy.py` script based on your latest requirements:

## 1. ✅ Extra Substrate Layers Only for XZ and YZ Projections

**Implementation**: The logic was already correctly implemented in the previous version:
- **XY projection**: Shows only the top substrate layer (single layer)
- **XZ and YZ projections**: Show 3 substrate layers to reveal the 3D stacking structure

**Code Logic**:
```python
# Use appropriate substrate layers based on projection
if proj == "xy":
    # For XY projection, use only top layer
    substrate_layer = _top_layer_indices(
        vis.atom_positions[0], substrate_indices, tol=top_layer_tol
    )
else:
    # For XZ and YZ projections, use multiple layers to show 3D structure
    substrate_layer = _multi_layer_indices(
        vis.atom_positions[0], substrate_indices, num_layers=3, tol=top_layer_tol
    )
```

## 2. ✅ Custom Frame Selection for Polygons and Molecule Snapshots

**New Features**:
- `--custom-frames`: Specify exact frame indices for polygon/convex hull display
- `--custom-mol-frames`: Specify exact frame indices for molecule snapshots

**New Command Line Arguments**:
```bash
--custom-frames "0,10,20,30"        # Custom frames for polygons
--custom-mol-frames "5,15"          # Custom frames for molecule snapshots
```

**Implementation Details**:
- Custom frames override the automatic spacing (`--samples` and `--num-mol-snapshots`)
- Frame indices are validated to ensure they're within the trajectory range
- Both parameters are optional - if not provided, the script uses the original behavior

**Usage Examples**:

### Using Custom Polygon Frames
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --custom-frames "0,5,10,15,20" \
    --projections xy,xz,yz \
    --out custom_polygons.png
```

### Using Custom Molecule Snapshot Frames
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --custom-mol-frames "3,7,12" \
    --projections xy,xz,yz \
    --out custom_molecules.png
```

### Using Both Custom Frame Types
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --custom-frames "0,8,16,24" \
    --custom-mol-frames "4,12,20" \
    --projections xy,xz,yz \
    --out custom_both.png
```

## 3. ✅ Square Subplots with Consistent Dimensions

**Problem Solved**: Previously, xy plots appeared square while xz and yz plots appeared rectangular due to different axis ranges.

**Solution**: Enhanced the `_calculate_common_limits()` function to ensure all subplots have the same aspect ratio:

**Key Changes**:
```python
def _calculate_common_limits(vis, substrate_indices, molecule_indices):
    # ... calculate individual ranges ...
    
    # Make all ranges equal to the maximum range for square subplots
    max_range = max(ranges.values())
    
    for coord in ['x', 'y', 'z']:
        current_center = (limits[coord][0] + limits[coord][1]) / 2
        half_max_range = max_range / 2
        limits[coord] = (current_center - half_max_range, current_center + half_max_range)
    
    return limits
```

**Benefits**:
- All subplots (xy, xz, yz) now have identical dimensions
- Square aspect ratio maintained for all projections
- Data is centered within each subplot
- Consistent visual appearance across all projections

## Complete Feature Summary

### Original Features (from previous improvements):
1. ✅ Polygon display control with `--samples=0`
2. ✅ Improved bond detection using 3D distances
3. ✅ Multiple substrate layers for xz/yz projections

### New Features (this update):
4. ✅ Custom frame selection for polygons (`--custom-frames`)
5. ✅ Custom frame selection for molecule snapshots (`--custom-mol-frames`)
6. ✅ Square subplots with consistent dimensions

## Testing Results

All improvements have been tested successfully:

### Test 1: Square Subplots
```bash
python visualize_top_layer_xy.py \
    --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --projections xy,xz,yz \
    --out test_square_subplots.png
```
✅ **Result**: All subplots now have identical square dimensions

### Test 2: Custom Frame Selection
```bash
python visualize_top_layer_xy.py \
    --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --custom-frames "0,5,10,15" \
    --custom-mol-frames "2,8" \
    --projections xy,xz,yz \
    --out test_custom_frames.png
```
✅ **Result**: Custom frame selection works correctly for both polygons and molecule snapshots

## Backward Compatibility

All changes maintain full backward compatibility:
- Existing commands continue to work exactly as before
- New parameters are optional
- Default behavior is unchanged when new parameters are not used

## Error Handling

Added robust error handling for custom frame indices:
- Validates that frame indices are within the trajectory range
- Provides clear error messages for invalid indices
- Gracefully handles empty or malformed input
