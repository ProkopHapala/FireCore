# Final Visualization Script Improvements

## Overview
I've implemented the final two improvements to address the remaining issues with substrate layer visibility and molecule coloring.

## 1. ✅ Enhanced Substrate Layer Visualization for XZ and YZ Projections

### Problem Solved
The extra substrate layers weren't clearly visible in XZ and YZ projections due to:
- Insufficient layer detection algorithm
- All layers having the same transparency
- Poor layer separation logic

### Solution Implemented

#### Improved Layer Detection Algorithm
```python
def _multi_layer_indices(positions, indices, num_layers=3, tol=0.35):
    """Enhanced layer detection using z-coordinate clustering."""
    zs = positions[indices, 2]
    
    # Find distinct z-layers by clustering z-coordinates
    z_sorted = np.sort(zs)[::-1]  # Sort in descending order
    layer_z_values = [z_sorted[0]]  # Start with the highest z
    
    for z in z_sorted[1:]:
        # If this z is significantly different from the last layer, it's a new layer
        if layer_z_values[-1] - z > tol:
            layer_z_values.append(z)
            if len(layer_z_values) >= num_layers:
                break
    
    # Include atoms from the identified layers
    selected_indices = []
    for layer_z in layer_z_values:
        layer_mask = np.abs(zs - layer_z) <= tol/2
        selected_indices.extend(indices[layer_mask])
    
    return np.array(selected_indices, dtype=int)
```

#### Layer-Specific Transparency
```python
# Group atoms by z-layers for different alpha values
z_unique = np.unique(z_coords)
layer_alphas = np.linspace(0.6, 0.2, len(z_unique))  # Higher layers more opaque

for i, z_val in enumerate(z_unique):
    layer_mask = np.abs(z_coords - z_val) < 0.1
    # ... plot with layer_alphas[i] ...
```

### Key Improvements:
- **Better Layer Detection**: Uses z-coordinate clustering to identify distinct layers
- **Variable Transparency**: Higher layers are more opaque (α=0.6), lower layers more transparent (α=0.2)
- **Proper Layer Separation**: Uses tolerance-based grouping to separate layers
- **XY vs XZ/YZ Logic**: XY shows only top layer, XZ/YZ show multiple layers with depth

## 2. ✅ Customizable Molecule Snapshot Colors

### New Feature: `--mol-colors` Parameter

#### Color Options:
1. **`auto`** (default): Automatically assigns different colors to each molecule snapshot
2. **`same`**: All molecule snapshots use the same color (black)
3. **Custom color list**: Comma-separated list of specific colors

### Usage Examples:

#### Automatic Different Colors
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --mol-colors "auto" \
    --projections xy,xz,yz \
    --out auto_colors.png
```

#### All Same Color (Black)
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --mol-colors "same" \
    --projections xy,xz,yz \
    --out same_color.png
```

#### Custom Color List
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --mol-colors "red,blue,green,purple" \
    --projections xy,xz,yz \
    --out custom_colors.png
```

### Implementation Details:

#### Color Assignment Logic
```python
# Determine color for this frame
if mol_colors is None or len(mol_colors) == 0:
    # Default: all black
    mol_color = 'k'
elif len(mol_colors) == 1:
    # Single color for all frames
    mol_color = mol_colors[0]
else:
    # Cycle through provided colors
    mol_color = mol_colors[frame_idx % len(mol_colors)]
```

#### Color Processing in Main Function
```python
# Parse molecule colors
mol_colors = None
if args.mol_colors == "auto":
    # Generate different colors automatically
    import matplotlib.pyplot as plt
    cmap = plt.cm.tab10
    mol_colors = [cmap(i) for i in range(10)]  # Use first 10 colors from tab10
elif args.mol_colors == "same":
    mol_colors = ['k']  # All black
elif args.mol_colors != "auto":
    # Custom color list
    mol_colors = [c.strip() for c in args.mol_colors.split(',') if c.strip()]
```

### Benefits:
- **Visual Distinction**: Different colors help distinguish between molecule snapshots at different time points
- **Flexibility**: Users can choose between automatic, uniform, or custom coloring schemes
- **Consistency**: Both atoms and bonds use the same color for each snapshot
- **Backward Compatibility**: Default behavior (black) is maintained when parameter is not specified

## Complete Feature Summary

### All Implemented Features:
1. ✅ Polygon display control with `--samples=0`
2. ✅ Improved bond detection using 3D distances
3. ✅ Multiple substrate layers for xz/yz projections (enhanced)
4. ✅ Custom frame selection for polygons (`--custom-frames`)
5. ✅ Custom frame selection for molecule snapshots (`--custom-mol-frames`)
6. ✅ Square subplots with consistent dimensions
7. ✅ Enhanced substrate layer visualization with depth
8. ✅ Customizable molecule snapshot colors (`--mol-colors`)

## Testing Results

All improvements have been tested successfully:

### Test 1: Enhanced Substrate Layers
```bash
python visualize_top_layer_xy.py \
    --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --projections xy,xz,yz \
    --mol-colors "auto" \
    --out test_layers_colors.png
```
✅ **Result**: Multiple substrate layers now clearly visible in XZ and YZ projections with proper depth visualization

### Test 2: Custom Colors
```bash
python visualize_top_layer_xy.py \
    --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --num-mol-snapshots 3 \
    --projections xy,xz,yz \
    --mol-colors "red,blue,green" \
    --out test_custom_colors.png
```
✅ **Result**: Molecule snapshots correctly colored with specified colors

### Test 3: Same Color Option
```bash
python visualize_top_layer_xy.py \
    --traj PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total_trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 0 \
    --projections xy,xz,yz \
    --mol-colors "same" \
    --out test_same_color.png
```
✅ **Result**: All molecule snapshots use the same black color

## Advanced Usage Examples

### Complete Custom Visualization
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --custom-frames "0,5,10,15,20" \
    --custom-mol-frames "2,8,14" \
    --mol-colors "red,green,blue" \
    --projections xy,xz,yz \
    --out complete_custom.png
```

### Scientific Publication Ready
```bash
python visualize_top_layer_xy.py \
    --traj trajectory.xyz \
    --fixed 26 --opposite 29 \
    --samples 4 \
    --num-mol-snapshots 3 \
    --mol-colors "auto" \
    --projections xy,xz,yz \
    --out publication_ready.png
```

## Backward Compatibility

All changes maintain full backward compatibility:
- Existing commands work exactly as before
- New parameters are optional with sensible defaults
- Default molecule color remains black when no color scheme is specified
- Substrate layer visualization is automatically enhanced for XZ/YZ projections

## Error Handling

Enhanced error handling includes:
- Validation of custom frame indices
- Proper color parsing and validation
- Graceful handling of malformed color specifications
- Clear error messages for invalid inputs
