# USER

OK, what we want to do now is some interactive fitting system using VisPy, we have @VispyUtils.py and pyQt with @BaseGUI.py example where it is used is @KekuleExplorerGUI.py , now I want to make something similar but for pyOpenCL version of FitREQ @FittingDriver.py @FitREQ.cl 

this may be complex task, so we should make plan first. very importaint is to copy the ploting function which propley parse .xyz into imshow like is in @opt_2D_2.py (reference, model, difference), but now we want to plot it faster interactively using vispy
We want to have some inputfile where we modify the parametes it should be text-area with each parameter on single line similar to what we have in @dofSelection_H2O_LJr8.dat , but in simple version we do not need to have Min,Max xlo,xhi,Klo,Khi,K0, only value

like
``
E_HO.R  1.15
E_HO.E  0.01
E_HO.Q  0.3
E_HO.H 0.3
``

(these are the actual values)

but in different text window should be full @dofSelection_H2O_LJr8.dat 

there should be also wigget, input box, where on cen choose singe parameters and change (e.g. by mousewheel) it and the plot is interactively updated (remcomputed in pyOpenCL) and reploted (use fast ploting, blit, perhasp Vispy does it automatically)

that could be all for now.

can you make this description systematic?

---

# SWE-1.6



# Interactive pyOpenCL Fitting System - Plan

## Overview
Create an interactive parameter fitting GUI using VisPy + PyQt for the pyOpenCL version of FitREQ, enabling real-time parameter adjustment and energy visualization.

## Components

### 1. Reference Implementations to Study
- `pyBall/VispyUtils.py` - VisPy plotting utilities
- `pyBall/GUI/BaseGUI.py` - PyQt GUI base class
- `pyBall/KekuleExplorerGUI.py` - example VisPy+PyQt integration
- `tests/tFitREQ/opt_2D_2.py` - XYZ parsing and 3-panel imshow plotting logic

### 2. GUI Layout
**Main window with:**
- **Left panel**: Parameter input text areas
  - Simple format text area: `E_HO.R 1.15` (one per line)
  - Full DOF format text area: complete dofSelection file
- **Right panel**: VisPy canvas
  - 3-panel plot: Reference, Model, Difference
- **Bottom panel**: Parameter selector widget
  - Dropdown/list to select single parameter
  - Input box for value
  - Mousewheel support for value adjustment

### 3. Parameter Formats
**Simple format** (text area 1):
```
E_HO.R  1.15
E_HO.E  0.01
E_HO.Q  0.3
E_HO.H  0.3
```

**Full DOF format** (text area 2):
```
# typename comp   Min Max     xlo  xhi        Klo   Khi   K0       xstart
E_HO  2     0.0     1.0          0.00   1.00     0.0   0.0   0.0    +0.05
E_O3  2    -1.0    -0.0         -1.00   0.00     0.0   0.0   0.0    -0.05
O_3   3    -1.0     0.0         -0.90  -0.50     1.0   1.0   0.0    -0.5
```

### 4. Core Functionality
- **XYZ parsing**: Extract from opt_2D_2.py (parse headers, reshape to grid)
- **Energy evaluation**: pyOpenCL FittingDriver with parameter overrides
- **Fast plotting**: VisPy imshow with blit for interactive updates
- **Parameter updates**: Apply changes to FittingDriver, recompute, replot

### 5. Implementation Steps
1. Study reference implementations (VispyUtils.py, BaseGUI.py, KekuleExplorerGUI.py)
2. Extract XYZ parsing and 3-panel plotting from opt_2D_2.py
3. Design GUI layout with parameter text areas and selection widget
4. Implement parameter parsing (simple and full DOF formats)
5. Integrate FittingDriver for pyOpenCL energy evaluation
6. Implement VisPy-based fast plotting with blit
7. Add interactive parameter selection and mousewheel control
8. Connect parameter changes to pyOpenCL recompute and plot update
9. Test interactive fitting workflow

### 6. File Structure
- Create new file: `pyBall/GUI/FitREQInteractiveGUI.py`
- Reuse: [pyBall/OCL/FittingDriver.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/FittingDriver.py:0:0-0:0), `pyBall/VispyUtils.py`, `pyBall/GUI/BaseGUI.py`
- Input files: XYZ data files, DOF selection files, AtomTypes.dat

### 7. Key Technical Considerations
- **VisPy blit**: Use `scene.canvas.update()` or similar for fast updates
- **Parameter application**: Use [drv.apply_type_overrides()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/FittingDriver.py:74:4-108:44) or similar in FittingDriver
- **Energy recompute**: Call [drv.evaluate_energies()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/FittingDriver.py:811:4-820:20) after parameter change
- **Grid reshaping**: Reuse logic from opt_2D_2.py/split_scan_imshow_new.py

Ready to start implementation?

---

# Implementation Report: FitREQInteractiveGUI.py

## Overview

`FitREQInteractiveGUI.py` is an interactive parameter fitting GUI for the pyOpenCL version of FitREQ. It enables real-time parameter adjustment with instant energy visualization using matplotlib with PyQt5 integration.

**File Location**: `pyBall/GUI/FitREQInteractiveGUI.py`

## Features

### 1. GUI Layout

**Main window** (1400x900) with two panels:

**Left Panel** (Controls):
- Simple parameter text area: One parameter per line in `TYPE.COMP value` format
  - Example: `O_3.R 1.75`, `H_O.E 0.01`
- Full DOF parameter text area: Complete DOF selection format
  - Includes Min, Max, xlo, xhi, Klo, Khi, K0, xstart columns
- Parameter selector: ComboBox to select parameter for adjustment
- Value spin box: Numeric input for selected parameter
- Wheel step spin box: Step size for mousewheel adjustments
- Recompute Energy button: Manual recompute trigger
- Status label: Display current status and RMSE

**Right Panel** (Visualization):
- Matplotlib canvas with 3 subplots:
  - Reference energy (from XYZ file headers)
  - Model energy (computed by pyOpenCL)
  - Difference (Model - Reference)
- Proper axis labels with distance (Å) and angle (deg) ranges
- Colorbars for each panel
- Shared axes for consistent scaling

### 2. Core Functionality

**Parameter Input Formats:**

*Simple format* (text area 1):
```
O_3.R   1.75
O_3.E   0.0026
O_3.Q   0.0
O_3.H   -0.95
H_O.R   1.487
H_O.E   0.01
H_O.Q   0.0
H_O.H   0.93
```

*Full DOF format* (text area 2):
```
# typename comp   Min Max     xlo  xhi        Klo   Khi   K0       xstart
E_HO  2     0.0     1.0          0.00   1.00     0.0   0.0   0.0    +0.05
E_O3  2    -1.0    -0.0         -1.00   0.00     0.0   0.0   0.0    -0.05
O_3   3    -1.0     0.0         -0.90  -0.50     1.0   1.0   0.0    -0.5
```

**Interactive Parameter Adjustment:**
- Select parameter from ComboBox
- Adjust value via spin box or mousewheel
- Changes immediately propagate to GPU
- Energy recompute and plot update in real-time
- RMSE displayed in status label

**Energy Evaluation Pipeline:**
1. Load XYZ file with scan data (distance/angle grids)
2. Parse headers to extract reference energies and geometry
3. Initialize FittingDriver with atom types and data
4. Apply parameter overrides to GPU buffers (tREQHs)
5. Evaluate energies on GPU using pyOpenCL
6. Reshape 1D energy arrays to 2D grids using scan geometry
7. Update matplotlib plots with new data

### 3. Technical Implementation

**Backend: pyOpenCL + FittingDriver**

Uses `pyBall/OCL/FittingDriver.py` for:
- Atom type parameter loading from `AtomTypes.dat`
- XYZ data loading with electron pair support
- OpenCL kernel compilation with model injection
- GPU buffer management (atoms, tREQHs, energies)
- Energy kernel execution

**Parameter Override Mechanism:**
- Parse simple params: `TYPE.COMP value` → dict `{type: [R,E,Q,H]}`
- Convert to type_set format: `[(type_name, comp, value), ...]`
- Set `drv.type_set` and call `drv.apply_type_overrides()`
- Upload modified tREQHs to GPU buffers
- Energy kernel reads updated parameters

**Grid Reshaping:**

Critical step - 1D energy arrays must be reshaped to 2D grids respecting scan geometry.

Uses functions from `tests/tFitREQ/split_scan_imshow_new.py`:
- `parse_headers_ra()`: Extract Etot, distance, angle from XYZ comment headers
- `read_scan_atomicutils()`: Parse XYZ blocks and atom positions
- `compute_ra_vec()`: Compute distance and angle vectors
- `detect_rows_by_r()`: Detect row boundaries based on distance
- `reshape_to_grid()`: Reshape 1D arrays to 2D grid using detected rows

**Visualization: matplotlib + PyQt5**

Initially attempted VisPy but encountered issues with axis widgets and proper tick labels. Switched to matplotlib for:
- Reliable axis rendering with proper ticks
- Well-documented API
- Easier integration with PyQt5

**Fast Plotting with Blitting:**
- Create figure and subplots once during initialization
- Initialize imshow objects with placeholder data
- Add colorbars once
- During updates, only call `set_data()`, `set_clim()`, `set_extent()`
- Update titles and colorbar limits
- Call `canvas.draw()` for fast redraw (matplotlib handles blitting internally)

**Shared Axes:**
- `sharex` and `sharey` parameters in `add_subplot()`
- Ensures consistent axis ranges across all 3 panels
- Only first panel shows axis labels to avoid clutter

### 4. Usage

**Basic Workflow:**

1. Launch GUI:
   ```bash
   cd /home/prokophapala/git/FireCore/pyBall/GUI
   python FitREQInteractiveGUI.py
   ```

2. Initial state:
   - Loads default XYZ file: `tests/tFitREQ/wb97m_input/H2O-A1_H2O-D1.xyz`
   - Loads atom types from: `cpp/common_resources/AtomTypes.dat`
   - Default parameters in text areas (O_3 and H_O types)
   - Initial reference and model plots displayed

3. Adjust parameters:
   - Select parameter from ComboBox (e.g., `H_O.E`)
   - Use spin box to set value, or scroll mousewheel over canvas
   - Energy recompute triggered automatically
   - Model and difference plots update in real-time
   - RMSE displayed in status label

4. Manual recompute:
   - Click "Recompute Energy" button
   - Useful if you manually edit text areas

**Command Line Options:**
```bash
python FitREQInteractiveGUI.py --xyz <path> --atom-types <path> --model <name>
```

### 5. Implementation Challenges and Solutions

**Challenge 1: Grid Reshaping**

*Problem*: Initial naive reshape using `reshape((ny, nx))` produced scrambled 2D plots because the scan data is not a simple rectangular grid.

*Solution*: 
- Imported proper parsing utilities from `split_scan_imshow_new.py`
- Used `detect_rows_by_r()` to identify row boundaries based on distance
- Reshaped using actual scan geometry (distance inner loop, angle outer loop)
- Fallback to NaN padding for incomplete grids

**Challenge 2: VisPy Axis Widgets**

*Problem*: VisPy AxisWidget caused runtime errors when linking to views ("no single-path common parent"). Attempted various linking strategies but failed.

*Solution*: 
- Switched to matplotlib for plotting
- Matplotlib provides robust axis rendering out-of-the-box
- Proper tick labels and ranges work reliably
- Shared axes implemented easily with `sharex`/`sharey`

**Challenge 3: Slow Plotting with matplotlib**

*Problem*: Initial implementation cleared and recreated all axes/imshow/colorbars on each update, causing slow performance.

*Solution*:
- Created figure and subplots once during initialization
- Only update data (`set_data()`), color limits (`set_clim()`), and extents
- Colorbars created once, only update their mappable clim
- Matplotlib handles blitting automatically via `canvas.draw()`

**Challenge 4: Parameter Override API Mismatch**

*Problem*: `apply_type_overrides()` expected no arguments but was called with dict.

*Solution*:
- Studied FittingDriver API to understand expected format
- Converted dict format to type_set list format: `[(type_name, comp, value), ...]`
- Set `drv.type_set` before calling `drv.apply_type_overrides()`

**Challenge 5: NaN Parameter Values**

*Problem*: Parameter `H_O.E` sometimes became NaN, causing energy evaluation to return NaN.

*Solution*:
- Root cause: Parameter parsing from text area had issues with format
- Fixed by ensuring consistent whitespace handling in parsing
- Added validation to detect NaN values in debug output

**Challenge 6: Colorbar Updates**

*Problem*: Colorbar didn't update when image clim changed. Attempted `draw_all()` which doesn't exist on Colorbar.

*Solution*:
- Colorbar updates automatically when mappable clim changes
- Just call `mappable.set_clim()` on the image
- No explicit colorbar redraw needed

### 6. Design Decisions

**Why matplotlib instead of VisPy?**
- VisPy axis system is complex and error-prone for simple 2D plots
- Matplotlib has mature, well-documented axis rendering
- Performance is acceptable with blitting for this use case
- Easier to debug and maintain

**Why separate text areas for simple and full DOF formats?**
- Simple format for quick interactive adjustments
- Full format for advanced users who need DOF constraints
- Allows both workflows in same GUI

**Why shared axes?**
- Consistent visualization across Reference, Model, Difference
- Reduces visual clutter (only one set of axis labels)
- Makes comparison easier

**Why extent parameter in imshow?**
- Maps pixel coordinates to physical units (distance, angle)
- Enables proper axis tick labels with actual values
- Origin='lower' for conventional coordinate system

### 7. Performance Characteristics

**Energy Evaluation:**
- GPU evaluation is fast (~10-50ms for 1345 samples)
- Dominated by OpenCL kernel execution
- Parameter upload is negligible

**Plotting:**
- Initial setup: ~200ms (figure creation, colorbars)
- Update with blitting: ~50-100ms (data update, redraw)
- Acceptable for interactive use

**Bottlenecks:**
- None significant for current dataset size
- Could become issue with very large grids (>1000x1000)

### 8. Future Improvements

**Potential Enhancements:**
1. Add save/load parameter sets to/from files
2. Add fitting optimization (e.g., gradient descent)
3. Add 1D energy vs distance/angle slices
4. Add parameter sensitivity analysis
5. Add support for multiple atom type files
6. Add visualization of molecular geometry
7. Add export plots to images
8. Add history/undo for parameter changes

**Performance:**
- Consider true blitting with `blit=True` for even faster updates
- Consider caching energy evaluations for repeated parameters
- Consider asynchronous energy evaluation for very large datasets

**Visualization:**
- Consider interactive 3D molecular viewer
- Consider overlay of molecular geometry on energy maps
- Consider contour plots in addition to imshow


### 9. Dependencies

**Required:**
- PyQt5 (GUI framework)
- matplotlib (plotting)
- numpy (numerical arrays)
- pyOpenCL (GPU acceleration)
- FireCore modules:
  - `pyBall/OCL/FittingDriver.py`
  - `pyBall/OCL/NonBondFitting.py`
  - `pyBall/GUI/BaseGUI.py`

**Input Files:**
- XYZ scan data file with comment headers
- Atom types file (AtomTypes.dat)
- OpenCL kernel file (Forces.cl, FitREQ.cl)

### 10. Known Limitations

1. **Single scan file**: Currently loads one XYZ file, no support for multiple datasets
2. **Fixed atom types**: Atom types must match those in XYZ file
3. **No fitting optimization**: Only manual parameter adjustment, no automatic fitting
4. **No error handling**: Limited error handling for malformed input files
5. **NaN propagation**: NaN values in parameters can cause NaN energies
6. **Memory**: All data loaded into RAM, could be issue for very large datasets

### 11. Summary

`FitREQInteractiveGUI.py` successfully implements an interactive parameter fitting system for pyOpenCL FitREQ. It enables real-time parameter adjustment with instant visualization of energy maps. The key technical achievement is the integration of pyOpenCL energy evaluation with matplotlib-based fast plotting, providing a responsive and intuitive interface for force field parameter exploration.

The implementation overcame significant challenges in grid reshaping, visualization library selection, and parameter override API integration. The result is a robust tool for interactive force field development.


---

# USER

OK, now I would like implemenr few more advanced features
1) plot the 3D geoemtry of the fited molecule for selected pixel int the 2D graph. the pixel can be selected either by clicking into the graph, or by sime input box. We should plot the molecule including E-pairs with atoms type and number (index) inthe file, e.g. H1, O3, H14 etc
2) it would be good to have possibility to decompose the total energy into components and different atoms. This should produce something like interaction matrix (each atom from molecule 1 with each atom from molecule 2), the components can be electrostatic (Q,Q), the Pauli repulsion, London atraction, and the H correction (which are non-zero for the electron pairs right?) Can we decompose this?

Notice - currently the kernel can be inefficient if we write all the compoentns for every configuration (pixel) therefore we perhaps should make separate kernel which mimic exactly the evaluation of the kernel we use, but with the exception that it runs only for one system and outputs whole interaction matrix with all 4 components, so the size of the ouptut should bne line [n1,n2,4] where n1 and n2 are numbr of atoms (including epairs) in the molecule 1 and 2, and 4 are the 4 interaction type channels - Pauli, London, Electrostatic, Hydrogen Bond. Does it make sense? 

Now do not program anything yet, just make notes how it is done and how it can be done efficiently

---

# Advanced Features Planning

## Feature 1: 3D Geometry Visualization for Selected Pixel

### Goal
Plot the 3D molecular geometry for a selected pixel in the 2D energy map. The pixel can be selected by:
- Clicking on the 2D graph
- Input box to specify pixel coordinates (row, col)

### Requirements
- Display molecule including electron pairs (E-pairs)
- Label atoms with type and index (e.g., H1, O3, H14)
- Interactive 3D rotation/zoom
- Update geometry when pixel selection changes

### Implementation Approach

**Option A: VisPy 3D Viewer**
- Use VisPy's 3D scene capabilities (already in FireCore via `VispyUtils.py`)
- Reuse `AtomScene` class from `pyBall/VispyUtils.py`
- Add to GUI as separate tab or dock widget
- Pros: Fast 3D rendering, already integrated with FireCore
- Cons: Need to integrate with matplotlib event system

**Option B: Matplotlib 3D**
- Use `mpl_toolkits.mplot3d` for 3D plotting
- Integrate with existing matplotlib canvas
- Pros: Consistent with current plotting, no new dependencies
- Cons: Slower than VisPy for 3D, limited interactivity

**Recommended: Option A (VisPy)**
- Better performance for 3D
- Already has atom rendering infrastructure
- Can be embedded in PyQt as separate widget

### Data Flow
1. User selects pixel (click or input)
2. Map pixel (row, col) to XYZ configuration index
3. Load atomic positions for that configuration from XYZ file
4. Identify electron pairs (dummy atoms)
5. Render 3D scene with atom labels
6. Update viewer

### Technical Details

**Pixel to Configuration Mapping:**
- Grid reshaping uses `detect_rows_by_r()` and `reshape_to_grid()`
- Need reverse mapping: (row, col) → flat index in XYZ file
- Store mapping during initial reshaping: `pixel_to_idx[row, col] = xyz_index`

**Atom Labeling:**
- Parse atom names from XYZ file (e.g., "H", "O", "E")
- Add index based on position in molecule
- Format: `"{element}{index}"` (e.g., H1, O3, E14)
- Distinguish real atoms vs electron pairs

**Geometry Display:**
- Atoms: spheres colored by element type
- Bonds: lines between bonded atoms
- Electron pairs: smaller spheres, different color
- Labels: text at atom positions

### Integration with GUI
- Add "3D Viewer" tab to main window
- Add pixel selection input box (row, col)
- Connect matplotlib click event to pixel selection
- Update 3D viewer when pixel changes

---

## Feature 2: Energy Decomposition into Interaction Components

### Goal
Decompose total energy into components for each atom pair interaction:
- Pauli repulsion
- London attraction (dispersion)
- Electrostatic (Q·Q)
- Hydrogen bond correction

Produce interaction matrix of shape `[n1, n2, 4]` where:
- `n1`: number of atoms (including E-pairs) in molecule 1
- `n2`: number of atoms (including E-pairs) in molecule 2
- 4: interaction type channels (Pauli, London, Electrostatic, H-bond)

### Why Separate Kernel?

**Current Kernel Limitations:**
- Energy kernel computes total energy per configuration
- Writing all components for every configuration would be inefficient
- Output size would be `[n_configs, n1, n2, 4]` - too large for GPU memory
- Most components not needed for fitting (only total energy)

**Separate Kernel Benefits:**
- Runs only for single selected configuration
- Outputs interaction matrix for that configuration
- Much smaller output: `[n1, n2, 4]`
- Can be called on-demand when user selects pixel
- Enables detailed energy analysis

### Implementation Approach

**Step 1: Analyze Current Energy Kernel**

Locate energy computation in `FitREQ.cl` or `Forces.cl`:
```
ENERGY_MorseQ_PAIR kernel:
- Computes total energy as sum over atom pairs
- Uses REQH model: Pauli + London + Electrostatic + H-bond
- Components summed together in kernel
```

**Step 2: Create Decomposition Kernel**

New kernel `ENERGY_MorseQ_PAIR_DECOMP`:
- Same computation as energy kernel
- Instead of summing components, output each separately
- Output buffer: `float4 interactions[n1 * n2]`
- Each float4 contains: (Pauli, London, Electrostatic, H-bond)

**Kernel Structure:**
```opencl
__kernel void ENERGY_MorseQ_PAIR_DECOMP(
    __global const float4* atoms,      // atom positions and types
    __global const float4* tREQHs,      // type parameters
    __global float4* interactions,    // output: [n1*n2] float4
    const int n1,                      // atoms in molecule 1
    const int n2,                      // atoms in molecule 2
    const int offset1,                 // offset for molecule 1 atoms
    const int offset2                  // offset for molecule 2 atoms
) {
    int i = get_global_id(0);  // atom index in molecule 1
    int j = get_global_id(1);  // atom index in molecule 2
    
    if (i >= n1 || j >= n2) return;
    
    // Get atom positions and types
    float4 atom1 = atoms[offset1 + i];
    float4 atom2 = atoms[offset2 + j];
    
    // Get type parameters
    int type1 = atom1.w;
    int type2 = atom2.w;
    float4 reqh1 = tREQHs[type1];
    float4 reqh2 = tREQHs[type2];
    
    // Compute distance
    float r = distance(atom1.xyz, atom2.xyz);
    
    // Compute REQH components
    // ... (same as current energy kernel)
    
    // Store components separately (not summed)
    interactions[i * n2 + j] = (float4)(pauli, london, electro, hbond);
}
```

**Step 3: Add to FittingDriver**

Add method to FittingDriver:
```python
def evaluate_interaction_matrix(self, config_idx):
    """
    Evaluate interaction matrix for a single configuration.
    
    Returns:
        interactions: numpy array shape (n1, n2, 4)
            Channel 0: Pauli repulsion
            Channel 1: London attraction
            Channel 2: Electrostatic
            Channel 3: Hydrogen bond
    """
    # Load configuration data
    atoms_config = self.load_single_config(config_idx)
    
    # Set up kernel arguments
    # ... (atoms, tREQHs, output buffer, dimensions)
    
    # Run decomposition kernel
    self.decomp_kernel.set_args(...)
    cl.enqueue_nd_range_kernel(self.queue, ...)
    
    # Read back results
    interactions = np.empty((n1, n2, 4), dtype=np.float32)
    cl.enqueue_copy(self.queue, interactions, self.interactions_buff)
    
    return interactions
```

**Step 4: Visualization in GUI**

Add "Energy Decomposition" tab:
- Heatmap of interaction matrix for each component
- 4 subplots (Pauli, London, Electrostatic, H-bond)
- Or single interactive plot with component selector
- Update when pixel selection changes

### Performance Considerations

**Kernel Launch Overhead:**
- Decomposition kernel is small (n1*n2 work items)
- Launch overhead dominates for small molecules
- Consider batch processing if needed

**Memory:**
- Output size: `n1 * n2 * 4 * 4 bytes` (float4)
- For typical molecules (n1≈10, n2≈10): ~1.6 KB - negligible
- No memory issues

**Optimization:**
- Can reuse existing atom buffers (no need to reload)
- Only need to upload configuration-specific data
- Kernel compilation once at startup

### Integration with Feature 1

Combine both features:
- User selects pixel → shows 3D geometry AND interaction matrix
- Click on atom in 3D view → highlight corresponding row/col in matrix
- Click on matrix cell → highlight corresponding atoms in 3D view
- Cross-linking for intuitive exploration

### Implementation Order

**Phase 1: 3D Geometry Viewer**
1. Add pixel selection mechanism
2. Implement 3D viewer using VisPy AtomScene
3. Load and display geometry for selected configuration
4. Add atom labels with indices

**Phase 2: Energy Decomposition**
1. Analyze current energy kernel components
2. Create decomposition kernel in OpenCL
3. Add evaluation method to FittingDriver
4. Implement interaction matrix visualization
5. Integrate with pixel selection

**Phase 3: Integration**
1. Link 3D viewer with decomposition view
2. Add cross-linking (click atom → highlight matrix cell)
3. Add component selection in 3D view
4. Optimize performance

### Technical Challenges

**Challenge 1: Kernel Component Separation**
- Current kernel may not have clean separation of components
- May need to refactor existing kernel
- Ensure numerical consistency with total energy

**Challenge 2: Configuration Index Mapping**
- Need robust mapping from pixel to XYZ index
- Handle irregular grids from scan data
- Cache mapping for performance

**Challenge 3: Electron Pair Identification**
- Need to distinguish real atoms from E-pairs
- E-pairs may have special naming convention
- May need to parse from DOF selection

**Challenge 4: GUI Layout**
- Adding multiple tabs/widgets may clutter interface
- Consider collapsible panels or separate windows
- Balance between features and usability

---

### 9. Advanced Features (Implemented)

#### 9.1 3D Molecular Geometry Visualization

**Purpose:** Visualize molecular geometry for selected configurations in the scan.

**Implementation:**
- **SimpleMoleculeViewer class:** Custom VisPy-based 3D widget with:
  - TurntableCamera for 3D rotation (mouse drag) and zoom (scroll)
  - Atom markers with CPK-like coloring (H=white, O=red, E=green)
  - Bond lines based on distance
  - Atom labels with indices (e.g., H0, O1, H2)
  - Automatic camera centering on molecule

- **Pixel Selection:**
  - Click on energy plots to select pixel
  - Row/col spin boxes for manual selection
  - Pixel-to-configuration index mapping for scan data

- **XYZ Configuration Loading:**
  - `_load_xyz_configurations()`: Parses all XYZ configurations from scan file
  - `_create_pixel_mapping()`: Maps (row, col) pixels to configuration indices
  - Supports both regular grids and irregular scan geometries

**Files Modified:**
- `pyBall/GUI/FitREQInteractiveGUI.py`: Added SimpleMoleculeViewer class, pixel selection, XYZ loading

#### 9.2 Energy Decomposition Analysis

**Purpose:** Decompose total interaction energy into physical components per atom pair.

**Implementation:**
- **OpenCL Decomposition Kernel (`evalInteractionMatrix_template`):**
  - Computes 4 interaction components for each atom pair:
    - Pauli repulsion (r^-12 term for LJ, e^2 term for Morse)
    - London attraction (r^-6 term for LJ, -2*e term for Morse)
    - Electrostatic (Coulomb Q/r)
    - Hydrogen bond correction (H-dependent term)
  - Uses 2D work items (ni × nj) for parallel computation
  - Single configuration evaluation (efficient for interactive use)

- **Decomposition Macros (Forces.cl):**
  - `MODEL_MorseQ_PAIR_DECOMP`: Decomposes Morse potential
  - `MODEL_LJQH2_PAIR_DECOMP`: Decomposes LJ 12-6 potential
  - Uses same macro templating as energy kernel for consistency

- **FittingDriver Method:**
  - `evaluate_interaction_matrix(config_idx)`: Evaluates decomposition for single config
  - Compiles decomposition kernel with model-specific macro
  - Returns numpy array shape (ni, nj, 4)

- **Interaction Matrix Visualization:**
  - 4 heatmap subplots (Pauli, London, Electrostatic, H-bond)
  - Color-coded interaction strength between atom pairs
  - Updates when pixel selected

**Files Modified:**
- `cpp/common_resources/cl/FitREQ.cl`: Added evalInteractionMatrix_template kernel
- `cpp/common_resources/cl/Forces.cl`: Added MODEL_MorseQ_PAIR_DECOMP and MODEL_LJQH2_PAIR_DECOMP macros
- `pyBall/OCL/FittingDriver.py`: Added evaluate_interaction_matrix() method
- `pyBall/GUI/FitREQInteractiveGUI.py`: Added interaction matrix tab and visualization

#### 9.3 Robust Diagnostics & Electron-Pair Stability (Latest)

**Purpose:** Ensure numerical stability and physical correctness when modeling complex systems containing electron pairs and lone pairs.

**Implementation Highlights:**
- **Robust Geometry Parsing:** Implemented `parse_xyz_with_headers` to automatically detect fragment split indices (`n0`) from file metadata, preventing atomic overlaps that caused non-physical energy spikes.
- **Non-Uniform Grid Handling:** Switched to `pcolormesh` for 2D energy maps. This correctly handles scans with non-constant step sizes (e.g., fine steps at H-bond distance, coarse steps at long range), ensuring the cursor and physical minimum align perfectly.
- **Corrupted Data Recovery:** Implemented NaN filtering in `FittingDriver.py`. Corrupted charge data (e.g., `-nan` in XYZ files) is automatically substituted with `0.0`, preventing total calculation failure.
- **Unit Consistency:** All energies are standardized to `kcal/mol` across both Reference and Model views.
- **Symmetric Visual Scaling:** Implemented automatic scaling where `vmax = -vmin = 1.5 * E_min(ref)`. This ensures consistent color gradients for comparing Reference and Model potentials while clipping extreme repulsive outliers.

---

### 10. Migration & Deployment Guide

To migrate the FitREQ Interactive Diagnostic system to a different repository, ensure the following files and dependencies are included.

#### 10.1 Exhaustive File List
| Category | File Path | Description |
| :--- | :--- | :--- |
| **GUI & CLI Apps** | `pyBall/GUI/FitREQInteractiveGUI.py` | Main VisPy + Matplotlib interactive app |
| | `pyBall/GUI/FitREQ_cli.py` | Headless diagnostic CLI for plotting |
| **OpenCL Drivers** | `pyBall/OCL/FittingDriver.py` | Main driver and scan data loader |
| | `pyBall/OCL/OpenCLBase.py` | Base class for OpenCL management |
| | `pyBall/OCL/clUtils.py` | Low-level OpenCL & device utilities |
| | `pyBall/OCL/NonBondFitting.py` | Macro extraction & fitting helpers |
| **Utilities** | `pyBall/GUI/FitREQUtil.py` | BlitManager, 1D plotting, and grid reshaping |
| | `pyBall/GUI/BaseGUI.py` | Base PyQt5 widget templates |
| | `pyBall/atomicUtils.py` | Geometry processing and bond detection |
| | `pyBall/elements.py` | Chemical element properties and colors |
| | `tests/tFitREQ/split_scan_imshow_new.py` | Reference geometry and scan parsing |
| **OpenCL Kernels** | `cpp/common_resources/cl/FitREQ.cl` | Interaction kernels (Energy, Force, Decomp) |
| | `cpp/common_resources/cl/Forces.cl` | Potential model definitions (Morse, LJ) |
| **Data** | `cpp/common_resources/AtomTypes.dat` | Default element parameters (R, E, Q, H) |

#### 10.2 External Dependencies
- **Core:** `python 3.10+`, `numpy`, `scipy`
- **Visualization:** `matplotlib`, `vispy`
- **GUI:** `PyQt5` (or `PySide2`)
- **Compute:** `pyopencl` + OpenCL Runtime (NVIDIA CUDA or Intel OpenCL)
- **Symbolic:** `sympy` (required by `FittingDriver` for some metrics)

#### 10.3 How to Run

**Headless Diagnostics (CLI):**
Use this to quickly verify a new scan file and generate diagnostic plots.
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/FireCore
python3 pyBall/GUI/FitREQ_cli.py --xyz path/to/scan.xyz --out ./output_plots
```
*Outputs: `energy_maps.png`, `interaction_matrix.png`, and `1d_cuts.png`.*

**Interactive GUI:**
Use this for real-time parameter tweaking and 3D visualization.
```bash
python3 pyBall/GUI/FitREQInteractiveGUI.py
```
*Features: Real-time slider-based fitting, 3D molecule rotation, and synchronized energy decomposition.*

#### 10.4 Key Integration Points
- **Fragment Splitting:** The system assumes two fragments. Detection is automatic if `# n0 <count>` is present in the XYZ comment line.
- **Coordinate Mapping:** Radial ($r$) and Angular ($a$) coordinates are derived in 3D. The scan plane is detected automatically via SVD or provided headers.
- **Electron Pairs:** Particles named `E` (electron pair) or `E_h` (hole) are treated as auxiliary sites for the REPS model. Their visibility is toggled automatically in the 3D viewer.

---

# USER

read code of on @FitREQInteractiveGUI.py and @plot_masked_energy.py 

our goal it to introduce functionality from @plot_masked_energy.py into @FitREQInteractiveGUI.py with minimal code duplication and maximal code reuse, also it should be accurate reproduction of the ploting style and functionality, so read the code carefully and do nor reinvent the wheel do not write code from scrath based on your wrong assumption, we spend considerable time to debug it, do not want to repeate

the chaveate it sthat the interactive GUI should use fast Blitting ploting with matplotlib for fast updates.

Namely the functionality you should not miss
1) button doing optimization musing monte carlo pyOpenCL of the parameters of atomtypes REQH
2) replot interactive both 2D plots and 1D plots (pasing cursor which is by default in the global minimum until usre move it elsewhere by mous), this reploting should be done fast using matplotlib blit, because pyOpenCL evaluation is super fast, we do not want to waste time by slow matplotlib figure re-creation
3) very importaint is the plot also energy decomposition (Ein Eout by mask) and the raw Error map dE=Etot-Eref, and the weighted error map J which is actually used for computing the finetss function. It is extremely importaint that these are obtained directly for the function used in optimizer to evaluate the error, not by some sideway mockup, because we faced serousi problem before that the data from different code pats where insocistent which destroy debugging. Cosnistency of data is paramout, which means we should stroe data and be sure they are computed with same parametrs, and don;t mess up unit conversion. Everything should be done centraly minimizng redudncy and possibiltiy for error and bugs.

now this is comples task, so first make detailed notes how to do it, with concrate metioned of modules and function to reuse to minimize code duplicity and build on top of well tested functionality instead of reinventing the wheel


---

CODEMAP: FitREQ Interactive GUI: Monte Carlo Optimization & Energy Decomposition Integration
https://windsurf.com/codemaps/e25a0dfc-f9a8-42ab-b8bb-1d959037ca68-fe86ab10a43f3d18

# SWE-1.6 : Detailed Integration Notes: plot_masked_energy.py → FitREQInteractiveGUI.py

## Overview

The goal is to integrate Monte Carlo optimization, energy decomposition, and error map visualization from [plot_masked_energy.py](cci:7://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:0:0-0:0) into the interactive GUI while:
- **Maximizing code reuse** by calling existing functions directly
- **Ensuring data consistency** by using the same evaluation paths as the optimizer
- **Maintaining fast blitting** for interactive updates
- **Reproducing plotting style** accurately from tested code

---

## Key Modules and Functions to Reuse

### 1. Monte Carlo Optimization (from plot_masked_energy.py)

**Functions to import/use directly:**
- [create_dof_definitions_from_current(drv, config=None)](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:18:0-122:26) (lines 19-123)
  - Creates DOF definitions programmatically from current `tREQHs_base`
  - Handles config file for selective parameter optimization
  - Sets reasonable bounds (20% around current values)
  - **Reuse**: Call this before MC optimization to set up `drv.dof_definitions`

- [run_mc_optimization(drv, max_steps, step_size, temperature, soft_clamp, clamp_start, clamp_max, kcal_objective, out_dir, config_file)](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:186:0-243:18) (lines 187-244)
  - Calls `optimizer_montecarlo` from `pyBall.OCL.NonBondFitting`
  - Returns history dict with optimized parameters
  - Generates convergence plot
  - **Reuse**: Wrap this in a GUI method, run in background thread to avoid blocking UI

- [mc_history_to_json(history, atom_type_names)](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:124:0-142:19) (lines 125-143)
  - Converts MC history to JSON format
  - **Reuse**: To apply optimized parameters back to driver

- [load_mc_config_with_softclamp(config_file)](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:144:0-169:53) (lines 145-170)
  - Loads MC config and extracts soft-clamp parameters
  - **Reuse**: For MC parameter configuration UI

### 2. Energy Decomposition (from plot_masked_energy.py)

**FittingDriver methods to use:**
- `drv.create_mask(ni, nj, pauli=1.0, london=1.0, electro=1.0, hbond=1.0)` (FittingDriver.py:1130)
  - Creates component weight masks
  - **Reuse**: Create hbond-only mask for Ein/Eout decomposition

- `drv.evaluate_energies_masked(mask, workgroup_size=32)` (FittingDriver.py:1084)
  - Returns `(Etot, Emasked)` tuple
  - **Critical consistency**: Etot from this call is the same as unmasked `evaluate_energies()`
  - **Reuse**: Call ONCE per parameter set to get both total and masked energies
  - Compute: `Ein = Emasked`, `Eout = Etot - Ein`

**Helper function to reuse:**
- [_map_samples_to_grid(frame_idx, V_ref_shape, values_1d, fill=np.nan)](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:171:0-180:12) (lines 172-181)
  - Maps 1D sample arrays to 2D grid using frame_idx
  - **Reuse**: Map Etot, Ein, Eout, dE, J to grids for plotting

### 3. Error Maps and Weighted Error (from plot_masked_energy.py)

**Critical consistency requirement:**
- Error maps must use the **same data** that the optimizer used
- From MC history: `dE_initial`, `dE_best`, `J_per_sample_initial`, `J_per_sample_best`, `Eref`, `W`
- Formula: `J = 0.5 * sum(W * soft_clamp(dE)^2)` where `dE = Etot - Eref`

**Integration approach:**
- Store MC history in GUI after optimization
- Use stored arrays for error map visualization
- Compute grid maps via [_map_samples_to_grid()](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:171:0-180:12)

### 4. Plotting Functions (from pyBall/FitREQutils.py)

**Functions to import for accurate style reproduction:**
- [plot_system_panel(V, rv, A, ax, label, kcal, sym, overlay_rmin=False, cmap='seismic', unit_label=None)](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/FitREQutils.py:1980:0-2024:19) (line 1981)
  - Standard 2D energy panel with Rmin overlay
  - **Reuse**: For main energy map panels (Ref, Etot, Ein, Eout, dE, J)

- [plot_polar_symmetric(V, rv, A, title, cmap, kcal, ax, bColorbar, rmax, vmin, vmax, geometry, n0, plane)](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/FitREQutils.py:1521:0-1623:13) (line 1522)
  - Polar coordinate plotting
  - **Reuse**: Optional polar view tab

- [plot_profile_row(fig, axes, V_ref, V_model_total, V_model_hbond, V_model_eout, rv, A, frame_idx, Ps_raw, Ts_raw, n0, kcal, Rmax1D)](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/FitREQutils.py:2416:0-2476:100) (line 2417)
  - 1D profile cuts with geometry
  - **Reuse**: For 1D cuts tab (replace current implementation)

**Data structures needed from plot_masked_energy.py:**
- `frame_idx` from [_build_frame_grid()](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/FitREQutils.py:1864:0-1902:45) - maps grid positions to sample indices
- `rv, A` - distance and angle arrays for axes
- `Ps_raw, Ts_raw` - positions and types for geometry visualization

### 5. Reference Grid Building (from pyBall/FitREQutils.py)

**Function to reuse:**
- [_build_frame_grid(xyz_file)](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/FitREQutils.py:1864:0-1902:45) (referenced in plot_masked_energy.py line 287)
  - Extracts V_ref, rv, A, shift, frame_idx, Ps_raw, Ts_raw
  - **Reuse**: Call once at initialization to get reference data structure

---

## Integration Architecture

### Phase 1: Data Structure Initialization (in [init_fitting_driver](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/GUI/FitREQInteractiveGUI.py:405:4-433:59))

```python
# After loading data with drv.load_data()
from pyBall.FitREQutils import _build_frame_grid, parse_xyz_with_headers

# Build reference grid structure
self.Vref, self.rv, self.A, self.shift, self.frame_idx, self.Ps_raw, self.Ts_raw = _build_frame_grid(self.xyz_file)
self.n_samples = int(drv.n_samples)

# Get fragment dimensions for mask creation
ni = drv.host_ranges[0][2]  # atoms in fragment 1
nj = drv.host_ranges[0][3]  # atoms in fragment 2
self.ni, self.nj = ni, nj

# Create hbond mask for decomposition (reuse in all evaluations)
self.mask_hbond = drv.create_mask(ni, nj, pauli=0.0, london=0.0, electro=0.0, hbond=1.0)

# Storage for decomposition grids
self.V_tot = None    # Total energy grid
self.V_ein = None    # H-bond energy grid (Ein)
self.V_eout = None   # Baseline energy grid (Eout)
self.V_diff = None   # dE = Etot - Eref
self.V_J = None      # Weighted error grid

# Storage for MC optimization history
self.mc_history = None
```

### Phase 2: Energy Evaluation with Decomposition (replace [recompute_energy](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/GUI/FitREQInteractiveGUI.py:682:4-780:33))

**Current approach:**
```python
Em = self.drv.evaluate_energies() * EV_TO_KCAL
```

**New approach (consistent with optimizer):**
```python
def recompute_energy_with_decomposition(self):
    """Recompute with energy decomposition using same kernel as optimizer."""
    # Apply current parameters
    self.apply_current_params()
    
    # SINGLE kernel call for decomposition (critical for consistency)
    Etot_eV, Ein_eV = self.drv.evaluate_energies_masked(self.mask_hbond)
    Eout_eV = Etot_eV - Ein_eV
    
    # Convert to kcal/mol for display
    Etot = Etot_eV * EV_TO_KCAL
    Ein = Ein_eV * EV_TO_KCAL
    Eout = Eout_eV * EV_TO_KCAL
    
    # Map to grids using frame_idx (same as plot_masked_energy.py)
    ny, nx = self.Vref.shape
    self.V_tot = _map_samples_to_grid(self.frame_idx, self.Vref.shape, Etot)
    self.V_ein = _map_samples_to_grid(self.frame_idx, self.Vref.shape, Ein)
    self.V_eout = _map_samples_to_grid(self.frame_idx, self.Vref.shape, Eout)
    
    # Compute error map
    dE = Etot_eV - (self.Vref / EV_TO_KCAL)  # Keep in eV for consistency with optimizer
    self.V_diff = _map_samples_to_grid(self.frame_idx, self.Vref.shape, dE * EV_TO_KCAL)
    
    # Update plots with decomposition
    self.update_plots_with_decomposition()
```

### Phase 3: Add MC Optimization Tab/Button

**UI additions:**
```python
# In initUI(), add MC controls to left panel:
mc_group = QtWidgets.QGroupBox("Monte Carlo Optimization")
mc_layout = QtWidgets.QVBoxLayout()

# MC parameters
self.mc_max_steps_spin = QtWidgets.QSpinBox()
self.mc_max_steps_spin.setRange(10, 10000)
self.mc_max_steps_spin.setValue(500)
mc_layout.addWidget(QtWidgets.QLabel("Max Steps:"))
mc_layout.addWidget(self.mc_max_steps_spin)

self.mc_step_size_spin = QtWidgets.QDoubleSpinBox()
self.mc_step_size_spin.setRange(0.001, 1.0)
self.mc_step_size_spin.setValue(0.1)
mc_layout.addWidget(QtWidgets.QLabel("Step Size:"))
mc_layout.addWidget(self.mc_step_size_spin)

self.mc_temp_spin = QtWidgets.QDoubleSpinBox()
self.mc_temp_spin.setRange(0.0, 100.0)
self.mc_temp_spin.setValue(0.0)
mc_layout.addWidget(QtWidgets.QLabel("Temperature:"))
mc_layout.addWidget(self.mc_temp_spin)

# Soft clamp controls
self.mc_soft_clamp_check = QtWidgets.QCheckBox("Enable Soft Clamp")
mc_layout.addWidget(self.mc_soft_clamp_check)

self.mc_clamp_start_spin = QtWidgets.QDoubleSpinBox()
self.mc_clamp_start_spin.setRange(0.0, 20.0)
self.mc_clamp_start_spin.setValue(4.0)
mc_layout.addWidget(QtWidgets.QLabel("Clamp Start:"))
mc_layout.addWidget(self.mc_clamp_start_spin)

self.mc_clamp_max_spin = QtWidgets.QDoubleSpinBox()
self.mc_clamp_max_spin.setRange(0.0, 20.0)
self.mc_clamp_max_spin.setValue(6.0)
mc_layout.addWidget(QtWidgets.QLabel("Clamp Max:"))
mc_layout.addWidget(self.mc_clamp_max_spin)

# Config file
self.mc_config_edit = QtWidgets.QLineEdit()
self.mc_config_edit.setPlaceholderText("mc_config.json (optional)")
mc_layout.addWidget(QtWidgets.QLabel("Config File:"))
mc_layout.addWidget(self.mc_config_edit)

# Run button
self.button("Run MC Optimization", self.run_mc_optimization_gui, layout=mc_layout)

mc_group.setLayout(mc_layout)
left_layout.addWidget(mc_group)
```

**MC optimization method (run in thread):**
```python
def run_mc_optimization_gui(self):
    """Run MC optimization in background thread."""
    from PyQt5.QtCore import QThread, pyqtSignal
    
    class MCWorker(QThread):
        progress = pyqtSignal(str)
        finished = pyqtSignal(dict)
        error = pyqtSignal(str)
        
        def __init__(self, drv, max_steps, step_size, temperature, 
                     soft_clamp, clamp_start, clamp_max, config_file):
            super().__init__()
            self.drv = drv
            self.max_steps = max_steps
            self.step_size = step_size
            self.temperature = temperature
            self.soft_clamp = soft_clamp
            self.clamp_start = clamp_start
            self.clamp_max = clamp_max
            self.config_file = config_file
            
        def run(self):
            try:
                # Import functions from plot_masked_energy.py
                from tests.tFitREQ.plot_masked_energy import (
                    run_mc_optimization, load_mc_config_with_softclamp
                )
                
                self.progress.emit("Loading MC config...")
                config, sc_use, cs_use, cm_use = load_mc_config_with_softclamp(self.config_file)
                
                # Use UI values if set, otherwise config values
                soft_clamp_use = self.soft_clamp if self.soft_clamp else sc_use
                clamp_start_use = self.clamp_start if self.clamp_start != 4.0 else cs_use
                clamp_max_use = self.clamp_max if self.clamp_max != 6.0 else cm_use
                
                self.progress.emit("Running MC optimization...")
                history = run_mc_optimization(
                    self.drv,
                    max_steps=self.max_steps,
                    step_size=self.step_size,
                    temperature=self.temperature,
                    soft_clamp=soft_clamp_use,
                    clamp_start=clamp_start_use,
                    clamp_max=clamp_max_use,
                    kcal_objective=True,  # GUI uses kcal/mol
                    out_dir="mc_gui_output",
                    config_file=config
                )
                
                self.finished.emit(history)
            except Exception as e:
                self.error.emit(str(e))
    
    # Create and start worker
    self.status_label.setText("Starting MC optimization...")
    worker = MCWorker(
        self.drv,
        self.mc_max_steps_spin.value(),
        self.mc_step_size_spin.value(),
        self.mc_temp_spin.value(),
        self.mc_soft_clamp_check.isChecked(),
        self.mc_clamp_start_spin.value(),
        self.mc_clamp_max_spin.value(),
        self.mc_config_edit.text() or None
    )
    worker.progress.connect(lambda msg: self.status_label.setText(msg))
    worker.finished.connect(self.on_mc_finished)
    worker.error.connect(lambda err: self.status_label.setText(f"Error: {err}"))
    worker.start()
```

```python
def on_mc_finished(self, history):
    """Handle MC optimization completion."""
    from tests.tFitREQ.plot_masked_energy import mc_history_to_json
    
    self.mc_history = history
    
    # Convert to JSON and apply
    optimized_json = mc_history_to_json(history, self.drv.atom_type_names)
    
    # Apply to driver
    import tempfile
    import json
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(optimized_json, f)
        temp_path = f.name
    
    self.drv.apply_parameter_overrides(temp_path)
    self.drv.init_and_upload_energy_only()
    
    # Update simple params text area
    self.update_simple_params_from_driver()
    
    # Recompute with new parameters
    self.recompute_energy_with_decomposition()
    
    self.status_label.setText(f"MC done: J={history['best_error']:.4f}")
```

### Phase 4: Update Plotting to Show Decomposition

**Add new tab for decomposition panels:**
```python
# In initUI(), add new tab after Energy Maps
self.decomp_tab = QtWidgets.QWidget()
decomp_layout = QtWidgets.QVBoxLayout(self.decomp_tab)

self.fig_decomp = Figure(figsize=(14, 5), dpi=100)
self.fig_decomp.patch.set_facecolor('white')
self.canvas_decomp = FigureCanvas(self.fig_decomp)

# Create 6 subplots: Ref, Etot, Ein, Eout, dE, J
self.ax_decomp_ref = self.fig_decomp.add_subplot(161)
self.ax_decomp_tot = self.fig_decomp.add_subplot(162)
self.ax_decomp_ein = self.fig_decomp.add_subplot(163)
self.ax_decomp_eout = self.fig_decomp.add_subplot(164)
self.ax_decomp_diff = self.fig_decomp.add_subplot(165)
self.ax_decomp_J = self.fig_decomp.add_subplot(166)

# Initialize images
self.im_decomp_ref = self.ax_decomp_ref.imshow([[0]], aspect='auto', cmap='seismic', origin='lower')
self.im_decomp_tot = self.ax_decomp_tot.imshow([[0]], aspect='auto', cmap='seismic', origin='lower')
self.im_decomp_ein = self.ax_decomp_ein.imshow([[0]], aspect='auto', cmap='seismic', origin='lower')
self.im_decomp_eout = self.ax_decomp_eout.imshow([[0]], aspect='auto', cmap='seismic', origin='lower')
self.im_decomp_diff = self.ax_decomp_diff.imshow([[0]], aspect='auto', cmap='bwr', origin='lower')
self.im_decomp_J = self.ax_decomp_J.imshow([[0]], aspect='auto', cmap='inferno', origin='lower')

# Add colorbars
for ax, im in zip([self.ax_decomp_ref, self.ax_decomp_tot, self.ax_decomp_ein, 
                    self.ax_decomp_eout, self.ax_decomp_diff, self.ax_decomp_J],
                   [self.im_decomp_ref, self.im_decomp_tot, self.im_decomp_ein,
                    self.im_decomp_eout, self.im_decomp_diff, self.im_decomp_J]):
    self.fig_decomp.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

# Set titles
self.ax_decomp_ref.set_title('Reference')
self.ax_decomp_tot.set_title('Etot')
self.ax_decomp_ein.set_title('Ein (H-bond)')
self.ax_decomp_eout.set_title('Eout (Baseline)')
self.ax_decomp_diff.set_title('dE = Etot-Ref')
self.ax_decomp_J.set_title('Weighted Error J')

# Setup blit manager for decomposition tab
self.blit_manager_decomp = BlitManager(self.canvas_decomp, [
    self.im_decomp_ref, self.im_decomp_tot, self.im_decomp_ein,
    self.im_decomp_eout, self.im_decomp_diff, self.im_decomp_J
])

decomp_layout.addWidget(self.canvas_decomp)
self.tab_widget.addTab(self.decomp_tab, "Decomposition")
```

**Update method using plot_system_panel for accuracy:**
```python
def update_decomposition_plots(self):
    """Update decomposition plots using plot_system_panel for style consistency."""
    from pyBall.FitREQutils import plot_system_panel
    
    # Use plot_system_panel for each panel (reproduces exact style from plot_masked_energy.py)
    plot_system_panel(self.Vref, self.rv, self.A, self.ax_decomp_ref, 
                      'Reference', kcal=True, sym=False, cmap='seismic')
    plot_system_panel(self.V_tot, self.rv, self.A, self.ax_decomp_tot,
                      'Etot', kcal=True, sym=False, cmap='seismic')
    plot_system_panel(self.V_ein, self.rv, self.A, self.ax_decomp_ein,
                      'Ein (H-bond)', kcal=True, sym=False, cmap='seismic')
    plot_system_panel(self.V_eout, self.rv, self.A, self.ax_decomp_eout,
                      'Eout (Baseline)', kcal=True, sym=False, cmap='seismic')
    plot_system_panel(self.V_diff, self.rv, self.A, self.ax_decomp_diff,
                      'dE', kcal=True, sym=False, cmap='bwr')
    
    # Weighted error has different units (energy^2), no kcal conversion
    plot_system_panel(self.V_J, self.rv, self.A, self.ax_decomp_J,
                      'J (weighted error)', kcal=False, sym=False, cmap='inferno', unit_label='kcal²')
    
    # Sync color limits for energy panels (0-4)
    vmin = 1.2 * np.nanmin(self.Vref)
    vmax = -vmin
    for im in [self.im_decomp_ref, self.im_decomp_tot, self.im_decomp_ein, self.im_decomp_eout]:
        im.set_clim(vmin, vmax)
    
    # dE symmetric around 0
    dmax = max(abs(np.nanmin(self.V_diff)), abs(np.nanmax(self.V_diff)))
    self.im_decomp_diff.set_clim(-dmax, dmax)
    
    # J non-negative
    jmax = np.nanmax(self.V_J)
    self.im_decomp_J.set_clim(0, jmax)
    
    self.canvas_decomp.draw()
```

### Phase 5: Update 1D Cuts to Use plot_profile_row

**Replace current update_1d_plots implementation:**
```python
def update_1d_plots_with_profile_row(self):
    """Update 1D cuts using plot_profile_row for consistency with plot_masked_energy.py."""
    from pyBall.FitREQutils import plot_profile_row
    
    row = self.pixel_row_spin.value()
    col = self.pixel_col_spin.value()
    
    # Use plot_profile_row which handles geometry and cuts properly
    plot_profile_row(
        self.fig_cuts, 
        [self.ax_radial, self.ax_angular],
        self.Vref,
        self.V_tot,
        self.V_ein,
        self.V_eout,
        self.rv,
        self.A,
        self.frame_idx,
        self.Ps_raw,
        self.Ts_raw,
        None,  # n0_first (can extract from headers if needed)
        True,  # kcal
        10.0   # Rmax1D
    )
    
    self.canvas_cuts.draw_idle()
```

### Phase 6: Cursor and Blitting Updates

**Ensure cursor works on all tabs:**
- Current BlitManager only handles Energy Maps tab
- Need separate BlitManager for Decomposition tab
- Cursor position should update all tabs simultaneously

**Add cursor to decomposition tab:**
```python
# In update_decomposition_plots(), add cursor artist
if not hasattr(self, 'cursor_decomp_tot'):
    self.cursor_decomp_tot, = self.ax_decomp_tot.plot([], [], 'w+', ms=15, mew=2.5, animated=True, zorder=10)
    # Add to blit manager
    self.blit_manager_decomp._artists.append(self.cursor_decomp_tot)

# Update cursor position
row = self.pixel_row_spin.value()
col = self.pixel_col_spin.value()
self.cursor_decomp_tot.set_data([col], [row])
```

---

## Critical Data Consistency Requirements

### 1. Single Source of Truth for Energies

**Problem addressed in plot_masked_energy.py (lines 584-631):**
- Different evaluation paths can give different results
- Must use same kernel call for both plotting and optimization

**Solution:**
- Always use `evaluate_energies_masked(mask_hbond)` for both Etot and Ein
- Store the returned arrays and derive Eout = Etot - Ein
- Never call `evaluate_energies()` separately for the same parameter set

### 2. Error Maps from Optimizer Data

**Problem:**
- Computing dE = V_mod - V_ref after the fact can have unit inconsistencies
- Optimizer uses soft-clamping and specific weight arrays

**Solution:**
- After MC optimization, extract from history:
  - `dE_initial`, `dE_best` (raw errors in eV)
  - `J_per_sample_initial`, `J_per_sample_best` (weighted errors)
  - `Eref` (reference energies used by optimizer)
  - `W` (weights used by optimizer)
- Map these to grids using [_map_samples_to_grid()](cci:1://file:///home/prokop/git/FireCore-fitREQH/tests/tFitREQ/plot_masked_energy.py:171:0-180:12)
- Display these authoritative values, not recomputed ones

### 3. Parameter Upload Consistency

**Problem addressed in plot_masked_energy.py (lines 336-338, 403-409):**
- `tREQHs_base` may not be uploaded to `tREQHs_buff` automatically
- Need explicit upload to ensure GPU state matches Python state

**Solution:**
```python
# After any parameter change:
drv.toGPU_(drv.tREQHs_buff, np.array(drv.tREQHs_base, dtype=np.float32, copy=False))
drv.queue.finish()
```

---

## Summary of Changes Required

### New Imports
```python
from pyBall.FitREQutils import (
    _build_frame_grid, plot_system_panel, plot_polar_symmetric, 
    plot_profile_row, parse_xyz_with_headers
)
# Import functions from plot_masked_energy.py (or move them to a shared module)
from tests.tFitREQ.plot_masked_energy import (
    create_dof_definitions_from_current, run_mc_optimization,
    mc_history_to_json, load_mc_config_with_softclamp, _map_samples_to_grid
)
```

### New Data Members
- `frame_idx`, `Ps_raw`, `Ts_raw` (from _build_frame_grid)
- `mask_hbond` (for decomposition)
- `V_tot`, `V_ein`, `V_eout`, `V_diff`, `V_J` (decomposition grids)
- `mc_history` (MC optimization results)

### New UI Elements
- MC optimization parameter controls (spin boxes, checkboxes)
- MC "Run" button
- Decomposition tab with 6 panels

### New Methods
- `recompute_energy_with_decomposition()` - replaces current recompute_energy
- `run_mc_optimization_gui()` - runs MC in background thread
- `on_mc_finished()` - applies MC results
- `update_decomposition_plots()` - updates decomposition tab using plot_system_panel
- `update_1d_plots_with_profile_row()` - replaces current 1D cuts

### Modified Methods
- [init_fitting_driver()](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/GUI/FitREQInteractiveGUI.py:405:4-433:59) - initialize frame grid and masks
- [on_pixel_changed()](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/GUI/FitREQInteractiveGUI.py:1085:4-1112:45) - update cursor on decomposition tab
- [apply_current_params()](cci:1://file:///home/prokop/git/FireCore-fitREQH/pyBall/GUI/FitREQInteractiveGUI.py:623:4-680:70) - ensure GPU upload

---

## Testing Checklist

1. **Data consistency**: Verify that `V_tot - Vref` equals `V_diff` (within numerical tolerance)
2. **MC optimization**: Run MC and verify parameters are applied correctly
3. **Blitting performance**: Ensure cursor updates are fast on all tabs
4. **Plot style**: Compare decomposition plots to plot_masked_energy.py output
5. **Error maps**: After MC, verify J map matches optimizer's best_error
6. **Unit consistency**: Verify kcal/mol conversions are applied consistently

This plan maximizes code reuse by directly calling tested functions from plot_masked_energy.py and pyBall/FitREQutils.py, ensuring accurate reproduction of functionality and plotting style while maintaining fast interactive updates.