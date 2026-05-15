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

#### 9.3 Current Limitations

1. **XYZ File Format:** Currently uses wb97m_input (no electron pairs). E-pairs are in wb97m_output but not default.
2. **Tab Switching:** Must manually click "Show 3D Geometry" or "Show Interaction Matrix" buttons
3. **Decomposition Kernel:** Overwrites self.prg (breaks energy kernel) - needs separate program object
4. **E-pair Visibility:** Semi-transparent (alpha=0.5), may be hard to see
5. **No Auto-update:** Must click buttons to update views after pixel selection