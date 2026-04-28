---
title: MoleculeEditor2D User Guide
---

# Overview
`pyBall/GUI/MoleculeEditor2D.py` provides a planar molecular editor built on PyQt5 and Matplotlib. It operates on `AtomicSystem` instances, enabling interactive atom placement, bond editing, selection manipulation, and topology export for other FireCore tools.

# Quick Start
- **Launch GUI**
```
python tests/tAttach/run_editor.py --molecule ../../cpp/common_resources/mol/porphirin.mol2
```
- **Command-line entry point**
```
python -m pyBall.GUI.MoleculeEditor2D path/to/molecule.mol2
```

# Core Components
- **`MoleculeDocument`** (`pyBall/GUI/MoleculeEditor2D.py`)
  - **`load(path)`**: Reads `.xyz`, `.mol`, or `.mol2` into `AtomicSystem`, normalizes element metadata, finds bonds, and caches neighbor lists.
  - **`save(path)`**: Persists the active system; format inferred from extension.
  - **`auto_bonds(Rcut=None, RvdwCut=None)`**: Rebuilds bond list via `atomicUtils.findBondsNP()` with configurable cutoffs.
  - **`add_atom(xy, element)` / `delete_atom(index)`**: Manage atomic coordinates and reindex bonds.
  - **`set_bond(i, j, order)` / `delete_bond(i, j)`**: Create or remove bonds while tracking bond order.
  - **`move_atom(index, xy)` / `move_atoms(indices, delta)`**: Translate atoms or selections in-plane.

- **`MoleculeCanvas`**
  - Matplotlib canvas handling draw/move/bond modes, selection visuals, keyboard shortcuts, and view framing.
  - **Interaction matrix**

    | Mode | LMB press/drag | LMB click on bond | RMB click | RMB drag | Wheel | Notes |
    |------|----------------|-------------------|-----------|----------|-------|-------|
    | **Draw** | Add atom at cursor; drag ignored | Cycle bond order 1→2→3→1 | Delete atom or bond under cursor | Rectangle selection (Shift=add, Ctrl=subtract) | Scroll zooms via `QDoubleSpinBox` or keyboard only | Selecting an existing atom with LMB primes it as first endpoint for bond creation when switching to Bond mode |
    | **Bond** | Start bond from picked atom; release on second atom to create/update bond with current order; clicking empty space cancels | Cycle bond order | Delete atom or bond | Rectangle selection | — | Pending atom highlighted through selection set |
    | **Move** | Drag entire selection if grabbed atom is already selected; otherwise drags single atom | Cycle bond order | Select atom/bond under cursor; toggles membership with Shift/Ctrl | Rectangle selection | — | Movement delta applied in document coordinates |

  - **Selection semantics**
    - **RMB single-click**: Select atom or bond; obeys modifier keys. Unmodified clicks replace current selection.
    - **RMB drag**: Creates rectangular marquee. Release finalizes according to replace/add/subtract mode determined by modifiers.
    - **Modifiers**: `Shift` adds to selection, `Ctrl` subtracts, none replaces. Applies to both clicks and drags.

  - **Keyboard controls**
    - **Arrow keys**: Translate selected atoms by `arrow_step` (default `0.1 Å`). No movement occurs without an active atom selection.
    - **Delete/Backspace**: Remove selected atoms; bonded entries are removed automatically. Selected bonds without endpoints in selection are deleted independently.
    - **`+` / `=` / numpad `+`**: Zoom in by multiplying `view_zoom` by `1.1`.
    - **`-` / numpad `-`**: Zoom out by dividing `view_zoom` by `1.1`.
    - **Focus**: `MoleculeEditor.keyPressEvent()` forwards key events even when UI widgets hold focus, so shortcuts always reach the canvas.

  - **View management**
    - **Zoom**: Controlled via keyboard shortcuts or the **Zoom** spin box. `view_zoom` is clamped to `[0.1, 10.0]` to avoid degenerate scales.
    - **Center**: **Center X/Y** spin boxes adjust `view_center`; updates immediately recenter the axes.
    - **Autoscale**: The **Autoscale view** button and `autoscale_view()` helper compute bounding box of all atoms, expand by configurable margin, and update view parameters.
    - **Canvas fill**: `subplots_adjust` and `ax.set_position([0,0,1,1])` remove Matplotlib padding so the drawing fills the available window area.

- **`ToolPanel`**
  - Groups controls into **Mode**, **Elements**, **Bond settings**, **Selection & move**, and **View** sections using `QGroupBox` + `QFormLayout`.
  - **Mode buttons** (`draw`, `bond`, `move`) toggle `MoleculeCanvas` mode.
  - **Element list** selects default element for new atoms.
  - **Bond settings** spin boxes configure `auto_bonds()` cutoffs (`Rcut`, `Rvdw ×`) and trigger recomputation.
  - **Selection & move** controls adjust `arrow_step` and pick radius multiplier.
  - **View** controls adjust zoom/center numerically and provide an “Autoscale view” button.

- **`MoleculeEditor`**
  - Main `QMainWindow` wiring canvas + tool sidebar, menus for file operations, and forwarding key events to canvas.

- **Utility entry point**
  - **`launch_editor(path=None)`**: Initializes `QApplication`, loads optional molecule, constructs `MoleculeEditor`, and enters the Qt event loop.

# Usage Tips
- **Auto bonds**: Tune `Rcut` (Å) and `Rvdw ×` multipliers before pressing **Auto bonds** to match desired connectivity.
- **Selections**: Hold **Shift** with RMB rectangle to add, **Ctrl** to subtract.
- **Bond order cycling**: LMB on a selected bond cycles 1→2→3→1 without reopening dialogs.
- **Keyboard nudging**: Reduce `Arrow step` for fine adjustments; keyboard moves entire selection.
- **Viewport**: Use zoom/center spin boxes for precise framing or press **Autoscale view** after large edits.

# Integration Points
- **`AtomicSystem`** (`pyBall/AtomicSystem.py`): All document operations reuse existing atomic data structures and neighbor generation.
- **`atomicUtils.findBondsNP()`** (`pyBall/atomicUtils.py`): Shared bond detection routine ensures consistency with other FireCore components.
- **Tests**: `tests/tAttach/run_editor.py` demonstrates launching the editor on sample molecules.

# Future Work
- [x] **Triangular grid snapping**: Introduce an optional triangular lattice overlay governing atom placement, drag endpoints, and rotation pivots. The grid spacing should derive from typical C–C bond lengths, with snapping tolerance adjustable from the tool panel. This enables rapid construction of hexagonal rings or 60°/30° motifs without manual tweaking. *Status*: Implemented via `MoleculeCanvas.set_grid_options()` with triangular/square patterns, snapping, and overlay controls.
- [x] **Selection rotation tools**: Implement rotation handles/UI allowing selected atoms to rotate about configurable axes. Default axis should match the active projection's normal (e.g., z-axis in `xy` view) with numeric entry for angle and arbitrary axis definition (two-point or vector input). *Status*: Delivered using `MODE_ROTATE`, manual/COG pivots, PageUp/PageDown steps, and `AtomicSystem.rotate_subset()`.
- [ ] **Projection plane switching**: Extend `MoleculeCanvas` to expose multiple orthographic projections (`xy`, `xz`, `yz`). Switching projections must remap controls (e.g., default rotation axis) and update labels while preserving the 3D positions stored in `AtomicSystem`.
- [ ] **Atomtatic orientation alignment** - Orient molecule along axis using selected atoms/bonds, or atumatic orientation using principal axes of inertia (PCA). 
- [x] **Whole-molecule selection shortcut**: Detect connectivity components in the current `MoleculeDocument` and add a single-click gesture (e.g., `Alt+LMB` or toolbar button) to select entire molecules/fragments. This will streamline translating or rotating complex assemblies.
- [ ] **Periodic boundary authoring**: Provide dialogs to define lattice vectors, toggle periodic boundary conditions, and replicate/delete atoms across cell boundaries. Edited geometry must remain compatible with downstream exporters that expect periodic metadata.
- [ ] **Substrate management**: Add helpers to import or parameterize crystalline substrates, expose lattice orientation controls, and lock substrate atoms while arranging adsorbates above the surface. Coordinate with the planned periodic tools for consistent lattice handling.
- [ ] **Force-field topology export**: Integrate with `pyBall/OCL/MMFF.py` and `pyBall/OCL/MMparams.py` to emit topology/parameter files after editing. The editor should validate valency, assign atom types, and generate data structures required for the OpenCL MMFF pipeline.
- [ ] **AFM imaging linkage**: Connect the editor with the AFM visualization module so current geometries can be forwarded for imaging simulation without intermediate file juggling. Support both live previews and batch exports.
- [ ] **Projective dynamics relaxation**: Implement an in-editor relaxation option inspired by `doc/ProjectiveDynamics.md`. Angle constraints (1–3 interactions) become auxiliary springs, enabling solution in a single step via direct sparse factorizations (Cholesky, LDLᵀ, or LU). Reuse the constraint assembly routines already documented, cache factorization for repeated solves, and fall back to iterative mat-vec if system sizes exceed direct-solve limits.
