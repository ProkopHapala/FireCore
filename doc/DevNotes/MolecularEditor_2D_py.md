
# (2D) Molecule Editor in python

## Scope & Goals
- **Purpose** Build a PyQt5 desktop tool for 2D planar molecular sketching/editing while maintaining reusable 3D data via [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78).
- **Primary Users** Surface-chemistry researchers prototyping adsorbed molecules.
- **Key Outcomes**
  - Top-down drawing canvas with 2D manipulation.
  - Integration with existing molecular data structures.
  - Export/import via formats already supported in [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78).

## Functional Requirements
- **Core Editing**
  - Add/remove atoms through mouse clicks.
  - Drag atoms in 2D while preserving `z` (either fixed or snapping to plane).
  - Change atom element (H/C/N/O/F initially).
  - Add/remove/change bond order (single/double/triple). Visual cues for bond type.
  - Toggle Draw vs Move mode. Draw controls topology; Move adjusts coordinates only.
- **Selection & Feedback**
  - Click-select atoms/bonds; highlight current selection.
  - Show atom index/element label near cursor or in status bar.
  - Undo/redo stack for edit actions.
- **Data Handling**
  - Use [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78) to represent structure; ensure updates to bonds/coords synchronize.
  - Auto-infer initial bonds via [atomicUtils.findBondsNP()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:211:0-230:61) if loading existing data.
  - Optional lattice handling for surface orientation (fix z-plane, allow global shift).
- **Import/Export**
  - Load `.xyz` (positions only) with [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78); rebuild bonds via [atomicUtils.findBondsNP()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:211:0-230:61). See more in `pyBall/doc/atomicUtils.md`.
  - Load `.mol`/`.mol2` preserving supplied bond topology through [AtomicSystem.__init__](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:18:4-78:25). See more in `pyBall/doc/AtomicSystem.md`.
  - Save back to `.xyz`/`.mol`/`.mol2` via [AtomicSystem.saveXYZ()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:80:4-86:144), [save_mol()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:89:4-90:63), [save_mol2()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:92:4-93:63).
- **Visualization**
  - Matplotlib-based drawing embedded in PyQt5 widget.
  - Depict bonds as lines with thickness/style per order.
  - Display atom colors from element dataset; optional size scaling by covalent radius.
- **UI Layout**
  - Central canvas (Matplotlib FigureCanvas).
  - Side panel: element palette, bond order controls, mode toggle, snap options.
  - Toolbar/status bar for undo/redo, file operations, coordinate readout.

## Non-Functional Requirements
- **Performance** Responsive for molecules up to ~500 atoms; use NumPy for coordinate updates.
- **Extensibility** Modular layering to allow future 3D view, surface lattice tools, force-field previews.
- **Reliability** Fail loudly on invalid operations; preserve internal consistency via assertions.
- **Portability** Linux desktop focus (PyQt5 compatibility).
- **Minimal Dependencies** Only PyQt5, NumPy, Matplotlib, existing project modules.

## System Architecture

- **UI Layer (PyQt5)**
  - MainWindow: orchestrates menus, sidebars, mode state.
  - CanvasWidget: Matplotlib FigureCanvas handling draw/move events.
  - SidePanel: element palette, bond selectors, mode toggles.
- **Controller Layer**
  - InteractionController: interprets mouse events based on current mode.
  - CommandManager: implements undo/redo via command pattern (add atom, move atom, etc.).
  - SelectionManager: tracks selected atoms/bonds, handles highlighting.
- **Model Layer**
  - [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78): authoritative molecule representation ([pyBall/AtomicSystem.py](cci:7://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:0:0-0:0)).
  - `AtomicViewModel`: wraps [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78) for 2D projection, caches drawing data.
  - `ActionValidators`: ensures valence/bond constraints before modifications.
- **Utilities**
  - `atomicUtils` for bond detection, neighbor lookup.
  - `elements` for element metadata (color, radii).
  - `plotUtils` for reference drawing routines (bond/atom plotting).

## Interaction Design

- **Mouse Controls**
  - Draw Mode:
    - Left-click empty space: add atom of selected element.
    - Drag from existing atom: add bond to new/existing atom; choose bond order via palette.
    - Click bond: cycle bond order or remove (with modifier key).
    - Right-click atom/bond: context menu for delete/change element/order.
  - Move Mode:
    - Click-drag atom: updates `x,y`; `z` remains constant (optionally lock to plane).
    - Drag selection box (with modifier) for multi-select; move group with translation.
- **Keyboard Shortcuts**
  - `Ctrl+Z / Ctrl+Y` undo/redo.
  - `M` toggle move mode; `D` toggle draw mode.
  - Numeric keys to set bond order (1/2/3).
- **Selection Feedback**
  - Highlight atoms/bonds using Matplotlib artists (color/linewidth change).
  - Status bar shows coordinates, selected atom index (from `AtomicSystem.aux_labels`).

## Data Flow
1. **Load File**: [AtomicSystem(fname=...)](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78) populates atoms/bonds; `atomicUtils` builds neighbors.
2. **Canvas Render**: `AtomicViewModel` projects 3D to 2D (use `apos[:, :2]`).
3. **User Action**: Controller updates [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78) (e.g., [delete_atoms](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:719:4-730:69), append to `apos`).
4. **Consistency**: After topology edits, refresh bonds via [atomicUtils.findBondsNP](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:211:0-230:61) or direct updates; call [AtomicSystem.neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23).
5. **Render Update**: Matplotlib artists refreshed with new positions or connectivity.

## AtomicSystem Integration

- **State Arrays**: `AtomicSystem.apos`, `AtomicSystem.enames`, `AtomicSystem.bonds`, and `AtomicSystem.ngs` hold coordinates, element labels, bond pairs, and neighbor maps. See more in `pyBall/doc/AtomicSystem.md`.
- **Topology Utilities**: [AtomicSystem.findBonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29) and [atomicUtils.neigh_bonds()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:261:0-268:17) maintain connectivity after edits. See more in `pyBall/doc/atomicUtils.md`.
- **Editing Hooks**: [AtomicSystem.add_atom()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:131:4-133:80), [add_bond()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:133:4-134:34), and [delete_atoms()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:719:4-730:69) map directly to UI actions. See more in `pyBall/doc/AtomicSystem.md`.
- **Coordinate Ops**: [AtomicSystem.shift()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:327:4-335:38) and [rotate_ax()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:337:4-341:54) support multi-atom transformations in move mode.

## Additional Helpers for Future Features

- **Ring Detection**: [atomicUtils.find_cycles()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:54:0-100:72) (preceded by [atomicUtils.preprocess_graph()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:41:0-52:17)) can identify aromatic and surface rings for template placement. See more in `pyBall/doc/atomicUtils.md`.
- **Angle & Dihedral Maintenance**: [atomicUtils.findAngles()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:335:0-349:21), [atomicUtils.findDihedral()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:352:0-386:24), and `AtomicSystem.findAngles()` / `AtomicSystem.findDihedral()` (see `pyBall/doc/AtomicSystem.md`) support enforcing planar constraints and monitoring torsions.
- **Bond-Length Preservation**: [AtomicSystem.store_bond_lengths()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:195:4-204:48) and [restore_bond_length()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:206:4-215:34) provide undo-safe geometric constraints during drag operations.
- **Neighbor Queries**: [atomicUtils.neigh_atoms()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:271:0-277:17) and [AtomicSystem.getNeighsOfType()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:172:4-174:69) enable valence-aware editing and future force-field previews.
- **Hydrogen Bond & Functional Group Detection**: [atomicUtils.findHBondsNP()](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:233:0-260:46) and [AtomicSystem.find_groups()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:157:4-167:69) support contextual selection tools. See more in `pyBall/doc/atomicUtils_long.md`.

## Reusable Modules & Specific Functionality

- **[pyBall/AtomicSystem.py](cci:7://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:0:0-0:0)**
  - Initialization/loading ([AtomicSystem.__init__](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:18:4-78:25)).
  - Bond management ([findBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:135:4-139:29), [delete_atoms](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:719:4-730:69), [saveXYZ](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:1099:0-1102:16)).
  - Neighbor discovery ([neighs()](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:150:4-154:23), [findBondsOfAtom](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:144:4-148:84)).
  - Coordinate manipulation ([shift](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:327:4-335:38), [rotate_ax](cci:1://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:337:4-341:54)).
- **[pyBall/atomicUtils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:0:0-0:0)**
  - Bond detection ([findBondsNP](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:211:0-230:61)).
  - Neighbor and angle utilities ([neigh_bonds](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:261:0-268:17), [findAngles](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:335:0-349:21)).
  - Element-based radius lookup ([getAtomRadiusNP](cci:1://file:///home/prokop/git/FireCore/pyBall/atomicUtils.py:208:0-209:55)).
  - Sampling helpers for label positioning (optional).
- **[pyBall/elements.py](cci:7://file:///home/prokop/git/FireCore/pyBall/elements.py:0:0-0:0)**
  - Element metadata (`ELEMENT_DICT`, [getColor](cci:1://file:///home/prokop/git/FireCore/pyBall/elements.py:154:0-159:18)); use for palette and rendering.
- **[pyBall/plotUtils.py](cci:7://file:///home/prokop/git/FireCore/pyBall/plotUtils.py:0:0-0:0)**
  - Reference for Matplotlib plotting of atoms/bonds ([plotAtoms](cci:1://file:///home/prokop/git/FireCore/pyBall/plotUtils.py:255:0-271:68), [plotBonds](cci:1://file:///home/prokop/git/FireCore/pyBall/plotUtils.py:273:0-295:20), [plotSystem](cci:1://file:///home/prokop/git/FireCore/pyBall/plotUtils.py:315:0-342:37)).
  - Provides color/size conventions already aligned with element data.
- **[pyBall/OCL/MMparams.py](cci:7://file:///home/prokop/git/FireCore/pyBall/OCL/MMparams.py:0:0-0:0)**
  - Potentially read `ElementTypes.dat` and `AtomTypes.dat` for extended properties (bond radii, valence). Useful for valence checks or bond-length suggestions.
- **Data Files**
  - `cpp/common_resources/ElementTypes.dat`, `AtomTypes.dat` for expanded parameters if needed for validation or display.

## External References
- **MolView** ([GitHub](https://github.com/molview/molview)): Web-based molecule editor; inspect interaction patterns and UI flows.
- **ChemDoodle** ([GitHub](https://github.com/zachcp/chemdoodle)): Desktop chemical drawing app; reference bond rendering and element palettes.
- **ChemDraw / MarvinSketch** (commercial benchmarks): Inform iconography, command layout, and advanced features (templates, ring builders).

## Open Questions & Future Enhancements
- **3D Adjustments**: Provide optional `z` slider for individual atoms to tweak adsorption height.
- **Snapping/Constraints**: Surface lattice alignment, bond angle snapping.
- **Template Library**: Quick insert of common fragments (benzene, methyl).
- **Force Field Integration**: Hook with `pyBall/OCL/MMFF.py` to preview energies or optimize geometry.
- **Export to GPU Pipelines**: Ensure compatibility with OpenCL workflows (requires setting `bMMFF=True` when necessary).

## Summary
Prepared comprehensive design for a PyQt5-based planar molecule editor using existing `pyBall` infrastructure (notably [AtomicSystem](cci:2://file:///home/prokop/git/FireCore/pyBall/AtomicSystem.py:16:0-890:78), `atomicUtils`, `elements`, and `plotUtils`). Document captures requirements, architecture, interactions, and reuse strategy, referencing relevant open-source tools for UI inspiration. No code changes performed yet.