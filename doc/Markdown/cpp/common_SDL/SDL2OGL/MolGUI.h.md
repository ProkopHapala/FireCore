# MolGUI.h

This header defines the `MolGUI` class, a 3D molecular editor and visualization app built on SDL2/OpenGL. It manages molecular data binding, rendering (atoms, bonds, ESP, AFM, isosurfaces), non-bonded probing, GUI panels, and interactive editing with a transformation gizmo.

## Includes

- `globals.h`, `macroUtils.h`, `testUtils.h`, `IO_utils.h`
- SDL/OpenGL: `<SDL2/SDL.h>`, `<SDL2/SDL_opengl.h>`
- Draw utilities: `Draw3D.h`, `SDL_utils.h`, `Draw3D_Molecular.h`, `MolecularDraw.h`
- Simulation/FF: `MolWorld_sp3.h`, `MarchingCubes.h`, `AtomsInGrid.h`
- GUI and tools: `GUI.h`, `Console.h`, `EditorGizmo.h`, `SimplexRuler.h`, `AppSDL2OGL_3D.h`
- Misc: `raytrace.h`, `MethodDict.h`, `DipoleMap.h`, `<chrono>`

---

## Free Functions (helpers)

- `plotNonBondLine(const NBFF&, Quat4d, double, Vec3d, Vec3d, int, Vec3d up, bool bForce)` — Sample LJ+Q along a segment and draw energy/force strip.
- `evalNonBondGrid2D(const NBFF&, Quat4d, double, Vec2i ns, double*& Egrid, Vec3d p0, Vec3d a, Vec3d b, Vec3d up, bool bForce)` — Evaluate a 2D grid of energy/force; returns value range.
- `drawDipoleMapGrid(DipoleMap&, Vec2d sc, bool radial, bool azimuthal)` — Plot precomputed dipole/field map in radial/azimuthal strips.

---

## Class `MolGUI`

High-level 3D molecular GUI application. Inherits from `AppSDL2OGL_3D` to get camera, event loop, and 3D drawing.

### Construction

- `MolGUI(int& id, int WIDTH, int HEIGHT, MolWorld_sp3* W = 0)` — Create app window, optionally bind a world.

### Core virtuals (override)

- `draw()` — 3D scene rendering (molecule, FF grids, AFM, ESP, helpers).
- `drawHUD()` — 2D overlays: labels, GUI panels, status.
- `eventHandling(const SDL_Event&)` — Dispatch input to modes and gizmo.
- `keyStateHandling(const Uint8* keys)` — Continuous key state processing.

### Event modes and mouse

- `eventMode_default(const SDL_Event&)` — Base interaction mode.
- `eventMode_scan(const SDL_Event&)` — AFM/scan-specific interactions.
- `eventMode_edit(const SDL_Event&)` — Edit mode (selection, gizmo).
- `mouse_default(const SDL_Event&)` — Default mouse handler.

### Binding / data setup

- `bindMolWorld(MolWorld_sp3* W)` — Bind an MM world object.
- `bindMolecule(int natoms, int nnode, int nbonds, int* atypes, Vec3d* apos, Vec3d* fapos, Quat4d* REQs, Vec3d* pipos, Vec3d* fpipos, Vec2i* bond2atom, Vec3d* pbcShifts)` — Bind raw arrays.
- `bindMolecule(const MolWorld_sp3* W)` — Bind directly from `MolWorld_sp3`.
- `unBindMolecule()` — Release bound arrays.

### Rendering and visualization

- `drawSystem(Vec3i ixyz = Vec3iZero)` — Draw atoms, bonds, labels, neighbors.
- `drawBuilder(Vec3i ixyz = Vec3iZero)` — Draw builder primitives/bonds.
- `drawPi0s(float sc)` — Visualize pi-orbital reference directions.
- `renderSurfAtoms(Vec3i nPBC, bool bPointCross, float qsc, float Rsc, float Rsub)` — Draw surface atoms; returns vertex count.
- `renderGridFF(double isoVal, int isoSurfRenderType, double colorScale)` — Old grid FF isosurface renderer.
- `renderGridFF_new(double isoVal, int isoSurfRenderType, double colorScale, Quat4d REQ)` — New grid FF isosurface renderer with probe REQ.
- `renderESP(Quat4d REQ)` — Render electrostatic potential surface via point sampling.
- `Fz2df(int nxy, int izmin, int izmax, const Quat4f* F, float* dfout)` — Convert Fz slice stack to frequency shift map.
- `renderAFM(int iz, int offset)` — Render AFM image layer(s).
- `renderAFM_trjs(int di)` — Render AFM tip relaxation trajectories.

### Scanning and analysis

- `makeAFM(int iz = -1)` — Build AFM data (from grid and probe params).
- `lattice_scan(int n1, int n2, const Mat3d& dlvec)` — Perform scan over lattice translations.
- `makeBondLengths0()` — Precompute equilibrium bond lengths for coloring.
- `makeBondColoring(Vec2i types, Vec2d lrange, double*& clr, bool bNew)` — Color bonds by deviation from length range.

### Non-bonded probing UI

- `tryPlotNonBond()` — Initialize/update non-bond plots.
- `plotNonBondLines()` — Draw non-bond energy/force along predefined lines.
- `plotNonBondGrid()` — Draw 2D non-bond grid as shaded strip.
- `plotNonBondGridAxis()` — Draw axes/guides for non-bond grid.
- `relaxNonBondParticles(double dt, double Fconv, int niter)` — Integrate attached probe particles to equilibrium.
- `drawParticles()` — Draw auxiliary particles used in probing.
- `drawDipoleMap()` — Draw dipole/field map overlays.
- `showBonds()` — Debug bond visualization helper.
- `showAtomGrid(char* s, int ia, bool bDraw)` — Print/draw atom grid occupancy for atom.
- `Vec3d showNonBond(char* s, Vec2i b, bool bDraw)` — Evaluate/report non-bond between two atoms.
- `printMSystem(int isys, int perAtom, int na, int nvec, bool bNg, bool bNgC, bool bPos)` — Dump multi-system buffers.

### GUI and commands

- `initGUI()` — Create and place panels/widgets.
- `updateGUI()` — Sync GUI state from runtime variables.
- `clearGUI(int n)` — Remove last N panels/widgets.
- `initWiggets()` — Populate core widgets (panels, lists, toggles).
- `nonBondGUI()` — Build non-bond plotting GUI (checkboxes, controls).
- `initMemberOffsets()` — Register offsets of members for scripting/GUI binding.
- `initCommands()` — Populate `MethodDict<MolGUI>` command table.
- `setPanelAction(int ipanel, const char* name)` — Bind a panel command by name via `panelActions`.

### Important members (selected)

- __Camera/Gizmo__: `EditorGizmo gizmo`, `SimplexRuler ruler`, `Vec3d rotation_center/axis`, `double rotation_step`.
- __Runtime__: `int perFrame`, `Quat4f qCamera0`, `double cameraMoveSpeed`, `bool bConsole`.
- __Modes/flags__: `bool bDoMM`, `bool bRunRelax`, `bool bConverged`, `bool useGizmo`, `bool bHexDrawing`, `bool bDrawHexGrid`.
- __Panels__: `GUI gui`, `Console console`, `GUIPanel* Qpanel`, `DropDownList* panel_Frags`, `GUIPanel* panel_iMO`, `GUIPanel* panel_AFM`.
- __Non-bond GUI__: `CheckBoxList* panel_NBPlot`, multiple `MultiPanel*` for edit/plot/type/grid.
- __Probe particles__: `Quat4d particle_REQ/REQ2`, `double particle_Lbond`, `std::vector<Quat4d> particles/particles2`, `std::vector<int> particlePivots`.
- __Actions/binding__: `Dict<Action> panelActions`, `unordered_map<string,int> member_offsets_double/int`, `MethodDict<MolGUI> actions`.
- __Visualization params__: `float textSize`, `double ForceViewScale`, `double mm_Rsc`, `double mm_Rsub`, many `bView*` toggles (atoms, bonds, charges, labels, spheres, forces, bonds lengths, pis, substrate, boxes), `int isoSurfRenderType`.
- __Data bindings__: sizes `natoms,nnode,nbonds`, arrays `atypes`, `apos`, `fapos`, `pipos`, `fpipos`, `REQs`, `bond2atom`, `pbcShifts`, neighbor arrays `Quat4i* neighs/neighCell`, multi-system buffers `M_neighs/M_neighCell/M_apos`, `Constrains* constrs`, `double* bL0s`.
- __Graphics handles__: `int fontTex/fontTex3D`, display lists like `ogl_afm`, `ogl_afm_trj`, `ogl_esp`, `ogl_sph`, `ogl_mol`, `ogl_isosurf`, `ogl_surfatoms`, `ogl_MO`, `ogl_nonBond`, `ogl_Hbonds`, `ogl_trj`, `ogl_surf_scan`.
- __Debug__: `std::vector<Quat4f> debug_ps/debug_fs/debug_REQs`, `std::vector<Vec2i> bondsToShow`, `Vec3d* bondsToShow_shifts`.
- __AFM__: `GridShape afm_scan_grid`, `GridShape afm_ff_grid`, buffers `Quat4f* afm_ff/afm_Fout/afm_PPpos/afm_ps0`, `int afm_iz`, `int afm_nconv`, `bool afm_bDf`.
- __Maps__: `bool bDipoleMap`, `DipoleMap dipoleMap`.
- __World__: `MolWorld_sp3* W` (bound molecular system).
- __Misc__: `Mat3d dlvec/dlvec2` (lattice shifts), `Vec2f mouse_pix`, `double subs_iso`, `double z0_scan`.

---

## Notes

- Rendering uses immediate-mode OpenGL for simplicity; display lists are used for heavier geometry (isosurfaces, spheres).
- Many booleans gate visualization features (labels, charges, ESP, AFM). GUI panels toggle these at runtime.
- Non-bond probing supports energy vs. force views, 2D maps, and relaxation of attached probe particles.
- AFM tools support generating df/Fz maps and visualizing tip trajectories.

---

## Event modes

- `eventMode_default(const SDL_Event&)` — Basic camera and selection behavior.
- `eventMode_scan(const SDL_Event&)` — Interactions specialized for AFM/scan operations.
- `eventMode_edit(const SDL_Event&)` — Editing with `EditorGizmo` (translate/scale/rotate selected points).
- `mouse_default(const SDL_Event&)` — Default mouse handling when no specialized mode is active.

The active mode is tracked by `gui_mode` (`Gui_Mode::{base,edit,scan}`) and may switch via GUI commands or key handling in `eventHandling()`.

---

## View toggles (selected)

- Atoms/bonds: `mm_bAtoms`, `bViewBonds`, `bViewAtomSpheres`, `bViewBondLenghts`.
- Labels/props: `bViewAtomLabels`, `bViewAtomTypes`, `bViewMolCharges`, `bViewHBondCharges`, `textSize`.
- Forces/helpers: `bViewAtomForces`, `ForceViewScale`, `bViewGroupBoxes`, `bViewAxis`, `bViewCell`.
- Advanced: `bViewPis`, `bViewSubstrate`, `isoSurfRenderType`, `bDebug_scanSurfFF`.

These toggles are commonly exposed via `panel_Edit`, `panel_NonBondPlot`, and related GUI panels.

### Full view toggle list (from source)

- Builder and scene framing: `bViewBuilder`, `bViewAxis`, `bViewCell`.
- Atoms/bonds core: `mm_bAtoms`, `bViewBonds`, `bViewAtomSpheres`.
- Labels/annotations: `bViewAtomLabels`, `bViewBondLabels`, `bViewAtomTypes`.
- Charges and properties: `bViewMolCharges`, `bViewHBondCharges`.
- Forces/diagnostics: `bViewAtomForces`, `bViewGroupBoxes`.
- Advanced visuals: `bViewPis`, `bViewSubstrate`.
- Isosurface/scan flags: `isoSurfRenderType`, `bDebug_scanSurfFF`.

Note: some visuals (e.g., isosurfaces, AFM/ESP) are gated by runtime data presence in `W` and grid flags like `W->bGridFF`.

---

## GUI panels and controls

The following panels are created in `initWiggets()` and `nonBondGUI()` in `cpp/common_SDL/SDL2OGL/MolGUI.h`.

### Top controls

- `TableView("lattice")` — editable 3x3 matrix bound to `W->builder.lvec.{a,b,c}`; supports inline numeric edits via `GUITextInput`.
- `GUIPanel("Zoom:")` — slider in [5, 50]; updates `zoom`.
- `DropDownList("Pick Mode:")` — entries: `pick_atoms`, `pick_bonds`, `pick_angles` (affects selection mode).
- `DropDownList("Fragments:")` — selects fragment; command: `W->selectFragment(i)`.
- `DropDownList("View Side")` — entries: `Top/Bottom/Front/Back/Left/Right`; sets `qCamera` to pre-defined quaternions and updates `cam.rot`.

### Panel: Edit

- `print.builder.atoms` — `W->builder.assignPiFragments()`, `addCappingNeighborsToFragments()`, `printAtoms()`.
- `print.nonB` — `W->ffl.print_nonbonded()`.
- `print.Aconf` — `W->builder.printAtomConfs()`.
- `Sel.All` / `Sel.Inv` — select all/invert in builder or world (depending on `bViewBuilder`).
- `Sel.Cap` — `W->builder.selectCaping()` and merges into `W->selection`.
- `Add.CapHs` — `builder.addAllCapsByPi(H)`; sets `bBuilderChanged` if atoms added.
- `toCOG` — centers selection (`W->center(true)`); syncs builder if active.
- `toPCAxy` — aligns selection to PCA axes in XY; syncs builder if active.
- `save.xyz` — saves as `.xyz` (default `out.xyz`); builder or world based on `bViewBuilder`.
- `save.mol:` — `W->updateBuilderFromFF()` then `builder.saveMol` (default `out.mol`).
- `VdwRim` — generates/samples `dipoleMap`, calls `W->hideEPairs()`, saves `dipoleMap.xyz`, sets `bDipoleMap=true`.
- `scanSurfFF` / `SurfFF_view` / `SurfFF_scan` — invokes `scanSurfFF(...)` for a picked atom with parameters read from this panel (view, scan index, scaling).

### Panel: Gizmo

- `On/Off` — toggles `useGizmo`.
- `Mode:Trans` — `gizmo.mTrans = 'm'` (translate).
- `Mode:Rot` — `gizmo.mTrans = 'r'` (rotate).
- `AutoCOG` — toggles `gizmoAutoPivotCOG`.
- `ArcPick rmin` — sets `gizmo.rotPickRminFrac` in [0.5, 0.95] (rotation annulus inner fraction).
- `Pos=SelCOG` — sets `gizmo.pose.pos = W->center(false)` for current selection (uses builder selection when `bViewBuilder`).
- `SelFragOfPick` — selects fragment of the picked atom in world and updates builder selection if active.

### Panel: Run (CheckBoxList)

- `NonBond` → `W->bNonBonded`.
- `NonBondNG` → `W->bNonBondNeighs`.
- `GridFF` → `W->bGridFF`.
- `tricubic` → `W->bTricubic`.

### Panel: NBPlot (CheckBoxList)

- `minima` → `bDrawParticles`.
- `lines` → `bDrawNonBondLines`.
- `gridXY` → `bDrawNonBondGrid`.
- `hideEp` → `hideEp` (temporarily hides electron-pair pseudo-atoms while evaluating grids).

### Panel: PlotNonBond (MultiPanel)

- `Mode` — integer mode selector (0..2).
- `Ezoom` — energy zoom; effective scale uses `pow(10., value)`.
- `Rplot` — line length for 1D scans.
- `dstep` — sampling step for 1D scans.
- `Rdamp` — damping/ramping of REQ blending.
- `Rcut` — cutoff.
- `findHb` — button; runs `W->findHbonds_PBC(Rc, 0.01, 30°)` with slider `Rc`.

### Panel: GridXY (MultiPanel)

- `n` — grid resolution per axis.
- `size` — half-extent; area spans `[-size, +size]` in X and Y.
- `vmin` — log10 of value cap for visualization (used to derive `vmax`).
- `z_cut` — extra Z offset or slice parameter for display.

### Panel: BondLenghs (MultiPanel)

- `types:` — text input (e.g., `Si-Si`); parsed by `params.parseBondAtomTypes()`.
- `min.BL:` — minimum bond length for coloring threshold; triggers recalculation.
- `max.BL:` — maximum bond length for coloring threshold; triggers recalculation.

Action: updates bond coloring via `makeBondColoring(types, {min,max}, bL0s, true)` and sets
`bViewBuilder=false`, `bViewBondLenghts=true`, `bViewBondLabels=false`, `bViewAtomLabels=false`.

### Panel: TestType (probe REQ) and PickedType (picked atom REQ)

- `RvdW` — van der Waals radius.
- `EvdW` — van der Waals epsilon.
- `Charge` — point charge.
- `Hbond` — hydrogen-bond scaling; internally combined as `w *= sqrt(EvdW)`.

### Panel: Mol. Orb. and AFM

- `Mol. Orb.` — slider `which_MO` in [-5, 5]; calls `renderOrbital(HOMO + which_MO)`.
- `AFM iz` (only if GPU world) — integer `afm_iz` in [0, 20]; calls `makeAFM(afm_iz)`.

---

## Continuous key-state controls

Handled in `keyStateHandling(const Uint8* keys)`; disabled when `bConsole` or `gui.bTextEvents` are true.

- Camera translation (world-relative):
  - `Left/Right` arrows → `cam.pos += cam.rot.a * (±cameraMoveSpeed)`.
  - `Up/Down` arrows → `cam.pos += cam.rot.b * (±cameraMoveSpeed)`.
- Molecule shift (numeric keypad):
  - `KP_4/KP_6` → `W->nbmol.shift({∓0.1, 0, 0})`.
  - `KP_8/KP_2` → `W->nbmol.shift({0, ±0.1, 0})`.
  - `KP_7/KP_9` → `W->nbmol.shift({0, 0, ±0.1})`.

Notes:
- Commented scancode handlers for `WASD/QE` and camera pitch/yaw exist in code as reference but are currently inactive.
- These continuous controls apply in all GUI modes (`base`, `edit`, `scan`).

---

## Discrete keybindings by mode

All modes: mouse wheel zooms in/out.

### Default mode (`eventMode_default`)

__Console__
- `` ` `` Toggle console

__Selection__
- `Delete` Delete selected atoms and clear selections

__Lattice vectors__
- `,` Add `dlvec`
- `.` Subtract `dlvec`
- `0` Add `dlvec2`
- `9` Subtract `dlvec2`

__AFM slice index__
- `;` Increment `afm_iz` and call `renderAFM()`
- `'` Decrement `afm_iz` and call `renderAFM()`
- `KP_MULTIPLY` Increment `afm_iz` and call `renderAFM()`
- `KP_DIVIDE` Decrement `afm_iz` and call `renderAFM()`
- `(` Increment `afm_iz` and call `renderAFM()`
- `)` Decrement `afm_iz` and call `renderAFM()`

__Save/Export__
- `s` Save `snapshot.xyz` (with atomic numbers)
- `p` Save PNG screenshot, SVGs, and tiled XYZs

__Charges/AFM__
- `c` Auto-assign charges via `W->autoCharges()`
- `v` Build AFM via `makeAFM()` (GPU world only)

__Replica/system__
- `[` Switch to previous system replica
- `]` Switch to next system replica

__Gizmo__
- `g` Toggle `useGizmo`

__View toggles__
- `a` Toggle atom spheres
- `l` Toggle atom labels
- `t` Toggle atom types
- `b` Toggle bond labels
- `q` Toggle molecule charges
- `h` Toggle H-bond charges
- `f` Toggle atom forces
- `w` Toggle substrate
- `i` Toggle bond-length coloring
- `j` Toggle builder view
- `k` Toggle fragment coloring

__Camera__
- `KP_0` Reset camera orientation to `qCamera0`

__Method__
- `m` Switch computational method

__Run/relax__
- `Space` Toggle `bRunRelax`; sync builder/FF as needed; re-render MO when stopping

__Actions__
- any other key Dispatch via `actions.actionDispatch(this, key)`

__Utilities__
- `o` Run `lattice_scan(20,20, …)` (heavy; debug/analysis)
- `u` Upload population from `population.xyz`

Mouse (with gizmo, builder view only):
- MMB sets gizmo pivot to nearest builder atom under cursor.
- When rotating with LMB and gizmo engaged, events are consumed to avoid deselection.

### Scan mode (`eventMode_scan`)

__Lattice vectors__
- `,` Add `dlvec`
- `.` Subtract `dlvec`

__Probe rotation__
- `a` Rotate tip around +Z
- `d` Rotate tip around −Z

__Probe scale__
- `w` Scale tip vector up
- `s` Scale tip vector down

__Scan height__
- `KP_PLUS` Increase `z0_scan`
- `KP_MINUS` Decrease `z0_scan`

__Run/relax__
- `Space` Toggle `bRunRelax`

### Edit mode (`eventMode_edit`)

__Lattice transform__
- `t` Apply `affineTransform(ff.apos, builder.lvec -> new_lvec)`, update builder PBC, then swap `lvec`/`new_lvec`

__Selection helper__
- `i` Call `selectShorterSegment(ray0, cam.rot.c)`

---

## Non-bonded probing workflow

1. Open non-bond GUI: `nonBondGUI()` builds `panel_NBPlot`, `panel_PlotNonBond`, and `panel_GridXY`. Enable desired plots (`bDrawNonBondLines`, `bDrawNonBondGrid`, `bDrawParticles`) and set parameters (Mode, Ezoom, Rplot, dstep, n, size, vmin, z_cut). Configure probe REQ in `TestType`/`PickedType` panels if needed.
2. Initialize/update plots with `tryPlotNonBond()`.
3. Line plots: call `plotNonBondLines()`; internally uses `plotNonBondLine()` helper over segments.
4. 2D grid: call `plotNonBondGrid()` (use `plotNonBondGridAxis()` for guides); evaluation uses `evalNonBondGrid2D()`.
5. Relax attached particles: `relaxNonBondParticles(double dt=0.2, double Fconv=1e-6, int niter=1000)`; visualize with `drawParticles()`.

Tip: Use `hideEp` to suppress electrostatic overlays when focusing on LJ-only behavior.

---

## AFM and ESP rendering

- Grid isosurfaces: `renderGridFF(isoVal, isoSurfRenderType, colorScale)` or `renderGridFF_new(isoVal, isoSurfRenderType, colorScale, REQ)` where `REQ` defines probe size/charge; `colorScale` tunes coloring intensity.
- ESP shells: `renderESP(REQ)` samples points on shells and colors by electrostatic potential.
- AFM synthesis: `makeAFM(int iz=-1)` prepares force slices; `Fz2df()` converts Fz to df; `renderAFM(iz, offset)` and `renderAFM_trjs(di)` display images and tip trajectories.

Parameter hints:
- `isoVal` — isosurface value (typ. small positive).
- `isoSurfRenderType` — selects shading/style (implementation-specific).
- `colorScale` — scales color mapping magnitude.

---

## See also

- Source header: `cpp/common_SDL/SDL2OGL/MolGUI.h`
- Related docs: [EditorGizmo.h.md](EditorGizmo.h.md), [Draw3D_Molecular.h.md](Draw3D_Molecular.h.md), [MolecularDraw.h.md](MolecularDraw.h.md)
