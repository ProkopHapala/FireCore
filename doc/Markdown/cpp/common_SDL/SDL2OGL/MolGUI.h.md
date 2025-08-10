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

---

## Non-bonded probing workflow

1. Configure probe parameters: `particle_REQ`, `particle_REQ2`, `particle_Lbond`. Attach to atoms via `particlePivots` if needed.
2. Open non-bond GUI (`nonBondGUI()` builds `panel_NonBondPlot`) and enable desired plots: `bDrawNonBondLines`, `bDrawNonBondGrid`, `bDrawParticles`.
3. Initialize/update with `tryPlotNonBond()`.
4. Line plots: call `plotNonBondLines()`; internally uses `plotNonBondLine()` helper over segments.
5. 2D grid: call `plotNonBondGrid()` (with optional `plotNonBondGridAxis()` for guides); evaluation uses `evalNonBondGrid2D()`.
6. Relax attached particles: `relaxNonBondParticles(double dt=0.2, double Fconv=1e-6, int niter=1000)`; visualize with `drawParticles()`.

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
