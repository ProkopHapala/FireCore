# Molecular Editor: Fragment Selection & Manipulation (Design Scratchpad)

Status: draft notes (2025-08-10)

This document captures the design for enabling whole-fragment selection and rigid-body manipulation (move/rotate; optional scale) in the molecular editor when simulation is off. It references existing code to maximize reuse and minimize duplication.

See also:
- `cpp/common/molecular/MMFFBuilder.h`
- `cpp/common_SDL/SDL2OGL/MolGUI.h`
- `cpp/common/molecular/MolWorld_sp3.h`
- `cpp/common_SDL/SDL2OGL/EditorGizmo.h`
- Docs: `MolGUI.h.md`, `EditorGizmo.h.md`, `MolecularDraw.h.md`, `Draw3D_Molecular.h.md`

## Goals

- Enable selecting entire fragments by click or box select in edit mode.
- Use `EditorGizmo` to move/rotate selected fragments as rigid bodies.
- Apply transforms to `MM::Builder` data when `bRunRelax == false` (simulation off).
- Reuse `startFragment()`/`finishFragment()` and existing selection/transform APIs.
- Avoid code duplication; keep selection and data ownership consistent.

## Relevant APIs and Data

- MMFFBuilder (`cpp/common/molecular/MMFFBuilder.h`)
  - `struct Fragment { Vec2i atomRange, confRange, bondRange, angRange, dihRange; Vec3d pos; Quat4d rot; uint32_t color; }`
  - `std::vector<Fragment> frags;`
  - `void startFragment()` / `int finishFragment(int i=-1)`.
  - Atom field `Atom.frag` (int) indicates fragment membership.
  - Selection (atoms): `std::unordered_set<int> selection;`
  - Transforms:
    - `void move_atoms(Vec3d d, int i0=0, int imax=-1)`
    - `void transform_atoms(Mat3d M, Vec3d orig_old=0, Vec3d orig_new=0, int i0=0, int imax=-1)`
    - `void rotate_atoms(double angle, Vec3d axis=Z, Vec3d orig_old=0, Vec3d orig_new=0, int i0=0, int imax=-1)`
  - Selection helpers:
    - `int selectRect(const Vec3d& p0, const Vec3d& p1, const Mat3d& rot)` — fills `builder.selection`.

- MolWorld_sp3 (`cpp/common/molecular/MolWorld_sp3.h`)
  - `int selectFragment(int ifrag)` — builds `selection` vector by `a.frag == ifrag`.
  - `int selectRect(...)` — world-space selection alternative.

- MolGUI (`cpp/common_SDL/SDL2OGL/MolGUI.h`)
  - State: `bool bRunRelax`, `bool bBuilderChanged`, `bool useGizmo`, `EditorGizmo gizmo`.
  - Event: `void eventMode_edit(const SDL_Event& event)`; calls `gizmo.onEvent(mouse_pix, event)` when `useGizmo`.
  - Selection routing: `void selectRect(const Vec3d& p0, const Vec3d& p1)` → `builder.selectRect(...)` when builder is viewed; otherwise `W->selectRect(...)`.
  - UI: `panel_Frags` calls `W->selectFragment(i)`.
  - Gizmo bind targets in `updateGUI()`: `gizmo.bindPoints(natoms, apos); gizmo.bindEdges(nbonds, bond2atom);` (these reflect current view).

- EditorGizmo (`cpp/common_SDL/SDL2OGL/EditorGizmo.h`)
  - Current internal selection: `unordered_map<int,int> selection` (point_index → group).
  - Transform application currently moves bound `points` directly:
    - `applyTranslation(const Vec3d& shift)` → `move(nsel, sel, points, shift)`.
    - `applyScaling(...)` present; rotation TBD.
  - Axis picking and drag logic available; pose/orig is `gizmo.pose`.

## Reuse Plan

- Use `builder.selection` to track selected atoms in builder view.
- Use atom field `Atom.frag` to extend selection from individual atoms to whole fragments.
- For actual transforms, prefer invoking `builder` methods over directly modifying bound `points` when `bRunRelax == false`.
- Leverage existing `panel_Frags` to select by fragment id explicitly.

## Missing Pieces / Proposed Additions (C++; not implemented yet)

- In `MM::Builder`:
  - Add selection extension helpers:
    - `int extendSelection_toFragments();`  // read current `selection` (atoms), collect their `frag`, then include all atoms with those fragment ids.
    - Alternatively: `int selectFragmentAtoms(const std::unordered_set<int>& frag_ids);`.
  - Add selection-based transforms:
    - `void move_selection(const Vec3d& d);`  // iterate `selection` and add `d` to `atoms[i].pos`.
    - `void rotate_selection(double ang, const Vec3d& axis, const Vec3d& origin);`
    - `void transform_selection(const Mat3d& M, const Vec3d& orig_old, const Vec3d& orig_new);`
  - Optional: per-fragment rigid transform using ranges for performance:
    - `void transform_fragment(int ifrag, const Mat3d& M, const Vec3d& orig_old, const Vec3d& orig_new);`  // apply to `[atomRange.x, atomRange.y)`.

- In `MolGUI`:
  - A fragment selection mode toggle (GUI + keybinding) that makes box/click selections extend to full fragments by calling `builder.extendSelection_toFragments()`.
  - A policy for where `EditorGizmo` writes changes:
    - When viewing builder and `!bRunRelax`: forward drag deltas to `builder.move_selection(...)` / `rotate_selection(...)`.
    - When viewing world and/or `bRunRelax`: keep existing behavior (or disable gizmo).

- In `EditorGizmo`:
  - Callback interface to apply transforms externally instead of modifying `points` directly:
    - `std::function<void(const Vec3d& shift)> onTranslateSelection;`
    - `std::function<void(const Vec3d& axis, double phi, const Vec3d& origin)> onRotateSelection;`
    - `std::function<void(const Vec3d& scale, const Vec3d& origin)> onScaleSelection;`
  - On drag end/increment, invoke callbacks; keep visual feedback by moving a temporary preview array or by re-fetching positions after callback.

## UI/UX Flow

- Modes:
  - __Edit Mode__: enable gizmo (`useGizmo=true`) when `!bRunRelax`.
  - __Fragment Select Toggle__: checkbox/button in GUI; keybinding (e.g., `F`).

- Selection:
  - Box select or click select atoms as today.
  - If Fragment Select Toggle is on: immediately call `builder.extendSelection_toFragments()` to include entire fragments.
  - `panel_Frags` continues to call `W->selectFragment(ifrag)`; when in builder view, a builder-level helper can mimic this via `Atom.frag`.

- Transform:
  - Drag gizmo in move/rotate modes.
  - Gizmo calls back into `MolGUI`/`builder` to apply transform to selected atoms as a rigid body.
  - Origin:
    - Default: use gizmo pose position (`gizmo.pose.pos`).
    - Optionally: compute center of mass (or average) of selected atoms to set gizmo origin on pick start.

## Data Ownership / Synchronization

- When `bRunRelax == false`, editing should modify builder (`MM::Builder`). Ensure current `natoms, apos` bound to gizmo reflect builder arrays.
- When toggling `bRunRelax` on, existing code updates world from builder if `bBuilderChanged` is true. Keep this contract.
- Avoid writing to `W->ff.apos` directly in this mode to prevent divergence.

## Edge Cases / Constraints

- __Cross-fragment bonds__: if present, rigid-body transforms of one fragment will strain bonds. Options:
  - Block such transforms (warn user), or
  - Auto-extend selection to include any bonded fragments.
- __Rotation__: `EditorGizmo` lacks rotation application at the moment. Implement with a well-defined origin and axis in builder space.
- __Scaling__: generally non-physical; keep off by default; consider only for debugging.

## Implementation Steps (planned)

1. Builder helpers:
   - `extendSelection_toFragments()` and `move_selection()`; optional `rotate_selection()`.
2. MolGUI integration:
   - Fragment selection toggle; call extension after selections.
   - Bind `EditorGizmo` callbacks to builder selection transforms when `!bRunRelax` and builder view.
3. EditorGizmo callbacks:
   - Add callable hooks; invoke on drag update/end.
   - For preview, keep using current bind to positions; since builder updates `apos`, redraw will reflect changes.
4. Rotation support:
   - Compute axis from gizmo axes and mouse; call `rotate_selection(...)` with origin `gizmo.pose.pos`.
5. Testing:
   - Load simple multi-fragment system; select by fragment; move; rotate; verify `Atom.frag` consistency; verify `bBuilderChanged` gating to world.

## Pseudocode Sketches (non-binding)

- Extend selection to fragments in builder:
```cpp
int Builder::extendSelection_toFragments(){
    std::unordered_set<int> frs;
    for(int ia : selection){ frs.insert(atoms[ia].frag); }
    for(int i=0;i<atoms.size();++i){ if(frs.contains(atoms[i].frag)) selection.insert(i); }
    return selection.size();
}
```

- Move selection in builder:
```cpp
void Builder::move_selection(const Vec3d& d){
    for(int ia : selection){ atoms[ia].pos.add(d); }
}
```

- EditorGizmo callback usage in MolGUI:
```cpp
// on init
gizmo.onTranslateSelection = [&](const Vec3d& d){ if(!bRunRelax) { W->builder.move_selection(d); bBuilderChanged=true; } };
```

## Open Questions

- Which object is bound to `natoms, apos` in `updateGUI()` under all view modes? Ensure it matches builder in edit mode.
- Should we add per-fragment cached COM/orientation (`Fragment.pos/rot`) updates on transform to keep metadata coherent?
- How to expose fragment selection toggle to the GUI layout concisely?

## TODOs (short list)

- [ ] Confirm `Atom.frag` is consistently assigned on load and after any atom reorder.
- [ ] Decide behavior for cross-fragment bonds.
- [ ] Add builder selection helpers and gizmo callbacks.
- [ ] Implement gizmo rotation application path.
- [ ] Wire `bBuilderChanged` flipping on transform.
