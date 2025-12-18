
# Session Summary (molgui_web refactor + performance + interaction fixes)

## Achieved / Implemented Changes

### 1) **Authoritative model migration groundwork**
- **`EditableMolecule` is now the authoritative molecular model** (AoS: `atoms[]` with `pos: Vec3`, `Z`, stable `id`).
- Legacy `MoleculeSystem` usage was pushed into a **packed render buffer** only.
- **[MoleculeRenderer](cci:2://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeRenderer.js:38:0-166:1) now renders from [PackedMolecule](cci:2://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeRenderer.js:4:0-36:1)** (typed arrays + capacity mgmt), synced from `EditableMolecule` via `exportToMoleculeSystem()` *only when dirty flags require it*.

### 2) **Legacy IO cleanup**
- **`IO.js` was removed**.
- XYZ export / parsing and mol2 parsing were **integrated into `EditableMolecule.js`** (so IO lives with the data model).
- [main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:0:0-0:0) / [GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0) were adjusted to use `EditableMolecule` directly.

### 3) **Math refactor: Vec3 + Mat3**
- Introduced **`web/common_js/Mat3.js`** for clean matrix operations and to eliminate ad-hoc array math.
- `CrystalUtils.js` and `PolymerUtils.js` were refactored to use:
  - `Vec3` for vectors
  - `Mat3` for rotations / alignment

### 4) **GUI boilerplate reduction**
- Expanded **[GUIutils.js](cci:7://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:0:0-0:0)** with helpers like:
  - [textArea()](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:134:4-150:5)
  - [setSelectOptions()](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:63:4-79:5)
  - plus existing helpers ([row](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:30:4-30:107), [btn](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:32:4-32:159), [num](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:34:4-40:5), [textInput](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:41:4-47:5), [selectList](cci:1://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:51:4-61:5), etc.)
- Large portions of [GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0) were rewritten to **use one-liners**, preserving behavior while cutting DOM boilerplate.

### 5) **Builder panels extracted from GUI.js**
- Created **[web/molgui_web/js/BuildersGUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/BuildersGUI.js:0:0-0:0)** and moved:
  - **Builder: Substrate**
  - **Builder: Polymers**
  out of [GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0).
- [GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0) delegates to [this.buildersGUI.addSubstrateSection(sidebar)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/BuildersGUI.js:11:4-149:5) and [addPolymersSection(sidebar)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/BuildersGUI.js:151:4-382:5).

### 6) **Rendering + CPU performance: on-demand rendering**
- The continuous [animate()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:239:4-246:5) loop was replaced by a **render-on-demand mechanism**:
  - [MolGUIApp.requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:20:4-28:5) schedules one RAF and renders once.
- Added ability to switch modes:
  - **[MolGUIApp.setContinuousRender(true/false)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:30:4-40:5)**
  - GUI checkbox: **“Continuous Render (Animate)”**
- Key insight: even if data is on GPU, **Three.js render loop costs CPU** (scene traversal, state sorting, driver calls). On-demand mode eliminates idle cost.

### 7) **Gizmo + selection fixes (TransformControls integration)**
Major issues were found and fixed:

- **Gizmo initially became non-draggable** after on-demand rendering.
  - Fixed by explicitly calling [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5) on TransformControls `change` / `dragging-changed`.

- **Selection got broken when gizmo was visible**
  - Root cause: selection code was stealing pointer events or gizmo helpers caused false positives.
  - Fixes:
    - raycast-based gizmo detection ([isPointerOnGizmo](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/Editor.js:74:4-93:5))
    - handle filtering by name (`X,Y,Z,XY,YZ,XZ,XYZ,E`)
    - explicit flag `gizmoPointerActive`
    - use pointerdown listener in **capture phase**
    - clear flags to prevent “stuck” states

End result:
- **Gizmo works**
- **Selection works even with gizmo visible**

### 8) **Reliable redraw after model changes (critical UX fix)**
On-demand rendering caused some operations (substrate generation etc.) to update internal state but not show visually.

Fix:
- Introduced [GUI.requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5) helper (calls [window.app.requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5)).
- Ensured *model-changing operations* call [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5) after [renderer.update()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeRenderer.js:52:4-69:5):
  - substrate/polymer generation + attachment (in [BuildersGUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/BuildersGUI.js:0:0-0:0))
  - load XYZ / clear scene / view changes / label/axes UI (in [GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0))
  - editor operations already request renders in key places.

### 9) **Visual / UX improvements**
- Selection highlight became **semi-transparent gray/white** (alpha blending) instead of solid yellow.
  - Shader updated (`atom.glslf` got `uAlpha`)
  - [MeshRenderer](cci:2://file:///home/prokop/git/FireCore/web/common_js/MeshRenderer.js:2:0-344:1) sets alpha + blending states for selection mesh.
- Orthographic camera near/far range widened to avoid clipping artifacts.
- Gizmo default-enabled behavior fixed, and size reduced.

### 10) **Crystallography / substrate generation (CIF + bonds + Miller slab cut)**
- Added a minimal CIF pipeline in `web/molgui_web/js/CrystalUtils.js`:
  - `parseCIF()` / `cifToCrystalData()` and `genCrystalFromCIF()`
  - **Optional symmetry expansion** via separate function `applySymmetryOpsFracSites()` (used only if requested)
- `BuildersGUI` was kept thin:
  - CIF presets (`*-sym.cif` and `*-nosym.cif`) + file load
  - checkbox to toggle symmetry application

- Added **bond building using known bond lengths**:
  - `MMParams` now parses `BondTypes.dat` and exposes bond-length lookup
  - `genReplicatedCell(..., buildBonds=true)` builds bonds efficiently using basis + 26 neighbor images and replicates them
  - Added GUI toggle **"Build bonds (BondTypes)"**

- Implemented **efficient slab cutting by Miller indices**:
  - New `genReplicatedCellSlab()` generates atoms only if `cmin <= dot(r,nHat) <= cmax` (Å)
  - Uses **per-cell overlap pruning** (skip whole unit cells if slab cannot intersect)
  - `genCrystalFromCIF()` / `genCrystalFromMPJson()` accept `params.slab={hkl,cmin,cmax}`
  - Builders GUI has `Slab cut (cmin/cmax along n)` controls
  - Existing `Orient by Miller (h k l) -> z` provides optional rotation to keep surface in XY

Current limitation:
- **Slab cut is currently not compatible with the fast bond replication path** (`buildBonds`).
- Rationale: slab removes atoms, so pre-replicated periodic bond patterns break. A slab-aware bond builder is needed.

---

## Key Insights / Learnings (future mol_gui development)

### **1) On-demand rendering requires explicit “invalidate” strategy**
When not continuously animating, the app must have a reliable policy:
- any user-visible change must call [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5).

Best practice recommendation:
- Keep **one central `invalidate()`** (aka [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5)) and call it from:
  - Editor actions
  - GUI actions
  - any async operations (file load, library import, etc.)

### **2) Dirty flags are necessary but not sufficient**
- Dirty flags prevent costly buffer exports / uploads.
- But on-demand rendering still needs a redraw trigger.
- Correct architecture is:
  - `EditableMolecule` sets dirty flags on mutation
  - [MoleculeRenderer.update()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeRenderer.js:52:4-69:5) syncs buffers only if dirty
  - [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5) draws a frame only if needed

### **3) TransformControls vs custom pointer handling is tricky**
- TransformControls can swallow pointer events and uses invisible helper objects.
- Robust integration requires:
  - detect true handle hits (not helper planes)
  - prevent selection from starting on gizmo interaction
  - do not rely solely on `gizmo.dragging` (can be late)
  - use capture phase if needed

### **4) Rendering “alone” still costs CPU in Three.js**
Even without geometry changes:
- scene graph traversal
- uniform updates
- WebGL calls + driver overhead
So continuous 60 FPS rendering is inherently CPU-active.

### **5) Keep packed GPU buffer separate from edit model**
This is working well:
- AoS editing model for correctness and features
- SoA/typed-array packed structure for GPU efficiency

---

# Current TODOs / Unfinished Work

## High priority (next session)

- **[EditableMolecule integration completion]** (`id: 21`, in_progress)
  - implement missing legacy attachment functions currently stubbed in `EditableMolecule.js`:
    - `attachGroupByMarker(...)`
    - `attachParsedByDirection(...)`
  - confirm all Editor/GUI workflows are fully on `EditableMolecule` (no old model assumptions).

- **Rendering mode switch (on-demand default)** (`id: 43`, done)
  - explicit GUI checkbox exists: **“Continuous Render (Animate)”**
  - default: unchecked => on-demand rendering

- **Finish GUIutils consolidation** (`id: 38`, done)

## Medium priority

- **Crystal/Polymer generators API cleanup** (`id: 24`, medium)
  - ensure generator APIs and GUI call sites consistently use `Vec3`/`Mat3` (Miller alignment etc.)
  - unify generator input conventions (origin, lvec, offsets)

- **Slab-aware bond builder** (`id: 54`, medium)
  - Enable bonds for `genReplicatedCellSlab()` (currently disabled)
  - Options:
    - per-atom bond search restricted to within-slab atoms + 26 neighbor images, or
    - compute basis bonds + neighbor offsets then only instantiate bonds whose endpoints survive the cut
  - Must avoid `O(N^2)` for large slabs.

- **Bond order propagation** (`id: 55`, medium)
  - `BondTypes.dat` contains bond order and spring constant; current crystal builder creates bonds but does not store order/type.
  - Decide how bond order/type should be represented in `EditableMolecule.Bond` and in rendering.

## Low priority / deferred

- **IO integration correctness audit** (`id: 22`, low)
  - XYZ/mol2 import/export is integrated into `EditableMolecule`, but not exhaustively tested.
  - address only if issues appear.

- **Bond/atom z-buffer ordering investigation** (`id: 28`, low)
  - determine how impostor depth is written (quad vs sphere depth)
  - document current behavior and options (depth prepass, sphere-depth shader)

---

# Files Touched (major)
- **[web/molgui_web/js/main.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:0:0-0:0)**
  - on-demand rendering + continuous toggle
  - [requestRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:18:4-20:5), [setContinuousRender()](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/main.js:30:4-40:5)
- **[web/molgui_web/js/Editor.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/Editor.js:0:0-0:0)**
  - gizmo/selection event routing fixes
  - requestRender hooks
- **[web/molgui_web/js/GUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/GUI.js:0:0-0:0)**
  - GUIutils refactor + requestRender helper
  - continuous render checkbox
- **[web/molgui_web/js/BuildersGUI.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/BuildersGUI.js:0:0-0:0)** *(new)*
  - builder panels extracted
  - ensured render invalidation after builder actions
- **[web/molgui_web/js/MoleculeRenderer.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeRenderer.js:0:0-0:0)**
  - uses packed buffer; dirty-gated export
- **`web/molgui_web/js/EditableMolecule.js`**
  - authoritative model + IO (XYZ/mol2)
- **[web/common_js/GUIutils.js](cci:7://file:///home/prokop/git/FireCore/web/common_js/GUIutils.js:0:0-0:0)**
  - helper expansions
- **`web/common_js/Mat3.js`** *(new)*
- **`web/molgui_web/js/CrystalUtils.js`, `web/molgui_web/js/PolymerUtils.js`**
  - refactor toward `Vec3`/`Mat3`
- **Shaders + rendering**
  - `web/common_resources/shaders/atom.glslf`
  - [web/common_js/MeshRenderer.js](cci:7://file:///home/prokop/git/FireCore/web/common_js/MeshRenderer.js:0:0-0:0)

---

## Session completion status
- **Stabilized app usability** after on-demand rendering:
  - gizmo works
  - selection works
  - builder operations redraw reliably
- **Core refactor direction is established**:
  - EditableMolecule authoritative
  - packed renderer buffer
  - math cleaned via Vec3/Mat3
  - GUI extraction + GUIutils consolidation underway
- **Crystallography features are usable**:
  - CIF import (preset + file)
  - optional symmetry expansion
  - bonds for full periodic crystals (BondTypes)
  - Miller slab cut (cmin/cmax) with optional rotation to Z
  - remaining blocker: slab bonds