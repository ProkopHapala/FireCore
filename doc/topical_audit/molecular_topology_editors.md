# Molecular Topology — Editors & Advanced Editing

> Cross-language audit of interactive molecular editors: GUI architectures, advanced editing operations, and consolidation roadmap. See [base topology](molecular_topology.md) for graph representations, [type assignment](molecular_topology_types.md) for typing, and [interactive codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18).

## Quick Navigation

| I want to... | Go to |
|--------------|-------|
| Understand graph representations | [molecular_topology.md](molecular_topology.md) |
| Understand type assignment | [molecular_topology_types.md](molecular_topology_types.md) |
| See which editor to use | [Editor Decision Tree](#which-editor-should-i-use) |
| See what needs consolidation | [Redundancies & Roadmap](#redundancies--consolidation-roadmap) |

---

## All GUIs Compared

| GUI | Technology | Backend | Status | Maintained? |
|-----|-----------|---------|--------|-------------|
| **KekuleExplorerGUI** | Python / Vispy / PyQt5 | `KekuleBackend` + `AtomicGraph` | **Active** | Yes |
| **MoleculeEditor2D** | Python / Matplotlib / PyQt5 | `AtomicSystem` directly | Experimental | **No** |
| **molgui_web** | JS / WebGL | `EditableMolecule` + `MMParams` | **Legacy** | No |
| **molgui_webgpu** | JS / WebGPU | `EditableMolecule` + `MMParams` + `MMFFLTopology` | **Active** | Yes |

---

## Which Editor Should I Use?

```
Need to edit molecules in Python?
├── Need 3D visualization + hex grid? → KekuleExplorerGUI
└── Need simple 2D sketching? → MoleculeEditor2D (experimental)

Need to edit molecules in browser?
├── New project / production → molgui_webgpu
└── Legacy WebGL code → migrate to molgui_webgpu
```

**Recommendation**: Use `KekuleExplorerGUI` (Python) or `molgui_webgpu` (JS). `MoleculeEditor2D` and `molgui_web` are effectively deprecated.

---

## Editor Architecture Patterns

### Pattern A: Direct Manipulation (KekuleExplorerGUI)
- GUI widgets directly call backend mutation methods
- Backend emits signals for GUI refresh
- No command queue / undo stack

### Pattern B: Document Model (MoleculeEditor2D)
- Document wrapper around AtomicSystem
- Explicit `auto_bonds()`, `auto_types()` calls
- Grid snapping modes (triangular, square, none)

### Pattern C: Versioned Object Graph (EditableMolecule)
- `topoVersion` counter invalidates cached indices
- `lockTopology()` / `unlockTopology()` for batch ops
- `id2atom` Map for stable object identity

```javascript
// EditableMolecule.js
lockTopology() {
    this.lockDepth++;
    if (this.lockDepth === 1) this._topoVersionLocked = this.topoVersion;
}

addAtom(x, y, z, Z) {
    this._assertUnlocked("addAtom");
    this._touchTopo();  // increments topoVersion
}
```

---

## Python Editors

### KekuleExplorerGUI (Active)

**File**: `pyBall/KekuleExplorerGUI.py` (1413 lines)

- Vispy-based 3D scene with PyQt5 widgets
- Hexagonal grid snapping for graphene-like structures
- Atom type selector panel
- Auto hydrogen capping
- Passivation groups (N, NH, CH, H, O, C=O, C-OH)
- Bond visualization
- XYZ export

**Backend**: `KekuleBackend.py` (1975 lines)
- `AtomicGraph` as authoritative state
- `AtomicSystem` for rendering/export (synced via `_sync_sys()`)
- `_target_sigma()`: Calculate target sigma bonds from valence + hybridization
- Ring detection: `detect_geometry_rings()` → `graph.detect_rings()`

**Key method**: `add_atom(pos, ename, atype, pin=None, parent=None, subtype='')`
- `pin`: hex grid node key `(rx, ry)`
- `parent`: parent atom (for H capping)
- `subtype`: 'C_sp2', 'H_cap', etc.

### MoleculeEditor2D (Deprecated)

**File**: `pyBall/GUI/MoleculeEditor2D.py` (1602 lines)

- Matplotlib-based 2D canvas
- Operates directly on `AtomicSystem`
- Modes: draw, move, bond, rotate
- Element canonicalization
- Bond order management
- **Status**: Less maintained than KekuleExplorerGUI. Consider for removal or archival.

---

## JavaScript Editors

### molgui_webgpu (Active)

**File**: `web/molgui_webgpu/main.js` + supporting modules

- WebGPU compute shaders for real-time MD
- XPBD/MMFFL force fields
- Parameter inspection
- Editing + simulation in same interface

**Key modules**:
- `EditableMolecule.js` — topology editing
- `MMParams.js` — parameter loading
- `MMFFLTopology.js` — topology → XPDB packing
- `MoleculeRenderer.js` — WebGPU rendering
- `GUI.js` — UI components

### molgui_web (Legacy)

**File**: `web/molgui_web/` directory

- Original WebGL implementation
- Same `EditableMolecule.js` and `MMParams.js` duplicated
- **Status**: Deprecated, plan migration to WebGPU

---

## Advanced Editing Operations

### Insert Ring / Remove Ring

| Language | Function | Backend |
|----------|----------|---------|
| **Python** | `KekuleBackend.add_ring()`, `remove_ring()` | `AtomicGraph` |
| **C++** | Not explicit | — |
| **JavaScript** | Not explicit | — |

**Python hex grid**: `KekuleBackend` defines hexagonal grid coordinates (`pin = (rx, ry)`). Adding a ring places atoms at hex vertices.

### Collapse Bond / Split by Bond

| Language | Function | File |
|----------|----------|------|
| **C++** | `MM::Builder::splitByBond()` | `MMFFBuilder.h` |
| **Python** | Not explicit | — |
| **JavaScript** | Not explicit | — |

**C++**: `splitByBond(ib)` removes bond `ib` and reassigns atoms to separate fragments.

### Fragment Manipulation

| Language | Class/Function | Features |
|----------|---------------|----------|
| **C++** | `Fragment` struct | `atomIds`, `bondIds`, `bounds`, `isClosed` |
| **JavaScript** | `Fragment` class | `atomIds`, `bondIds`, `bounds`, `updateBounds()` |
| **Python** | Not explicit | — |

---

## Geometry Transforms

| Language | Translate | Rotate | Orient |
|----------|-----------|--------|--------|
| **C++** | `move_atoms(dshift)` | `rotate_atoms(angle, axis)` | `orient_atoms(fw, up)` |
| **Python** | `AtomicSystem.translate()` | `AtomicSystem.rotate()` | — |
| **JavaScript** | `EditableMolecule.translateAtoms()` | `EditableMolecule.rotateAtoms()` | — |

**Critical bug fixed** (see memory): JS `translateAtoms` and `rotateAtoms` were using atom IDs as array indices instead of looking up `id2atom` Map. Fixed by using `mol.getAtomIndex(id)`.

---

## Redundancies & Consolidation Roadmap

### Critical: EditableMolecule.js Duplicate

| File | Lines | Status |
|------|-------|--------|
| `web/molgui_web/js/EditableMolecule.js` | 1057 | **Legacy** |
| `web/molgui_webgpu/EditableMolecule.js` | 1057 | **Active** |

**Action**: Move shared code to `web/common_js/EditableMolecule.js`. Update both molgui_web and molgui_webgpu to import from common location. Remove duplicate files.

### Critical: MMParams.js Versions

| File | Lines | Features | Status |
|------|-------|----------|--------|
| `web/molgui_web/js/MMParams.js` | 466 | Basic parsing | **Legacy** |
| `web/molgui_webgpu/MMParams.js` | 524 | Logger, verbosity, QEq | **Active** |

**Action**: Use WebGPU version as canonical. Update molgui_web to import from molgui_webgpu or common location. Deprecate WebGL version.

### Moderate: Selection Systems

| Language | Implementation | Status |
|----------|---------------|--------|
| **C++** | `std::unordered_set<int>` in MMFFBuilder | Keep |
| **Python** | Set-based in AtomicSystem | Keep |
| **JavaScript** | `Selection.js` with banks | Keep |

**Rationale**: Different approaches for different use cases; no consolidation needed.

### Moderate: Graph Representations

| Language | Representations | Recommendation |
|----------|----------------|---------------|
| **Python** | `AtomicGraph` (object) + `AtomicSystem` (array) | Standardize on one |
| **C++** | Integer-indexed arrays only | Keep |
| **JavaScript** | Object ID maps only | Keep |

**Python decision needed**: `AtomicGraph` for editing, `AtomicSystem` for FF conversion. Document when to use each.

### Low: Parameter File Parsing

Three separate parsers for same `.dat` format (C++, Python, JS).

**Options**:
1. Keep separate (language-specific needs)
2. Add JSON export from one canonical parser
3. Document parity and keep as-is

---

## Crystal Building — Cross-Language Gap Analysis

### Overview

Crystal building is a **critical gap** in the Python ecosystem. The JavaScript `CrystalUtils.js` (1188 lines) provides a complete crystallography toolkit with no Python equivalent. This blocks the goal of a unified Python VisPy GUI.

### JavaScript Implementation (Active, Complete)

**`web/molgui_webgpu/CrystalUtils.js`** (1188 lines) — the canonical crystal builder:

| Function | Lines | Purpose |
|----------|-------|---------|
| `latticeVectorsFromParams(params)` | 229-255 | Build 3 lattice vectors from (a, b, c, α, β, γ) |
| `fracToCart(p, lvec, out)` | 257-264 | Fractional → Cartesian coordinate conversion |
| `cartToFrac(p, lvec, out)` | 266-279 | Cartesian → Fractional coordinate conversion |
| `parseLatticeText(txt)` | 281-293 | Parse 3-line text input as lattice vectors |
| `parseSitesTextXYZ(txt, mode)` | 308-327 | Parse atomic sites (fractional or Cartesian) |
| `parseSymOpsText(txt)` | 382-391 | Parse symmetry operation strings |
| `parseSymOpXYZ(opStr)` | 485-498 | Parse single sym op "x,y,z" → matrix + translation |
| `parseCIF(cifText)` | 418-466 | Full CIF parser (tags + loops) |
| `cifToCrystalData(cifText)` | 540-580 | Extract lattice, sites, symmetry from CIF |
| `applySymmetryOpsFracSites(sites, symOps)` | 500-528 | Apply symmetry operations to fractional sites |
| `dedupFracSitesByTolA(sites, lvec, tol)` | 104-157 | Grid-bucket deduplication of fractional sites |
| `dedupMolAtomsByTolA(mol, tol)` | 159-227 | Dedup atoms in EditableMolecule by distance |
| `cellDataFromFracSites(lvec, sites)` | 582-602 | Convert fractional sites → flat arrays (pos, types, charges) |
| `genReplicatedCell(params)` | 928-1026 | Replicate unit cell (bulk) with optional bond building |
| `genReplicatedCellSlab(params)` | 667-750 | Replicate with slab cut (HKL normal, cmin/cmax range) |
| `genReplicatedCellCutPlanes(params)` | 752-864 | Replicate with arbitrary plane cuts (multiple planes) |
| `reciprocalLattice(lvec)` | 1063-1083 | Compute reciprocal lattice vectors |
| `expandPlaneTemplates(templates, cSym)` | 11-29 | Expand {100}, {110}, {111} plane families |
| `_computeBasisBonds(lvec, basisPos, types, mmParams)` | 877-926 | Compute bonds within unit cell + 27 neighbor offsets |

**`web/molgui_webgpu/Nanocrystals.js`** (631 lines):
- `buildCrystalFromCIFText(cifText, opts)` — high-level CIF → crystal pipeline
- `buildPlanesFromTemplates(lvec, templates, ...)` — plane template expansion for slabs
- `loadMMParams(args)` — load FF parameters for bond building

**`web/molgui_webgpu/BuildersGUI.js`** (1130 lines):
- Unit cell editor GUI (lattice vectors, atoms, symmetry, dedup controls)
- CIF file loading → fill editor fields
- Crystal generation (bulk, slab, cut planes) with replication counts
- Preset crystal structures (NaCl, CaF2, diamond, etc.)
- Preview rendering with cell box, atom spheres, bond lines

### Python Implementation (Missing)

| Capability | JS Status | Python Status |
|-----------|-----------|---------------|
| Lattice vectors from params | ✅ Complete | ❌ **Missing** |
| Frac↔Cart conversion | ✅ Complete | ❌ **Missing** |
| CIF parsing | ✅ Complete | ❌ **Missing** (no Python CIF parser in pyBall) |
| Symmetry operations | ✅ Complete | ❌ **Missing** |
| Site deduplication | ✅ Grid-bucket | ❌ **Missing** |
| Bulk cell replication | ✅ Complete | ❌ **Missing** |
| Slab cutting (HKL) | ✅ Complete | ❌ **Missing** |
| Plane cutting | ✅ Complete | ❌ **Missing** |
| Bond building across cells | ✅ Complete | ❌ **Missing** |
| Reciprocal lattice | ✅ Complete | ❌ **Missing** (exists in `doc/Julia/EwaldGrid.jl` only) |
| Nanocrystal generation | ✅ Complete | ❌ **Missing** |
| Preset crystals | ✅ Complete | ❌ **Missing** |

**Note:** `KekuleBackend.py` (1975 lines) builds hexagonal graphene grids only — not general crystallography.

### C++ Implementation (Partial)

`MMFFBuilder` (`cpp/common/molecular/MMFFBuilderBase.h`, 808 lines) has:
- `insertAtoms(n, pos, types)` — add atoms to builder
- `autoBonds(rCut)` — distance-based bond finding
- No CIF parsing, no lattice vectors, no symmetry, no slab cutting

### Recommended Python Port Plan

Create `pyBall/crystal/CrystalBuilder.py` with these functions (mapping from JS):

```python
# Core crystallography
def lattice_vectors_from_params(a, b=None, c=None, alpha=90, beta=90, gamma=90) -> np.ndarray(3,3)
def frac_to_cart(frac, lvec) -> np.ndarray(N,3)
def cart_to_frac(cart, lvec) -> np.ndarray(N,3)
def reciprocal_lattice(lvec) -> np.ndarray(3,3)

# CIF parsing
def parse_cif(cif_text) -> dict  # tags + loops
def cif_to_crystal_data(cif_text) -> dict  # lattice, sites, sym_ops

# Symmetry
def parse_sym_op_xyz(op_str) -> (matrix, translation)
def apply_symmetry_ops_frac_sites(sites, sym_ops, tol=1e-6) -> list

# Deduplication
def dedup_frac_sites_by_tol(sites, lvec, tol=0.1) -> list

# Cell generation
def gen_replicated_cell(lvec, basis_pos, basis_types, n_rep=(1,1,1), ...) -> AtomicGraph
ndef gen_replicated_cell_slab(lvec, basis_pos, basis_types, n_rep, n_hat, cmin, cmax, ...) -> AtomicGraph
def gen_replicated_cell_cut_planes(lvec, basis_pos, basis_types, n_rep, planes, ...) -> AtomicGraph

# Bond building
def compute_basis_bonds(lvec, basis_pos, basis_types, mm_params, bond_tol=0.2) -> list
```

**Key design decisions:**
- Use NumPy arrays (not Vec3 objects) for positions — vectorized, cache-friendly
- Return `AtomicGraph` instances (not raw arrays) for GUI integration
- Use `mm_params` from existing Python FF parameter loading (`pyBall/OCL/UFFbuilder.py`)
- Fail loudly on invalid inputs (matching JS error behavior)

---

## Consolidation Checklist

- [ ] **Port `CrystalUtils.js` → Python `CrystalBuilder.py`** (CRITICAL — blocks unified GUI)
- [ ] Move `EditableMolecule.js` to `web/common_js/`
- [ ] Consolidate `MMParams.js` to single version
- [ ] Archive or remove `MoleculeEditor2D.py`
- [ ] Document `AtomicGraph` vs `AtomicSystem` usage in Python
- [ ] Plan molgui_web → molgui_webgpu migration
- [ ] Add cross-language parity tests for topology operations
- [ ] Port VSEPR capping from `EditableMolecule.js` → Python
- [ ] Create `crystal_building.md` audit document
- [ ] Build unified VisPy GUI combining KekuleExplorer + crystal builder + VSEPR

---

## File Index

### Python
- `pyBall/KekuleExplorerGUI.py` (1413 lines) — active 3D editor
- `pyBall/KekuleBackend.py` (1975 lines) — hex grid backend (graphene only, NOT general crystallography)
- `pyBall/GUI/MoleculeEditor2D.py` (1602 lines) — deprecated 2D editor
- `pyBall/AtomicGraph.py` (392 lines) — object-based molecular graph
- `pyBall/AtomicSystem.py` (1314 lines) — array-based molecular system

### JavaScript (Active)
- `web/molgui_webgpu/CrystalUtils.js` (1188 lines) — **crystal builder (no Python equivalent)**
- `web/molgui_webgpu/Nanocrystals.js` (631 lines) — nanocrystal generation from CIF
- `web/molgui_webgpu/BuildersGUI.js` (1130 lines) — crystal builder GUI
- `web/molgui_webgpu/EditableMolecule.js` (1057 lines) — topology editor
- `web/molgui_webgpu/MMFFLTopology.js` (826 lines) — XPDB packing
- `web/molgui_webgpu/MMParams.js` (524 lines) — parameters
- `web/molgui_webgpu/MoleculeRenderer.js` — WebGPU rendering
- `web/molgui_webgpu/GUI.js` — UI components
- `web/molgui_webgpu/Selection.js` — selection system

### JavaScript (Legacy)
- `web/molgui_web/js/EditableMolecule.js` (1057 lines) — **[DUPLICATE]**
- `web/molgui_web/js/MMParams.js` (466 lines) — **[LEGACY]**

### C++
- `cpp/common/molecular/MMFFBuilderBase.h` (808 lines) — FF builder with autoBonds (no crystallography)
- `cpp/apps/MolecularEditor/MolGUI.h` — C++ SDL-based GUI

---

## See Also

- [Base Topology](molecular_topology.md) — graph representations, bond finding, rings, bridges
- [Type Assignment](molecular_topology_types.md) — atom type assignment, parameter files
- [GUI Feature Audit](gui_audit.md) — detailed visualization & editor feature matrices, consolidation plan
- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Interactive Codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18) — visual navigation
- Canonical default parameters: `tests/tUFF/data_UFF/{ElementTypes,AtomTypes,BondTypes,AngleTypes,DihedralTypes}.dat`
