# Molecular Topology — Base Graph & Connectivity

> Cross-language audit of molecular graph representations, bond finding, neighbor lists, ring/bridge detection, and hybridization geometry. See [interactive codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18) for visual navigation.

## Quick Navigation

| I want to... | Go to |
|--------------|-------|
| Assign atom types (C_sp2, C_sp3, etc.) from bonding | [molecular_topology_types.md](molecular_topology_types.md) |
| Build/edit molecules interactively or compare GUIs | [molecular_topology_editors.md](molecular_topology_editors.md) |
| Understand the two-phase model (topology vs FF) | [Conceptual Overview](#conceptual-overview) below |
| Compare graph representations across C++/Python/JS | [Graph Representations](#graph-representations) |
| Find rings in a molecule | [Ring Detection](#ring-detection) |
| Find critical bonds (bridges) | [Bridge Finding](#bridge-finding) |

---

## Conceptual Overview: The Two-Phase Model

Every molecular simulation follows the same two-phase pattern:

```
Phase 1: Topology Building/Editing        Phase 2: Force Field Evaluation
├─ Dynamic structures                        ├─ Dense, flat arrays
├─ Object graphs (atoms, bonds, rings)       ├─ Neighbor tables, GPU buffers
├─ Type assignment (sp1/sp2/sp3)             ├─ Bond/angle parameter lookup
└─ Algorithmically complex                   └─ Performance-critical
```

**The conversion boundary** (Topology → FF) is the critical bottleneck. Each language implements its own:
- **C++**: `MMFFBuilder::toMMFFsp3_loc()` (implied in Builder)
- **Python**: `MMFF.toMMFFsp3_loc()` [@/home/prokop/git/FireCore/pyBall/OCL/MMFF.py:377](pyBall/OCL/MMFF.py:377)
- **JavaScript**: `MMFFLTopology.buildXPDBInputsFromMol()` [@/home/prokop/git/FireCore/web/molgui_webgpu/MMFFLTopology.js:31](web/molgui_webgpu/MMFFLTopology.js:31)

---

## Graph Representations

All three languages maintain the same duality: a **dynamic graph** for editing and a **static array** for force field evaluation.

| Language | Dynamic Graph (Editing) | Static Array (FF Evaluation) | Conversion |
|----------|------------------------|------------------------------|------------|
| **C++** | `MM::Builder` (integer-indexed arrays) | `MMFFsp3_loc` flat arrays | Builder → FF arrays |
| **Python** | `AtomicGraph` (object-based) + `AtomicSystem` (array-based) | `MMFF` packed arrays | `toMMFFsp3_loc()` |
| **JavaScript** | `EditableMolecule` (object ID maps) | `XPDB` packed buffers | `buildXPDBInputsFromMol()` |

### C++: Integer-Indexed Arrays

`MMFFBuilderBase.h` defines flat arrays:
- `atoms[i]` — position, type, conf indices
- `bonds[j]` — atom indices a/b, order, type
- `confs[k]` — `AtomConf` with `nSigma`, `nPi`, `nEPair`, `capping`

```cpp
struct AtomConf {
    int nSigma;   // sigma bonds
    int nPi;      // pi orbitals
    int nEPair;   // electron pairs
    int capping;  // capping atom index
    // Methods: addNeigh, fillIn_npi_sp3, updateNtot
};
```

Key feature: `AtomConf` separates bonding topology from atom identity. One atom can have multiple confs.

### Python: Dual Representations

Python has **both** object-based and array-based representations, with no clear guidance on when to use each:

**`AtomicGraph`** (object-based, stable IDs):
- `Atom`, `Bond`, `Ring` objects with `id` (no renumbering on deletion)
- `atom.bonds` → list of `Bond` objects
- Used by: `KekuleBackend`, `KekuleExplorerGUI`

**`AtomicSystem`** (array-based):
- `apos[natoms, 3]`, `atypes[natoms]`, `bonds[nbonds, 2]`
- `neighs` list built from bonds
- Used by: `MoleculeEditor2D`, `MMFF.toMMFFsp3_loc()`

**Recommendation needed**: Standardize on one canonical Python representation. `AtomicGraph` is better for editing; `AtomicSystem` is better for FF conversion.

### JavaScript: Object ID Maps

`EditableMolecule` uses stable object IDs with lazy index resolution:
- `atoms[]` array + `id2atom` Map
- `Bond.aId`, `Bond.bId` (stable IDs) → `Bond.a`, `Bond.b` (array indices cached via `topoVersion`)
- `lockTopology()` / `unlockTopology()` for batch mutations

```javascript
class Bond {
    aId: number;  // Stable ID
    bId: number;
    a: number;     // Cached array index (lazy, invalidated by topoVersion)
    b: number;
    topoVersionCached: number;
}
```

---

## Bond Finding

Distance-based bond detection from atomic positions using covalent radii.

| Language | Function | File |
|----------|----------|------|
| **C++** | `MM::BuilderBase::autoBonds()` | `MMFFBuilderBase.h:518` |
| **Python** | `atomicUtils.findBondsNP()` | `atomicUtils.py` |
| **Python** | `AtomicSystem.findBonds()` | `AtomicSystem.py` |
| **JavaScript** | `EditableMolecule.recalculateBonds()` | `EditableMolecule.js` |

**Algorithm**: For each atom pair (i,j), compute `d = |ri - rj|`. If `d < (Rcov_i + Rcov_j) * bondFactor`, create bond. Default `bondFactor = 1.1`.

**C++ details**:
- Uses `MMFFparams` covalent radii via `ElementType.Rcov`
- `autoBonds()` can be called after `insertAtoms()` to auto-generate connectivity

**Python details**:
- `findBondsNP()` uses NumPy vectorized distance computation
- `findHBonds()` for hydrogen bonds (shorter cutoff)

**JavaScript details**:
- `MMParams.bondCutoff2(z1, z2, ...)` computes squared cutoff
- `MMParams.bondLengthEstimate(zA, zB)` estimates l0 from radii

---

## Neighbor Lists

Building adjacency from bonds for force field evaluation.

| Language | Function | Purpose |
|----------|----------|---------|
| **C++** | `BuilderBase::makeNeighs()` | Build neighbor list from bonds |
| **Python** | `AtomicSystem.neighs()` | List of neighbor atom indices per atom |
| **Python** | `MMFF.make_back_neighs()` | Reverse neighbor lookup for force recoil |
| **JavaScript** | `MMFFLTopology.buildMMFFLTopology()` | Build `bondsAdj` adjacency array |

**Key concept: Back-neighbors**. In a directed neighbor list (i → j), the back-neighbor is the index of i in j's neighbor list. Required for force accumulation (when j feels force from i, i feels equal and opposite from j).

---

## Hybridization & VSEPR Geometry

Determining sp1/sp2/sp3 from bonding topology.

| Language | Concept | Implementation |
|----------|---------|---------------|
| **C++** | `AtomConf.npi`, `AtomConf.nEPair` | `MMFFBuilderBase.h:60` |
| **Python** | `npi_list`, `nep_list` | `MMFF.initAtomProperties()` `pyBall/OCL/MMFF.py:51` |
| **JavaScript** | `npiList`, pi-orbital directions | `MMFFLTopology.computePiOrientations()` `MMFFLTopology.js:248` |

**The VSEPR pipeline** (same in all languages):
1. Count sigma bonds (`nSigma` = number of bonded neighbors)
2. Count pi orbitals (`nPi` = from atom type / conjugation)
3. Count electron pairs (`nEPair` = `valence - nSigma - nPi`)
4. Total domains = `nSigma + nPi + nEPair` → geometry:
   - 2 domains → linear (sp1)
   - 3 domains → trigonal planar (sp2)
   - 4 domains → tetrahedral (sp3)

**JavaScript VSEPR helpers** (for geometry editing):
- `missingDirsVSEPR(vs, nMissing, totalDomains)` — compute missing bond directions
- `orthonormalBasisFromDir(dir)` — build orthonormal frame

---

## Ring Detection

| Language | Algorithm | File | Lines |
|----------|-----------|------|-------|
| **Python** | DFS cycle detection | `AtomicGraph.py` | ~392 |
| **C++** | Not implemented in core (external) | — | — |
| **JavaScript** | Not implemented | — | — |

**Python `AtomicGraph.detect_rings(max_ring_size=8)`**:
- Build adjacency dict from bond graph
- DFS from each start atom, tracking path
- Cycle found when neighbor == start and path length >= 3
- Creates `Ring` objects with ordered `atoms` and `bonds`

See codemap trace [6] for algorithm walkthrough.

---

## Bridge Finding

Critical bonds whose removal disconnects the graph.

| Language | Algorithm | File | Lines |
|----------|-----------|------|-------|
| **C++** | Tarjan's DFS | `LimitedGraph.h` | 181 |
| **Python** | Not implemented | — | — |
| **JavaScript** | Not implemented | — | — |

**C++ `LimitedGraph<T,M>::bridge()`**:
- Template class with compile-time max neighbors `M`
- `addEdge()`, `bridgeUtil()` (recursive DFS)
- Tracks `disc[u]` (discovery time) and `low[u]` (lowest reachable)
- Bridge condition: `low[v] > disc[u]`

See codemap trace [7] for algorithm walkthrough.

---

## Selection & Connected Components

| Language | API | File |
|----------|-----|------|
| **C++** | `std::unordered_set<int> selection` | `MMFFBuilder.h` |
| **Python** | `grow_selection()`, `shrink_selection()`, `select_all_connected()` | `AtomicSystem.py` |
| **JavaScript** | `Selection` class with banks | `web/common_js/Selection.js` |

**Python selection operations**:
- `grow_selection(mask, n=1)`: Expand by n bonded shells
- `shrink_selection(mask, n=1)`: Contract by n shells
- `select_all_connected(i0)`: BFS/DFS from atom i0

**JavaScript `Selection`**:
- `add/remove/toggle` with tombstones (`-1` in vec)
- `SelectionBanks` for multiple named selections
- `selectByPredicate(range, predicate)`: Functional selection

---

## Cross-Language Feature Matrix

| Feature | C++ | Python | JS |
|---------|-----|--------|-----|
| Add/remove atoms | `insertAtom`, `removeAtom` | `AtomicGraph.add_atom` | `EditableMolecule.addAtom` |
| Add/remove bonds | `insertBond` | `AtomicGraph.add_bond` | `EditableMolecule.addBond` |
| Auto-detect bonds | `autoBonds()` | `findBondsNP()` | `recalculateBonds()` |
| Neighbor lists | Built into Builder | `neighs()` | `bondsAdj` array |
| Back-neighbors | Yes | `make_back_neighs()` | Not explicit |
| Ring detection | No | `detect_rings()` | No |
| Bridge finding | `LimitedGraph::bridge` | No | No |
| Selection/grow | `std::unordered_set` | `grow_selection()` | `Selection.js` |
| Hybridization calc | `AtomConf.npi/nEPair` | `initAtomProperties()` | `computePiOrientations()` |
| **Crystal building** | ❌ Not in builder | ❌ **Missing** | ✅ `CrystalUtils.js` (1188 lines) |
| **CIF parsing** | ❌ | ❌ **Missing** | ✅ `parseCIF()`, `cifToCrystalData()` |
| **Lattice vectors** | ❌ | ❌ **Missing** | ✅ `latticeVectorsFromParams()` |
| **Symmetry operations** | ❌ | ❌ **Missing** | ✅ `applySymmetryOpsFracSites()` |
| **Cell replication** | ❌ | ❌ **Missing** | ✅ `genReplicatedCell()`, `genReplicatedCellSlab()` |
| **Slab cutting (HKL)** | ❌ | ❌ **Missing** | ✅ `genReplicatedCellSlab()` with reciprocal lattice |
| **Plane cutting** | ❌ | ❌ **Missing** | ✅ `genReplicatedCellCutPlanes()` |
| **Site deduplication** | ❌ | ❌ **Missing** | ✅ `dedupFracSitesByTolA()` (grid-bucket) |
| **Bonds across cells** | `autoBonds()` (single cell) | ❌ **Missing** | ✅ `_computeBasisBonds()` (27 neighbor offsets) |

**Critical gap:** Crystal building exists only in JavaScript. Python has no CIF parsing, no lattice vectors, no symmetry operations, no cell replication. This blocks the unified Python VisPy GUI goal. See [molecular_topology_editors.md](molecular_topology_editors.md) § "Crystal Building — Cross-Language Gap Analysis" for detailed function inventory and port plan.

---

## File Index

### C++
- `cpp/common/molecular/MMFFBuilderBase.h` (808 lines) — base topology structures
- `cpp/common/molecular/MMFFBuilder.h` (837 lines) — advanced topology ops
- `cpp/common/dataStructures/LimitedGraph.h` (181 lines) — graph + Tarjan bridges

### Python
- `pyBall/atomicUtils.py` (2186 lines) — bond/angle/dihedral utilities
- `pyBall/AtomicSystem.py` (1314 lines) — array-based atomic system
- `pyBall/AtomicGraph.py` (392 lines) — object-graph representation
- `pyBall/OCL/MMFF.py` (1107 lines) — topology → FF conversion
- **GAP:** No crystal building module exists in Python

### JavaScript
- `web/common_js/Selection.js` (168 lines) — generic selection class
- `web/common_js/MeshBuilder.js` (1688 lines) — mesh + selection building
- `web/molgui_webgpu/EditableMolecule.js` (1057 lines) — molecular topology editor
- `web/molgui_webgpu/MMFFLTopology.js` (826 lines) — topology → XPDB conversion
- `web/molgui_webgpu/MMParams.js` (524 lines) — parameter loading
- `web/molgui_webgpu/CrystalUtils.js` (1188 lines) — **crystal builder (no Python equivalent)**
- `web/molgui_webgpu/Nanocrystals.js` (631 lines) — nanocrystal generation from CIF

---

## See Also

- [Type Assignment](molecular_topology_types.md) — atom type assignment, parameter files, MMFF/UFF
- [Editors & Editing](molecular_topology_editors.md) — GUIs, advanced editing, crystal building gap analysis, consolidation roadmap
- [GUI Feature Audit](gui_audit.md) — detailed visualization & editor feature matrices
- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Interactive Codemap](https://windsurf.com/codemaps/692593e6-1efe-495f-bbf6-2ad291a285c9-fe86ab10a43f3d18) — visual navigation of topology code across languages
- `doc/topical_audit/intramolecular_forcefields.md` — force field evaluation phase
- Canonical default parameters: `tests/tUFF/data_UFF/{ElementTypes,AtomTypes,BondTypes,AngleTypes,DihedralTypes}.dat`
