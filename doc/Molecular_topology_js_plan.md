
# TODO

- [x] Use **explicit hydrogens** as normal atoms (with neighbors etc.).
- [x] Selection/reference policy: store **atom IDs** everywhere; only use indices as transient dense-array slots.
- [x] Implement `EditableMolecule` core (Atom/Bond/Fragment/Bounds) with swap-remove and stable IDs.
- [x] Implement **direct export** `EditableMolecule.exportToMoleculeSystem(ms)` to fill existing packed render system (pos/types/bonds/selection).
- [x] Wire editor/GUI operations to mutate `EditableMolecule`, then export to packed system for rendering.

---

# Status update (Dec 2025)

## What was implemented

- `EditableMolecule` is now the authoritative topology/geometry model:
  - stable atom IDs + swap-remove deletion
  - bonds store `aId/bId` as source of truth with cached indices via `topoVersion`
- Export boundary exists and is used by renderer:
  - `EditableMolecule.exportToMoleculeSystem(ms)` feeds the packed renderer buffer
- Attachment and polymer assembly (ported from `tests/tAttach`):
  - marker attach: `EditableMolecule.attachGroupByMarker(...)`
    - supports distinct marker pairs on backbone vs group (`groupMarkerX/groupMarkerY`)
  - direction attach: `EditableMolecule.attachParsedByDirection(...)`
  - polymer builder: `PolymerUtils.assemblePolymerFromTokens(...)`
    - repeat-vector sign is chosen (±`lvec[1]`) to minimize head↔tail distance (prevents long wrap-around join bonds)
- IO/debug support:
  - `EditableMolecule.parseMol2/parseXYZ`
  - `EditableMolecule.toXYZString/toMol2String` + GUI download buttons

## Still missing / planned next

- Connected components / fragments computation (beyond the current `Fragment` scaffolding)
- Ring/cycle detection + bridge detection
- Group/functional substitution beyond marker-based anchors (chemically-aware replacement, valence checks)

# Best design (after reading the ring/chain editor doc + [MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0))

Your doc is implicitly asking for **two different responsibilities** for the same molecule:

- **Editable / semantic representation** (IDs, graph edits, rings/chains, bridge queries, merging)
- **Packed / numeric buffers** for the renderer / GPU / simulation (fast iteration + GPU upload + WebGL simulation)

Trying to make *one* structure do both is what caused the premature SoA/typed-array design trap.

So the clean design is:

## Design decision: **Editable graph is AoS (Atom/Bond classes), export/pack is SoA (into existing MoleculeSystem)**
- **During editing**: work with [Atom](cci:2://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:49:0-66:1)/[Bond](cci:2://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:197:0-220:1) objects + stable IDs + adjacency.
- **When rendering/simulating**: export/pack directly into the already-existing `MoleculeSystem` buffers that `MoleculeRenderer` consumes.

In JavaScript, all `number` values are already **double precision**, so “use doubles internally” is automatically satisfied unless you explicitly put positions into `Float32Array`. So: keep `Atom.pos` as [Vec3](cci:2://file:///home/prokop/git/FireCore/web/common_js/Vec3.js:0:0-197:1) (numbers), and only downcast when exporting to GPU.

---

# Core requirements distilled from your doc

From [Molecular_topology_editor_with_explicit_rings.md](cci:7://file:///home/prokop/git/FireCore/doc/Molecular_topology_editor_with_explicit_rings.md:0:0-0:0) the “must support later” operations are:

- **Dynamic topology edits**:
  - add/remove atoms/bonds
  - merge/collapse atoms (graph contraction)
  - split/rotate subgraph around a bond if it’s a bridge
- **Higher-level structures**:
  - rings/chains as “super-particles” with `(center, radius)` and optionally AABB
  - atoms/bonds may belong to multiple rings, but atom should store only *one* “primary structure” for UI
- **Efficient adjacency**:
  - organic typical degree ≤ 4, but allow up to ~8

From [MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0) the key idea to keep is:
- dense arrays + swap-remove + adjacency stored **per atom** as small fixed-ish lists of incident bonds
- stable-ish `id` separate from `index`

---

# Proposed JS architecture

## 1) `EditableMolecule` (graph + chemistry + geometry)
This is your new “real” molecule editor core.

### `class Atom`
- **Identity**
  - `id: int` immutable unique ID (for selection, undo, external refs)
  - `i: int` current dense-array index (changes on swap-remove)
- **Chemistry**
  - `Z: int` atomic number (or symbol)
  - `charge: number` optional
  - `flags: int` (bitfield: isCap, isSelected, isFrozen, etc.)
- **Geometry**
  - `pos: Vec3` (store as object, no conversions)
- **Topology**
  - `bonds: int[]` list of incident bond indices (small array)
  - hydrogens are **normal atoms** (no special capping representation)
- **Structure membership**
  - `frag: int` primary fragment id (or `-1`)
  - `fragSlot: int` index within fragment ordering (optional)

### `class Bond`
- **Identity**
  - `id: int`, `i: int`
- **Endpoints (ID + cached index)**
  - `aId: int`, `bId: int` are the **source of truth**
  - `a: int`, `b: int` are cached **atom indices** for fast tight loops when topology is not dirty
- **Chemistry**
  - `order: number` (allow aromatic)
  - `type: int` (single/double/etc) if you prefer
- **Topology flags computed later**
  - `isBridge: bool`
  - `isRingEdge: bool`
- **Helper**
  - `other(ai)` returns the other endpoint

### `class Fragment` (unified ring/chain)
Matches your “no Ring vs Chain inheritance” preference.
- `id: int`
- `isClosed: bool`
- `atoms: int[]` ordered atom indices (node atoms)
- `bonds: int[]` ordered edge indices (optional but very useful)
- `bounds: Bounds` (stores both sphere + aabb cheaply)
- `updateBounds(mol)` recompute center/radius/aabb from member atom positions

### `class Bounds` (generic “bucket”)
Keep both AABB and sphere together (cheap, avoids “switching between two”).
- `min,max` (Vec3)
- `center` (Vec3)
- `radius` (number)
- `reset()`, `addPoint(p)`, `finalize()`
- intersection methods: `sphereSphere`, `aabbAabb` (later)

### `class EditableMolecule`
- `atoms: Atom[]` dense
- `bonds: Bond[]` dense
- `fragments: Fragment[]`
- `id2atom: Map<int,int>` id→index
- `id2bond: Map<int,int>` id→index
- `lastAtomId: int`, `lastBondId: int` (monotonic counters)
- dirty flags:
  - `dirtyTopo`, `dirtyGeom`, `dirtyFrags`, `dirtyExport`

#### ID generation policy (important)

- IDs are generated by incrementing counters (`++lastAtomId`, `++lastBondId`).
- IDs are **never reused** after deletion.
- If an atom/bond is deleted and recreated, it always gets a **new** ID.
- IDs are not random hashes (no collision probability).

#### Cached indices validity (performance contract)

To support `Bond{aId,bId}` + cached `Bond{a,b}` without constant recomputation:

- `id2atom` maps **atomId -> atomIndex** (atoms are stored canonically in dense array `atoms[]`)
- `id2atom` is equivalent of C++ `std::unordered_map` (hash map) but we can implement it either as:
  - JS `Map<int,int>` (direct equivalent, simple)
  - or as an `Array<int>`/"sparse array" mapping: `id2i[id]=index` (often faster in V8 if IDs are dense monotonic)
- maintain `topoVersion: int` (increment on any topology change: add/remove atom/bond, merge, etc.)
- each bond stores `topoVersionCached: int` (or `cacheVersion`)
- `bond.ensureIndices(mol)` updates `bond.a/b` from `aId/bId` using `id2atom` only if `bond.topoVersionCached !== mol.topoVersion`

This gives:
- **Fast path**: numeric loops over indices with no map lookups when nothing changed
- **Correctness**: IDs remain stable across swap-remove; cached indices are lazily refreshed

#### Topology lock mechanism (for index-only fast operations)

Some performance-sensitive code paths should operate purely on indices (no `Map` lookup) under a guarantee that indexing will not change.

Proposed mechanism:

- `lockTopology()` increments `lockDepth` and returns a `topoVersionLocked = topoVersion` snapshot
- while `lockDepth>0`, **topology edits are forbidden** (add/remove atom/bond, merge, swap-remove)
  - such edits should throw (fail loudly)
- `unlockTopology()` decrements `lockDepth`
- index-based kernels assert `topoVersion === topoVersionLocked` (or just `lockDepth>0`)

This makes "fast-index" code safe and explicit.

#### Swap-remove done correctly (important)
If you remove atom at index `k` by swapping with last:
- swap `atoms[k] = atoms[last]`, set moved atom’s `i=k`
- update `id2atom` for moved atom
- **fix all bonds incident to moved atom** to replace endpoint `last→k`
- remove all bonds incident to removed atom (each bond also swap-removed, which requires updating adjacency lists)

This is exactly the “cost” noted in your doc and in the earlier draft: **swap-remove in a graph requires endpoint fixups**. But it’s manageable and fast because degree is small.

---

## 2) Export to existing `MoleculeSystem` (renderer/GPU boundary)

There should be **no extra intermediate container** between the editor graph and the renderer.

Instead, `EditableMolecule` provides an explicit boundary function:

- `exportToMoleculeSystem(ms)`
  - fills `ms.pos`, `ms.types`, `ms.bonds` (and `ms.nAtoms`, `ms.neighborList`)
  - this is where we accept any unavoidable conversions (e.g. `Vec3` -> packed arrays)

### Key point
**Export is the only place** where we touch packed numeric buffers / strides. Everything else uses `Atom.pos` ([Vec3](cci:2://file:///home/prokop/git/FireCore/web/common_js/Vec3.js:0:0-197:1)) directly.

Export can be:
- full rebuild on demand (fine for <1000 atoms)
- incremental later (optional)

---

# Why this is the best fit for your stated goals

## Programming convenience
- Topology edits become straightforward and local (swap-remove + small fixups).
- Stable IDs are first-class for undo/redo, selection, and external tools.
- Rings/chains are separate objects (fragments) with cached bounds → fast picking + collision tests.

## Speed
- Editing workloads are dominated by **small-degree graph operations**, not bulk numeric loops.
- When you *do* need bulk loops (rendering, WebGL forcefield), you pack once and run on packed buffers (or GPU).
- You avoid the slow pattern you hated: `Float32Array stride` + creating [Vec3](cci:2://file:///home/prokop/git/FireCore/web/common_js/Vec3.js:0:0-197:1) for every operation.

## Precision
- Internal geometry is `number` (double).
- GPU export is float32 (almost unavoidable in WebGL), but that’s a conscious conversion at a boundary.

---

# What to do with current [web/molgui_web/js/MoleculeSystem.js](cci:7://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSystem.js:0:0-0:0)

Right now it is SoA (`Float32Array pos`, `Uint8 types`, `bonds: [i,j]`), and it’s already entangled with renderer expectations.

I recommend:

- **Keep the existing [MoleculeSystem](cci:2://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSystem.js:0:0-799:1) as the buffer container** used by `MoleculeRenderer`.
- Add a new file/class (e.g. `EditableMolecule.js`) implementing the design above.
- Add export/import functions (no new intermediate class):
  - `EditableMolecule.exportToMoleculeSystem(ms)`
  - optionally `EditableMolecule.importFromMoleculeSystem(ms)` for loading an already-packed system into the editor graph

This lets you evolve editor features without breaking rendering.

---

# Two clarification questions (to avoid a wrong fork)

1) **Do you want to represent capping atoms (H) explicitly as atoms**, or keep them implicit (`nH`/valence-derived) until export?  
- I’d start explicit (simpler topology + ring fusion), and later add “implicit H mode” for UX.

2) For selection/GUI: do you want to select by **atom ID** everywhere (recommended), even if internally arrays are index-based?  
- If yes: selection should store `id`s, not indices.

---

# Status
- **Completed:** read and extracted requirements from the topology editor doc + compared to [MMFFBuilderBase.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/MMFFBuilderBase.h:0:0-0:0).
- **Proposed:** a concrete architecture (Editable AoS + Packed SoA) aligned with your workflow and future algorithms.
- **Next step:** once you confirm the two clarifications above, you should switch to **Code mode** and I’ll implement the minimal class skeleton + packing adapter without touching unrelated code.