# `EditableMolecule.js`

## Motivation / core ideas

`EditableMolecule` is the authoritative **editable molecular graph** used by MolGUI Web.

It is designed around a few practical constraints:

- **Stable identity under dynamic edits**: atoms/bonds have monotonic `id` values; dense-array indices can change due to swap-remove deletions.
- **Fast tight loops when possible**: topology operations can run on dense indices (`atoms[]`, `bonds[]`), while `id2atom`/`id2bond` provide robust lookup when topology changes.
- **Fail loudly**: operations that require invariants (e.g. locked topology) throw on violations.
- **Chemistry from `MMParams`**: passivation and typing use `AtomTypes.dat`/`ElementTypes.dat` via `MMParams`, not hard-coded tables.
- **Adjacency is stored incrementally**: `Atom.bonds[]` is the adjacency list; it is maintained by `addBond/removeBond*` and is used for selection queries and local chemistry logic.

This file also contains small geometry helpers (Vec3/Mat3 frames), a VSEPR-like completion routine for missing directions, and parsers/exporters for XYZ/MOL2.

---

## Data classes

- **`Atom`**: editable atom node with stable `id`, current dense index `i`, element `Z`, optional atom type index `atype`, position `pos`, and adjacency `bonds` (incident bond indices).
- **`Bond`**: editable bond edge with stable `id`, dense index `i`, endpoints by stable IDs (`aId/bId`) plus cached dense indices (`a/b`) guarded by `topoVersionCached`.
- **`Bounds`**: combined AABB + bounding sphere accumulator used for cheap spatial bounds of atom groups/fragments.
- **`Fragment`**: ring/chain-like grouping (`isClosed`) storing member `atomIds`/`bondIds` and cached `bounds`.

---

## Internal helpers (file-scope)

- **`_bondCut2()`**: compute squared bond cutoff distance using covalent radii from `MMParams` (with stats counters).
- **`_maxRcovFromMolAtoms()`**: find maximum covalent radius in the molecule (used to set AABB overlap margins).
- **`_ensureMMTypeIndex()`**: lazily build atomType name↔index mapping on `mmParams` (monkey-patched fields).
- **`_atomTypeNameToIndex()` / `_atomTypeIndexToName()`**: resolve atomType names and indices using the lazy mapping.
- **`_resolveTypeOrElementToAtomType()`**: interpret a token as atomType or element and return `{name, atype, iZ}`.
- **`_parseCountSet()`**: parse integer count sets like `{1,2}` for selection queries.
- **`_compileTokenSetToMatcher()`**: compile a token-set like `N|C` or `C_3|O_OH` into a fast predicate.
- **`_compileSelectQuery()`**: parse the AND-only query DSL into compiled matchers + count constraints.
- **`_getAtomTypeForAtom()`**: resolve the effective atom type for an atom (prefers `atom.atype`, falls back to element name).
- **`_bondLengthEst()`**: estimate bond length from covalent radii.
- **`_orthonormalBasisFromDir()`**: build a stable orthonormal basis around a direction.
- **`_missingDirsVSEPR()`**: VSEPR-like completion of missing directions given existing neighbor directions and desired domain count.

---

## `EditableMolecule` (core editable model)

### Lifecycle / state

- **`constructor()`**: initializes dense atom/bond storage, stable ID maps, selection set, dirty flags, and topology versioning.
- **`nAtoms`**: number of atoms (`atoms.length`).
- **`_touchTopo()`**: mark topology dirty and increment `topoVersion`.
- **`_touchGeom()`**: mark geometry dirty.
- **`_assertUnlocked()`**: enforce topology lock (fail loudly).
- **`lockTopology()` / `unlockTopology()` / `assertLocked()`**: explicit lock to guarantee indices don’t change during index-only kernels.

### ID/index lookup

- **`getAtomIndex(id)`**: stable atom ID → current dense index.
- **`getBondIndex(id)`**: stable bond ID → current dense index.

### Atom typing / chemistry helpers

- **`setAtomTypeByName(id, typeName, mmParams)`**: assign atom type by name (or element token) and update `Z` if provided.
- **`addCappingAtoms(mmParams, cap='H', opts)`**: add missing sigma-bond caps using atom-type valence + VSEPR-like direction completion.
- **`addExplicitEPairs(mmParams, opts)`**: add explicit lone-pair dummy atoms (`Z==200`) using `AtomTypes.dat` (`nepair`, `epair_name`).

### Geometry edits

- **`addAtom()` / `addAtomZ()`**: append a new atom (stable ID) to the dense array.
- **`setAtomPosById(id,x,y,z)`**: update atom position and mark geometry dirty.

### Topology edits (bonds)

- **`addBond(aId,bId,order,type)`**: create a bond and update both endpoint adjacency lists.
- **`removeBondByIndex(ib)`**: swap-remove a bond and patch adjacency lists.
- **`removeBondById(id)`**: remove a bond by stable ID.

### Selection

- **`selectAtom(id, mode)`**: update `selection` set (`replace`/`add`/`subtract`).
- **`select(idOrIndex, mode)`**: compatibility wrapper accepting either stable ID or dense index.
- **`clearSelection()`**: clear selection.

### Selection by chemical environment (query DSL)

- **`compileSelectQuery(q, mmParams)`**: compile query string once into matchers and count constraints.
- **`applySelectQuery(compiled, opts)`**: apply compiled query using `Atom.bonds[]` adjacency; supports `mode=replace|add|subtract`.

### Deletion / clearing

- **`deleteSelectedAtoms()`**: delete all selected atoms by stable ID.
- **`removeAtomByIndex(i)`**: swap-remove an atom and remove all incident bonds.
- **`removeAtomById(id)`**: remove an atom by stable ID.
- **`clear()`**: remove all atoms, bonds, fragments, and selection.

### Import / append

- **`addAtomsFromArrays(pos3, types1)`**: bulk add atoms from packed arrays.
- **`updateNeighborList()`**: no-op placeholder (adjacency is maintained incrementally).
- **`appendParsedSystem(parsed, opts)`**: append a parsed system (XYZ/MOL2-like) with optional translation/rotation.

### Group attachment

- **`attachGroupByMarker(groupParsed, markerX, markerY, opts)`**: attach a parsed group onto a backbone marker pair (iterative marker consumption).
- **`attachParsedByDirection(capAtom, groupParsed, params)`**: attach a parsed group using a cap/back bond direction + up vector + optional twist.

### Bond rebuild (geometry-derived)

- **`recalculateBonds(mmParams, opts)`**: brute-force O(N²) bond rebuild using covalent cutoffs.
- **`recalculateBondsBucketNeighbors(mmParams, bucketGraph, opts)`**: bucket-neighbor O(N)ish rebuild for large systems.
- **`recalculateBondsBucketAllPairsAABB(mmParams, bucketGraph, opts)`**: AABB-tested bucket all-pairs rebuild.

### Renderer boundary

- **`exportToMoleculeSystem(ms)`**: pack current atoms/bonds/selection into the renderer’s `MoleculeSystem` buffers.

### File I/O

- **`toXYZString(opts)`**: export to XYZ (optionally includes lattice vectors and charges).
- **`toMol2String(opts)`**: export to MOL2 (optionally includes lattice vectors).
- **`parseXYZ(text)`**: parse XYZ into `{pos, types, bonds, lvec}`.
- **`parseMol2(text)`**: parse MOL2 into `{pos, types, bonds, lvec}`.

### Static symbol/element utilities

- **`normalizeSymbol(s)`**: normalize element symbols (`cl`→`Cl`).
- **`symbolToZ(sym)`**: element symbol → atomic number.
- **`asZ(x)`**: accept number or symbol and return atomic number.
- **`SYMBOL_TO_Z` / `Z_TO_SYMBOL`**: small built-in element tables for common symbols.
