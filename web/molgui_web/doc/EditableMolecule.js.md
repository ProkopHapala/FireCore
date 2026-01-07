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

## refactoring

Here’s a fresh high-level overview of the EditableMolecule ecosystem (one-line bullets, purpose-first, no args):

Core [EditableMolecule](cci:2://file:///home/prokop/git/FireCore/web/molgui_web/js/EditableMolecule.js:240:0-1188:1) (base chem/geom/topology)
- Atoms/bonds/fragments management, add/remove/update, selection set, transforms (translate/rotate), VSEPR helpers (missingDirsVSEPR/orthonormalBasisFromDir), bond rebuilds, lattice replicate, attach/remove by IDs, dirty flags.
- Parsing/IO, selection, and high-level attach/polymer are intentionally stubbed here—extensions install the working methods.

`MoleculeIO` (IO and symbol helpers)
- [exportAsParsed](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:4:0-22:1) / [appendParsedSystem](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:24:0-68:1): convert to/from flat arrays (pos/types/bonds/lvec) for interoperability.
- [toXYZString](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:70:0-93:1) / [toMol2String](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:95:0-132:1), [parseXYZ](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:180:0-202:1) / [parseMol2](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:134:0-178:1): serialize/parse common formats.
- Symbol utilities: [normalizeSymbol](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:204:0-210:1), [symbolToZ](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:212:0-218:1), [asZ](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:220:0-225:1), `SYMBOL_TO_Z`, `Z_TO_SYMBOL`.
- Installer: [installMoleculeIOMethods(EditableMolecule)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeIO.js:238:0-254:1) wires the above as instance + static helpers.

`MoleculeSelection` (selection DSL)
- [parseCountSet](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:2:0-16:1), [compileTokenSetToMatcher](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:18:0-49:1), [compileSelectQuerySpec](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:51:0-87:1), [applySelectQuery](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:89:0-134:1): token/neighbor-count DSL to build and apply selections.
- Installer: [installMoleculeSelectionMethods(EditableMolecule)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:136:0-164:1) adds instance methods and overrides static stubs so [compileSelectQuery](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/EditableMolecule.js:590:4-593:5)/[applySelectQuery](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeSelection.js:89:0-134:1) work after install; logs when invoked.

`MoleculeUtils` (high-level operations)
- Polymer/build: [assemblePolymerFromTokens](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:5:0-64:1).
- Attach/merge: [appendMolecule](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:110:0-115:1), [attachMoleculeByMarker](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:117:0-122:1), [attachMoleculeByDirection](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:124:0-129:1) (also parsed variants), replication: [replicateMolecule](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:168:0-214:1).
- Helper frames and marker finders kept local for attach workflows.
- Installer: [installMoleculeUtilsMethods(EditableMolecule)](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/MoleculeUtils.js:351:0-384:1) adds these as instance helpers.

Install order/rationale
- Core keeps stubs to avoid hard dependencies; extensions install real methods at startup.
- Call once after imports (as done in GUI.js):
  ```
  installMoleculeIOMethods(EditableMolecule);
  installMoleculeUtilsMethods(EditableMolecule);
  installMoleculeSelectionMethods(EditableMolecule);
  ```
- Installer logs indicate when wiring occurred; ensure installs run before any static [compileSelectQuery](cci:1://file:///home/prokop/git/FireCore/web/molgui_web/js/EditableMolecule.js:590:4-593:5) or parse calls.

Context: parsed vs molecule-level
- “Parsed” refers to raw array structures (pos/types/bonds/lvec) from format parsers; used for fast transforms/removals before converting to live molecules.
- Molecule-level attach/merge wraps parsed flow by exporting from source molecules, preserving API compatibility.

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

---

## Architectural notes (modular add-ons)

- **Dependency direction**: `MMParams` → `EditableMolecule` → optional add-ons (`MoleculeIO`, `MoleculeSelection`, `MoleculeUtils`). Core stays independent of add-ons to avoid cycles.
- **MoleculeIO**: imports `EditableMolecule` and exposes `installMoleculeIOMethods(cls=EditableMolecule)` to patch instance helpers (`exportAsParsed`, `appendParsedSystem`, `toXYZString`, `toMol2String`) and static symbol helpers (`normalizeSymbol`, `symbolToZ`, `asZ`). Pure IO functions remain exported. @web/molgui_web/js/MoleculeIO.js#1-249
- **MoleculeSelection**: imports `EditableMolecule` and exposes `installMoleculeSelectionMethods(cls=EditableMolecule)` to add instance methods that compile/apply queries using `mol.mmParams` and `cls.asZ`. DSL helpers stay exported. @web/molgui_web/js/MoleculeSelection.js#1-150
- **MoleculeUtils**: imports `EditableMolecule` and exposes `installMoleculeUtilsMethods(cls=EditableMolecule)` to attach `assemblePolymerFromTokens` as static and instance helper. @web/molgui_web/js/MoleculeUtils.js#1-74
- **EditableMolecule.mmParams**: constructor now has an optional `mmParams` property so add-ons can read it without threading args. @web/molgui_web/js/EditableMolecule.js#241-256

---

Here’s a one-line purpose + context for every function/method in `EditableMolecule.js`:

- **_bondCut2** – Computes squared bond cutoff from covalent radii (or default) to decide if two atoms should be bonded; used by all bond-rebuild routines. @web/molgui_web/js/EditableMolecule.js#4-21  
- **_maxRcovFromMolAtoms** – Finds the largest covalent radius in the molecule to size search margins for spatial bucket bond builds. @web/molgui_web/js/EditableMolecule.js#23-33  
- **Atom** – Data holder for a single atom (id, type, charge, bonds, position) used throughout EditableMolecule state. @web/molgui_web/js/EditableMolecule.js#35-48  
- **_ensureMMTypeIndex** – Builds cached name↔index maps for mmParams atom types so later lookups are O(1). @web/molgui_web/js/EditableMolecule.js#50-60  
- **_atomTypeNameToIndex** – Resolves atom-type name to index via mmParams caches; used by type parsing and selection. @web/molgui_web/js/EditableMolecule.js#62-66  
- **_atomTypeIndexToName** – Inverse lookup from index to atom-type name; used when back-computing types for atoms. @web/molgui_web/js/EditableMolecule.js#68-72  
- **_resolveTypeOrElementToAtomType** – Accepts an atom-type name or element symbol and returns standardized atom-type info; used when setting atom types or building caps. @web/molgui_web/js/EditableMolecule.js#74-98  
- **_parseCountSet** – Parses count-set strings like `{1,2}` for select queries’ neighbor-count constraints. @web/molgui_web/js/EditableMolecule.js#100-113  
- **_compileTokenSetToMatcher** – Turns token sets (elements or atom-types) into matcher functions for selection queries. @web/molgui_web/js/EditableMolecule.js#115-145  
- **_compileSelectQuery** – Compiles a selection query string into matchers and neighbor-count constraints for later application. @web/molgui_web/js/EditableMolecule.js#147-182  
- **_getAtomTypeForAtom** – Resolves the effective atom-type object for an atom (by atype or element), feeding VSEPR and capping logic. @web/molgui_web/js/EditableMolecule.js#184-197  
- **_bondLengthEst** – Estimates bond length from covalent radii (or fallback); drives placement distances for caps/e-pairs. @web/molgui_web/js/EditableMolecule.js#199-208  
- **_orthonormalBasisFromDir** – Builds an orthonormal pair perpendicular to a direction; helper for VSEPR geometry. @web/molgui_web/js/EditableMolecule.js#210-218  
- **_missingDirsVSEPR** – Computes ideal missing bond/lone-pair directions from existing vectors using VSEPR; used to place caps and explicit electron pairs. @web/molgui_web/js/EditableMolecule.js#220-326  
- **Bond** – Bond record (endpoints, order, cached indices) with helpers to resolve atom indices against the molecule. @web/molgui_web/js/EditableMolecule.js#328-354  
- **Bond.other** – Returns the opposite atom index on a bond; utility for neighbor traversals. @web/molgui_web/js/EditableMolecule.js#343-344  
- **Bond.ensureIndices** – Refreshes cached atom indices from ids when topology changes; used before most bond operations. @web/molgui_web/js/EditableMolecule.js#345-353  
- **`Bounds.reset/addPoint/finalize/intersectsSphere/intersectsAABB`** – Maintains bounding boxes/spheres for fragments to speed overlap checks. @web/molgui_web/js/EditableMolecule.js#356-406  
- **Fragment (ctor)** – Tracks a fragment’s atom/bond ids and bounds; used in fragment decomposition. @web/molgui_web/js/EditableMolecule.js#408-416  
- **Fragment.updateBounds** – Recomputes fragment bounds from current atom positions. @web/molgui_web/js/EditableMolecule.js#417-425  
- **EditableMolecule (ctor)** – Initializes molecule state (atoms/bonds/fragments, selection, lattice, dirty flags). @web/molgui_web/js/EditableMolecule.js#428-456  
- **EditableMolecule.nAtoms** – Convenience getter for atom count. @web/molgui_web/js/EditableMolecule.js#458  
- **EditableMolecule._touchTopo** – Marks topology-dependent caches dirty and bumps topoVersion after edits. @web/molgui_web/js/EditableMolecule.js#460-465  
- **EditableMolecule._touchGeom** – Marks geometry/export dirty after coordinate changes. @web/molgui_web/js/EditableMolecule.js#467-471  
- **EditableMolecule._assertUnlocked** – Guards topology mutations unless lockDepth is zero. @web/molgui_web/js/EditableMolecule.js#473-475  
- **EditableMolecule.lockTopology/unlockTopology/assertLocked** – Manage a topo lock and detect accidental edits during locked operations. @web/molgui_web/js/EditableMolecule.js#477-492  
- **EditableMolecule.getAtomIndex/getBondIndex** – Map ids to current array indices. @web/molgui_web/js/EditableMolecule.js#494-502  
- **EditableMolecule.addAtom** – Adds an atom (legacy signature) and marks topology dirty; base for most builders/importers. @web/molgui_web/js/EditableMolecule.js#504-513  
- **EditableMolecule.setAtomTypeByName** – Applies an atom type to an atom (and optionally adjusts Z) using mmParams; used before VSEPR/capping. @web/molgui_web/js/EditableMolecule.js#515-524  
- **EditableMolecule.addCappingAtoms** – Adds hydrogens or other caps to undercoordinated atoms using VSEPR directions and optional bonding; uses selection by default. @web/molgui_web/js/EditableMolecule.js#526-602  
- **EditableMolecule.addExplicitEPairs** – Inserts explicit lone-pair pseudo-atoms in VSEPR directions where needed. @web/molgui_web/js/EditableMolecule.js#604-671  
- **addAtomZ** – Convenience to add an atom by Z (aliases addAtom). @web/molgui_web/js/EditableMolecule.js#673  
- **setAtomPosById** – Updates an atom’s coordinates and marks geometry dirty. @web/molgui_web/js/EditableMolecule.js#675-681  
- **addBond** – Creates a bond between two atom ids, wiring adjacency lists and caching indices. @web/molgui_web/js/EditableMolecule.js#683-700  
- **`selectAtom/select`** – Mutate selection by id or index-compatible input; drives downstream ops like delete or capping. @web/molgui_web/js/EditableMolecule.js#702-718  
- **translateAtoms** – Translates a set of atoms by a Vec3 and marks geometry dirty. @web/molgui_web/js/EditableMolecule.js#720-729  
- **rotateAtoms** – Rotates atoms around an axis/center using Mat3; used for interactive transforms. @web/molgui_web/js/EditableMolecule.js#731-758  
- **`clearSelection/selectAll`** – Reset or fill selection set, marking export dirty. @web/molgui_web/js/EditableMolecule.js#761-775  
- **EditableMolecule.compileSelectQuery** – Static wrapper to compile selection queries with mmParams. @web/molgui_web/js/EditableMolecule.js#777-780  
- **EditableMolecule.applySelectQuery** – Executes a compiled selection query against the molecule, updating selection set. @web/molgui_web/js/EditableMolecule.js#782-826  
- **EditableMolecule.deleteSelectedAtoms** – Removes all atoms in selection (and their bonds). @web/molgui_web/js/EditableMolecule.js#828-835  
- **`_removeBondIndexFromAtom/_replaceBondIndexInAtom`** – Internal helpers to maintain per-atom bond index lists during removals/moves. @web/molgui_web/js/EditableMolecule.js#837-876  
- **`removeBondByIndex/removeBondById`** – Delete bonds safely, compacting bond array and adjacency. @web/molgui_web/js/EditableMolecule.js#844-883  
- **`removeAtomByIndex/removeAtomById`** – Delete atoms (and incident bonds), compacting arrays and selection. @web/molgui_web/js/EditableMolecule.js#884-910  
- **clear** – Wipes all atoms/bonds/fragments/selection, resetting topo state. @web/molgui_web/js/EditableMolecule.js#912-921  
- **addAtomsFromArrays** – Bulk-add atoms from flat position and type arrays; used by importers. @web/molgui_web/js/EditableMolecule.js#923-933  
- **updateNeighborList** – No-op placeholder (adjacency is maintained incrementally). @web/molgui_web/js/EditableMolecule.js#935  
- **attachGroupByMarker** – Attaches a parsed group to backbone markers (X/Y) via computed rotation/translation, deleting markers after bonding; used for modular assembly. @web/molgui_web/js/EditableMolecule.js#937-979  
- **attachParsedByDirection** – Attaches a parsed fragment to a “cap” atom following specified forward/up references and optional twist; used for directional grafting. @web/molgui_web/js/EditableMolecule.js#981-1041  
- **EditableMolecule._getParsedPos** – Utility to read a Vec3 from flat parsed positions. @web/molgui_web/js/EditableMolecule.js#1043-1046  
- **EditableMolecule._findMarkerPairsMol** – Finds backbone marker pairs (X–Y plus anchor) in the current molecule for group attachment. @web/molgui_web/js/EditableMolecule.js#1048-1065  
- **EditableMolecule._findMarkerPairsParsed** – Same as above but on parsed group data. @web/molgui_web/js/EditableMolecule.js#1068-1084  
- **EditableMolecule._buildFrame** – Builds an orthonormal frame from forward/up vectors; foundational for attachment transforms. @web/molgui_web/js/EditableMolecule.js#1087-1095  
- **EditableMolecule._rotateFrameAroundForward** – Rotates a frame around its forward axis by an angle; used for twist in attachments. @web/molgui_web/js/EditableMolecule.js#1098-1106  
- **EditableMolecule._computeMarkerAttachRotation** – Computes rotation aligning group markers to backbone markers for attachment. @web/molgui_web/js/EditableMolecule.js#1109-1132  
- **EditableMolecule._transformParsed** – Applies rotation/translation to parsed group coordinates relative to an anchor. @web/molgui_web/js/EditableMolecule.js#1135-1149  
- **EditableMolecule._removeAtomsFromParsed** – Produces a copy of parsed data with specific atoms removed and indices remapped; used when discarding markers. @web/molgui_web/js/EditableMolecule.js#1152-1177  
- **EditableMolecule.replicate** – Replicates the current system over an n×m×k lattice (updating lvec) to build supercells. @web/molgui_web/js/EditableMolecule.js#1180-1225  
- **EditableMolecule.exportAsParsed** – Exports molecule to flat arrays (pos, types, bonds, lvec) for interoperability. @web/molgui_web/js/EditableMolecule.js#1227-1243  
- **EditableMolecule.appendParsedSystem** – Imports a parsed system into this molecule with optional rotation/translation, returning new atom ids; used by attachments and replication. @web/molgui_web/js/EditableMolecule.js#1246-1288  
- **EditableMolecule.recalculateBonds** – Rebuilds all bonds by brute-force pairwise distance with covalent radii cutoffs. @web/molgui_web/js/EditableMolecule.js#1291-1323  
- **EditableMolecule.recalculateBondsBucketNeighbors** – Faster bond rebuild using neighbor buckets (precomputed neighbors) to limit pair checks. @web/molgui_web/js/EditableMolecule.js#1325-1372  
- **recalculateBondsBucketAllPairsAABB** – Bond rebuild using bucket AABB overlap pruning plus covalent-radius margins. @web/molgui_web/js/EditableMolecule.js#1374-1424  
- **exportToMoleculeSystem** – Copies current editable molecule into a `MoleculeSystem` structure (positions, types, bonds, selection). @web/molgui_web/js/EditableMolecule.js#1426-1463  
- **toXYZString** – Serializes the molecule to XYZ format (with optional charges and lattice vectors). @web/molgui_web/js/EditableMolecule.js#1465-1487  
- **toMol2String** – Serializes to MOL2 format including bonds and optional lattice vectors. @web/molgui_web/js/EditableMolecule.js#1489-1525  
- **normalizeSymbol** – Normalizes element symbols to proper capitalization. @web/molgui_web/js/EditableMolecule.js#1527-1532  
- **symbolToZ** – Converts element symbol to atomic number via lookup. @web/molgui_web/js/EditableMolecule.js#1534-1539  
- **asZ** – Coerces number or symbol to atomic number; used throughout parsing. @web/molgui_web/js/EditableMolecule.js#1541-1545  
- **parseMol2** – Parses MOL2 text into flat arrays (pos/types/bonds/lvec) for import. @web/molgui_web/js/EditableMolecule.js#1547-1594  
- **parseXYZ** – Parses XYZ text into positions/types arrays. @web/molgui_web/js/EditableMolecule.js#1597-1617  
- **`EditableMolecule.SYMBOL_TO_Z` / `Z_TO_SYMBOL`** – Minimal symbol↔Z lookup tables used by parsing/serialization. @web/molgui_web/js/EditableMolecule.js#1620-1629