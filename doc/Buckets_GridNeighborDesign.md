# Design notes for flat Buckets + Grid/AABB neighbor search in MolGUI/CrystalUtils

# Status checklist (current JS implementation)

- [x] `web/common_js/Buckets.js`: JS bucket graph for storing membership + neighbor buckets (array-based, not typed arrays)
- [x] `web/common_js/GridIndexer.js`: grid helpers for mapping positions to cells (if/when used)
- [x] `web/common_js/BucketAABBs.js`: AABB helpers (kept for compatibility; newer AABB utilities exist elsewhere)
- [x] Bucket-based bond rebuild modes in `web/molgui_web/js/EditableMolecule.js`
- [x] Bucket visualization in MolGUI (wireframe boxes/cells) + debug atom→bucket lines
- [x] Robustness for topology edits (atom deletions): bucket storage can temporarily switch from indices to stable atom IDs
- [x] Post-process atom dedup on an existing molecule: `CrystalUtils.dedupMolAtomsByTolA(mol,tolA,opts)`
- [ ] Dedup pipeline refactor (generate → bucketize → dedup → rebuild → bond)
- [ ] Add correctness checks for dedup (charge consistency, deterministic ordering)
- [ ] Perf notes / benchmarking (n=1k, 10k, 100k)

# Implemented notes and “don’t forget” items (important for future work)

## Index vs ID storage (topology robustness)

- Buckets normally store **atom indices** for speed.
- Any topology operation that changes atom ordering (e.g. `swap-with-last` deletion) invalidates stored indices.
- To keep buckets general/reusable (even when not rebuildable from metadata), we use a **temporary flip**:
  - `BucketGraph.toIds(mol)` before topology change
  - perform topology change (delete atoms)
  - `BucketGraph.toInds(mol)` after topology change (drops deleted atoms)

Justification:
- **Fast path stays fast** (indices).
- **Robustness is explicit** around topology edits.
- Avoids forcing a full rebuild of bucket membership from metadata (not always possible for general buckets).

## Bucket bounds and visualization correctness

- Bucket visualizations use per-bucket `pmin/pmax`.
- After topology changes (and also after membership edits), these bounds must be recomputed:
  - `BucketGraph.recalcBounds(mol)` recomputes AABBs from current member atom positions.

## Empty bucket pruning

- After deletions, some buckets become empty (all their member atoms got deleted).
- We prune them to avoid rendering a full empty grid and to keep pair loops tight:
  - `BucketGraph.pruneEmptyBuckets()` removes empty buckets and remaps `bucket.neigh[]`.
  - Note: pruning changes bucket indexing; any dense-cell → bucket index mapping must be treated as invalid.

## MolGUI debug UI / rendering hooks

- `autoUpdateBuckets` (default `true`): after topology edits (deleteSelection) refresh bucket debug overlays.
- `showBucketAtomLines` (default `false`): show line segments from each atom to its bucket center (computed from AABB bounds).
- All bucket debug refresh work is centralized in `window.app.refreshBucketDebug()`.

## Node test scripts in `web/molgui_web/scripts`

- `web/molgui_web/package.json` sets `"type":"module"`, so in that directory tree Node treats **`.js` as ES modules**.
- Therefore `test_*.js` and `test_*.mjs` are redundant there; keep only one style.
- Dedup test: `web/molgui_web/scripts/test_dedupMolAtomsByTolA.js` (or `.mjs`) can be run via `node`.

# Flat Buckets + Grid/AABB Neighbor Search (JS port)

## Motivation

We want to accelerate multiple operations which are essentially neighborhood queries:

- **De-duplication** (merge atoms closer than `tol`)
- **Bond construction** (find neighbors within type-dependent cutoff)
- (Later) other local queries (nonbonded lists, contact search, selection, etc.)

A naive implementation is often `O(n^2)` in number of atoms. For practical systems (typically **~10^3 atoms**, sometimes **up to 10^5**) this becomes slow.

The desired approach is a **flat bucketing** data structure:

- `n = k * m` where
  - `m` = number of buckets
  - `k` = typical atoms per bucket (goal ~ **8–16**, e.g. crystal unit cell, benzene ring)

If we can restrict comparisons to a small neighborhood of buckets, we can reduce work to approximately:

- **Grid case (best):** `O(n * k)` (27 neighbor buckets)
- **General bucket case:** `O(m^2)` bucket overlaps (AABB-pair) and inside overlap `O(k^2)`

The key requirement is to have a **single reusable bucket membership storage**, and a **pluggable neighbor enumeration**:

- **Grid neighbor enumeration** when buckets correspond to a spatial grid (fast)
- **AABB-pair enumeration** when buckets are arbitrary groups (still useful, avoids per-atom global search)

We explicitly avoid BVH/trees; only flat buckets.

## Design decisions

### 1) Keep core `Buckets` as “storage only”

The core `Buckets` class (ported from C++ `cpp/common/dataStructures/Buckets.h`) should:

- store bi-directional mapping between objects and buckets
- be compact and cache-friendly
- not include geometry, neighbor logic, or special cases

**Reasoning**:

- We want the same structure to serve many use cases (grid buckets, chemical groups, unit cells, fragments).
- Geometry and neighbor enumeration vary; storage does not.

### 2) Only two neighbor providers: Grid and AABB-pair

Neighbor enumeration should be implemented by *providers*:

- **Grid neighbor provider**: enumerates 27 neighboring buckets using integer index math
- **AABB-pair provider**: enumerates candidate bucket pairs using AABB overlap

**Grid is a special case of AABB** (each grid cell has an implicit AABB), but we do not need to materialize AABBs unless useful for visualization/debug.

### 3) Support both dense and sparse grid indexers

Real systems often occupy only a small fraction of a bounding box (vacuum, polymers on surface). A dense grid of `nx*ny*nz` can contain many empty cells.

We therefore support two *grid indexers* (both used with the same grid neighbor enumeration):

- **DenseGridIndexer**: allocates all cells, fastest for dense small grids
- **SparseGridIndexer**: allocates only occupied cells, avoids empty cells

Both indexers map `(ix,iy,iz)` to `icell`, but differ in storage:

- Dense: direct indexing `icell = ix + nx*(iy + ny*iz)`
- Sparse: `Map(key(ix,iy,iz) -> icell)` + arrays `icell -> (ix,iy,iz)`

### 4) JS implementation uses arrays for flexibility

The current JS implementation favors **plain JS arrays** (for bucket membership and neighbor lists) rather than a strict typed-array port.

Justification:
- Easier to support dynamic membership edits, prune buckets, and keep it “general-purpose reusable”.
- Typed arrays can be reintroduced later for specific hot paths once algorithms stabilize.

## Target problem sizes (web MolGUI)

- Typical: **n ~ 1,000 atoms**
- Possible: **n ~ 100,000 atoms**
- Space scale: up to ~ **400×400×400 Å**
- Cell size: typical **5–10 Å** → in worst box that is ~40³ = 64k logical grid coordinates

This motivates sparse grid option for vacuum-heavy systems.

# Core classes (pseudocode)

## `Buckets` (membership storage)

Port of `cpp/common/dataStructures/Buckets.h`.

### Data members

- `ncell`            : number of buckets (cells)
- `nobj`             : number of objects
- `nobjSize`         : allocated capacity for objects
- `cellNs[ncell]`    : number of objects in each cell
- `cellI0s[ncell]`   : prefix offsets into `cell2obj`
- `cell2obj[nobj]`   : contiguous object lists per cell
- `obj2cell[nobj]`   : mapping object -> cell index (`-1` means unassigned)

### Update pipeline

```text
updateCells(obj2cell):
  cellNs[:] = 0
  for i in 0..nobj-1:
    ic = obj2cell[i]
    if ic>=0: cellNs[ic]++

  # prefix sum
  ntot=0
  for ic in 0..ncell-1:
    cellI0s[ic] = ntot
    ntot += cellNs[ic]
    cellNs[ic] = 0          # reset for fill phase

  # fill
  for i in 0..nobj-1:
    ic = obj2cell[i]
    if ic<0: continue
    j = cellI0s[ic] + cellNs[ic]
    cell2obj[j] = i
    cellNs[ic]++
```

### Query

```text
range = inCellRange(ic):
  i0 = cellI0s[ic]
  n  = cellNs[ic]
  return [i0, n]
```

## `DenseGridIndexer`

### Data members

- `p0` (Vec3) origin
- `h`  (float) cell size
- `nx,ny,nz`

### Mapping

```text
ix = floor((x - p0.x)/h)
iy = floor((y - p0.y)/h)
iz = floor((z - p0.z)/h)
if outside -> -1
icell = ix + nx*(iy + ny*iz)
```

### Neighbor enumeration

```text
for dx,dy,dz in {-1,0,1}:
  jx=ix+dx; jy=iy+dy; jz=iz+dz
  if inside:
     jcell = jx + nx*(jy + ny*jz)
     yield jcell
```

## `SparseGridIndexer`

### Data members

- same `p0,h`
- `Map key -> icell`
- arrays `cellX[icell], cellY[icell], cellZ[icell]`

### Keying

We need a stable integer key for `(ix,iy,iz)`.

Options:

- string key: `"ix,iy,iz"` (simple but slower)
- packed 32-bit key: use a bias and bit packing (fast)

For robustness, start with string key, upgrade later.

### Build from objects

```text
map.clear(); ncell=0
for each object i:
  (ix,iy,iz) = floor((p[i]-p0)/h)
  key = pack(ix,iy,iz)
  if key not in map:
    map[key] = ncell
    cellX[ncell]=ix; cellY[ncell]=iy; cellZ[ncell]=iz
    ncell++
  obj2cell[i] = map[key]
```

### Neighbor enumeration

```text
(icell -> ix,iy,iz) from cellX/Y/Z
for dx,dy,dz in {-1,0,1}:
  keyN = pack(ix+dx,iy+dy,iz+dz)
  jcell = map.get(keyN)
  if jcell exists: yield jcell
```

## `AABBPairProvider`

For arbitrary buckets (chemical groups, rings, etc.). Requires AABBs per bucket:

- `pmin[icell]`, `pmax[icell]`

### Enumerate candidate bucket pairs

```text
for i in 0..ncell-1:
  for j in i+1..ncell-1:
    if overlapAABB(pmin[i],pmax[i], pmin[j],pmax[j], margin):
       yield (i,j)
```

Then algorithms nest `cell2obj` loops over buckets `(i,j)`.

---

# Algorithms using Buckets + neighbor providers

## Deduplication (grid neighbor)

Inputs:
- `pos[i]` positions
- `type[i]` element/type
- `q[i]` charge (optional)
- `tol` merge radius

```text
h = tol
build GridIndexer (dense or sparse)
obj2cell[i] = grid.cellOfPos(pos[i])
Buckets.updateCells(obj2cell)

keep[i] = i
for i in 0..n-1:
  ic = obj2cell[i]
  for each neighbor cell jc around ic:
    for each j in Buckets.objects(jc):
      if j<=i: continue
      if type[j]!=type[i]: continue
      if |q[j]-q[i]|>qtol: error
      if dist2(pos[i],pos[j]) < tol^2:
          merge j -> i (or to min index)

compact arrays using keep[] mapping
```

Notes:
- In practice use union-find or “keep smallest index” scheme.

## Bond construction (grid neighbor)

Inputs:
- final deduped atoms
- bond cutoff `rcut(type_i,type_j)` from mmParams

```text
h = maxBondCutoff
build grid buckets
for each atom i:
  for neighbor cells:
    for j in objects:
      if j<=i: continue
      if dist2 < rcut^2(type_i,type_j): addBond(i,j)
```

## General bucket interactions (AABB-pair)

Useful when buckets are groups not defined by spatial grid.

```text
compute bucket AABBs from member atom positions
for each overlapping pair (A,B):
  for i in bucket A:
    for j in bucket B:
      do expensive check
```

# Integration points in web code (MolGUI / CrystalUtils)

## Concrete integration points (current)

- `web/common_js/Buckets.js`
  - provides `Bucket` and `BucketGraph` used by bond rebuild and visualization
  - provides `buildCrystalCellBucketsFromMol()` used by crystal builders to create cell-buckets + neighbor lists
- `web/molgui_web/js/BuildersGUI.js`
  - stores last generated buckets in `window.app.lastBucketGraph`
  - triggers `window.app.updateBucketOverlay()` (and now `refreshBucketDebug()` can be used)
- `web/molgui_web/js/Editor.js`
  - wraps `deleteSelection()` with `toIds()` before deletion and `toInds()` after
  - triggers `refreshBucketDebug()` after deletion when `autoUpdateBuckets` is enabled
- `web/molgui_web/js/main.js`
  - owns the overlay objects and the `refreshBucketDebug()` method
- `web/molgui_web/js/GUI.js`
  - provides checkboxes toggling bucket debug flags

## Proposed pipeline (replace current “dedup during replication”)

1. **Generate atoms only** (replication, symmetry, plane cuts)
2. **Build buckets**
3. **Dedup once** (grid provider)
4. **Rebuild final arrays**
5. **Build bonds once** (grid provider)

This makes code paths simpler and avoids per-atom `Map` hashing inside replication loops.

# API summary (new / important functions and classes)

## `web/common_js/Buckets.js`

- `class Bucket`: per-bucket storage of member atoms + neighbor buckets + bounds (`pmin/pmax`).
- `class BucketGraph`: container of buckets + optional metadata; supports topology-robust conversion and helpers.
- `BucketGraph.toIds(mol)`: convert bucket membership from atom indices to stable atom IDs (before topology change).
- `BucketGraph.toInds(mol)`: convert bucket membership from atom IDs back to indices (after topology change); drops deleted atoms.
- `BucketGraph.recalcBounds(mol)`: recompute each bucket AABB (`pmin/pmax`) from current member atom positions.
- `BucketGraph.pruneEmptyBuckets()`: remove buckets with no atoms; remap `bucket.neigh[]`.
- `BucketGraph.getBucketCenterFromBounds(ib,out)`: compute bucket center as `(pmin+pmax)/2` (used for debug lines).
- `buildCrystalCellBucketsFromMol(mol, na,nb,nc, lvec, origin)`: build buckets by crystal unit cell index and precompute neighbor bucket lists.
- `buildWireframeCellVerts(A,B,C,O, out, i0)`: build wireframe line vertices for a parallelepiped cell (visualization).
- `buildWireframeAABBVerts(pmin,pmax, out, i0)`: build wireframe line vertices for an AABB (visualization).

## `web/molgui_web/js/Editor.js`

- `Editor.deleteSelection()`: wraps atom deletion with bucket `toIds()` / `toInds()` and triggers bucket debug refresh when enabled.
- `Editor.recalculateBonds()`: supports bucket-based modes; ensures buckets are in index mode before bucket bond search.

## `web/molgui_web/js/CrystalUtils.js`

- `dedupFracSitesByTolA(sites,lvec,tolA)`: deduplicate fractional unit-cell sites using a grid neighbor check.
- `dedupMolAtomsByTolA(mol,tolA,opts)`: post-process deduplicate atoms in an existing `EditableMolecule` (default `tolA=0.1`), with `bPrint` and `bError` options.

## `web/molgui_web/js/main.js`

- `MolGUIApp.refreshBucketDebug()`: recompute bucket indices, prune empties, recompute bounds, update bucket overlay + optional atom→bucket lines.
- `MolGUIApp.updateBucketOverlay()`: rebuild bucket wireframe visualization (cells or AABBs).
- `MolGUIApp.autoUpdateBuckets`: flag enabling automatic refresh after topology edits.
- `MolGUIApp.showBucketAtomLines`: flag enabling atom→bucket debug lines.

## `web/molgui_web/js/GUI.js`

- `Auto-update buckets` checkbox: toggles `window.app.autoUpdateBuckets`.
- `Show atom→bucket lines` checkbox: toggles `window.app.showBucketAtomLines` and refreshes debug overlay.



