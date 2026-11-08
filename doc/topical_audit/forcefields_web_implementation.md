# Web-Based Force Field Implementations

FireCore includes browser-based implementations of force field kernels using WebGL and WebGPU. These are used for interactive molecular visualization, educational demos, and rapid prototyping without requiring a local GPU compute stack.

**Related Windsurf Codemaps:**
- No dedicated WebGL/WebGPU forcefield codemaps are currently listed in `topical_audit.md`. The WebGPU visualization codemap (`WebGPU Molecular Visualization & Physics Simulation`) would cover shader-based force evaluation if available.
- Cross-reference: [FireCore Force Field Navigation](https://windsurf.com/codemaps/d550a435-7c6f-47b1-aeb5-efc9b098564f-fe86ab10a43f3d18) traces the C++/OpenCL reference implementations that the WebGPU shaders mirror.

---

## 1. WebGL Implementation (`web/molgui_web/`)

### Architecture

The WebGL implementation uses standard HTML5 Canvas with WebGL 1.0/2.0 contexts. Force field evaluation is performed in JavaScript on the CPU or via fragment shaders for simple grid-based visualizations.

### Implementation Files

- **`web/molgui_web/js/ProjectiveDynamics.js`** — JavaScript port of the ProjectiveDynamics solver:
  - Implements the same Jacobi/Gauss-Seidel iteration logic as `ProjectiveDynamics_d.h`.
  - Uses dense matrices (Float32Array) and simple loops.
  - Suitable for small systems (<500 atoms) in educational demos.
  - Lacks sparse matrix support; limited to truss-like structures.

### Key Physics

The JS implementation mirrors the C++ math exactly:

$$\mathbf{A} \mathbf{x} = \mathbf{b}, \quad \mathbf{A} = \mathbf{K} + \frac{\mathbf{M}}{h^2} + \mathbf{D}$$

but with dense matrices stored as flat `Float32Array`. For a system with $N$ points and $M$ bonds, the stiffness matrix $\mathbf{K}$ is $3N \times 3N$ and assembled by:

```javascript
for (let b = 0; b < nBonds; b++) {
    let i = bonds[b].i, j = bonds[b].j;
    let k = bonds[b].k;
    // Add k to diagonal blocks K[ii], K[jj]
    // Add -k to off-diagonal blocks K[ij], K[ji]
}
```

### Performance Considerations

- **Single-threaded**: JavaScript runs on the main UI thread. For systems >200 atoms, the solver blocks the event loop.
- **No GPU acceleration**: WebGL 1.0 does not support compute shaders or atomic operations. The fragment shader pipeline is ill-suited to force field reduction.
- **Memory overhead**: Dense matrix storage limits systems to ~1000 atoms before browser memory limits are hit.

---

## 2. WebGPU Implementation (`web/molgui_webgpu/`)

### Architecture

WebGPU is the successor to WebGL and provides compute shaders, workgroups, and explicit memory management—making it viable for production-grade force field kernels in the browser.

### Implementation Files

- **`web/molgui_webgpu/`** — Main WebGPU molecular GUI:
  - `Draw3D_webgpu.js` — WebGPU rendering pipeline for atoms and bonds.
  - `ModularGL.py` / `ModularGL.js` — Modular shader generation system.
  - `CrystalUtils.js` — Lattice and PBC handling for periodic substrates.
- **`web/common_resources/shaders/`** — Shared shader code:
  - `compute_force.wgsl` — WGSL compute shader for non-bonded force evaluation.
  - `compute_grid_sample.wgsl` — WGSL shader for GridFF-style grid sampling.

### Key Physics

**Compute Shaders for Force Evaluation**:

WebGPU compute shaders use the same physics as the OpenCL kernels, translated to WGSL:

```wgsl
@compute @workgroup_size(64)
fn eval_nonbonded(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i: u32 = gid.x;
    if (i >= nAtoms) { return; }

    var f: vec3<f32> = vec3(0.0);
    for (var j: u32 = 0; j < nAtoms; j = j + 1) {
        if (i == j) { continue; }
        let dp: vec3<f32> = pos[j] - pos[i];
        let r2: f32 = dot(dp, dp);
        let ir2: f32 = 1.0 / r2;
        let u2: f32 = (REQ[i].x * REQ[j].x) * ir2;
        let u6: f32 = u2 * u2 * u2;
        let vdW: f32 = u6 * sqrt(REQ[i].y * REQ[j].y);
        let fmag: f32 = -12.0 * (u6 - 1.0) * vdW * ir2;
        f = f + dp * fmag;
    }
    force[i] = f;
}
```

This is directly analogous to `getLJQH()` in `relax_multi.cl:164`.

**Grid Sampling in WGSL**:

The GridFF interpolation logic is also ported to WGSL. The `sample3D` kernel computes B-spline weights and samples the grid texture:

```wgsl
fn basis(u: f32) -> vec4<f32> {
    let t: f32 = 1.0 - u;
    return vec4(
        t*t*t / 6.0,
        (3.0*u*u*(u - 2.0) + 4.0) / 6.0,
        (3.0*u*(1.0 + u - u*u) + 1.0) / 6.0,
        u*u*u / 6.0
    );
}
```

This matches the OpenCL `basis()` function in `GridFF.cl:66` exactly.

### Performance Considerations

- **Workgroup Size**: WebGPU limits workgroup size to 256 threads (vs. 1024 in OpenCL). The WGSL kernels are tuned to 64 threads for portability across devices.
- **Storage Buffers vs. Textures**: GridFF grids are stored as 3D textures (not buffers) for hardware-accelerated trilinear interpolation. This is ~2× faster than manual B-spline evaluation in shader code.
- **Memory Limits**: WebGPU enforces a 1 GB buffer limit per device. For large grids ($512^3$ floats = 512 MB), only one channel (Pauli, London, or Coulomb) can be resident at a time.
- **No Pinned Memory**: Unlike OpenCL, WebGPU cannot share buffers with the host CPU. Data must be explicitly copied via `mapAsync()`, adding ~1 ms latency per frame for small systems.

---

## 3. Common Resources (`web/common_js/`)

### Spatial Data Structures

- **`web/common_js/Buckets.js`** — Generic `Bucket` + `BucketGraph` classes for spatial partitioning. Provides `addAtomIndex(i, pos, radius)` with VdW-aware AABB bounds, `overlapAABB()`, `findOverlapNeighbors()`, `exportFlat()` for GPU-ready typed arrays, `recalcBounds(mol, radius)`, `pruneEmptyBuckets()`, ID/index conversion for topology robustness. Also exports standalone `aabbOverlap3D()` and `dist2ToAabb()`.
- **`web/common_js/BucketGrid3D.js`** — Uniform 3D grid partitioning for crystal unit cells. `buildCrystalCellBucketsFromMol()`, `crystalCellKey()`, `buildWireframeCellVerts()`, `buildWireframeAABBVerts()`. Uses `BucketGraph` from `Buckets.js`.
- **`web/common_js/Buckets_SoA.js`** — Structure-of-Arrays variant for SIMD-friendly layout.
- **`web/common_js/BucketAABBs.js`** — AABB construction and overlap tests for rigid-body packing.

### Collision Workgroups & Nonbonded Solver

- **`web/common_js/CollisionWorkgroups.js`** — Surface atom selection and spatial workgroup partitioning (uses `BucketGraph` from `Buckets.js` internally):
  - `surfaceAtomIndices(mol)` — selects H atoms + undercoordinated heavy atoms (<4 heavy-atom neighbors)
  - `buildCollisionWorkgroups({ pos, mol, radius, groupCap })` — farthest-pivot clustering into bounded-size groups; AABB bounds computed via radius-aware `Bucket.addAtomIndex`; exports flat typed arrays via `BucketGraph.exportFlat()`
  - `buildExclIcol_1_2_3(...)` — sorted 1-2/1-3 exclusion lists in icol space (matches `getNonBond_ex2` kernel pattern)
  - **CPU analog (atlas relax):** `pyBall.nanocrystal_pipeline.mask_bulk_nonbond` zeros `REQs` of 4-coordinated heavies (same atom set). Spatial AABB groups are for sparse GPU eval; L0/L1 is O(N²) with that mask. See [`ChemAtlas_MMFF_relax.md`](../Topics/FTIR_Nanocrystals/ChemAtlas_MMFF_relax.md).

- **`web/common_js/Nonbonded.js`** — CPU nonbonded force evaluation with workgroup acceleration (imports `aabbOverlap3D`, `dist2ToAabb` from `Buckets.js`):
  - `collisionPair(dp, R_min, E_min, R_cut, R_cut2)` — JS port of `getSR_x2_smooth` parabolic collision potential (harmonic repulsive + quadratic smoothing), designed for Projective Dynamics integration
  - `ljqh_pair(dp, REQ)` — full Lennard-Jones + Coulomb + H-bond pair force
  - `computeNonbondByGroups(...)` — workgroup-accelerated solver: AABB overlap for group discovery, sorted exclusion pointer advance, supports both full-LJ and collision modes
  - `computeNonbondBruteForceKernelStyle(...)` — O(N²) reference for parity checking
  - Verified parity: exact match (machine precision) with margin ≥ max VdW diameter

- **`scripts/debug_nanocrystal_nonbond_groups.mjs`** — CLI for parity testing and collision pair export

**Physics:** The collision potential splits the interaction into a harmonic region (`r < R_cut`, linear constraint for PD solver) and a smoothing region (`R_cut < r < R_cut2`, external force), matching the design in `ToDo_FastCollision_3.md`. Parameters derived from VdW: `R_min = Ri + Rj`, `E_min` from `EvdW` mixing.

### Topology Export Pipeline

- **`web/common_js/MolIO.js`** — Generic molecule I/O: `loadMMParamsFromDir()`, `loadMolFromMol2()`, `applyPositions()`, `stripHHBonds()` (H···H is steric, not covalent), `bondsForVisualization()`, `getCrystalBondsFromFiles()`
- **`web/common_js/exportFF.js`** — Topology NPZ builder: orchestrates MMFFL topology + collision workgroups + exclusions, writes enriched NPZ with group arrays via `LinearizedTopologyNpz.js`
- **`web/molgui_webgpu/LinearizedTopologyNpz.js`** — Packs typed arrays into numpy `.npz` format

### Visualization

- **`web/common_js/nanocrystalViewer.html`** — p5.js WEBGL viewer with workgroup coloring, AABB rendering, and collision pair visualization (intra-group green, inter-group orange)

---

## 4. Status and Recommendations

| Implementation | Language | GPU Compute | Max Atoms | Status | Recommendation |
|--------------|----------|-------------|-----------|--------|----------------|
| **WebGL / JS** | JavaScript | No | ~500 | Legacy | Deprecated; use WebGPU |
| **WebGPU / WGSL** | WGSL + JS | Yes (Compute) | ~50,000 | Active | Recommended for browser demos |
| **OpenCL (C++)** | C / OpenCL | Yes | Unlimited | Production | Primary for research |
| **PyOpenCL** | Python + OpenCL | Yes | Unlimited | Production | Debugging and prototyping |

### Porting Strategy

The WebGPU shaders are maintained as **direct translations** of the OpenCL kernels. When updating the physics in `relax_multi.cl` or `GridFF.cl`, the corresponding WGSL files should be updated in tandem. The math is kept identical to ensure parity.

Key parity checkpoints:
- `evalLJQH()` in `relax_multi.cl:164` ↔ `eval_nonbonded.wgsl`
- `basis()` in `GridFF.cl:66` ↔ `basis()` in `compute_grid_sample.wgsl`
- `fe3d_pbc()` in `GridFF.cl:108` ↔ `sample_grid_3d()` in WGSL

---

## See Also

- [Topical Audit Index](topical_audit.md) — priority ranking, dependency graph, missing topics
- [Forcefields Overview](forcefields_overview.md) — high-level taxonomy of all force field classes
- [Intramolecular Forcefields](intramolecular_forcefields.md) — UFF, MMFFsp3, ProjectiveDynamics, XPBD, RigidBody
- [Non-Bonding Forcefields](nonbonding_forcefields.md) — NBFF, exclusion schemes, FMM, PME
- [Surface Interactions](surface_interactions.md) — GridFF, FoldedAtomicFunctions, Ewald2D
- [GUI Feature Audit](gui_audit.md) — visualization & editor feature matrices, VisPy consolidation plan
- [Molecular Topology Editors](molecular_topology_editors.md) — WebGPU molecular editor (EditableMolecule.js)

---

*Last updated: 2026-06-23*
