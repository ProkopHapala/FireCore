# Linearized Topology (MMFFL) — Design & Progress

**Scope:** Build a **spring-only** force field topology (1–2, 1–3, 1–4 interactions as harmonic sticks) in **JavaScript (Node headless)**, with visual debug of stick classes **K₁₂ / K₁₃ / K₁₄**, for later use in vibration Hessian assembly and GPU Jacobi solvers.

**Related:** `web/molgui_webgpu/MMFFLTopology.js`, `web/molgui_web/js/ProjectiveDynamics.js`, `pyBall/OCL/MMFFL.py`, `pyBall/OCL/LFFSolver.py`, `intramolecular_forcefields.md`.

**Status:** **Implemented** (Job 2, 2026-06-15). Headless batch + K₁₄ + NPZ export + HTML viewers. See daily report below.

**Fixtures:** `Parallel_agent_fixtures.guide.md` — bootstrap before starting.

---

## Design goals

1. **JavaScript-first batch path:** `node scripts/build_linearized_topology.mjs` — faster iteration than Python MMFFL for large catalogs of generated NCs.
2. **Complete bonded linearization for diamond/sp³:**
   - **K₁₂ (red):** primary σ bonds (1st neighbors).
   - **K₁₃ (green):** angle-derived springs between **2nd neighbors** (graph distance 2, “meta” on hexagon).
   - **K₁₄ (blue):** dihedral/torsion-derived springs between **4th neighbors** (graph distance 4, “para” on hexagon — opposite site of the six-ring).
3. **Explicit stiffness channels:** Store `(i, j, l0, K, tag)` with `tag ∈ {bond, angle, dihedral}`; allow global scalars `K12`, `K13`, `K14` overrides for debugging.
4. **Bounded neighbor lists:** Know **max sticks per atom** for GPU packing (`MAX_NEIGHBORS` in kernels).
5. **Visual debug with minimal code:** Auto-generate a standalone HTML viewer (p5.js WEBGL or Three.js CDN) — mouse orbit, color-coded sticks, no full molgui dependency.
6. **Parity with Python:** Same topology on adamantane / minimal diamond as `MMFFL.build_topology()` within numerical tolerance.

---

## Physics: what each stick class means

### K₁₂ — 1–2 bonds (existing)

Harmonic: \(E = \frac{1}{2} k (r - r_0)^2\) between bonded pairs.

**Source:** `BondTypes.dat` via `mmParams.getBondParams(ta, tb)`.

**Implementation:** `buildMMFFLTopology` primary pass (`bondsAdj1`) — **done**.

### K₁₃ — 1–3 angles as 2nd-neighbor springs (existing)

For angle **j–i–k**, add spring between **j** and **k** with:

\[
l_0 = \sqrt{r_{ij}^2 + r_{ik}^2 - 2 r_{ij} r_{ik} \cos\theta_0}
\]

**Source:** `AtomTypes.dat` → `Ass`, `Kss` on central type at **i**.

**Implementation:** `buildAngleBonds()` in `MMFFLTopology.js` — **done** (tag `'angle'`).

**Graph distance:** endpoints **j**, **k** are 2nd neighbors of each other through **i** (distance 2 in bond graph).

### K₁₄ — 1–4 torsion as distance spring (**implemented**)

For dihedral **a–b–c–d**, add spring between **a** and **d** (endpoints of 1–4 interaction) with equilibrium distance \(l_0\) derived from:

- Bond lengths \(r_{ab}, r_{bc}, r_{cd}\)
- Equilibrium angles \(\theta_{abc}, \theta_{bcd}\) from types at **b**, **c**
- Equilibrium dihedral \(\phi_0\) from `DihedralTypes.dat` (or UFF fallback)

**Implementation:** `buildDihedralBonds()` + `enumerateDihedralQuadruples()` in `MMFFLTopology.js`; `MMParams.dihedralEndDistanceL0()` / `convertDihedralToDistance()` in `MMParams.js`.

1. Enumerate dihedral quadruples **a–b–c–d** (central bond **b–c**), dedupe undirected **(a,d)** pairs.
2. \(l_0(a,d)\) from equilibrium geometry: bond lengths \(r_{ab}, r_{bc}, r_{cd}\), angles \(\theta_{abc}, \theta_{bcd}\) (`Ass` / `AngleTypes.dat`), dihedral \(\phi_0\) (`DihedralTypes.dat` or **UFF sp³–sp³** fallback, \(\phi_0=60°\)).
3. \(k_r\) from numerical \(dl/d\phi\) at \(\phi_0\): \(k_r = k_\phi / (dl/d\phi)^2\).

**Graph distance note:** MMFF “1–4” dihedral endpoints **a**, **d** are **3 bonds apart** (path a–b–c–d). This is distinct from hexagon “para” pairs at bond-graph distance **4** (histogram bucket still useful for crystal inventory).

**Diamond / hexagon analogy (user intent):**

On a {111} hexagonal ring of Si atoms:

| Label | Graph dist along ring | Role |
|-------|----------------------|------|
| ortho | 1 (1–2) | bond |
| meta | 2 (1–3) | angle spring |
| para | 3 on ring = 4 in bond graph along chair path | 1–4 spring |

In bulk diamond, graph-distance-4 pairs are the “opposite” sites of local six-membered motifs — exactly the partners needed to replace torsion stiffness without explicit dihedral forces.

---

## Maximum neighbors per atom (packing budget)

For **interior** sp³ atom in diamond lattice:

| Class | Count touching atom **i** | Notes |
|-------|---------------------------|-------|
| K₁₂ bonds | **4** | Tetrahedral |
| K₁₃ angle springs | **0 centered on i** (springs between pairs of i's neighbors); as **endpoint**: up to **6** (each neighbor pair (j,k) gives spring j–k; i is neighbor of both j and k for 3 pairs × 2?) — exact count: each of 4 neighbors participates in C(3,2)=3 angle centers → i endpoint of angle springs to 2nd neighbors: up to **4×3/2?** |

**Careful enumeration needed.** For interior Si:

- Springs **centered** on i as angle hub: C(4,2) = **6** virtual bonds between neighbor pairs (none use i as endpoint).
- Springs where **i** is endpoint of K₁₃: for each neighbor j, j has 3 other neighbors; i–k springs where k is 2nd neighbor of i — typically **4 × 3 = 12** directed, **~6** unique undirected 2nd-neighbor contacts per atom.

| K₁₄ dihedral springs | Estimate **12–24** directed entries; **~6–12** unique 4th-neighbor contacts per interior atom (to be counted exactly on G1/G2 crystals) |

**Combined per-atom neighbor list** (for Jacobi gather — all springs incident on i):

| Component | Conservative upper bound (interior diamond) |
|-----------|---------------------------------------------|
| K₁₂ | 4 |
| K₁₃ | 12 |
| K₁₄ | 12 |
| **Total** | **~28** |

**Current code limits:**

| Location | `MAX_NEIGHBORS` / `nMaxBonded` |
|----------|-------------------------------|
| `LFF.cl` | **8** — insufficient for full K₁₂+K₁₃+K₁₄ |
| `MMFFLTopology.buildXPDBInputsFromMol` | default **16** — OK for K₁₂+K₁₃ only |
| `XPDB_new` | **4** bonded ports — different paradigm |

**Design decision:** New vibration/LFF kernel should use **`MAX_NEIGHBORS = 32`** (power of 2) or separate **three lists** (bondIdx₁₂, bondIdx₁₃, bondIdx₁₄) with fixed small max each (4+6+8) to save registers.

---

## JavaScript headless batch architecture

```
gen_nanocrystals.mjs  →  minimal.xyz / .mol2
        ↓
build_linearized_topology.mjs   (NEW)
  - load MMParams (ElementTypes, AtomTypes, BondTypes, AngleTypes, DihedralTypes)
  - EditableMolecule.fromMol2 / fromXYZ
  - buildMMFFLTopology(mol, mm, { add_angle, add_dihedral, K12, K13, K14 })
  - write:
      topology.json     # [{i,j,l0,K,tag}, ...]
      topology.mol2       # optional: dummy sticks as pseudo-bonds for Jmol
      debug_viewer.html   # auto-generated p5/three viewer
```

**Reuse (do not rewrite):**

- `buildMMFFLTopology`, `packBondArrays`, `lawOfCosines`
- `MMParams.resolveTypeNameTable` (Si → `Si` / `Si3`)
- Type assignment: mirror `MMFFL.assign_type_names` with `type_source: 'table'`

**New code:**

- `buildDihedralBonds(mol, mmParams, bondsAdj1, opts, outLinear)` — K₁₄
- `enumerateGraphDistances(mol, maxDist)` — validate 1–2 / 1–3 / 1–4 classification
- `exportTopologyMol2(topo)` — color by tag or separate bond types
- `exportP5Viewer(topo, htmlPath)` — ~100–150 lines standalone HTML

---

## Visual debug viewer (low-effort option)

**Recommended:** Single static HTML + **Three.js** (CDN) or **p5.js WEBGL**:

- Atoms: small spheres (Si gray, H white).
- Sticks:
  - K₁₂ → **red**
  - K₁₃ → **green**
  - K₁₄ → **blue**
- OrbitControls (Three) or p5 `mouseDragged` rotation.
- Input: embed `topology.json` or pass via `fetch` from same directory.

**Alternative:** Export multi-segment `.xyz` with fake atoms at stick midpoints (hacky) — avoid.

**Alternative:** Use existing Jmol with `topology.mol2` if bond orders/types encode stick class — workable but less clear than RGB sticks.

**Deliverable:** Opening `debug_viewer.html` in browser shows adamantane-sized crystal with correct stick counts without starting molgui.

---

## Python parity (secondary)

Keep `pyBall/OCL/MMFFL.py` as reference implementation:

- After JS milestone M-L2, run `test_TiledJacobi_molecules.py --dump_topo_debug` vs Node `topology.json` diff.
- Python remains for XPDB MD; **batch NC catalog uses Node**.

---

## Planned work packages

### WP-L1 — Inventory & graph-distance validator

- On G1/G2 geometries from generation milestone: count graph distances 1,2,3,4 from each atom.
- Print histogram; confirm hexagon para = distance 4.

### WP-L2 — Extend `buildMMFFLTopology` with `add_dihedral` / K₁₄

- Load `DihedralTypes.dat` in `MMParams.js` if not already exposed.
- Implement `buildDihedralBonds` (law-of-geometry for \(l_0\)).
- Tags: `['dihedral', [a,b,c,d]]`.

### WP-L3 — Headless batch script

- `scripts/build_linearized_topology.mjs`
- CLI: `--input`, `--out-dir`, `--K12`, `--K13`, `--K14`, `--add-dihedral 0|1`, `--html-viewer 0|1`

### WP-L4 — HTML stick viewer generator

- Function writes self-contained HTML; document in this file.

### WP-L5 — Packing format for GPU solver

- Export `neighs[natoms, MAX_N]`, `KLs[natoms, MAX_N, 2]`, `stick_class[natoms, MAX_N]` (0=empty, 1/2/3).
- Assert no overflow; fail with atom id if `n_sticks > MAX_N`.

### WP-L6 — Adamantane parity

- Compare spring count and \((l_0, K)\) lists JS vs Python MMFFL.

---

## Testing milestones

| Milestone | Criterion |
|-----------|-----------|
| **M-L0** | Adamantane: K₁₂ count = number of σ bonds; K₁₃ count = number of angle triples (each undirected 2nd-neighbor pair once) |
| **M-L1** | G1 diamond primitive: every interior Si has 4 red sticks; K₁₃ springs only between 2nd neighbors (verify graph dist = 2) |
| **M-L2** | K₁₄ springs only between graph-distance-4 pairs; count matches manual enumeration on hexagon ring |
| **M-L3** | `debug_viewer.html` rotates; colors match legend; no WebGPU/molgui required |
| **M-L4** | `topology.json` round-trips: Node builder → Python reader → same adjacency |
| **M-L5** | `MAX_NEIGHBORS` budget documented and satisfied for G2 capped SiNC (largest atom degree &lt; limit) |
| **M-L6** | Optional: analytical spring Hessian for K₁₂+K₁₃+K₁₄ matches finite-difference on LFF energy (small molecule) |

---

## Relationship to vibration solver

Linearized topology defines **sparse K** (Hessian of harmonic network):

- Off-diagonal 3×3 blocks only between stick endpoints.
- Assembly is **O(n_sticks)** — trivial on CPU; no C++ MMFF FD needed for LFF-based vibrations.
- Full MMFF angle/dihedral FD remains reference for **parity** until LFF spectrum is trusted.

---

## Open questions

1. **Dihedral \(l_0\) formula:** Closed form from \((r_{ab}, r_{bc}, r_{cd}, \theta_0, \phi_0)\) vs tabulated from ideal diamond geometry?
2. **Stiffness mapping \(k_\phi \to k_r\):** Use finite-difference calibration on adamantane single dihedral?
3. **Merge vs split neighbor arrays:** One `MAX_NEIGHBORS=32` list vs three specialized arrays for GPU register pressure?
4. **Hydrogen caps:** Include H in K₁₃/K₁₄ or treat as massless caps (topology only on heavy-atom subgraph)?

---

## Parallel agent contract (Job 2)

Full rules: **`Parallel_agent_fixtures.guide.md`**.

| | |
|--|--|
| **Fixture root** | `tests/tSiNCs/fixtures/vibration_parallel/` (gitignored) |
| **Read** | `structures/adamantane.mol2`, `structures/si_G1_caps_only.mol2`, `structures/si_G2_facet111_caps_only.mol2` |
| **May edit** | `MMFFLTopology.js`, `MMParams.js`, `LinearizedTopologyNpz.js`, `LinearizedTopologyViewer.js`, `web/common_js/npzWrite.js`, `scripts/build_linearized_topology.mjs` |
| **May write** | `tests/tSiNCs/linearized/` (`.npz`, optional `.json`, debug HTML) |
| **Must NOT edit** | `gen_nanocrystals.mjs`, `EditableMolecule.js`, `FTIR.py`, `LFF.cl`, `LFFSolver.py`, C++ MMFF |
| **Parity (read-only)** | `pyBall/OCL/MMFFL.py` |

---

*Last updated: 2026-06-15 — Job 2 complete; module split + direct NPZ from Node.*

---

## Daily report — 2026-06-15 (Job 2 complete)

### Summary

Implemented the full **JavaScript headless** linearized-topology pipeline: K₁₂ + K₁₃ + K₁₄ spring lists, GPU packing, **`.npz`** export for `LFFSolver` / numpy, and two **Three.js** debug viewers (embedded JSON + NPZ loader). Validated on fixture structures adamantane, `si_G1_caps_only`, `si_G2_facet111_caps_only`.

**Refactor (same day):** topology math, NPZ I/O, and HTML viewers split into separate modules; `.npz` written **directly from Node** (typed arrays → binary, no JSON intermediate, no Python subprocess). JSON export is opt-in (`--export-json 1`) for debugging only.

### How it works (end-to-end)

```
.mol2 / .xyz
    → scripts/build_linearized_topology.mjs
        → MMParams (Element/Atom/Bond/Angle/Dihedral .dat)
        → EditableMolecule + bonds
        → buildMMFFLTopology({ add_angle, add_dihedral, K12, K13, K14 })   [MMFFLTopology.js]
        → packLinearTopologyForGPU({ maxNeighbors: 48 })                     [MMFFLTopology.js]
        → writeTopologyNpzFile(…)                                          [LinearizedTopologyNpz.js → npzWrite.js]
        → exportStickViewer* (optional HTML)                               [LinearizedTopologyViewer.js]
        → outputs under tests/tSiNCs/linearized/ (local, not for git)
```

| Output | Purpose |
|--------|---------|
| `{base}_topology.npz` | **Primary** — LFFSolver-ready arrays (`neigh_idx`, `KLs`, `pos`, `Z`, …); written directly from Node |
| `{base}_debug_viewer.html` | Self-contained viewer; topology **embedded** in HTML |
| `{base}_npz_viewer.html` | Same viewer UI; loads `.npz` via file picker or HTTP `fetch` |
| `{base}_topology.json` | **Optional** (`--export-json 1`) — meta + spring list only (no packed arrays) |

**CLI:**

```bash
node scripts/build_linearized_topology.mjs \
  --input tests/tSiNCs/fixtures/vibration_parallel/structures/adamantane.mol2 \
  --out-dir tests/tSiNCs/linearized \
  [--add-angle 1] [--add-dihedral 1] \
  [--K12 0] [--K13 0] [--K14 0] \
  [--max-neighbors 48] \
  [--html-viewer 1] [--export-npz 1] [--export-json 0] [--npz-viewer 1] \
  [--graph-dist 1]
```

### Module layout (post-refactor)

| Module | Role |
|--------|------|
| `web/molgui_webgpu/MMFFLTopology.js` | **Topology only** — build, pack, spring list, graph distances |
| `web/molgui_webgpu/LinearizedTopologyNpz.js` | Pack typed arrays + `writeTopologyNpzFile()` |
| `web/common_js/npzWrite.js` | Minimal numpy `.npy`/`.npz` encoder (Node zlib, no deps) |
| `web/molgui_webgpu/LinearizedTopologyViewer.js` | Three.js HTML generators (`exportStickViewer*`) |
| `scripts/build_linearized_topology.mjs` | Headless CLI orchestration |

**Removed:** `scripts/topology_to_npz.py` (obsolete — was JSON bundle → Python → `.npz`).

### Core library changes

**`web/molgui_webgpu/MMFFLTopology.js`** (topology math only)

| Export | Role |
|--------|------|
| `buildMMFFLTopology` | Extended: `add_dihedral`, `K12`/`K13`/`K14` overrides, `bonds_dihedral` |
| `buildDihedralBonds` | K₁₄ springs (internal) |
| `enumerateGraphDistances` | BFS histogram + per-atom max distance |
| `packLinearTopologyForGPU` | `neighs`, `KLs`, `stick_class`; fails loud if degree > `maxNeighbors` |
| `exportTopologySpringList` | `[{i,j,l0,K,tag}, …]` |

**`web/molgui_webgpu/LinearizedTopologyViewer.js`** (debug viewers — not topology)

| Export | Role |
|--------|------|
| `exportStickViewerHTML` | Embedded-data Three.js viewer |
| `exportStickViewerNpzLoaderHTML` | NPZ file-picker / HTTP viewer |

**`web/molgui_webgpu/LinearizedTopologyNpz.js`** + **`web/common_js/npzWrite.js`**

| Export | Role |
|--------|------|
| `buildTopologyNpzArrays` | Typed arrays from `topo` + `packing` + `mol` |
| `writeTopologyNpzFile` | Direct `.npz` on disk (no Python) |
| `encodeNpzCompressed` / `writeNpzCompressed` | Low-level numpy ZIP writer |

**`web/molgui_webgpu/MMParams.js`**

| Export | Role |
|--------|------|
| `parseDihedralTypes` | Load `DihedralTypes.dat` |
| `getDihedralParams` / `_uffDihedralFallback` | Torsion params or UFF sp³–sp³ |
| `getAngleParamsOptional` | Fallback to `Ass` without log spam |
| `dihedralEndDistanceL0` | 3D equilibrium \|a−d\| |
| `convertDihedralToDistance` | \(k_\phi \to k_r\) via \(dl/d\phi\) |

### NPZ schema (`*_topology.npz`)

Compatible with `pyBall/OCL/LFFSolver.upload_state(neighs=…, KLs=…)`:

| Key | Shape | dtype |
|-----|-------|-------|
| `pos` | (N, 3) | float64 |
| `Z` | (N,) | int32 |
| `neigh_idx` | (N, M) | int32, **−1 padded** |
| `neigh_count` | (N,) | int32 |
| `bond_l0`, `bond_k` | (N, M) | float32 |
| `KLs` | (N, M, 2) | float32 |
| `stick_class` | (N, M) | int32 (0 empty, 1 bond, 2 angle, 3 dihedral) |
| `max_neighbors`, `natoms`, `n_bond`, `n_angle`, `n_dihedral` | (1,) | int32 |

Written **directly in Node** by `LinearizedTopologyNpz.writeTopologyNpzFile()` → `npzWrite.encodeNpzCompressed()`. No `elements` / `source` object arrays in JS output (viewer uses `Z`; LFFSolver needs numeric arrays only).

### HTML viewers

- **Colors:** K₁₂ red, K₁₃ green, K₁₄ blue.
- **Legend:** clickable checkboxes — toggle each stick class independently.
- **Atom sizes:** H smaller (r≈0.07), C/Si larger (0.15 / 0.18).
- **`debug_viewer.html`:** open directly (`file://` OK); data baked in.
- **`npz_viewer.html`:** use **file picker** on `file://`; or serve directory (`python3 -m http.server`) for auto-fetch of `{base}_topology.npz`.

Viewer **logic lives in** `LinearizedTopologyViewer.js` (`exportStickViewer*`); per-structure HTML under `tests/tSiNCs/linearized/` are **generated artifacts**, not source.

### Validation results (fixtures)

| Structure | Atoms | K₁₂ | K₁₃ | K₁₄ | max degree | budget |
|-----------|-------|-----|-----|-----|------------|--------|
| adamantane | 26 | 28 | 60 | 96 | 20 | 48 |
| si_G1_caps_only | 214 | 287 | 720 | 1380 | 42 | 48 |
| si_G2_facet111 | 609 | 924 | 2478 | 4872 | 42 | 48 |

Graph-distance checks: K₁₂ → dist **1**; K₁₃ → dist **2**; K₁₄ (dihedral) → dist **3** (MMFF 1–4).

**Python parity (adamantane K₁₂+K₁₃):** counts and \((l_0, K)\) match `pyBall/OCL/MMFFL.py`; atom indices differ between Python `AtomicSystem` and JS `parseMol2` loaders (isomorphic topology). **K₁₄ not yet in Python MMFFL.**

### Milestone status

| Milestone | Status |
|-----------|--------|
| M-L0 | **Pass** (adamantane) |
| M-L1 | **Pass** (G1/G2 graph-dist validation) |
| M-L2 | **Pass** (K₁₄ at graph dist 3; see note above) |
| M-L3 | **Pass** (HTML viewers + toggles) |
| M-L4 | **Partial** (counts/values OK; index order differs mol2 loaders) |
| M-L5 | **Pass** (max degree 42 < 48 on G2) |
| M-L6 | **Open** (analytical LFF Hessian vs FD; integration job) |

### Git — add these reusable sources (not generated outputs)

| Path | Action |
|------|--------|
| `scripts/build_linearized_topology.mjs` | **new** — headless batch CLI |
| `web/common_js/npzWrite.js` | **new** — numpy `.npz` writer (Node) |
| `web/molgui_webgpu/LinearizedTopologyNpz.js` | **new** — topology → `.npz` |
| `web/molgui_webgpu/LinearizedTopologyViewer.js` | **new** — Three.js HTML viewers |
| `web/molgui_webgpu/MMFFLTopology.js` | **modified** — K₁₄, packing, spring list |
| `web/molgui_webgpu/MMParams.js` | **modified** — dihedral types + distance linearization |

**Removed:** `scripts/topology_to_npz.py`.

**Do not commit:** `tests/tSiNCs/linearized/*` (`.json`, `.npz`, `*_debug_viewer.html`, `*_npz_viewer.html`), `tests/tSiNCs/fixtures/vibration_parallel/` (local fixtures).

Optional later: commit a **single** generic `scripts/linearized_topology_npz_viewer.template.html` if we want a static viewer without re-running the batch script; today the reusable viewer is **`exportStickViewerNpzLoaderHTML()`** in `LinearizedTopologyViewer.js`.

### Remaining / integration

1. Port K₁₄ to `pyBall/OCL/MMFFL.py` for full JS↔Python parity.
2. Wire `*_topology.npz` into vibration Hessian assembly (integration with sparse solver job).
3. Resolve mol2 atom-order parity or standardize on `.xyz` + `recalculateBonds` for golden tests.
4. `LFF.cl` `MAX_NEIGHBORS` still 8 — vibration kernel fork (`vib_jacobi.cl`) should use ≥48 per this packing.
