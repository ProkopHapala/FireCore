# Nanocrystal Generation — Design & Progress

**Scope:** Structure generation for Si/C nanocrystals (`tests/tSiNCs/nanocrystals.mjs generate`) with **correct tetrahedral passivation** and a **minimal reproducible test matrix** before large-scale sampling.

> **Canonical scripts (2026-07):** `tests/tSiNCs/` — see [`tests/tSiNCs/AGENTS.md`](../../../tests/tSiNCs/AGENTS.md). Legacy `scripts/` copies deprecated.

**Related:** [`gen_nanocrystals.chat.md`](gen_nanocrystals.chat.md), [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../topical_audit/Nanocrystal_Vibrations.md), `web/molgui_webgpu/EditableMolecule.js` (`addCappingAtoms`, `missingDirsVSEPR`).

**Status:** **Implemented** (2026-06-15) — M-G0–M-G2 PASS; JS/Python unified; cross-check harness green.

**Fixtures:** `Parallel_agent_fixtures.guide.md` — bootstrap before starting.

---

## Design goals

1. **Minimal geometry first:** Every bulk Si/C atom has four σ-neighbors in a tetrahedral arrangement (or explicit H caps completing tetrahedra on surfaces). Test on the smallest structures that expose the bug, not on 500-atom random shapes.
2. **Inspectable output:** Each test case writes `.xyz` (and optionally `.mol2`) suitable for Jmol / OVITO visual inspection, with a one-line geometry sanity summary (Si–H distances, H–Si–H angles, clash count).
3. **Defect catalog:** Document and exercise every shape/defect knob in `gen_nanocrystals.mjs` so passivation can be validated under controlled perturbations.
4. **Separation of concerns:** Generation produces **topologically correct** capped lattices; optional MMFF relax is a **downstream** step (not mixed into defect logic until baseline passivation is fixed).
5. **Fail loudly:** If valence cannot be satisfied (unsupported `totalDomains`, zero-length bond vectors), abort with atom id and local environment — no silent skip.

---

## Root-cause analysis (why H looks wrong today)

| Factor | Mechanism | Severity |
|--------|-----------|----------|
| **No post-cap relax** | Raw VSEPR placement only; Python vibration tests always `MMFF.run()` before Hessian | High (visual) |
| **Distorted surface Si–Si–Si angles** | Facet cuts break bulk 109.5°; VSEPR completes tetrahedron **locally** relative to distorted bonds → H may not align with facet outward normal | Medium |
| **`insertBridge` after capping** | Carbon-default geometry (`hDist: 1.3`, `upOffsetFactor: 0.5`), atom promoted to Si; adds non-tetrahedral SiH₂ | High when `insertProb > 0` |
| **`collapseBridgeAt`** | Removes bridge + H, bonds neighbors — changes topology after caps | Medium |
| **No explicit `atype` assignment** | Falls back to `ElementTypes` → `AtomTypes['Si']` (usually OK, fragile) | Low |
| **Prune before/after cap** | `minSiDegree=2` removes Si with &lt;2 Si neighbors; can leave exotic undercoordination before cap | Medium on small cuts |
| **Bond recalc after cap** | May add spurious H–H or long Si–Si bonds if `bondFactor`/`defaultRcut` too loose | Medium |

**Key insight:** `missingDirsVSEPR()` is mathematically correct for ideal VSEPR completion. Failures are mostly **pipeline ordering** (bridge ops, no relax, surface distortion), not a wrong cross-product formula.

---

## Minimal test geometries (ordered by complexity)

Use these before any random `--samples` batch.

| ID | System | Atoms (approx) | Purpose |
|----|--------|----------------|---------|
| **G0** | Single SiH₄-like unit (manual `.xyz`) | 5 | Reference tetrahedron |
| **G1** | `diamond_primitive.xyz` 1×1×1, no cut, H-cap all undercoordinated | ~10–20 | Bulk corner / single cell |
| **G2** | `Si-sym.cif` replicate `nx=ny=nz=1`, `--centered 0`, single `{111}` plane, large `planeCScale` (gentle cut) | ~30–80 | One facet, few surface types |
| **G3** | G2 with `--insertProb 0 --collapseProb 0` (caps only) vs defaults | compare | Isolate bridge defects |
| **G4** | Adamantane (`cpp/common_resources/mol/adamantane.mol2`) converted to Si analogue OR C reference | ~26 | Known tetrahedral cage; cross-check cap logic on organic |
| **G5** | `diamond_tersoff_2x2x2.xyz` or `2×2×2` supercell from CIF, spherical/no plane cut | ~50–100 | Regular interior + corners |

**Acceptance metrics (automated, per structure):**

- For each Si/C with 4 σ-domains: all four X–Si–Y angles in **[100°, 120°]** (bulk) or document surface exception.
- Si–H bond lengths in **[1.40, 1.55] Å** (UFF/MM estimate ~1.48 Å).
- No pair distance &lt; **1.8 Å** (H–H clash) unless flagged as intentional defect test.
- Surface Si: `nDang` caps match `4 - n_real_sigma_neighbors`.

---

## Defect & shape modification catalog

All knobs in `tests/tSiNCs/nanocrystals.mjs generate` (or deprecated `tests/tSiNCs/gen_nanocrystals.mjs`):

### A. Shape / size (Miller cutting)

| Option | Effect |
|--------|--------|
| `--nx-range`, `--ny-range`, `--nz-range` | Supercell replication before cut |
| `--centered 0\|1` | Centered vs corner replication |
| `--planeTemplates` | e.g. `a111`, `a100`, `a110` — Miller families via `CrystalUtils.expandPlaneTemplates` |
| `--planeSymC` | Symmetry expansion factor for plane set |
| `--planeMode ang\|frac` | Plane offset definition |
| `--planeC0`, `--planeCScale`, `--planeCJitter` | Cut position / randomness (`cmin`, `cmax` from average lattice extent) |
| `--cif` | Source cell (`Si-sym.cif` default; can use `diamond_primitive.cif`) |

**Test matrix suggestion:** fix replication, sweep `{111}` only → `{111}+{100}` → jitter on/off.

### B. Topology cleanup

| Option | Effect |
|--------|--------|
| `--minSiDegree` | Prune Si with fewer than N **Si** neighbors (default 2) |
| `--pruneMaxIter` | Iterative prune passes |
| `--bondFactor`, `--defaultRcut` | Bond detection thresholds in `recalculateBonds` |

### C. Passivation

| Option | Effect |
|--------|--------|
| `--caps H\|0\|none` | `addCappingAtoms(mm, cap, { onlySelection: false, bBond: true })` |
| `--resolveClashes 0\|1` | Rotate clashing H caps around parent heavy atom (default **1**); use **0** for geometry validation |
| `--outwardBias` | Blend cap directions toward outward COG (0–1) |
| `--capHHBonds 0\|1` | **Optional:** add explicit H–H bonds for cap pairs closer than `--capHHBondDist` (default **0** = off) |
| `--capHHBondDist` | Distance threshold (Å) for cap H–H bonding (default **1.8**) |

**Cap H–H bonds (vibration):** `recalculateBonds` uses covalent H–H cutoff ≈ 0.96 Å, so CH₂ H–H (~1.77 Å) and inter-surface cap clashes are **not** bonded automatically. With `--capHHBonds 1`, the generator adds topological H–H bonds after all topology edits (post-cap, post-bridge, **after** last `recalculateBonds`). Stiffness comes from `BondTypes.dat` (`H H 1 0.740 10.000`) and can be overridden later for the vibration/LFF pipeline — goal is to replace expensive non-bonded H–H repulsion with bonded terms. Python sphere mode writes pair list to `*.hh_bonds.json`; planes mode forwards flag to Node.

### D. Surface bridge defects (post-cap)

| Option | Effect |
|--------|--------|
| `--insertProb` | On surface Si–Si bonds: `insertBridge` → SiH₂-like bridge, Z→14 |
| `--collapseProb` | On surface SiH₂-like: `collapseBridgeAt` (remove bridge + H, bond neighbors) |
| `--collapseAll` | Legacy `collapseAllBridges` (carbon path) |
| `--requireH2` | Filter for collapse candidates |

**Important:** Bridge section deliberately **does not** recalc bonds after collapse (preserves R1–R2). Insert still uses geometry defaults from `MoleculeUtils.insertBridge`.

### E. Sampling / output

| Option | Effect |
|--------|--------|
| `--samples`, `--seed`, `--maxFiles` | Stochastic batch |
| `--stackedXyzOut` | Multi-frame XYZ for animation |
| `--statsCsv` | Counts: SiH, SiH₂, SiH₃, bare, bridge, pseudo-energy |

### F. Not yet in generator (future)

- Explicit outward **facet normal** bias for H (chemically preferred over pure VSEPR).
- Post-cap **MMFF/UFF relax** hook (Python subprocess or separate script).
- **Adamantane-shaped** Si cage template without random plane jitter.

---

## Planned work packages

### WP-G1 — Diagnostic harness (headless Node)

**Deliverable:** `tests/tSiNCs/test_nanocrystal_geometry.mjs` (or flags on existing script)

- Input: fixed CLI presets (G0–G5 table above).
- Output: `.xyz` + `geometry_report.json` (angles, distances, clash count).
- No MMFF dependency — pure geometry audit.

### WP-G2 — Isolate hydrogen placement

1. Run G2/G3 with `insertProb=collapseProb=0`.
2. Compare VSEPR output to Python builder in `test_nanocrystal_sparse_hessian.py` (`build_diamond_nanoparticle` H placement) on same topology.
3. If VSEPR matches Python but both look bad → problem is **missing relax**, not JS cap logic.
4. If JS differs → debug `missingDirsVSEPR` edge cases (`nb=1`, `nb=2` on high-curvature facets).

### WP-G3 — Bridge defect regression

- Fixed seed, `insertProb=0.2`, `collapseProb=0.3` vs zeros.
- Document which surface motifs break tetrahedra (SiH₂ bridge insertion geometry).

### WP-G4 — Explicit typing + optional outward bias

- Before cap: `setAtomTypeByName(id, 'Si')` / `'Si3'` for all Si.
- Optional: bias last cap direction toward **outward** = negative mean surface normal (estimate from undercoordinated Si only).

### WP-G5 — Minimal crystal library

- Commit reference `.xyz` for G1, G2 (caps-only), G3 variants under `tests/tSiNCs/geometry/` (or `OUT_nanocrystals_ref/`).
- Used by linearized-topology and solver milestones.

---

## Testing milestones

| Milestone | Criterion | Status |
|-----------|-----------|--------|
| **M-G0** | G0 manual SiH₄ passes all angle/length checks | **PASS** |
| **M-G1** | G1 primitive cell: every Si has 4 neighbors (Si or H), tetrahedral angles ±5° | **PASS** |
| **M-G2** | G2 single-facet cut, caps-only: no H–H clashes &lt; 1.8 Å; visual Jmol OK | **PASS** |
| **M-G3** | G3 proves bridge ops are **opt-in** defects; caps-only path is default for spectroscopy | **PASS** (default `collapseProb=0`; G3 preset exercises bridges) |
| **M-G4** | Automated `geometry_report.json` in CI (Node, no GPU) | harness exists; CI hook pending |
| **M-G5** | Optional MMFF relax reduces max force &lt; 1e-3 eV/Å on G2 (Python wrapper) | pending |

---

## Dependencies & handoff

- **→ Linearized topology:** Needs G2 `.mol2`/`.xyz` with clean tetrahedral graph.
- **→ Sparse solver:** Needs relaxed structures eventually; solver milestone can start on adamantane / diamond primitive before G2 is perfect.

---

## Open questions

1. Should spectroscopy pipeline **require** relaxed structures (Python), or should generator call relax internally?
2. Prefer **C** (diamond) or **Si** for first minimal crystals? (Si is target; C has more reference data.)
3. Is `minSiDegree=2` correct for small facets, or should minimal tests use `minSiDegree=1` with aggressive prune disabled?

---

## Parallel agent contract (Job 1)

Full rules: **`Parallel_agent_fixtures.guide.md`**.

| | |
|--|--|
| **Fixture root** | `tests/tSiNCs/fixtures/vibration_parallel/` (gitignored; run `python3 tests/tSiNCs/bootstrap_vibration_parallel_fixtures.py`) |
| **Read** | `structures/diamond_primitive.xyz`, `structures/si_G1_caps_only.*`, `structures/si_G2_facet111_caps_only.*` |
| **May edit** | `tests/tSiNCs/gen_nanocrystals.mjs`, `tests/tSiNCs/test_nanocrystal_geometry.mjs`, `scripts/mol2_to_xyz.mjs`, `EditableMolecule.js` / `MoleculeUtils.js` (capping/bridges only) |
| **May write** | `tests/tSiNCs/geometry/` (new outputs; do not overwrite fixture basenames) |
| **Must NOT edit** | `FTIR.py`, `OCL/*`, `MMFFLTopology.js`, `LFF.cl`, `vib_jacobi.cl`, `test_vibration*` |

---

*Last updated: 2026-06-15 — planning + parallel contract.*

---

## Daily report — 2026-06-15 (implementation wrap-up)

### Summary

Nanocrystal generation is **working end-to-end** for both **Si** and **C** diamond/sphalerite lattices. Two CLIs share the same tetrahedral passivation physics; a headless geometry harness and a Python cross-check script verify parity. Milestones **M-G0, M-G1, M-G2** pass; bridge defects are opt-in (`collapseProb` default **0**).

### What was implemented

#### 1. Root-cause fixes (passivation geometry)

| Bug | Fix |
|-----|-----|
| Si–H length ~1.62 Å (hardcoded) | `capBondLength()` / `BondTypes.dat` → **1.46 Å** (Si), **1.095 Å** (C) |
| H–H angles ~53° on nb=2 caps | Replaced `d = v ± 0.5·p` with **`missingDirsVSEPR` / `missing_dirs_vsepr_tetra`** (tetrahedral completion) |
| Wrong geometry audit | H–Si–H angle measured **at the heavy atom**, not at the H neighbor |
| `outwardBias` Si-only | Uses **all non-H atoms** for outward COG estimate (C caps work) |
| `resolveCapHClashes` C failure | `parentHeavyIndex` finds **any heavy parent**, not Si-only |
| Python Si sphere: no caps | `default_rcut_heavy`: **C 1.75 Å**, **Si 2.55 Å** (Si–Si ~2.35 Å) |
| Primitive CIF symOps error | `--applySymmetry 0` for `diamond_primitive.cif` / `Si_primitive.cif` |

#### 2. Dual generator architecture

**JS** — `tests/tSiNCs/gen_nanocrystals.mjs` (full feature set, native)

- **Miller plane cut** (`--cutMode planes`, default): CIF → replicate → `CrystalUtils` plane filter → prune → cap → optional bridges.
- **Spherical cut** (`--cutMode sphere`): replicate cell, center, filter `r² ≤ R²`.
- **Element-agnostic**: heavy element from CIF (`heavyZ`); same VSEPR for C and Si.
- **New flags**: `--sphereR`, `--sphereNrep`, `--rcutHeavy`, `--minHeavyDegree`, `--resolveClashes`, `--outwardBias`.
- Outputs **`.mol2` + `.xyz`** per sample.

**Python** — `tests/tSiNCs/gen_nanocrystals.py` + `pyBall/nanocrystal_gen.py`

- **Sphere cut**: native (`build_spherical_nanoparticle`) — parity with JS sphere mode.
- **Miller planes + bridges**: delegates to Node (`gen_nanocrystals.mjs` subprocess).
- **`test_nanocrystal_sparse_hessian.py`** refactored to import `pyBall.nanocrystal_gen` (single VSEPR source).

#### 3. Shared cap pipeline (both languages)

```
replicate lattice → cut (planes | sphere) → recalculateBonds
  → prune undercoordinated heavy (min degree 2)
  → assignCrystalAtomTypes
  → addCappingAtoms (VSEPR + BondTypes lengths + optional outwardBias)
  → [optional] resolveCapHClashes (default on; --resolveClashes 0 for geometry validation)
  → [optional] insertBridge / collapseBridge (planes only, opt-in)
  → recalculateBonds (last pass)
  → [optional] addCapHHBonds (--capHHBonds 1; after last recalc so bonds persist)
  → write .mol2 / .xyz
```

**VSEPR reference:** `web/molgui_webgpu/EditableMolecule.js` `missingDirsVSEPR` ↔ `pyBall/nanocrystal_gen.py` `missing_dirs_vsepr_tetra`.

**Clash resolve note:** Default clash resolver rotates H around its parent to separate H–H pairs &lt; 1.8 Å. This fixes dense facet cuts but can distort local tetrahedral angles when many **inter-surface** H clashes exist (e.g. 42 on C sphere R=6). For geometry validation / cross-check, use `--resolveClashes 0`; for production facet cuts with downstream MMFF relax, keep default **on** (G2: 0 clashes with resolve on).

#### 4. Test & validation harnesses

| Script | Role |
|--------|------|
| `tests/tSiNCs/test_nanocrystal_geometry.mjs` | Presets **G0–G5**; writes `geometry_report.json` + `.xyz`; no MMFF |
| `tests/tSiNCs/crosscheck_nanocrystal_generators.py` | Py vs JS parity: C/Si sphere R=6, Si planes G2, C planes G1, Si bridges demo |
| `scripts/mol2_to_xyz.mjs` | Thin mol2→xyz converter for cross-check |

**Cross-check results** (all PASS, `tests/tSiNCs/crosscheck/crosscheck_report.json`):

| Case | Atoms (py = js) | Key metric |
|------|-----------------|------------|
| C sphere R=6 | 270 | H–H ≈ 109.47° both; radial RMSD ≈ 2.3×10⁻⁵ Å |
| Si sphere R=6 | 80 | identical counts |
| Si planes G2 | 609 | py delegates to node; identical |
| C planes G1 | 71 (js) | H–H ≈ 109.47° |
| Si bridges demo | 615 | 18 collapsed, 20 inserted |

**Geometry milestones:**

- **M-G0** PASS — manual SiH₄ reference
- **M-G1** PASS — primitive cell caps
- **M-G2** PASS — `node tests/tSiNCs/test_nanocrystal_geometry.mjs --preset G2` → Si–H 1.460 Å, 0 clashes

#### 5. How to run

```bash
# Geometry audit (G0–G5)
node tests/tSiNCs/test_nanocrystal_geometry.mjs --preset G2

# Full py/js cross-check
python3 tests/tSiNCs/crosscheck_nanocrystal_generators.py

# JS — Si facet, caps only
node tests/tSiNCs/gen_nanocrystals.mjs --cutMode planes --nx-range 2,2 --ny-range 2,2 --nz-range 2,2 \
  --planeTemplates a111 --planeCScale 0.40 --caps H --insertProb 0 --collapseProb 0 \
  --outDir OUT_nanocrystals --prefix si_nc

# JS — C sphere
node tests/tSiNCs/gen_nanocrystals.mjs --cutMode sphere --sphereR 6 --sphereNrep 5 \
  --cif cpp/common_resources/crystals/diamond_primitive.cif --applySymmetry 0 \
  --caps H --resolveClashes 0 --outDir OUT_nanocrystals --prefix C_sphere

# Python — Si sphere (native)
python3 tests/tSiNCs/gen_nanocrystals.py --cutMode sphere --element Si \
  --sphere-r 6 --sphere-nrep 5 --outDir OUT_nanocrystals_py --prefix Si_sphere

# Vibration-oriented: explicit cap H-H bonds (skip NB H-H repulsion in LFF later)
node tests/tSiNCs/gen_nanocrystals.mjs --cutMode planes --nx-range 2,2 --ny-range 2,2 --nz-range 2,2 \
  --planeTemplates a111 --planeCScale 0.40 --caps H --capHHBonds 1 --capHHBondDist 1.8 \
  --outDir OUT_nanocrystals --prefix si_nc_hhbond
```

#### 6. Output locations (generated artifacts — not for git)

| Directory | Contents |
|-----------|----------|
| `tests/tSiNCs/crosscheck/` | Cross-check `.xyz`, `.mol2`, `crosscheck_report.json` |
| `tests/tSiNCs/geometry/` | G0–G5 harness outputs, `*_geometry_report.json` |
| `OUT_nanocrystals/` | Default JS CLI output |
| `OUT_nanocrystals_py/` | Default Python CLI output |

Canonical comparison files: `tests/tSiNCs/crosscheck/C_sphere_R6_{py,js}.xyz`, `Si_sphere_R6_{py,js}.xyz`, `Si_planes_G2_{py,js}.xyz`, `Si_bridges_demo.xyz`.

### Feature parity matrix (final)

| Feature | JS | Python |
|---------|----|--------|
| Miller plane cut | native | Node delegate |
| Spherical cut | native | native |
| VSEPR H caps | ✓ | ✓ |
| BondTypes bond lengths | ✓ | ✓ |
| outwardBias | ✓ | ✓ |
| Clash resolve (optional) | ✓ | ✓ |
| Cap H–H bonds (optional) | ✓ | ✓ (sphere: `*.hh_bonds.json`; planes: Node) |
| Bridges insert/collapse | ✓ (planes) | Node delegate |
| C + Si | ✓ | ✓ |
| `.mol2` export | ✓ | via Node (planes) / xyz only (sphere) |

### Known limitations / next steps

1. **Python Miller cut** is not a native port of `CrystalUtils` — still calls Node for planes/bridges.
2. **Clash resolve** vs **tetrahedral angles**: trade-off on dense spheres; MMFF relax remains downstream (M-G5).
3. **Fixture regen**: `diamond_nc_R6_init.xyz` in vibration fixtures predates VSEPR fix; re-bootstrap when Hessian pipeline needs updated init geometry.
4. **M-G3–M-G5**: bridge regression documented (G3 preset); CI hook for `geometry_report.json` not yet wired.
5. **Cap H–H bond stiffness:** topology bonds use default `BondTypes.dat` H–H (`l0=0.74 Å`, `k=10`); for vibration, tune `l0`/`k` per actual cap H–H distance (~1.77 Å on CH₂) when wiring into LFF/MMFFL.

#### 7. Cap H–H bonds (2026-06-15 add-on)

**Motivation:** Dense capped nanocrystals have many H–H pairs in the 1.7–1.8 Å range (same-parent CH₂ + inter-surface clashes). Standard `recalculateBonds` does not create these (H–H covalent cutoff ≈ 0.96 Å), so vibration/Hessian paths need non-bonded H–H terms. Optional `--capHHBonds 1` adds explicit bonds so repulsion can be modeled as **bonded** interactions (faster sparse solvers, no NB pair list for those pairs).

**Implementation:**

| Layer | Function | Notes |
|-------|----------|-------|
| JS | `EditableMolecule.addCapHHBonds(maxDist)` | Scans H–H pairs; skips already-bonded; called **after** last `recalculateBonds` in `gen_nanocrystals.mjs` |
| Python | `find_cap_hh_pairs`, `save_hh_bonds_json` | Sphere CLI writes `*.hh_bonds.json` sidecar (0-based index pairs) |

**Verified (C sphere R=6, `--resolveClashes 0`):** 381 → 417 bonds; 42 H–H bonds total (36 newly added; 6 already present from sub-0.96 Å detection). Si G2 facet with clash resolve: 0 new H–H bonds (geometry clean).

**Code:** `web/molgui_webgpu/EditableMolecule.js`, `pyBall/nanocrystal_gen.py`, flags on `tests/tSiNCs/gen_nanocrystals.{mjs,py}`.

### Files to add to git (reusable code only)

**Core modules & generators**

| Path | Role |
|------|------|
| `tests/tSiNCs/gen_nanocrystals.mjs` | JS generator (modified: sphere mode, C/Si, resolveClashes, outwardBias) |
| `tests/tSiNCs/gen_nanocrystals.py` | Python CLI (sphere native; planes → Node) |
| `pyBall/nanocrystal_gen.py` | Shared Python builder: VSEPR, sphere cut, xyz I/O |
| `web/molgui_webgpu/EditableMolecule.js` | Cap pipeline + `addCapHHBonds` (modified) |

**Test & utility scripts**

| Path | Role |
|------|------|
| `tests/tSiNCs/test_nanocrystal_geometry.mjs` | G0–G5 geometry harness |
| `scripts/mol2_to_xyz.mjs` | mol2 → xyz converter |
| `tests/tSiNCs/crosscheck_nanocrystal_generators.py` | Py/JS parity cross-check |
| `tests/tMMFF/test_nanocrystal_sparse_hessian.py` | Refactored to use `nanocrystal_gen` |

**Suggested `git add` command**

```bash
git add \
  tests/tSiNCs/gen_nanocrystals.mjs \
  tests/tSiNCs/gen_nanocrystals.py \
  tests/tSiNCs/test_nanocrystal_geometry.mjs \
  scripts/mol2_to_xyz.mjs \
  pyBall/nanocrystal_gen.py \
  web/molgui_webgpu/EditableMolecule.js \
  tests/tSiNCs/crosscheck_nanocrystal_generators.py \
  tests/tMMFF/test_nanocrystal_sparse_hessian.py
```

**Do not add** (generated outputs): `tests/tSiNCs/crosscheck/*.{xyz,mol2,json}`, `tests/tSiNCs/geometry/**/*.{xyz,mol2,json}`, `OUT_nanocrystals*/`.

**Optional** (reusable viewers from related linearized-topology work, not nanocrystal-gen core): `tests/tSiNCs/linearized/adamantane_npz_viewer.html`, `tests/tSiNCs/linearized/si_G2_facet111_caps_only_debug_viewer.html` — add separately if NPZ/topology inspection is in scope.

---

*Report closed: 2026-06-15. Cap H–H bonds documented same day.*
