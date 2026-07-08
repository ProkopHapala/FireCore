# USER

plese read 
@doc/Topics/FTIR_Nanocrystals/Nanocrystal_generation.progress.md @doc/Topics/FTIR_Nanocrystals/Linearized_topology.progress.md @doc/Topics/FTIR_Nanocrystals/Nanocrystal_vibration_sparse.progress.md 


OK, now I realized the solution of quitelarge nancrystals is quite fast, so we should not spend much time on optimizing the diagonalization od dynamical matrix and just use exiring methods. eigh (dense)

we should rather focus on making full end-to-end pipleline form generation of some ensambles of crystals to theis fiagonalization and producing the spectra.


Spectra can be in fist approximation produced just as histogram of eignestate (we can then blur that histogram by convolution or something), we will sum spectra for many nanocrystals.

For debugging it would be good to have option to plot geoemtry of nanocrystas in some projection (orthographic) with arbitrary rotation matrix into .svg then we can plot for every nanocrystal generate the initial generated geometry, the relaxed geometry, and the spectrum. But this is just debug run, in production run we will not plot every crystal as it takes longer than clauclation of the spectra. We will acumulate the spectrum.

We should have some system how to store (cache) relaxed crystal geometry, topology and hessian in .npz files for fast load (parsing .xyz or other text format is slow, we wnat high throguhput system in the end, but with some diganostic tools for debuggin)

we should also think about how to make robust and flexible system for generation. We should choose some rnges of crystal sizes, make them more prolonged (disk like, tube like) be choosing different cutoff shapes for generation, and some probability of indroducing defects - colapsing -CH2- bridges, inserging -CH2- bridges. Any other ideas?

now think about it, make notes and plan what we will use from what is already implemented, what we need to implemt new and how the whole system should work including the data managemnt.

do not program yet


---

# Cursor Composer 2.5

# End-to-End Nanocrystal Spectrum Pipeline — Plan & Notes

**Docs:** [`README.md`](README.md) · [`doc/topical_audit/Nanocrystal_Vibrations.md`](../../topical_audit/Nanocrystal_Vibrations.md) · [`tests/tSiNCs/README.md`](../../../tests/tSiNCs/README.md)

> **Canonical scripts (2026-07):** `tests/tSiNCs/nanocrystals.mjs` (`generate`, `ensemble`, …). Legacy `scripts/` copies are deprecated. Run examples below with `tests/tSiNCs/` paths.

**Strategic pivot:** Dense diagonalization (`np.linalg.eigh` + mode/histogram spectrum) is fast enough for the target crystal sizes (~500–1500 atoms). **Do not invest further in sparse/GPU frequency-domain solvers** until ensemble throughput is blocked by `eigh` itself. The bottleneck to watch is **Hessian assembly** (MMFF FD), not diagonalization.

## 1. What we already have (reuse as-is)

| Layer | Status | Key artifacts |
|-------|--------|---------------|
| **Generation** | Done (M-G0–G2) | `tests/tSiNCs/nanocrystals.mjs generate`, `pyBall/nanocrystal_gen.py`, G0–G5 harness, bridge/cap knobs |
| **Linearized topology** | Done (M-L0–L5) | `tests/tSiNCs/nanocrystals.mjs topology` → `*_topology.npz` (K₁₂/K₁₃/K₁₄ packed) |
| **MMFF relax + Hessian** | Working | `MMFF.run()`, `getHessian3Nx3N`, `getHessianSparseBlocks` |
| **Dense vibration** | Working | `FTIR.project_rigid_modes`, `FTIR.vibration_spectrum_from_modes`, `test_vibration_spectra.py` |
| **Fixtures/bootstrap** | Working | `tests/tSiNCs/bootstrap_vibration_parallel_fixtures.py`, benchmark NPZ guide |
| **Size calibration** | Documented | `Vibration_benchmark_npz.guide.md`: dense `eigh` &lt;1 s up to ~558 atoms / 1674 DOF |

**Deferred (not on critical path):** GPU Jacobi (`vib_jacobi.cl`), `mechanical_greens_probing_sparse`, frequency scanning — all documented in `Sparse_vibration_solver.progress.md` but not needed for v1 production.

## 2. Target architecture

```mermaid
flowchart TB
    subgraph config [Ensemble config]
        C[ensemble.yaml: size ranges, cut shapes, defect probs, seed]
    end

    subgraph gen [Generation]
        G[nanocrystals.mjs generate / py] --> INIT[init crystal]
    end

    subgraph relax [Relaxation]
        INIT --> R[MMFF.run bonds+angles]
        R --> REL[relaxed geometry]
    end

    subgraph hess [Hessian - choose path]
        REL --> H_A[Path A: MMFF FD dense H]
        REL --> TOPO[build_linearized_topology.mjs]
        TOPO --> H_B[Path B: analytical LFF K assembly]
    end

    subgraph vib [Vibration - dense only]
        H_A --> EIGH[np.linalg.eigh on mass-weighted H]
        H_B --> EIGH
        EIGH --> MODES[omega_i modes]
        MODES --> HIST[histogram + optional Gaussian blur]
    end

    subgraph accum [Ensemble accumulation]
        HIST --> SUM[running sum spectrum]
        SUM --> OUT[ensemble_spectrum.npz]
    end

    subgraph cache [NPZ cache layer]
        INIT -.-> NPZ1[init.npz]
        REL -.-> NPZ2[relaxed.npz]
        TOPO -.-> NPZ3[topology.npz]
        H_A -.-> NPZ4[hessian.npz]
        HIST -.-> NPZ5[spectrum.npz]
    end

    subgraph debug [Debug mode only]
        INIT --> SVG1[init.svg]
        REL --> SVG2[relaxed.svg]
        HIST --> PLOT[spectrum.png]
    end

    C --> G
```

**Two Hessian paths (sequential adoption, not parallel investment):**

| Path | Physics | Hessian cost | When to use |
|------|---------|--------------|-------------|
| **A — MMFF FD** | Full angles/dihedrals | O(N) FD calls, expensive | Validation, first ensemble runs, parity |
| **B — LFF analytical** | Spring network K₁₂+K₁₃+K₁₄ | O(n_sticks) assembly, cheap | Production throughput once modes match Path A on 2–3 reference crystals |

Diagonalization is **always dense `eigh`** for v1. Path B only changes how `K` is built, not how it is solved.

## 3. Spectrum model (v1 approximation)

### Per-crystal spectrum

1. Build `K`, diagonal mass `M`, positions `pos`.
2. `K_proj = FTIR.project_rigid_modes(K, M, pos)`.
3. Mass-weight: `K_mw = M^{-1/2} K_proj M^{-1/2}`.
4. `eigvals, eigvecs = np.linalg.eigh(K_mw)` → `omega_i = sqrt(max(eigval, 0))`.
5. **Histogram:** bin `omega_i` into uniform grid `[ω_min, ω_max]` with width `dω` (e.g. 1–5 cm⁻¹ equivalent).
6. **Optional broadening:** convolve bin counts (or a delta-comb) with Gaussian kernel `σ` (instrument linewidth / thermal broadening).
7. Store per-crystal `spectrum_hist`, `omega_centers`, raw `omegas_modes` (for debugging).

This is simpler than `vibration_spectrum_from_modes` (no dipole weighting). Good first approximation; IR activity weighting `w_i ∝ |∂μ/∂Q_i|²` is phase 2.

### Ensemble accumulation

```
S_total(ω) = Σ_k S_k(ω)     # sum over crystals k
```

Optionally normalize by crystal count or by total mode count. Running accumulator in float64 avoids storing every per-crystal histogram in production.

**Existing code to wrap:** `FTIR.project_rigid_modes`, mass matrix from `FTIR.get_mass_matrix`, `np.linalg.eigh`. Add thin `modes_to_histogram(omegas, grid, sigma)` — do not reuse per-ω Green probing.

---

## 4. NPZ cache system (high-throughput data management)

Text formats (`.xyz`, `.mol2`) are for human inspection and generator I/O only. **All pipeline stages read/write NPZ.**

### Directory layout

```
{run_root}/
  ensemble_meta.json          # config snapshot, code version, param hash
  manifest.jsonl              # one JSON line per crystal: id, status, paths, timings, checksums
  crystals/
    {crystal_id}/
      meta.json                 # generation params (reproducible)
      01_init.npz
      02_relaxed.npz
      03_topology.npz           # optional until LFF path
      04_hessian.npz
      05_spectrum.npz
      debug/                    # only if --debug
        init.svg, relaxed.svg, spectrum.svg
  accumulated/
    spectrum_running.npz        # updated every N crystals
    spectrum_final.npz
```

### NPZ schemas

**`01_init.npz`**
| Key | Shape | Notes |
|-----|-------|-------|
| `pos` | (N,3) float64 | Generated geometry |
| `Z` | (N,) int32 | Atomic numbers |
| `elements` | (N,) str | Optional |
| `bonds_ij` | (B,2) int32 | Optional bond list |
| `gen_params` | scalar str | JSON blob: seed, cutMode, nx/ny/nz, planeTemplates, insertProb, … |

**`02_relaxed.npz`**
| Key | Notes |
|-----|-------|
| `pos` | Relaxed |
| `fmax` | Final max force |
| `n_steps` | MMFF convergence |
| `mmff_switches` | JSON for reproducibility |

**`03_topology.npz`**
Reuse schema from `Linearized_topology.progress.md` (`neigh_idx`, `KLs`, `stick_class`, …) — same as `{base}_topology.npz`.

**`04_hessian.npz`**
| Key | Notes |
|-----|-------|
| `H` | Full dense (3N×3N) — primary for v1 |
| `H_projected` | After rigid-mode shift |
| `M` | Diagonal mass |
| `pos` | At Hessian geometry |
| `source` | `"mmff_fd"` or `"lff_analytical"` |

For N &gt; ~600 atoms: still store dense `H` if `eigh` is fast; only omit when calibration says otherwise.

**`05_spectrum.npz`**
| Key | Notes |
|-----|-------|
| `omega_centers` | Histogram bin centers |
| `hist` | Binned mode counts |
| `omegas_modes` | Raw eigenfrequencies |
| `sigma` | Broadening used |

### Cache policy

| Operation | Behavior |
|-----------|----------|
| **Resume** | If `04_hessian.npz` exists and `meta.hash` matches, skip relax+hessian |
| **Invalidate** | Generation param change → invalidate downstream stages for that crystal |
| **Id** | `{seed:08d}_{index:06d}` or content hash of `gen_params` |
| **Manifest** | Append-only JSONL for crash recovery and timing stats |

Extend `bootstrap_vibration_parallel_fixtures.py` pattern into `scripts/nanocrystal_cache.py` (read/write helpers, manifest) — not the fixtures themselves.

## 5. Generation ensemble design

### Shape / size sampling (existing knobs)

| Goal | Mechanism (already in `gen_nanocrystals.mjs`) |
|------|--------------------------------------------------|
| **Size distribution** | `--nx-range`, `--ny-range`, `--nz-range` sampled independently |
| **Disk-like** | `nx >> ny ≈ nz` + plane cuts ⊥ long axis (e.g. `a100` on short dims) |
| **Rod / tube-like** | `ny, nz >> nx` or use **sphere** with ellipsoidal scaling (see new below) |
| **Faceted NC** | `--cutMode planes`, `--planeTemplates`, `--planeCScale`, `--planeCJitter` |
| **Spherical NC** | `--cutMode sphere`, `--sphereR`, `--sphereNrep` |
| **Surface defects** | `--insertProb` (add SiH₂ bridge), `--collapseProb` (remove bridge) |
| **Clean baseline** | `insertProb=0`, `collapseProb=0` (default for spectroscopy validation) |

### Recommended ensemble config (YAML sketch)

```yaml
ensemble:
  n_crystals: 1000
  element: Si
  seed_base: 42

  size:
    nx: { dist: uniform_int, min: 2, max: 5 }
    ny: { dist: uniform_int, min: 2, max: 5 }
    nz: { dist: uniform_int, min: 2, max: 5 }

  shape_class:  # categorical mixture
    - { weight: 0.4, cutMode: sphere, sphereR: { uniform: [4, 8] } }
    - { weight: 0.3, cutMode: planes, planeTemplates: [a111], aspect: disk }  # nx>ny,nz
    - { weight: 0.2, cutMode: planes, planeTemplates: [a100,a110], aspect: rod }
    - { weight: 0.1, cutMode: planes, planeTemplates: [a111,a100], multi_facet: true }

  defects:
  insertProb: { dist: bernoulli, p: 0.15 }   # per-bond probability
  collapseProb: { dist: bernoulli, p: 0.10 }

  passivation:
    caps: H
    capHHBonds: 1
    resolveClashes: 1
```

### Additional generation ideas (prioritized)

| Idea | Value | Effort | Notes |
|------|-------|--------|-------|
| **Ellipsoidal cut** | Disk/tube without awkward replication | Medium | `r² = (x/Rx)²+(y/Ry)²+(z/Rz)² ≤ 1` filter; not in generator yet |
| **Log-normal size** | Realistic NC size distribution | Low | Sample `sphereR` or `nx*ny*nz` from log-normal |
| **Stratified defects** | Control defect density bins | Low | Bucket by `n_bridges_inserted` from `statsCsv` |
| **Cap chemistry** | Surface chemistry variation | Medium | `--caps` variants (H vs D for isotope broadening) |
| **Isotope disorder** | Lifetime broadening in spectrum | Low | Random `M` scaling on H (cache in hessian stage) |
| **Twin / stacking fault** | Realistic faceting | High | Asymmetric plane offsets (`planeC0` per facet) |
| **Dangling bonds (no cap)** | Surface trap states | Low | `--caps 0` on subset — probably separate experiment |
| **Oxidized surface** | SiO₂-like passivation | High | New cap element + types — future |

**Recommendation:** v1 uses **mixture of sphere + single-facet + disk-like replication** with **low defect probability** (0–20%). Validate caps-only subset first, then ramp defects.

## 6. Debug visualization (optional, slow path)

**Production:** accumulate spectrum only — no geometry plots.

**Debug mode (`--debug` or `--debug-max N`):**

| Output | Content |
|--------|---------|
| `init.svg` | Orthographic projection of generated geometry |
| `relaxed.svg` | Same view after MMFF |
| `spectrum.svg` | Line plot or bar histogram of `omega_i` |

**Projection API (new, small module):**

```python
# pyBall/nanocrystal_viz.py (planned)
def project_orthographic(pos, Z, R_view, width=400, height=400) -> svg_str
# R_view: (3,3) rotation matrix (world → view)
# Atoms: circles colored by Z, radius by element
# Bonds: optional lines from bond list
```

Reuse conventions from `plotUtils.export_pov` (camera as `look_at`, `camera_up`, `camera_right`) but emit SVG directly — no POV-Ray dependency. Rotation matrix `R_view` sampled from ensemble config or fixed for comparability.

**Spectrum debug:** matplotlib `spectrum.png` or minimal SVG polyline — use `plotUtils` patterns from existing tests.

## 7. Orchestration: new components vs reuse

### Reuse (consolidated into unified CLI — 2026-06-22)

| Component | Role | Consolidated into |
|-----------|------|-------------------|
| `gen_nanocrystals.mjs` | Structure generation | `tests/tSiNCs/nanocrystals.mjs generate` (deprecated standalone in same folder) |
| `run_nanocrystal_ensemble.mjs` | Ensemble orchestration | `tests/tSiNCs/nanocrystals.mjs ensemble` (deprecated) |
| `build_linearized_topology.mjs` | Topology NPZ | `tests/tSiNCs/nanocrystals.mjs topology` (deprecated) |
| `test_nanocrystal_geometry.mjs` | Geometry audit | `tests/tSiNCs/nanocrystals.mjs audit` (deprecated) |
| `debug_nanocrystal_nonbond_groups.mjs` | Nonbond debug | `tests/tSiNCs/nanocrystals.mjs nonbond` (deprecated) |
| `MMFF.init/run/getHessian3Nx3N` | Relax + Hessian (Path A) | Python subprocess (unchanged) |
| `FTIR.project_rigid_modes`, `get_mass_matrix` | Standard vibration prep | Python subprocess (unchanged) |
| `np.linalg.eigh` | Mode frequencies | Python subprocess (unchanged) |
| `bootstrap_vibration_parallel_fixtures.py` | Pattern for artifact generation | Unchanged |

### Module updates (shared functions moved from scripts → existing modules)

| Module | Functions added | Source script (deprecated) |
|--------|----------------|---------------------------|
| `web/molgui_webgpu/MoleculeUtils.js` | `atomDist`, `atomAngleDeg`, `neighborIndices`, `bondedPairSet` | `test_nanocrystal_geometry.mjs` |
| `web/molgui_webgpu/Nanocrystals.js` | `isBulkSi`, `auditGeometry`, `oneLineSummary`, `buildSiH4`, `sampleValue`, `sampleCutSpec`, `crystalId`, `appendManifest`, `pyEnv`, `runCmd`, `runPy`, `writeTimingReport`, `genArgsFromConfig`, `genParamsFromAtlasEntry` | `test_nanocrystal_geometry.mjs`, `run_nanocrystal_ensemble.mjs` |
| `web/common_js/nanocrystalSvg.js` | `writeCrystalCompareSvgs`, `writeCrystalViewerJson`, `enrichViewerJsonWithTopology`, `installViewer` | `run_nanocrystal_ensemble.mjs` |
| `web/common_js/MolIO.js` | `loadMol` (generalized loader for .mol2 and .xyz) | `build_linearized_topology.mjs` |
| `web/common_js/exportFF.js` | `buildPrimaryAdj`, `validateSpringGraphDistances`, `buildTopologyFull` | `build_linearized_topology.mjs` |

### Build new (integration work packages)

| WP | Deliverable | Priority | Status |
|----|-------------|----------|--------|
| **WP-I1** | `ensemble.yaml` schema + param sampler | P0 | Done (JSON config) |
| **WP-I2** | `pyBall/nanocrystal_pipeline.py` — stage runner with cache | P0 | Done (JS orchestrator) |
| **WP-I3** | `modes_to_histogram()` + running accumulator | P0 | Done (Python) |
| **WP-I4** | NPZ cache I/O + `manifest.jsonl` | P0 | Done |
| **WP-I5** | `pyBall/nanocrystal_viz.py` — orthographic SVG | P1 | Done (JS `nanocrystalSvg.js`) |
| **WP-I6** | CLI: `scripts/run_nanocrystal_ensemble.py` | P0 | **Done** — replaced by unified `scripts/nanocrystals.mjs` CLI |
| **WP-I7** | LFF analytical Hessian assembler (Python, from topology NPZ) | P1 (throughput) | Pending |
| **WP-I8** | Ellipsoidal cut in `gen_nanocrystals.mjs` | P2 | Pending |

### Explicitly not building (v1)

- GPU vibration solver
- Per-frequency sparse probing
- Per-crystal HTML stick viewers in production
- Committing generated artifacts to git

## 8. Pipeline stages (per crystal)

```
1. SAMPLE    ensemble.yaml → gen_params
2. GENERATE  node nanocrystals.mjs generate … → parse to 01_init.npz (skip if cached)
3. RELAX     MMFF bonds+angles, no NB → 02_relaxed.npz
4. TOPOLOGY  node nanocrystals.mjs topology → 03_topology.npz (optional for Path A)
5. HESSIAN   Path A: getHessian3Nx3N → 04_hessian.npz
             Path B: assemble K from 03_topology.npz → 04_hessian.npz
6. DIAG      eigh → omegas_modes
7. SPECTRUM  histogram → 05_spectrum.npz; add to accumulator
8. DEBUG     if enabled: SVG + spectrum plot
```

**MMFF relax settings (match existing tests):** bonds + angles on, non-bonded off, `Fconv ~ 1e-3–1e-4 eV/Å`, ~1000–2000 steps. M-G5 (post-cap relax) is **required** for this pipeline — raw generated `.xyz` is not the vibration input.

**Cap H–H bonds:** use `--capHHBonds 1` in generation so MMFF/LFF does not need expensive NB H–H repulsion.

## 9. Testing milestones (integration phase)

| ID | Criterion |
|----|-----------|
| **M-P0** | Single G2 crystal: generate → relax → MMFF H → eigh → histogram; matches `test_vibration_spectra.py` peak positions ± broadening |
| **M-P1** | Cache round-trip: delete stage 5, resume from `04_hessian.npz` only |
| **M-P2** | Ensemble N=10: manifest complete, `spectrum_final.npz` = sum of per-crystal hists |
| **M-P3** | Debug mode: 3 SVGs + spectrum plot for 1 crystal |
| **M-P4** | Defect on/off: ensemble with `insertProb=0.2` runs without crash; defect counts in manifest |
| **M-P5** | Path B parity: LFF analytical ω_i vs MMFF FD on adamantane + G2 (±5% tolerance) |
| **M-P6** | Throughput: 100 crystals caps-only, report timing breakdown (gen / relax / hessian / eigh) |

## 10. Expected bottlenecks (where to optimize later)

| Stage | ~609-atom G2 | Notes |
|-------|--------------|-------|
| Generation | Fast (Node) | Not bottleneck |
| MMFF relax | Seconds | Acceptable |
| **MMFF FD Hessian** | **Slowest** | First optimization target → Path B (LFF analytical) |
| Dense `eigh` | &lt;0.1 s | **Not worth optimizing** per your observation |
| Histogram + accumulate | Negligible | — |
| Debug SVG | Per-crystal | Disable in production |

## 11. Suggested execution order

```
Phase 0 — Validate single-crystal MMFF path (M-P0)
    ↓
Phase 1 — Cache layer + ensemble runner, caps-only, N=10–50 (M-P1, M-P2)
    ↓
Phase 2 — Debug viz + defect sampling (M-P3, M-P4)
    ↓
Phase 3 — Shape mixture (disk/sphere/facet) + large ensemble N=500+
    ↓
Phase 4 — LFF analytical Hessian for throughput (M-P5), only if MMFF FD dominates wall time
```

**Parallel work from the three completed jobs:** no further parallel agent splits needed. This is a **fourth integration job** owning `pyBall/nanocrystal_pipeline.py`, `scripts/run_nanocrystal_ensemble.py`, and cache helpers — read-only on `gen_nanocrystals.mjs`, `MMFFLTopology.js`, `FTIR.py` (only add histogram helpers, do not rework sparse solver).

## 12. Open decisions (need your call before coding)

1. **Element:** Si-only first, or Si/C mixture in ensemble?
2. **Hessian v1:** MMFF FD only until M-P0 passes, or parallel Path B from the start?
3. **Spectrum grid:** fixed `ω ∈ [0, 3500] cm⁻¹`, `dω = 2 cm⁻¹`, `σ = 10 cm⁻¹` — or tie to experimental reference?
4. **IR weighting:** plain mode histogram (your v1) vs dipole-weighted (`vibration_spectrum_from_modes`) for v2?
5. **Artifact root:** `OUT_nanocrystal_ensemble/` (gitignored) vs configurable `--cache-dir`?

I can turn this into `doc/Topics/FTIR_Nanocrystals/Nanocrystal_spectrum_ensemble.progress.md` with agent contract and milestones when you want to start implementation. No code until you confirm the open decisions above (or tell me to use sensible defaults).

---

# USER

important consideration
we have part of it in javascript (generation of cryustasl, and javascript is faster than python, can we save .npz and generate .svg images direcly from there, so we do not wast time loading the crystal geoemtry in python? We can do all the steps (generation of crustal, generation of linearized topology) using javascrib based script, and only the eigensolver do in python (loading the .npz loaded from js)

about additiond geometry ensamble vadiation ideas
* we do not ned isotopic effects, log size distirubtion, dagnilmb bonds, or oxidized surface for now. But we should have broades variability of cutting planes, and curring cylinger and cutting casula, and w should think how to let user choose some variation ranges from which we then randomly toss these paramters how to cut the crystals before passivation

this must be very well documented

the stage time estimates/bottlenec you just estimated, you do not have solid suppoer for that, we should add benchmarking/timing measurements in the pipleline and report them

Decisions:

1) we start with Carbon (diamond) for now, Si letf for later when eveything is tested
2) what you ask here? if you should run the hessian ? (yes, how you can get spectrum without it), or if you should optimize this? (no we should finish pipeline first)
3) spectrum grid  (for histogram) is given by eighevalues, we use eigh which give use all eigenvalues, so we should consider all for them adjust grid accorignly (with some margin)
4) we atart with just plain historgram, but it is true that later we should do some weighting.

most rigorous way of such weighting is to consider some probe vectors (excitation vectros) and then our eigenvectror should be respons to it.

the proper equation is (I.e - H)v=b, where b is probe (impulse), e is eigenvalue, v is eigenvector, or something like that. What I;m saying we should not do just simple 

<v|b> to get the proper aimplitide of the response, but perhaps |<v|H|b>|^2 or somthing like that? Think it throu this is just branstroming, not rigorous.
Now to do this efficiently and quicly we can jsut multiply the matrixes? dot(dot(V,H),b) or something like this?

5) peopre data management is on you, but yet, we should have option to specify work folder or data folder/output folder

now do not program yet, write updated plan and rigorous analysis of the problem considering my imputs and decisions

---

# Cursosr Composer 2.5

I've written the updated plan to [`doc/Topics/FTIR_Nanocrystals/Nanocrystal_pipeline.progress.md`](doc/Topics/FTIR_Nanocrystals/Nanocrystal_pipeline.progress.md) (correct spelling; your open file was `prograss.md`).

## Architecture: JS front, Python back

**Node owns everything that touches the molecular graph:**

- Generation → `01_init.npz`
- Linearized topology on relaxed geometry → `03_topology.npz`
- Analytical LFF **K** assembly → `04_hessian.npz` (target)
- Orthographic **SVG** (init + relaxed) — no Python geometry load

**Python subprocesses (NPZ in/out only, no `.xyz` on the hot path):**

- `nc_relax.py` — MMFF relax → `02_relaxed.npz` (C++ MMFF is Python-only today)
- `nc_solve_spectrum.py` — `project_rigid_modes` + `eigh` + histogram → `05_spectrum.npz`

Orchestrator: `scripts/run_nanocrystal_ensemble.mjs` spawns Python, writes `manifest.jsonl`, respects `--work-dir` / `--data-dir` / `--output-dir`.

## Cut ensemble (documented variation ranges)

v1 adds **`cylinder`** cut; **casual** = weighted mixture in `ensemble.json`, not a new algorithm.

| Mode | Sampled parameters |
|------|-------------------|
| `planes` | `planeTemplates`, `planeSymC`, `planeCScale`, `planeCJitter`, `nx/ny/nz` |
| `sphere` | `sphereR`, `sphereNrep` |
| `cylinder` *(new)* | `cylAxis`, `cylR`, `cylHalfLen` + replication |
| `casual` | `cut_mixture[]` with weights + nested ranges |

Out of scope for now: isotopes, log-normal size, dangling bonds, oxidation.

Full schema contract → `ensemble_config.schema.md` (to write at implementation).

## Spectrum grid (your decision #3)

Per crystal, after `eigh`:

```
omega_min = min(ω_i) − margin
omega_max = max(ω_i) + margin
margin = margin_frac × span
```

Histogram on that span; ensemble accumulation **extends and rebins** when a new crystal exceeds bounds (or store raw `omegas_modes` and histogram once at the end — both documented).

## Probe weighting — rigorous analysis (your decision #4)

**Correct equation** (frequency domain):

\[
A(\omega)\, u = b, \quad A(\omega) = K - \omega^2 M + i\eta\omega M
\]

**b** is a **force probe**, not an eigenvalue. Modal expansion with \(K v_i = \omega_i^2 M v_i\):

\[
u(\omega) = \sum_i \frac{v_i^T b}{\omega_i^2 - \omega^2 - i\eta\omega}\, v_i
\]

**Intensity:** \(I(\omega) \propto \sum_i |v_i^T b|^2 / |\omega_i^2 - \omega^2 - i\eta\omega|^2\)

This is what `FTIR.vibration_spectrum_from_modes` already does (with the mass-weighted inner product).

**Not correct:** \(|v_i^T K b|^2\) or \(|v_i^T H b|^2\) — since \(v_i^T K b = \omega_i^2\, v_i^T M b\), that overweights high modes and is not the Green's function response.

**Efficient after one `eigh`:** one matmul `f_tilde = V^T @ (f * m_inv_sqrt)` per probe; no further **K** or **H** multiplies. v1 = unit weights (plain histogram); v2 = bin \(|f_{\tilde{i}}|^2\) at \(\omega_i\) instead of 1.

## Hessian (your decision #2)

Yes — spectrum needs **K**. v1: finish pipeline with MMFF FD in Python if JS analytical K is not ready; do not optimize until **M-P6** timing report exists.

## Timing (no unsupported claims)

Every stage logs `timing_ms` to `manifest.jsonl`. End of run → `timing_report.json` with min/median/max per stage. **M-P6** requires a 100-crystal C run before any bottleneck decisions.

## Milestones summary

| Phase | Goal |
|-------|------|
| 0 | NPZ plumbing, single C crystal end-to-end |
| 1 | Orchestrator + cache + N=10 ensemble |
| 2 | `cut_mixture` + cylinder + schema docs |
| 3 | JS analytical K + debug SVG |
| 4 | Large ensemble + timing report |
| 5 | v2 probe weighting |
| 6 | Si (swap CIF/element only) |

No code until you want to start Phase 0. I can also merge/rename the typo `prograss.md` file if you want a single canonical doc.

---

# USER

Good, now implement it you can code now. And thest it for few crystal, I want to see the plots, svg structures and spectra, and the roganization of the output dirs

notice I expect that slowest would be optimization, possibly at some point we want to do it on GPU using projective dynamics or compute shaders, we hav it maybe already implemented in javacript, check it out in 
/home/prokop/git/FireCore/web/molgui_webgpu
 but we do not want to change out backed now, just make notes

then we will implemen this, fist using python MMFF backedn, only then we can think about moving it to js completely

---

# Nanocrystal Spectrum Ensemble Pipeline — Design & Progress

**Scope:** End-to-end pipeline from **stochastic generation** of C (diamond) nanocrystal ensembles → relaxation → Hessian → dense `eigh` → **histogram spectra** summed over many crystals.

**Related:** `Nanocrystal_generation.progress.md`, `Linearized_topology.progress.md`, `Sparse_vibration_solver.progress.md`, `Vibration_benchmark_npz.guide.md`.

**Status:** **Pass 1 implemented** (2026-06-16) — 3-crystal C sphere test green in `OUT_nc_ensemble_v2/`.
**Consolidation:** **Complete** (2026-06-22) — 5 scripts merged into unified `scripts/nanocrystals.mjs` CLI; shared functions moved to existing modules. Old scripts deprecated, not deleted.

**Canonical doc:** `Nanocrystal_pipeline.progress.md` (not `Nanocrystal_pipeline.prograss.md`).

**Strategic decisions (locked):**

| # | Decision |
|---|----------|
| 1 | **Carbon (diamond) first**; Si deferred until C pipeline validated |
| 2 | **Hessian required** for spectrum; **do not optimize** Hessian/eigensolver until full pipeline works |
| 3 | **Histogram grid** derived from eigenvalues per crystal (with margin), not a fixed global grid |
| 4 | **v1:** plain mode-count histogram; **v2:** probe-weighted response (see §6) |
| 5 | **Configurable** `--work-dir` / `--data-dir` / `--output-dir` (see §8) |
| 6 | **JS-first** for generation, topology, NPZ I/O, SVG; **Python-only** for MMFF relax + `eigh` + spectrum (see §2) |
| 7 | **No** isotope disorder, log-normal size, dangling bonds, oxidized surfaces in v1 |
| 8 | **Timing:** measured in pipeline, not guessed — no performance claims without data (see §9) |

## 1. Language split: why JS owns the front half

JavaScript (Node headless) is already faster for molecular graph operations and is where generation + MMFFL topology live. **Avoid round-tripping geometry through Python text parsers** (`.xyz`/`.mol2`).

| Stage | Runtime | Rationale |
|-------|---------|-----------|
| Ensemble param sampling | **JS** | Same process as generator; no IPC |
| Crystal generation | **JS** `nanocrystals.mjs generate` | Native; Miller/sphere cuts, caps, bridges |
| Write `01_init.npz` | **JS** `npzWrite.js` | Already used by `LinearizedTopologyNpz.js` |
| Debug SVG (init geometry) | **JS** (`nanocrystalSvg.js`) | Orthographic projection without Python load |
| MMFF relaxation | **Python subprocess** (`nc_relax.py`) | **v1 only** — C++ MMFF via `pyBall.MMFF`; GPU/JS relax deferred (§15) |
| Write `02_relaxed.npz` | **Python** | Output of relax stage only |
| Linearized topology (K₁₂/K₁₃/K₁₄) | **JS** `nanocrystals.mjs topology` | On **relaxed** geometry; writes `03_topology.npz` |
| Analytical LFF Hessian **K** assembly | **JS** (new) or Python v1 | Target: JS exports dense **K** in `04_hessian.npz`; v1 may use MMFF FD in Python |
| Debug SVG (relaxed) | **JS** | Reads `02_relaxed.npz` directly |
| Rigid-mode projection + `eigh` | **Python** | `numpy.linalg.eigh`; loads NPZ only |
| Histogram + ensemble accumulate | **Python** | Loads `04_hessian.npz` or mode list |
| Debug SVG (spectrum) | **JS** or Python | Minimal line plot; JS preferred for consistency |
| Timing / manifest | **Both** | JS stages log to shared `manifest.jsonl` |

**Python never parses `.xyz` for production throughput.** Subprocess boundaries exchange **NPZ + JSON meta only**.

```mermaid
flowchart LR
    subgraph js [Node.js orchestrator]
        CFG[ensemble.json] --> GEN[nanocrystals.mjs generate]
        GEN --> NPZ1[01_init.npz]
        GEN --> SVG1[init.svg]
        NPZ1 --> TOPO[nanocrystals.mjs topology]
        TOPO --> NPZ3[03_topology.npz]
        NPZ2 --> TOPO
        NPZ3 --> KJS[assemble K from springs]
        KJS --> NPZ4[04_hessian.npz]
        NPZ2 --> SVG2[relaxed.svg]
    end

    subgraph py [Python subprocesses]
        NPZ1 --> RELAX[MMFF relax]
        RELAX --> NPZ2[02_relaxed.npz]
        NPZ4 --> EIGH[eigh + histogram]
        EIGH --> NPZ5[05_spectrum.npz]
    end

    NPZ5 --> ACC[ensemble accumulator]
```

**Orchestrator:** `tests/tSiNCs/nanocrystals.mjs ensemble` — loops crystals, calls Python via `spawn` for relax and solve steps, owns cache invalidation and timing. Replaces deprecated `tests/tSiNCs/run_nanocrystal_ensemble.mjs`.

**Python entry points (thin):**

| Script | Input | Output |
|--------|-------|--------|
| `scripts/nc_relax.py` | `01_init.npz` | `02_relaxed.npz` |
| `scripts/nc_solve_spectrum.py` | `04_hessian.npz` | `05_spectrum.npz` |

Optional v1 fallback: `nc_build_hessian_mmff.py` if analytical K in JS is not ready — still NPZ in/out.

## 2. Pipeline stages (per crystal)

```
STAGE 0  SAMPLE     ensemble.json → gen_params (JSON blob, stored in meta)
STAGE 1  GENERATE   JS nanocrystals.mjs generate → 01_init.npz + optional init.svg
STAGE 2  RELAX      Python MMFF (bonds+angles, NB off) → 02_relaxed.npz
STAGE 3  TOPOLOGY   JS nanocrystals.mjs topology on relaxed coords → 03_topology.npz
STAGE 4  HESSIAN    JS analytical K from 03 + pos → 04_hessian.npz
                     (v1 fallback: Python MMFF FD → 04_hessian.npz)
STAGE 5  DIAG+SPEC  Python: project_rigid_modes, eigh, histogram → 05_spectrum.npz
STAGE 6  ACCUM      Add to running ensemble sum (05 → accumulated/)
STAGE 7  DEBUG      Optional SVG/plots (init, relaxed, spectrum) — off in production
```

**Element (v1):** C diamond from `diamond_primitive.cif` (`--applySymmetry 0`), H caps, `--capHHBonds 1` recommended.

**Passivation before relax:** always capped; bridge defects opt-in via `insertProb` / `collapseProb`.

## 3. Cut-shape ensemble — user-controlled variation ranges

All cut parameters are **sampled from ranges** declared in `ensemble.json`. The generator already supports `planes` and `sphere`; v1 adds **`cylinder`** and documents **casual** (mixed/random) sampling.

### 3.1 Cut modes

| Mode | `cutMode` | Parameters (sampled from ranges) | Shape intent |
|------|-----------|----------------------------------|--------------|
| **Miller facets** | `planes` | `nx, ny, nz`, `planeTemplates[]`, `planeSymC`, `planeCScale`, `planeC0`, `planeCJitter`, `planeMode` | Faceted NC; multi-plane = broad facet variability |
| **Sphere** | `sphere` | `sphereR`, `sphereNrep` | Compact, roughly isotropic |
| **Cylinder** | `cylinder` *(new)* | `cylAxis` (x\|y\|z), `cylR`, `cylHalfLen`, `sphereNrep` or `nx,ny,nz` | Tube / rod |
| **Casual** | `casual` *(meta)* | Weighted mixture of above + random `planeTemplates` subset | User declares mixture weights |

**Cylinder definition (to implement):** replicate cell, center origin, keep atom if `(r_⊥ ≤ R) ∧ (|z_axis| ≤ L/2)` in axis-aligned frame (optional random rotation of axis in v2).

**Broad plane variability (existing knobs, document ranges):**

| Parameter | Effect | Example range |
|-----------|--------|---------------|
| `planeTemplates` | Miller families | `["a111"]`, `["a100"]`, `["a110"]`, `["a111","a100"]`, `["a111","a100","a110"]` |
| `planeSymC` | Symmetry expansion count | `[1, 12]` |
| `planeCScale` | Cut depth vs lattice | `[0.25, 0.65]` |
| `planeCJitter` | Random offset fraction | `[0, 0.15]` |
| `planeC0` | Base offset | `[0, 0.2]` |
| `nx, ny, nz` | Supercell extent | independent `[1, 5]` each |

**Casual sampling:** not a separate algorithm — a **mixture model** in config:

```json
"cut_mixture": [
  { "weight": 0.25, "cutMode": "sphere",   "sphereR":   { "uniform": [4.0, 8.0] } },
  { "weight": 0.25, "cutMode": "cylinder", "cylR":      { "uniform": [3.0, 6.0] },
                                      "cylHalfLen": { "uniform": [4.0, 10.0] } },
  { "weight": 0.50, "cutMode": "planes",   "planeTemplates": { "choice": [
      ["a111"], ["a100"], ["a110"], ["a111","a100"], ["a111","a100","a110"]
    ]}, "planeCScale": { "uniform": [0.30, 0.55] }, "planeCJitter": { "uniform": [0.0, 0.12] } }
]
```

Each crystal: draw `u ~ Uniform(0,1)`, pick component by cumulative weights, then sample all nested ranges.

### 3.2 Defects (v1)

| Parameter | Range | Notes |
|-----------|-------|-------|
| `insertProb` | `[0, 0.2]` | SiH₂-like bridge insertion (C analogue in generator) |
| `collapseProb` | `[0, 0.15]` | Bridge collapse |

Default validation subset: both **0** (caps-only).

### 3.3 `ensemble.json` schema (documented contract)

Top-level keys:

| Key | Type | Description |
|-----|------|-------------|
| `schema_version` | int | Currently `1` |
| `element` | str | `"C"` for v1 |
| `n_crystals` | int | Ensemble size |
| `seed_base` | int | Master seed; crystal `i` uses `seed_base + i` |
| `cif` | str | Path to unit cell |
| `replication` | object | Default `nx/ny/nz` ranges (used when cut mode needs them) |
| `cut_mixture` | array | Weighted cut specs (§3.1) |
| `passivation` | object | `caps`, `capHHBonds`, `resolveClashes`, `outwardBias` |
| `defects` | object | `insertProb`, `collapseProb` ranges |
| `spectrum` | object | `margin_frac`, `min_bins`, `sigma` (Gaussian blur, optional) |
| `probes` | object | v2: probe definitions (§6) |
| `debug` | object | `enabled`, `max_crystals`, `svg_rotation` |

Full schema with examples → **`doc/Topics/FTIR_Nanocrystals/ensemble_config.schema.md`** (to write at implementation).

## 4. NPZ cache & folder layout

User-configurable roots:

| CLI flag | Purpose |
|----------|---------|
| `--work-dir` | Per-run scratch, subprocess temp |
| `--data-dir` | Persistent crystal cache (resume) |
| `--output-dir` | Final ensemble spectra, timing reports |

Default: `--data-dir OUT_nc_ensemble/data`, `--output-dir OUT_nc_ensemble/out`.

```
{data-dir}/
  ensemble_meta.json           # frozen copy of ensemble.json + git hash
  manifest.jsonl               # one line per crystal per stage completion
  crystals/{crystal_id}/
    meta.json                  # gen_params + stage checksums
    01_init.npz
    02_relaxed.npz
    03_topology.npz
    04_hessian.npz
    05_spectrum.npz
    debug/                     # if debug.enabled
      init.svg, relaxed.svg, spectrum.svg

{output-dir}/
  accumulated/
    spectrum_running.npz       # updated every N crystals
    spectrum_final.npz
  timing_report.json           # measured stats (§9)
  timing_report.md             # human-readable summary
```

### 4.1 NPZ schemas

**`01_init.npz`** (JS writes)

| Key | dtype | shape |
|-----|-------|-------|
| `pos` | float64 | (N, 3) |
| `Z` | int32 | (N,) |
| `neigh_idx`, `neigh_count`, `bonds_ij` | int32 | optional bond list from mol |
| `gen_params` | UTF-8 str | JSON |

**`02_relaxed.npz`** (Python writes)

| Key | Notes |
|-----|-------|
| `pos` | Relaxed coordinates |
| `Z` | Unchanged |
| `fmax`, `n_steps`, `converged` | Relax diagnostics |
| `timing_ms` | Stage wall time |

**`03_topology.npz`** (JS writes — same as `Linearized_topology.progress.md`)

`pos`, `Z`, `neigh_idx`, `KLs`, `stick_class`, `neigh_count`, `max_neighbors`, counts.

**`04_hessian.npz`**

| Key | Notes |
|-----|-------|
| `K` | Dense (3N×3N) stiffness |
| `M` | Diagonal mass (3N×3N diagonal or (N,) atom masses) |
| `pos` | Geometry at which K built |
| `Z` | |
| `source` | `"lff_analytical"` or `"mmff_fd"` |
| `timing_ms` | |

**`05_spectrum.npz`**

| Key | Notes |
|-----|-------|
| `omegas_modes` | All eigenfrequencies from `eigh` |
| `omega_centers`, `hist` | Binned histogram |
| `grid_meta` | JSON: `omega_min`, `omega_max`, `margin_frac`, `n_bins` |
| `timing_ms` | |

### 4.2 Cache / resume rules

| Rule | Behavior |
|------|----------|
| Crystal ID | `{seed:08d}_{index:06d}` |
| Skip stage | If output NPZ exists and `meta.json` `gen_params_hash` matches |
| Invalidate downstream | Changing `gen_params` invalidates stages 2–5 |
| Text export | Optional `--export-xyz 1` for human inspection only; not on hot path |

## 5. Spectrum v1: plain histogram from eigenvalues

### 5.1 Per-crystal algorithm (Python)

1. Load `K`, `M`, `pos` from `04_hessian.npz`.
2. `K_proj = project_rigid_modes(K, M, pos)` (same as `FTIR.project_rigid_modes`).
3. Mass-weight: `K_mw = M^{-1/2} K_proj M^{-1/2}`.
4. `w, V = np.linalg.eigh(K_mw)` → `omega_i = sqrt(max(w_i, 0))`.
5. **Adaptive grid:**
   - `omega_min = max(0, min(omega_i) - margin)`
   - `omega_max = max(omega_i) + margin`
   - `margin = margin_frac * (omega_max - omega_min)` (config, default e.g. 0.02)
   - `n_bins = max(min_bins, ceil((omega_max - omega_min) / domega_target))` OR one bin per mode cluster — document choice: **uniform bin width from span / min_bins**
6. Histogram: count modes per bin (each mode = δ(ω − ω_i) → bin index).
7. Optional Gaussian blur: convolve bin counts with `σ` (instrument broadening).

### 5.2 Ensemble accumulation

When crystal k has grid `[omega_centers_k, hist_k]`:

- Maintain running global `[omega_min, omega_max]`; **extend** when new crystal exceeds bounds (rebin accumulated sum onto new grid — document rebinning as linear interpolation of bin centers).
- `hist_total += hist_k_rebinned`.

Alternative (simpler v1): store **raw `omegas_modes` only** per crystal in cache; histogram **once** when finalizing ensemble (memory cost O(N_crystals × 3N) — acceptable for moderate ensembles).

**Do not use fixed `[0, 3500]` grid** — grid follows data.

## 6. Spectrum v2: probe-weighted response (rigorous analysis)

### 6.1 Correct frequency-domain equation

Undamped harmonic system with driving force:

\[
M\ddot{u} + K u = f e^{i\omega t}
\]

In frequency domain:

\[
A(\omega)\, u(\omega) = b, \qquad A(\omega) = K - \omega^2 M + i\eta\omega M
\]

Here **b** is the **force probe** (generalized force vector), not an eigenvalue λ. The user’s `(Iλ - H)v = b` mixes notation — the standard form uses **K** (stiffness), not Hamiltonian H, and **ω²** enters through the dynamical matrix, not a single eigenvalue λ per solve.

### 6.2 Modal expansion (correct weighting)

Generalized eigenproblem (physical modes):

\[
K v_i = \omega_i^2 M v_i
\]

Assume **M-orthonormal** modes: `v_i^T M v_j = δ_ij`.

Force projection onto mode i:

\[
f_i = v_i^T b
\]

Displacement response:

\[
u(\omega) = \sum_i \frac{f_i}{\omega_i^2 - \omega^2 - i\eta\omega}\, v_i
\]

**Observable intensity** (displacement / mechanical):

\[
I_{\mathrm{mech}}(\omega) = \|u(\omega)\|^2 = \sum_i \frac{|f_i|^2}{|\omega_i^2 - \omega^2 - i\eta\omega|^2}
\]

(+ cross terms if computing vector response; intensity uses squared magnitude of scalar projection or dipole — see below.)

This is **exactly** what `FTIR.vibration_spectrum_from_modes` implements (mass-weighted variant: absorb `M^{1/2}` into `V` and `f`).

### 6.3 Why NOT `|<v|K|b>|²` or `|<v|H|b>|²`

Using `v_i^T K b` as coupling is **wrong for linear response**:

\[
v_i^T K b = \omega_i^2\, (v_i^T M b) \neq v_i^T b
\]

This injects an extra factor ω_i² — **overweights high-frequency modes** and is not the mechanical Green’s function. Similarly `|v_i^T H b|²` has no standard interpretation in this problem unless H is a specifically defined operator (e.g. dipole moment surface).

**Correct coupling for IR:** replace scalar `f_i = v_i^T b` with **dipole derivative**:

\[
\mu(\omega) = \sum_i \frac{(\partial\mu/\partial Q_i)\, q_i(\omega)}{...}
\]

or in linear approximation: effective probe `b_μ` such that induced dipole `μ(ω) = Γ u(ω)`; intensity `|μ|²`.

For **symmetric stretch** probes (existing in `test_vibration_spectra.py`): `b_n = q_n * ê` at atom n — still enters as **`v_i^T b`**, not `v_i^T K b`.

### 6.4 Efficient computation after `eigh` (no per-ω solves)

Given `V = [v_1 | … | v_N]` and force vector `b`:

```python
# Mass-weighted convention (matches FTIR.vibration_spectrum_from_modes):
m_inv_sqrt = 1.0 / np.repeat(np.sqrt(m_atom), 3)
f_mw = f * m_inv_sqrt
f_tilde = V.T @ f_mw          # shape (N_modes,)  — ONE matmul O(N²)
# For each omega (or histogram of mode contributions):
# weight_i(omega) = |f_tilde_i|² / |omega_i² - omega² - i*eta*omega|²
```

**No matrix H or K needed after diagonalization** for mode-sum spectrum. One `V.T @ f_mw` per probe direction.

Multiple probes (x, y, z symmetric): three matmuls or store `F_tilde = V.T @ (F * m_inv_sqrt)` with F shape (N_modes, n_probes).

**v2 histogram variant:** instead of counting each mode equally, bin **|f_tilde_i|²** (or sum over probes) at ω = ω_i — weighted δ-comb, then blur.

### 6.5 Relation to v1 plain histogram

| v1 | v2 |
|----|-----|
| Each mode weight 1 | Each mode weight `\|f_tilde_i\|^2` (or `\|μ_i\|^2`) |
| Density of states-like | IR / mechanical response to specified probes |
| No probe config | `ensemble.json` → `probes[]` |

## 7. JS components to add

| Module | Role |
|--------|------|
| `scripts/run_nanocrystal_ensemble.mjs` | Main orchestrator, cache, timing, Python spawn |
| `web/molgui_webgpu/NanocrystalNpz.js` | `writeCrystalNpz` / `readCrystalNpz` (init + relaxed coords) |
| `web/molgui_webgpu/NanocrystalSvg.js` | Orthographic SVG: `pos, Z, bonds, R_view[3×3]` |
| `web/molgui_webgpu/LFFHessian.js` | Assemble dense K from `neigh_idx`, `KLs`, `pos` (3×3 blocks per spring) |
| `scripts/nanocrystals.mjs generate` | Add `cutMode: cylinder` |
| `scripts/nc_relax.py` | MMFF relax NPZ → NPZ |
| `scripts/nc_solve_spectrum.py` | `eigh` + histogram NPZ → NPZ |

**Reuse unchanged:** `npzWrite.js`, `LinearizedTopologyNpz.js`, `exportFF.js` (topology), `Nanocrystals.js` (generation, extended).

## 8. What we reuse vs build

| Existing | Use in pipeline |
|----------|-----------------|
| `nanocrystals.mjs generate` | Generation + caps + defects |
| `nanocrystals.mjs topology` | Stage 3 topology NPZ |
| `LinearizedTopologyNpz.js` / `npzWrite.js` | All JS NPZ writes |
| `MMFF.run`, `getHessian3Nx3N` | Relax; optional v1 FD Hessian |
| `FTIR.project_rigid_modes` | Python solve stage |
| `FTIR.vibration_spectrum_from_modes` | Reference for v2 probe weighting |
| `bootstrap_vibration_parallel_fixtures.py` | Pattern for NPZ export |

| **Deferred** | Reason |
|--------------|--------|
| GPU / sparse frequency scan | Not needed — dense `eigh` sufficient for now |
| Si element | After C validated |
| Isotope, oxidation, dangling bonds | Out of v1 scope |

## 9. Timing & benchmarking (mandatory)

**No stage timing claims in docs without measured data.** Every stage records wall time:

```json
{ "crystal_id": "...", "stage": "relax", "timing_ms": 1234.5, "natoms": 270, "status": "ok" }
```

Appended to `manifest.jsonl`. At ensemble end:

**`timing_report.json`:**

```json
{
  "n_crystals": 100,
  "stages": {
    "generate":   { "min_ms": 0, "median_ms": 0, "max_ms": 0, "total_ms": 0 },
    "relax":        { ... },
    "topology":     { ... },
    "hessian":      { ... },
    "eigh":         { ... },
    "histogram":    { ... },
    "svg_debug":    { ... }
  },
  "notes": "Auto-generated; do not cite as hardware-independent benchmarks."
}
```

Orchestrator prints one-line summary per crystal and final table. **M-P6** milestone requires this report from a 100-crystal caps-only C run before any optimization decisions.

## 10. Debug vs production

| | Production | Debug |
|---|------------|-------|
| SVG geometry | off | init + relaxed (same `R_view`) |
| SVG spectrum | off | per-crystal line histogram |
| Cache per-crystal `05_spectrum.npz` | optional | on |
| HTML stick viewer | off | manual, not in loop |

Debug controlled by `ensemble.json` → `debug.enabled`, `debug.max_crystals`.

**SVG (JS):** orthographic: `p_view = (R_view @ pos.T).T` → use (x, y) for 2D; circles by covalent radius / Z color.

## 11. Milestones

| ID | Criterion |
|----|-----------|
| **M-P0** | Single C NC: JS gen → NPZ → Python relax → JS topology → Python `eigh` → histogram; no `.xyz` on hot path |
| **M-P1** | Resume: skip completed stages via cache |
| **M-P2** | Ensemble N=10: accumulated spectrum = sum of per-crystal hists (rebin documented) |
| **M-P3** | Debug: init/relaxed SVG + spectrum plot for 1 crystal |
| **M-P4** | `cut_mixture` with planes + sphere + cylinder samples without crash |
| **M-P5** | LFF analytical K (JS) vs MMFF FD K (Python): ω_i agree within tolerance on adamantane + one NC |
| **M-P6** | N=100 caps-only C: **`timing_report.json` published** |
| **M-P7** | v2: one probe direction; weighted hist matches `vibration_spectrum_from_modes` peak locations |

## 12. Execution order

```
Phase 0 — NPZ plumbing + single-crystal M-P0 (MMFF FD Hessian OK in Python)
Phase 1 — JS orchestrator + ensemble.json schema + cache (M-P1, M-P2)
Phase 2 — cut_mixture + cylinder cut + docs (M-P4)
Phase 3 — JS analytical K + debug SVG (M-P3, M-P5)
Phase 4 — Large ensemble + timing report (M-P6)
Phase 5 — Probe weighting v2 (M-P7)
Phase 6 — Si port (change element + CIF only, after C stable)
```

## 13. Agent contract (integration job)

| May edit | Must not edit |
|----------|---------------|
| `scripts/run_nanocrystal_ensemble.mjs` (new) | `LFF.cl`, GPU kernels |
| `web/molgui_webgpu/NanocrystalNpz.js`, `NanocrystalSvg.js`, `LFFHessian.js` (new) | Core cap logic unless bugfix |
| `scripts/gen_nanocrystals.mjs` (cylinder cut only) | |
| `scripts/nc_relax.py`, `scripts/nc_solve_spectrum.py` (new) | |
| `doc/Topics/FTIR_Nanocrystals/ensemble_config.schema.md` (new) | |
| `pyBall/FTIR.py` (thin histogram helpers only) | Sparse solver rework |

**Outputs:** gitignore `{data-dir}`, `{output-dir}`; do not commit generated crystals.

## 14. Implementation pass 1 — what to build and test (next coding session)

**Goal:** Few-crystal **debug run** with visible artifacts: directory layout, init/relaxed **SVG**, per-crystal **spectrum SVG/PNG**, accumulated ensemble spectrum, and **`timing_report.json`** with **measured** stage times (no guessed bottlenecks).

### 14.1 Scope for pass 1 (minimal viable pipeline)

| In scope | Out of scope (pass 1) |
|----------|------------------------|
| C diamond, `cutMode: sphere` + one `planes` preset | `cylinder` cut, full `cut_mixture` |
| `n_crystals = 3–5`, `debug.enabled = true` | N=100 production run |
| JS: gen → `01_init.npz` + `init.svg` | JS analytical K assembly |
| Python: MMFF relax → `02_relaxed.npz` | GPU relax |
| JS: topology → `03_topology.npz` | |
| Python: MMFF FD Hessian → `04_hessian.npz` | LFF analytical K |
| Python: `eigh` + histogram → `05_spectrum.npz` + `spectrum.svg` | Probe weighting v2 |
| JS: `relaxed.svg` from `02_relaxed.npz` | |
| Accumulated `spectrum_final.npz` + `ensemble_spectrum.svg` in `--output-dir` | |

**Hot path:** no `.xyz` parse in Python for generation/topology. Python may write optional `*.xyz` sidecar when `--export-xyz 1`.

### 14.2 Files to create (pass 1)

| File | Role |
|------|------|
| `scripts/run_nanocrystal_ensemble.mjs` | Orchestrator: sample params, spawn Python, cache, manifest, timing |
| `scripts/ensemble.example.json` | Minimal C ensemble config (3 crystals, debug on) |
| `web/molgui_webgpu/NanocrystalNpz.js` | Read/write crystal NPZ (`pos`, `Z`, `bonds_ij`, `gen_params`) |
| `web/molgui_webgpu/NanocrystalSvg.js` | Orthographic SVG from NPZ |
| `scripts/nc_relax.py` | Load `01_init.npz` → MMFF relax → `02_relaxed.npz` |
| `scripts/nc_build_hessian_mmff.py` | Load `02_relaxed.npz` → `getHessian3Nx3N` → `04_hessian.npz` |
| `scripts/nc_solve_spectrum.py` | Load `04_hessian.npz` → `eigh` + histogram → `05_spectrum.npz` |
| `scripts/nc_accumulate_spectrum.py` | Merge per-crystal hists → `spectrum_final.npz` (or inline in orchestrator) |
| `doc/Topics/FTIR_Nanocrystals/ensemble_config.schema.md` | Documented JSON schema |

**Extend (thin):** `build_linearized_topology.mjs` — accept NPZ path instead of mol2 when `--input-npz` set.

### 14.3 Expected output directory layout (debug run)

```bash
node scripts/run_nanocrystal_ensemble.mjs \
  --config scripts/ensemble.example.json \
  --data-dir OUT_nc_ensemble/data \
  --output-dir OUT_nc_ensemble/out \
  --work-dir OUT_nc_ensemble/work
```

```
OUT_nc_ensemble/
  data/
    ensemble_meta.json
    manifest.jsonl
    crystals/
      00000042_000000/
        meta.json
        01_init.npz
        02_relaxed.npz
        03_topology.npz
        04_hessian.npz
        05_spectrum.npz
        debug/
          init.svg
          relaxed.svg
          spectrum.svg
      00000043_000001/
        ...
  out/
    accumulated/
      spectrum_running.npz
      spectrum_final.npz
    plots/
      ensemble_spectrum.svg      # sum over crystals
      ensemble_spectrum.png      # optional matplotlib if Python writes
    timing_report.json
    timing_report.md
  work/                          # subprocess temp (gitignored)
```

User inspects: `debug/*.svg` per crystal, `out/plots/ensemble_spectrum.svg`, `timing_report.md`.

### 14.4 MMFF relax settings (pass 1 — match existing NC tests)

From `test_nanocrystal_sparse_hessian.py`:

```python
MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, Dihedrals=+1, SurfAtoms=-1, GridFF=-1, ...)
MMFF.run(nstepMax=2000, dt=0.05, Fconv=1e-4, damping=0.1, ialg=2, omp=False)
```

**Init from NPZ:** reconstruct MMFF state from `pos`, `Z`, and bond list in `01_init.npz` (or export minimal mol2 to temp only if MMFF.init requires file — prefer in-memory if API allows; if not, write `{work-dir}/_tmp_{id}.mol2` and delete after relax).

**Fail loudly:** if `fmax > Fconv` after `nstepMax`, record `converged: false` in NPZ and manifest; do not silently continue unless `--allow-unconverged 1`.

### 14.5 First test command (acceptance)

```bash
# From repo root; requires tests/tMMFF data symlink + MMFF .so built
node scripts/run_nanocrystal_ensemble.mjs \
  --config scripts/ensemble.example.json \
  --data-dir OUT_nc_ensemble/data \
  --output-dir OUT_nc_ensemble/out \
  --n-crystals 3
```

**Pass when:**

1. All 3 crystals have complete stage NPZs under `data/crystals/`.
2. Each has `debug/init.svg`, `debug/relaxed.svg`, `debug/spectrum.svg`.
3. `out/plots/ensemble_spectrum.svg` exists.
4. `out/timing_report.json` has non-zero `timing_ms` for every stage that ran.
5. `manifest.jsonl` has 3×(stages) lines with `status: ok` (or documented failures).

## 15. Relaxation bottleneck & future GPU path (notes only — do not wire in pass 1)

### 15.1 Hypothesis

**Relaxation is expected to dominate wall time** for capped nanocrystals (many bonds/angles, hundreds–thousands of MMFF FD force evaluations). **Diagonalization (`eigh`) and histogram are expected to be fast** at target sizes — but pass 1 **must measure** this via `timing_report.json`; no doc claims without data.

### 15.2 v1 policy: Python MMFF only

| Rule | Rationale |
|------|-----------|
| Pass 1 uses **`pyBall.MMFF`** (C++ backend) for relax + FD Hessian | Same physics as validated `test_vibration_spectra.py` / `test_nanocrystal_sparse_hessian.py` |
| **Do not** call WebGPU / XPDB / PD from orchestrator in pass 1 | Different force field (LFF springs vs full MMFF angles); parity not proven on NC caps |
| **Do not** change C++ MMFF or `LFF.cl` for this job | Integration only |

After pass 1 timing report: if `relax` ≫ `hessian` ≫ `eigh`, optimize relax first — but only after pipeline is end-to-end green.

### 15.3 Existing GPU / JS relaxation inventory (`web/molgui_webgpu`)

| Component | Path | What it does | Headless? | NC pipeline fit |
|-----------|------|--------------|-----------|-----------------|
| **XPDB WebGPU** | `XPDB_WebGPU.js`, `xpdb.wgsl`, `run_xpdb_webgpu_headless.mjs` | MMFFL linear springs + optional collisions; Jacobi on GPU | **Yes** — `headless_init.mjs` + WebGPU in Node | Topology from `buildXPDBInputsFromMol` / MMFFL; **not full MMFF angles** |
| **XPDB CPU** | `main.js` `xpdbRelaxSteps`, `XPDBTopology.js` | Same springs, CPU Jacobi | Browser / molgui | Same caveat |
| **PD Jacobi GPU** | `ScriptRunner.js` `relaxJacobiGPU` → `pdSimulation` | Projective dynamics bond constraints on GPU | **Browser** (needs `MolGUIApp`, WebGPU canvas) | Bond-only; angles via virtual bonds in PD |
| **PD Jacobi CPU** | `ScriptRunner.js` `relaxJacobiCPU` → `runCpuJacobi` | CPU reference for PD | Browser | Debug / parity |
| **PD core** | `web/molgui_web/js/ProjectiveDynamics.js` `PDSimulation` | Constraint Jacobi, `setTopology(mol, mmParams)` | Imported by molgui_webgpu | Shared math with XPDB family |
| **pyOpenCL LFF** | `pyBall/OCL/LFFSolver.py`, `LFF.cl` | GPU spring dynamics (`lff_jacobi`) | Python + OpenCL | Relax MD, not static minimization; `MAX_NEIGHBORS=8` too small for full K₁₂+K₁₃+K₁₄ |
| **MMFFL topology** | `MMFFLTopology.js`, `build_linearized_topology.mjs` | K₁₂/K₁₃/K₁₄ packing for GPU | Node headless | Already in pipeline stage 3 |

**Headless WebGPU entry point:** `run_xpdb_webgpu_headless.mjs` — loads XYZ, builds MMFFL topology, runs `XPDB_WebGPU` steps; parity checked vs `test_XPDB_new_dump_headless.py` on pentacene.

### 15.4 Future relax acceleration options (ordered, after pass 1 timing)

| Priority | Approach | Pros | Cons / gate |
|----------|----------|------|-------------|
| **F1** | Tune Python MMFF (`nstepMax`, `Fconv`, `omp=True`) | Zero new code | Still CPU C++; may be enough for medium N |
| **F2** | **XPDB WebGPU headless** with MMFFL topology from stage 3 | Already headless; fast Jacobi; same springs as vibration LFF path | Relaxed geometry ≠ MMFF-minimized; **must validate** ω spectrum shift vs MMFF relax on 2–3 crystals |
| **F3** | Port `relaxJacobiGPU` / `PDSimulation` to Node headless (like XPDB) | Reuses PD infrastructure | Browser-coupled today; bond-topology setup from NPZ needed |
| **F4** | pyOpenCL `LFFSolver` energy minimization mode | Matches LFF Hessian path | `MAX_NEIGHBORS=8`; need `vib_jacobi`/48-neighbor fork; OpenCL device |
| **F5** | Full MMFF on GPU | Exact physics | **Not implemented** — C++ MMFF is CPU |

**Recommended future path (when timing demands it):**

1. Keep **MMFF Python relax** as **gold reference**.
2. Add optional `--relax-backend mmff|xpdb_webgpu` only after **F2 parity**: compare `02_relaxed.npz` RMSD and resulting `omegas_modes` vs MMFF on same `01_init.npz`.
3. Do **not** replace MMFF relax until parity milestone **M-P8** (new): `max_rel_pos_RMSD < ε` and `median |Δω|/ω < δ` on test set.

### 15.5 Why not GPU relax in pass 1

- Vibration pipeline validation used **MMFF-relaxed** geometries; switching relax changes modes unpredictably.
- XPDB/PD use **linearized springs** (or PD constraints), not MMFF angle/dihedral terms.
- WebGPU PD path is **not** headless-batch-ready without molgui app bootstrap (unlike `run_xpdb_webgpu_headless.mjs`).
- User directive: **finish Python MMFF pipeline first**, then consider JS/GPU relax.

## 16. Updated execution order (implementation)

```
Pass 1 (NOW) — M-P0 + M-P3 subset:
  NPZ plumbing, nc_relax.py (MMFF), nc_build_hessian_mmff.py, nc_solve_spectrum.py
  orchestrator, 3-crystal debug run, SVG + timing_report.json

Pass 2 — M-P1, M-P2:
  cache resume, N=10 ensemble, accumulated spectrum

Pass 3 — M-P4:
  cut_mixture, cylinder cut, ensemble_config.schema.md

Pass 4 — M-P5:
  JS LFFHessian.js analytical K (optional MMFF FD bypass)

Pass 5 — M-P6:
  N=100 timing report → decide relax optimization (§15)

Pass 6 — M-P8 (if relax slow):
  XPDB WebGPU relax option + parity vs MMFF

Pass 7 — M-P7, Si port
```

## 17. New milestone

| ID | Criterion |
|----|-----------|
| **M-P8** | Optional GPU/LFF relax: XPDB WebGPU vs MMFF — pos RMSD and ω spectrum within documented tolerance on 3 C NCs |

*Last updated: 2026-06-16 — implementation pass 1 spec, relax/GPU notes, Python MMFF-first policy, debug output layout.*

Updated notes are in [`doc/Topics/FTIR_Nanocrystals/Nanocrystal_pipeline.progress.md`](doc/Topics/FTIR_Nanocrystals/Nanocrystal_pipeline.progress.md) (§14–17). No code written.

## What changed in the plan

### Pass 1 — next coding session (when you say go)

Minimal pipeline for **3–5 C crystals** with full debug artifacts:

| Stage | Runtime |
|-------|---------|
| Generate + `01_init.npz` + `init.svg` | JS |
| MMFF relax → `02_relaxed.npz` | **Python only** |
| Topology → `03_topology.npz` | JS |
| MMFF FD Hessian → `04_hessian.npz` | Python |
| `eigh` + histogram → `05_spectrum.npz` + `spectrum.svg` | Python |
| `relaxed.svg` | JS |
| Ensemble sum → `out/plots/ensemble_spectrum.svg` | JS or Python |
| **`timing_report.json`** | Measured, not estimated |

Output layout is spelled out under `OUT_nc_ensemble/{data,out,work}/` in §14.3.

### Relaxation: Python MMFF first (your directive)

- **v1:** `nc_relax.py` → `MMFF.run(...)` (same settings as `test_nanocrystal_sparse_hessian.py`).
- **No GPU / XPDB / PD** in pass 1.
- Relax is **hypothesized** to be slowest; pass 1 **must prove that** via `timing_report.json`.

### GPU relax inventory (notes only — §15)

| Existing | Headless? | Notes |
|----------|-----------|-------|
| `run_xpdb_webgpu_headless.mjs` + `XPDB_WebGPU.js` | **Yes** | MMFFL springs, not full MMFF |
| `ScriptRunner.relaxJacobiGPU` + `PDSimulation` | Browser/molgui | Bond/PD constraints |
| `pyBall/OCL/LFFSolver.py` | Python+OpenCL | MD dynamics, `MAX_NEIGHBORS=8` |

**Future path (after pass 1 timing):** optional `--relax-backend xpdb_webgpu` only after new milestone **M-P8** (RMSD + ω parity vs MMFF on 3 crystals). MMFF stays the reference.

### Deferred to later passes

- Cylinder cut, full `cut_mixture`
- JS analytical K (`LFFHessian.js`)
- Probe weighting v2
- N=100 production run

## Canonical doc

Use `Nanocrystal_pipeline.progress.md` (not `Nanocrystal_pipeline.prograss.md`).

When you want implementation, say so and I'll start with Pass 1: orchestrator + 3-crystal debug run with SVGs, spectra, and timing report.

---

## 18. Pass 1 implemented (2026-06-16)

```bash
node scripts/run_nanocrystal_ensemble.mjs \
  --config scripts/ensemble.example.json \
  --data-dir OUT_nc_ensemble/data \
  --output-dir OUT_nc_ensemble/out \
  --n-crystals 3
```

Test run: `OUT_nc_ensemble_v2/` — 3 C sphere NCs, all stages + debug SVGs + `timing_report.md`.

Measured medians (crystals 2–3, ~164 atoms): relax ~123 ms, hessian ~172 ms, spectrum ~146 ms, generate ~59 ms, topology ~8 ms.

---

## 19. Wrap-up report — Pass 1 consolidation (2026-06-16)

Session goal: **end-to-end C diamond nanocrystal ensemble → MMFF relax → Hessian → dense `eigh` → accumulated IR-style spectra**, with **debug geometry** and a **local HTML viewer**. Planning doc §14–17 is largely implemented; below is what actually landed, how to run it, and what to commit.

### 19.1 Big picture (current architecture)

```mermaid
flowchart LR
  CFG[ensemble.example.json] --> ORCH[run_nanocrystal_ensemble.mjs]
  ORCH --> GEN[nanocrystals.mjs generate]
  GEN --> INIT[01_init.npz + init.xyz + init.mol2]
  INIT --> RELAX[pyBall.nanocrystal_pipeline relax]
  RELAX --> TOPO[exportFF.js]
  TOPO --> HESS[pyBall.nanocrystal_pipeline hessian]
  HESS --> SPEC[pyBall.nanocrystal_pipeline spectrum]
  SPEC --> ACC[pyBall.nanocrystal_pipeline accumulate]
  ORCH --> DBG[SVG + crystal JSON + viewer.html]
  ORCH --> ATLAS[--atlas atlas_shapes.json]
```

**Division of labour**

| Layer | Runtime | Responsibility |
|-------|---------|----------------|
| `scripts/run_nanocrystal_ensemble.mjs` | Node | Single orchestrator: sampling, cache paths, spawn gen/topology/Python, timing, debug artifacts, atlas mode |
| `scripts/nanocrystals.mjs generate` | Node | Crystal generation (consolidated; called by orchestrator) |
| `web/common_js/npzIO.js` | Node | NPZ read/write, `molToCrystalArrays`, `writeCrystalNpz`, XYZ/JSON crystal I/O |
| `web/common_js/MolIO.js` | Node | mol2 load, visualization bonds, `loadMMParamsFromDir`, `loadMolFromMol2`, `applyPositions`, `bondsForVisualization` |
| `web/common_js/exportFF.js` | Node | `buildTopologyNpz` (linearized topology NPZ) |
| `web/common_js/nanocrystalSvg.js` | Node | Cartesian-view SVG (single + init\|relaxed compare), atlas HTML table |
| `web/common_js/nanocrystalViewer.html` | Browser | Self-contained p5.js orthographic viewer (embedded JSON, `file://` safe) |
| `pyBall/nanocrystal_pipeline.py` | Python | `relax`, `hessian`, `spectrum`, `accumulate` subcommands |

**Strategic choices retained from §1–17**

- Dense `eigh` + histogram (v1); probe-weighted polarization spectra (partial v2).
- MMFF Python reference for relax + FD Hessian (no GPU relax in v1).
- NPZ stage cache under `data/crystals/<id>/` for throughput.
- Debug SVG/viewer optional; production runs use `--no-debug`.

### 19.2 Pipeline stages & on-disk layout

**Per crystal** (`<data-dir>/crystals/<crystal_id>/`):

| File | Stage | Producer |
|------|-------|----------|
| `01_init.npz` | generate | JS (`writeCrystalNpz`) |
| `init.mol2`, `init.xyz` | generate | JS |
| `meta.json` | generate | JS (`gen_params`, `bonds_ij` for viz) |
| `02_relaxed.npz`, `relaxed.xyz` | relax | Python MMFF |
| `03_topology.npz` | topology | JS (`buildTopologyNpz`) |
| `04_hessian.npz` | hessian | Python MMFF FD |
| `05_spectrum.npz` | spectrum | Python `eigh` + mode metadata |
| `*.status.json` | relax/hessian/spectrum | Python timing/status |

**Ensemble output** (`<output-dir>/`):

| Path | Content |
|------|---------|
| `accumulated/spectrum_histogram.npz` | Naive equal-weight mode histogram (cm⁻¹) |
| `accumulated/spectrum_pol_{x,y,z}.npz` | Probe-weighted spectra \|e·u_mode\|² |
| `plots/spectrum_histogram.png` | Line plot (not bar chart) |
| `plots/spectrum_pol_{x,y,z}.png` | Polarization line plots |
| `timing_report.md` | Per-stage min/median/max ms |
| `crystals/<id>/compare_{111,100,001}.svg` | Debug init vs relaxed (if debug on) |
| `crystals/<id>/crystal_{init,relaxed}.json` | Viewer geometry |
| `viewer.html` | Self-contained 3D viewer (debug crystals only) |
| `viewer_manifest.json` | Optional manifest for HTTP reload |

**Atlas mode** (`--atlas scripts/atlas_shapes.json`):

- Output: `<output-dir>/atlas/` — per-shape dirs, `index.html`, `viewer.html`, `compare.svg`, `gen_params.json`.
- No full ensemble spectrum; shape showcase only.
- Verified output tree: `OUT_nc_atlas/atlas/`.

**Test runs (local, do not commit)**

- `OUT_nc_ensemble_v2/` — 8-crystal ensemble with spectra + polarization plots.
- `OUT_nc_atlas/` — diamond shape atlas.

### 19.3 Key fixes & lessons (physics + tooling)

| Issue | Root cause | Fix |
|-------|------------|-----|
| H–C–H ~53° on dense spheres | `resolveClashes: 1` rotates H caps, breaks tetrahedral passivation | `ensemble.example.json` → `resolveClashes: 0`; orchestrator always passes flag |
| Absurd rigid modes in histogram | ω ≈ 1000 internal units treated as physical | `filter_vibrational_modes` before binning |
| Spectrum axis wrong units | Internal MMFF ω not cm⁻¹ | `omega_internal_to_cm1()` in Python pipeline |
| Histogram looked like bars | Matplotlib `bar` | Line plots (`plot_spectrum_line`) |
| Atlas index broken thumbnails | `index.html` used `atlas/<id>/…` prefix while already inside `atlas/` | Relative paths `sphere_R4p5/compare.svg` |
| Viewer CORS on `file://` | `fetch()` blocked | `installViewer()` embeds all crystal JSON in `viewer.html` |
| p5.js bond/axis flicker | WebGL `line()` + depth buffer at extreme zoom; ortho near/far tweaks unreliable in p5 | Reverted to `ortho(…, 0.1, 2000)` + `depthTest` off for bonds/axes; **good enough, stop here** |

### 19.4 Consolidation (what replaced scattered scripts)

**Merged / removed during consolidation** (do not resurrect):

- `web/molgui_webgpu/NanocrystalNpz.js`, `NanocrystalSvg.js` → `web/common_js/npzIO.js`, `nanocrystalSvg.js`
- `npzRead.js`, `npzWrite.js` → `npzIO.js`
- `scripts/build_topology_npz.mjs` → `exportFF.js` + `MolIO.js`
- `scripts/nc_relax.py`, `nc_build_hessian_mmff.py`, `nc_solve_spectrum.py` → `pyBall/nanocrystal_pipeline.py`

**Still reused unchanged**

- `tests/tSiNCs/nanocrystals.mjs generate`, `pyBall/nanocrystal_gen.py` (deprecated standalone scripts in `tests/tSiNCs/`)
- `tests/tSiNCs/nanocrystals.mjs topology` (logic in `exportFF.js` for headless NPZ write)
- `tests/tMMFF/test_nanocrystal_sparse_hessian.py` (MMFF reference)

### 19.5 Files to commit (recommended)

All currently **untracked** (`git status ??`):

| Path | Role |
|------|------|
| `tests/tSiNCs/nanocrystals.mjs` | **Main entrypoint** (unified CLI) |
| `tests/tSiNCs/run_nanocrystal_ensemble.mjs` | Deprecated wrapper |
| `tests/tSiNCs/ensemble.example.json` | Default ensemble config (8 spheres, debug views, `resolveClashes: 0`) |
| `tests/tSiNCs/atlas_shapes.json` | 8 C diamond shape presets for `--atlas` |
| `tests/tSiNCs/make_small_symmetric_nc.mjs` | Symmetric Si/C sphere batch (&lt;100 atoms) |
| `pyBall/nanocrystal_pipeline.py` | Python stages: relax, hessian, spectrum, accumulate |
| `web/common_js/npzIO.js` | Shared NPZ + crystal I/O |
| `web/common_js/exportFF.js` | Topology build + NPZ export |
| `web/common_js/MolIO.js` | Molecule I/O + bond extraction for viz |
| `web/common_js/nanocrystalSvg.js` | SVG debug + atlas HTML helpers |
| `web/common_js/nanocrystalViewer.html` | p5 viewer template (copied/embedded into outputs) |

**Suggested commit message**

```
Add nanocrystal ensemble pipeline (JS orchestrator + Python MMFF/eigh stages).

Single entrypoint run_nanocrystal_ensemble.mjs with NPZ caching, probe-weighted
spectra, debug SVG/viewer, and shape atlas mode. Consolidates prior scattered
nc_*.py and NanocrystalNpz/Svg JS into common_js + nanocrystal_pipeline.py.
```

**Do not commit**

- `OUT_nc_*` output trees (ensemble/atlas runs)
- Editor swap files, duplicate `classify_md copy*.py` under `scripts/`
- Large binary NPZ/plots unless explicitly desired as fixtures

**Optional follow-up doc edits** (separate commit or same):

- Update §7 table in this file: replace planned paths (`web/molgui_webgpu/…`, `scripts/nc_*.py`) with actual `web/common_js/…` and `pyBall/nanocrystal_pipeline.py`.
- Cross-link `Nanocrystal_generation.progress.md` § on `resolveClashes` (already documented there).

### 19.6 How to run

**Ensemble (8 crystals, full pipeline)**

```bash
node tests/tSiNCs/nanocrystals.mjs ensemble \
  --config tests/tSiNCs/ensemble.example.json \
  --data-dir tests/tSiNCs/OUT_nc_ensemble_v2/data \
  --output-dir tests/tSiNCs/OUT_nc_ensemble_v2/out \
  --work-dir tests/tSiNCs/OUT_nc_ensemble_v2/work
```

**Atlas only**

```bash
node tests/tSiNCs/nanocrystals.mjs ensemble \
  --atlas tests/tSiNCs/atlas_shapes.json \
  --output-dir tests/tSiNCs/OUT_nc_atlas
```

Open `tests/tSiNCs/OUT_nc_atlas/atlas/viewer.html` or `…/out/viewer.html` directly in browser (embedded data).

**Python stages standalone** (for debugging)

```bash
python3 -m pyBall.nanocrystal_pipeline relax --init-xyz … --out-npz …
python3 -m pyBall.nanocrystal_pipeline hessian --relaxed-xyz … --out-npz …
python3 -m pyBall.nanocrystal_pipeline spectrum --hessian-npz … --out-npz …
python3 -m pyBall.nanocrystal_pipeline accumulate --inputs … --out-dir …
```

Requires `tests/tMMFF/data` → `cpp/common_resources` symlink (handled by `ensure_mmff_data()`).

### 19.7 Deferred (unchanged from §16–17)

| Item | Notes |
|------|-------|
| Formal `ensemble_config.schema.md` | JSON works; schema doc not written |
| `cut_mixture` multi-weight sampling | Only first-weight sphere/planes in pass 1 |
| Cylinder cut | Not wired in orchestrator CLI |
| JS analytical Hessian (`LFFHessian.js`) | Python MMFF FD remains reference |
| GPU / XPDB relax (M-P8) | After N=100 timing report |
| Continuous probe-weighted spectrum v2 | Partial: `probe_weight_*` in `05_spectrum.npz`, line histograms only |
| Si element port | After C ensemble validated |
| Better 3D viewer | p5.js `line()` limits; consider Three.js only if viewer becomes production-critical |

### 19.8 Milestone status

| ID | Status |
|----|--------|
| **M-P1** (3–5 crystals, all stages, debug SVG, timing) | **Done** — `OUT_nc_ensemble_v2`, timing medians documented §18 |
| **M-P2** (ensemble accumulate + cm⁻¹ plots) | **Done** — histogram + pol x/y/z |
| **M-P3** (shape atlas) | **Done** — `atlas_shapes.json`, `OUT_nc_atlas` |
| **M-P4** (cut_mixture, cylinder, defects) | **Not started** |
| **M-P5** (JS analytical K) | **Deferred** |
| **M-P6** (N=100 timing) | **Not run** |
| **M-P8** (GPU relax parity) | **Deferred** |

*Pass 1 wrap-up complete. Next session: commit listed files, run N=100 timing (M-P6), or extend cut/defect sampling — not viewer polish.*

### 19.9 Physics path clarification — MMFF vs linearized (not mixed)

**Short answer:** Relax, Hessian, and spectrum use **full MMFF** (C++ backend via `pyBall.MMFF`) with **bonds + angles** enabled. They do **not** use the linearized spring topology for vibrations in pass 1.

| Stage | What runs | Bonds | Angles | Uses `03_topology.npz`? |
|-------|-----------|-------|--------|-------------------------|
| **Relax** (`02_relaxed`) | `MMFF.run()` | yes (`MMFF=+1`) | yes (`Angles=+1`) | no |
| **Hessian** (`04_hessian`) | `MMFF.getHessian3Nx3N` (finite-difference on full MMFF energy) | yes | yes | no |
| **Spectrum** (`05_spectrum`) | `eigh(K_projected)` from `04_hessian.npz` | — | — | no |
| **Topology** (`03_topology`) | JS `buildMMFFLTopology` → linearized sticks | K₁₂ only | **off** (`add_angle: false`, `add_dihedral: false`) | written, **not read** by spectrum |

So stage 3 is **cache-only / forward-looking** (planned “Path B”: assemble dense K from K₁₂/K₁₃/K₁₄ NPZ in JS or GPU). The doc §14 originally listed two paths:

- **Path A (implemented):** MMFF FD Hessian → `eigh` — **this is what the ensemble runs today.**
- **Path B (deferred):** `03_topology.npz` → analytical spring Hessian (`LFFHessian.js` / XPDB) — **not wired into `nanocrystal_pipeline.py`.**

**Why stage 3 still runs:** prepares linearized topology artifacts for future JS/GPU solver work (`Linearized_topology.progress.md`) and keeps the pipeline layout stable. It does **not** change relax or spectra until Path B is hooked up.

**Visualization bonds** (SVG/viewer) come from **mol2 topology** (`bondsForVisualization`), not from MMFFL linearized springs — unrelated to the vibration Hessian.

**If we want to avoid confusion:** either (a) skip stage 3 in production with `--no-topology` until Path B exists, or (b) enable `add_angle: true` in `buildTopologyNpz` so `03_topology.npz` at least matches the linearized design doc — still unused by spectrum until Path B lands.

---

Your instinct is right — the repo has **both** systems, but pass 1 **does not mix them** for vibrations.

### What actually drives relax + spectra today

**Full MMFF via Python → C++** (`pyBall.MMFF`), with **bonds and angles on**:

```103:103:pyBall/nanocrystal_pipeline.py
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
```

- **Relax** → `MMFF.run()` on `init.xyz`
- **Hessian** → `MMFF.getHessian3Nx3N` (finite differences on the **full** MMFF energy, bonds + angles)
- **Spectrum** → `eigh(K_projected)` from **`04_hessian.npz` only**

No linearized springs, no replacement of angles by 1–3 sticks on this path.

### What stage 3 (`03_topology.npz`) is

**Linearized MMFFL topology** (springs only) — built in JS via `buildMMFFLTopology`, but currently with **angles and dihedrals turned off**:

```79:86:web/common_js/exportFF.js
    const topo = buildMMFFLTopology(mol, mmParams, {
        type_source: 'table',
        add_angle: addAngle,
        add_dihedral: addDihedral,
        ...
        K12: 0, K13: 0, K14: 0,
    });
```

(`add_angle` / `add_dihedral` default to `false` in the orchestrator call.)

That file is **written to the cache but not read** by relax, hessian, or spectrum. It’s prep for the planned **Path B**: assemble K from K₁₂/K₁₃/K₁₄ sticks (JS/GPU / `LFFHessian.js`), documented in §14 but **not implemented** in the ensemble yet.

### Summary table

| Component | MMFF (bonds+angles) | Linearized springs |
|-----------|---------------------|-------------------|
| Relax | ✅ | ❌ |
| Hessian / spectrum | ✅ | ❌ |
| `03_topology.npz` | ❌ | ✅ (K₁₂ only today) |
| SVG/viewer bonds | mol2 bonds (viz only) | ❌ |

So we’re **not** accidentally using linearized topology for spectra — we’re **running both in parallel**, but only MMFF FD feeds `eigh`. Stage 3 is misleading until Path B is wired or skipped.

I added **§19.9** to `Nanocrystal_pipeline.progress.md` with this clarification. Sensible next cleanup: add `--no-topology` for production runs, or enable `add_angle: true` so `03_topology.npz` matches the linearized design while Path B is still deferred.

---

## 20. Script consolidation plan — unified `nanocrystals.mjs`

> **Relocated (2026-07-08):** all nanocrystal orchestration scripts now live under `tests/tSiNCs/`. Legacy `scripts/` copies remain as deprecated stubs. Historical paths below refer to the consolidation event; use `tests/tSiNCs/` for current runs.

**Goal:** Merge 6 scattered nanocrystal scripts into a single unified CLI with subcommands, moving reusable logic to shared modules. Old scripts are **deprecated but not deleted** until the unified script passes equivalent tests.

### 20.1 Current scripts — inventory & overlap analysis

| Script | Lines | Role | Reusable logic trapped inside? |
|--------|-------|------|-------------------------------|
| `scripts/gen_nanocrystals.mjs` | 278 | CLI generation (mol2+xyz output) | **DEPRECATED** — logic in `Nanocrystals.js`, CLI in `nanocrystals.mjs generate` |
| `scripts/run_nanocrystal_ensemble.mjs` | 479 | Ensemble orchestrator: gen→relax→topo→hessian→spectrum→accumulate | **Yes** — `mulberry32` (dup), `sampleValue`, `sampleCutSpec`, `crystalId`, `buildGenCli`, `writeTimingReport`, `writeCrystalCompareSvgs`, `writeCrystalViewerJson`, `installViewer`, `enrichViewerJsonWithTopology`, `appendManifest`, `runCmd`, `runPy` |
| `scripts/build_linearized_topology.mjs` | 262 | Standalone topology builder (mol2/xyz→topology.npz + HTML viewer) | **Yes** — `loadMolecule`, `buildPrimaryAdj`, `validateSpringGraphDistances` |
| `scripts/test_nanocrystal_geometry.mjs` | 264 | Geometry audit (bond lengths, angles, clashes, valence) | **Yes** — `auditGeometry`, `neighborIndices`, `bondedPairSet`, `angleDeg`, `dist`, `buildSiH4`, `oneLineSummary` |
| `scripts/debug_nanocrystal_nonbond_groups.mjs` | 123 | Nonbond parity check (brute vs grouped) + collision pair export | Minimal — mostly orchestration of `Nonbonded.js` / `CollisionWorkgroups.js` |
| `web/web_tests/test_ring_detection.mjs` | 160 | Ring detection test (adamantane + diamond NC defects) | No — already thin, imports from modules |

### 20.2 Key overlaps identified

1. **Generation called 3 ways:** `gen_nanocrystals.mjs` (standalone CLI), `run_nanocrystal_ensemble.mjs` (spawns `gen_nanocrystals.mjs` as subprocess then re-parses mol2), `test_ring_detection.mjs` (calls `generateNanocrystal()` directly in-process). The subprocess round-trip in the orchestrator is wasteful — should call `generateNanocrystal()` directly.

2. **Topology built 2 ways:** `build_linearized_topology.mjs` (full standalone with graph validation, HTML viewer, JSON export) vs `exportFF.js` `buildTopologyNpz()` (simpler, adds collision workgroups). The standalone has richer output; `exportFF.js` is stripped down.

3. **`mulberry32` duplicated:** `Nanocrystals.js` exports it; `run_nanocrystal_ensemble.mjs` has its own copy.

4. **SVG/viewer helpers trapped in orchestrator:** `writeCrystalCompareSvgs`, `writeCrystalViewerJson`, `installViewer`, `enrichViewerJsonWithTopology` — all reusable, should be in `nanocrystalSvg.js`.

5. **Geometry audit not reusable:** `auditGeometry()` is valuable for post-generation QC but trapped in a script. The ensemble pipeline could use it as a post-generation validation step.

6. **Ring detection not in pipeline:** `test_ring_detection.mjs` tests ring detection but it's not part of the ensemble. Ring analysis could be a debug stage.

### 20.3 Proposed unified CLI: `scripts/nanocrystal.mjs`

```
node scripts/nanocrystal.mjs <subcommand> [options]

Subcommands:
  generate   — single crystal generation (mol2+xyz output)
  ensemble   — batch ensemble pipeline (gen→relax→topo→hessian→spectrum→accumulate)
  topology   — linearized topology build (mol2/xyz/npz → topology.npz + viewer)
  rings      — ring detection on a crystal (SVG output)
  audit      — geometry audit (bond lengths, angles, clashes, valence)
  nonbond    — nonbond group debug + collision pair export
```

### 20.4 Dual CLI + JSON config (all subcommands)

**Every subcommand accepts both CLI flags and JSON config:**

- CLI flags: convenient for quick single-crystal operations, scripting, debugging
- `--config path.json`: JSON file with full parameter set, enables fine-grained control and reproducibility
- CLI flags **override** JSON values when both are present (CLI > JSON > defaults)
- JSON schema is subcommand-specific but shares common sections (cut, passivation, defects, debug)

**Example:**

```bash
# CLI only
node scripts/nanocrystal.mjs generate --cif C_diamond_sym.cif --cutMode sphere --sphereR 5 --caps H

# JSON config
node scripts/nanocrystal.mjs generate --config my_crystal.json

# JSON + CLI override (override sphereR from JSON)
node scripts/nanocrystal.mjs generate --config my_crystal.json --sphereR 8

# Ensemble with JSON
node scripts/nanocrystal.mjs ensemble --config ensemble.json --n-crystals 100 --data-dir OUT/data --output-dir OUT/out
```

**JSON config structure (shared sections):**

```json
{
  "cif": "cpp/common_resources/crystals/C_diamond_sym.cif",
  "cut": { "cutMode": "sphere", "sphereR": 5.0, "sphereNrep": 5 },
  "replication": { "nx": 2, "ny": 2, "nz": 2 },
  "passivation": { "caps": "H", "resolveClashes": false, "outwardBias": 0.35 },
  "defects": { "insertProb": 0.0, "collapseProb": 0.0 },
  "debug": { "enabled": true, "views": ["111", "100", "001"] }
}
```

### 20.5 Module consolidation — what moves where

| Module | What moves here | Status |
|--------|----------------|--------|
| `web/molgui_webgpu/Nanocrystals.js` | Already has generation logic + `mulberry32` + `defaultArgs`. Remove duplicate `mulberry32` from orchestrator. | ✅ Done |
| `web/common_js/nanocrystalSvg.js` | Move `writeCrystalCompareSvgs`, `writeCrystalViewerJson`, `installViewer`, `enrichViewerJsonWithTopology` from orchestrator. Already has `bondsForSvgMol`, `exportCrystalSvg`, `exportCrystalCompareSvgViews`, `atlasIndexHtml`. | **TODO** |
| `web/common_js/exportFF.js` | Already has `buildTopologyNpz`. Add `validateSpringGraphDistances` from `build_linearized_topology.mjs` as optional validation step. | **TODO** |
| **NEW: `web/common_js/NanocrystalAudit.js`** | Move `auditGeometry`, `neighborIndices`, `bondedPairSet`, `angleDeg`, `dist`, `buildSiH4`, `oneLineSummary` from `test_nanocrystal_geometry.mjs`. | **TODO** |
| **NEW: `web/common_js/NanocrystalEnsemble.js`** | Move `sampleValue`, `sampleCutSpec`, `crystalId`, `buildGenParams`, `writeTimingReport`, `appendManifest`, `runCmd`, `runPy` from orchestrator. These are reusable ensemble utilities. | **TODO** |

### 20.6 What stays in scripts

| Script | What remains |
|--------|-------------|
| `scripts/nanocrystal.mjs` (new unified) | Subcommand dispatch + CLI/JSON parsing per subcommand. ~50-80 lines per subcommand (thin glue). |
| `scripts/gen_nanocrystals.mjs` | **DEPRECATED** — header note: "Absorbed into `nanocrystal.mjs generate`. Do not delete until `nanocrystal.mjs generate` passes equivalent tests." |
| `scripts/run_nanocrystal_ensemble.mjs` | **DEPRECATED** — header note: "Absorbed into `nanocrystal.mjs ensemble`. Do not delete until tested." |
| `scripts/build_linearized_topology.mjs` | **DEPRECATED** — header note: "Absorbed into `nanocrystal.mjs topology`. Do not delete until tested." |
| `scripts/test_nanocrystal_geometry.mjs` | **DEPRECATED** — header note: "Absorbed into `nanocrystal.mjs audit`. Do not delete until tested." |
| `scripts/debug_nanocrystal_nonbond_groups.mjs` | **DEPRECATED** — header note: "Absorbed into `nanocrystal.mjs nonbond`. Do not delete until tested. Collision group export is important for downstream simulation." |
| `web/web_tests/test_ring_detection.mjs` | **Keep for now** — already thin test script. Ring detection also available via `nanocrystal.mjs rings`. |

### 20.7 Key design decisions

1. **No subprocess for generation** — call `generateNanocrystal()` directly in-process. The ensemble orchestrator currently spawns `gen_nanocrystals.mjs` then re-reads the mol2 — wasteful. The unified `ensemble` subcommand calls `generateNanocrystal()` directly like `test_ring_detection.mjs` does.

2. **Python stages stay as subprocess** — `nc_relax`, `nc_build_hessian`, `nc_solve_spectrum` via `pyBall.nanocrystal_pipeline` are C++ backed, can't run in Node.

3. **Dual CLI + JSON for all subcommands** — JSON provides fine-grained control (especially for ensemble configs with `cut_mixture` arrays, probability distributions, etc.); CLI flags are convenient for quick operations. CLI overrides JSON when both present.

4. **Pipeline stages as composable functions** — each stage (generate, relax, topology, hessian, spectrum, accumulate) is a function in a shared module. The `ensemble` subcommand chains them. Individual subcommands call just one.

5. **Debug output unified** — SVG, viewer HTML, and timing report all go through `nanocrystalSvg.js` functions. No SVG/viewer code in the CLI script.

6. **Collision groups and nonbond pair export preserved** — `nanocrystal.mjs nonbond` subcommand keeps the collision group debug + pair export functionality. This is important for downstream simulation (PD collision constraints, GPU kernel pair lists).

7. **Deprecated scripts get header notes** — each old script gets a clear deprecation header:
   ```javascript
   /// @deprecated — Absorbed into scripts/nanocrystal.mjs <subcommand>.
   /// Do not delete until nanocrystal.mjs passes equivalent tests.
   /// Kept for backward compatibility during transition.
   ```

### 20.8 Subcommand → shared module mapping

| Subcommand | Primary module calls | Python stages |
|------------|---------------------|---------------|
| `generate` | `Nanocrystals.js: generateNanocrystal()`, `MolIO.js: loadMMParamsFromDir()`, `MoleculeIO.js: toMol2String/toXYZString` | None |
| `ensemble` | `NanocrystalEnsemble.js: sampleValue/sampleCutSpec/crystalId/runCmd/runPy/writeTimingReport`, `Nanocrystals.js: generateNanocrystal()`, `exportFF.js: buildTopologyNpz()`, `nanocrystalSvg.js: writeCrystalCompareSvgs/installViewer` | `pyBall.nanocrystal_pipeline: relax/hessian/spectrum/accumulate` |
| `topology` | `exportFF.js: buildTopologyNpz()`, `MMFFLTopology.js: buildMMFFLTopology/packLinearTopologyForGPU`, `LinearizedTopologyNpz.js: writeTopologyNpzFile`, `LinearizedTopologyViewer.js: exportStickViewerHTML` | None |
| `rings` | `Graph.js: runRingsOnMol()`, `nanocrystalSvg.js: bondsForSvgMol/exportCrystalSvg`, `Nanocrystals.js: generateNanocrystal` (if generating) | None |
| `audit` | `NanocrystalAudit.js: auditGeometry/oneLineSummary`, `MolIO.js: loadMolFromMol2` | None |
| `nonbond` | `CollisionWorkgroups.js: buildCollisionWorkgroups/buildExclIcol_1_2_3`, `Nonbonded.js: computeNonbondBruteForceKernelStyle/computeNonbondByGroups/maxAbsDiff`, `MolIO.js`, `npzIO.js` | None |

### 20.9 Implementation order

1. **Create `NanocrystalAudit.js`** — extract `auditGeometry` + helpers from `test_nanocrystal_geometry.mjs`
2. **Create `NanocrystalEnsemble.js`** — extract ensemble utilities from `run_nanocrystal_ensemble.mjs`
3. **Extend `nanocrystalSvg.js`** — move SVG/viewer helpers from orchestrator
4. **Extend `exportFF.js`** — add `validateSpringGraphDistances` as optional validation
5. **Create `scripts/nanocrystal.mjs`** — unified CLI with all 6 subcommands
6. **Add deprecation headers** to 5 old scripts
7. **Test each subcommand** against old script output for parity
8. **Delete old scripts** only after all subcommands pass equivalent tests

### 20.10 Testing strategy for consolidation

| Test | How |
|------|-----|
| `generate` parity | Run `gen_nanocrystals.mjs` and `nanocrystal.mjs generate` with same args → compare mol2/xyz output |
| `ensemble` parity | Run `run_nanocrystal_ensemble.mjs` and `nanocrystal.mjs ensemble` with same config → compare NPZ outputs + timing report |
| `topology` parity | Run `build_linearized_topology.mjs` and `nanocrystal.mjs topology` with same input → compare topology.npz |
| `audit` parity | Run `test_nanocrystal_geometry.mjs` and `nanocrystal.mjs audit` with same preset → compare geometry_report.json |
| `nonbond` parity | Run `debug_nanocrystal_nonbond_groups.mjs` and `nanocrystal.mjs nonbond` with same args → compare force diff + pair output |
| `rings` parity | Run `test_ring_detection.mjs` standalone and `nanocrystal.mjs rings` with same input → compare SVG output |

*Last updated: 2026-06-22 — §20 consolidation plan (review only, no code yet).*