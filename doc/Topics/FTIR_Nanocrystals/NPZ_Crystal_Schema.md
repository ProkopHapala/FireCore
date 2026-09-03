# NPZ Crystal & Topology Schema (v1.2)

**Contract document** for pipeline writers and consumers (JS export, Python relax/Hessian/spectrum, C++/Python viewers). Implementations must stay compatible with [`web/common_js/npzIO.js`](../../../web/common_js/npzIO.js), [`pyBall/io/crystal_npz.py`](../../../pyBall/io/crystal_npz.py), and [`cpp/common/io/`](../../../cpp/common/io/).

**Endianness:** little-endian (`<f8`, `<f4`, `<i4`, `<i8`). **Compression:** deflate (zip method 8), same as NumPy `np.savez_compressed`.

**Pipeline guide:** [`Nanocrystal_NPZ_Pipeline.guide.md`](Nanocrystal_NPZ_Pipeline.guide.md).

---

## Pipeline stages and filenames

| Stage | Canonical files | Topology mutates? | Coordinates |
|-------|-----------------|-------------------|-------------|
| **1 — Generate** | `01_crystal.npz`, `03_topology.npz`, `meta.json` | **Defines** connectivity + FF params | Init geometry |
| **2 — Relax** | `02_relaxed.npz`, `relaxed.xyz` (opt.) | **No** | Updated |
| **3 — Hessian** | `04_hessian.npz` | **No** (reads `03`) | Uses `02` `pos` |
| **4 — Spectrum** | `05_spectrum.npz` | — | — |

`03_topology.npz` may contain a `pos` snapshot from export time; it is **not** required to equal `02_relaxed.pos`. Hessian uses **exported spring tables** evaluated at relaxed coordinates.

---

## File kinds (consumer index)

| Kind | Typical filename | Show in viewer grid? | Primary consumers |
|------|------------------|----------------------|-------------------|
| Crystal geometry | `01_crystal.npz`, `02_relaxed.npz` | **yes** | Viewers, relax |
| Topology + FF + surface | `03_topology.npz` | **yes** | Viewers, Hessian |
| Hessian | `04_hessian.npz` | **no** | Spectrum stage |
| Spectrum | `05_spectrum.npz` | **no** | Plots, ensemble |
| Status | `*.status.json` | — | Humans, CI |

---

## 1. Crystal geometry NPZ (`01_*`, `02_*`)

**Writers:** `writeCrystalNpz()` (JS), `nanocrystal_pipeline.relax` (Python).  
**Readers:** `readCrystalNpz()` (JS), `load_crystal_npz()` (Python), `CrystalNpz.h` (C++).

| Key | dtype | shape | Required | Description |
|-----|-------|-------|----------|-------------|
| `pos` | float64 | (N, 3) | yes | Cartesian coordinates (Å) |
| `Z` | int32 | (N,) | yes | Atomic numbers |
| `natoms` | int32 | scalar | recommended | N (= `len(Z)`) |
| `bonds_ij` | int32 | (Nb, 2) | **yes** | Covalent bond pairs, **0-based** atom indices (not distance-inferred) |
| `gen_params` | uint8 | (L,) | no | UTF-8 JSON generation metadata |
| `timing_ms` | float64 | scalar | no | Stage wall time (ms) |
| `fmax` | float64 | scalar | no | Post-relax max \|force\| (eV/Å), `02` only |
| `converged` | bool (`\|b1`) | scalar | no | Relax converged flag, `02` only |
| `n_steps` | int64 (`<i8`) | scalar | no | MMFF relax step count, `02` only |

**Bonds:** required. Writers (`writeCrystalNpz`, `relax --init-mol2/--init-xyz/--init-npz`) fail loud if topology is missing. Distance inference drops 5-ring closers (~2.5 Å C–C, longer Si–Si). Atlas also writes `01_init.npz` (same schema; alias of `01_crystal.npz`).

**Relax contract:** `02_relaxed.npz` must preserve `Z` and `bonds_ij` from init mol2/xyz/npz; only `pos` (and relax metadata) change.

---

## 2. Topology + linearized FF + surface groups (`03_topology.npz`)

**Writers:** `buildTopologyNpzFull()` / `exportNanocrystalBundle()` (JS).  
**Readers:** `load_topology_npz()` (Python), `TopologyNpz.h` + `CrystalNpz.h` (C++), `readNpzFile` (JS).

**Immutability:** treat all keys below as **fixed parameters** after export unless a topology edit (add/remove atom/bond) occurs.

### 2a. MMFFL spring packing (Hessian input)

Harmonic sticks: \(E = \frac{1}{2} k (r - l_0)^2\) between endpoint atoms in each slot. `stick_class` distinguishes primary bonds vs angle/dihedral replacement sticks (K₁₂/K₁₃/K₁₄).

| Key | dtype | shape | Description |
|-----|-------|-------|-------------|
| `pos` | float64 | (N, 3) | Geometry **at export** (reference snapshot) |
| `Z` | int32 | (N,) | Atomic numbers (must match crystal) |
| `neigh_idx` | int32 | (N, M) | Neighbor atom index per slot (−1 = empty) |
| `neigh_count` | int32 | (N,) | Active neighbor slots per atom |
| `bond_l0` | float32 | (N, M) | Equilibrium distance per slot (Å) |
| `bond_k` | float32 | (N, M) | Stiffness per slot (eV/Å²) |
| `KLs` | float32 | (N, M, 2) | Packed [l0, K] (duplicate of bond_*) |
| `stick_class` | int32 | (N, M) | 0=empty, 1=bond, 2=angle, 3=dihedral |
| `max_neighbors` | int32 | scalar | M (typically 48) |
| `n_bond` | int32 | scalar | Primary bond spring count |
| `n_angle` | int32 | scalar | Angle replacement stick count |
| `n_dihedral` | int32 | scalar | Dihedral replacement stick count |

**Hessian assembly:** iterate slots with `stick_class > 0`, unique pairs `i < j`, use `bond_k[i,s]`, `bond_l0[i,s]` with **current** `pos` from `02_relaxed.npz`. Implementation: `FTIR.build_hessian_from_linear_topology`.

### 2b. Surface collision workgroups (viewer AABB + future nonbond)

From [`web/common_js/CollisionWorkgroups.js`](../../../web/common_js/CollisionWorkgroups.js):

| Key | dtype | shape | Description |
|-----|-------|-------|-------------|
| `radius` | float64 | (N,) | VdW radius per atom (Å) |
| `icol` | int32 | (N,) | Collision group id (−1 = bulk) |
| `icolGroup` | int32 | (N,) | Display / color group index |
| `group_atoms` | int32 | (G, group_cap) | Atom indices per surface group (−1 pad) |
| `group_nAtoms` | int32 | (G,) | Atom count per group |
| `group_bbox_min` | float64 | (G, 3) | AABB min corner per group |
| `group_bbox_max` | float64 | (G, 3) | AABB max corner per group |
| `n_groups` | int32 | scalar | G |
| `group_cap` | int32 | scalar | Max atoms per group (default 32) |
| `excl_icol` | int32 | (N, excl_max) | Excluded collision partner group ids |
| `excl_count` | int32 | (N,) | Exclusion list length |

**Viewer rule:** draw AABB wireframes from `group_bbox_min` / `group_bbox_max` when `icolGroup` present (`[g]` in C++ MolBrowser).

### 2c. Optional metadata

| Key | dtype | Description |
|-----|-------|-------------|
| `schema_version` | int32 | Start at `1` |
| `source_mol2` | uint8 | UTF-8 path string |
| `defects_json` | uint8 | UTF-8 JSON: `{nCollapsed, nInserted, nFused, …}` |

---

## 3. Hessian NPZ (`04_hessian.npz`)

**Writer:** `nanocrystal_pipeline build_hessian_bundle`.  
**Readers:** `nanocrystal_pipeline spectrum`; analysis scripts.

| Key | dtype | shape | Description |
|-----|-------|-------|-------------|
| `K` | float64 | (3N, 3N) | Stiffness matrix (Hessian of harmonic network) |
| `K_projected` | float64 | (3N, 3N) | `K` after rigid-mode shift (`FTIR.project_rigid_modes`) |
| `M` | float64 | (3N, 3N) | Diagonal mass matrix (amu), block 3×3 per atom |
| `pos` | float64 | (N, 3) | Coordinates used for assembly (**from `02_relaxed`**) |
| `Z` | int32 | (N,) | From `01_crystal` / topology (authoritative) |
| `natoms` | int32 | scalar | N |
| `bonds_ij` | int32 | (Nb, 2) | Copied from init crystal when present |
| `source` | unicode | scalar | `topology_linear` (default) or `mmff_fd` (diagnostic) |
| `timing_ms` | float64 | scalar | Assembly time |
| `parity_max_rel_diff` | float64 | scalar | ‖K_topo−K_mmff‖/‖K_mmff‖ when `--compare-mmff`; else NaN |
| `topology_path` | unicode | scalar | Provenance path to `03_topology.npz` |

**Cached from `03_topology.npz` (surface / nonbond, for downstream):**  
`group_bbox_min`, `group_bbox_max`, `icolGroup`, `icol`, `radius`, `group_atoms`, `group_nAtoms`, `n_groups`, `group_cap`, `excl_icol`, `excl_count` — same dtypes as §2b.

**Spectrum stage** uses `K_projected` and `M` by default.

---

## 4. Spectrum NPZ (`05_spectrum.npz`)

**Writer:** `nanocrystal_pipeline solve_spectrum`.

| Key | dtype | shape | Description |
|-----|-------|-------|-------------|
| `omega_centers` | float64 | (n_bins,) | Histogram bin centers (cm⁻¹) |
| `hist` | float64 | (n_bins,) | Mode-count histogram (or weighted if extended) |
| `omegas_modes` | float64 | (3N,) | All eigenfrequencies (√(ω²) internal units → stored as ω) |
| `omegas_modes_vib` | float64 | (n_vib,) | Physical modes after rigid-artifact filter |
| `probe_weight_x` | float64 | (3N,) | \|e·u\|² probe weight per mode (x polarization) |
| `probe_weight_y` | float64 | (3N,) | y polarization |
| `probe_weight_z` | float64 | (3N,) | z polarization |
| `grid_meta` | unicode | scalar | JSON: bin range, margin, `spectrum_kind` |
| `timing_ms` | float64 | scalar | Diagonalization + histogram time |
| `sigma_bins` | float64 | scalar | Gaussian blur σ (bins); 0 = none |
| `units` | unicode | scalar | `cm-1` |

**Note:** eigenvectors are **not** persisted in v1.2 (computed internally for probe weights). Future versions may add `eigenvectors` (3N×3N) for mode visualization.

---

## 5. Status JSON (`<stage>.status.json`)

Plain JSON beside stage NPZ. Fields vary by stage; always includes `status: "ok"` on success.

| Stage | Example keys |
|-------|----------------|
| relax | `natoms`, `fmax`, `converged`, `n_steps`, `timing_ms` |
| hessian | `hessian_source`, `ndof`, `n_sticks`, `parity_max_rel_diff`, `topology_npz`, `relaxed_npz` |
| spectrum | `n_modes`, `omega_min_cm1`, `omega_max_cm1`, `spectrum_kind_hist` |

---

## 6. Standalone `.npy` fast path

Uncompressed single arrays for mmap-heavy viewers:

- `pos.npy`, `Z.npy`, `bonds_ij.npy`

Full `.npz` remains canonical for pipeline stages.

---

## 7. Viewer listing rules

When browsing a pipeline output directory, **list only geometry/topology stages** in the thumbnail grid:

| Pattern | Show in browser? |
|---------|------------------|
| `01_*.npz`, `02_*.npz`, `03_*.npz` | **yes** |
| `*hessian*.npz`, `*spectrum*.npz`, `04_*.npz`, `05_*.npz` | **no** |
| `pos.npy` + companions | yes (`Z.npy` / `bonds_ij.npy` hidden when merged) |

C++: `BrowserView::filterNonViewerNpzFiles()`. Python: `VispyMolBrowser` mirrors this.

**C++ NPZ reader:** decodes `<f8`, `<f4`, `<i8`, `<i4`, `|u1`, `|b1`; **skips** `<U*>` / `|S*` with warning — see [`CPP_MolecularBrowser_NPZ.md`](CPP_MolecularBrowser_NPZ.md).

---

## 8. Validation

**Python crystal:**

```python
from pyBall.io.crystal_npz import load_crystal_npz, load_topology_npz, validate_topology_crystal_parity
c = load_crystal_npz("02_relaxed.npz")
t = load_topology_npz("03_topology.npz", full=True)
validate_topology_crystal_parity({"pos": c["pos"], "Z": load_crystal_npz("01_crystal.npz")["Z"], "bonds_ij": c.get("bonds_ij")}, t)
```

**Node crystal:**

```bash
node -e "
import { readCrystalNpz } from './web/common_js/npzIO.js';
import fs from 'fs';
const c = readCrystalNpz(fs, 'path/to/01_crystal.npz');
console.log(c.natoms, c.Z.length);
"
```

---

## Changelog

| Version | Date | Change |
|---------|------|--------|
| 1.2 | 2026-07-07 | Full `04_hessian` / `05_spectrum` tables; topology immutability; pipeline stages; consumer rules |
| 1.1 | 2026-07-07 | Relaxed `converged`/`n_steps` dtypes; viewer listing; C++ skip unicode |
| 1 | 2026-07-07 | Initial contract for parallel agent work |
