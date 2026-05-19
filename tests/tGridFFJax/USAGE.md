# GridFF-JAX — How to use it for ANY substrate

This is the operational guide for the `pyBall.gridff_jax` framework and its
benchmark driver `tests/tGridFFJax/benchmark_substrate_6d.py`. It is written
so that any new student, collaborator, or future-you can take a substrate
+ a quantum-chemistry / MLIP teacher and produce a fitted GridFF — without
reading any code.

If you are looking for the **research results** (CHONH2, PTCDA, NaCl
findings), see `runs/SUMMARY_6D_apples.md` instead.

---

## Status: v1 — what's done in this cleanup pass

These flaws (from `tests/tGridFFJax/Flaws_Findings.md`, both sections) are
fixed:

- Positional `FeatureToggles(...)` everywhere → keyword args (the
  silent-`use_reactive_grid=True` bug is gone)
- Substrate-class owns the recipe: `metal` gets London damping + pairwise C6
  + raw 1/r⁶ Stage-2; `ionic` gets none of those by default. The driver no
  longer hardcodes a metal recipe onto ionic substrates.
- QEQ + image-charge energy double-counting fixed in both NumPy and JAX
  paths. The model now sums `energy + image_energy` without overlap.
- `save_apples()` persists `adsorbate.charges` (was silently dropped on
  `--load-apples`, which rewrote the Hamiltonian).
- Chunked teacher evaluation preserves `backend / model_path / tile /
  interaction_energy / forces_available` metadata across chunks.
- `SurfaceXYZBackend` preserves the full lattice (off-diagonals included);
  no more silent diagonalisation of skewed cells.
- Ionic substrate without non-zero point charges now raises (hard error).
- Gas-phase molecule reference cache now keyed on geometry + charges +
  anchor, not just `adsorbate.name`. ASE `Atoms` carries `initial_charges`
  for charge-aware calculators.
- Cube backend potentials are now converted **Hartree → eV**, not divided
  by Bohr³. Density still gets the e/Bohr³ → e/Å³ conversion.
- Cube backend uses ASE's full periodic table (no hardcoded subset).
- `np.fromstring` (deprecated) replaced with `np.fromiter` in cube reader.
- Cache signature now fingerprints teacher kind, model path, tiling,
  interaction_energy mode, density backend kind and file paths.
- Stratified split honors zero `--val-fraction` / `--test-fraction`.
- `GenericASECalcBackend` no longer swallows constructor `TypeError`;
  `teacher_tile` accepts `'auto'` as a string.
- `PrecomputedTeacherBackend` validates force shape `(n_pose, n_atom, 3)`
  against the pose batch's adsorbate atom count.
- `--prefer-jax` is `BooleanOptionalAction` (you can `--no-prefer-jax`).
- All hardcoded workstation paths (CHGCAR / LOCPOT / MACE model) now come
  from `GRIDFF_JAX_CHGCAR / GRIDFF_JAX_LOCPOT / GRIDFF_JAX_MODEL_PATH`
  environment variables.
- FireCore export now writes the dispersion grid + per-element C₆ + offset
  alongside the PLQ grid, and warns loudly when the fitted student
  includes a C₆ channel that PLQ alone can't reproduce.
- A complete `HybridParameters_full.npz` checkpoint is written for every
  run (chi, hardness, image plane, reservoir terms, all per-element
  parameters), so a fit is now fully reproducible from disk.
- Export failure is fail-loud (the run aborts). Set
  `GRIDFF_JAX_EXPORT_LENIENT=1` only for debugging.
- `TrainingConfig.batch_size` and `HybridModelConfig.use_qeq/use_image/
  use_reactive_grid` are documented as informational only — runtime
  behaviour is driven by `FeatureToggles`.

The full design direction for the next step (bundle-driven plug-and-play
workflow, `ionic_dft_volumetric` recipe, validated manifest, hash-stable
input identity) is in `tests/tGridFFJax/new_PLAN.md`.

## 1. What the framework does

Input:
1. **A substrate** — a slab of metal (Ag, Au, Cu, …) or an ionic crystal
   (NaCl, KBr, CsCl, …), in either VASP CHGCAR/LOCPOT, Gaussian-cube, or
   point-charge xyz form.
2. **A teacher** — interaction energies (and ideally forces) of a probe
   adsorbate (H₂O, CHONH₂, PTCDA, …) at many poses above the substrate.
   The teacher may be any DFT code (VASP, Psi4, Gaussian, Q-Chem, …) or
   any MLIP (MACE, ANI, NEP, AIMNet2, schnetpack, …).

Output:
1. A **`Bspline_PLQd.npy`** grid (Pauli / London / Coulomb channels) that
   FireCore's C++/OpenCL runtime reads at MD time. Same format as the
   existing point-charge GridFFs in `cpp/common_resources/<surface>/`.
2. A 11-figure diagnostic pack showing energy/force parity, per-z RMSE,
   per-element error breakdown, lateral 2D slice, and channel decomposition.
3. An `apples_dataset.npz` you can re-use to re-fit with different
   `--fit-mode`, weights, etc., without re-running the slow teacher.

---

## 2. The decision tree

```
                ┌─── Conductor (Ag, Au, Cu, Pt, …)?
                │       → --substrate-class metal
                │       → density: VASP CHGCAR + LOCPOT  (--density-kind vasp_volumetric)
                │       → builder_mode = metal_density_plq
                │       → physics: PLQ from electron density + LOCPOT image plane
                │
   Substrate ──┼─── Ionic crystal (NaCl, KBr, CsCl, …)?
                │       → --substrate-class ionic
                │       → density: point-charge xyz  (--density-kind surface_xyz, default)
                │                    or DFT cube      (--density-kind cube)
                │       → builder_mode = parity_core
                │       → physics: PLQ + Ewald-slab Coulomb (no work function)
                │
                └─── Other (semiconductor, oxide, …)?
                        Phase B — open a ticket; you'll need a new
                        substrate-class entry in pyBall/gridff_jax/substrate_classes.py.
```

```
   Teacher ──┬─── Run an MLIP at runtime?
              │       → --teacher-kind generic_calc
              │       → --teacher-module <python_module> (e.g. mace.calculators)
              │       → --teacher-callable <Class>       (e.g. MACECalculator)
              │       → --model-path /path/to/model
              │       Works for: MACE, ANI (TorchANI), NEP, AIMNet2,
              │       schnetpack, MattersSim, Allegro, MACE-MP-0 …
              │       Anything with an ASE-compatible calculator.
              │
              ├─── MAD-SURF MACE (coin metals ONLY)?
              │       → --teacher-kind madsurf
              │       (this is a thin specialisation of generic_calc)
              │
              └─── Pre-computed DFT (Psi4 / VASP / Gaussian / …) on disk?
                      → --teacher-kind precomputed
                      → --teacher-npz /path/to/dataset.npz
                      Generate the npz with `dft_to_npz_template.py`
                      (see section 5 below).
```

---

## 3. Required Python packages

| Package | Why | Source |
|---|---|---|
| `numpy`, `scipy`, `matplotlib` | core | pip |
| `jax`, `jaxlib`, `optax` | substrate-field assembly + Stage 1 fit | pip |
| `ase` | atomic structure + universal calculator interface | pip |
| `mace-torch` (optional) | only for --teacher-kind madsurf or MACE | `pip install mace-torch` |
| `torchani` (optional) | only for --teacher-callable TorchANI calc | pip |
| any other MLIP package | if you want runtime MLIP teacher | their docs |

For pre-computed DFT (`--teacher-kind precomputed`) you need **none of the
MLIP packages** — just numpy. Convenient when running on machines that don't
have CUDA / torch.

---

## 4. Common recipes (copy-paste-runnable)

### 4.1 Ag(111) + CHONH₂ — the original metal recipe

```bash
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/benchmark_substrate_6d.py \
    --substrate-class metal \
    --chgcar /path/to/Ag/CHGCAR \
    --locpot /path/to/Ag/LOCPOT \
    --adsorbate CHONH2 \
    --teacher-kind madsurf \
    --model-path tests/tGridFFJax/mad-surf_data/models/full_dataset_config_weights/MACE_model.model \
    --device cuda \
    --n-u 5 --n-v 5 --n-z 10 --n-orient 3 --n-yaw 2 \
    --out-dir tests/tGridFFJax/runs/ag6d_chonh2
```

### 4.2 Ag(111) + ANY molecule, ANY MLIP teacher

Replace `--teacher-kind madsurf …` with:

```bash
    --teacher-kind generic_calc \
    --teacher-module mace.calculators \
    --teacher-callable MACECalculator \
    --model-path /path/to/your.model \
```

For TorchANI:

```bash
    --teacher-kind generic_calc \
    --teacher-module torchani \
    --teacher-callable ase   \   # torchani.models.ANI2x().ase()
    # (you may need to wrap your specific model — see section 7)
```

### 4.3 NaCl + H₂O — ionic from DFT (the user's student's case)

**Step 1**: Compute interaction energies (and ideally forces) for the
6D pose grid in your DFT code of choice. Use the template script:

```bash
# Edit dft_to_npz_template.py: fill in the run_psi4() body to match your group's
# psi4 invocation (or use --extractor coulomb_only / vasp_template).
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/dft_to_npz_template.py \
    --substrate-xyz cpp/common_resources/xyz/NaCl_1x1_L1_jax.xyz \
    --adsorbate H2O \
    --extractor psi4_template \
    --out tests/tGridFFJax/runs/nacl_psi4_dataset.npz \
    --n-u 4 --n-v 4 --n-z 8 --n-orient 1 --n-yaw 1
```

This generates the SAME 6D poses the benchmark driver will use later,
runs your DFT extractor on each pose, and writes the npz.

**Step 2**: Feed the npz to the benchmark driver:

```bash
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/benchmark_substrate_6d.py \
    --substrate-class ionic \
    --density-xyz cpp/common_resources/xyz/NaCl_1x1_L1_jax.xyz \
    --adsorbate H2O \
    --teacher-kind precomputed \
    --teacher-npz tests/tGridFFJax/runs/nacl_psi4_dataset.npz \
    --n-u 4 --n-v 4 --n-z 8 --n-orient 1 --n-yaw 1 \
    --out-dir tests/tGridFFJax/runs/nacl_ionic_psi4
```

⚠ **The pose-grid CLI args (n-u, n-v, n-z, n-orient, n-yaw, z-min, z-max,
z-bias, random-fraction, seed) MUST MATCH between the two scripts**,
otherwise the precomputed teacher will not find the queried poses and
will raise a clear `LookupError`.

**Step 3** (optional): compare against the existing point-charge GridFF:

```bash
python3 tests/tGridFFJax/compare_nacl_grids.py \
    --ref cpp/common_resources/NaCl_1x1_L1/Bspline_PLQd.npy \
    --new tests/tGridFFJax/runs/nacl_ionic_psi4/Bspline_PLQd.npy \
    --out tests/tGridFFJax/runs/nacl_ionic_psi4/compare_vs_pointcharge.png
```

### 4.4 ANY substrate + ANY DFT — the general pattern

```
  1. Pick a substrate. Put its xyz (with point charges if ionic) OR its
     DFT density (CHGCAR/LOCPOT, or two .cube files) in a known path.

  2. Decide the substrate class:
        metal  → --substrate-class metal --chgcar <X> --locpot <Y>
        ionic  → --substrate-class ionic --density-xyz <Z>
                  (or --density-kind cube + --chgcar density.cube --locpot pot.cube)

  3. Decide the teacher source. If you have DFT energies in any format,
     write them into the npz schema (see dft_to_npz_template.py) and use
     --teacher-kind precomputed. If you have an MLIP ready, use
     --teacher-kind generic_calc.

  4. Run the benchmark driver. Output goes under --out-dir.
```

---

## 5. The pre-computed teacher npz schema

`PrecomputedTeacherBackend` reads a single `.npz` with these arrays:

| Array | Shape | Dtype | Units | Required? | Meaning |
|---|---|---|---|---|---|
| `pose_params` | (n_pose, 7) | float64 | dimensionless | YES | (u, v, z, qw, qx, qy, qz) — fractional lateral + Å height + unit quaternion |
| `energies` | (n_pose,) | float64 | eV | YES | Interaction energy E_complex − E_slab − E_mol |
| `forces` | (n_pose, n_atom, 3) | float64 | eV/Å | NO | Forces on the molecule's atoms only |
| `positions` | (n_pose, n_atom, 3) | float64 | Å | NO | Adsorbate positions in lab frame (saved for traceability) |
| `adsorbate_symbols` | (n_atom,) | str/object | — | NO | Sanity check vs the pose batch |
| `metadata` | 0-d object | dict | — | NO | Provenance: {extractor, software, version, …} |

If `forces` is absent or all-zero, the driver auto-switches `--fit-mode
joint → energy`. **The framework is fully usable with energy-only DFT
data** (the common Psi4 / Gaussian / Q-Chem single-point case).

The simplest possible npz that works:

```python
import numpy as np
np.savez(
    "my_teacher.npz",
    pose_params=pose_params,   # (n_pose, 7) float64
    energies=energies,          # (n_pose,)   float64, eV
)
```

---

## 6. The substrate-class abstraction (for developers)

Defined in `pyBall/gridff_jax/substrate_classes.py`. Adding a new class
(e.g. `semiconductor`, `oxide`, `2D_material`) is one dict entry there:

```python
SUBSTRATE_CLASS_DEFAULTS = {
    ...,
    "your_new_class": {
        "description": "<one-line>",
        "builder_mode": "<existing or new builder>",
        "toggles": FeatureToggles(...),
        "density_kind_default": "<vasp_volumetric|cube|surface_xyz>",
        "teacher_tile_default": (NX, NY),
        "use_qeq": ..., "use_image": ..., "use_reactive_grid": ...,
    },
}
```

If your substrate needs new substrate-field physics, add a new
`builder_mode` branch in `pyBall/gridff_jax/substrate_fields.py`'s
`build_substrate_grids()`.

The four existing builder modes:

| `builder_mode` | Used by | Density input | Coulomb physics |
|---|---|---|---|
| `metal_density_plq` | metal class | CHGCAR + LOCPOT | LOCPOT image plane |
| `parity_core` | ionic class | point-charge xyz | Ewald slab |
| `surrogate_density` | (dev/test) | CHGCAR only | Poisson |
| `metal_dft_plq` | (fallback) | density or pairwise | LOCPOT or pairwise |

---

## 7. Plugging in your own MLIP teacher

The `generic_calc` teacher backend wraps any ASE-compatible calculator.
ASE-compatible means it implements `get_potential_energy()` and
`get_forces()` on an `ase.Atoms` object. For most modern MLIPs this is
true out of the box.

Recipes for common MLIPs:

| MLIP | `--teacher-module` | `--teacher-callable` | Notes |
|---|---|---|---|
| MACE | `mace.calculators` | `MACECalculator` | (MAD-SURF is a special case) |
| MACE-MP-0 | `mace.calculators` | `MACECalculator` | Foundation model; supply `--model-path` |
| AIMNet2 | `aimnet.calculators` | `AIMNet2Calculator` | install: `pip install aimnet` |
| TorchANI | (custom wrapper) | (custom wrapper) | needs ~10 line wrapper; ask if blocked |
| NEP | `nep` (calorine pkg) | `Nep` | calorine wraps it as an ASE calc |
| AllergoEX | `allegro` | `AllegroCalculator` | the original Allegro implementation |
| schnetpack | `schnetpack.interfaces.ase_interface` | `SpkCalculator` | needs the spk torch checkpoint |

For anything more exotic, write a 20-line wrapper exposing `__init__(model_paths,
device, default_dtype, …)` and an ASE-compatible interface, then point
`--teacher-module` + `--teacher-callable` at it.

---

## 8. The framework is generic — what we do NOT hardcode

This list exists so you can verify nothing about the framework is
specific to a single substrate/molecule:

- ✅ Substrate identity (metal/ionic chosen at CLI; defaults loaded from
  `substrate_classes.py`).
- ✅ Adsorbate identity (`--adsorbate H2O|CHONH2|<path/to/xyz>` or any
  XYZ file).
- ✅ DFT code (VASP CHGCAR/LOCPOT, Gaussian-cube, or precomputed
  energies in npz — three backends; adding more is one file each).
- ✅ MLIP teacher (any ASE calc via `generic_calc`; specific MACE shim
  exists for backwards compat).
- ✅ Pose grid (n_u, n_v, n_z, n_orient, n_yaw, z range, tilt, random
  fraction, seed — all CLI).
- ✅ Fit mode (energy-only or joint E+F; force weight tunable; outlier
  filter threshold tunable).
- ✅ Output channels (PLQ + reactive + image; toggled by class +
  FeatureToggles).

The only hardcoded assumption is the **physics layer** — Pauli is
exponential repulsion, London is r⁻⁶, Coulomb is 1/r. If you want a
fundamentally different functional form (e.g. anisotropic dispersion),
that's a substrate-field change, not a config change.

---

## 9. Verifying your install

```bash
# 1. Unit tests
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/test_cube_backend.py
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/test_precomputed_teacher.py
# Both should print "All tests PASSED".

# 2. Regression smoke (~30 s on CPU, uses cached Ag dataset)
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/benchmark_substrate_6d.py \
    --substrate-class metal --adsorbate CHONH2 \
    --load-apples tests/tGridFFJax/runs/ag6d_chonh2/apples_dataset.npz \
    --fit-mode joint --force-weight 1.0 \
    --out-dir /tmp/ag6d_regression \
    --device cpu --student-chunk 64
# Train RMSE should be 147.46 meV (bit-identical to the reference run).

# 3. Ionic smoke (synthetic teacher; ~30 s)
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/make_nacl_synthetic_teacher.py \
    --density-xyz cpp/common_resources/xyz/NaCl_1x1_L1_jax.xyz \
    --out /tmp/nacl_synth.npz \
    --n-u 4 --n-v 4 --n-z 8 --n-orient 1 --n-yaw 1
JAX_PLATFORMS=cpu python3 tests/tGridFFJax/benchmark_substrate_6d.py \
    --substrate-class ionic \
    --density-xyz cpp/common_resources/xyz/NaCl_1x1_L1_jax.xyz \
    --adsorbate H2O \
    --teacher-kind precomputed --teacher-npz /tmp/nacl_synth.npz \
    --n-u 4 --n-v 4 --n-z 8 --n-orient 1 --n-yaw 1 \
    --out-dir /tmp/nacl_smoke --device cpu --student-chunk 64
# Should complete in 1–2 min with a sensible Bspline_PLQd.npy export.
```

If all three pass, you can drive any substrate/molecule/teacher combination
with confidence.

---

## 10. Known limitations and design notes

Document these so users and developers don't get surprised:

- **Cube resampling is nearest-neighbour, not trilinear.** `_resample_zyx`
  in `density_backends/cube.py` (and the matching VASP volumetric backend)
  uses integer index arrays from `linspace(...).astype(int)`. This aliases
  the field when you resize the grid; for force-field fitting, aliasing on
  a potential creates artificial corrugation. **Recommendation**: do not
  use `--grid-shape` to coarsen a fine input; resample upstream with a
  proper trilinear/cubic filter, or supply the cube at the resolution you
  actually want.
- **`parity_core` builder is 3D-periodic FFT Coulomb**, not a strict 2D
  slab Ewald. `_poisson_opencl_style` uses `fftn` over x/y/z and drops the
  `k=0` mode. The label "Ewald slab" is a misnomer carried from the
  OpenCL C++ implementation; effectively this assumes full 3D periodicity
  with implicit charge-neutrality. For most lab-scale ionic-slab work
  this is fine (the z-vacuum is large enough that interaction across the
  periodic image of the slab is negligible), but it is NOT the
  textbook Yeh-Berkowitz slab correction.
- **Auto teacher-tiling uses Cartesian extents along lattice-vector
  norms.** `_auto_tile` does `ceil(extent / |a|) + 1`. This ignores the
  reciprocal-cell projection and the MLIP's neighbour-list cutoff. For
  cubic-ish cells (Ag, NaCl in this repo) it is fine; for sharply skewed
  cells or models with long-range cutoffs you may need to set
  `--teacher-tile NX,NY` explicitly to avoid undersized supercells.
- **QEQ has no conditioning checks or charge-bound assertions.** With
  default `--substrate-class metal|ionic` we do not enable QEQ, but if
  you turn it on (`use_ct_qeq=True`), low-hardness / near-plane
  geometries can produce numerically-finite-but-unphysical charges. A
  proper fix adds condition-number logging during fit; for now, treat
  QEQ as experimental.
- **`HybridModelConfig.use_qeq / use_image / use_reactive_grid` are
  informational.** The runtime energy assembly reads `FeatureToggles`.
  These hybrid-model fields are still set by `apply_substrate_class()`
  so downstream exporters and provenance logs reflect intent; they are
  not the source of truth at MD time.
- **`TrainingConfig.batch_size` is currently unused.** The JAX optimiser
  packs the full split into one JIT'd loss. For large 6D datasets that
  can OOM the device; tighten `--max-abs-teacher-eV` and/or shrink the
  pose grid. Field is kept for forward compat.
- **No double-counting audit between PLQ, fixed-charge Coulomb, QEQ,
  image-charge and raw 1/r⁶ channels.** Several of these describe the
  same polarisation/dispersion stabilisation. The code exposes all
  channels but does not enforce orthogonality. When enabling more than
  one of these per-system, run a held-out ablation to confirm each term
  is adding distinct physics rather than absorbing the same residual.

## 11. Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `LookupError: pose not found in npz` | pose-grid args differ between `dft_to_npz_template.py` and `benchmark_substrate_6d.py` | use the SAME `--n-u --n-v --n-z --n-orient --n-yaw --z-min --z-max --seed --random-fraction` for both |
| `metal_density_plq builder requires rho_zyx from CHGCAR` | you chose `--substrate-class metal` but your density backend doesn't supply density | switch to `--density-kind vasp_volumetric` with CHGCAR, or `--density-kind cube` with a density cube |
| `surface_xyz backend requires grid configuration` | density backend called without the GridConfig (very old code path) | now fixed in this repo; pull latest. |
| `Unknown metal '<X>' for C6 grid` warning | the C₆ TS-damped grid only knows the symbols listed in `substrate_fields.py:_ts_c6_pair`; unknown elements give zero dispersion grid | benign for ionic (we use raw 1/r⁶ in the fit instead); add your element to that table for production runs |
| Teacher RMSE > 100 meV | dataset variance is dominated by tilted poses near the Pauli wall (z ≲ 2.2 Å). Filter with `--max-abs-teacher-eV 1.0` (or smaller) and `--z-min 2.2` to focus on the physisorption regime | see `runs/SUMMARY_6D_apples.md` for the full Ag diagnosis |
| Driver auto-switched to fit-mode energy | your npz has no `forces` array; this is correct behaviour | If you DO have forces, store them under the key `forces` in the npz |
| `--substrate-class ionic but no/zero atomic charges` | the substrate xyz column 5 (charges) is missing or all zero | edit the xyz to include physical point charges (e.g. ±0.7 for Na/Cl) — ionic class is forbidden from running on a neutral substrate |
| `FireCore export failed` (run aborts) | export crashed and we now fail-loud (was warn-and-continue) | set `GRIDFF_JAX_EXPORT_LENIENT=1` to downgrade to a warning for debugging; otherwise inspect the traceback — usually a missing grid metadata field |
| `precomputed npz forces have shape X, but the pose batch adsorbate has N atoms` | npz was generated with a different molecule than the one passed to `--adsorbate` | regenerate the npz; verify `adsorbate_symbols` matches; for energy-only data, just omit the `forces` array |

---

## 12. Quick-reference: file layout

```
pyBall/gridff_jax/
├── substrate_classes.py            ← metal / ionic factory (this PR)
├── density_backends/
│   ├── vasp_volumetric.py          ← CHGCAR + LOCPOT
│   ├── cube.py                     ← Gaussian-cube (this PR)
│   ├── surface_xyz.py              ← point-charge xyz
│   └── pseudo_density.py / ml_density.py  ← (R&D paths)
├── teacher_backends/
│   ├── madsurf.py                  ← MAD-SURF MACE (coin-metal only)
│   ├── generic_calc.py             ← any ASE calc (this PR)
│   ├── precomputed.py              ← npz on disk (this PR)
│   └── synthetic.py                ← model-vs-perturbed-self (dev only)
├── substrate_fields.py             ← physics: P/L/Q grid builders
├── hybrid_energy/model.py          ← HybridGridFFModel: channel assembly
├── fit/                            ← Stage 1 PLQ + Stage 2 C₆ lstsq
├── diagnostics/                    ← compare_methods + parity stats
├── pose_sampling/                  ← 6D pose batch builder
└── export/firecore.py              ← Bspline_PLQd.npy writer

tests/tGridFFJax/
├── benchmark_substrate_6d.py       ← single substrate-agnostic driver
├── benchmark_ag_6d.py              ← legacy metal-only driver (kept for ref.)
├── dft_to_npz_template.py          ← turn DFT outputs into the npz schema
├── make_nacl_synthetic_teacher.py  ← synthetic-Coulomb NaCl example
├── compare_nacl_grids.py           ← per-channel grid comparison
├── benchmark_nacl_parity.py        ← OpenCL vs JAX parity (existing)
├── test_cube_backend.py            ← cube reader unit test
├── test_precomputed_teacher.py     ← precomputed teacher unit tests
└── runs/                           ← (gitignored) all output, persists across reboots
    └── SUMMARY_6D_apples.md        ← research results + analysis
```
