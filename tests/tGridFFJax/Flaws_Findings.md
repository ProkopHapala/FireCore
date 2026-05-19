**Findings**
- **Critical:** `FeatureToggles` is often constructed positionally, and the argument order is wrong. In `plq` mode, this silently sets `use_reactive_grid=True`, so “strict PLQ” runs are not strict PLQ. Same bug appears in [fit_hybrid_gridff.py](/home/niel/git/FireCore/tests/tGridFFJax/fit_hybrid_gridff.py:90), [validate_hybrid_gridff.py](/home/niel/git/FireCore/tests/tGridFFJax/validate_hybrid_gridff.py:73), [analyze_ag_residuals.py](/home/niel/git/FireCore/tests/tGridFFJax/analyze_ag_residuals.py:64), and [run_nacl_parity.py](/home/niel/git/FireCore/tests/tGridFFJax/run_nacl_parity.py:54). Use keyword args only.

- **Critical:** The new “substrate-agnostic” driver still forces a metal recipe onto ionic substrates: London damping, pairwise C6, raw `1/r^6`, REQ-PLQ fitting, and C6 fitting are always enabled in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:295). For NaCl this mixes metal vdW assumptions into the ionic parity path.

- **Critical:** FireCore export does not include the fitted C6 channel or `energy_offset`. The 6D benchmark fits `PLQ + raw 1/r^6 + offset`, but [export_firecore.py](/home/niel/git/FireCore/pyBall/gridff_jax/export/firecore.py:24) only writes `Pauli/London/Coulomb` plus limited sidecars. The exported grid cannot reproduce the reported fitted student if C6/offset matters.

- **High:** Teacher cache invalidation is unsafe. `pose_signature()` hashes only pose/position data, not teacher kind, model path, interaction-energy mode, tiling, density backend, or config. `load_cached_teacher()` then reuses stale labels in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:463). Changing the teacher can silently reuse old data.

- **High:** Cube potential units are wrong. `read_cube()` always converts values as density `e/Bohr^3 -> e/A^3` in [cube.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/cube.py:137), but the same reader is used for potential cubes in [cube.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/cube.py:181). A potential cube should not be divided by volume.

- **High:** `TrainingConfig.batch_size` is effectively dead. The JAX optimizer packs whole splits into one JIT loss in [optimize.py](/home/niel/git/FireCore/pyBall/gridff_jax/fit/optimize.py:414). Large 6D datasets can blow memory, and the config field gives false confidence.

- **High:** Split logic ignores zero fractions. `stratified_split_indices()` uses `max(1, ...)` for val/test in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:371), so `--val-fraction 0 --test-fraction 0` still creates val/test entries and can leave tiny strata with no train samples.

- **Medium:** `GenericASECalcBackend` swallows constructor `TypeError` and retries after dropping kwargs, including `model_paths`, in [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:118). That can instantiate the wrong calculator or hide a real constructor error.

- **Medium:** `teacher_tile` docs say `"auto"` is accepted, but `GenericASECalcBackend.teacher_tile` casts `t[0]`/`t[1]` to int in [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:137). A JSON config with `"auto"` crashes.

- **Medium:** Precomputed teacher force shape is under-validated. If forces exist, output shape follows the NPZ force shape, not the pose adsorbate shape, in [precomputed.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/precomputed.py:160). Missing `adsorbate_symbols` can let wrong force tensors reach fitting.

- **Medium:** `HybridModelConfig.use_qeq/use_image/use_reactive_grid` is mostly config theater. The model energy path reads `FeatureToggles`, not those hybrid flags, except `use_req_plq`. Search shows actual energy switches at [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:503). This invites config drift.

- **Medium:** Hard-coded machine paths are everywhere, e.g. [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:105), [build_ag_dataset.py](/home/niel/git/FireCore/tests/tGridFFJax/build_ag_dataset.py:23), and [prepare_caf2_gridff_jobs.py](/home/niel/git/FireCore/tests/tGridFFJax/prepare_caf2_gridff_jobs.py:10). This is not portable or reproducible off your workstation.

- **Scientific flaw:** The latest 6D summary itself shows the core model is structurally under-expressive: lateral shape is locked to density-derived PLQ and fails badly off-symmetry, with anti-correlation for CHONH2 and weak 2D corrugation for PTCDA. See [SUMMARY_6D_apples.md](/home/niel/git/FireCore/tests/tGridFFJax/runs/SUMMARY_6D_apples.md:159). This is not a tuning problem.

**Fix First**
1. Replace all positional `FeatureToggles(...)` with keyword construction.
2. Make substrate-class defaults own C6/damping/fitting choices too, not just builder/toggles.
3. Extend export or stop claiming C6-fitted runs are FireCore-exportable.
4. Add cache signatures over full config + teacher identity.
5. Fix cube potential units and add a unit test for density cube vs potential cube separately.

---

## Additional Critical Review Findings

These are not repeats of the findings above. They are second-pass issues from the physics, caching, geometry, teacher-backend, and JAX/NumPy paths.

**Failure Path Sketch**

```text
wrong cell/charges/cache metadata
        ↓
teacher labels or student features change silently
        ↓
fit still converges numerically
        ↓
reported meV metrics look precise but describe the wrong physics
```

**More Findings**
- **Critical:** QEQ image energy appears to be counted twice. `solve_qeq_with_reservoir()` adds the image matrix into the interaction matrix before computing `energy` in [qeq.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/qeq.py:67), then separately returns `image_energy` in [qeq.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/qeq.py:98). The model sets `e_qeq = qeq["energy"]`, `e_image = qeq["image_energy"]`, and adds both to the total in [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:491) and [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:540). Same pattern exists in the JAX path at [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:637) and [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:692). This can over-stabilize charged adsorbates near the metal plane and corrupt z-forces.

- **Critical:** `--load-apples` can silently remove molecule partial charges. `save_apples()` stores symbols, positions, and anchor, but not `adsorbate.charges`, in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:460). The reload path rebuilds `AdsorbateDefinition` without charges in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:1020). The energy model uses `poses.adsorbate.charges` to decide whether fixed-charge Coulomb is active in [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:747). Result: re-fitting an apples dataset for H2O, CO, or CHONH2 can change the Hamiltonian even though the cached poses and teacher labels look identical.

- **Critical:** Chunked teacher evaluation drops teacher identity and force-availability metadata. For `n > chunk`, [evaluate_teacher()](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:413) rebuilds a new `TeacherResult` with only timing metadata in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:429). Backend, model path, tile, interaction-energy mode, and `forces_available` are lost. Later, force-mode detection falls back to `np.any(forces != 0.0)` in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:1145), which can misclassify valid zero/symmetric forces as missing forces.

- **High:** Non-orthogonal cells are broken in multiple places. The canonical helper converts Cartesian to fractional with `np.linalg.solve(cell.T, cart.T).T` in [utils.py](/home/niel/git/FireCore/pyBall/gridff_jax/utils.py:66), but `wrap_positions_into_cell()` uses `shifted @ cell_inv.T` and reconstructs with `frac @ cell.T` in [pbc.py](/home/niel/git/FireCore/pyBall/gridff_jax/pbc.py:165). That is invisible for diagonal cells and wrong for skewed cells. Any high-index surface, rotated supercell, or non-orthogonal adsorbate scan can be wrapped to the wrong physical position.

- **High:** `SurfaceXYZBackend` destroys the lattice geometry. It reads `atoms.lvec`, then converts it to diagonal lengths and sets `cell = np.diag(lengths)` in [surface_xyz.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/surface_xyz.py:71). Off-diagonal lattice vectors are discarded. For stepped surfaces, reconstructed slabs, non-cubic NaCl cuts, or rotated supercells, the Coulomb/Morse fields are built in a different cell from the input structure.

- **High:** Ionic charge fields can silently become zero. `SurfaceXYZBackend` takes `charges = np.asarray(atoms.qs)` in [surface_xyz.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/surface_xyz.py:57), but there is no validation that the XYZ actually contains nonzero charges. The resulting `DensityData` is accepted in [surface_xyz.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/surface_xyz.py:91). For an ionic benchmark, missing XYZ charges should be a hard error, not a silent neutral substrate.

- **High:** The parity-core Coulomb path is a 3D periodic FFT potential, not a slab electrostatics model. `_poisson_opencl_style()` uses `fftn` over x/y/z and drops the `k=0` mode in [substrate_fields.py](/home/niel/git/FireCore/pyBall/gridff_jax/substrate_fields.py:109). It is enabled for `builder_mode == "parity_core"` in [substrate_fields.py](/home/niel/git/FireCore/pyBall/gridff_jax/substrate_fields.py:664). That is useful for code parity, but physically it implies full 3D periodicity and hidden neutralization behavior, not 2D slab Ewald/open-z electrostatics.

- **High:** Automatic teacher tiling underestimates image interactions for skewed cells and long-range models. `generic_calc`, `madsurf`, and `pbc` use Cartesian x/y extents divided by cell-vector norms in [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:51), [madsurf.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/madsurf.py:16), and [pbc.py](/home/niel/git/FireCore/pyBall/gridff_jax/pbc.py:440). This ignores reciprocal-cell projection, cutoffs, electrostatics, and skew. The safe criterion should be based on minimum image distance plus calculator cutoff.

- **High:** Gas-phase molecule reference caching is keyed only by adsorbate name. Both [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:192) and [madsurf.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/madsurf.py:145) cache the reference under `poses.adsorbate.name`. Same name but different geometry, charges, protonation, conformer, or anchor will reuse the wrong gas-phase energy and forces.

- **High:** Teacher ASE `Atoms` objects ignore stored adsorbate/substrate charges. The framework carries charges in `AdsorbateDefinition` and `DensityData`, but `GenericASECalcBackend` builds `Atoms(...)` without setting initial charges in [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:199) and [generic_calc.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/generic_calc.py:210). The same pattern exists in [madsurf.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/madsurf.py:152) and [madsurf.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/madsurf.py:210). For charge-aware calculators or DFT wrappers, the teacher may not be evaluating the charged state the student assumes.

- **Medium:** The synthetic teacher does not test the configured framework. It instantiates `HybridGridFFModel(density=density, reactive_elements=elements, prefer_jax=False)` without passing `cfg.toggles`, `cfg.grid`, or `cfg.hybrid_model` in [synthetic.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/synthetic.py:27). Therefore synthetic end-to-end tests can pass while substrate-class settings, JAX execution, QEQ toggles, C6 channels, or interpolation choices are broken.

- **Medium:** Cube resampling is nearest-neighbor while described as trilinear-style. `_resample_zyx()` builds integer index arrays from `linspace(...).astype(int)` in [cube.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/cube.py:150). This aliases densities and potentials during grid-size changes. For force-field fitting, aliasing in a potential grid directly creates artificial force corrugation.

- **Medium:** The cube backend hard-codes a partial periodic table. `_PERIODIC_TABLE` in [cube.py](/home/niel/git/FireCore/pyBall/gridff_jax/density_backends/cube.py:63) omits many elements. This is unnecessary because the project already has element utilities elsewhere. A scientific workflow should not fail on a valid cube just because an element was absent from a handwritten map.

- **Medium:** `fit_params.json` is not a complete model checkpoint. The saved payload includes PLQ, C6, offset, and REQ terms in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:1215), but `HybridParameters` contains additional active parameters such as static charges, QEQ electronegativities/hardnesses, image plane, image damping, sample shifts, and reservoir terms in [model.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/model.py:37). If those are ever fitted or changed from defaults, the run cannot be exactly reconstructed from `fit_params.json`.

- **Medium:** FireCore export failure is swallowed. The main benchmark wraps `export_firecore_artifacts()` in `try/except Exception` and prints only a warning in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:1234). For a workflow whose output is meant to feed FireCore, export failure should fail the run or at least set an explicit failed status in the summary.

- **Medium:** `--prefer-jax` cannot be turned off from the CLI. It is declared as `action="store_true", default=True` in [benchmark_substrate_6d.py](/home/niel/git/FireCore/tests/tGridFFJax/benchmark_substrate_6d.py:222). There is no `--no-prefer-jax`, so users cannot force the NumPy path for debugging numerical parity, CPU-only reproducibility, or isolating JAX compilation issues.

- **Medium:** Precomputed-teacher lookup is exact-row matching with quadratic scaling. Each query pose compares against every cached pose in [precomputed.py](/home/niel/git/FireCore/pyBall/gridff_jax/teacher_backends/precomputed.py:165). This is fragile to floating-point serialization and slow for large 6D grids. Pose IDs or hashed rounded pose keys would be safer.

- **Scientific flaw:** The image-plane QEQ model has no conditioning checks or physical charge bounds. Both NumPy and JAX solve dense linear systems directly in [qeq.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/qeq.py:89) and [qeq.py](/home/niel/git/FireCore/pyBall/gridff_jax/hybrid_energy/qeq.py:187). Low hardness, near-plane geometries, or badly fitted reservoir terms can produce huge unphysical charges while still returning a finite energy. At minimum, log condition numbers and charge extrema during fitting.

- **Scientific flaw:** The current model mixes empirical density-derived PLQ, optional point-charge Coulomb, optional QEQ, optional image charge, and optional raw `1/r^6` without an explicit double-counting audit. Several of these terms can describe the same polarization/dispersion stabilization. The code exposes all channels, but there is no ablation constraint or orthogonality test that proves each term is adding distinct physics.

**Most Dangerous Next Bugs To Fix**
1. QEQ/image double counting.
2. Apples reload dropping adsorbate charges.
3. Full-cache signatures and metadata preservation through chunking.
4. Non-orthogonal-cell handling in wrapping, surface XYZ, and tiling.
5. Hard validation that ionic substrates contain real, nonzero charges.
6. Complete reproducibility payload: fitted parameters, toggles, density metadata, teacher identity, and export status.
