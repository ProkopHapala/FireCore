# GridFF JAX TODO

## Phase 1: NaCl Numerical Parity
- [ ] Match current OpenCL/GridFF `Q` generation more tightly, including the `makeCoulombEwald_slab()` path, before treating NaCl parity as closed.
- [x] Build a real one-to-one NaCl benchmark for `total`, `morse`, and `coulomb` scans with timings and linear/log error plots.
- [ ] Extend the NaCl benchmark to `pauli` and `london` component scans once the main `PL/Q` gates are stable.
- [ ] Compare JAX-generated NaCl scan outputs against the existing LAMMPS reference workflow used in `tests/tGridFF`.

## Phase 2: Ag(111) Strict PLQ
- [x] Keep Ag phase-1 strict `PLQ` only: fixed charges, rigid slab, no `CT`, no image term, no reactive channels.
- [x] Switch the focused Ag molecular benchmark to explicit `primary/stress/full` windows with fixed charges by default.
  Result: `tests/tGridFFJax/benchmark_ag_zscan.py` now fits only on `2.0-5.6 A`, reports `1.8-2.0 A` stress metrics separately, and saves a `constraint_report` for the final REQ fit.
- [ ] Use `LOCPOT` as the primary Ag `Q` source and keep Poisson-from-`CHGCAR` as a cross-check.
- [x] Rebuild Ag `P/L` from the parity-grade GridFF pairwise substrate-field builder, not the earlier surrogate-density path.
- [x] Run the single-atom probe ladder on Ag for `H`, `C`, and `O` over top/bridge/hollow sites.
- [x] Fit dense strict-`PLQ` molecule benchmarks against MAD-SURF for `H`, `CO`, and `H2O`.
- [x] Test whether a fitted grid-sampling `z` shift materially improves Ag strict-`PLQ`.
  Result: it does not dominate the Ag error; fitted `sample_shift_z` is near zero on atomic probes and only becomes large when the model tries to compensate missing molecule-scale physics.
- [x] Test whether fitting molecule-side static charge scales improves strict-`PLQ`.
  Result: it improves energy RMSE for `H`, `CO`, and `H2O`, but does not materially reduce the large force RMSE on `CO`/`H2O`.
- [x] Fix MAD-SURF interaction-force labeling for rigid polyatomics.
  Result: subtracting the isolated-molecule force rotated into the current pose reduced the strict-`PLQ` aggregate Ag force RMSE from about `2.74 eV/A` to about `0.71 eV/A` on the `654`-pose `H/CO/H2O` benchmark.
- [x] Add the strict-`PLQ` Ag curriculum runner with probe initialization and deterministic orientation sampling.
  Result: Ag fitting now uses a cleaner phase-1 workflow without ad hoc `z`-shift compensation.
- [ ] Calibrate Ag substrate `P/L` builder parameters (`alpha_morse`, `RvdW`, `EvdW`) against the saved `H/C/O` probe dataset.
  Current status: override support, chunk-safe deterministic scan diagnostics, and a real GPU sweep runner exist. The latest deterministic scans show the dominant strict-`PLQ` error is still near-surface lateral corrugation for `H` and `CO`, so this remains the next justified phase-1 improvement.
- [x] Run deterministic fixed-path validation against MAD-SURF for strict `PLQ`.
  Result: full `H/CO/H2O` line, skew-line, rotation, and plane benchmarks now exist with real MAD-SURF labels, raw `*_trace.npz` dumps, linear/log error plots, and per-path timing. `H` and `CO` remain biased positive in energy, especially on the `xy` plane near `z=2.2 A`, while `H2O` is substantially closer to the teacher.
- [x] Run a focused fixed-orientation CHON `z`-scan benchmark on Ag and compare coarse-fit vs dense holdout.
  Result: `tests/tGridFFJax/benchmark_ag_zscan.py` now benchmarks `CHONH2` on `top/bridge/hollow` with only `z` variation. The real `121`-point-per-site run with MAD-SURF on CPU and JAX on `gpu:0` gave aggregate dense-scan errors of about `0.459 eV` and `0.722 eV/A`, while train and holdout errors stayed very similar. That means the remaining error is not primarily coarse-ladder interpolation; it is short-range force-shape error in the Ag strict-`PLQ` model.
- [x] Replace free molecule-side `P/L` amplitudes with a stricter `REQ -> PLQ` molecule coupling on the focused CHON `z` scan.
  Result: the real rerun at `/tmp/ag_zscan_chon_121_reqplq_cpu_teacher` improved the dense CHON benchmark to about `0.411 eV` and `0.715 eV/A`, and improved the low-`z` weighted objective from `8.65` to `8.51`. This confirms that molecule-side radius control matters. The remaining near-surface force error is still present, but it is now less likely to be a pure adsorbate-parameterization problem.
- [x] Run the fixed-charge phase-2 CHON molecular benchmark with deterministic primary/stress windows and a held-out tilted orientation.
  Result: the real run at `/tmp/ag_transfer_phase2_real/CHONH2` reached primary-window `0.0867 eV` and `0.2857 eV/A`, with held-out tilted primary errors `0.0787 eV` and `0.3762 eV/A`. This is a real improvement, but the fit is still `constraint_limited`.
- [x] Add the CHON-to-transfer molecular suite for `CO` and `H2O` under the same strict-`PLQ` protocol.
  Result: `tests/tGridFFJax/run_ag_transfer_suite.py` now runs `CHONH2`, `CO`, and `H2O` with the same builder and windowing. Real results in `/tmp/ag_transfer_phase2_real` show `H2O` transfers very well, `CHONH2` is good but constraint-limited, and `CO` still has a sizeable primary energy residual (`0.487 eV`) despite a good primary force RMSE (`0.256 eV/A`).
- [ ] Resolve the remaining constraint-limited REQ saturation in the phase-2 molecular fits.
  Current status: `CHONH2`, `CO`, and `H2O` all finish with at least one `REQ` parameter within `5%` of its allowed range. `CHONH2` saturates `C/H` radius or `H` energy scale; `CO` saturates `O` radius and both `C/O` energy scales. This is the next molecule-side fitting issue to fix before treating the phase-2 protocol as physically relaxed.
- [x] Run and interpret the H/C/N/O atomic-anchor diagnostics under the same phase-2 protocol.
  Result: the real run at `/tmp/ag_anchor_phase2_real` completed. All four anchors are still `constraint_limited`. Primary-window energy RMSE is about `0.451 eV` for `H`, `0.643 eV` for `C`, `0.413 eV` for `N`, and `0.866 eV` for `O`. This strongly suggests the remaining `CO` energy residual is not a random molecule-only artifact; the current strict-`PLQ` substrate plus molecule-side REQ fit is still too weak in the short-range `C/O` channels.
- [x] Check whether the CHON MAD-SURF teacher path is stable on CUDA for the focused `z` scan.
  Result: the CHON `z`-scan teacher hangs before `teacher_done` on `--device cuda`, even on a tiny `27`-pose benchmark, while the CPU path completes. Treat this as a MAD-SURF/MACE CUDA-path issue for this scan chemistry and keep the focused CHON benchmark on CPU-teacher mode until that path is understood.
- [x] Re-run the coarse Ag substrate rescale sweep after the stricter `REQ -> PLQ` molecule update.
  Result: the focused sweep at `/tmp/ag_zscan_calibration_reqplq_norm` still failed to beat the updated unscaled CHON baseline before it was stopped. Best-so-far objective stayed around `8.86`, worse than the updated baseline `8.51`. This means the remaining CHON error is not fixed by a simple global Ag rescale of `alpha_morse`, `RvdW`, and `EvdW`.

## Phase 3: Residual Analysis
- [ ] Produce residual maps versus MAD-SURF by `z`, site, lateral position, and orientation after the best strict-`PLQ` fit.
- [x] Produce exact path-wise residual diagnostics for fixed lines, rotations, and planes.
  Result: the path benchmark now provides direct teacher-vs-student error on exact `z`, `x`, `y`, `xy`, `xyz`, yaw, tilt, `xy` plane, and `xz` plane scans. This removed the ambiguity from mixed-pose `z` plots.
- [ ] Use those residuals to decide the next missing physics, rather than adding `CT` or image terms by default.
  Current leading hypothesis after the corrected-force and CHON `REQ -> PLQ` runs: the remaining strict-`PLQ` error is concentrated mainly in the short-range Ag substrate field itself. The adsorbate-side factorization is improved, but the near-surface force is still too soft and simple global Ag rescaling is not enough.
