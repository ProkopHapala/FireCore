---
name: visual-debugging
description: Use when creating diagnostic plots or visualizations for debugging
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/*debug*.py"
    - "**/test_*.sh"
    - "**/run_*.sh"
    - "**/*benchmark*.py"
---

## Shared Utilities

Before writing ad-hoc debugging/plotting code, check these existing modules:

**Python Visualization:**
- `pyBall/plotUtils.py` - matplotlib utilities: 1D/2D function plots, geometry visualization, field slices, scan profiles, derivative plots
- `pyBall/VispyUtils.py` - GPU-accelerated 3D visualization: AtomScene class for interactive molecular viewing, bond visualization, force vectors

**Python Testing/Diagnostics:**
- `pyBall/DFTB/TestUtils.py` - RMS error computation, checkpoint management (save/load/compare), grid generation, eigenvec printing
- `pyBall/atomicUtils.py` - Atomic utilities: normalize, findAllBonds, graph preprocessing, adjacency lists

**C++ Testing/Diagnostics:**
- `cpp/common/testUtils.h` - Print arrays/vectors/matrices, compareVecs, derivative checking (checkDeriv, checkDeriv3d), timing (StopWatch), error macros (TEST_ERROR_PROC_N, SPEED_TEST_FUNC)

## Test Artifacts

- **Structured outputs:** Group all debugging, benchmarking, and testing outputs into organized, numbered directories (e.g., `tests/003_case_name/`). Do not clutter root directories. Explicitly report their location.
- **Foreground execution:** Run tests synchronously in the foreground with full output. Never hide output or use background commands (`&`, `| tail`, `| head`, or silent redirects). Full `stdout` must be visible.

## Visual Review

- **Python tests:** Generate diagnostic plots using `matplotlib` saved as `.png` files. Use shared helpers like `plotUtils.py` (e.g., `plot_scan_profile`, `plot_field_slice`, `plotGeometryWithForces`).
- **Optional plotting:** Make plotting optional via flags (e.g., `--noPlot`, `--saveFig`). Isolate `plt.show()` strictly to the CLI/main entry point.
- **Report paths:** Always report the exact paths/folders of generated plots.

## Diagnostics

- **Numerical range sanity:** Strategically place checks throughout calculations to ensure values are within reasonable limits and are not `NaN`, infinity, or unexpected zeros.
- **Checkpointing:** Use `pyBall/DFTB/TestUtils.py` checkpoint functions (`save_checkpoint`, `load_checkpoint`, `compare_checkpoint`) for parity testing and reproducible debugging.
- **RMS error:** Use `compute_rms_error` from `TestUtils.py` for array comparisons.

## Consolidation Principle

- **Reuse over reinvent:** Before writing new debug/plot/test functions, search existing utility modules. Generalize existing functions if they almost fit your needs.
- **Separate concerns:** Keep compute algorithms separate from plotting/diagnostics. Move ad-hoc plotting code from test scripts into shared utilities.
- **Zero-copy buffers:** For Python-C++ interop, use `np.ctypeslib.as_array` pattern (see `python_native_bindings` skill) instead of copying data.
