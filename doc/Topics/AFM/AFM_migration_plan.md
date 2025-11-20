# AFM PyOpenCL Migration Plan (2025 Update)

This plan describes how to replace the `libOCL_GridFF.so` C/C++ host layer with a pure Python implementation that relies on NumPy for data preparation and PyOpenCL for GPU execution. It expands previous drafts by synchronizing with the verified documentation in @/doc/Topics/AFM/AFM.md and enumerating every C/C++ function, class, and data dependency that must be re-created in Python.

---

## 1. Scope and Goals

- **Goal:** deliver a Python package (working name `pyBall.pyocl_dft`) that exposes the same API currently provided by `pyBall.DFT.oclfft`, while driving the existing OpenCL kernels (`myprog.cl`, `GridFF.cl`, `relax.cl`) directly from Python.
- **Out of scope:** rewriting the OpenCL kernels themselves; re-deriving Fireball DFT routines; changing physics models.
- **Success criteria:**
  1. All Python tests under `tests/tDFT*` and AFM tutorials run using the new package without loading `libOCL_GridFF.so`.
  2. Density grids, potentials, relaxation trajectories, and frequency-shift outputs match C++ results within numerical tolerance (define per dataset in §6).
  3. Documentation and examples default to the Python host while keeping the legacy C++ path as fallback until feature parity is proven.

---

## 2. Summary of the Current Workflow

The AFM workflow has five computational stages (@/doc/Topics/AFM/AFM.md#37-112):

1. **Fireball DFT** (Fortran) produces wavefunction coefficients `wfcoef` and density matrices.
2. **Density projection** uses `projectAtomsDens`, `projectAtomsDens0`, or `projectDenmat` to populate FFT-aligned grids (`float2` buffers) via `projectOrbDenToGrid_texture`, `projectAtomDenToGrid_texture`, or `projectDenmatToGrid` kernels (@/doc/Topics/AFM/AFM.md#41-54).
3. **Potential assembly** performs FFT-based Hartree/Pauli/VdW operations (`poisson`, `convolve`, `evalLJC_QZs`) (@/doc/Topics/AFM/AFM.md#55-90).
4. **Probe relaxation** runs `relaxStrokesTilted` to locate equilibrium tip positions (@/doc/Topics/AFM/AFM.md#88-100).
5. **Frequency-shift post-processing** applies `convolveZ` weights to relaxed force traces (@/doc/Topics/AFM/AFM.md#103-111).

Each stage is currently orchestrated in C++ via `OCL_DFT` and `OCL_PP` classes, exported in `OCL_GridFF.cpp` and invoked through `pyBall.DFT.oclfft.py` (summarized in Appendix A of AFM.md @/doc/Topics/AFM/AFM.md#214-233). The migration replicates those orchestrations in Python.

---

## 3. Target Python Package Layout

Create a new directory `pyBall/pyocl_dft/` (parallel to `pyBall/DFT/`) with the following modules:

| Module | Responsibility | C++ Source to Replace |
|---|---|---|
| `context.py` | Device/context management, kernel compilation, buffer registry | `OCL_GridFF.cpp` (init/release), `OCL_DFT::makeMyKernels`, `OCL_PP::makeKernels_PP` |
| `fft.py` | clFFT/pyclfft planning, FFT/IFFT helpers, buffer allocation | `OCL_GridFF.cpp::initFFT`, `OCL_DFT::runFFT`, buffer bookkeeping |
| `io.py` | Loads/saves (`saveToBin`, `loadWfBasis`, `saveToXsf*`) using NumPy | `OCL_GridFF.cpp`, `OCL_DFT::loadWfBasis` |
| `density.py` | Implements `project_dens_GPU`, `project_denmat_GPU`, `projectAtomsDens0` using PyOpenCL | `OCL_DFT::projectAtomsDens`, `projectAtomsDens0`, `projectDenmat`, helper preparation routines |
| `potentials.py` | Wraps `poisson`, `convolve`, `gradient`, `evalLJC_QZs` | `OCL_GridFF.cpp` functions & `OCL_PP::evalLJC_QZs` |
| `relax.py` | Probe-particle relaxation, `getFEinStrokes`, parameter utilities | `OCL_PP::relaxStrokesTilted`, `getFEinStrokes` |
| `df.py` | Implements Giessibl frequency-shift convolution using `convolveZ` | `relax.cl::convolveZ` orchestration (currently C++ utilities) |
| `assets.py` | Discovers `Fdata/basis/`, validates grid descriptors, manages `acumCoef` presets | New |
| `tests/` (package) | Python-level regression harness mirroring `run.sh` scripts | New |

Provide a top-level facade `pyBall/pyocl_dft/__init__.py` exposing functions with the same names/signatures as `pyBall.DFT.oclfft` so existing user scripts import the new module transparently.

---

## 4. Detailed Migration Tasks

### 4.1 Stage 0 – Preparatory Work

1. **Dependency inventory**
   - Confirm availability of `pyopencl`, `pyclfft` (or alternative FFT bindings). Pin versions in `requirements.txt` or project docs.
   - Audit `Fdata/basis/` contents; document mandatory elements (per AFM.md §2.1 @/doc/Topics/AFM/AFM.md#115-123).
2. **Kernel packaging**
   - Copy `myprog.cl`, `GridFF.cl`, `relax.cl` to a runtime-accessible directory; ensure `context.py` loads them relative to repository root or an environment variable to mimic `cl_src_dir` semantics.
3. **API compatibility layer**
   - Design decorators to match current `ctypes` signatures (e.g., `poisson(ibuff_in:int, ibuff_out:int, dcell)`), raising `NotImplementedError` until each feature is ported.

### 4.2 Stage 1 – Context & Resource Management (C++: `OCL_GridFF.cpp`, `OCL_DFT` core)

| Task | Python Target | Notes |
|---|---|---|
| Context initialisation | `context.py:init()` | Create PyOpenCL context/queue, compile kernels, mirror `OCL_DFT::makeMyKernels()`; manage logging/verbosity. |
| Buffer bookkeeping | `context.BufferTable` class | Store metadata (`shape`, `dtype`, `format`) for each logical buffer; replace global arrays `buffers[]`, `Ns[]`. |
| Upload/download helpers | `context.upload`, `context.download` | Support float/float2/double conversions (cf. `upload_d`) and ensure host/device offsets follow memory rules (check retrieved memory `e22a0816...`). |
| Texture/image handling | `context.create_image3d` etc. | Mirror `newFFTimage`, respecting CL image formats required by `evalLJC_QZs_toImg`. |
| Grid descriptors | `assets.make_grid_descriptor` | Provide `pos0`, `dA/B/C`, ensure compatibility with kernels (per AFM.md supporting assets). |

### 4.3 Stage 2 – FFT Infrastructure (C++: `initFFT`, `runFFT`, `convolution` scaffolding)

1. Implement `fft.PlanCache` to encapsulate clFFT plan creation keyed by grid dimensions.
2. Provide `fft.run_fft(buffer_id, direction)` to enqueue forward/inverse FFTs, matching the existing three-buffer setup (`inputA`, `inputB`, `outputC`).
3. Add helpers to resize FFT buffers when grid dimensions change (mimic `newFFTbuffer` logic).

### 4.4 Stage 3 – Density Projection Pipeline (C++: `OCL_DFT` methods)

| Python Function | Replaces | Key Steps | Reference |
|---|---|---|---|
| `density.project_atoms_dens()` | `OCL_DFT::projectAtomsDens` | 1) prepare `float4` atom arrays & coefficient matrices; 2) set `acumCoef=[0.0, 2.0]` (default) or alternative; 3) call `projectOrbDenToGrid_texture`. | AFM.md Step 2 @/doc/Topics/AFM/AFM.md#41-54 |
| `density.project_atoms_dens0()` | `OCL_DFT::projectAtomsDens0` | Upload neutral atom templates, set `acumCoef=[1.0, -1.0]` for difference densities. | AFM.md Step 2 bullet on difference density |
| `density.project_denmat()` | `OCL_DFT::projectDenmat` | Re-implement `atoms2box` chunking in Python (likely using NumPy slicing) before calling `projectDenmatToGrid`. | AFM.md Step 2 alternative workflow |
| `density.prepare_basis()` | `OCL_DFT::loadWfBasis`, `convCoefs`, `assignAtomDensCoefs` | Parse `.wf1/.wf2`, perform 1D interpolation (NumPy/scipy), reorder orbitals, compute accumulators. | Appendix A, Supporting Assets |

Additionally, mirror supportive structures:

- **Coefficient conversion**: reproduce the `convCoefs` logic for s/p decomposition; include assertions for expected orbital counts to "fail loudly" per user rules.
- **Atom boxing**: implement `atoms2box` in Python for chunk-wise processing; optional but needed for large systems.

### 4.5 Stage 4 – Potential & Energy Kernels (C++: `OCL_GridFF.cpp` functions)

| Python Function | Responsibilities | Replacement Details |
|---|---|---|
| `potentials.poisson()` | FFT density, apply `poissonW`, inverse FFT | Follow C++ sequence; parameterize `dcell` (float4) handling. |
| `potentials.convolve()` | FFT two buffers, call `mul`, inverse FFT | Ensure optional scaling factors match legacy behavior. |
| `potentials.gradient()` | Launch `gradient` kernel | Provide Python wrapper for gradient magnitude vector fields. |
| `potentials.eval_ljc_qzs()` | Launch `evalLJC_QZs` | Accept atom arrays, LJ parameters, charges; handle `float4` output (xyz force + energy). |
| `potentials.eval_ljc_qzs_to_img()` | Launch image variant | Create or reuse 3D image; ensure sampler setup matches existing usage. |

If the project chooses to expose `make_MorseFF` (currently only through PyOpenCL demos), provide optional wrapper using the same infrastructure.

### 4.6 Stage 5 – Probe Relaxation & Frequency Shift (C++: `OCL_PP`)

1. `relax.relax_strokes_tilted()` – replicate FIRE-based iteration (`update_FIRE` logic lives in kernel) by setting kernel arguments, including scan grid, step sizes, damping parameters.
2. `relax.get_fe_in_strokes()` – capture force/energy along predefined trajectories without relaxation.
3. `df.convolve_z()` – wrap the `convolveZ` kernel to produce frequency-shift values using weight arrays; provide convenience functions to generate Giessibl weights.
4. Provide Python utilities to configure scan grids mirroring `makeStartPointGrid` and `setGridShapePP`.

### 4.7 Stage 6 – I/O and Utilities

- Recreate `saveToBin`, `loadFromBin`, `saveToXsf`, and friends using NumPy and shared helper functions (`io.py`).
- Implement logging/verbosity toggles matching `setVerbosity`, `setErrorCheck`.
- Ensure Python functions reuse the "Supporting Assets" rules from AFM.md (basis availability, buffer formats, etc.).

---

## 5. Validation & Testing Strategy

1. **Unit tests** for each Python wrapper using small synthetic systems (e.g., single hydrogen atom) to confirm kernel argument wiring.
2. **Integration tests** replicating `tests/tDFT_pentacene/run.py` and related harnesses. Update each directory’s `run.sh` to allow `PYTHON_AFMPATH=pyBall.pyocl_dft` (feature flag) to switch between implementations.
3. **Numerical parity checks**:
   - Compare density grid norms and pointwise differences after `project_atoms_dens`.
   - Compare Hartree potentials (`poisson`), VdW grids (`eval_ljc_qzs`), and relaxed positions.
   - Validate `convolveZ` outputs against saved reference traces.
4. **Performance benchmarking** to ensure PyOpenCL overhead is acceptable; gather timings similar to existing `GridFF_cl::make_MorseFF()` prints.

---

## 6. Risks and Mitigation

| Risk | Impact | Mitigation |
|---|---|---|
| Missing basis data | Python path fails to load orbital textures | Add pre-flight checks in `assets.py`; provide descriptive error messages referencing @/doc/Topics/AFM/AFM.md#115-123. |
| FFT API differences | Incorrect scaling/orientation | Cross-check with C++ output; wrap clFFT plan creation to mimic stride/order; add regression tests. |
| Memory alignment | Kernel reads wrong offsets | Mirror buffer-offset logic noted in internal memory reminders (host/device offsets) and validate with assertions before kernel launch. |
| Numeric drift | Breaks regression comparisons | Keep computations in `float32` to match GPU; where double precision is required, document conversions. |

---

## 6.1 Critical Implementation Notes

### Reuse of `OpenCLBase`

All new PyOpenCL solver classes introduced by this migration must inherit from `pyBall.OCL.OpenCLBase.OpenCLBase`. The base class already encapsulates:

- Device selection (`select_device`) and queue creation with logging hooks.
- Kernel compilation (`load_program`) plus header parsing utilities used across existing OpenCL components.
- Buffer bookkeeping (`create_buffer`, `check_buf`, `toGPU`, `fromGPU`).

Adopting it ensures the new Python host remains aligned with current GPU tooling (e.g., atomic MMFF solvers). When implementing modules such as `pyocl_dft/context.py`, wrap the functionality in subclasses that extend `OpenCLBase` rather than rebuilding context logic from scratch.

### FFT Initialisation Patterns

The project already integrates clFFT/pyclfft in several locations:

- `tests/tMMFF/run_test_GridFF_ocl.py`
- `pyBall/tests/ocl_GridFF.py`
- `pyBall/OCL/GridFF.py`
- `pyBall/OCL/clUtils.py`

These files demonstrate correct plan creation, buffer registration, and execution order (FFT → kernel → inverse FFT) using the shared utility functions. The PyOpenCL AFM migration should follow the same conventions—initialise FFT plans via the helper routines in `clUtils.py` (or refactor shared helpers if necessary) and reuse the `GridFF_cl` approach for organising kernels and buffers. Documenting this reuse avoids subtle bugs such as incorrect strides or plan lifetime mismatches.

### Fireball DFT Integration

The Python AFM host continues to obtain molecular orbital data from the Fireball Fortran library via `pyBall.FireCore`. Critical entry points:

- **Initialisation sequence** (`FireCore.initialize`, `FireCore.preinit`, `FireCore.init`): sets verbosity, prepares species tables, and loads atomic coordinates. Example usage is shown in `tests/pyFireball/relax_molecules.py`.
- **SCF driver** (`FireCore.SCF`, `FireCore.assembleH`, `FireCore.solveH`, `FireCore.updateCharges`): run after `initialize` to populate the wavefunction coefficients (`wfcoef`) and density matrix.
- **Coefficient access**:
  - `FireCore.get_wfcoef(ikp, wfcoef)` retrieves the full molecular-orbital coefficient matrix as a NumPy array (shape `norbitals × norbitals`).
  - `FireCore.getPointer_wfcoef` exposes the raw pointer when zero-copy access is required (ensure lifetime management is understood before using).
- **Density exports**:
  - `FireCore.setupGrid(Ecut, g0, ngrid, dCell)` configures the real-space grid used by Fireball’s own projection routines and returns `(ngrid, dCell, lvs)` descriptors that should be forwarded to the OpenCL host.
  - `FireCore.getGridDens(ngrid=..., Cden=1.0, Cden0=0.0)` outputs an electron density grid directly from Fireball if we want to bypass the GPU projection for validation.
  - `FireCore.dens2points(points, f_den=1.0, f_den0=0.0)` evaluates the density at arbitrary coordinates, useful for spot checks.
  - **Grid sizing & ordering:** Fireball’s `setupGrid` will snap the requested lattice to FFT-friendly dimensions (e.g. returning `48×52×52` instead of the nominal `50×50×50`). Always persist the exact `ngrid` reported by the Fortran driver and forward it unchanged to the GPU pipeline. The Fortran buffers are column-major; when exposing them to NumPy reverse the shape (`ngrid[::-1]`) or allocate Fortran-contiguous arrays so that `z` remains the fastest axis when plotting or projecting.
- **Reference scripts**: explore `tests/pyFireball/relax_molecules.py` for end-to-end initialisation and relaxation, and `pyBall/FireCore.py` for the authoritative list of bindings to `fortran/MAIN/libFireCore.f90`.

During the migration, design the PyOpenCL pipeline so that the Fireball SCF step can hand off either (a) the MO coefficient matrix for GPU projection (`projectOrbDenToGrid_texture`) or (b) the already-sampled density grid for comparison. Record grid metadata (`pos0`, `dA/B/C`, `ngrid`) immediately after `setupGrid` so both paths remain consistent.

---

## 7. Roadmap & Milestones

1. **Milestone A – Fireball Data Extraction**
   - Verify Fireball SCF workflow via `pyBall.FireCore` to obtain wavefunction coefficients (`get_wfcoef`) and real-space density grids (`getGridDens`).
   - Capture grid descriptors from `setupGrid` and persist both the density matrix (2D NumPy array) and reference density (3D NumPy array) to disk for GPU comparison.
   - Provide `tests/pyocl_dft/test_firecore_data.py` CLI script that loads a sample molecule (from `tests/pyFireball`), runs the Fireball calculation, saves the coefficient/density files, and plots diagnostic slices using `matplotlib`.
2. **Milestone B – Skeleton & Context**
   - Implement `context.py`, `fft.py`, stub wrappers.
   - Validate kernel compilation and simple FFT round-trip.
   - Produce a standalone smoke-test script (now superseded by milestone A artefacts).
3. **Milestone C – Density Projection**
   - Port basis loading, coefficient conversion, `project_atoms_dens`.
   - Achieve parity on simple systems; add unit tests.
   - Deliver `tests/pyocl_dft/test_density.py` CLI script that loads Fireball-produced coefficient/density files, runs the Python host projection, saves the complex grid to `.bin`/`.xsf`, and plots a density slice using `matplotlib`.
4. **Milestone D – Potentials**
   - Implement FFT-based Hartree/Pauli/VdW wrappers.
   - Validate against C++ outputs for pentacene/CO tests.
   - Provide `tests/pyocl_dft/test_potentials.py` CLI script that loads density files, performs Hartree/Pauli convolutions, saves resulting grids, and plots radial profiles with `matplotlib`.
5. **Milestone E – Relaxation & df**
   - Port probe relaxation, stroke utilities, `convolveZ`.
   - End-to-end AFM image produced via Python host.
   - Ship `tests/pyocl_dft/test_relaxation.py` CLI script that relaxes a minimal scan grid, stores trajectories, and plots force traces/frequency shifts using `matplotlib`.
6. **Milestone F – Integration & Documentation**
   - Update tutorials, ensure `run.sh` scripts can switch hosts.
   - Collect benchmarks, document known differences.
   - Create `tests/pyocl_dft/test_end_to_end.py` CLI script orchestrating the full pipeline (density → potentials → relaxation → df), writing summary artifacts and generating figure panels for documentation.

---

## 8. Deliverables Checklist

- [ ] `pyBall/pyocl_dft/` package with modules described in §3.
- [ ] Unit and integration tests passing with Python host enabled.
- [ ] Updated documentation referencing the new Python path (AFM.md already compatible).
- [ ] Migration guide for developers (how to switch, debug, extend).
- [ ] Performance comparison report (Python vs C++ host).

Once all boxes are checked and parity is confirmed, the project can decide whether to deprecate `libOCL_GridFF.so` or keep it as a fallback for legacy workflows.

