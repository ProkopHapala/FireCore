# GridFF B-spline Fitting: CG vs Damped MD

## Scope
- Implemented faster fitting for 3D B-spline grid potentials using Conjugate Gradient (CG) alongside existing damped MD/Jacobi.
- Applies to Morse-only fitting (no Coulomb/Ewald) in the OpenCL GridFF path.

## Key Implementations
- **Python driver**: `pyBall/OCL/GridFF.py`
  - `fit3D_CG(...)`: Matrix-free CG using gather-only kernels (`BsplineConv3D`) to avoid atomics and write collisions.
  - Stall detection: stops when force norm or residual stops improving for `nconf` steps (mirrors MD `fit3D`).
  - Residual logging: true fitting error `|Ref - A·Gs|` computed via `BsplineConv3D` with coefficients `[1,-1]`.
  - GPU dot product: `dot_gpu` uses `dot_wg` kernel with optional grid-stride loop (can launch fewer threads than `ntot`).
  - Buffers added: `pGs`, `rGs` for CG state.
- **OpenCL kernels**: `cpp/common_resources/cl/GridFF.cl`
  - `dot_wg`: workgroup-local reduction for dot products (no atomics); grid-stride loop; `reqd_work_group_size(64,1,1)`.
  - Existing gather kernels reused: `BsplineConv3D`, `addMul`, `setLinear`, `setMul`.
- **Test harness**: `pyBall/tests/ocl_GridFF_new.py`
  - Job `MorseFit`: fits Pauli/London without electrostatics/FFT; optional `use_CG` flag.
  - Plotting: `plotTrjs` saves convergence figures (`--save-fig`, `--fig-path`).
- **CLI wrapper**: `tests/tMMFF/run_test_GridFF_ocl_new.py`
  - Arguments: `--job`, `--use-CG`, `--nmaxiter`, `--nPerStep`, `--damp`, `--save-fig`, `--fig-path`.
  - Job `CG` aliases `MorseFit` with CG enabled.

## Design Decisions
- **Matrix-free CG**: Use repeated gather convolutions (`A` and `A^T`) to avoid atomics; kernels stay contention-free.
- **True residual reporting**: Residual computed explicitly as `b - A·x`; avoids sign mistakes and reflects fitting error.
- **Stall-based stopping**: Both MD and CG stop after no improvement for `nconf` steps (force and residual separately in CG).
- **GPU dot**: Local reduction avoids atomics; grid-stride loop retained for flexibility, though `nG` is typically `roundup(ntot,64)`.
- **Defaults preserved**: MD behavior unchanged; CG added as opt-in path.

## BsplineConv3D coefficient conventions
- `BsplineConv3D(queue, ..., src0, src1, dst, coeffs)` computes `dst = coeffs[0]*src0 + coeffs[1]*(A·src1)` where `A` is the Bspline convolution operator.
- For the forward apply `A·Gs`, we pass `coeffs = [1, 0]` and `src0 = dGs` (ignored) so the output is purely `A·src1`.
- For the true residual `b - A·Gs`, we call with `src0 = Ref` (b), `src1 = Gs`, and `coeffs = [1, -1]`, so `dst = 1*Ref + (-1)*(A·Gs) = Ref - A·Gs`.
- Sign matters: earlier a flipped sign hid the real fitting error; the current convention matches the intended residual.

## Performance (NaCl_1x1_L3, Morse-only)
- Grid size: 40×40×200 (nxyz=320000). Workgroup 4×4×4 for Bspline conv; dot uses wg=64.
- **CG**: <150 iterations to converge (force ~1e-11, residual ~3.6e-6) — significantly faster and more accurate.
- **Damped MD**: ~3000 iterations (force ~1e-7, residual ~1.4e-6); slower convergence.

## How to Run
```bash
# CG Morse-only fit with saved plot
cd tests/tMMFF
python run_test_GridFF_ocl_new.py --job MorseFit --use-CG 1 --nmaxiter 1000 --nPerStep 50 --save-fig

# MD baseline
python run_test_GridFF_ocl_new.py --job MorseFit --use-CG 0 --nmaxiter 10000 --damp 0.15 --save-fig
```
- Outputs: fitted Bspline coefficients and optional convergence figure under `./data/<name>/`.

## Relevant Paths
- Python CG/MD implementation: `pyBall/OCL/GridFF.py`
- OpenCL kernels: `cpp/common_resources/cl/GridFF.cl`
- Test job and plotting: `pyBall/tests/ocl_GridFF_new.py`
- CLI entrypoint: `tests/tMMFF/run_test_GridFF_ocl_new.py`

## Notes / TODO
- RepeatedKernelRetrieval warnings are benign; could cache kernel handles to silence.
- For larger grids, consider second-stage reduction of `partial` on GPU to avoid CPU summation.
- Extend CG to Coulomb/Ewald path once correctness validated for Morse-only.
