
# GridFF CaF₂ z-scan tutorial & developer notes

This document explains the physics, workflow, CLI, pitfalls, and debugging history behind `run_test_GridFF_CaF2.py`. It is both a user tutorial and a developer reference for the CaF₂ slab PLQ (Pauli–London–Coulomb) GridFF pipeline.

## 1) Physical model (what is computed)
- **Surface PLQ grids (GPU)**: OpenCL kernels (`GridFF.cl`) generate three fields on a 3D grid over the periodic CaF₂ slab:
  - `VPaul` (Pauli repulsion), `VLond` (London dispersion), `VCoul` (electrostatic). They are stored in `Bspline_PLQd_*.npy` with shape `(nx, ny, nz, 3)` in the grid frame.
- **Probe combination (Python)**: For a probe with Morse parameters `(R0, E0, a)` and charge `Q`, `pyBall.tests.utils.getPLQH` computes coefficients `cP, cL, cQ` and multiplies the precomputed grids to form
  - `VNonElec = VPaul * cP + VLond * cL`
  - `VTotal   = VNonElec + VCoul * cQ`
- **Exponent sharing**: The Morse exponent `alpha` must be consistent between the surface build (OpenCL) and the probe combination (Python). The CaF₂ workflow defaults to `alpha=2.0` for both.
- **Prefactor split (surface vs probe)**: From `GridFF.cl` + `getPLQH` the pair prefactors are
  - `cL ∝ E_surface * E_probe * exp(alpha * (Rsurf + Rprobe))`
  - `cP ∝ E_surface * E_probe * exp(2*alpha * (Rsurf + Rprobe))`
  Surface radii enter during grid build; probe radii enter during combination. This is detailed in the CaF₂ debugging note (@doc/Topics/OnSurfaceAssembly/ElectrostaticContinuumEmbeding_report.md lines 393-457).
- **Smearing**: Optional Gaussian smearing `sigma` (Å) applied in reciprocal-space Poisson. Large `sigma` damps short-range Coulomb and helps avoid divergence for attractive pairs.

## 2) Key files and roles
- `tests/tMMFF/run_test_GridFF_CaF2.py`: Main script (z-scans, plots, report).
- `pyBall/tests/ocl_GridFF_new.py`: GPU driver that builds PLQ grids; accepts `alpha_morse`.
- `cpp/common_resources/AtomTypes.dat`: Surface/probe RvdW, EvdW, aliases (`Ca+2`, `F-`, `H`, `O`).
- `cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz`: Slab with charges and lattice vectors.
- `pyBall/tests/utils.py::getPLQH`: Probe PLQH coefficients.
- OpenCL kernels: `cpp/common_resources/cl/GridFF.cl` (`make_MorseFF`, `poissonW`, samplers).

## 3) CLI reference (script arguments)
Run **from `tests/tMMFF`** (script enforces cwd).

- Geometry & grid
  - `--src_xyz` path (default CaF₂ slab)
  - `--dgx --dgy --dgz` grid spacing (Å)
  - `--job` GPU job (PLQ / MorseFit / PLQ_lin / Ewald / Morse)
  - `--save_name` GPU output layout selector
- GPU fit controls
  - `--use_CG` (0/1), `--nmaxiter`, `--nPerStep`, `--damp`
- Probe defaults for combined plots (not scan-specific)
  - `--R0 --E0 --a --Q`
- Smearing & scan combos
  - `--sigma` main run sigma
  - `--sigma_scan` comma list for comparative scans
  - `--probe_specs` comma list `name:q` (e.g., `H:0.2,O:-0.4`)
  - `--probe_a_scan` override alpha for scan probes (default uses `--a`)
  - `--probe_e0_scale` scale EvdW for scan probes
  - `--probe_h_e0 --probe_o_e0` EvdW overrides per probe
  - `--probe_h_r0 --probe_o_r0` RvdW overrides per probe
- Scan window
  - `--scan_zmin --scan_zmax` (Å above target atom)
- Plot windows
  - `--ylim_total min,max` fixed y-lim for total panel
  - `--ylim_comp min,max` fixed y-lim for component panels
  - `--auto_ylim_span` half-span if y-lim not set (default 1 eV)
  - `--auto_total_center` center for total auto-lim (default 0)
  - `--auto_comp_center` center for component auto-lims (default min)
- Slices/diagnostics
  - `--slice_z` z-index for 2D component plots (default auto near top)
  - `--z_above_top` xy cut height above top atom (negative = auto-best-contrast)

## 4) Data flow & outputs

1) **GPU grid build (per sigma)**: `Bspline_PLQd_{sigma}.npy` (shape `(nx,ny,nz,3)`, components ZYX→XYZ transposed) plus XSF slices: `_VPaul.xsf`, `_VLond.xsf`, `_VCoul.xsf` via `export_plq_xsf()`.
2) **2D potential maps** (diagnostics per sigma):
   - `plq_components_gpu_{sigma}.png`: xy component maps at chosen z.
   - `plq_xy_gpu_{sigma}.png`: xy maps for Pauli/London/Coulomb/Total at `z_cut_idx`.
   - `plq_xz_gpu_{sigma}.png`: x–z cuts through cell center.
   - `plq_surface_cuts` (if invoked) combines xy+ xz panels for Pauli/London/Coulomb/Total.
   - `plq_linecuts_gpu_{sigma}.png`: 1D lines along x,y,z through the center.
3) **Z-scans per target/probe/sigma**: `.npz` with (`h`, `VTotal`, `VCoul`, `VPaul`, `VLond`, `VNonElec`, `ix/iy`, `target_z`, probe params).
4) **Comparative z-scan plots**: 4-panel (Total, Coulomb, Pauli, Non-electrostatic) per target/probe across sigmas.
5) **Markdown summary**: `gaussian_smearing_report.md` (one per run) with probe/grid params, prefactors, value ranges, plot paths.

### Directory walkthrough (default path)
Outputs land in `tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/`:
- `Bspline_PLQd_*.npy` — raw GPU PLQ grids per sigma.
- `gpu_sigma_*_VPaul.xsf|VLond.xsf|VCoul.xsf` — XSF slices for visualization.
- `plq_components_gpu_*` / `plq_xy_gpu_*` / `plq_xz_gpu_*` / `plq_linecuts_gpu_*` — 2D and 1D diagnostics.
- `z_scans/` — per-pair `.npz` and comparative PNGs (`zscan_compare_*`).
- `gaussian_smearing_report.md` — textual summary of the run.

### Grid generation workflow (how the grid is built)
1) `run_test_GridFF_CaF2.py` calls `gff_ocl.test_gridFF_ocl(...)` with the slab xyz, grid spacing (`dgx,dgy,dgz`), sigma, and `alpha_morse`.
2) Inside `pyBall/tests/ocl_GridFF_new.py`, the OpenCL solver:
   - Allocates a grid with origin `g0=(-Lx/2,-Ly/2,z0)` and shape matching `dg` and lattice vectors.
   - Runs `make_MorseFF` to accumulate Pauli/London from surface atoms (using surface RvdW/EvdW and `alpha_morse`).
   - Runs `poissonW`/`makeCoulombEwald_slab` to compute Coulomb with optional Gaussian smearing `sigma`.
   - Saves the PLQ grid to `Bspline_PLQd.npy` (transposed to XYZ order) and convergence plots if enabled.
3) Back in the script, per-sigma grids are re-saved with sigma tag, exported to XSF, and plotted as 2D/1D diagnostics before any z-scan extraction.

### Using 2D potentials
- **Visual inspection**: Use `plq_xy_gpu_*` to check lateral structure at a chosen height (auto-selected near top layer if `--z_above_top < 0`).
- **Identify contrast plane**: `choose_reasonable_cut_z` scans z to maximize range in combined PLQ; the chosen `z_cut_idx` is reused for xy/xz/linecuts.
- **Linecuts**: `plq_linecuts_gpu_*` help spot asymmetries or origin mistakes—expect periodicity aligned with lattice vectors.
- **Surface cuts**: `plot_plq_surface_cuts` (callable from code) produces paired xy/xz views for each component and total.

## 5) Known pitfalls and fixes (history-informed)
- **Grid origin**: PLQ grids use centered origin `g0=(-Lx/2,-Ly/2,z0)`. Sampling must subtract `g0` when mapping atom xy → grid indices (`coord_to_index`). Fixed in script; avoid reverting.
- **Alpha consistency**: Use the same `alpha` for surface build and probe combination (`--a` and `--probe_a_scan`). Default is 2.0. Mismatch weakens/over-strengthens Pauli.
- **Weak Pauli**: Original `alpha=1.5` plus small EvdW (Ca ~1.3e-3 eV, F ~2.17e-3 eV) led to collapse for Ca–O with smearing. Raising `alpha` to 2.0 fixes short-range repulsion without inflating EvdW.
- **Aliases**: AtomTypes support `Ca+2` and `F-`. Script tries alias map `{'Ca':'Ca+2','F':'F-'}` for surface lookup; probes also accept aliases. Ensure `AtomTypes.dat` has those rows.
- **Run location**: Must run from `tests/tMMFF` so relative paths to resources and kernel outputs stay valid.
- **Y-limits/centering**: Use `--ylim_total/comp` for fixed ranges or `--auto_total_center/auto_comp_center` to recenter around 0 or minima; prevents plots dominated by far tails.
- **Large EvdW**: Over-inflating EvdW can swamp Coulomb; prefer tuning `alpha` first. If scaling EvdW, note pair prefactor growth is exponential in `(Rsurf+Rprobe)`.

## 6) Quickstart (baseline)
```bash
cd tests/tMMFF
python3 run_test_GridFF_CaF2.py \
  --sigma_scan 0.0,1.0 \
  --probe_specs H:0.2,O:-0.4 \
  --a 2.0 \
  --ylim_total -1,1 --ylim_comp -1,1 \
  --auto_total_center 0.0 --auto_comp_center 0.0
```
Outputs in `tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/` with z-scan plots under `z_scans/`.

## 7) Tuning examples
- **Stronger Pauli only**: `--probe_a_scan 2.5` (keep EvdW default) to raise short-range wall without altering dispersion.
- **Softer Coulomb**: add smearing `--sigma_scan 0.0,0.5,1.0` to compare unsmeared vs smeared.
- **Probe-specific overrides**:
  - Increase H Pauli depth only: `--probe_h_e0 0.003 --probe_e0_scale 1.0`.
  - Larger probe radius for O: `--probe_o_r0 1.9` (shifts wall outward; prefactors updated accordingly).
- **Narrow plot around well**: `--ylim_total -0.5,0.5 --ylim_comp -0.5,0.5`.

## 8) Interpreting the plots
- Panel 1 (Total): Look for correct sign near contact (repulsive) and reasonable zero-crossing height for attractive pairs (e.g., Ca–O crosses ~0.7 Å above Ca at sigma=1.0 in the fixed case).
- Panel 2 (Coulomb): Should reflect probe charge sign; far field decays smoothly. Smearing reduces magnitude near contact.
- Panel 3 (Pauli): Should dominate short range; raising `alpha` steepens the rise.
- Panel 4 (Non-electrostatic = Pauli+London): Shows net short-range balance; London is attractive, Pauli repulsive.
- Use `gaussian_smearing_report.md` to check numeric ranges and prefactors printed by `print_pair_diagnostics()` / `summarize_scan_case()`.

## 9) Developer reference (where things happen)
- **Surface grid build**: `pyBall/tests/ocl_GridFF_new.py::test_gridFF_ocl()` builds PLQ with `alpha_morse`; saves `Bspline_PLQd.npy` (transposed to XYZ).
- **Kernel physics**: `cpp/common_resources/cl/GridFF.cl::make_MorseFF` uses `exp(-alpha*(r - R_surface))`, so surface radius is baked into exponent; Coulomb via `poissonW` with optional Gaussian.
- **Probe coefficients**: `pyBall/tests/utils.py::getPLQH` uses probe `(R0, E0, a)` to generate `cP, cL, cQ` and prints them.
- **Scan extraction**: `extract_vertical_scan` samples GPU grids at `(ix, iy)` and builds `VNonElec`/`VTotal`; `coord_to_index` accounts for centered `g0`.
- **Aliasing**: Surface type lookup tries direct key then alias map; probes require the name to exist in `AtomTypes.dat` (can be aliases).
- **Outputs**: `export_plq_xsf` writes component XSFs for visualization; `.npz` z-scans include metadata (`ix, iy, target_z, probe params`).

## 10) Troubleshooting checklist
- **KeyError for atom type**: Ensure `AtomTypes.dat` has `Ca+2`/`F-`; probe names match; aliases enabled.
- **Nonphysical sign**: Verify run from `tests/tMMFF`; confirm centered indexing (`coord_to_index`) not modified; check probe charge sign in `--probe_specs`.
- **Pauli too weak/strong**: Adjust `--probe_a_scan` first; secondarily scale EvdW via `--probe_e0_scale` or per-probe overrides.
- **Plots off-scale**: Set `--ylim_total/comp` or use auto-centering spans.
- **Smearing effect unclear**: Compare multiple `sigma` values in one run; inspect `gaussian_smearing_report.md` ranges.

## 11) Suggested validation steps
1) Run baseline quickstart (above) and confirm Coulomb sign per pair (F–H attractive mid-range, F–O repulsive, Ca–H repulsive, Ca–O repulsive at contact with zero-cross ~0.7 Å when smeared).
2) Check far-field tail ~0 near largest z; verify zero-cross with `find_zero_crossing` output in logs.
3) If modifying `AtomTypes.dat`, re-run and confirm `print_pair_diagnostics` lines: inspect `pair_cP`, `pair_cL`, and sample potentials at `r0_sum±{0,0.5,1.0}`.

## 12) Reproducibility notes
- Set `OMP_NUM_THREADS` as desired; GPU selection prefers NVIDIA via `pyBall/OCL/OpenCLBase.py` helper.
- Script enforces cwd to avoid stale paths; outputs are versioned by sigma tag.
- Keep `alpha` single-sourced (CLI) to avoid drift between surface build and probe combination.

## 13) Proven fixes (from CaF₂ debugging note)
- Centered grid origin in scans → correct Coulomb sign.
- Exposed `alpha_morse` and defaulted to 2.0 → strengthened Pauli, avoided Ca–O collapse even with sigma=1.0.
- Added y-limit auto-centering options → clearer plots near wells (±1 eV typical).
- Added alias handling for `Ca+2`, `F-` in surface/probe lookup.

## 14) Minimal workflow template
```bash
cd tests/tMMFF
python3 run_test_GridFF_CaF2.py \
  --src_xyz ../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz \
  --sigma 0.0 --sigma_scan 0.0,0.5,1.0 \
  --probe_specs H:0.2,O:-0.4 \
  --a 2.0 \
  --scan_zmin 0.0 --scan_zmax 6.0 \
  --ylim_total -1,1 --ylim_comp -1,1 \
  --auto_total_center 0.0 --auto_comp_center 0.0
```
After the run, open `data/<slab>/z_scans/zscan_compare_*` PNGs and `gaussian_smearing_report.md` for numeric ranges.
