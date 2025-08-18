# H‑bond Map Plotting and Metrics — Documentation

Script: `tests/tFitREQ/plot_Hbond_ref.py`

This tool loads hydrogen‑bond angle–distance energy maps from multiple methods, visualizes them, extracts minima lines, and computes a soft‑clamp RMSE metric versus a reference. It also sweeps over many systems to summarize accuracy across methods.


## Dataset layout
Under `--root`:
- `data/<method>/*.dat` — evaluated methods
- `ref/<method>/*.dat` — reference DFT methods (e.g., `wb97m`, `b3lyp`)
- Display labels are mapped by `LABELS` (e.g., `wb97m → DFT-wb97m`).

Each `.dat` is triplets: `angle  distance  energy` (comments/blank lines ignored). Grids are regularized internally.


## Quick start
- 2D map (single system):
  ```bash
  python tests/tFitREQ/plot_Hbond_ref.py \
    --methods DFT-wb97m gfn2-xtb pm7 \
    --file 'H2O-D1_CH2O-A1-z' --common-scale
  ```
- 1D minima lines + repulsive filter + RMSE print:
  ```bash
  python tests/tFitREQ/plot_Hbond_ref.py \
    --line-plot --methods DFT-wb97m gfn2-xtb pm7 \
    --file 'H2O-D1_CH2O-A1-z' --min-rmax 2.8 --print-metrics
  ```
- Metric sweep (first 30 systems), show + save:
  ```bash
  python tests/tFitREQ/plot_Hbond_ref.py \
    --metric-plot softclamp --multi 30 \
    --methods DFT-wb97m gfn2-xtb pm7 \
    --clamp-start 0 --clamp-limit 1 --metric-show
  ```
Outputs are saved under `--root/--outdir/` (default `plots/`).


## CLI options
- General I/O
  - `--root PATH` dataset root containing `data/` and `ref/` (default inside script).
  - `--outdir NAME` subfolder for saved figures (default `plots`).
  - `--methods NAMES...` subset of display names to include; else uses `METHODS` order or discovered.
  - `--list` list discovered methods and exit.
  - `--limit N` limit number of methods per figure.
  - `--file STR` select per‑method file by substring or glob.
  - `--file-index I` fallback index if `--file` doesn’t match.
- Rendering
  - `--cmap NAME` colormap (default `seismic`).
  - `--panel W H` per‑panel size in inches (default `3 5`).
  - `--non-uniform` use `pcolormesh` with physical axes; otherwise `imshow` with tick relabels.
  - `--common-scale` force one symmetric vmin/vmax across methods.
  - `--no-share-scale` disable shared scale (deprecated; prefer `--common-scale`).
- Line plots
  - `--line-plot` plot `r_min(angle)` and `E_min(angle)` instead of 2D maps.
  - `--min-rmax R` if `r_min > R`, drop that angle (set to NaN) to ignore repulsive/no‑minimum cases.
- Metrics
  - `--ref-method NAME` reference display name (default `DFT-wb97m`).
  - `--clamp-start y1` soft‑clamp start (default `0.0`).
  - `--clamp-limit y2` soft‑clamp limit (default `1.0`).
  - `--print-metrics` print metrics to console.
- Multi‑system sweep
  - `--multi N` iterate first N files aligned by index across methods; saves figures.
  - `--metric-plot softclamp` (alias `--metric_plot`) compute per‑system soft‑clamp RMSE and plot:
    - `metrics_softclamp_lines.png` — per‑system lines, one line per method.
    - `metrics_softclamp_agg.png` — aggregate mean bar chart across systems.
  - `--metric-show` show the metric figures (in addition to saving).


## How comparisons are normalized
For each method’s 2D grid `g` we apply an asymptotic baseline shift before any comparison:
- Find min of the last distance row: `ref = min(g[last_row,:])` over finite entries.
- Use `gs = g - ref` for plotting and metrics.
This zeroes the far‑distance baseline and makes methods comparable.


## 2D plotting behavior
- `imshow` (uniform pixel spacing) with axis tick relabels to physical values.
- `pcolormesh` if `--non-uniform` for physical axes.
- Color scaling
  - Per‑panel symmetric around 0 using that panel’s min (clipped to ≤ 0).
  - `--common-scale` uses one symmetric range for all methods based on the minimal panel min.


## 1D minima lines
- For each angle column, pick the distance index of the energy minimum → `r_min(angle)` and `E_min(angle)`.
- `--min-rmax R` filters angles where the found minimum sits beyond R (treated as no minimum; set to NaN).
- Metrics in line mode: RMSE vs reference after interpolating onto reference angles.


## Soft‑clamp RMSE metric
Compare shifted grids `gs_model` vs `gs_ref`:
- Compute `diff = gs_model - gs_ref`.
- For `diff > y1`, soft‑clamp towards `y2` with a smooth saturating transform.
- RMSE over finite entries (NaNs ignored). Emphasizes attractive region while limiting penalty for extreme repulsion.


## Multi‑system metric sweep
- Align systems by file index across methods (first N files per method).
- If the reference method is missing for a system, that system’s metrics become NaN for all compared methods.
- Produces a per‑system line plot (x=system labels taken from filenames; y=soft‑clamp RMSE) and an aggregate bar chart of mean RMSE per method.


## Functions (names only, 1–2 lines each)
- `label`
  - Map a folder name (e.g., `wb97m`) to a display label (e.g., `DFT-wb97m`).
- `pick_file`
  - Choose a `.dat` file for a method by substring/glob or by index fallback.
- `_ref_shift`
  - Shift a 2D grid by `min(last_row)` and return the shifted grid and its post‑shift global minimum.
- `load_shifted`
  - Load (angle, distance, grid) from files and apply `_ref_shift()` for all items.
- `softclamp_rmse`
  - Soft‑clamped RMSE between two same‑shape grids using `y1/y2` thresholds; NaNs ignored.
- `rmse`
  - Basic RMSE between arrays with shape check and finite‑mask handling.
- `extract_min_curves`
  - Get `r_min(angle)` and `E_min(angle)` from a shifted grid; optional `rmax` filter sets out‑of‑range minima to NaN.
- `plot_min_lines`
  - Plot `r_min(angle)` and `E_min(angle)` for multiple methods; returns the curves for metrics.
- `align_and_rmse`
  - Interpolate a curve onto the reference angle grid, then compute RMSE.
- `load_energy_map_gnuplot`
  - Parse triplet text into regular grids; handles comments/NaNs robustly.
- `plot_energy_map`
  - Quick single 2D map visualization (standalone helper).
- `plot_all_methods`
  - Multi‑panel 2D visualization with shared/individual scaling; returns loaded/shifted data for metrics.


## Extending
- New metrics (e.g., L1, weighted RMSE) can clone `softclamp_rmse()` and add a `--metric-plot <name>` branch.
- Line processing: add smoothing, angle masking, or custom minima criteria inside `extract_min_curves()`.
- Export metrics to CSV/JSON alongside figures for downstream analysis.


## Notes & limitations
- Multi‑system sweep relies on index alignment across methods; if counts differ, some entries become NaN.
- Interpolation for 1D RMSE assumes angles are monotonic.
- Shared color scale assumes symmetric ranges around zero after baseline shift.
