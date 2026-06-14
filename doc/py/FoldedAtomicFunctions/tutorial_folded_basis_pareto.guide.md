---
"description": "Folded basis Pareto search tutorial"
---

# Folded Atomic Functions: Pareto Search Tutorial

This tutorial explains the small toolkit in `doc/py/FoldedAtomicFunctions/` for exploring and fitting compressed 1D z-basis sets (cutoff-polynomial / composite bases) used in folded surface potentials. It covers the scripts, their inputs/outputs, where results live, and the math behind the scoring.

## Files and roles

- **`Optimize_basis_batch.py`** – batch driver that sweeps many basis definitions, fits to sample potentials, and appends results to a JSONL DB.
- **`Optimize_basis.py`** – single-run optimizer for a specific basis config (handy for quick tests/debugging).
- **`plot_pareto_basis.py`** – reads a results DB, computes the Pareto front (RMSE vs nBasis), and plots both the scatter and the basis/fits for each Pareto point.
- **`results_db.py`** – thin JSONL-backed key/value store; includes `pareto_front()` helper.
- **`basis_utils.py`** – basis generators and utilities (polynomial, cutoff-polynomial, composite bases; metrics like nBasis / max power).
- **`plot_db_utils.py`** – helpers to decode basis keys (`basis_key_to_struct`) and to plot basis definitions from DB entries.
- **`results.jsonl` (or similar)** – JSON lines database produced by batch runs. Each line is a record keyed by a string encoding `(z0, [(z_cut_i, [powers_i]) ...])`.
- **`Morse*.json` sample sets** – reference potentials (z grids + samples) that the batch optimizer fits against.

## Quickstart: run a batch and plot the Pareto front

```bash
cd doc/py/FoldedAtomicFunctions
# Run a batch sweep (writes/updates results.jsonl by default)
python Optimize_basis_batch.py --samples Morse1.json --results results.jsonl

# Plot Pareto front and per-basis fits from the DB
python plot_pareto_basis.py --db results.jsonl --samples Morse1.json
# Outputs: pareto_scatter.png, pf_<idx>_* images in the current directory
```

Typical outputs:
- `results.jsonl` (appended/updated DB)
- `pareto_scatter.png` (all points in grey, Pareto front in red)
- `pf_<idx>_nbas<N>_basis.png` (basis function shapes)
- `pf_<idx>_nbas<N>_fits.png` or `_basis_fits.png` (recon vs reference samples)

## Inputs

- **Sample file** (`Morse1.json`, etc.): contains `z_grid`, `Y_samples` (reference potentials), optional weights, and Morse params for labeling. Loaded by `load_morse_samples_json()`.
- **Basis definitions**: encoded as lists of `(z_cut, [powers])`, e.g. `[(4.8, [2,5,7,10,13]), (8.6, [2,4,7])]`. Stored as string keys in the DB: `"z0,[(z_cut,[powers]), ...]"`.
- **Config knobs (batch)**: degree ranges, number of functions, z-cut ranges, seeds; see CLI flags in `Optimize_basis_batch.py`.

## Outputs

- **DB records (JSONL)**: one line per basis key, containing at least `rmse`, `timestamp`, `samples` name, and sometimes extra metrics. Example (from `results copy.jsonl`):
  ```json
  {"1.0,[(4.8, [2, 5, 7, 10, 13]), (8.6, [2, 4, 7])]": {"rmse": 1.9075e-05, "timestamp": "2025-06-18T15:14:18Z", "samples": "Morse1.json"}}
  ```
- **Pareto plots**: `pareto_scatter.png` plus per-Pareto-point basis/fits PNGs.
- **Basis/fits visuals**: overlaid reconstructions vs. reference samples with coefficient labels.

## How the Pareto front is computed

- Metric: \((x, y) = (n_{\text{basis}}, \text{RMSE})\).
- `ResultsDB.pareto_front(metric_fn, unique_x=True)` sorts by \(x\) and keeps records where \(y\) is strictly improving (monotone lower envelope). `unique_x=True` keeps only the best RMSE per nBasis bucket before the sweep.

**RMSE definition** (per sample set):
\[
\text{RMSE} = \sqrt{\frac{1}{M}\sum_{j=1}^M w_j\,\|\hat{V}_j - V_j\|_2^2}
\]
where \(V_j(z)\) is the reference potential sample, \(\hat{V}_j(z)\) is the basis reconstruction, and \(w_j\) are optional per-sample weights.

**Reconstruction**:
Given basis matrix \(\Phi \in \mathbb{R}^{P \times N_z}\) and sample vector \(V \in \mathbb{R}^{N_z}\), coefficients \(c \in \mathbb{R}^P\) are found by weighted least squares
\[
\min_c \|W (\Phi^T c - V)\|_2^2 \quad\Rightarrow\quad c = \arg\min_c \|W\Phi^T c - W V\|_2^2
\]
Then \(\hat{V} = \Phi^T c\). In the scripts this is done column-wise for multiple samples: `coeffs = fit_coefficients(Y_samples, Phi, weights)`.

## Basis parameterization (cutoff-polynomial composite)

A basis is a list of blocks \((z_{\text{cut}}, [p_1, p_2, \dots])\). For a given block:
- Define \(h(z) = \max(0, z_{\text{cut}} - z)^2\) (optionally on scaled z).
- Basis rows: \(h(z)^{p_k}\) for each power \(p_k\) in the list.
- Composite basis concatenates all blocks’ rows.

`basis_utils.construct_composite_cutoff_basis(z_grid, basis_def, z0)` builds \(\Phi\) and labels; `get_basis_metrics(basis_def)` returns `(nBasis, nZcut, nMaxPow, sumPow)`.

## Typical workflow (end-to-end)

1) **Prepare samples**: ensure `Morse*.json` is present (or your own sample set in the same format).
2) **Run batch search**:
   ```bash
   python Optimize_basis_batch.py --samples Morse1.json --results results.jsonl
   ```
   This sweeps many basis candidates, appending RMSE entries to `results.jsonl`.
3) **Inspect Pareto front**:
   ```bash
   python plot_pareto_basis.py --db results.jsonl --samples Morse1.json
   ```
   Check `pareto_scatter.png` and the `pf_*` images to see the best trade-offs.
4) **Pick a basis**: choose a Pareto point balancing nBasis and RMSE; the key string encodes \(z_0\) and the block structure.
5) **Reuse elsewhere**: plug the chosen basis definition into folded-surface fitting code (e.g., pass it to any caller that accepts `basis_def` or convert it into kernel parameters).

## Where results live

- Default DB: `doc/py/FoldedAtomicFunctions/results.jsonl` (or `results copy.jsonl` if you duplicated it).
- Plots: written to the current working directory where you ran the script (often the same folder).
- Sample sets: `Morse1.json` (extend with your own if needed).

## Tips and pitfalls

- Run from the folder so relative imports work: `cd doc/py/FoldedAtomicFunctions`.
- Use `unique_x=True` (default in `plot_pareto_basis.py`) to get one Pareto point per nBasis bucket.
- Large batches can be slow; constrain degree/cut ranges to a few hundred candidates for quick iterations.
- Keep `results.jsonl` under version control only if you want the historical sweep; otherwise treat it as scratch.
- Always inspect `pareto_scatter.png` plus the `pf_*` basis/fits to avoid overfitting weird shapes.
