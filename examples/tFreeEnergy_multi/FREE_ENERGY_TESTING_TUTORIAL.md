# Free Energy Testing Tutorial (Entropic Spring)

Document date: **2026-02-26**  
Update this date whenever this workflow is cleaned/refactored.

## 1) What this tests

This folder tests TI/JE/BOTH free-energy workflows for entropic spring systems (`N=10,20,30`) and compares against an analytic reference.

Main scripts:
- `run_ES.py`: single run
- `run_ES_optimization.py`: parameter scan
- `run_ES_params_optimization.sh`: parameter scan wrapper using JSON config (`es_optimization_config.json`)
- `collect_ES_summary.py`: rebuild summary tables/plots from finished runs
- `analyze_free_energy_accuracy.py`: metric calculation against reference

## 2) Recommended workflow (parameter optimization)

Run from `examples/tFreeEnergy_multi`:

```bash
bash run_ES_params_optimization.sh
```

With explicit config:

```bash
bash run_ES_params_optimization.sh es_optimization_config.json
```

Reset previous output before rerun:

```bash
bash run_ES_params_optimization.sh --reset es_optimization_config.json
```

Outputs are written under `bench_ES/` (or `output_root` from config):
- `runs/` per-case logs and data
- `summary/all_runs.csv`
- `summary/best_by_method.csv`
- `summary/best_pareto.csv`
- `summary/plots/summary_dashboard.html`

## 3) Concrete parameter recipes

### A) Fast smoke test (TI only)

```bash
python3 run_ES.py \
  --mode TI \
  --nSys 16 \
  --xyz_name "../tMMFF/data/entropic_spring_10.xyz" \
  --nLambda 30 \
  --nMDsteps 3000 \
  --nEQsteps 300 \
  --Fconv 1e-6 \
  --constraints constraints_ES.txt \
  --K 2.0 \
  --dt 0.05 \
  -T 300 \
  --t_damp 150
```

### B) Balanced optimization batch (TI+JE+BOTH)

```bash
python3 run_ES_optimization.py \
  --N-list 10,20,30 \
  --modes TI,JE,BOTH \
  --nSys-list 64,100 \
  --nLambda-list 80,120 \
  --nEQsteps-list 1000 \
  --nMDsteps-list 100000,200000 \
  --je-k-list 2.0 \
  --output-root bench_ES \
  --dt 0.05 \
  --temperature 300 \
  --t_damp 150 \
  --skip-existing
```

### C) Rebuild summary only (no rerun)

```bash
python3 collect_ES_summary.py --output-root bench_ES
```

## 4) Important: Required conditions:
- Non-bonded interaction is OFF.
- Angle terms are OFF.
- Bond terms are hardcoded in OpenCL kernel to:
  - equilibrium length `1.198`
  - stiffness `40`
  - see `cpp/common_resources/cl/relax_multi.cl:346`
- Debug initializer switch is OFF:
  - `bool initial = false`
  - see `cpp/common/molecular/MolWorld_sp3_multi.h:2016`

If any of these assumptions is changed, results are not directly comparable to current entropic-spring reference benchmarking.

## 5) Quick post-check

After a run, verify:
- data file exists: `entropic_spring_*_free_energy.dat`
- run metrics exist: `runs/.../metrics.json`
- summary was generated: `summary/all_runs.csv`
- dashboard opens: `summary/plots/summary_dashboard.html`
