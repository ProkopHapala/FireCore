# Si ~1 nm passivation examples

Small silicon nanocrystals (≈50–130 atoms) with **different surface terminations** — full **01→05 NPZ pipeline** outputs for viewer QA and vibration benchmarks.

**Contracts:** [`Nanocrystal_NPZ_Pipeline.guide.md`](../../../../doc/Topics/FTIR_Nanocrystals/Nanocrystal_NPZ_Pipeline.guide.md) · [`NPZ_Crystal_Schema.md`](../../../../doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md)

Each subfolder contains:

| File | Stage | Role |
|------|-------|------|
| `01_crystal.npz` | Generate | Init `pos`, `Z`, `bonds_ij` |
| `03_topology.npz` | Generate | Fixed MMFFL springs + surface AABBs |
| `meta.json` | Generate | Defect counts, gen params |
| `02_relaxed.npz` | Relax | Minimized `pos`; same `Z`/`bonds_ij` |
| `04_hessian.npz` | Hessian | `K`, `K_projected`, `M` at relaxed coords |
| `05_spectrum.npz` | Spectrum | Mode frequencies + histogram |
| `spectrum.png` | Spectrum | Line plot |
| `*.status.json` | — | Timing / diagnostics |

## Gallery

| Folder | Atoms | Passivation / surface |
|--------|------:|------------------------|
| `01_sphere_11A_H_standard` | 68 | 11 Å sphere — H monolayer (SiH + SiH2), no defects |
| `02_sphere_13A_H_standard` | 118 | 13 Å sphere — H monolayer, no defects |
| `03_sphere_13A_H_bridge_insert` | 133 | 13 Å sphere — H + 5 inserted `-SiH2-` bridges |
| `04_sphere_13A_H_bridge_collapse` | 106 | 13 Å sphere — H + 4 collapsed bridges (fused rings) |
| `05_sphere_13A_H_mixed_defects` | 124 | 13 Å sphere — H + insert + collapse |
| `06_sphere_13A_bare_surface` | 60 | 13 Å sphere — **no H caps** (bare dangling Si) |
| `07_faceted_111_H` | 125 | {111} faceted G1 cut — mostly SiH-like |
| `08_faceted_100_H` | 204 | {100} faceted G1 cut — more SiH2-like sites |
| `09_sphere_13A_H_fuse` | 112 | 13 Å sphere — **SiH2+SiH2 clash fusion** |

Surface group counts and defect tallies are in each `meta.json` and aggregated in `index.json`.

## Regenerate init bundles (01 + 03 only)

```bash
node tests/tSiNCs/generate_si_passivation_examples.mjs
```

## Run relax + vibration pipeline (all crystals)

```bash
cd /path/to/FireCore
export LSAN_OPTIONS=detect_leaks=0
ASAN="$(g++ -print-file-name=libasan.so)" && [ -f "$ASAN" ] && export LD_PRELOAD="$ASAN"

for d in tests/tSiNCs/fixtures/si_1nm_passivation/*/; do
  python3 -m pyBall.nanocrystal_pipeline relax --init-npz "$d/01_crystal.npz" \
    --out-npz "$d/02_relaxed.npz" --out-xyz "$d/relaxed.xyz" --allow-unconverged
  python3 -m pyBall.nanocrystal_pipeline hessian --relaxed-npz "$d/02_relaxed.npz" \
    --topology-npz "$d/03_topology.npz" --out-npz "$d/04_hessian.npz"
  python3 -m pyBall.nanocrystal_pipeline spectrum --hessian-npz "$d/04_hessian.npz" \
    --out-npz "$d/05_spectrum.npz" --out-plot "$d/spectrum.png"
done
```

**Do not** rebuild `03_topology.npz` after relax — topology params are fixed at export.

## Viewers

- C++: `./tests/tSiNCs/run_cpp_mol_browser.sh tests/tSiNCs/fixtures/si_1nm_passivation`
- Python: `./tests/tSiNCs/run_vispy_mol_browser.sh` (same path)

Load `01_crystal.npz` or `02_relaxed.npz` for geometry; companion `03_topology.npz` enables bond sticks and AABB overlay. Do not open `04_*` / `05_*` in the viewer grid.
