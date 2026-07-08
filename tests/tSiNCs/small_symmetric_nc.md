---
type: guide
title: Small symmetric nanocrystals
description: H-passivated Si and C diamond sphere batch under 100 atoms
tags: [nanocrystal, generation]
---

# Small symmetric nanocrystals (&lt;100 atoms)

Deterministic centred sphere cuts — no defects, `resolveClashes: 0`. Outputs land in **`tests/tSiNCs/OUT_small_nc/`** (relative to `TEST_DIR`, not repo root).

## Files

- `tests/tSiNCs/make_small_symmetric_nc.mjs` — generator
- `tests/tSiNCs/small_symmetric_nc.json` — material list and radii
- `tests/tSiNCs/OUT_small_nc/index.md` — atom-count table after run

## Run

```bash
node tests/tSiNCs/make_small_symmetric_nc.mjs
# from anywhere:
node /path/to/FireCore/tests/tSiNCs/make_small_symmetric_nc.mjs
```

Custom config: pass JSON path as first argument; `outDir` in JSON is relative to `tests/tSiNCs/` unless absolute.

## Sizes

Edit `materials[].cuts[].sphereR`. Atom counts jump in discrete lattice shells — check `index.md` after each run.

| Material | CIF | `applySymmetry` |
|----------|-----|-----------------|
| Si | `cpp/common_resources/crystals/Si-sym.cif` | 1 |
| C diamond | `cpp/common_resources/crystals/diamond_primitive.cif` | 0 |

No Python / MMFF — JS generation only (`Nanocrystals.js` + `cpp/common_resources`).
