---
description: Silicon nanocrystal generator (scripts/gen_nanocrystals.mjs)
---

# Silicon Nanocrystal Generator

Reference for developers and students on how to generate silicon nanocrystals with configurable cutting planes, pruning, capping, and surface bridge operations (collapse/insert) using `scripts/gen_nanocrystals.mjs`.

## What it does

- Builds a replicated Si supercell from CIF (`cpp/common_resources/crystals/Si-sym.cif` by default).
- Applies plane cuts to shape the nanocrystal (e.g., {111} and {100}).
- Recalculates bonds using MM parameters (ElementTypes/AtomTypes/BondTypes/AngleTypes).
- Prunes undercoordinated Si (iteratively).
- Adds caps (H) to dangling bonds.
- Optionally performs surface bridge insertion (adds SiH2 bridges) and collapse (removes SiH2 bridges, bonds the neighbors).
- Outputs per-sample `.mol2` and optional stacked multi-frame `.xyz` for Jmol/visual inspection.

## CLI (key options)

```
node scripts/gen_nanocrystals.mjs [options]
```

- Geometry & planes:
  - `--cif <path>`: CIF file (default Si-sym.cif)
  - `--nx-range a,b` `--ny-range a,b` `--nz-range a,b`: replication ranges (single value allowed)
  - `--centered 0|1`: center replication
  - `--planeTemplates a111,a100,...`: plane families
  - `--planeSymC <float>`: symmetry scaling for templates
  - `--planeMode ang|frac`: plane definition mode
  - `--planeCScale <float>` `--planeCJitter <float>`: scale/jitter plane offsets

- Bonding & pruning:
  - `--bondFactor <float>` `--defaultRcut <float>`: bonding thresholds
  - `--minSiDegree <int>` `--pruneMaxIter <int>`: prune undercoordinated Si iteratively

- Caps:
  - `--caps H|0|none`: add caps (H) or disable

- Surface bridges:
  - `--collapseProb <0..1>`: probability to collapse a selected surface SiH2 bridge (remove bridge, bond neighbors)
  - `--insertProb <0..1>`: probability to insert a SiH2 bridge on each surface Si–Si bond candidate
  - `--collapseAll 0|1`: legacy carbon bridge collapse

- Energetics (simple counting model):
  - `--E_SiH2` `--E_SiH3` `--E_bare` `--E_bridge` `--muH`

- Output:
  - `--outDir <path>`: output directory
  - `--prefix <str>`: file prefix
  - `--samples <int>`: number of samples (or use `--maxFiles`)
  - `--statsCsv <path>`: append stats (counts, E)
  - `--stackedXyzOut <path>`: write all frames into one multi-frame XYZ
  - `--seed <int>`: deterministic RNG

## Current surface-bridge logic (Si)

- Surface filter: `nSi < 4` (targets surface atoms).
- Insertion (surface-only):
  - Enumerate surface Si–Si bonds; with probability `insertProb`, call `insertBridge` (adds SiH2) and promote the new atom to Si (Z=14).
- Collapse:
  - Select surface Si with 2 heavy neighbors and ≥2 H (SiH2-like) and collapse with probability `collapseProb` (remove bridge + H, bond neighbors).
- No bond rebuild after bridge ops (to preserve newly added R1–R2 bonds).

## Usage examples

Fixed shape, random passivation:
```
node scripts/gen_nanocrystals.mjs \
  --samples 10 \
  --nx-range 2,2 --ny-range 2,2 --nz-range 2,2 \
  --planeTemplates a111,a100 --planeCScale 0.45 --planeCJitter 0 \
  --caps H --minSiDegree 2 \
  --collapseProb 0.3 --insertProb 0.2 \
  --seed 2 \
  --outDir OUT_fixedshape \
  --prefix fixed_si \
  --statsCsv OUT_fixedshape/stats.csv \
  --stackedXyzOut OUT_fixedshape/stacked.xyz
```

Smaller crystals for inspection:
```
node scripts/gen_nanocrystals.mjs \
  --samples 10 \
  --nx-range 1,2 --ny-range 1,2 --nz-range 1,2 \
  --planeTemplates a111,a100 --planeCScale 0.45 --planeCJitter 0.20 \
  --caps H --minSiDegree 2 \
  --collapseProb 0.3 --insertProb 0.2 \
  --seed 2 \
  --outDir OUT_nanocrystals_small10 \
  --prefix small_si \
  --stackedXyzOut OUT_nanocrystals_small10/stacked.xyz
```

## Troubleshooting / gotchas

- **Do not recalc bonds after bridge ops**: recalculating removes the newly formed R1–R2 bond; we removed those recalc calls in the bridge section.
- If inserts stay zero: check the debug log `surface Si-Si candidate bonds: N`. If N>0 but inserts=0, raise `insertProb` or move insertion earlier (before heavy pruning/capping). If N=0, loosen the surface filter or pruning.
- The script promotes inserted atoms to Si (Z=14) after `insertBridge` (which defaults to carbon geometry but works for SiH2 here).

## Development notes / lessons learned

- Root cause of missing R1–R2 bonds was bond recalculation after collapse; fixed by removing post-op recalc.
- Insertion initially failed due to overly tight filters and try/catch masking; simplified to deterministic candidate loop with surface-only filter and direct `insertBridge` calls.
- Debug logging (candidate count) is essential to see if the surface actually exposes insertable Si–Si bonds.
- Keep randomness explicit (`--seed`) when comparing passivation variants on fixed shapes.
