# NPZ viewer fixtures

Small committed NPZ bundles for **C++** and **Python Vispy** molecular browsers (fast smoke tests).

**Schema:** [`doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md`](../../../../doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md)  
**Full pipeline gallery:** [`../si_1nm_passivation/README.md`](../si_1nm_passivation/README.md) (01→05 per crystal)

## Files

| File | Atoms | Description |
|------|------:|-------------|
| `diamond_sphere_R6_caps.npz` | 270 | C diamond sphere R=6 Å, H caps only (crystal geometry) |
| `diamond_sphere_R6_defects_crystal.npz` | 288 | C sphere R=6 with insert+collapse defects |
| `diamond_sphere_R6_defects_topology.npz` | 288 | Full MMFFL topology + surface group AABBs |
| `si_planes_G2_defects_topology.npz` | 289 | Si {111} plane cut G2 (nx=2), defects, topology only |
| `*_meta.json` | — | Human-readable defect counts and gen_params |

## View

```bash
./tests/tSiNCs/run_cpp_mol_browser.sh tests/tSiNCs/fixtures/npz_viewer
./tests/tSiNCs/run_vispy_mol_browser.sh tests/tSiNCs/fixtures/npz_viewer
```

## Regenerate

```bash
# C caps (crystal only)
node scripts/export_nanocrystal_bundle.mjs \
  --cif cpp/common_resources/crystals/diamond_primitive.cif \
  --sphere 6 --insertProb 0 --collapseProb 0 --seed 1001 \
  --out /tmp/nc_caps --id diamond_sphere_R6_caps --crystalOnly 1
cp /tmp/nc_caps/01_crystal.npz tests/tSiNCs/fixtures/npz_viewer/diamond_sphere_R6_caps.npz

# C defects (crystal + topology)
node scripts/export_nanocrystal_bundle.mjs \
  --cif cpp/common_resources/crystals/diamond_primitive.cif \
  --sphere 6 --insertProb 0.15 --collapseProb 0.10 --seed 1 \
  --out /tmp/nc_defects --id diamond_sphere_R6_defects
cp /tmp/nc_defects/01_crystal.npz tests/tSiNCs/fixtures/npz_viewer/diamond_sphere_R6_defects_crystal.npz
cp /tmp/nc_defects/03_topology.npz tests/tSiNCs/fixtures/npz_viewer/diamond_sphere_R6_defects_topology.npz
cp /tmp/nc_defects/meta.json tests/tSiNCs/fixtures/npz_viewer/diamond_sphere_R6_defects_meta.json

# Si planes G2 defects (topology)
node scripts/export_nanocrystal_bundle.mjs \
  --cif cpp/common_resources/crystals/Si-sym.cif \
  --planes a111 --nx 2 --planeCScale 0.34 --planeCJitter 0 \
  --insertProb 0.15 --collapseProb 0.10 --seed 3 \
  --out /tmp/nc_si --id si_planes_G2_defects
cp /tmp/nc_si/03_topology.npz tests/tSiNCs/fixtures/npz_viewer/si_planes_G2_defects_topology.npz
cp /tmp/nc_si/meta.json tests/tSiNCs/fixtures/npz_viewer/si_planes_G2_defects_meta.json
```

## Validation

```bash
node tests/tSiNCs/test_nanocrystal_defect_export.mjs
./tests/tSiNCs/test_cpp_npz_load.sh
PYTHONPATH=. pytest tests/tSiNCs/test_crystal_npz_load.py -q
```
