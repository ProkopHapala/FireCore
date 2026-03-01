## Tutorial

### Why we are doing this (purpose & motivation)
We need training/diagnostic samples to **fit a force field against DFT references** for a molecule on a stepped NaCl surface. The workflow builds a **5D sampling grid** (u,v,y,tilt,phi) that covers plausible adsorption poses. The coarse+1D scan modes are **diagnostic**: they generate small, interpretable lines you can quickly visualize and compare to DFT to spot issues (geometry drift, symmetry bugs, wrong tilts, etc.).

### Physical & mathematical background (quick recap)
- **SDF isolines**: We compute a signed-distance field to the step and take isolines at prescribed `u` levels. Each isoline is parameterized by `v` between fixed x-boundaries.
- **Conformal streamline tracking**: `u` walks along streamlines; `v` samples along each streamline segment, keeping topology consistent between levels.
- **Two-corner lattice**: The step can expose **Na** or **Cl** at the corner. We build both UV grids (Na-corner, Cl-corner) and **blend x,z across y** (cosine by default) to represent one lattice constant along the step direction.
- **Local “up” from streamline**: Tilt is measured from the **local connector along u** (streamline tangent) at fixed v, not from the global z-axis. This yields a per-(u,v) local frame; phi then rotates around the tilted axis.
- **Triangle-fan tilts + symmetric phi**: Tilts sample a spherical cap; phi is reduced by symmetry (m-fold) and can be limited by `phi_max`.

### What this does
You want to sample adsorption **poses** of a molecule on a NaCl step edge in **5D**:
- **Position**:
  - `u` = distance shell index (SDF isoline)
  - `v` = index along an isoline segment between fixed x-boundaries
  - `y` = position along the step direction (one unit cell)
- **Orientation**:
  - `tilt` = direction of molecule axis on a spherical cap (from triangle-fan sampling)
  - `phi` = rotation around that axis (symmetry-reduced)

Key points:
- The 2D (x,z) grid is generated twice (Na-corner and Cl-corner) and then **blended across y**.
- “Tilt from z-axis” is **not used anymore** in pose generation:
  - tilt is measured from the **local streamline direction** (the connector along `u` at fixed `v`), so each (u,v) has its own “up”.

### Main script for combined sampling: [sample_driver.py](cci:7://file:///home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/sample_driver.py:0:0-0:0)
Run it from:
[doc/py/WrapedSamplinggGrid/](cci:9://file:///home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid:0:0-0:0)

#### Mode: counts (no files)
Prints grid sizes, coarse grid sizes, and estimated scan totals:
```bash
python sample_driver.py --mode counts
```

#### Mode: scan (1D scans + plots + many xyz files)
This generates:
- a PNG showing scan line overlays (`--png ...` or default)
- a PNG showing local up-vectors: `<out_prefix>_up_vectors.png`
- **separate XYZ for every 1D scan line**, i.e. one file per fixed combination of the other parameters

Example (scan all axes, using your coarse defaults):
```bash
python sample_driver.py --mode scan --scan_axes u,v,y,tilt,phi --png scan_lines.png
```

To limit file explosion, reduce coarse lists, e.g.:
```bash
python sample_driver.py --mode scan --scan_axes u \
  --coarse_v 0,8 --coarse_y 0,4 --coarse_tilt 0 --coarse_phi 0,4
```

#### Mode: random (random samples)
```bash
python sample_driver.py --mode random --Nrandom 500 --seed 1
```

### Angular-only visualization / orientation movie: [tilt_sampling_fan.py](cci:7://file:///home/prokophapala/git/FireCore/doc/py/WrapedSamplinggGrid/tilt_sampling_fan.py:0:0-0:0)
Shows the fan/disk/dome/phi sampling and can generate a pure orientation movie:
```bash
python tilt_sampling_fan.py --fans 4 --subdiv 2 --theta_max 45 --phi_samples 8 --phi_max 180
```

### Output naming conventions
- Scan-line XYZ files are named like:
  - `samples_u_line_v=0_y=4_tilt=0_phi=4.xyz`
Meaning: “scan along axis u, with fixed v=0, y=4, tilt=0, phi=4”.

### Controlling u-levels with `u_list`
- By default, `sample_driver.py` builds `u` levels from `u_min/u_max/u_levels` (linear/exponential schedules).
- To **pin exact isolines**, pass `--u_list` as a comma-separated list (Å). Example:
  ```bash
  python sample_driver.py --mode scan --u_list 0.1,0.5,1.0,2.0 --scan_axes u
  ```
- When `u_list` is provided, it overrides `u_levels`/`u_min`/`u_max` and is used for both Na and Cl grids before y-blending.

### Practical workflow for students
1) Start with:
```bash
python sample_driver.py --mode counts
```
2) Decide coarse indices (small numbers first).
3) Generate one axis scan (e.g. u) to validate:
```bash
python sample_driver.py --mode scan --scan_axes u --coarse_v 0 --coarse_y 0 --coarse_tilt 0 --coarse_phi 0
```
4) Inspect:
- `*_scan_lines.png`
- `*_up_vectors.png`
- a couple of XYZ files in a viewer
5) Scale up coarse lists and add other axes scans.

## Status
- **Task 1 implemented**: y blending between Na/Cl xz grids (cosine blend).
- **Task 2 implemented**: tilt measured from local streamline connector, with debug plot.
- **Tests run successfully**.

If you want the blend mode selectable from CLI, I’ll add `--y_blend {cos,linear}` next.