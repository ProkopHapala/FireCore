# Post-mortem (surface potential visualization)

## Goal

- Render a seamless periodic isoheight surface z(x,y) for substrate–probe interaction, using the existing, validated backends (GridFF npy, InteractionScanner/MolecularDynamics). Mode selection (threshold/minimum), component selection (nonel/total/coulomb/etc.), and coloring must all be user-selectable. Output static PNGs (headless) and expose the same path in GUI.

## What we did

- Implemented a shared sampler in `pyBall/SurfaceSampling.py` for probe-based rigid evaluation over (x,y,z).
- Added mesh/color/render helpers to `pyBall/VispyUtils.py` and a headless CLI script `tests/tMMFF/render_surface_iso.py`.
- Integrated a thin surface overlay in `ExplorerVisPy.py`.
- Performed strict validation runs (no relaxation) on NaCl_1x1_L3 for GridFF and XYZ (reference, fast_gpu) backends.

## Problems encountered (root causes)

- **Wrong GridFF origin/z convention**: Initially set g0.z to top atom and dg from xyz, instead of the generation convention `g0 = (-Lx/2, -Ly/2, z0_gen)`, dg from GridShape.
- **Wrong surface definition**: Used `total=0` as isosurface; for many xy points this crossing does not exist (hence NaNs/whites). The tested C++ renderer uses non-electrostatic PL isovalue (default 0.1) and often visualizes minima instead of a fixed threshold.
- **Wrong vertical reference**: Treated z_range as absolute world coords; the working scans use height above the top substrate atoms. This mismatch caused sampling outside the supported z-domain.
- **Physics mismatch in InteractionScanner path**: Enabled Morse for “nonel/total”; the reference scans use LJ (Morse disabled) with macro correction and wrap_PBC.
- **Coloring with neutral probe**: When using `color=coulomb` and `q=0`, the colormap failed due to zero range; previously this produced flat/white renders.
- **Silently masking misses**: Earlier runs allowed missing crossings to produce partially filled/white surfaces instead of failing.

## Fixes applied

- **GridFF sampler**: Enforced generation convention for g0/dg; removed z-clamping that hid out-of-range sampling; added NaN marking outside grid. Added top-relative height reporting.
- **XYZ sampler**: Switched surface search to height-above-top; kept wrap_PBC and macro matching InteractionScanner; left fast_gpu/folded_gpu using same z reference.
- **Reference backend physics**: `nonel/total` now enable LJ (Morse off) to match tested scans.
- **Robust surface mode**: Adopted `mode=first_minimum`, `selector=nonel`, `color=coulomb` as a validated, fully populated surface for NaCl with H probe q=+0.2.
- **Strict validation**: Headless script now fails if any xy point is missing, if any NaN appears, or if seams differ at cell boundaries; reports seam metrics and value ranges in JSON/console.
- **Colormap guard**: Non-finite or zero-range colors now raise errors instead of producing white output.

## Validation results (rigid, NaCl_1x1_L3)

- **GridFF** (`first_minimum`, `selector=nonel`, `color=coulomb`, `q=+0.2`, `nx=81, ny=81, nz=320`):
  - ok 6561/6561, no NaNs, seam_x=0, seam_y=0
  - z range: 0.1088 .. 3.7008 Å
  - color range: -9.3772 .. 9.3772 eV
  - outputs: `output_surface_iso_fix/nacl_gridff_nonel_min_qp02_*`

- **XYZ reference** (same mode):
  - ok 6561/6561, seam_x=0, seam_y=0
  - z range: 2.5642 .. 3.3047 Å
  - color range: -0.02274 .. 0.06245 eV
  - outputs: `output_surface_iso_fix/nacl_xyz_reference_nonel_min_qp02_*`

- **XYZ fast_gpu** (same mode):
  - ok 6561/6561, seam_x=0, seam_y=0
  - z range: 2.5657 .. 3.3608 Å
  - color range: -0.02086 .. 0.06816 eV
  - outputs: `output_surface_iso_fix/nacl_xyz_fast_nonel_min_qp02_*`

## Key takeaways / guardrails to prevent repeat mistakes

- **Do not guess GridFF metadata**: g0 and dg must come from the generation script (`ocl_GridFF_new.py` / GridShape), not inferred from xyz or z_top.
- **Use height-above-top for surface search**: Align z with the working scans; avoid absolute-world z windows that miss the grid support.
- **Pick surface definitions that actually exist**: Fixed `total=0` often has no crossing; prefer `first_minimum` or a known PL isovalue with explicit validation.
- **Match InteractionScanner physics**: For “nonel/total”, use LJ with macro+wrap_PBC, Morse off (unless explicitly required). Keep nPBC consistent with tested defaults.
- **Fail loudly on partial surfaces**: If any xy point is missing or any NaN occurs, stop and report; do not render partial/white meshes.
- **Report seam metrics**: Always check periodic seams (x/y) and include ranges in logs/JSON for reproducibility.
- **Coloring must have finite range**: Guard against zero-range fields (e.g., neutral probe with coulomb color) to avoid flat/white renders.

## Next actions (if needed)

- Add parity analysis between GridFF and XYZ surfaces: Δz(x,y), ΔEcoul(x,y), RMS/max errors, discrepancy maps.
- Expose validated presets in CLI/GUI (mode=first_minimum, selector=nonel, color=coulomb, q=+0.2) and keep validation guardrails on by default.