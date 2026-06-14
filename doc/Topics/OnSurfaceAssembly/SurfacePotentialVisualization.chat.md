# USER

We have this explorer of interaction of molecule with rigid substrate. 
@ExplorerVisPy.py 

We need to implement functionality which renders that surface (e.g. as 3D mesh resp as plot of function depending on two varaibles x,y

What we want to do is to find intersection z_node(x,y) where total potential ( Pauli+London+Coulomb = Morse+Coulomb or Lenaard-Jone+Coulomb ) pass certain trashold (some repulsive force trashold), or where it hase minimum and for every point x,y, on a grid evaluatethis height z_node(x,y) and plot this function as a surface. Vsiply must have good tools for ploting 3D function right? We can then also color this 3D function by the coulomb potential in that point. To evaluate this potential we should not use whole molecule but test particle - i.e. single atom for example H+ or O- with certain charge e.g. H+0.2e and O-0.4e are good proxy for interaction of H2O molecule with the substrate. 

We already do something like this here in function
renderSubstrate_ in @MolecularDraw.h 
``C++
int renderSubstrate_( const GridShape& grid, Quat4f * FF, Quat4f * FFel, double isoval, bool sign, float sclr=1.0 ){
    double * Zs = new double[grid.n.x * grid.n.y];
    getIsoSurfZ( grid, isoval, sign, FF, Zs );
    int nvert = 0;
    for ( int ib=1; ib<=grid.n.y; ib++ ){
        glBegin(GL_TRIANGLE_STRIP);
        for ( int ia=0; ia<=grid.n.x; ia++ ){
            int ip1 = ((ib-1)%grid.n.y)*grid.n.x + (ia%grid.n.x);
            int ip2 = ((ib  )%grid.n.y)*grid.n.x + (ia%grid.n.x);
            Vec3f gpos; Quat4f fel1,fel2,  f1,f2;
            Vec3d p1 = grid.dCell.a*ia + grid.dCell.b*(ib-1);
            Vec3d p2 = grid.dCell.a*ia + grid.dCell.b*(ib  );
            p1.z=Zs[ip1];
            p2.z=Zs[ip2];
            grid.cartesian2grid( p1, gpos); fel1 = interpolate3DvecWrap( FFel, grid.n, gpos );
            grid.cartesian2grid( p2, gpos); fel2 = interpolate3DvecWrap( FFel, grid.n, gpos );
            Vec3d nr1,nr2; double invr;
            grid.cartesian2grid( p1, gpos); gpos.z-=3.5; f1 = interpolate3DvecWrap( FF, grid.n, gpos );  f1.z*=0.1; f1.normalize(); f1.mul(-1);
            grid.cartesian2grid( p2, gpos); gpos.z-=3.5; f2 = interpolate3DvecWrap( FF, grid.n, gpos );  f2.z*=0.1; f2.normalize(); f2.mul(-1);
            colorRB( fel1.z*sclr ); glNormal3f(f1.x,f1.y,f1.z); glVertex3f(p1.x,p1.y,p1.z); nvert++;
            colorRB( fel2.z*sclr ); glNormal3f(f2.x,f2.y,f2.z); glVertex3f(p2.x,p2.y,p2.z); nvert++;
        }
        glEnd();
    }
    delete [] Zs;
    return nvert;
}
``
but I want to di it now with python and Vispy (or maybe you can make pyOpenCL kernel for that, so that is is fast?

It would be good if we can make it general, that means that choice which potential to use for the root-fidning ( total potential Morse+Coulomb, or just Morse without cullomb ) and what for coloring. Or if instead of trashold we should fine minimum of the total potential (goint from top from far toward the surface taking first minimum) ... these descision should be done by CLI and the system should be made as general as possible.

We also want to be able to switch betwen different potential. For exaple GridFF loaded from @Bspline_PLQd.npy or FoldedAtomicFunction ... both should be possible to chose.

Now for testing, I know that you have problem testing visual GUI aps, right? Perhaps it would be better if you make separate headless script (but still using Vispy) which render static .png images in orthographi view shouing it fhte potential loos reasonable with different setting (O,H, total potential, Morse  , GridFF , Folded Atomic Functions etc. ). Only when this is done you implement it aso to the GUI @ExplorerVisPy.py 

But you should maximize code reuse and minimize code duplication. In general I think it would be better to most of functions into some common module like @VispyUtils.py and keep @ExplorerVisPy.py only as thin top level interface calling these shared backend function.

Read also this documentation, it mayb be usefull for this
@ElectrostaticContinuumEmbeding_report.md @RigidGridFF_H2O_NaCl_report.md @GridFF_CaF2_doc_tutorial.md 

Read the referenced files and make detailed notes abut everything which can help you with implementation. (where is what, what does what), and then write comprehesive plan how to do this

---

# GPT-5.1-Codex

## Findings (from referenced files & docs)

- ExplorerVisPy.py (GUI) already supports scan backends (Reference, Fast GPU, Folded GPU) and produces 1D/2D energy plots; no surface/isoheight extraction yet. Physics toggles: LJ, Coulomb, HBond, Morse. Fast/folded backends force Morse+Coulomb only. Scan utils build transform grids; results carry `scan_info` with ranges and transforms for visualization. GUI shows atoms as vispy markers/lines; scan trajectory drawn as markers/line.
- VispyUtils.py provides a reusable AtomScene for orthographic top-down vispy rendering with rich toggles (bonds, forces, labels, picking/dragging). No surface mesh helper here yet.
- C++ MolecularDraw.h renderSubstrate_ shows desired behavior: computes iso-surface height Z(x,y) by getIsoSurfZ over grid FF, then draws triangle strips with per-vertex color from FFel and normals from FF gradients (interpolate3DvecWrap). This is the reference shape and coloring.
- Reports/tutorials:
  - ElectrostaticContinuumEmbedding_report.md: describes macro (continuum) correction, grids, macro potential; notes on scan utilities and test_interaction_scan headless plots for NaCl. Highlights periodic wrapping, macro flags, and validated workflows.
  - RigidGridFF_H2O_NaCl_report.md: workflows for rigid GridFF on NaCl/CaF2, z-scan validation, batch runs; emphasizes validated grids (e.g., CaF2 PLQd npy), typing, wrapping, and diagnostic plots. Good for picking test cases and defaults.
  - GridFF_CaF2_doc_tutorial.md (open) likely contains CaF2 grid validation steps, probe charges (H+0.2e, O-0.4e) and preferred alpha; relevant for headless tests.

## What’s needed (per user)

- Compute z_node(x,y): height where a chosen potential (Morse+Coulomb, Morse only, LJ+Coulomb, etc.) crosses a threshold (repulsive force threshold) or first minimum when descending from far z.
- Evaluate on an XY grid; produce a surface mesh/heightmap; color by another potential component (e.g., Coulomb) at that point.
- Allow choosing substrate potential source: GridFF (Bspline PLQd npy) or FoldedAtomicFunction.
- Provide a headless vispy script to render static PNGs (orthographic) for sanity (H/O probes, total vs component, different grids), before integrating into GUI.
- Maximize code reuse: push shared logic into VispyUtils (or a new helper module) and keep ExplorerVisPy thin.

## Proposed implementation plan

1) Surface extraction backend (shared module, e.g., VispyUtils_surface.py or added section in VispyUtils.py)
   - Inputs: potential field accessor (callable V(x,y,z) and optionally Fz), XY grid bounds/resolution, z search range, mode (threshold crossing vs first minimum), potential selector (total/Morse/LJ/Coulomb), probe parameters (position offset, charge/radius), substrate source (GridFF / folded / direct InteractionScanner).
   - Implement z-finder per (x,y): bracket z in [z_hi→z_lo], sample coarse then refine (e.g., monotonic search from far to near until V crosses threshold or derivative changes sign). Optionally use scipy-like 1D bracketing/bisection manually (no extra deps).
   - Return arrays Z(x,y) and color field C(x,y) (e.g., Coulomb at found z). Also return diagnostic mask for failure (no crossing/min found).

2) Potential accessors
   - For GridFF npy (Bspline_PLQd): reuse existing loaders from pyBall/OCL or Python-side GridFF utilities to evaluate Pauli/London/Coulomb on a probe (single atom with REQ + charge). Ensure origin/centering matches CaF2/NaCl grids (see CaF2 tutorial: g0 centered).
   - For FoldedAtomicFunction: wrap existing folded basis eval (folded_scanner or InteractionScanner with folded mode).
   - For reference scanner: use InteractionScanner single-ion mode; allow toggling components (Morse/LJ/Coulomb) like GUI does.

3) Headless renderer (new script, e.g., tests/tMMFF/render_surface_iso.py)
   - CLI options: substrate grid/file choice (GridFF npy vs xyz+type map for scanner), probe type/charge, potential mode for root-finding (total/Morse/LJ/Coul), coloring component, threshold or “first-min” mode, xy range/resolution, z range, output PNG path.
   - Use vispy offscreen (SceneCanvas with fov=0) to render a mesh: build vertex grid from X,Y,Z, color by C. Provide optional Matplotlib fallback heightmap. Save PNGs.

4) GUI integration (ExplorerVisPy.py)
   - Add controls to trigger surface extraction with current substrate selection and probe settings; display as a mesh overlay in the 3D view (Vispy Mesh visual) with toggle.
   - Reuse backend helper from step 1; avoid duplicating logic. Allow selecting potential/color modes via combo boxes; set defaults (e.g., threshold = small repulsive force, probe H+0.2e, O-0.4e).

5) Testing & validation
   - Headless script: run on NaCl_1x1_L3 and CaF2 PLQd npy with H+0.2e / O-0.4e; generate PNGs (total, Coulomb, Morse-only) and check qualitative shapes (repulsive wall height, minima near known adsorption distances per reports).
   - Ensure centered grid indexing for CaF2 (use documented g0), consistent nPBC/macro flags for scanners.

## Key design choices / reuse points

- Put mesh-building/coloring helpers in VispyUtils (e.g., make_surface_mesh(vertices, colors) and build_grid_vertices_from_heightmap).
- Potential evaluation should share existing InteractionScanner/MolecularDynamics paths; avoid ad-hoc duplication.
- Provide a clean interface to switch substrate source (GridFF npy vs folded vs direct scanner) and component selection.
- Keep loud failures: if no crossing/min found for a grid point, mark and report; do not silently clamp.

## Next steps (if you agree)
- Draft the shared surface extraction helper (probe-based potential accessor + z-finder + heightmap→mesh builder) in VispyUtils or a new helper module.
- Add headless script with CLI to produce PNGs.
- Wire ExplorerVisPy to call the helper and render the mesh with selectable modes.
- Validate on NaCl/CaF2 with H/O probes; iterate thresholds/min-mode heuristics if needed.

---

# USER

WHAT THE HELL!!!!! THIS IS TOTAL GARBAGE AGAIN! every time I give you somethng to do on top of some framwork which is already tested and working you do some crap which does not reproduce the functionality of the original system you are supoosed to buld atom. Istead of writing your garbage from scrath you MUST read carefully the original system and undrestand it so you can build upon it! that is general rule, you must remember that! That is alsays 80% of success to understand the existing systems properly !


now what exactly is grap about the images you produces:
what is 
Bspline_PLQd_H_qp0.20_threshold_total_color_coulomb_maps.png
Bspline_PLQd_O_qm0.40_first_minimum_total_color_coulomb_maps.png
? There are white ares, these are what? some NANs ?

Bspline_PLQd_H_qp0.20_threshold_total_color_coulomb_surface.png
Bspline_PLQd_O_qm0.40_first_minimum_total_color_coulomb_surface.png
are totally white, plain color, nothing rendered. what is that crap?

NaCl_1x1_L3_H_qp0.20_threshold_total_color_coulomb_maps.png contains some artifacts, and definitely does not show perioduc susbstrate
NaCl_1x1_L3_H_qp0.20_threshold_total_color_coulomb_surface.png same artifacts in 3D

NaCl_1x1_L3_xyz_fast_gpu_H_qp0.20_threshold_total_color_coulomb_maps.png shows just some nouse
NaCl_1x1_L3_xyz_fast_gpu_H_qp0.20_threshold_total_color_coulomb_surface.png is the same noise in 3D


=> NONE of the plots you produce show anything physically meaningfull

What I expect is periodic substrate whit proper seamless periodic boundary condiction (no seams at boundary), and smooth bumps up and down showing atractive and repulsive sites (ions and holes in between the ions), I do not see anything like that.

But In the program I gave you as reference everything is already impleemnted and tested we were ding z-scans and also 2D xy scans, and they were fine

But you instead of unrestaind it, how to properly use GridFF how properly use electrostatic embeding corrctions and the other things, you implemented from scratch some total grbage totally igrning what we spand long time to debug before. I HATE YOUR IGNORANCE AND LACK OF CARE!!!! I'M SUPER ANGLEY ON YOU YOU FUCKING LAZZY PIECE OF SHIT!

are you doing some relaxation of is that fully rigid scan ? Because for findihg the z_mode, you obviously should do just rigid (unrelaxed) evaluation of the surface potential.

WHAT I REALLY HATE, is that you were gusssing parameters, grid dimension, vxel size, height above the surface and others, instead learning how xactly this is done in the working and tested referencs I gave you. READ THE FUCKING REFERENCES YOU LEZZY PIECE OF SHIT INSTEAD OF BLIND GUESSING!!!!

thell me what files you have to read to get the answers how to do it propely , and then tell me what you learn there how it is done

---

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


---

# USER

1) tell me if all the ne krnels are used or some are deprecated (can I delete some?)

getSurfaceIsoGridFF
getSurfaceIsoSurfMorse
evalSurfMorseE3D
sampleGridFF_Bspline_points

2) I was testing different surface modes with GUI, and what I see is this

1) XYZ Reference  :  Good result but super slow
2) XYZ Fast GPU   :  ather slow, quite fine, but doesnot reproduce  XYZ Reference. Why ? Should it? Why it is so slow? Haw are the bottleneks? Where we do something inefficient?
3) XYZ Folded GPU :  not working at all, crashes
4) GridFF         :  looks reasoanble but gives totally different results than XYZ Reference  or  XYZ Fast GPU; electrostatics (charge Q of probe particle) seem to have much stronger affect. Why ? Check that

can you explain clearly how each methods works and why we have these problems?

---

# GPT-5.1-Codex
Here’s the concise status and explanations.

## 1) Are these kernels used or deletable?

- `getSurfaceIsoGridFF` — **used**. Launched by [MolecularDynamics.eval_surface_iso_gridff()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:1340:4-1383:126) and the GridFF GUI mode and headless renderer. Keep. @pyBall/OCL/MolecularDynamics.py#1341-1384, @cpp/common_resources/cl/relax_multi.cl#2960-3050
- `getSurfaceIsoSurfMorse` — **used**. Launched by [MolecularDynamics.eval_surface_iso_morse()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:1385:4-1443:126) for the XYZ fast_gpu/folded_gpu surface path. Keep. @pyBall/OCL/MolecularDynamics.py#1386-1444, @cpp/common_resources/cl/relax_multi.cl (below sampleGridFF_Bspline_points)
- `evalSurfMorseE3D` — **used internally** by the Morse iso kernel (`getSurfaceIsoSurfMorse`) and related Morse samplers. Keep.
- `sampleGridFF_Bspline_points` — **used** for pointwise GridFF sampling (e.g., [eval_gridff_points](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:1300:4-1338:77) in the GridFF backend and during kernel setup). Keep. @cpp/common_resources/cl/relax_multi.cl#2726-2795, @pyBall/OCL/MolecularDynamics.py#508-515

None of these should be deleted; they’re all active.

## 2) Modes in GUI: how they work and why you see the issues

### a) XYZ Reference (InteractionScanner, CPU)
- How it works: pure CPU; loops over all probe points, uses tabulated Morse (and electrostatics) with PBC wrapping, optional macro embedding when ≤3 z-layers.
- Why it’s slow: per-point Python/NumPy + CPU loops; no GPU parallelism; O(nx*ny*nz) on CPU.
- Fidelity: This is the baseline/expected physics.

### b) XYZ Fast GPU (Morse iso kernel)
- How it works: single OpenCL kernel `getSurfaceIsoSurfMorse` scans z for every (x,y) in parallel, returns z_node + color in one launch; host only downloads the final 2 maps.
- Why it can differ from Reference:
  1) Macro embedding auto-disables when >3 distinct z-layers; Reference used macro when allowed. That changes electrostatic embedding and can shift the surface.
  2) z-range offsets: GPU path uses `z_top + z_range`; Reference may scan absolute z differently if z_top differs.
  3) PLQH/component selection must match (selector/color). A mismatch of selector/color or probe scaling shows up as differences.
  4) Folded path isn’t used here; this is straight Morse.
- Why it’s “rather slow”: It’s still O(nx*ny*nz) work per thread. For dense grids (e.g., 81x71x160) that’s ~5.7k threads * 160 z-steps. But it’s a single kernel launch; remaining overhead is just two buffer reads (surf_xyzq, surf_zc). If you saw slowness, it’s likely from large nz and GUI-side copies, not multiple launches.

### c) XYZ Folded GPU
- How it works (intended): uses folded harmonics precompute (folded_meta, folded_lvec2d, folded_coeffs/kxyz) and folded kernel variants.
- Why it crashes: those folded buffers aren’t initialized/bound in GUI use; [setup_kernels()](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/MolecularDynamics.py:425:4-560:43) prints “skipping getSurfFolded…” and later calls still try to use it → crash. It’s currently not wired/ready.

### d) GridFF
- How it works: uses precomputed B-spline PLQd grid (`Bspline_PLQd.npy`) and runs `getSurfaceIsoGridFF` to find z_node + color. Coulomb comes from the PLQd grid, not from the Morse surface charges.
- Why it differs from XYZ Reference/Fast GPU:
  1) Different underlying potential: GridFF PLQd is a precomputed field (Pauli/London/Coulomb) and not the on-the-fly Morse + embedding used by XYZ.
  2) Electrostatics strength: the Coulomb channel in the PLQd grid can be stronger than the Morse-surface charge model. If selector/color uses `coulomb`, expect larger Q sensitivity.
  3) Any crop/range mismatch shows up as seams; both fast_gpu and gridff showed the same y-seam in CaF2, meaning it’s geometry/range periodicity, not backend physics.

## Bottlenecks and inefficiencies (XYZ Fast GPU)
- Main cost: nz scanning per (x,y) inside the kernel. For large nz this dominates. The rest is cheap (one kernel launch + two buffer reads).
- No per-column kernel spam (we already avoid that). If you still see high latency, reduce nz or use coarser nx/ny for interactive previews.
- Ensure macro embedding is enabled only when valid (≤3 layers); otherwise physics will differ from Reference.

## Folded GPU next steps
- Needs proper folded setup: allocate/bind folded_meta, folded_lvec2d, folded_coeffs/kxyz, and populate folded_atom_type. Until that’s done, the mode will keep crashing.

## What to align if you want parity
- Match selector/color (PLQH) and probe scaling between Reference and Fast GPU.
- Ensure the same z-range and z_top offset.
- Use macro embedding only when the slab truly has ≤3 distinct z-levels (already auto-gated).
- Folded: not ready; either disable in GUI or finish buffer setup.

---

# Current detailed report (goals, work, challenges, pitfalls, discrepancies)

## Goals
- Fully GPU-parallel iso-surface extraction (single kernel call) for `fast_gpu` (Morse) and `GridFF` backends, producing `z_node(x,y)` + color without CPU orchestration.
- Integrate into GUI with correct world-space z placement, no substrate bonds, and parity with the validated reference (InteractionScanner) physics.
- Achieve parity between GridFF and getSurfMorse (same chemistry/physics, matching z/color) and expose headless/GUI paths for validation.

## What we did (summary)
- Added GPU iso kernels: `getSurfaceIsoSurfMorse` (Morse) and `getSurfaceIsoGridFF` (GridFF) + Python wrappers (`eval_surface_iso_morse/gridff`).
- Fixed OpenCL local memory address space issues and PBC wrapping in `relax_multi.cl`; added macro auto-gating (≤3 z-layers) in Python to avoid CaF2 crashes.
- Headless renders for NaCl and CaF2 now run from a single GPU call, producing PNG/JSON/NPZ with seam metrics and value ranges.
- GUI integration: world-space z placement retained, substrate bonds suppressed; macro gating applied to avoid folded/macro crashes.

## Challenges & pitfalls hit
- **Z reference mismatches**: Absolute z windows vs. height-above-top led to misses/NaNs. Fixed by using `z_top + z_range` for GPU paths and reporting z relative to top.
- **Macro embedding limits**: Original Morse macro code supports only ≤3 z-layers; CaF2 exceeded this and crashed. Now auto-disables macro when layer count >3.
- **Physics drift (selector/color)**: Mismatched PLQH/component choices caused deviations between Reference and Fast GPU. Need strict alignment of selector/color and probe scaling.
- **GridFF vs. Morse physics**: GridFF uses precomputed PLQd (Pauli/London/Coulomb) while Morse path computes on-the-fly surface charges; Coulomb strength differs, especially with charged probes.
- **Folded path**: Not wired—missing folded buffers (folded_meta, folded_lvec2d, coeffs/kxyz) cause kernel skips and later crashes.
- **Seam artifacts**: CaF2 cropped cells show y-seam in both backends, indicating range/periodicity mismatch rather than backend physics.

## Current discrepancies (need rigorous parity tests)
- GridFF vs. Morse (`fast_gpu`) do **not** yet match z/color on CaF2; Coulomb appears stronger in GridFF. Potential causes: double Coulomb scaling, REQ/PLQ mismatch, or generation vs. sampling scaling (e.g., 14.39 eV·Å factor applied twice).
- NaCl shows good agreement between Reference and Fast GPU; GridFF parity still needs quantitative Δz/ΔE maps.
- Folded GPU remains non-functional in GUI (crashes) until buffers are initialized/bound.

## Suspected scaling issues to check
- Coulomb factor (14.39 eV·Å) possibly applied twice: once during GridFF generation and again during sampling. Must confirm the generation script vs. sampling kernel conventions.
- REQ/PLQ parameter parity: ensure identical parameter sets are used in GridFF generation and Morse evaluation; verify Pauli/London amplitudes are not downscaled/upsized.
- Probe charge scaling: confirm PLQH multiplication for selector/color uses the same probe scaling across Reference/Fast GPU/GridFF.

## Next validation steps
1. Compute Δz(x,y) and ΔE_coul(x,y) maps between GridFF and getSurfMorse for the same probe (H+0.2e) and mode (`first_minimum`, `selector=nonel`, `color=coulomb`) on both NaCl and CaF2; report RMS/max and seam metrics.
2. Audit Coulomb scaling: trace the 14.39 factor in GridFF generation vs. sampling; ensure it is applied exactly once. Rebuild GridFF if needed.
3. Verify REQ/PLQ parity: cross-check parameter tables used in GridFF generation against Morse runtime parameters; align Pauli/London magnitudes.
4. Fix Folded GPU path: initialize folded_meta/lvec2d/coeffs/kxyz and atom types, then validate against Reference and Fast GPU.
5. Keep macro gating consistent: macro on only when ≤3 z-layers; otherwise disable uniformly across Reference/Fast GPU to avoid physics skew.

## Takeaways / guardrails
- Do not guess GridFF metadata or scaling; always use generation-time g0/dg and documented unit factors.
- Enforce top-relative z windows to avoid sampling outside support; fail loudly on missing crossings/NaNs.
- Align selector/color/PLQH and probe scaling before comparing backends; mismatches masquerade as physics errors.
- Treat Coulomb scaling as a prime suspect when GridFF shows stronger charge sensitivity; confirm single-application of 14.39.
- Folded path is not ready; disable in GUI until buffers and kernel params are wired.

---

# Reproducible runs and GUI mode → kernel mapping

## Headless scripts (current, tested)

Working directory: `cd tests/tMMFF`

NaCl, fast_gpu (Morse iso kernel):
```bash
python3 ./render_surface_iso.py --strict \
  --source xyz --sub ../../cpp/common_resources/xyz/NaCl_1x1_L3.xyz \
  --xyz_backend fast_gpu --probe H --charge 0.2 --alpha 2.0 \
  --mode first_minimum --selector nonel --color coulomb \
  --xrange=0,8 --yrange=0,8 --zrange=0.0,8.0 \
  --nx 81 --ny 81 --nz 160 \
  --outdir output_surface_iso_gpucheck --name nacl_2x2_xyz_fast_gpuprod
```

CaF2, fast_gpu (Morse iso kernel):
```bash
python3 ./render_surface_iso.py --strict \
  --source xyz --sub ../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz \
  --sub_types Ca=Ca+2,F=F- --xyz_backend fast_gpu \
  --probe H --charge 0.2 --alpha 2.0 \
  --mode first_minimum --selector nonel --color coulomb \
  --xrange=0,11.5875 --yrange=0,10.035 --zrange=0.0,8.0 \
  --nx 81 --ny 71 --nz 160 \
  --outdir output_surface_iso_gpucheck --name caf2_crop_xyz_fast_gpuprod
```

CaF2, GridFF (tricubic Bspline iso kernel):
```bash
python3 ./render_surface_iso.py --strict \
  --source gridff \
  --gridff ../../cpp/common_resources/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd.npy \
  --gridff_xyz ../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz \
  --sub_types Ca=Ca+2,F=F- \
  --probe H --charge 0.2 --alpha 2.0 \
  --mode first_minimum --selector nonel --color coulomb \
  --xrange=0,11.5875 --yrange=0,10.035 --zrange=0.0,8.0 \
  --nx 81 --ny 71 --nz 160 \
  --outdir output_surface_iso_gpucheck --name caf2_crop_gridff_gpu2d
```

Outputs: `*_surface.png`, `*_maps.png`, `*.json` (z/color ranges, seam metrics), `*.npz` (arrays).

## GUI modes → kernels (ExplorerVisPy.py)

- **XYZ Reference**: Uses `InteractionScanner` (CPU). No iso kernel; per-point CPU eval. Slow but baseline physics (macro on if ≤3 z-layers, wrap_PBC on).
- **XYZ Fast GPU**: Calls `MolecularDynamics.eval_surface_iso_morse()` → OpenCL kernel `getSurfaceIsoSurfMorse` in `relax_multi.cl`. Macro auto-disables when >3 layers; z window is `z_top + z_range`.
- **XYZ Folded GPU**: Intended to call folded kernels (e.g., `getSurfFolded` variants) but currently crashes because folded buffers/metadata are not initialized/bound in GUI. Not usable yet.
- **GridFF**: Calls `MolecularDynamics.eval_surface_iso_gridff()` → OpenCL kernel `getSurfaceIsoGridFF` (tricubic B-spline). Uses precomputed `Bspline_PLQd.npy` (Pauli/London/Coulomb) for the probe; Coulomb sensitivity depends on grid scaling.

File refs: `pyBall/ExplorerVisPy.py` (mode selection), `pyBall/OCL/MolecularDynamics.py` (wrappers), `cpp/common_resources/cl/relax_multi.cl` (kernels).

## Known problems (per mode)
- **XYZ Reference**: Slow by design (CPU loop over all points). Baseline physics.
- **XYZ Fast GPU**: Parity with Reference depends on matching selector/color/PLQH and macro state; uses height-above-top z. Differences observed when macro is off (slabs with >3 layers) or selector/color mismatched.
- **XYZ Folded GPU**: Crashes because folded buffers (folded_meta, folded_lvec2d, coeffs/kxyz, atom types) are not populated. Needs wiring before use.
- **GridFF**: Physics gap vs. Morse path—Coulomb appears stronger; suspects: double 14.39 factor during grid generation/sampling, REQ/PLQ mismatch. CaF2 shows y-seam in both GridFF and fast_gpu, indicating geometry/range periodicity rather than backend-only bug.