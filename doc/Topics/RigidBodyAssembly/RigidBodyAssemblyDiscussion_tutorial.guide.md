## 1) What the module does (scientific overview)
- We sample a large set of rigid-body poses of a molecule inside a periodic lattice, enforcing a 6-fold hexagonal “flower” symmetry. The evaluated lattice tiling is controlled by `--nPBC_test` (e.g., 1→3×3, 2→5×5).
- **Rotation Sampling**: Can be full 3D (Super-Fibonacci sequence), perfectly flat 2D in-plane rotations, or 2D in-plane rotations with small orthogonal tilts.
- **Translation Sampling**: In-plane shifts are tiled over the fractional coordinates of the unit cell vectors. We support a compact triangular shift region (`ca + cb <= 1`) to eliminate redundancy in highly symmetric cells.
- Before GPU evaluation, we compute each pose’s molecular thickness (z-span) on the host and discard any pose with z-span ≥ `--zspan_max`. This prunes the search while keeping dense sampling.
- GPU scoring (`evaluate_packing_3d`) returns **two metrics** per pose: clash sum (quadratic overlaps with early abort at `--penalty`) and **minimum inter-molecular atom distance**. Both are later filtered on the host.
- A composite score = clash + `--zpenalty` * z-span is used for ranking on host.
- Exports/plots are strictly filtered: `(zspan < --zspan_max)`, `(clash < --clash_max)`, and `(min_dist >= --dist_min)`. No silent fallbacks.

## 2) Key components

### Host harness: `AssemblyOCL`
- Loads and caches kernels `emit_configuration_xyz` and `evaluate_packing_3d`.
- `evaluate_packing`: flattens transforms, launches one workgroup per pose, returns **both clash sums and minimum distances**.
- `emit_configuration`: emits all transformed atom positions for a chosen pose (used to write XYZ frames).

### Transform generation utilities & Diagnostics
- `generate_rotations(mode, nrot, tilt_range, n_tilt)`: Generates rotation matrices based on mode (`full3d`, `inplane`, `tilt`). In `tilt` mode, rotations are ordered yaw-fast for smooth inspection.
- `custom_generate_transform_buffer(...)`: Builds rotations, in-plane shifts with optional compact masking, applies 6-fold symmetry, tiles lattice replicas according to `--nPBC_test`, flattens to `(n_confs, nmols, 4, 4)`. Export/plot replicas are independently controlled via `--nPBC_xyz` but cannot exceed what was evaluated.
- `plot_translations`, `plot_rotations`, `make_toy_molecule_movie`, `make_sym_toy_movie`: Diagnostic tools to visually verify that the translational and rotational sampling spaces are dense, physically correct, and properly symmetrized.

### CLI driver: `tests/tMMFF/test_assembly.py`
- Args: `--preset`, `--mol`, `--cell`, `--nrot`, `--rot_mode`, `--n_tilt`, `--tilt_range`, `--nshift`, `--shift_range`, `--shift_region`, `--shift_sum_max`, `--penalty`, `--zpenalty`, `--zspan_max`, `--clash_max`, `--dist_min`, `--export_max`, `--nPBC_test`, `--nPBC_xyz`, `--xyz_simple`, `--radius`, `--wg`, `--device`, `--dump`, `--plot_trans`, `--plot_rot`, `--plot_toy`, `--plot_toy_sym`, `--plot_best_k`, `--z_highlight`, `--top_k`, `--outdir`.
- Presets: small (HCCOCCH), helicene, triptycene.
- Workflow: build atoms → generate transforms → prefilter by z-span → GPU scoring (clash + min-dist) → host filtering → exports/plots.
- Outputs (when `--dump` or plotting flags are used):
  - **Double Pareto plot**: left = clash vs z-span; right = min-dist vs z-span. Dense black = all configs; red = exported; blue line = Pareto fronts; dashed guides = thresholds.
  - **XYZ movie** of all exported configs with metadata. `--nPBC_xyz` chooses a centered lattice block for export (must be ≤ `--nPBC_test`); capped at `--export_max`; element names simplified if `--xyz_simple`.
  - **Top-K Best-config 2D PNGs** using the same replica selection as the movie, with depth-fog alpha and red highlights for highest atoms; bonds and supercell bounds overlaid.
  - **Sampling Diagnostics**: translation plot, rotation plot, toy orientation movie, symmetrized toy movie.
- Default output dir: `tests/tMMFF/assembly_out_<preset>/`.

## 3) How to run (examples)

### Triptycene preset (GPU) with Tilt Mode and Diagnostics
Perfect for 2D molecular adlayers where molecules should mostly lie flat but need slight tilting to interlock bulky groups (like triptycene paddles). We use compact triangular translation limits (`--shift_region triangle`) to prevent testing redundant shifts, and cap the outputs to prevent massive dumps from dense grids.
```bash
python tests/tMMFF/test_assembly.py --preset triptycene --rot_mode tilt \
  --nrot 32 --n_tilt 7 --tilt_range 0.25 \
  --nshift 20 --shift_region triangle --shift_sum_max 1.0 \
  --zspan_max 40.0 --clash_max 2 --export_max 100 --nPBC_test 1 --nPBC_xyz 0 \
  --plot_rot --plot_trans --plot_best_k 3 --z_highlight 0.4 --dump
```
- Generates base orientations and compact shifts, visualizes the sampling density and symmetry, runs fast GPU evaluation, and reports the top 3 best packings as height-aware PNGs. The movie output only contains the central cell + 6-fold symmetry (`--nPBC_xyz 0`) and is capped at 100 frames.

### Helicine with tighter shifts and stricter clash (Full 3D)
```bash
python tests/tMMFF/test_assembly.py --preset helicene --dump \
  --rot_mode full3d --nrot 800 --nshift 12 --shift_range 0.4 \
  --zspan_max 10 --clash_max 5
```

### Custom molecule/lattice with In-Plane rotation
```bash
python tests/tMMFF/test_assembly.py \
  --mol path/to/molecule.xyz \
  --cell "lvs 32.7 0 0  16.35 28.319 0  0 0 40" \
  --rot_mode inplane --nrot 36 --nshift 10 --dump
```

## 4) Interpreting outputs
- **Diagnostic Plots**: 
  - `trans_sampling_*.png`: Verifies that translation shifts densely cover the chosen unit cell region (e.g. triangle vs square).
  - `rot_sampling_*.png`: Shows how the primary axes of the molecule explore the rotational sphere. For `tilt` mode, you will see a dense band around the equator for X/Y and a tight cluster at the poles for Z.
  - `toy_rotations_*.xyz`: Load in a viewer to watch a simple 4-atom triad sweep through the exact sampled orientations systematically (yaw sweeps first).
  - `toy_rotations_sym_*.xyz`: Verifies that the 6-fold flower symmetry (and optionally the 3x3 lattice tiling) is correctly applied to the generated orientations.
- **Double Pareto plot**: clash vs z-span (minimize both); min-dist vs z-span (maximize min-dist, minimize z-span). Dense black = all; red = exported; blue = Pareto front; dashed = thresholds.
- **XYZ movie**: Each frame has a comment line:
  `# idx=<orig_conf> rank=<order> clash=<c> mindist=<dmin> zspan=<z> transform: r00 r01 ... r22 tx ty tz`
  Replica selection uses the centered lattice block defined by `--nPBC_xyz` (must be ≤ `--nPBC_test`).
- **Best-config PNGs (`assembly_best_*_<rank>.png`)**: Overhead 2D view. Atoms have a depth-fog effect (lower atoms are more transparent/grey). Atoms near the absolute top (within `--z_highlight` limit) are colored solid red to simulate AFM highest-point visibility. Uses the same replica selection as the XYZ.

## 5) Scientific/technical notes
- **Penalty function**: Quadratic in overlap depth. Early exit at `--penalty`.
- **Sampling quality & Modes**: 
  - `full3d`: Near-uniform SO(3) coverage via Super-Fibonacci. Good for bulk packing.
  - `inplane`: Strict rotation around Z. Good for perfectly flat surfaces.
  - `tilt`: Combines Z-rotation with small X/Y tilts. **Highly recommended for surface science**, as realistic adlayers often have slight out-of-plane tilts to minimize steric hindrance.
- **Performance**: Z-span prefiltering is done on the host. GPU evaluation is extremely fast (sub-0.1 s for hundreds–thousands of configs). Adjust `--nrot/--nshift` upward freely; prune with `--zspan_max` and `--shift_region`.

## 6) Suggested experiments for students
1) Compare `--rot_mode inplane` vs `--rot_mode tilt`. How much does a slight tilt (e.g. `--tilt_range 0.1`) improve the best clash score?
2) Inspect `trans_sampling_*.png` and `rot_sampling_*.png` to gain an intuition for how `--nrot` and `--nshift` translate to physical coverage density.
3) Use `--plot_toy_sym` and `--toy_sym_replicas 54` to visually inspect how the dense periodic boundary conditions are applied.
4) Tighten `--clash_max` or raise `--dist_min` to see how the Pareto fronts change.
5) Modify radii via `--radius` to study sensitivity of packing to vdW assumptions.
