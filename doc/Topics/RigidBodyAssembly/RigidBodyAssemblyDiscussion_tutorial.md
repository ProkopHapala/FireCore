## 1) What the module does (scientific overview)
- We sample a large set of rigid-body poses of a molecule inside a periodic lattice, enforcing a 6-fold hexagonal “flower” symmetry and a 3×3 lattice tiling (54 replicas total per pose). Rotations are drawn from a low-discrepancy Super-Fibonacci sequence; in-plane shifts are tiled over the cell vectors. @pyBall/OCL/Assembly.py#89-151
- Before GPU evaluation, we compute each pose’s molecular thickness (z-span) on the host and discard any pose with z-span ≥ `--zspan_max` (default 10 Å). This prunes the search while keeping dense sampling. @tests/tMMFF/test_assembly.py#152-182
- Remaining poses are scored on GPU by the OpenCL kernel (`evaluate_packing_3d`): it sums quadratic overlaps `(r_sum - dist)^2` for atom pairs whose distance is below the sum of their radii, with early abort if the running score exceeds `--penalty`. (Kernel defined in [cl/Assembly.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/cl/Assembly.cl:0:0-0:0); host launch shown at @pyBall/OCL/Assembly.py#22-53)
- A composite score = clash + `--zpenalty` * z-span is used for ranking on host. @tests/tMMFF/test_assembly.py#182-207
- We then export and plot only “reasonable” poses: `zspan < --zspan_max` and `clash < --clash_max` (default 10). @tests/tMMFF/test_assembly.py#194-234

## 2) Key components

### Host harness: [AssemblyOCL](cci:2://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:5:0-76:29)
- Loads and caches kernels `emit_configuration_xyz` and `evaluate_packing_3d`. @pyBall/OCL/Assembly.py#6-77
- [evaluate_packing](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:21:4-52:27): flattens transforms, launches one workgroup per pose, returns clash scores. @pyBall/OCL/Assembly.py#22-53
- [emit_configuration](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:54:4-76:29): emits all transformed atom positions for a chosen pose (used to write XYZ frames). @pyBall/OCL/Assembly.py#55-77

### Transform generation utilities
- [super_fibonacci_rotations(N)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:88:0-101:41): low-discrepancy SO(3) sampling. @pyBall/OCL/Assembly.py#89-103
- [generate_transform_buffer(lattice_a, lattice_b, n_rot, n_shift)](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:103:0-150:49): builds rotations, in-plane shifts, applies 6-fold symmetry, tiles 3×3 lattice, flattens to (n_confs, 54, 4, 4). @pyBall/OCL/Assembly.py#104-151
- [pack_transforms](cci:1://file:///home/prokophapala/git/FireCore/pyBall/OCL/Assembly.py:78:0-86:14): packs 3×3 R and 3D shift into 4×4 float32 buffers for the kernel. @pyBall/OCL/Assembly.py#79-88

### CLI driver: [tests/tMMFF/test_assembly.py](cci:7://file:///home/prokophapala/git/FireCore/tests/tMMFF/test_assembly.py:0:0-0:0)
- Args: `--preset`, `--mol`, `--cell`, `--nrot`, `--nshift`, `--shift_range`, `--penalty`, `--zpenalty`, `--zspan_max`, `--clash_max`, `--radius`, `--wg`, `--device`, `--dump`, `--top_k`, `--outdir`. @tests/tMMFF/test_assembly.py#39-58
- Presets: small (HCCOCCH), helicene, triptycene; set lattice if not provided. @tests/tMMFF/test_assembly.py#61-92
- Builds base atoms with radii (vdW from elements table or override). @tests/tMMFF/test_assembly.py#16-32
- Generates transforms (simple or full symmetry), then prefilters by z-span. @tests/tMMFF/test_assembly.py#97-168
- Runs GPU scoring on the filtered set. @tests/tMMFF/test_assembly.py#170-179
- Filters exports by `(zspan < zspan_max) & (clash < clash_max)`; sorts by composite score. @tests/tMMFF/test_assembly.py#182-205
- Outputs (when `--dump`):
  - Pareto plot (clash vs z-span) only for exported configs. @tests/tMMFF/test_assembly.py#213-228
  - Multi-frame XYZ movie of all exported configs with metadata comment line: original index, clash, z-span, rotation matrix, shift. @tests/tMMFF/test_assembly.py#230-245
  - Best-config 2D plot with atoms (dots) and bonds (thin) plus cell overlays. @tests/tMMFF/test_assembly.py#284-296
- Default output dir: `tests/tMMFF/assembly_out_<preset>/`. @tests/tMMFF/test_assembly.py#209-212

## 3) How to run (examples)

### Triptycene preset (GPU)
```bash
python tests/tMMFF/test_assembly.py --preset triptycene --dump
```
- Generates 50k poses (nrot=500, nshift=10), prefilters to ~900 by z-span, GPU run ~0.01 s, exports ~87 configs meeting clash/zspan thresholds, writes outputs into `tests/tMMFF/assembly_out_triptycene/`.

### Helicine with tighter shifts and stricter clash
```bash
python tests/tMMFF/test_assembly.py --preset helicene --dump \
  --nrot 800 --nshift 12 --shift_range 0.4 \
  --zspan_max 10 --clash_max 5
```

### Custom molecule/lattice
```bash
python tests/tMMFF/test_assembly.py \
  --mol path/to/molecule.xyz \
  --cell "lvs 32.7 0 0  16.35 28.319 0  0 0 40" \
  --nrot 600 --nshift 10 --dump
```

## 4) Interpreting outputs
- **Pareto plot**: Each point is a kept config (after clash/zspan thresholds); color = composite score (clash + zpenalty*zspan). The red X marks the best composite score. @tests/tMMFF/test_assembly.py#213-228
- **XYZ movie**: Each frame has a comment line:
  ```
  # idx=<orig_conf> rank=<order> clash=<c> zspan=<z> transform: r00 r01 ... r22 tx ty tz
  ```
  Load in VMD/OVITO/Avogadro to step through plausible packings. @tests/tMMFF/test_assembly.py#230-245
- **Best-config PNG**: Overhead view with atoms as dots, bonds as thin lines, blue unit cell, dashed black 3×3 supercell. @tests/tMMFF/test_assembly.py#284-296

## 5) Scientific/technical notes
- **Penalty function**: Quadratic in overlap depth; scales with how far atoms interpenetrate, not just count of clashes. Early exit at `--penalty`. (Kernel in [cl/Assembly.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/OCL/cl/Assembly.cl:0:0-0:0); host launch @pyBall/OCL/Assembly.py#22-53)
- **Sampling quality**: Super-Fibonacci rotations give near-uniform SO(3) coverage; dense in-plane shifts plus 6-fold symmetry and 3×3 tiling ensure all replicas are considered simultaneously.
- **Performance**: After z-span prefiltering, GPU evaluation is extremely fast (sub-0.1 s for hundreds–thousands of configs). Adjust `--nrot/--nshift` upward freely; prune with `--zspan_max`.

## 6) Suggested experiments for students
1) Vary `--shift_range` (e.g., 0.2 vs 0.5) and observe how compactness affects surviving configs.
2) Tighten `--clash_max` (e.g., 5 → 2) to see how the Pareto front shrinks.
3) Increase `--nrot` and check convergence of best clash/z-span.
4) Modify radii via `--radius` to study sensitivity of packing to vdW assumptions.
5) Inspect the XYZ movie in a viewer and correlate frames with the Pareto positions.
