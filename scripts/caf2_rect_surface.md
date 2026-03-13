# CaF₂ Rectangular Slab Builder (scripts/caf2_rect_surface.py)

A concise tutorial + developer note on how the CaF₂ slab rectangularization script works, how to use it, and how to extend it to other substrates.

## What the script does
- Loads a skew (hex-like) CaF₂ slab from `.xyz` with an `lvs` comment.
- Constructs an **exact rectangular superlattice** using an integer lattice transform (no geometric distortion).
- Selects the **fundamental domain** with a strict half-open rule (periodic, duplicate-free, no missing atoms).
- **Shifts the origin** to the largest empty gap in fractional space so borders avoid atom rows (cleaner edges while staying periodic).
- Trims by detected **layer groups** and replicates to user-requested supercell size.
- Writes an `.xyz` with `lvs` comment and a `.png` preview.

## Core math (periodic, duplicate-free)
1) **Start from fractional coords** `f = r @ inv(L0)` in the original skew cell `L0 = [a, b, c]`.
2) **Choose an integer superlattice** that is orthogonal in-plane. For CaF₂ the best rectangular choice is:
   - `A = a`
   - `B = -a + 2 b`
   - Integer matrix (in the original (a,b) basis):
     ```
     M = [[ 1, 0],
          [-1, 2]]
     det(M) = 2  (area doubles)
     ```
   - Full 3×3 with z untouched: `M = [[1,0,0], [-1,2,0], [0,0,1]]`.
   - New lattice: `L_rect = M @ L0` (rectangular in-plane, exact)
3) **Transform fractional coords into the new supercell**:
   - `g = f @ inv(M)`
4) **Half-open fundamental domain** (exact periodic tiling): keep atoms with
   - `0 <= g_x < 1` and `0 <= g_y < 1` (small tolerance)
   - This guarantees no overlaps and no missing atoms when tiling.
5) **Origin shift to largest empty gaps**:
   - Find largest gap on the circle for `g_x` and `g_y` separately; shift origin to the mid-gap.
   - This moves borders away from atom rows while preserving periodicity.
6) **Layer detection**:
   - Cluster z by tolerance (default 0.12 Å) into planes, then group planes by a median-gap threshold (default 1.5×).
   - Keeps requested top/bottom layer groups without breaking stoichiometry.

## Usage (examples)
From repo root:

```bash
# 2×3 in-plane repeats, keep 4 top layer groups (CaF2_6L_Ni3)
python3 scripts/caf2_rect_surface.py \
    --in_xyz cpp/common_resources/Substrates/CaF2_6L_Ni3.xyz \
    --nx 2 --nz 3 --layers 4

# Slimmer slab: 2×1 repeats, keep top 2 layer groups
python3 scripts/caf2_rect_surface.py \
    --in_xyz cpp/common_resources/Substrates/CaF2_6L_Ni3.xyz \
    --nx 2 --nz 1 --layers 2
```

Outputs land in `cpp/common_resources/Substrates/generated_rect/`:
- `...rect_nx*_nz*_L*_*.xyz` (with `lvs`)
- `...rect_nx*_nz*_L*_*.png` (top/side scatter preview)

Key options:
- `--nx`, `--nz` : in-plane repeats
- `--layers` : number of detected layer groups to keep; `--layers_from top|bottom`
- `--vacuum` : override c-axis vacuum; default keeps original vacuum amount
- `--coeff_search` : search radius for rectangular superlattice (default 4)
- `--no_charges` : drop charges in output xyz
- `--out_dir`, `--out_prefix` : control output placement/naming

## Why edges looked jagged before
A naive approach rotates atoms and clips by Cartesian bounds. That is only approximate to the true fundamental domain; it can leave apparent holes on one side and duplicates on the other. The integer-superlattice + half-open fractional window fixes this; the gap-centered origin improves border aesthetics without breaking periodicity.

## Notes for developers
- Core helpers are in `scripts/caf2_rect_surface.py`:
  - `find_rectangular_supercell` → integer transform
  - `build_supercell_from_transform` → exact fractional selection, gap shift
  - `cluster_z_planes` / `group_planes_to_layers` → layer grouping
- Periodicity proof sketch:
  - Any atom `f` in the original cell maps to `g = f @ inv(M)`; adding an original lattice vector shifts `g` by an integer. Half-open selection `0 <= g < 1` keeps exactly one representative per lattice orbit.
- Origin shift is cosmetic but periodic-safe: shifting by a fractional vector just redefines the origin; we re-wrap into `[0,1)` afterward.

## Extending to other substrates
- The same method applies if you can find an integer in-plane transform that yields orthogonality (or any desired metric). Steps:
  1) Obtain in-plane lattice vectors `(a,b)` from the input `lvs`.
  2) Search small integer combinations `[ [p,q], [r,s] ]` for near-orthogonal rows with nonzero area and positive determinant.
  3) Use that `M` in `build_supercell_from_transform` (or mirror the current search strategy).
  4) Keep the half-open selection and gap-based origin shift.
- If no small orthogonal integer combo exists, you can still pick a nearly rectangular metric, but then the cell will not be exactly commensurate; avoid this for production periodic work.

## Challenges encountered
- **Jagged edges / missing atoms**: caused by Cartesian clipping after rotation. Fixed by exact integer-lattice transform + fractional half-open selection.
- **Duplicate overlaps on tiling**: solved by strict half-open domain.
- **Edge aesthetics**: solved by shifting origin to largest empty gaps in fractional space.
- **Layer stoichiometry**: we group z-planes by gap detection, then trim whole layer groups from top or bottom.

## Outputs to inspect (examples generated)
- `cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz3_L4_top.xyz/.png`
- `cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz/.png`

These are strictly periodic under the rectangular superlattice and have borders shifted off atom rows.

---

# GPU GridFF / PLQ generation notes (CaF₂)

## What we generated
- Full GPU PLQ (Pauli, London, Coulomb) grids for `CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz` using `tests/tMMFF/run_test_GridFF_CaF2.py`.
- Outputs:
  - `tests/tMMFF/data/CaF2_6L_Ni3_rect_nx2_nz1_L2_top/Bspline_PLQd.npy`
  - Plots: `plq_components_gpu.png`, `plq_xy_gpu.png`, `plq_xz_gpu.png`, `plq_linecuts_gpu.png`
  - XSF exports: `gpu_VPaul.xsf`, `gpu_VLond.xsf`, `gpu_VCoul.xsf` (same folder)

## Key fixes and caveats
- **vdW params missing for Ca**: Added a temporary `Ca` entry to `cpp/common_resources/AtomTypes.dat` by copying Na (placeholder until true Ca/F parameters are available). Without this, Pauli/London were zeroed.
- **Parameter source**: GPU helper now loads REQ from `AtomTypes.dat` (not `ElementTypes.dat`).
- **Grid spacing**: Passed `dgx/dgy/dgz` through Morse and Coulomb so the grid stays FFT-friendly and consistent.
- **Slices near surface**: Auto-selects a z-cut ~0.8–3 Å above the topmost atom based on max total variation; prints min/max/range to stdout.
- **Plots with physical units**: Separate `plq_xy_gpu.png` and `plq_xz_gpu.png` use physical extents (Å) and `aspect='equal'` for readability.

## Probe PLQ coefficients (test particle)
- Computed via `getPLQH(R0, E0, a, Q, H)` in `pyBall/tests/utils.py`.
- Defaults (H-like): `R0=1.443 Å`, `E0=0.00190802 eV`, `a=1.5`, `Q=0.4`, `H=0.0`.
- Derived coefficients: `cL = exp(a*R0)*E0`, `cP = exp(a*R0)*cL`, `Q` passes through. Adjust CLI `--R0 --E0 --a --Q` to model another probe.

## Observed magnitudes (near-surface cut, iz≈45)
- VPaul range ≈ 0.56 eV, VLond ≈ 0.24 eV, VCoul ≈ 11.5 eV, VTotal ≈ 4.53 eV.
- Global ranges: VPaul up to ~13 eV, VLond ~2 eV span, VCoul ~277 eV span, VTotal ~110 eV span.

## Known limitations / TODO
- Ca/F vdW parameters are placeholders (copied from Na); needs proper Ca²⁺/F⁻ data for production.
- XSF axis order follows `saveXSF` transpose; visual conventions may appear swapped in some viewers (arrays are correct in Python).
- Orientation cross-check (xs/yz swap suspicion) still pending a formal note—current lvs are diagonal with z normal to slab.

## How to rerun
```bash
cd tests/tMMFF
python3 run_test_GridFF_CaF2.py \
  --dgx 23.175/225.0 --dgy 20.070/200.0 --dgz 48.472/382.0 \
  --job PLQ --use_CG 1
```
Important switches:
- `--z_above_top -1` (default) → auto-pick a near-surface cut with highest total variation.
- `--slice_z` → force a specific z index if desired.
- `--R0/E0/a/Q` → set probe coefficients.
