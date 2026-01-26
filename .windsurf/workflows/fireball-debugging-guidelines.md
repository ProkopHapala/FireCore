---
description: Evidence-based Fireball Fortran↔OpenCL parity and debugging
---

## Purpose
Use this workflow when porting or verifying Fireball terms (S/T/Vna/dipole/Vca/Ewald/Vxc) from Fortran to OpenCL. It encodes the proven strategy, guardrails, and pitfalls.

## Rules (must follow)
1) **Fortran is source of truth.** If outputs differ, OpenCL is wrong.
2) **Fortran instrumentation:**
   - Add only gated, clearly commented debug prints inside Fortran loops (e.g., `! DEBUG OCL PARITY START/END`).
   - Gate by verbosity/atom indices; never change physics/logic.
   - Avoid new modules/C-bindings just to expose scalars.
3) **Neighbor identity:** use (iatom, ineigh, mbeta) as unique key; do not collapse by (iatom, jatom); respect periodic shifts.
4) **Table encoding:** SK path only for SK tables (S/T/Vna…); generic recover+rotate for non-SK (dipole_*). Drive by per-pair flag.
5) **Rotations:** match Fortran epsilon/twister inputs per term (e.g., vna_atom may use absolute r2). Don’t reuse a generic frame blindly.
6) **SFIRE overwrite:** replicate 3c overwrite of reverse slots (transpose/conjugate), not just additive accumulation.
7) **Density/charges:** when Fortran uses Qin (itheory=1), compare against Qin-weighted results; make the comparison target explicit.
8) **Structural vs numeric:** keep structural/self-consistency checks separate from float parity; expect exact zeros only in structural checks.

## Workflow
1) **Setup references**
   - Export Fortran sparse data (neigh, s_mat, rho/rho_off, etc.) and identify required tables.
   - If needed, insert gated Fortran prints for inputs/outputs around the target loop; mark with `! DEBUG OCL PARITY` comments.

2) **Scan microtests**
   - For each table, run `scanHamPiece` vs OCL scan across multiple distances/angles with matching rotation flags.
   - Confirm SK vs generic path matches Fortran (dipole_* are generic).

3) **Structural integrity**
   - Verify neighbor indices (neigh_self/neigh_back) exactly (integer compare).
   - Reconstruct pair lists using (iatom, ineigh, mbeta); ensure periodic shifts are preserved.

4) **Block-level parity**
   - Compare 2c/3c blocks from OCL vs Fortran exports using the exact neighbor slots; handle SFIRE overwrite for 3c.
   - Keep structural checks separate from numeric summaries.

5) **Dense reassembly**
   - Rebuild dense matrices from blocks and compare to Fortran dense references (respect orientation nu,mu vs mu,nu; SFIRE overwrites).

6) **End-to-end term**
   - After inputs match (S/dip/rho), validate composed terms (EwaldSR/LR, Vca, future Vxc) against Fortran outputs or known zeros.

7) **Debug drill when mismatches appear**
   - Bisect: scan → neighbor-slot blocks → dense assembly → composed term.
   - Use worst-diff reporting; if needed, add/adjust gated prints in Fortran (never alter logic).

## Pitfalls to avoid
- Forcing SK on non-SK tables (dipole_*).
- Ignoring mbeta/periodic shifts; collapsing by (iatom, jatom).
- Mixing structural checks into numeric summary; zero diffs become misleading.
- Guessing rotation frames; always copy Fortran’s epsilon inputs.
- Comparing avg_rho to neutral weights when Fortran uses Qin.
- Skipping scan tests and jumping to end-to-end.

## When to stop and fix prerequisites
- If any scan test fails, fix before neighbor-slot tests.
- If neighbor mapping fails, do not proceed to numeric checks.
- If inputs (dip/S/rho) don’t match, do not attempt EwaldLR/Vxc integration.
