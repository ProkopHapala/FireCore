# Fireball assemble Hamiltonian reimplementation overview

Here’s a survey of the non–DensityToGrid markdowns under [/fortran/doc](cci:7://file:///home/prokop/git/FireCore/fortran/doc:0:0-0:0), emphasizing how each can guide a McWEDA-focused PyOpenCL rewrite:

1. **CodeStructure.md** — quick reference to Fireball option flags (`itheory`, `itheory_xc`, `iks`, etc.) and module naming conventions. Useful to know which assembler paths must be replicated for the McWEDA branch.

2. **CodeStructure_discussion.md** — describes the “core vs. services” distinction and lists assembler/interaction file families (1c/2c/3c, DOGS, SNXC, etc.). Helps decide which Fortran pieces to mirror first in PyOpenCL (e.g., `assemble_2c`, `assemble_olsxc_off`, `trescentros`).

3. **FireballOCL_implementation_progress.md** — living progress log (dated 2025‑12‑21) detailing parity status for each Hamiltonian term, sign fixes, pending tasks (e.g., 3c L0 mismatch, XC not yet ported). Directly informs what remains for a full McWEDA rewrite.

4. **FireballOCL_kernel_opt.md** — deep back-and-forth about GPU kernel batching, pair sorting, contraction strategies, and fused assembly kernels. It covers data layouts, shared-memory considerations, and is invaluable for designing performant PyOpenCL kernels.

5. **Fireball_API.md** — describes the Python `pyBall.FireCore` entry points, SCF workflow, export helpers (`get_HS_sparse`, `scanHamPiece2c/3c`). Provides the high-level hooks the PyOpenCL code must integrate with.

6. **Fireball_Hamiltonian_Assembly.md** — thorough breakdown of `assemble_mcweda`: term definitions, when each assembler runs (Kscf=1 vs every SCF step), pseudocode for 2c/3c, McWEDA math, and complexity estimates. Core blueprint for matching Fortran logic.

7. **Fireball_code_structure.md** — extended documentation of directories, modules, SCF flow, and data structures (e.g., `interactions.f90`, `integrals.f90`). Clarifies where to source reference behavior and data layouts for a GPU rewrite.

8. **Hamiltonian_Mapping.md** — concise mapping table tying physical terms to Fortran variables, interpolation routines, and SCF workflow pseudocode (geometry vs charge-dependent stages). Handy checklist to ensure no term is omitted in PyOpenCL.

9. **Fireball_general.md / Fireball_general_gemini.md** — high-level explanations of the method, history (Harris → SNXC → McWEDA), accuracy considerations, and conceptual differences from Kohn-Sham DFT. Good for background/context but less actionable technically.

10. **Fireball_equations.md** (not quoted above but listed) likely contains explicit mathematical derivations—worth consulting when implementing kernels that need the exact formulas.

11. **Fireball_Review.md** (also listed) reportedly summarizes the codebase for reviewers; check for architectural notes or pitfalls.

12. **FireballOCL_implementation_progress.md** + **FireballOCL_kernel_opt.md** together function as the “lab book” for OpenCL work—no separate log was found.

13. **ToDo.md** — short list urging definition of a minimal `itheory`/`idogs` path and a prototype Hamiltonian, plus misc. cleanup ideas.

Alignment with current needs:
- **McWEDA specifics**: [Fireball_Hamiltonian_Assembly.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_Hamiltonian_Assembly.md:0:0-0:0), [Hamiltonian_Mapping.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Hamiltonian_Mapping.md:0:0-0:0), and the code-structure docs outline exactly which 2c/3c + XC routines to port, when they run, and their data dependencies.
- **PyOpenCL implementation guidance**: [FireballOCL_implementation_progress.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/FireballOCL_implementation_progress.md:0:0-0:0) and [FireballOCL_kernel_opt.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/FireballOCL_kernel_opt.md:0:0-0:0) capture the state of GPU parity, kernel batching strategies, and remaining gaps (notably 3c and XC).
- **API integration/tests**: [Fireball_API.md](cci:7://file:///home/prokop/git/FireCore/fortran/doc/Fireball_API.md:0:0-0:0) plus the existing verification scripts ([verify_C2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_C2.py:0:0-0:0), etc.) show how to hook the GPU path into Python and validate term-by-term parity.

No additional lab notebook was found beyond these markdowns and the progress/optimization files.

## Summary of `libFireCore.f90` and `pyBall/FireCore.py`:

### Dimension/metadata helpers (non-redundant, all used during sparse export)
1. **`firecore_get_HS_dims`** – returns global dimensions (`natoms`, `norbitals`, `nspecies`, neighbor cap, orbital-cap per atom, `nsh_max`, `ME2c_max`, tensor sizes for `mu`, `xl` beta count, number of species present in Fdata, integerized total electrons). Wrapped by [fc.get_HS_dims()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:438:0-466:5) returning a [FireballDims](cci:2://file:///home/prokop/git/FireCore/pyBall/FireCore.py:355:0-370:62) object. Every test uses it to size buffers. @fortran/MAIN/libFireCore.f90#725-791, @pyBall/FireCore.py#355-372.

2. **`firecore_get_HS_neighs`** – exports per-atom descriptors: `num_orb`, `degelec`, atom species `iatyp`, shell layout (`lssh`, `nssh`, `mu`, `nu`, `mvalue`), Fdata species list `nzx`, neighbor counts and adjacency (`neighn`, `neigh_j`), neighbor block indices (`neigh_b`), and lattice vectors `xl`. Wrapped inside [fc.get_HS_sparse()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:501:0-506:15) (it calls `_get_HS_structures` first to fill these arrays before fetching matrix blocks). Needed once per assembled geometry. @libFireCore.f90 lines continue after 799; wrapper at @pyBall/FireCore.py lines beyond 400 (class [FireballData](cci:2://file:///home/prokop/git/FireCore/pyBall/FireCore.py:372:0-434:70) instantiation).

3. **`firecore_get_HS_sparse`** – copies the actual block matrices `h_mat` and `s_mat` (shape `[iatom, ineigh, inu, imu]` in Fortran order) into provided outputs that match the dims. Python wrapper [fc.get_HS_sparse(dims)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:501:0-506:15) allocates [FireballData](cci:2://file:///home/prokop/git/FireCore/pyBall/FireCore.py:372:0-434:70) buffers and fills them. Used in both [verify_C2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_C2.py:0:0-0:0) and (indirectly) the scan tests when they need reference data. Critical for reconstructing dense matrices. @libFireCore.f90 (following the neighbor export) plus @pyBall/FireCore.py around lines 420+.

4. **`firecore_get_HS_k`** – applies the Fortran `ktransform` to produce dense Bloch-space matrices `H(k)`/`S(k)` given a 3-vector. Not used in the cited tests (they work at Γ using real-space blocks), but handy for band plotting. @libFireCore.f90#912-932, wrapper near @pyBall/FireCore.py lines ~520.

5. **Simple getters** `firecore_get_nspecies`, `firecore_get_nsh_max`, `firecore_get_ME2c_max` – convenience dimension access when callers don’t want the big struct. The scan helpers reuse them to size outputs. Minimal overlap with [get_HS_dims](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:438:0-466:5) but kept for quick queries. @libFireCore.f90#934-958, wrappers near @pyBall/FireCore.py earlier.

### 2-center scan exports
6. **`firecore_scanHamPiece2c` / `_batch`** – call `doscentros` once or over multiple `dR` vectors. They assemble the chosen interaction block (e.g., overlap, kinetic, individual Vna components) optionally applying rotation matrices. Python wrappers [scanHamPiece2c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:595:0-600:17) and [scanHamPiece2c_batch](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:602:0-609:17) allocate `[norb,norb]` or `[npoints,norb,norb]` arrays and pass flattened coordinates. Used in [tests/pyFireball/verify_scan_2c.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_2c.py:0:0-0:0) to sweep radial/ angular distances and compare against PyOpenCL tables. @libFireCore.f90#960-1025 and wrappers @pyBall/FireCore.py#593-610.

### 3-center scan exports
7. **`firecore_scanHamPiece3c` / `_batch`** – full three-center integrals (bcna, bcna-ca, density overlaps, etc.) evaluated by `trescentros` including Legendre reconstruction and Slater-Koster rotation when requested. Used for verifying bcna parity. Wrappers mirror the 2c pattern and are used in [verify_scan_3c.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_3c.py:0:0-0:0) when comparing full rotated matrices. @libFireCore.f90#1027-1328, wrappers @pyBall/FireCore.py#612-695.

8. **`firecore_scanHamPiece3c_raw` / `_batch`** – exports the underlying 2D spline values `bcna_0{1..5}` without angular sums or rotations. Output shape `(5, ME3c_max [, npoints])`. Useful for verifying interpolation alone (tests [verify_scan_3c.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_3c.py:0:0-0:0) uses [scanHamPiece3c_raw](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:626:0-641:14) to debug L0 mismatch). @libFireCore.f90#1081-1207, wrappers @pyBall/FireCore.py#627-662.

9. **`firecore_export_bcna_table`** – dumps an entire raw bcna grid slice (`itheta`, `iME`) into a caller-provided buffer, along with grid spacing `(hx, hy)` and sizes `(nx, ny)`. This is more of a user-facing table export than a runtime scan; seldom used in tests but allows offline inspection of Fdata. Wrapper returns reshaped arrays. @libFireCore.f90#1209-1271, wrapper @pyBall/FireCore.py#618-685.

### Other small helpers
10. **`firecore_scanHamPiece3c_batch`** (full rotated) already covered; no redundant alternative except raw vs full.  
11. **`firecore_scanHamPiece3c` (single)** vs `_batch` – same computation, just vectorized vs scalar. Tests mostly prefer single-call per geometry, but the batched variant is necessary for high-throughput scanning or GPU parity checks.  
12. **`firecore_scanHamPiece2c_batch`** – same reasoning.

### Overlap / redundant considerations
- `firecore_get_HS_dims` subsumes `get_nspecies`, `get_nsh_max`, `get_ME2c_max`; the latter are merely shortcuts.
- [scanHamPiece3c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:619:0-624:20) vs [scanHamPiece3c_batch](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:685:0-693:20) vs raw variants are complementary (full vs raw vs batched). Not redundant because raw exports different data (pre-rotation, per-L channel).
- `firecore_get_HS_neighs` and `firecore_get_HS_sparse` must both be called sequentially; the first conveys topology and species mapping required to interpret block matrices the second delivers. Neither replaces the other.

### How tests use them
1. **[tests/pyFireball/verify_C2.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_C2.py:0:0-0:0)**
   - Calls [fc.get_HS_dims()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:438:0-466:5) to size buffers.
   - Invokes [fc.get_HS_sparse(dims)](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:501:0-506:15) (which internally called [get_HS_neighs](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:470:0-497:15) and [get_HS_sparse](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:501:0-506:15)).
   - Uses the exported neighbor lists to reconstruct dense matrices for comparison with PyOpenCL. It also toggles [set_options](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:589:0-590:97) to isolate individual Hamiltonian pieces for each call. File @tests/pyFireball/verify_C2.py#97-165.

2. **[tests/pyFireball/verify_scan_2c.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_2c.py:0:0-0:0)**
   - Uses [fc.scanHamPiece2c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:595:0-600:17) to generate reference blocks for various distances/angles and compares them to PyOpenCL interpolation results. It also leverages `firecore_get_nspecies`, `firecore_get_ME2c_max`, etc., implicitly via [fc._get_norb](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:535:0-542:23) when sizing arrays. File lines @tests/pyFireball/verify_scan_2c.py#1-200 (partial view shows [_extract_components](cci:1://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_2c.py:11:0-32:31); later sections loop over interactions calling [fc.scanHamPiece2c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:595:0-600:17)).

3. **[tests/pyFireball/verify_scan_3c.py](cci:7://file:///home/prokop/git/FireCore/tests/pyFireball/verify_scan_3c.py:0:0-0:0)**
   - Calls both [fc.scanHamPiece3c_raw](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:626:0-641:14) (for raw table parity) and [fc.scanHamPiece3c](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:619:0-624:20)/`_batch` to validate full rotated bcna blocks. It also uses the batched raw variant to profile entire grids. File @tests/pyFireball/verify_scan_3c.py (not reread but previously summarized).

### Proper usage sequence (typical workflow)
1. Initialize Fortran (preinit → set lattice → init).
2. Assemble Hamiltonian (or SCF) to populate `h_mat`/`s_mat`.
3. Call [get_HS_dims()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:438:0-466:5) once to know sizes.
4. Call [get_HS_sparse()](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:501:0-506:15) which internally uses [get_HS_neighs](cci:1://file:///home/prokop/git/FireCore/pyBall/FireCore.py:470:0-497:15) to fill structural arrays, then fetches H/S (and optionally `rho`) blocks.
5. Build Python-side dense matrices or pass neighbor-lists to GPU assembly.
6. For diagnostics: use `scanHamPiece2c/3c` to produce reference blocks for specific geometries, or raw scans for integral-table validation.
7. Use `scanHamPiece*_batch` when you need many points efficiently, typically for verifying PyOpenCL vectorized kernels.

No stored memories influenced this response.

## Summary of DensityToGrid documentation

Summary of DensityToGrid documentation
--------------------------------------
1. **SparseDenmatToGriod_Fotran_to_pyOpenCL.md** — concise “how-to” for exporting Fireball sparse H/S/ρ blocks directly into pyOpenCL without flattening, plus current GPU task-builder status and TODOs for real Fireball radial basis replacement and optional CSR improvements. Sections detail buffer shapes, host plans, kernel indexing, and testing workflow @/home/prokop/git/FireCore/fortran/doc/DensityToGrid/SparseDenmatToGriod_Fotran_to_pyOpenCL.md#1-87.
2. **DensityProjection_TDD_Plan.md** — a milestone-style roadmap (geometry harness → CPU DM/Grid prototypes → GPU DM → GPU projection → tests). It lists references into existing Fortran/Python APIs and test scripts, defines required data exports, algorithms, and validation strategy @/home/prokop/git/FireCore/fortran/doc/DensityToGrid/DensityProjection_TDD_Plan.md#1-94.
3. **DensToGrid_Opt.md** — the longest “lab notebook”: back-and-forth design discussions covering interpolation choices (cubic B-splines vs texture), hierarchical task building, work-group batching, DM construction strategies, and load balancing. It effectively captures evolving decisions and rationale for kernel architecture and is the primary narrative record @/home/prokop/git/FireCore/fortran/doc/DensityToGrid/DensToGrid_Opt.md#1-630 (continues beyond line 630 with further refinements).
4. **DensityToGridPartitioning.md** — focuses on spatial partitioning, load balancing, and a 2D benzene prototype to reason about block workloads. It also introduces hierarchical “macro-grid + block” filtering and provides pseudocode plus plotting ideas for overlap statistics @/home/prokop/git/FireCore/fortran/doc/DensityToGrid/DensityToGridPartitioning.md#1-415.
5. **GridCoreRefinement.md** — explores advanced ideas (PAW-like blunted basis, overlapping grids, pyramidal plugs, custom quadrature weights, optimized finite-support basis functions). It is broader R&D brainstorming rather than immediate implementation guidance @/home/prokop/git/FireCore/fortran/doc/DensityToGrid/GridCoreRefinement.md#1-400 (continues further).

Collectively, **DensToGrid_Opt.md** + **DensityProjection_TDD_Plan.md** serve as the most actionable references for rewriting the Hamiltonian assembly/projection logic in pyOpenCL, while the other files capture specialized research threads.

McWEDA Fortran assembler context
--------------------------------
[assemble_mcweda.f90](cci:7://file:///home/prokop/git/FireCore/fortran/ASSEMBLERS/assemble_mcweda.f90:0:0-0:0) orchestrates the full Fireball SCF assembly: neighbor-map refresh, Ewald energy, 1c/2c/3c contributions (including VNA, VNL, OLS/SN XC, DOGS terms), and finally `buildh`. Any pyOpenCL rewrite must replicate the staged sequence (neighbors → per-term assembly → build H) or at least ensure equivalent data exports from existing Fortran routines @/home/prokop/git/FireCore/fortran/ASSEMBLERS/assemble_mcweda.f90#1-366.

Existing pyBall verification harnesses
--------------------------------------
1. **verify_H2.py** and **verify_C2.py** compare Fortran vs PyOpenCL overlap, kinetic, Vna, and Vnl matrices by reconstructing dense matrices from sparse exports and by invoking `pyBall.FireballOCL.OCL_Hamiltonian`. They demonstrate how to fetch neighbors via `fc.get_HS_sparse` and feed them into the Python assembler @tests references:  
   - @/home/prokop/git/FireCore/tests/pyFireball/verify_H2.py#1-133  
   - @/home/prokop/git/FireCore/tests/pyFireball/verify_C2.py#1-177
2. **verify_scan_2c.py** sweeps radial/angular separations for 2-center tables, extracting individual orbital components to spot discrepancies. It shows how to call `fc.scanHamPiece2c` vs OCL tables, useful for regression once kernels move to pyOpenCL @/home/prokop/git/FireCore/tests/pyFireball/verify_scan_2c.py#1-234.
3. **verify_scan_3c.py** (and **verify_Fdata_read_3c.py**) target 3-center bcna tables, either via full Legendre reconstruction or raw table access, again cross-validating Fortran exports with the Python loader @/home/prokop/git/FireCore/tests/pyFireball/verify_scan_3c.py#1-235 and @/home/prokop/git/FireCore/tests/pyFireball/verify_Fdata_read_3c.py#1-85.

PyOpenCL assembly components
---------------------------
- **OCL_Hamiltonian.py** and **FdataParser.py** underpin the current PyOpenCL-side Hamiltonian builder. They already parse Fdata, prepare splines/tables, and expose `assemble_full`, `assemble_2c`, `scanHamPiece2c/3c`, etc., which the tests above leverage.
- **Grid.py** (plus `pyBall/pyocl_dft/density.py`) implement task builders and projection scaffolding aligned with the documented strategies, including GPU and CPU task-building, B-spline packing, and tiled kernels.

Lab-book / progress notes
-------------------------
Besides the markdown files above, no dedicated “labbook” file was found elsewhere in the repo. The conversational structure inside **DensToGrid_Opt.md**, **DensityToGridPartitioning.md**, and **GridCoreRefinement.md** effectively functions as the running progress log of design decisions for Density-to-Grid and Hamiltonian porting work.

Next steps (if desired)
-----------------------
1. Normalize the documentation: extract key actionable requirements from DensToGrid_Opt.md into an ordered implementation checklist alongside the TDD plan.
2. Ensure PyOpenCL components (`OCL_Hamiltonian`, `GridProjector`, etc.) cover every term present in `assemble_mcweda` so future rewrites don’t miss contributions.
3. Expand the test suite by wiring the DensityProjection TDD milestones into pytest to keep GPU vs Fortran parity as kernels evolve.

Memories applied
----------------
No stored user memories were relevant to this response.