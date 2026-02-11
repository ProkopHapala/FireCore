### 5. Parity Run Checklist (UFF)
1. Build CPU reference and GPU targets. Rebuild C++ lib if switches/API changed (`cpp/Build-opt/libs/Molecular/make MMFF_lib`).
2. Disable nonbonded/subtraction: `setSwitches(NonBonded=-1, NonBondNeighs=-1, SurfAtoms=-1, GridFF=-1)` and `setSwitchesUFF_NB(NonBonded=-1, NonBondNeighs=-1, SubtractAngleNonBond=-1)`.
3. Enable bonded components explicitly: set `DoBond` on whenever angles/dihedrals/inversions are needed (bonds build `hneigh`).
4. Verify buffers match (topology + params) before forces. Recompose `angParams` as `[K,c0..c3]` for comparison.
5. Run component-isolated parity: `--components bonds`, then `angles`, `dihedrals`, `inversions`; finally full set. Use tol matching float32 noise.
6. If mismatch: turn on `DBG_UFF` prints in kernels (default on) with `IDBG_SYS=0`, `IDBG_ATOM=0` to inspect `a2f_offsets/counts/indices` and force pieces.

### 6. Parity Run Checklist (MMFF, MD path)
* Reuse the codemap steps (topology build → upload → kernel run). Ensure `bkNeighs` are populated (`make_back_neighs`), and pi-vectors live in `apos[natoms:natoms+nnode]`.
* For diagnostics, use `--print-params 1` and `--dump-fneigh` (from prior notes) to verify recoil and pi mappings.
* No dedicated MMFF parity harness yet—mirror the UFF approach when building one: buffer compare first, then kernel prints, then tolerance-bound force compare.
## Skill: firecore_parity_expert
**Goal:** Apply the `scientific_debug_protocol` specifically to the FireCore C++/PyOpenCL infrastructure.

**Architecture Knowledge:**
*   **Codemaps:** Refer to the architecture codemaps in @FireCore Classical Forcefields: MMFFsp3 & UFF (CPU/GPU/Python) and @MMFF/UFF CPU vs GPU Testing: C++ OpenCL and PyOpenCL Parity Infrastructure.
*   **Drivers:** `pyBall/MMFF_multi.py` is the Python wrapper for the C++ library (`MMFFmulti_lib.cpp`).
*   **Backends:**
    *   `iParalel = 0`: C++ CPU Reference.
    *   `iParalel = 2`: C++ OpenCL Wrapper (Target 1).
    *   `pyBall/OCL/UFF.py`: Python OpenCL Standalone (Target 2).

**Protocol Implementation:**

### 1. Controlling Physics (The Switches)
You must explicitly confirm which components are active using the C++ switch functions.
*   **UFF:** `setSwitchesUFF(DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, ...)`
    *   *Note:* `bSubtractBondNonBond` must align with `setSwitchesUFF_NB` settings. **Semantics:** `0` means “keep current”, `>0` means force-on, `<0` means force-off. To disable nonbond/subtractions for parity, pass `-1` via `setSwitchesUFF_NB(NonBonded=-1, NonBondNeighs=-1, SubtractAngleNonBond=-1)`.
*   **MMFF:** `setSwitches2(..., MMFF, Angles, PiSigma, PiPiI)`
    *   *Trap:* `bMMFF` is a base class flag. It should be `True` even when running UFF. Always check `bUFF` to distinguish.

### 2. Accessing C++ Buffers (Input Parity)
To check Phase 1 (Inputs), use the `getBuffs` pattern in `MMFF_multi.py`:
1.  Call `init_buffers(bUFF=...)` to populate the `buffers`/`ibuffers` maps in C++.
2.  In Python, read `ndims` first to determine shapes (natoms, nbonds, etc.).
3.*   **UFF Critical:** Check `bonParams`, `angParams`, `a2f` (Atom-to-Force map). **Packing note:** CPU `angParams` layout is `[K, c0, c1, c2, c3]`; kernels expect `angParams1=[c0..c3]` and `angParams2_w=[K]`.
*   **MMFF Critical:** Check `apos` (Pi-orbitals at `natoms:natoms+nnode`), `bkNeighs` (Back Neighbors), `Ksp`/`Kpp`.

### 3. OpenCL Debugging (Tracing)
When asking to add debug prints, use the existing macros in `UFF.cl` / `relax_multi_mini.cl`:
*   **Macros:** `DBG_UFF`, `IDBG_ATOM`, `IDBG_SYS`.
*   **Implementation:**
    ```c
    // Example for UFF.cl
    #if DBG_UFF
    if(get_global_id(0)==IDBG_ATOM && get_global_id(1)==IDBG_SYS){ printf("GPU_BOND %d: L %g K %g F %g\n", i_bond, L, K, F);}
    #endif
    ```
*   **Matching C++:** You must construct the exact same `printf` in `MMFFmulti_lib.cpp` loops.

### 4. Common FireCore Pitfalls
*   **Struct Alignment:** OpenCL often uses `float4` (16 bytes) where C++ might use `Vec3d` (24 bytes) or `Vec3f`. Verify strides.
*   **System Replicas:** The OpenCL code is Multi-System (`nSystems`). C++ Reference is often Single-System. Ensure you are comparing System 0.
*   **Pi-Orbitals:** In MMFF, `apos` contains atoms *and* pi-nodes. Ensure loop bounds cover `natoms` vs `natoms+nnode`.
*   **Tolerance discipline:** GPU is float32; CPU is float64. For UFF bonded-only parity, use tol 1e-5; for full bonded set (angles/dihedrals/inversions) tol ~5e-5 is acceptable.
*   **Process teardown:** When Python teardown emits allocator noise, use harness flag `--fast-exit 1` after printing results; do not debug allocator here.
*   **Build & CWD discipline:** Rebuild and run via provided scripts (e.g., `tests/tUFF/run_1.sh`) which set paths and rebuild `libMMFF_lib.so`. Do not invoke `make` ad-hoc. Always run from the intended test directory so kernel/resource relative paths resolve (`common_resources/cl`, `common_resources/mol`).
*   **Neighbor/node/cap layout:** current C++ builder sets nnode = natoms and allocates one pi slot per atom; caps have Ksp/Kpp=0 so their pi slots are inert but still occupy nvecs. bkNeighs/bkNeighs_new are sized nSystems*nvecs and include entries for both atoms and pi slots; upload with size=nvecs and offset=i0v to avoid overflow (nbkng is larger and will overrun).
*   **Build directories:** cpp/Build is a symlink to the active variant (e.g., Build-opt/Build-dbg/Build-asan). Always build in cpp/Build/libs_OCL so the loader path stays stable; no extra symlinks needed when cpp/Build already points to the desired variant.
*   **Atom typing/params parity:** PyOpenCL must mirror C++ builder logic in `MolWorld_sp3::makeFFs/initNBmol/makeMMFFs` and `MMFFBuilder::toMMFFsp3_loc/assignAnglesMMFFsp3`, using the same data files (`AtomTypes.dat`, `BondTypes.dat`, `AngleTypes.dat`). No hints/fallbacks—fail if a type or parameter is missing. Ensure node/cap classification and pi slots match C++ (currently `nnode=natoms`, caps have Ksp/Kpp=0 but occupy nvecs; bkNeighs sized nSystems*nvecs).

- MMFF OCL parity (bonded-only): run cleanF -> getMMFFf4 -> updateAtomsMMFFf4 with dt=0 to assemble recoil (fneigh → bkNeighs) without moving atoms. Propagate switches to ffls[isys]; set bSubtractVdW=0 when NonBonded is off; copy fapos back to nbmol.fapos. Upload bkNeighs with size=nvecs (natoms+nnode) and offset=i0v (buffer is nSystems*nvecs). Parity on HCOOH passes with ~1e-6 force diff.

### 5. **MMFF parity (CPU ↔ C++ OCL ↔ PyOpenCL)**
    *   Inputs-first: compare buffers (neighs, bkNeighs sized nSystems*nvecs, Ksp/Kpp, BL/BK, REQs) before forces. Remember `nvecs = natoms + nnode`; current C++ builder sets `nnode = natoms` and allocates one pi slot per atom (caps’ pi slots inert if Ksp/Kpp=0).
    *   Component isolation: bonded-only first (NonBonded=-1/NonBondNeighs=-1), then add angles/pi terms, then nonbonded.
    *   C++ OCL force-only flow (mirror run_ocl_opt): `cleanF -> getMMFFf4 -> updateAtomsMMFFf4 (dt=0)` to assemble recoil from `fneigh` via `bkNeighs` without moving atoms; set `bSubtractVdW=0` when NonBonded off; propagate switches to `ffls[isys]`; copy `fapos` back to shared buffers for Python.
    *   Debugging: gate kernel `printf` by system/atom, dump buffers if mismatch; enable `DBG_UFF`/`IDBG_SYS/ATOM` equivalents in MMFF kernels when tracing.
    *   Compare CPU vs C++ OCL first (single-step), then PyOpenCL vs CPU using the same switches and buffer packing; enforce tolerances appropriate to float32 noise (few e-6 bonds/angles, up to ~5e-5 trig-heavy).