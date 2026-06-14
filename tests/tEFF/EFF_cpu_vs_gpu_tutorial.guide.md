

Here’s a concise user tutorial to use eFF on GPU, regenerate electron .xyz files, run scans, and perform detailed debug tests. I’ve reviewed [tests/tEFF/EFF_cpu_vs_gpu_discussion.md](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/EFF_cpu_vs_gpu_discussion.md:0:0-0:0) and [tests/tEFF/check_EFF_forces_cpu_vs_gpu.md](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/check_EFF_forces_cpu_vs_gpu.md:0:0-0:0) and the scripts in `tests/tEFF/`.

## 0) Prereqs
- Env: `PYTHONPATH=/home/prokop/git/FireCore`
- Build CPU lib once: `cd cpp/Build && cmake .. && make -j` (or rely on regen scripts that auto-build).
- GPU: PyOpenCL + OpenCL driver working; small-system kernel uses one workgroup per system.

## 1) Generate electron-containing .xyz from atom-only .xyz
Tool: [tests/tEFF/regen_scans.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/regen_scans.py:0:0-0:0) (writes `na,ne,core` header via `bOutCoreHeader=True` so GPU loader parses it).

Examples (atom-only bases live in [tests/tEFF/export/scan_data/](cci:9://file:///home/prokop/git/FireCore/tests/tEFF/export/scan_data:0:0-0:0)):
```bash
cd tests/tEFF
# Regenerate all variants (spins_a/spins_fc/spins_e, pairs_a/pairs_fc/pairs_e) from atom-only bases
PYTHONPATH=../.. python3 regen_scans.py --backup --in-dir export/scan_data --out-dir export/scan_data
```
Notes:
- Skips files with `__` (already electron-augmented).
- Skips [distscan_Oe.xyz](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/export/scan_data/distscan_Oe.xyz:0:0-0:0) (not atom-only, contains pseudo-atom E).
- Uses `nstepMax=0` to avoid electron collapse; `bChangeEsize=True` to init sizes.
- Modes encoded: `coreMode='a'` → all-electron; `coreMode='f'` → frozen core; `_e` variants generate all electrons explicitly.

Output files: e.g. [distscan_CH4__spins_a.xyz](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/export/scan_data/distscan_CH4__spins_a.xyz:0:0-0:0), [distscan_CH4__pairs_fc.xyz](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/export/scan_data/distscan_CH4__pairs_fc.xyz:0:0-0:0), etc.

## 2) Run scans and produce GPU vs CPU plots (1D/2D, rigid or relaxed)
- **Rigid (no relaxation):** `tests/tEFF/batch_eval_plot_cpu_gpu_maps.sh` → [eval_plot_cpu_gpu_maps.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/eval_plot_cpu_gpu_maps.py:0:0-0:0) computes CPU/GPU energies directly.
  ```bash
  cd tests/tEFF
  OUTDIR=export/plots_cpu_gpu_regen PYTHONPATH=../.. OFFLOAD_CORE=0 bash batch_eval_plot_cpu_gpu_maps.sh
  ```
  Outputs: `*_cpu.npy`, `*_gpu.npy`, `*__1d.png`/`*__2d.png`, `*__stats.txt`.

- **Relaxed 1D scans (fixed ions, electron relax, parity check):**
  - Runner: `tests/tEFF/run_relaxed_scans_1d.sh` (H2/CH4/H2O distance scans; dt=0.01, damping=0.1, steps=2000, fixed ions).
  - Plotter: `tests/tEFF/plot_scan_parity.py` (CPU dotted lw=1.5, GPU solid lw=0.5; per-component diffs).
  ```bash
  cd tests/tEFF
  bash run_relaxed_scans_1d.sh
  ```
  Outputs: `export/scan_parity_*/*{Es5_cpu.npy,Es5_gpu.npy,scan_parity.png}` plus reports.

- **Relaxed 2D scans (spins/pairs, frozen core):**
  - Runner: `tests/tEFF/run_relaxed_scans_2d.sh` auto-discovers `*scan*__spins_fc.xyz` / `*scan*__pairs_fc.xyz`, relaxes via `run_relax_parity_protocol.py` (dt=0.01, damping=0.1, steps=2000), then plots maps using `eval_plot_cpu_gpu_maps.py --from-es5-*`.
  ```bash
  cd tests/tEFF
  bash run_relaxed_scans_2d.sh
  ```
  Outputs: `export/relaxed_*/*{Es5_cpu.npy,Es5_gpu.npy,__stats.txt,__1d.png/__2d.png}`.

Files consumed: any `*__*.xyz` with `na,ne,core` header in [export/scan_data/](cci:9://file:///home/prokop/git/FireCore/tests/tEFF/export/scan_data:0:0-0:0).

## 3) Detailed debug tests with prints (few problematic examples)
Key scripts in `tests/tEFF/`:
- [test_ocl_vs_cpu.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/test_ocl_vs_cpu.py:0:0-0:0): CPU vs GPU force comparison on small systems; uses `localMD` small-system kernel; tol defaults to 1e-4. Enable debug prints via kernel compile flags in `eFF.cl` (DBG toggles) and CPU verbosity in [eFF.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/eFF.h:0:0-0:0) ([setVerbosity](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:101:0-107:1)).
- [debug_singlepoint_cpu_gpu.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/debug_singlepoint_cpu_gpu.py:0:0-0:0): single geometry energies CPU vs GPU (can pass dbg_pair/idbg_sys/offload_core).
- [eval_plot_cpu_gpu_maps.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/eval_plot_cpu_gpu_maps.py:0:0-0:0): batch energy eval for scans; supports `--worker cpu|gpu`, `--out-npy`, `--nloc`, `--device`.

Typical debug flow:
```bash
# Single geometry CPU vs GPU with optional pairwise debug
PYTHONPATH=../.. python3 debug_singlepoint_cpu_gpu.py \
    --xyz export/scan_data/distscan_CH4__spins_fc.xyz \
    --dbg_pair 1 --idbg_sys 0 --offload_core 0
```
For force parity:
```bash
PYTHONPATH=../.. python3 test_ocl_vs_cpu.py --xyz H2O_fixcore.xyz --tol 1e-4
```
(see [check_EFF_forces_cpu_vs_gpu.md](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/check_EFF_forces_cpu_vs_gpu.md:0:0-0:0) for debug print formats and known fixes).

## 4) Relevant functions and kernels (what each does)

### C++ (CPU) – files: [cpp/common/molecular/eFF.h](cci:7://file:///home/prokop/git/FireCore/cpp/common/molecular/eFF.h:0:0-0:0), [cpp/libs/Molecular/eFF_lib.cpp](cci:7://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:0:0-0:0)
- Energy pieces:
  - `evalKinetic` (electron kinetic)
  - `evalEE` (e–e Coulomb + Pauli)
  - `evalAE` (e–atom Coulomb + Pauli; core Coulomb if `bCoreCoul`)
  - `evalAE_ECP` (ECP branch if `bUseECPs`)
  - `evalAA` (atom–atom Coulomb)
  - `evalCoreCorrection` (adds frozen-core self terms if enabled)
  - [eval](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:109:0-109:34) orchestrates based on flags
- Modes/toggles:
  - `setCoreMode(char coreMode)`: `'a'` all-electron, `'f'` frozen core, `'e'` ECP
  - Flags on `eFF`: `bEvalKinetic`, `bEvalEE`, `bEvalCoulomb`, `bEvalPauli`, `bEvalAE`, `bEvalAECoulomb`, `bEvalAEPauli`, `bEvalAA`, `bEvalCoreCorect`, `bUseECPs`, `iECPmodel`
- Input/build:
  - `from_xyz_line`, `loadFromFile_xyz`: parse `na,ne,core` .xyz with e+/e-/e2 labels
  - [preAllocateXYZ](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:588:0-636:1), [processXYZ](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:638:0-757:1) (generate electrons from atom-only; optional valence/core)
  - [processXYZ_e](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:761:0-826:1) (explicit electrons with `na,ne,core` header)
  - [builder2EFFstatic](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:522:0-586:1) (construct electrons, add core e2 when requested)
- Scripts:
  - [run_process_xyz.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/run_process_xyz.py:0:0-0:0), [run_process_xyz_e.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/run_process_xyz_e.py:0:0-0:0), `run.sh`

### OpenCL (GPU) – files: `cpp/common_resources/cl/eFF.cl`, `pyBall/OCL/eFF_ocl.py`
- Kernel: `localMD` small-system WG-local kernel (one WG per system; ntot≤64). Does kinetic, e–e, e–ion, ion–ion. Modes:
  - ECP: not implemented (Gaussian Coulomb+Pauli only)
  - Frozen core: via atom params (`Zcore_eff`/`sP`); `bFrozenCore` arg passed but handled through params
  - Pairing: electrons are individual entries; `e2` in XYZ duplicated into ± spins in loader
- Host harness (`eFF_ocl.py`):
  - Loads multi-geometry `na,ne,core` XYZ; recognizes `e2`/`e+`/`e-`
  - Allocates batch buffers; one WG per system
  - `eval_energies_localmd` / `relax_systems` launch kernels
  - Supports optional core-constant offload and core map (recent work)
- Debug aids:
  - Kernel printf toggles for pairwise terms; CPU verbosity via [setVerbosity](cci:1://file:///home/prokop/git/FireCore/cpp/libs/Molecular/eFF_lib.cpp:101:0-107:1)
  - Force output buffer (`fout`) for comparisons ([test_ocl_vs_cpu.py](cci:7://file:///home/prokop/git/FireCore/tests/tEFF/test_ocl_vs_cpu.py:0:0-0:0))

## 5) Quick recipes

### Regenerate all electron scans (no relaxation, correct headers)
```bash
cd tests/tEFF
PYTHONPATH=../.. python3 regen_scans.py --backup --in-dir export/scan_data --out-dir export/scan_data
```

### Generate jittered single-config inputs (H2 / CH4) for trajectory parity
Tools: `tests/tEFF/export/h2_jitter.py`, `tests/tEFF/export/ch4_jitter.py` (sample names; use the provided jitter scripts next to the XYZs). They write `*jpos*` jittered XYZs used by long-parity tests.

Example:
```bash
cd tests/tEFF/export
python3 h2_jitter.py   # writes h2_jpos001_js002.xyz etc.
python3 ch4_jitter.py  # writes ch4_jpos001_js002.xyz etc.
```
Then run long-parity (CPU C++ ialg=3 vs GPU localMD) via `run_relax_parity_protocol.py` with `--long-parity --fix-atoms` pointing to the jittered file.

### Run full batch plots/stats (GPU vs CPU, no offload)
```bash
cd tests/tEFF
OUTDIR=export/plots_cpu_gpu_regen PYTHONPATH=../.. OFFLOAD_CORE=0 bash batch_eval_plot_cpu_gpu_maps.sh
```

### Single-point CPU vs GPU with debug pairs
```bash
PYTHONPATH=../.. python3 debug_singlepoint_cpu_gpu.py \
    --xyz export/scan_data/distscan_CH4__spins_fc.xyz \
    --dbg_pair 1 --idbg_sys 0 --offload_core 0
```

### Force parity test
```bash
PYTHONPATH=../.. python3 test_ocl_vs_cpu.py --xyz H2O_fixcore.xyz --tol 1e-4
```