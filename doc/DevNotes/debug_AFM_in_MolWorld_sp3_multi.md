# USER

I'm running AFM simulations in @MolWorld_sp3_QMMM_multi.h 
@MolGUI.h 

but there is something strange, it looks like there is only repulsion not attraction, what it can be.
Morse potential should be attractive, pauli potnetial should be repulsive, and electrostatic should depand on the charges. But the results I see there seems like the molecule is just repulsive. Or maybe there are some non-sense / corrupted params of atoms usd for AFM imaging. Plese do not try to code anything, just review the code related to AFM imaging and tell what could be worng. 

Also we plot some grid with some property ploted on it. For diagnostics it would be good if we can plot always one of the 3 component (Pauli, vdW, electrostatics) separately on that grid, to avoid confusion (without relaxation). But we need some way how to switch between the view (df relaxed, Pauli, vdW, eelectrostatics - just the force field z) Pehaps some GUI combo/selection box?@GUI.h 

can you figure out how to do it? Do not code yet.

---

# SWE-1.6

## Issue: AFM Force Field Uses Combined LJ Potential Instead of Separate Components

**Root cause in [relax_multi.cl](cci:7://file:///home/prokop/git/FireCore-master/cpp/common_resources/cl/relax_multi.cl:0:0-0:0) PPAFM_makeFF kernel (lines 2963-2967):**

The AFM kernel uses `getLJQH()` which computes a **combined Lennard-Jones + Coulomb** potential:
```cpp
fe += getLJQH(dp, REQ, 1.0f);  // Line 2963
```

The `getLJQH` function (lines 151-164) computes:
- **Lennard-Jones**: `E = E0*((R0/r)^6 - 2*(R0/r)^12)` 
  - The `(R0/r)^12` term is repulsive (Pauli-like)
  - The `(R0/r)^6` term is attractive (London dispersion)
- **Coulomb**: `Ec = Q1*Q2/sqrt(r^2+R2damp)`

**Problems:**
1. **No separate Pauli/London/Coulomb components** - The kernel doesn't use the GridFF infrastructure that stores these separately (`Bspline_PLQ` with Pauli, London, Coulomb as Vec3d)
2. **Hardcoded tip parameters** in `OCL_MM.h` line 152: `tipParams{1.661f, 0.0091063f, 0.0f, 0.0f}` modify the REQ values (lines 2955-2956), which may be inappropriate for your system
3. **Potential parameter corruption** - If the molecule's REQ parameters (R0, E0, Q) are wrong or corrupted, the LJ potential could appear purely repulsive

**Why it looks repulsive:**
- If `REQ.y` (E0 - well depth) is too small or negative, the attractive `(R0/r)^6` term becomes negligible
- If `REQ.x` (R0 - equilibrium distance) is too small, the repulsive `(R0/r)^12` term dominates
- The tip parameters multiply `REQ.yzw` by `0.0091063f`, which could make the attractive term tiny

## Solution: Add GUI View Mode Selection

**To add a dropdown for switching between df-relaxed, Pauli, vdW, and electrostatic views:**

1. **Add enum/view mode variable** in [MolGUI.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:0:0-0:0):
```cpp
enum class AFMViewMode { df_relaxed, Pauli, vdW, electrostatic };
AFMViewMode afm_view_mode = AFMViewMode::df_relaxed;
```

2. **Add DropDownList widget** in [MolGUI::initWiggets()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:519:0-730:1) (around line 718 where panel_AFM is added):
```cpp
DropDownList* panel_AFM_view = ((DropDownList*)gui.addPanel( new DropDownList("AFM View:",5,ylay.x0,5+100, 3) ) ));
panel_AFM_view->addItem("df_relaxed")
              ->addItem("Pauli")
              ->addItem("vdW")
              ->addItem("electrostatic");
panel_AFM_view->setCommand( [&](GUIAbstractPanel* p){ 
    afm_view_mode = (AFMViewMode)((DropDownList*)p)->iSelected;
    makeAFM( afm_iz );  // Re-render with new mode
    return 0; });
```

3. **Modify kernel to output separate components** - The `PPAFM_makeFF` kernel would need to be modified to:
   - Accept a view mode parameter
   - Use GridFF's `Bspline_PLQ` infrastructure to evaluate Pauli, London, Coulomb separately
   - Output only the selected component instead of the combined LJQH

4. **Modify [renderAFM()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:2012:0-2050:1)** to select which data channel to visualize based on `afm_view_mode`

**Key insight:** The GridFF class already has the infrastructure for separate Pauli/London/Coulomb components (`Bspline_PLQ`), but the AFM kernel doesn't use it. You'd need to either:
- Modify the AFM kernel to use GridFF's B-spline interpolation of separate components
- Or add separate output channels to the existing kernel for diagnostic purposes

Would you like me to proceed with implementing the GUI selection and kernel modifications?


---

# USER

I see, so this is not problem, lets use PPAFM_makeFF and getLJQH, that is fine, lets set electrostatics to zero (Q=0), and then we can switch only between df (relaxed), and the potential
it would be also good to have slider choosing the height of that imaging plane.

Also look of the positioning of the imaging grid is correct, it seems to me we have some wrong placement/alignment of the origin, so when the molecule move, and the grid move to cover molecule, we are ot really sampling the field aboove the molecule. 

now figure out updated infor about this. Also write explicitly relevant files and relevant function for AFM imaging

---

# SWE-1.6

Based on my analysis, here are the relevant files and the grid positioning bug I found:

## Relevant Files for AFM Imaging

**Kernel files:**
- [/home/prokop/git/FireCore-master/cpp/common_resources/cl/relax_multi.cl](cci:7://file:///home/prokop/git/FireCore-master/cpp/common_resources/cl/relax_multi.cl:0:0-0:0)
  - [PPAFM_makeFF()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:1004:4-1079:5) (lines 2893-2988) - Evaluates force field on 3D grid using `getLJQH()`
  - [PPAFM_scan()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:1081:4-1147:5) (lines 3033-3120) - Performs df-relaxed scan by relaxing tip positions

**Host-side OpenCL interface:**
- [/home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:0:0-0:0)
  - [PPAFM_makeFF()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:1004:4-1079:5) (line 1005) - Sets up and runs the force field kernel
  - [PPAFM_scan()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:1081:4-1147:5) (line 1082) - Sets up and runs the scan kernel
  - Contains tip parameters: `tipParams{1.661f, 0.0091063f, 0.0f, 0.0f}` (line 152)

**Molecular world interface:**
- [/home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h:0:0-0:0)
  - [evalAFMscan()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h:1346:0-1405:1) (lines 1347-1406) - Main scan function that auto-calculates grid bounds
  - [evalAFM_FF()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h:1407:0-1426:1) (lines 1408-1427) - Evaluates force field on grid

**GUI:**
- [/home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:0:0-0:0)
  - `afm_scan_grid` (line 324) - Grid shape definition
  - [makeAFM()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:2080:0-2091:1) (line 2081) - Main AFM rendering function
  - [renderAFM()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:2012:0-2050:1) (line 2013) - Renders AFM slices
  - [initWiggets()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:519:0-730:1) (line 720) - GUI setup including `panel_AFM`

## Critical Bug: Grid Cell Vector Error

**Location:** [MolWorld_sp3_multi.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h:0:0-0:0) line 1371

```cpp
scan.cell.a=Vec3d{pmax.x-pmin.x,0.0,0.0};
scan.cell.b=Vec3d{pmax.x-pmin.x,0.0,0.0};  // BUG: Should be {0.0, pmax.y-pmin.y, 0.0}
```

**Problem:** The `scan.cell.b` vector is set to `{pmax.x-pmin.x, 0.0, 0.0}` instead of `{0.0, pmax.y-pmin.y, 0.0}`. This means:
- The grid's y-direction cell vector points in x-direction with the same length as x
- The grid is completely collapsed in the actual y-direction
- When the molecule moves, the grid doesn't properly cover the y-extent of the molecule

**Fix:** Change line 1371 to:
```cpp
scan.cell.b=Vec3d{0.0, pmax.y-pmin.y, 0.0};
```

## Grid Origin Calculation

**Location:** [MolWorld_sp3_multi.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common/molecular/MolWorld_sp3_multi.h:0:0-0:0) line 1369

```cpp
scan.pos0=Vec3d{pmin.x,pmin.y,ztop+4.0+5.0+1.0};  // ztop + 10.0 Å above molecule
```

This sets the grid origin at:
- x: `pmin.x` (molecule min x - margin)
- y: `pmin.y` (molecule min y - margin)
- z: `ztop + 10.0` (10 Å above highest atom)

The z-offset seems reasonable for AFM imaging, but the bug in `cell.b` would cause the grid to be misaligned regardless.

## Summary of Required Changes

1. **Fix grid cell vector bug** in `MolWorld_sp3_multi.h:1371`
2. **Set Q=0** in [PPAFM_makeFF](cci:1://file:///home/prokop/git/FireCore-master/cpp/common/OpenCL/OCL_MM.h:1004:4-1079:5) kernel by passing `Qs=float4Zero` from host
3. **Add view mode selection** (df-relaxed vs potential) - this would require:
   - Add enum/variable in [MolGUI.h](cci:7://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:0:0-0:0)
   - Add DropDownList widget in [initWiggets()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:519:0-730:1)
   - Modify kernel to optionally skip relaxation and just sample potential
4. **Add height slider** - add GUIPanel for z-offset in [initWiggets()](cci:1://file:///home/prokop/git/FireCore-master/cpp/common_SDL/SDL2OGL/MolGUI.h:519:0-730:1)

Would you like me to proceed with implementing these fixes?

---

# Additional Fixes Implemented

## Bug: Wrong System Replica Selected in AFM Imaging

**Problem:** When running with multiple replicas (`-m 50`), the AFM scan was imaging a different replica than the one being moved. The kernel was reading atoms from an incorrect system index.

**Root cause:** In `OCL_MM.h:PPAFM_makeFF()`, the `nDOFs.z` field (used by the kernel as the system index `iSys`) was never set before the kernel call. It retained whatever stale value was left by a previous kernel operation.

**Fix:** In `OCL_MM.h:1051-1053`, explicitly set `nDOFs` fields before kernel argument setup:
```cpp
nDOFs.x = nAtoms;  // natoms per system (kernel uses ns.x for loop and i0a)
nDOFs.y = nnode;   // nnode per system  (kernel uses ns.y for i0v offset)
nDOFs.z = isys;    // system replica index -> i0a=isys*natoms, i0v=isys*(natoms+nnode)
```

The kernel computes atom buffer offsets as:
```cpp
const int i0a = iSys * natoms;
const int i0v = iSys * (natoms + nnode);
```

**Location:** `cpp/common/OpenCL/OCL_MM.h` lines 1051-1053

## OpenCL Device Selection

**Problem:** No way to select different OpenCL devices at runtime for testing resource limits on different hardware (e.g., CPU/PoCL vs NVIDIA GPU).

**Fix:** Added CLI flags for device selection in `MolGUIapp_QMMM_multi.cpp`:
- `--cl_list` — list available OpenCL devices with indices
- `--cl_device <index>` — select specific device by index

**Implementation:**
1. Added `OCL_device_picker.h` include
2. Added LambdaDict handlers for `--cl_list` and `--cl_device`
3. Added `ocl_device_idx` member to `MolWorld_sp3_multi` (default -1 = auto-select NVIDIA)
4. Modified `MolWorld_sp3_multi::init()` to use `ocl_device_idx` instead of auto-detecting NVIDIA

**Usage:**
```bash
./MolGUIapp_QMMM_multi --cl_list
./MolGUIapp_QMMM_multi --cl_device 1 -m 1 -x common_resources/xyz/PTCDA ...
```

**Locations:**
- `cpp/apps_OCL/MolecularEditorOCL/MolGUIapp_QMMM_multi.cpp` lines 7, 33-49
- `cpp/common/molecular/MolWorld_sp3_multi.h` line 108, lines 308-310

## OpenCL C Version Detection

**Problem:** Hardcoded `-cl-std=CL2.0` build option failed on CPU/PoCL devices which only support OpenCL C 1.2.

**Fix:** Added runtime detection of device OpenCL C version and use appropriate build option:
1. Query `CL_DEVICE_OPENCL_C_VERSION` during `OCLsystem::init()`
2. Parse version string (e.g., "OpenCL C 1.2 PoCL" → "CL1.2")
3. Use detected version in `buildProgram()` instead of hardcoded CL2.0

**Implementation:**
```cpp
char version_raw[256];
clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(version_raw), version_raw, NULL);
int major=1, minor=2;
sscanf(version_raw, "OpenCL C %i.%i", &major, &minor);
sprintf(cl_c_version, "CL%i.%i", major, minor);
// Later in buildProgram():
sprintf(build_opts, "-I. -cl-std=%s", cl_c_version);
```

**Locations:**
- `cpp/common/OpenCL/OCL.h` line 214 (field), lines 400-407 (detection), lines 506-509 (usage)

## Kernel OpenCL 1.2 Compatibility (PoCL)

**Problem:** `relax_multi.cl` kernel code used OpenCL 2.0 features incompatible with OpenCL 1.2 (PoCL):
1. Address space qualifier errors: passing `__local` pointers to functions expecting `__private` (or vice versa)
2. Macro redefinition warnings

**Fixes:**

1. **Address space qualifiers in `choose_inds_pbc` functions:**
   - Added explicit `__local` qualifier to function parameters
   - Removed const-overloaded versions (OpenCL 1.2 doesn't support overloading by const-ness)
   - Kept single non-const `__local int4*` version which works for both const and non-const arguments

2. **Macro redefinition warnings:**
   - Added `#undef NATOM_LOC` and `#undef NNODE_LOC` before redefining these macros in subsequent kernels

**Location:** `cpp/common_resources/cl/relax_multi.cl` lines 2193-2209, lines 3576-3579, lines 3982-3985

