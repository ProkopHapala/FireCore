# USER

OK, now we should test it properly. I saw in my prevous test were bing problems with proper collisions between atoms, atoms were not coliding as they should because they were not mapped properly into clusters or the mapping of local indexes (ghost atoms, neibhbor exclusion) was not do ne probperly

for that prupose we should add several debug potions into the kernels
1) detailed debug prints gated by #define macros and compile-flags
2) add option to selectively turn on/of the collision and bond evaluation selectively. 
3) then we should run synthetic tests e.g. for H2O checking if atoms which should colide really colide and which should not colide does not colide. We should make clearly defined list of interaction which should be seen in deteild debug print writen e.g. as python dictionary or python list of list. And then we should have post-processing function which analyse the debug prints from the kernels (only one iteration) which analyse it. We should think deeply about how to make such regorous and general testing system, that is quete essential and we can reuse that in other projects

the other think we should do is visial debugging where user (me) can pick and pull atoms, this will allow me to test many features quicky and in flaxyble way and gain some intiution what is the problem. I see we already have such functionality in @test_rigid_XPBD_molecules.py and in @XPTB_utils.py but the matplotlib bakcend is slow and problematic. And I wouldprefer to migrate to pyQt and Vispy for fster full-featured 3D interface. Look how we did it in @MolecularPlacerVisPy.py @SequencePlacerVisPy.py 

The mouse picking in 3D may be tricky (especially when using zoom and camera rotation) , so maybe we should test that separately firt. But we should also have alternative pseudo-2D mode where we still use Vispy D rendered but with camera being in top view so we operate only in x-y plane. Note that by any means we should use othographi projection, not perspective, as it is easier to navigate and also the mouse picking works better in that view.

now read the referenced files and make detailed notes what to implement and how to run the test then, Make detailed notes and plan.

---

# Gemini-3-pro

This is an excellent engineering approach. **"If you can't measure it, you can't fix it."**

To solve the complexity of Local Memory caching, Ghost mapping, and Exclusions, we need a **Data-Driven Verification System**.

Here is the comprehensive plan to build a rigorous testing suite for your `make_local_topology` and `solve_collisions_rigid` kernels.

### Phase 1: The "Interaction Truth" (Python Definition)

We define the expected behavior in Python **before** we run the simulation. This serves as the ground truth.

**Scenario: Two Water Molecules (H2O)**
*   **Cluster 0:** Molecule A (O1, H1, H2)
*   **Cluster 1:** Molecule B (O2, H3, H4) - placed close enough to interact.

```python
# verification_data.py

# Format: [Atom_ID_A, Atom_ID_B, Expected_Action]
# Actions: "EXCLUDE_BOND", "EXCLUDE_ANGLE", "COLLIDE", "IGNORE_FAR"

INTERACTION_TRUTH = [
    # --- Molecule A (Intra-molecular) ---
    (0, 1, "EXCLUDE_BOND"),  # O-H bond
    (0, 2, "EXCLUDE_BOND"),  # O-H bond
    (1, 2, "EXCLUDE_ANGLE"), # H-H angle (2nd neighbor)
    
    # --- Molecule B (Intra-molecular) ---
    (3, 4, "EXCLUDE_BOND"),
    (3, 5, "EXCLUDE_BOND"),
    (4, 5, "EXCLUDE_ANGLE"),

    # --- Inter-molecular (A vs B) ---
    (0, 3, "COLLIDE"),       # O-O collision
    (1, 3, "COLLIDE"),       # H-O collision
    (2, 4, "COLLIDE"),       # H-H collision
]

# We will use this dictionary to parse the GPU output logs.
```

### Phase 2: Kernel Instrumentation (OpenCL Macros)

We need a standardized logging format that Python can parse. We will use a macro system that can be completely compiled out for production.

**Header File (`debug_utils.cl`):**
```c
// Only enable prints if this flag is set during compilation
#ifdef ENABLE_DEBUG_PRINTS

// Filter: Only print for specific atoms to save buffer space
#define DEBUG_GID_START 0
#define DEBUG_GID_END   8 

#define LOG_TOPOLOGY(gid, lid, n_ghosts) \
    if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) \
        printf("TOPOLOGY: GID=%d LID=%d n_ghosts=%d\n", gid, lid, n_ghosts);

#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type) \
    if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) \
        printf("MAP: GID=%d L_IDX=%d G_IDX=%d TYPE=%s\n", gid, local_idx, mapped_global_idx, type);

#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action) \
    if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) \
        printf("COLL: MeG=%d OtherL=%d OtherG=%d Dist=%.4f Action=%s\n", \
               my_global, other_local, other_global, dist, action);

#else
// Empty macros for production
#define LOG_TOPOLOGY(gid, lid, n_ghosts)
#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type)
#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action)
#endif
```

### Phase 3: Instrumenting the Kernels

We inject these macros into your key kernels.

#### 1. `make_local_topology` (Checking the Map)
This kernel determines *who* gets loaded into Local Memory.

```c
__kernel void make_local_topology(...) {
    // ... calculate ghost indices ...
    
    if (lid == 0) {
        // Log the Ghosts found for this cluster
        for (int k=0; k<n_found_ghosts; k++) {
            LOG_MAPPING(get_group_id(0)*64, 64+k, ghost_indices[k], "GHOST");
        }
    }
    
    // Log exclusions mapping
    // "I am Atom GID, I am excluding Bond neighbor Global X, which maps to Local Y"
    int4 bond_globals = bonds[gid];
    // ... mapping logic ...
    LOG_MAPPING(gid, local_excl_idx, bond_globals.x, "EXCL_MAP");
}
```

#### 2. `solve_collisions_rigid` (Checking the Logic)
This is where the rubber meets the road.

```c
__kernel void solve_collisions_rigid(...) {
    // ... Load Data ...
    
    // Iterate Local Memory
    for (int j=0; j < total_local_atoms; j++) {
        if (j == lid) continue;

        // Reconstruct Global ID of 'j' for debugging
        int j_global = (j < 64) ? (group_id*64 + j) : ghost_global_ids[j-64];

        // Check Exclusion
        bool is_excluded = (j == my_local_exclusion.x || ...);
        
        float dist = distance(my_pos, l_pos[j]);
        
        if (is_excluded) {
            LOG_COLLISION_CHECK(gid, gid, j, j_global, dist, "SKIP_EXCL");
        } else if (dist < RADIUS_SUM) {
            LOG_COLLISION_CHECK(gid, gid, j, j_global, dist, "COLLIDE");
        } else {
            LOG_COLLISION_CHECK(gid, gid, j, j_global, dist, "TOO_FAR");
        }
    }
}
```

### Phase 4: The Python Test Runner

This script runs the test and validates the output automatically.

```python
import pyopencl as cl
import numpy as np
import re

# 1. SETUP
ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

# Define 2 Water Molecules (Positions in Angstroms)
# Cluster 0: Centered at 0,0,0
pos_data = np.array([
    [0.0, 0.0, 0.0, 1.0],  # O1 (Node)
    [1.0, 0.0, 0.0, 1.0],  # H1 (Cap)
    [0.0, 1.0, 0.0, 1.0],  # H2 (Cap)
    # ... padding to 64 ...
    [3.0, 0.0, 0.0, 1.0],  # O2 (Node) - Cluster 1, close to O1
    [4.0, 0.0, 0.0, 1.0],  # H3
    [3.0, 1.0, 0.0, 1.0],  # H4
    # ... padding ...
], dtype=np.float32)

# Compile with DEBUG flag
options = ["-DENABLE_DEBUG_PRINTS", "-I."]
prg = cl.Program(ctx, kernel_source).build(options=options)

# 2. RUN ONE STEP
# We only need to run the topology and collision kernels once to verify logic
# ... (Buffer creation and enqueue commands) ...
print("Running Kernels...")
queue.finish()

# 3. ANALYZE OUTPUT (Stdout capture logic needed here if not using Jupyter)
# Assuming we captured stdout to a list of strings called 'logs'

logs = [
    "COLL: MeG=0 OtherL=1 OtherG=1 Dist=1.00 Action=SKIP_EXCL",
    "COLL: MeG=0 OtherL=65 OtherG=3 Dist=3.00 Action=COLLIDE",
    # ...
]

print("-" * 30)
print("ANALYSIS REPORT")
print("-" * 30)

regex_coll = re.compile(r"COLL: MeG=(\d+) OtherL=(\d+) OtherG=(\d+) .* Action=(\w+)")

for line in logs:
    match = regex_coll.search(line)
    if match:
        me_g, other_l, other_g, action = match.groups()
        me_g = int(me_g)
        other_g = int(other_g)
        
        # Check against TRUTH
        expected = "UNKNOWN"
        for (a, b, exp) in INTERACTION_TRUTH:
            if (a == me_g and b == other_g) or (a == other_g and b == me_g):
                expected = exp
                break
        
        # Validation Logic
        if expected == "EXCLUDE_BOND" or expected == "EXCLUDE_ANGLE":
            if action != "SKIP_EXCL":
                print(f"FAIL: Atom {me_g} vs {other_g}. Expected {expected}, got {action}")
            else:
                print(f"PASS: {me_g} vs {other_g} correctly excluded.")
                
        elif expected == "COLLIDE":
            if action != "COLLIDE":
                print(f"FAIL: Atom {me_g} vs {other_g}. Expected COLLIDE, got {action}")
            else:
                print(f"PASS: {me_g} vs {other_g} correctly collided.")
```

### Summary of What This Achieves

1.  **Ghost Validation:** If `Action=COLLIDE` appears for an inter-molecular pair (0 vs 3), it proves that `make_local_topology` correctly found Cluster 1 as a neighbor of Cluster 0 and mapped Atom 3 into Cluster 0's Local Memory.
2.  **Exclusion Validation:** If `Action=SKIP_EXCL` appears for pair (0 vs 1), it proves your register-based exclusion logic works and that global indices were correctly converted to local indices.
3.  **Local Memory Integrity:** If `OtherG` (Global ID reconstructed from debug print) matches the expected Global ID, it proves your Local Memory loading loop is aligned correctly.

This setup moves you from "guessing why the physics explodes" to "knowing exactly which index mapping failed."

---

# USER

now we should figure out how to test this thing rigorously.
plese read @RRsp3_XPBD_verification_strategy.chat.md 

and also@Test_system_for_agentic_loops.chat.md 

We have two ways of testing 

we should first ansure momentum conservation (both anguar and linear), for that we need perhaps asymetric molecules, and we should try different kernels, that is the first thing

also it is complicated by the presence of capping atoms, as momentum exchnage and assnebly with capping atoms works differently then with normal atoms, and we often have problem with that. 

What also complicate momentum conservation that we often "hide" (store) some bomentum to rotational DOFs of atoms (angular momentum of qwuaternions) so we must consider total sum, and see if momentum is conserved in long run (no drift)

we also have many kernels, force(/impuls)-based, and position-based we should idensity which works and which not.

The other problem is the non-covalent interaction but we will adress that only after we resolve the purely covalent interaction.

We must create robust test system on two levels
1) automatic check (momentum conservation, topology, index maping etc.)
2) results which can be evaluated by AI (e.g. some set of numbers our output form the simulation/debug prints)
3) visual ouput which human user have to evaluate and judge (plots, interactive GUI)

Also we should consolidate our test script. There is no reason why we have so many test scripts, we need in principle two:
1) one headless script (it can produce .xyz or other trajectroy or other artifacts, e.g. also HTML for user) 
2) interactive GUI where user can pick and rag atoms
all the otehr should be consolidated into this two using CLI options.

We should not have add-hoc single purpose functions direcly in test scripts. Istead we should write shared module like RRsp3_utils.py which agregates all utility function for testing in reusable form

we should also define series of test systems of incresing size and complexity, but some of them must be asyetric. For easier debuging it would be good to have option to fix atom on 2D plane at least at start because it is easier to visualize and control. (we can do 3D calculation but with atoms constrained to 2D plane, perhaps we already have such feature)

Another thing is gated prints from inside kernels, we need to have robust system to do that, how to switch it on/of and differne verosaity level or picked atoms(workgroups) and componets to moniotr 

## Diagnostics / Debugging Tricks (Fail Loudly)

### NaN padding as an out-of-bounds / invalid-access sentinel
- Keep **padding atoms** positions as `NaN` on purpose in device buffers.
- This makes any accidental use of padding indices immediately visible (propagates NaNs into debug overlays / asserts).
- Important semantic rule: **`invm==0` currently implies "invalid/padding" in `RRsp3.upload_state()`**, because it is used to decide which atoms to overwrite as NaN.

### Do NOT use `invm=0` to implement pinning of real atoms
- Pinning/dragging must be done via `fixmask` (position DOF constraints), not by changing `invm`.
- Otherwise the pinned atom becomes NaN in `cl_pos` (by design), poisoning simulation state.

### Reset momentum when constraints change
- The solver uses heavy-ball momentum buffers (`cl_dpos_mom`, `cl_dquat_mom`).
- Whenever we introduce a discontinuous constraint/state change (pin/unpin, drag start/end, or a jump in a dragged target), call `RRsp3.reset_momentum()`.
- Without this, stale momentum can cause sudden jumps/divergence after constraint changes.

### Invariant checks (crash early, print structured stats)
- After each step/download, assert:
  - **Real atoms must be finite**: `np.isfinite(pos[real]).all()`
  - If violated: print indices of first bad atoms + per-cluster min/max/COG summary, then `raise`.
- Keep these checks in the GUI too (not only headless) so interactive debugging immediately stops at the first corruption.

### Build option caching (prevent silent recompilation spam)
- Cache the last-used `build_options` tuple and only recreate `RRsp3` object when options actually change.
- Without this, toggling a checkbox off → every step recompiles the entire OpenCL program, causing massive lag and masking real bugs.

### Bounding box visibility by default
- Bboxes are the primary diagnostic for group/cluster integrity. Keep them **on by default** with high-contrast color (e.g. dark red, 90% opacity, width=2.0).
- If bboxes are hidden, you won't notice when clusters drift apart or collapse.

### Real-time drag signal architecture
- `AtomScene` must emit a **continuous** `sig_atom_moved(idx, pos3)` during drag, not just press/release events.
- Without this, the parent GUI has no way to propagate mouse position into the simulation state.

### VISPY-SEG-ERR as a NaN propagation canary
- `VispyUtils._validated_segs()` prints `[VISPY-SEG-ERR]` when bond/line segments contain non-finite endpoints.
- If this fires after a code change, it is almost certainly because a real atom became NaN (not a Vispy bug).
- Do not suppress it — trace back to the root cause (usually `invm=0` on a real atom, or solver divergence).

### Debug verbosity levels (0..3) with per-cluster stats
- Level 0: silent
- Level 1: basic events (drag start/end, pin toggle)
- Level 2: per-cluster COG/min/max dump on drag events
- Level 3: per-cluster dump on every solver step (noisy, for deep debugging)

now read relevant fioles and make a definite plan 

---

# Kimi 2.6

# RRsp3 Testing Framework — Definite Plan

## Current State Assessment

**Existing test scripts (6 files, significant duplication):**
- [test_RRsp3_momentum.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_momentum.py:0:0-0:0) — Only tests H2O (symmetric), checks `dP` and `dL` from delta buffers, no long-term drift test
- [test_RRsp3_smoke.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_smoke.py:0:0-0:0) — Topology smoke test with local index range checks, ghost downloads
- [test_RRsp3_debug.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_debug.py:0:0-0:0) + [test_RRsp3_debug_runner.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_debug_runner.py:0:0-0:0) — Collision/exclusion verification with kernel printf parsing (good pattern, but ad-hoc)
- [test_RRsp3_convergence.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_convergence.py:0:0-0:0) — Parameter sweep (dt, relaxation, stiffness), complex monolithic script
- [test_RRsp3_vispy.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_vispy.py:0:0-0:0) — Interactive Vispy GUI with pick/drag, molecule loading, ~800 lines

**Key gaps identified:**
- Momentum test uses only H2O (symmetric → angular momentum tests are trivial/ambiguous)
- No systematic testing across the 5 port kernel variants (`current`, `orig`, `substep`, `shapematch`, `eigen`)
- Capping atoms (H) are not nodes → their momentum is handled via `bkSlots` recoil gather; no explicit validation that cap-atom momentum is conserved
- Debug printf system exists but is coarse (global GID range only, no per-component or per-workgroup targeting)
- No 2D-constrained mode for simplified debugging
- No structured output format for AI evaluation

## Phase 1: Shared Utility Module (`RRsp3_utils.py`)

**Location:** `pyBall/RigidAtomFF/RRsp3/RRsp3_utils.py`

**Functions to consolidate (extracted from existing test scripts):**

| Function | Source | Purpose |
|----------|--------|---------|
| [make_ports_from_neighs()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_vispy.py:73:0-86:28) | All test scripts | Build port_local + Kflat from positions and neighs |
| [reorder_nodes_first()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/XPDB_new/test_rigid_XPBD_molecules.py:33:0-55:31) | test_RRsp3_vispy.py | Reorder atoms: nodes first, caps second |
| [load_molecule_any_xyz()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_vispy.py:60:0-70:35) | test_RRsp3_vispy.py | Load XYZ, find bonds via AtomicSystem, reorder |
| `make_test_system()` | New | Factory for test molecules (H2O, CH3OH, CH2NH, etc.) |
| `compute_momentum_conservation()` | test_RRsp3_momentum.py | Compute `dP` and `dL` from downloaded deltas |
| `compute_total_momentum()` | New | Full system momentum including rotational DOFs |
| `parse_kernel_debug_logs()` | test_RRsp3_debug.py | Parse COLL/TOPOLOGY/MAP printf lines |
| `validate_interactions()` | test_RRsp3_debug.py | Compare parsed logs against INTERACTION_TRUTH |
| [write_xyz_frame()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_momentum.py:30:0-35:58) | test_RRsp3_momentum.py | Trajectory output |
| [check_local_ranges()](cci:1://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_smoke.py:34:0-53:147) | test_RRsp3_smoke.py | Validate local index bounds |

**New functions needed:**
- `make_asymmetric_molecule(name)` — CH3OH (methanol), CH2NH (formaldimine), CH3CH2OH fragments
- `compute_total_angular_momentum(pos, vel, quat, omega, m)` — Include both orbital and spin angular momentum
- `perturb_state_planar(pos, quat, rng, plane='xy')` — 2D perturbation for planar debugging
- `build_interaction_truth(molecules, shift)` — Auto-generate expected COLLIDE/EXCLUDE pairs

## Phase 2: Kernel Debug Print Enhancement

**Current state:** [RRsp3.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3.cl:0:0-0:0) lines 5-20 have basic `LOG_COLLISION_CHECK`, `LOG_TOPOLOGY`, `LOG_MAPPING` gated by `ENABLE_DEBUG_PRINTS` and hardcoded `DEBUG_GID_START/END`.

**Target design:** Verbosity-level + component-bitmask + workgroup-targeting system

```c
// In RRsp3.cl, replace lines 5-20
#ifdef ENABLE_DEBUG_PRINTS

#ifndef DEBUG_VERBOSITY
#define DEBUG_VERBOSITY 1   // 0=none, 1=errors only, 2=summary, 3=per-atom
#endif

#ifndef DEBUG_COMPONENTS
#define DEBUG_COMPONENTS 0xFFFF  // bitmask: 1=collision, 2=port, 4=topology, 8=correction
#endif

#ifndef DEBUG_TARGET_WG
#define DEBUG_TARGET_WG -1   // -1 = all workgroups
#endif

#define DBG_COLL  1
#define DBG_PORT  2
#define DBG_TOPO  4
#define DBG_CORR  8

#define DBG_ENABLED(comp) (DEBUG_VERBOSITY >= 3 && (DEBUG_COMPONENTS & comp))
#define DBG_WG_MATCH (DEBUG_TARGET_WG < 0 || get_group_id(0) == DEBUG_TARGET_WG)
#define DBG_ATOM_MATCH(gid) (DEBUG_TARGET_ATOM_START <= gid && gid < DEBUG_TARGET_ATOM_END)
```

**Compile-time options (passed via `build_options` in RRsp3.py):**
```python
# Example usage in test script
build_opts = [
    "-DENABLE_DEBUG_PRINTS",
    "-DDEBUG_VERBOSITY=3",
    "-DDEBUG_COMPONENTS=7",      # collision + port + topology
    "-DDEBUG_TARGET_WG=0",         # only workgroup 0
    "-DDEBUG_TARGET_ATOM_START=0",
    "-DDEBUG_TARGET_ATOM_END=8",
    "-DENABLE_COLL=1",
    "-DENABLE_PORT=0",           # selectively disable for isolation testing
]
```

**Why this matters:** When debugging capping atom momentum issues, we need to see only the port/gather logic without collision noise. When debugging ghost mapping, we need topology only.

## Phase 3: Test System Hierarchy

Define **5 test systems** of increasing complexity, all with explicit node/cap split:

| Level | System | Atoms | Nodes | Caps | Asymmetric? | Purpose |
|-------|--------|-------|-------|------|-------------|---------|
| L0 | Single H2O | 3 | 1 (O) | 2 (H) | No | Smoke test, basic momentum |
| L1 | Two H2O (close) | 6 | 2 | 4 | No | Collision + exclusion, ghost mapping |
| L2 | Single CH3OH (methanol) | 6 | 2 (C,O) | 4 (H×4) | **Yes** | Asymmetric inertia, proper angular momentum test |
| L3 | Two CH3OH (close) | 12 | 4 | 8 | **Yes** | Full covalent + non-covalent with asymmetric bodies |
| L4 | CH2NH + H2O mixed | 7 | 3 | 4 | **Yes** | Mixed topology, different node types |

**Why asymmetric molecules are critical:**
- H2O has O on the bisector of H-H; pure rotation around any axis through O leaves the H-H pair invariant → angular momentum conservation tests are insensitive
- CH3OH has C-O bond with methyl group on one side and OH on the other → rotation produces clearly distinguishable cap motions → any error in torque propagation is immediately visible in `dL`

**2D planar constraint:**
- Add `fixmask` bit interpretation: bits 1,2,4 fix x,y,z
- For 2D mode: set `fixmask |= 4` for all atoms (constrain to z=0)
- The 3D kernel still runs, but atoms cannot move in z → effectively 2D
- Easier to visualize and debug because you see the full configuration from top view

## Phase 4: Consolidated Headless Test Script

**File:** `test_RRsp3_headless.py` (replaces: momentum, smoke, debug, convergence, debug_runner)

**CLI routing via subcommands:**

```bash
# Momentum conservation test
python test_RRsp3_headless.py momentum --system L2 --kernel all --steps 100 --dt 0.1

# Topology verification (collision/exclusion mapping)
python test_RRsp3_headless.py topology --system L1 --shift 1.5

# Smoke test (basic step completion, NaN checks)
python test_RRsp3_headless.py smoke --system L0

# Convergence sweep
python test_RRsp3_headless.py converge --system L2 --dts 0.05,0.1,0.2 --relaxations 0.3,0.5,0.7

# Full suite
python test_RRsp3_headless.py suite --level L3 --kernels current,substep,shapematch,eigen
```

**Output format (AI-evaluable):**
```json
{
  "test_type": "momentum",
  "system": "CH3OH_x2",
  "kernel": "shapematch",
  "steps": 100,
  "results": {
    "dP_max": [1.2e-6, 3.4e-6, 2.1e-6],
    "dL_max": [5.6e-5, 1.2e-4, 8.9e-5],
    "dP_drift_final": 2.1e-6,
    "dL_drift_final": 9.8e-5,
    "pass": true
  },
  "artifacts": {
    "trajectory": "/tmp/rrsp3_CH3OH_x2_shapematch.xyz",
    "plot": "/tmp/rrsp3_momentum_drift.png"
  }
}
```

## Phase 5: Consolidated Interactive GUI Script

**File:** `test_RRsp3_gui.py` (replaces: vispy test)

**CLI options:**
```bash
python test_RRsp3_gui.py --system L2 --kernel current --dt 0.1 --relax 0.5
python test_RRsp3_gui.py --xyz /path/to/mol.xyz --kernel eigen --2d
```

**Features (from existing vispy test + improvements):**
- Orthographic projection (already partially there, enforce strictly)
- 2D top-down mode: constrain camera to z-axis view, disable rotation
- Pick and drag atoms with mouse
- Run/stop relaxation toggle
- Step-by-step advance (single kernel step)
- Show/hide: bboxes, ghost atoms, port vectors
- Real-time momentum monitor overlay (dP, dL plotted as bars)
- Kernel switcher dropdown (current/orig/substep/shapematch/eigen)

**Why consolidate:** The existing [test_RRsp3_vispy.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_vispy.py:0:0-0:0) is ~800 lines with heavy duplication of molecule setup logic. By importing from `RRsp3_utils.py`, the GUI script becomes ~300 lines of visualization-specific code only.

## Phase 6: Momentum Conservation — Detailed Implementation

**Current test logic (test_RRsp3_momentum.py):**
- Downloads `dpos_coll`, `dpos_node`, `drot_node`, `dpos_neigh`
- Gathers `dx_port` via `bkSlots`
- Computes `dP = sum(m * dx)` and `dL = sum(r × m*dx) + sum(I * dtheta)`
- Checks `|dP| < epsP` and `|dL| < epsL` per-step

**Problems with current approach:**
1. Only tests delta buffers, not full state evolution (no drift detection)
2. H2O is too symmetric
3. Does not test different kernels
4. Capping atoms: the test manually gathers `dx_port` via `bkSlots` on CPU, but this duplicates kernel logic → if there's a bug in the kernel gather, the test won't catch it because it reproduces the same gather logic

**Improved approach:**

```python
def test_momentum_conservation(system, kernel, steps, dt, relaxation):
    sim = build_system(system)
    sim.reset_momentum()
    
    # Reference state for drift tracking
    P_hist = []
    L_hist = []
    
    for step in range(steps):
        sim.step_cluster(..., port_kernel=kernel)
        
        # Compute total momentum from FULL STATE (not deltas)
        pos, quat = sim.download_pos_quat()
        P, L = compute_total_momentum(pos, quat, sim.masses)
        P_hist.append(P)
        L_hist.append(L)
    
    # Check 1: Per-step invariance (like current test, but more rigorous)
    dP = np.diff(P_hist, axis=0)
    dL = np.diff(L_hist, axis=0)
    
    # Check 2: Long-term drift (cumulative error)
    P_drift = np.max(P_hist, axis=0) - np.min(P_hist, axis=0)
    L_drift = np.max(L_hist, axis=0) - np.min(L_hist, axis=0)
    
    # Check 3: Cap-atom specific (verify caps move correctly when nodes rotate)
    # For asymmetric molecules, cap displacements should reflect node rotation
    # not just linear translation
```

**Cap-atom momentum validation:**
For CH3OH with 2 nodes (C, O) and 4 caps (H×3 on C, H on O):
- If C rotates, the 3 H caps on C should recoil in a way that conserves total L
- The test should verify that `sum(m_cap * v_cap)` from the gathered recoils balances `I_node * omega_node`

## Phase 7: Topology/Index Mapping Verification

**Reuses the existing [test_RRsp3_debug.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/test_RRsp3_debug.py:0:0-0:0) pattern but generalized:**

```python
# In RRsp3_utils.py
def build_interaction_truth(molecules, packed, shift):
    """Auto-generate expected interaction table from molecular topology."""
    truth = {}
    natoms_per_mol = [len(m['elems']) for m in molecules]
    
    # Intra-molecular exclusions
    for mol_idx, mol in enumerate(molecules):
        base = sum(natoms_per_mol[:mol_idx])
        for (i, j) in mol['bonds']:
            truth[(base+i, base+j)] = 'EXCLUDE_BOND'
        # 2nd neighbors (angles)
        # ... build from bond graph ...
    
    # Inter-molecular collisions
    # For each pair of molecules, compute if any atom pairs are within collision distance
    # ... using shift parameter ...
    
    return truth
```

**Validation pipeline:**
1. Build truth table from Python topology
2. Compile kernel with `ENABLE_DEBUG_PRINTS` + `DEBUG_COMPONENTS=1` (collision only)
3. Run one step with `ENABLE_PORT=0` (isolate collision logic)
4. Parse stdout for `COLL:` lines
5. Assert every truth-table pair appears with correct action

## Phase 8: Non-Covalent Interactions (Deferred)

Address **only after** covalent port dynamics is fully validated. The plan:
1. Validate port-only (collision disabled): `ENABLE_COLL=0, ENABLE_PORT=1`
2. Validate collision-only (ports disabled): `ENABLE_COLL=1, ENABLE_PORT=0`
3. Validate combined: both enabled

This isolation prevents debugging two complex systems simultaneously.

## Implementation Order

| Step | Task | Files touched | Estimated lines |
|------|------|-------------|-----------------|
| 1 | Create `RRsp3_utils.py` with extracted utilities + new molecule builders | New file | ~400 |
| 2 | Enhance kernel debug macros in [RRsp3.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3.cl:0:0-0:0) | [RRsp3.cl](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/RRsp3/RRsp3.cl:0:0-0:0) lines 5-20 | ~40 |
| 3 | Add `make_ch3oh_geometry()` and `make_ch2nh_geometry()` to [XPTB_utils.py](cci:7://file:///home/prokophapala/git/FireCore/pyBall/RigidAtomFF/shared/XPTB_utils.py:0:0-0:0) or `RRsp3_utils.py` | `RRsp3_utils.py` | ~60 |
| 4 | Build `test_RRsp3_headless.py` with `momentum`, `topology`, `smoke`, `converge` subcommands | New file | ~500 |
| 5 | Port existing vispy test to `test_RRsp3_gui.py` using `RRsp3_utils.py` | New file | ~350 |
| 6 | Write momentum test for all 5 kernels × 5 test systems | `test_RRsp3_headless.py` | ~200 |
| 7 | Run full suite, fix failures, document which kernel passes which test | Chat/audit.md | — |
| 8 | Clean up deprecated test scripts | Delete old files | — |

## Kernel Pass/Fail Expectations (Hypothesis)

Based on the code analysis:

| Kernel | Momentum Conserved | Notes |
|--------|-------------------|-------|
| `current` (compute_ports_cluster_rigid) | **Should pass** | Full XPBD with mass-weighted denominator |
| `orig` (compute_ports_cluster_rigid_orig) | **Should pass** | Earlier version, similar structure |
| `substep_optimized` | **May fail** | Massless rotation solver, angular momentum not conserved by design |
| `shapematch` | **May fail** | Massless Procrustes solve, no inertia in rotation |
| `eigen` | **May fail** | Massless quaternion eigensolve |

**Critical insight:** The "massless" solvers (substep, shapematch, eigen) solve for optimal geometry without considering inertia. They will violate angular momentum conservation by design. The test framework must:
1. Confirm they fail the strict conservation test
2. Measure *how much* they violate it (quantify the error)
3. Potentially accept them if the error is bounded and the gain in convergence speed is worth it

## Immediate Next Actions

1. **Switch to Code mode** so I can begin implementation
2. I will start with `RRsp3_utils.py` (the shared utility module) — this unblocks everything else
3. Then enhance kernel debug macros
4. Then build the headless test script with momentum and topology subcommands
5. Run against all 5 kernels to establish baseline pass/fail matrix

**Do you want me to proceed with implementation? Switch to Code mode and I will start with `RRsp3_utils.py`.**