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