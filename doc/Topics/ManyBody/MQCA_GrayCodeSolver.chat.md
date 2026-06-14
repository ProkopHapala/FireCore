# MQCA Ground-State Solver: Gray Code Brute-Force Search

---

# USER

What we want to do is to solve molecualr cellular automata MQCA
This is basically lattice of sites with some onsite energy and neighbor capacitve coupling described by Ising-like hamiltonina 

$$
H = \sum_i n_i[ \epsilon_i + \sum_{i<j} W_{ij} n_j ]
$$

n_i can be only 0 or 1 (binary number)

Now I already have some solver of such system by mote carlo for large systems, take a look here

look in @hubbard.cl @HubbardSolver.py @HubbardSolver.md @HubbardSolver_new.md 

but Now I want to do something even simpler also using pyOpenCL and @OpenCLBase.py 

The new solver should be for smaller system, max 16 active sites so each site can be evaluated by its own thread within workgroup (workgroupsize=16), every site can be either on or off => 2^16 combination 65536, that is doable by bruteforece sequnetially which significantly siplify the algortihm. Our goal is basically just to find the groundstate. If we do it sequentially we do not need reduction. Perhaps is it possible to desing such path over all binary combinations that we always switch only single site, so we need to evaluate only the change of the energy (i.e. update only single site), not the total energy (summing over all sites). The threads within workgroup will compute coupling of that one site with the others.

We assume that matrix Wij is sparse and there is max 8 non-zero neigbrs for each site with which it is coupled by Wij. Althoguh in pracitce this is on grid, we want to make it general so we store the coupling values and indexes if neighbors, it is up to python driver to set this neibhor indexes and coupling Wij to represent given geometry of system (typically cluster) 

The cluster can have also some non-active (fixed) neighbors points which are influencing active sites by coulomb coupling Wij, but since they are fixed we can precacluate their effect \sum_{i<j} W_{ij} n_j and simply adjust \epsilon_i - which is equvalent to constant bias postential. 

Now please read the provided files and make comprehnsive desing document

---

To make sure that we iterate over all binary combination 

What you want is a Gray code sequence.

In Gray code, each successive 16-bit value differs from the previous one by exactly one bit, so iterating through all 65536 combinations only flips a single bit each step.

For integer i, the Gray-coded value is:

G(i)=i\oplus\left(i\gg 1\right)

Example in C/C++:

```cpp
for(uint32_t i=0; i<(1u<<16); i++){
    uint32_t g = i ^ (i >> 1);
    // use g as the current 16-bit pattern
}
```

The sequence starts like:

```text
0000
0001
0011
0010
0110
0111
0101
0100
...
```

Each step changes only one bit.

If you also want to know which bit changed between successive Gray-code values, you can do:

```cpp
uint32_t prev = 0;
for(uint32_t i=1; i<(1u<<16); i++){
    uint32_t g    = i ^ (i >> 1);
    uint32_t diff = g ^ prev;   // exactly one bit set
    int bitIndex  = __builtin_ctz(diff);
    prev = g;
}
```

Here bitIndex tells you which bit flipped.


---

OK, now lets think about how to use it. We want to use massive paralelism of GPU, therefore we need to simulated many (like thousancs) of the autmata each described by its own copy of the system and processed by different workgroup.

The application for which we want to use it is to find parameters resp geometries of the molecular clustes (which are encoded in the values of Wij and e_i) for which the cluster perform some binary operations. Therefore we want to run some genertica algorithm where we vary the shapes of clusters and theri parameters, and always compute binary locical tablle. e.g. table for XOR, OR, AND, NAND, always 2 inpts and 1 output, each such table have 4 possible combinations. These 4 combinations can be computed in parallel, but maybe it is easier if we compute them sequentially.

What is  most importaint is to realize what are inputs ad outputs for our MQCA?

Inputs are fixed values on some sites. For example we can modify on-site energy of certain site by W*1 or W*0 indicating that neighboring site is on, or that neighboring site is of. This meas we consider that there is some boundary fixed site to which our active sites are coupled. 

The ouput is simply that we compute the ground-state for our hamiltonina (energy minimum - configuration with minimum energy) and read value of one selected site which we designated as output. Choice which sites are input and which are output is up to user. 

For testing we must make nice visialization where we plot input sites (those are fixed and not included in simulation, they enter simulations only via precomputed biases on the neighboring active sites), output sites, and internal=active sites (output site is one selected internal active sites)
We plot the ground-state for all combinations of the two inputs (e.g. as black/write, red/blue values of the active sites), and print final logical table, if it agree with given logical function ( XOR, OR, AND, NAND) - probably it is better to check all of them, to increase the change we get useful device.


Now we can split the task into two sub-tasks of increasing complexity

1) simple is to scan over paramters for fixed shape of the cluster. We fix value of epsilon0  (same for all active sites, unless bias for fixed sites is applied), without loss of generality lets say  e=1 . Then assuming square lattice we have two types of neighbor - alaong cartesian direction W1 and diagonally W2. We can vary W1 and W2 independnetly as two axis - this will basically plot an image showing for which combination of W1 W2 given device produce meaningfull logical table. We can get like phase diagram where same device would produce AND for some values and OR ofr other values, perhaps, maybe ...

2)  The more complicated task is to vary the geometry of the cluster.  Perhaps we need some claver stochastic algorithm which will start with some cluster as map of active sites on a grid, and will mover sites on the surface, but will always make sure the whole cluster is connected - active sites must be all adjecent to at least one other site. The input and output sites can be fixed for simplicity. Now this generated shapes should be replicated and simulated for different values of W1 W2 in parallel, se basically we created multiple low-res 2D imaged  as in (1), stack them all together and sent to GPU to evaluate them in parallel (to find ground-state) and then we check which desing implement which logical operation, if any

---

please think it through and write detailed notes more systematically into @MQCA_GrayCodeSolver.md 




# Plan (Kimi K2.5)


---

# MQCA Ground-State Solver: Gray Code Brute-Force Search

## Overview

This document describes a specialized OpenCL solver for Molecular Quantum Cellular Automata (MQCA) systems using **Gray code traversal** for efficient brute-force ground-state search. This solver is designed for **small systems** (max 16 active sites) where exhaustive enumeration of all 2^16 = 65536 configurations is tractable.

### Key Design Goals

1. **Max 16 active sites** - small enough for exhaustive enumeration
2. **Workgroup size = 16** - one thread per site within a workgroup
3. **Gray code iteration** - only 1 bit flips between consecutive states, enabling incremental energy updates
4. **Sparse coupling matrix** - max 8 neighbors per site
5. **Support for fixed (non-active) neighbors** - precalculated bias contribution to on-site energies


## Physical Model

### Hamiltonian

The system is described by an Ising-like Hamiltonian:

```
H = Σ_i n_i [ ε_i + Σ_{j<i} W_{ij} n_j ]
```

Where:
- `n_i ∈ {0, 1}` - occupancy of site i (binary)
- `ε_i` - on-site energy (including bias from fixed neighbors)
- `W_{ij}` - capacitive coupling between sites i and j (sparse, max 8 non-zero per site)

### Fixed Neighbors (Bias)

For fixed (non-active) neighbors k, their effect is precalculated:

```
ε_i^eff = ε_i + Σ_{k∈fixed} W_{ik} n_k
```

This is done on the Python side before uploading to GPU.

## Algorithm: Gray Code Traversal

### Why Gray Code?

Gray code ensures that consecutive 16-bit values differ by exactly one bit. This enables **incremental energy evaluation**:

- Instead of summing over all sites for each configuration
- We only calculate the energy **change** from flipping one site
- Thread i computes the coupling of site i to the current configuration

### Gray Code Formula

```c
// G(i) = i ^ (i >> 1)
uint16_t gray_code(uint16_t i) {
    return i ^ (i >> 1);
}
```

### Detecting Which Bit Flipped

```c
uint16_t prev_g = 0;
for (uint16_t i = 1; i < (1u << nSite); i++) {
    uint16_t g    = i ^ (i >> 1);     // Current Gray code
    uint16_t diff = g ^ prev_g;        // Exactly one bit set
    int bit_idx   = __builtin_ctz(diff);  // Index of flipped bit (0-15)
    prev_g = g;
    // Now evaluate energy change from flipping bit_idx
}
```

## OpenCL Kernel Design

### Kernel: `mqca_graycode_groundstate`

**Purpose**: Find the ground state energy and configuration for a single system instance.

**Parallelization Strategy**:
- Each **work-group** handles one independent system instance
- Within a work-group, **16 threads** collaborate on energy evaluation
- Thread `lid` computes the contribution from site `lid` to the total energy change

**Memory Layout**:

| Buffer | Type | Shape | Description |
|--------|------|-------|-------------|
| `Esite` | `__global float*` | `[nInstances, nSite]` | On-site energies ε_i (including fixed neighbor bias) |
| `W_val` | `__global float*` | `[nSite, maxNeigh]` | Sparse coupling values W_{ij} |
| `W_idx` | `__global int*` | `[nSite, maxNeigh]` | Neighbor indices j for each W_{ij} |
| `nNeigh` | `__global int*` | `[nSite]` | Number of neighbors per site |
| `E_out` | `__global float*` | `[nInstances]` | Output: ground state energy |
| `occ_out` | `__global ushort*` | `[nInstances]` | Output: ground state occupancy (16-bit mask) |

### Kernel Signature

```c
#define MAX_SITES     16
#define MAX_NEIGHBORS 8
#define N_COMBINATIONS (1 << MAX_SITES)  // 65536 for 16 sites

__kernel void mqca_graycode_groundstate(
    const int nSite,                    // 1: Number of active sites (<= 16)
    __global const float* Esite,        // 2: [nInstances, nSite] On-site energies
    __global const float* W_val,        // 3: [nSite, maxNeigh] Coupling values
    __global const int*   W_idx,        // 4: [nSite, maxNeigh] Neighbor indices
    __global const int*   nNeigh,       // 5: [nSite] Number of neighbors per site
    const int maxNeigh,                 // 6: Max neighbors (<= 8)
    __global float*       E_out,        // 7: [nInstances] Output: ground state energy
    __global ushort*      occ_out       // 8: [nInstances] Output: ground state occupancy
);
```

### Algorithm Flow

```c
// Per work-group (one system instance)
const int iinst  = get_group_id(0);     // Instance index
const int lid    = get_local_id(0);     // Thread 0-15, handles site lid
const int offset = iinst * nSite;

// Local memory for current configuration
__local ushort occ_mask;                // Current 16-bit occupancy mask
__local float dE_contrib[16];             // Each thread's contribution to dE

// Initialize
if (lid == 0) {
    occ_mask = 0;                       // Start with all sites unoccupied
    E_min    = 0.0f;                    // Initial energy (all n_i = 0)
    occ_min  = 0;                       // Initial ground state
}
barrier(CLK_LOCAL_MEM_FENCE);

// Compute initial total energy E_current
float E_current = 0.0f;
// ... parallel reduction over sites ...

// Gray code traversal
ushort prev_g = 0;
for (ushort i = 1; i < (1u << nSite); i++) {
    ushort g    = i ^ (i >> 1);         // Current Gray code
    ushort diff = g ^ prev_g;           // Bit that changed
    int flip_idx = __builtin_ctz(diff); // Index 0-15
    
    // All threads compute their site's contribution to dE
    if (lid == flip_idx) {
        // Energy change from flipping site flip_idx
        // dE = (1 - 2*n_i) * [ε_{flip_idx} + Σ_{j∈neigh} W_{flip_idx,j} * n_j]
        float site_E = Esite[offset + flip_idx];
        
        // Add coupling to occupied neighbors
        int iw0 = flip_idx * maxNeigh;
        for (int k = 0; k < nNeigh[flip_idx]; k++) {
            int j = W_idx[iw0 + k];
            if ((occ_mask >> j) & 1u) {
                site_E += W_val[iw0 + k];
            }
        }
        
        // If currently occupied, flipping removes energy; if unoccupied, adds energy
        int n_flip = (occ_mask >> flip_idx) & 1u;
        dE_contrib[lid] = (n_flip ? -1.0f : 1.0f) * site_E;
    } else {
        dE_contrib[lid] = 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Thread 0 aggregates dE and updates state
    if (lid == 0) {
        float dE = dE_contrib[0] + dE_contrib[1] + ... + dE_contrib[15]; // Reduction
        E_current += dE;
        occ_mask ^= (1u << flip_idx);   // Flip the bit
        
        if (E_current < E_min) {
            E_min   = E_current;
            occ_min = occ_mask;
        }
    }
    prev_g = g;
    barrier(CLK_LOCAL_MEM_FENCE);
}

// Write output
if (lid == 0) {
    E_out[iinst]   = E_min;
    occ_out[iinst] = occ_min;
}
```

## Python Interface

### Class: `MQCAGraySolver`

Extends `OpenCLBase` to provide high-level interface.

```python
class MQCAGraySolver(OpenCLBase):
    """
    MQCA Ground State Solver using Gray code brute-force search.
    
    Suitable for small systems (max 16 sites) where exhaustive
    enumeration of all 2^16 = 65536 states is tractable.
    """
    
    MAX_SITES = 16
    MAX_NEIGHBORS = 8
    WORKGROUP_SIZE = 16
    
    def __init__(self, device_index=0):
        super().__init__(nloc=self.WORKGROUP_SIZE, device_index=device_index)
        # Load kernel from cl/mqca_gray.cl
        
    def solve(self, Esite, W_sparse, nSite=None):
        """
        Find ground state for one or more system instances.
        
        Parameters
        ----------
        Esite : np.ndarray, shape (nInstances, nSite) or (nSite,)
            On-site energies (including bias from fixed neighbors)
        W_sparse : tuple (W_val, W_idx, nNeigh)
            - W_val: (nSite, maxNeigh) coupling values
            - W_idx: (nSite, maxNeigh) neighbor indices
            - nNeigh: (nSite,) number of neighbors per site
            
        Returns
        -------
        E_min : np.ndarray, shape (nInstances,)
            Ground state energies
        occ_min : np.ndarray, shape (nInstances,)
            Ground state occupancies (16-bit unsigned integers)
        """
        
    def solve_batch(self, Esite_list, W_sparse_list):
        """
        Solve multiple independent systems in parallel.
        Each system can have different site energies and couplings.
        """
```

### Buffer Management

```python
def _allocate_buffers(self, nInstances, nSite, maxNeigh):
    """Allocate GPU buffers for the solver."""
    sz_f = np.dtype(np.float32).itemsize
    sz_i = np.dtype(np.int32).itemsize
    sz_u = np.dtype(np.uint16).itemsize
    
    buffs = {
        "Esite":    sz_f * nInstances * nSite,
        "W_val":    sz_f * nSite * maxNeigh,
        "W_idx":    sz_i * nSite * maxNeigh,
        "nNeigh":   sz_i * nSite,
        "E_out":    sz_f * nInstances,
        "occ_out":  sz_u * nInstances,
    }
    self.try_make_buffers(buffs)
```

---

## File Structure

```
pyBall/OCL/
├── cl/
│   ├── hubbard.cl              # Existing kernels
│   └── mqca_gray.cl            # NEW: Gray code MQCA solver
├── OpenCLBase.py               # Base class (existing)
├── HubbardSolver.py            # Existing solver
└── MQCAGraySolver.py           # NEW: Gray code solver
```

---

## Comparison with Existing Solvers

| Feature | `solve_minBrute_fly` | `solve_local_updates` | `mqca_graycode_groundstate` |
|---------|---------------------|----------------------|----------------------------|
| Max sites | 32 | 256 | 16 |
| Workgroup size | 32 | 16 | 16 |
| Algorithm | Brute force | Monte Carlo | Gray code traversal |
| Energy update | Full recalculation | Local delta | Incremental (1 site) |
| States evaluated | 2^nSites | Sampled subset | All 2^nSites |
| Parallelism | Per tip position | Per tip position | Per system instance |
| Reduction | Workgroup reduction | Workgroup reduction | Within workgroup |
| Guaranteed ground state | Yes | No | Yes |

### When to Use Each Solver

- **`solve_minBrute_fly`**: Medium systems (16-32 sites), many tip positions, full enumeration needed
- **`solve_local_updates`**: Large systems (>32 sites), approximate solution acceptable
- **`mqca_graycode_groundstate`**: Small systems (<=16 sites), exact ground state required, independent systems

---

## Kernel Implementation Details

### Key OpenCL Features Used

1. **`__builtin_ctz()`**: Count trailing zeros to find flipped bit index
2. **`barrier(CLK_LOCAL_MEM_FENCE)`**: Synchronize threads within workgroup
3. **`__local` memory**: Store current occupancy mask and partial dE contributions
4. **16-bit ushort**: Compact representation of occupancy (fits in GPU registers)

### Performance Considerations

1. **Memory access pattern**: 
   - `Esite` accessed coalesced (instance-major layout)
   - `W_val`/`W_idx` same for all instances (can use `__constant` or `__local` cache)
   
2. **Register pressure**:
   - Minimize live variables in the 65536-iteration loop
   - Keep `E_current`, `occ_mask`, `prev_g` in registers
   
3. **Occupancy**:
   - 16 threads/workgroup ensures high occupancy on modern GPUs
   - Multiple workgroups per compute unit for latency hiding

### Numerical Precision

- All energies stored as `float32` (sufficient for eV-scale energies)
- Accumulation order is deterministic (Gray code sequence is fixed)
- No floating-point reduction across workgroups (each workgroup independent)

## Usage Example

```python
import numpy as np
from pyBall.OCL.MQCAGraySolver import MQCAGraySolver

# Setup: 8 sites in a ring
def make_ring_couplings(n, R, W0):
    """Create nearest-neighbor couplings for a ring."""
    W_val = np.zeros((n, 8), dtype=np.float32)
    W_idx = np.zeros((n, 8), dtype=np.int32)
    nNeigh = np.full(n, 2, dtype=np.int32)
    
    for i in range(n):
        j1 = (i + 1) % n
        j2 = (i - 1) % n
        W_idx[i, 0] = j1
        W_idx[i, 1] = j2
        W_val[i, 0] = W0
        W_val[i, 1] = W0
    
    return (W_val, W_idx, nNeigh)

# Create solver
solver = MQCAGraySolver()

# System parameters
nSite = 8
Esite = np.random.uniform(-1.0, 1.0, nSite).astype(np.float32)
W_sparse = make_ring_couplings(nSite, R=5.0, W0=0.5)

# Solve
E_min, occ_min = solver.solve(Esite, W_sparse)

print(f"Ground state energy: {E_min:.6f} eV")
print(f"Ground state occupancy: {occ_min:016b}")
print(f"Occupied sites: {[i for i in range(nSite) if (occ_min >> i) & 1]}")

# Batch solve multiple instances
nInstances = 100
Esite_batch = np.random.uniform(-1.0, 1.0, (nInstances, nSite)).astype(np.float32)
E_mins, occ_mins = solver.solve(Esite_batch, W_sparse)
```

## Future Extensions

### Possible Enhancements

1. **Multi-site Gray code**: For systems with >16 sites, use 32-bit or 64-bit Gray codes
2. **Hybrid approach**: Use Gray code for small sub-clusters, Monte Carlo for global
3. **Thermal averaging**: Extend to compute Boltzmann-weighted averages over all states
4. **Energy spectrum**: Track all energy levels (not just ground state) for small systems
5. **Bitwise optimization**: Use OpenCL vector types (`ushort16`) for parallel neighbor checks

### Thermal Averaging Extension

```c
// Track partition function Z and average energy
__local float Z_local;
__local float E_avg_local;

if (lid == 0) {
    Z_local = 0.0f;
    E_avg_local = 0.0f;
}

// In the loop:
float w = exp(-beta * E_current);
if (lid == 0) {
    Z_local += w;
    E_avg_local += w * E_current;
}
```

## Summary

The Gray code MQCA solver provides:
- **Guaranteed exact ground state** for systems up to 16 sites
- **Efficient incremental evaluation** (only 1 site changes per step)
- **Simple, deterministic algorithm** suitable for small systems
- **Parallel evaluation** of multiple independent instances
- **Minimal memory footprint** (16-bit masks, sparse coupling)

This complements the existing Monte Carlo solvers by handling the small-system regime where exhaustive enumeration is both feasible and desirable for accuracy.

---

# Part 2: Application - Logic Gate Design via Genetic Algorithm

## Overview

The MQCA Gray code solver is designed for **massively parallel evaluation** of thousands of small molecular clusters to find designs that implement binary logic operations (AND, OR, XOR, NAND, etc.). This section describes the application layer: how to define inputs/outputs, evaluate logic tables, and search the design space.

## Core Concept: MQCA as Logic Device

### Inputs
- **Fixed boundary sites** (not part of active simulation)
- Input encoded by biasing neighboring active sites:
  - Input = 1: Apply bias `+W` to adjacent active site (simulating occupied neighbor)
  - Input = 0: Apply bias `0` or `-W` (simulating unoccupied neighbor)
- Two input bits (A, B) → 4 combinations: (0,0), (0,1), (1,0), (1,1)

### Outputs
- One designated **active site** serves as output
- After ground-state computation, read its occupancy `n_out ∈ {0, 1}`
- This gives the logical result for each input combination

### Logic Table
For 2 inputs, we compute ground state for all 4 combinations:

| Input A | Input B | Ground State | Output n_out |
|---------|---------|--------------|--------------|
| 0 | 0 | `occ_00` | `n_out(00)` |
| 0 | 1 | `occ_01` | `n_out(01)` |
| 1 | 0 | `occ_10` | `n_out(10)` |
| 1 | 1 | `occ_11` | `n_out(11)` |

The 4 output bits define which logic function is implemented:
- `0000` = FALSE
- `0001` = AND (only 11 → 1)
- `0111` = OR (only 00 → 0)
- `0110` = XOR (01,10 → 1; 00,11 → 0)
- `1110` = NAND
- `1111` = TRUE
- etc.

## System Architecture

### Site Types

```
[Input A] ──┐
            ├──→ [Active cluster: nSite active sites] ──→ [Output site (one of active)]
[Input B] ──┘
             ↑
    Fixed boundary sites (not in simulation)
    Apply bias ε_bias = W * n_input to adjacent active sites
```

| Site Type | Role in Simulation | Representation |
|-----------|-------------------|----------------|
| Input | Fixed boundary, not simulated | Bias applied to neighbors |
| Active | Full ground-state computation | Binary variable `n_i ∈ {0,1}` |
| Output | One selected active site | Read `n_out` from ground state |

### Bias Calculation (Python side)

```python
def apply_input_bias(Esite_base, input_sites, input_values, W_coupling):
    """
    Modify on-site energies to encode input bits.
    
    Parameters
    ----------
    Esite_base : (nSite,) base on-site energies
    input_sites : list of (site_idx, neighbor_idx, W) tuples
        Defines which active sites neighbor input pads
    input_values : (2,) binary values for inputs A, B
    W_coupling : float, coupling strength to input pads
    
    Returns
    -------
    Esite : (nSite,) modified energies including input bias
    """
    Esite = Esite_base.copy()
    for inp_idx, (site_idx, neighbor_idx, W) in enumerate(input_sites):
        input_val = input_values[inp_idx]  # 0 or 1
        Esite[site_idx] += W * input_val   # Bias from fixed input neighbor
    return Esite
```

## Task 1: Parameter Scan for Fixed Geometry

### Goal
For a fixed cluster shape, find which coupling parameters (W1, W2) produce useful logic gates.

### Geometry Setup
- Square lattice arrangement of active sites
- Two coupling types:
  - **W1**: Cartesian neighbors (horizontal/vertical)
  - **W2**: Diagonal neighbors
- Fixed ε0 = 1.0 for all active sites (without input biases)

### Parallelization Strategy

```python
# W1-W2 parameter grid
nW1, nW2 = 100, 100
W1_vals = np.linspace(-2.0, 2.0, nW1)
W2_vals = np.linspace(-2.0, 2.0, nW2)
W1_grid, W2_grid = np.meshgrid(W1_vals, W2_vals)

# Each (W1, W2) point → 4 input combinations
nInstances = nW1 * nW2 * 4  # Total parallel evaluations
```

### GPU Batch Structure

| Instance Index | W1 | W2 | Input A | Input B | Purpose |
|---------------|----|----|---------|---------|---------|
| 0 | W1[0] | W2[0] | 0 | 0 | Logic table entry |
| 1 | W1[0] | W2[0] | 0 | 1 | Logic table entry |
| 2 | W1[0] | W2[0] | 1 | 0 | Logic table entry |
| 3 | W1[0] | W2[0] | 1 | 1 | Logic table entry |
| 4 | W1[0] | W2[1] | 0 | 0 | Next (W1,W2) point |
| ... | ... | ... | ... | ... | ... |

Total: `nW1 * nW2 * 4` independent ground-state computations.

### Post-Processing: Logic Identification

```python
def identify_logic(outputs_4):
    """
    outputs_4: [n_out(00), n_out(01), n_out(10), n_out(11)]
    Returns: logic_type, is_valid
    """
    out_bits = ''.join(str(int(o)) for o in outputs_4)
    logic_table = {
        '0000': 'FALSE', '0001': 'AND', '0010': 'A<B', '0011': 'A',
        '0100': 'A>B', '0101': 'B', '0110': 'XOR', '0111': 'OR',
        '1000': 'NOR', '1001': 'XNOR', '1010': 'NOT_B', '1011': 'A≥B',
        '1100': 'NOT_A', '1101': 'A≤B', '1110': 'NAND', '1111': 'TRUE'
    }
    return logic_table.get(out_bits, 'UNKNOWN')

# Aggregate results into phase diagram
phase_diagram = np.zeros((nW2, nW1), dtype=int)  # 0-15 encoding logic type
for i in range(nW1):
    for j in range(nW2):
        outputs = results[i, j, :]  # 4 outputs for this (W1,W2)
        phase_diagram[j, i] = encode_logic(outputs)
```

### Visualization

```python
def plot_phase_diagram(W1_vals, W2_vals, phase_diagram, cluster_shape):
    """Plot W1-W2 phase diagram showing which logic gates appear."""
    fig, ax = plt.subplots()
    im = ax.imshow(phase_diagram, origin='lower', extent=[W1_vals[0], W1_vals[-1], W2_vals[0], W2_vals[-1]],
                   cmap='tab20', vmin=0, vmax=15)
    ax.set_xlabel('W1 (Cartesian coupling)')
    ax.set_ylabel('W2 (Diagonal coupling)')
    ax.set_title(f'Logic Phase Diagram: {cluster_shape}')
    
    # Add colorbar with logic labels
    cbar = plt.colorbar(im, ax=ax, ticks=range(16))
    cbar.ax.set_yticklabels([logic_names[i] for i in range(16)])
    
    plt.savefig('phase_diagram.png', dpi=150)

def visualize_cluster_groundstates(cluster_pos, input_sites, output_site, 
                                     ground_states_4, outputs_4, logic_name):
    """
    Visualize the cluster and its 4 ground state configurations.
    
    Parameters
    ----------
    cluster_pos : (nSite, 2) positions of active sites
    input_sites : list of (pos, idx) for input pads
    output_site : int, index of output site
    ground_states_4 : (4, nSite) occupancy for each input combination
    outputs_4 : (4,) output values
    logic_name : str, identified logic function
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    titles = ['Input (0,0)', 'Input (0,1)', 'Input (1,0)', 'Input (1,1)']
    
    for ax, occ, out, title in zip(axes.flat, ground_states_4, outputs_4, titles):
        # Plot active sites
        for i, (x, y) in enumerate(cluster_pos):
            color = 'red' if occ[i] else 'blue'
            marker = 's' if i == output_site else 'o'
            size = 200 if i == output_site else 100
            ax.scatter(x, y, c=color, marker=marker, s=size, edgecolors='black')
        
        # Plot input sites (fixed, not in simulation)
        for inp_pos, inp_idx in input_sites:
            ax.scatter(*inp_pos, c='green', marker='^', s=300, edgecolors='black')
        
        # Draw couplings
        for i in range(len(cluster_pos)):
            for j in range(i+1, len(cluster_pos)):
                if W[i,j] != 0:
                    ax.plot([cluster_pos[i,0], cluster_pos[j,0]], 
                           [cluster_pos[i,1], cluster_pos[j,1]], 'k-', alpha=0.3)
        
        ax.set_title(f'{title} → Output: {int(out)}')
        ax.set_aspect('equal')
        ax.axis('off')
    
    plt.suptitle(f'Logic: {logic_name}')
    plt.tight_layout()
    plt.savefig(f'cluster_logic_{logic_name}.png', dpi=150)
```

## Task 2: Geometry Optimization via Genetic Algorithm

### Goal
Find cluster geometries (arrangements of active sites on a grid) that implement target logic functions, while maintaining connectivity.

### Representation

```python
class ClusterGenome:
    """
    Genome encoding a cluster geometry.
    
    Fixed elements (not in genome):
    - Input site positions (2 fixed locations)
    - Output site index (must be one of the active sites)
    
    Variable elements (genome):
    - Active site positions on a grid
    - nActive: number of active sites (1-16)
    - positions: (nActive, 2) integer grid coordinates
    
    Constraints:
    - All active sites must be connected (path exists between any pair)
    - Input sites must neighbor at least one active site
    - Output site is one of the active sites
    """
    
    def __init__(self, nActive, positions, input_neighbors, output_idx):
        self.nActive = nActive                    # 1-16
        self.positions = positions              # (nActive, 2) int
        self.input_neighbors = input_neighbors  # Which active sites neighbor inputs
        self.output_idx = output_idx            # 0 to nActive-1
        
    def is_valid(self):
        """Check connectivity constraint."""
        # Build adjacency graph
        adj = self._build_adjacency()
        # BFS from site 0, check all reachable
        visited = self._bfs(adj, 0)
        return len(visited) == self.nActive
```

### Genetic Algorithm Structure

```python
def run_genetic_algorithm(n_generations=100, population_size=1000,
                         target_logic='XOR', W1_range=(-2,2), W2_range=(-2,2)):
    """
    Genetic algorithm to find cluster geometries that implement target logic.
    
    Fitness function combines:
    1. Logic correctness: does it implement target function?
    2. Robustness: how large is the (W1,W2) region where it works?
    3. Simplicity: fewer sites is better
    """
    
    # Initialize population
    population = [random_cluster() for _ in range(population_size)]
    
    for gen in range(n_generations):
        # === EVALUATION PHASE ===
        # Each individual needs to be tested at multiple (W1,W2) points
        # → Massive batch evaluation on GPU
        
        all_instances = []
        individual_map = []  # Track which instances belong to which individual
        
        for ind_idx, ind in enumerate(population):
            # Sample W1,W2 points for this individual
            W1_samples = np.random.uniform(*W1_range, 10)
            W2_samples = np.random.uniform(*W2_range, 10)
            
            for W1 in W1_samples:
                for W2 in W2_samples:
                    # 4 input combinations per (W1,W2) point
                    for inp_A in [0, 1]:
                        for inp_B in [0, 1]:
                            # Build system instance
                            Esite = ind.build_Esite(W1, W2, inp_A, inp_B)
                            W_sparse = ind.build_W_sparse(W1, W2)
                            all_instances.append((Esite, W_sparse))
                            individual_map.append(ind_idx)
        
        # === GPU BATCH EVALUATION ===
        # Send all instances to GPU in one batch
        E_mins, occ_mins = solver.solve_batch(all_instances)
        
        # === FITNESS COMPUTATION ===
        fitness = np.zeros(len(population))
        for ind_idx in range(len(population)):
            instances = [i for i, idx in enumerate(individual_map) if idx == ind_idx]
            outputs = occ_mins[instances]  # Extract output bits
            
            # Check if target logic is implemented
            correct_fraction = compute_logic_match(outputs, target_logic)
            
            # Bonus for larger valid (W1,W2) region
            robustness = compute_robustness(outputs, W1_samples, W2_samples)
            
            # Penalty for complexity
            complexity_penalty = 0.1 * population[ind_idx].nActive
            
            fitness[ind_idx] = correct_fraction + 0.5 * robustness - complexity_penalty
        
        # === SELECTION, CROSSOVER, MUTATION ===
        parents = tournament_selection(population, fitness, n_parents=population_size//2)
        offspring = crossover_and_mutate(parents)
        population = parents + offspring
        
        # Logging
        best_idx = np.argmax(fitness)
        print(f"Gen {gen}: Best fitness = {fitness[best_idx]:.3f}, "
              f"Logic = {identify_logic(population[best_idx])}, "
              f"nSites = {population[best_idx].nActive}")
    
    return population[best_idx]
```

### Mutation Operators

```python
def mutate_cluster(cluster, mutation_rate=0.1):
    """
    Apply random mutations while preserving connectivity.
    
    Mutation types:
    1. Move: Shift one active site to adjacent empty grid position
    2. Add: Insert new active site adjacent to existing cluster
    3. Remove: Delete active site (must not disconnect cluster)
    4. Swap output: Change which active site is the output
    """
    mutant = cluster.copy()
    
    if random.random() < mutation_rate:
        mutation_type = random.choice(['move', 'add', 'remove', 'swap_output'])
        
        if mutation_type == 'move':
            # Move random site to adjacent empty position
            site_idx = random.randint(0, mutant.nActive-1)
            old_pos = mutant.positions[site_idx]
            neighbors = get_adjacent_positions(old_pos)
            valid_neighbors = [p for p in neighbors if p not in mutant.positions]
            if valid_neighbors:
                mutant.positions[site_idx] = random.choice(valid_neighbors)
                
        elif mutation_type == 'add' and mutant.nActive < 16:
            # Add site at boundary
            boundary = get_boundary_positions(mutant.positions)
            new_pos = random.choice(boundary)
            mutant.positions = np.vstack([mutant.positions, new_pos])
            mutant.nActive += 1
            
        elif mutation_type == 'remove' and mutant.nActive > 2:
            # Remove site (check connectivity first)
            removable = [i for i in range(mutant.nActive) 
                          if can_remove_without_disconnect(i, mutant)]
            if removable:
                remove_idx = random.choice(removable)
                mutant.positions = np.delete(mutant.positions, remove_idx, axis=0)
                mutant.nActive -= 1
                # Adjust output_idx if needed
                if mutant.output_idx == remove_idx:
                    mutant.output_idx = 0
                elif mutant.output_idx > remove_idx:
                    mutant.output_idx -= 1
                    
        elif mutation_type == 'swap_output':
            mutant.output_idx = random.randint(0, mutant.nActive-1)
    
    return mutant if mutant.is_valid() else cluster
```

### Massive Parallelism Architecture

```
GPU Work Distribution:

┌─────────────────────────────────────────────────────────────┐
│  Grid of Compute Units (e.g., 80 CUs on modern GPU)        │
│                                                             │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐          │
│  │ WorkGroup 0 │  │ WorkGroup 1 │  │ WorkGroup 2 │  ...     │
│  │             │  │             │  │             │          │
│  │ Instance 0  │  │ Instance 16 │  │ Instance 32 │          │
│  │ (W1=0.5,    │  │ (W1=0.5,    │  │ (W1=0.6,    │          │
│  │  W2=-1.2,   │  │  W2=-1.1,   │  │  W2=-1.2,   │          │
│  │  inp=(0,0)) │  │  inp=(0,1)) │  │  inp=(0,0)) │          │
│  │             │  │             │  │             │          │
│  │ 16 threads  │  │ 16 threads  │  │ 16 threads  │          │
│  │ → Gray code │  │ → Gray code │  │ → Gray code │          │
│  │   search    │  │   search    │  │   search    │          │
│  └─────────────┘  └─────────────┘  └─────────────┘          │
│                                                             │
│  Total: ~10,000+ instances running in parallel              │
│  Each: 65536 Gray code steps × 16 threads                   │
└─────────────────────────────────────────────────────────────┘
```

Batch size calculation:
- 1000 individuals
- 10 × 10 (W1, W2) samples each
- 4 input combinations per sample
- **Total: 400,000 instances per generation**

With 65536 iterations per instance at ~1-2 ns per iteration, this is tractable on modern GPUs.

### Fitness Landscape Visualization

```python
def visualize_fitness_landscape(population, fitness, gen):
    """Plot population in PCA-reduced geometry space, colored by fitness."""
    # Extract features: nActive, centroid positions, aspect ratio, etc.
    features = np.array([extract_features(ind) for ind in population])
    
    # PCA to 2D
    pca = PCA(n_components=2)
    coords = pca.fit_transform(features)
    
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(coords[:, 0], coords[:, 1], c=fitness, cmap='viridis', s=50)
    plt.colorbar(scatter, label='Fitness')
    plt.title(f'Generation {gen}: Fitness Landscape in Geometry Space')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.savefig(f'fitness_landscape_gen{gen}.png', dpi=150)
    plt.close()
```

## Implementation Roadmap

### Phase 1: Basic Parameter Scan (Week 1-2)
1. Implement kernel `mqca_graycode_groundstate`
2. Python wrapper `MQCAGraySolver`
3. W1-W2 parameter scan for fixed cluster shapes
4. Phase diagram visualization
5. Logic identification and validation

### Phase 2: Single Cluster Visualization (Week 2-3)
1. Visualizer for cluster geometry with input/output sites
2. Ground-state configuration plots for 4 input combinations
3. Logic table display
4. Animation of Gray code traversal (optional, for debugging)

### Phase 3: Genetic Algorithm (Week 3-4)
1. Cluster genome representation
2. Validity checking (connectivity)
3. Mutation operators
4. Batch evaluation pipeline
5. Fitness function design
6. Selection and breeding

### Phase 4: Analysis Tools (Week 4-5)
1. Fitness landscape visualization
2. Convergence diagnostics
3. Best solution analysis
4. Robustness quantification
5. Export to simulation formats

## Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Input encoding | Fixed boundary sites → bias | Natural for MQCA physics |
| Output read | Single active site occupancy | Direct binary readout |
| Logic test | All 4 input combinations | Exhaustive verification |
| W1-W2 sampling | Random or grid | Grid for phase diagrams, random for GA |
| GA population | 1000-5000 | Balance diversity vs. evaluation cost |
| Genome | Grid positions + output index | Simple, intuitive, easy to mutate |
| Constraint | Connectivity required | Physical realism |
| Batch size | 400K+ instances | Maximize GPU utilization |

## Success Metrics

1. **Logic correctness**: Device implements target logic for ≥90% of (W1,W2) samples
2. **Robustness**: Valid region in W1-W2 space spans ≥0.5×0.5 range
3. **Efficiency**: Fewer than 8 active sites preferred
4. **Novelty**: Non-trivial geometries (not just simple chains)

---

## File Structure (Proposed)

```
pyBall/OCL/
├── cl/
│   ├── hubbard.cl              # Existing kernels
│   ├── mqca_gray.cl            # NEW: Gray code MQCA solver
│   └── mqca_logic.cl           # NEW: Logic-specific kernels (optional)
├── MQCAGraySolver.py           # Base solver
├── MQCALogic.py                # NEW: Logic gate evaluation, input encoding
├── MQCAGenetic.py              # NEW: Genetic algorithm infrastructure
└── MQCAVis.py                  # NEW: Visualization tools

tests/tMQCA/
├── test_gray_solver.py         # Unit tests for Gray code solver
├── test_logic_gates.py         # Test logic identification
├── test_genetic.py             # Test GA components
├── demo_phase_diagram.py       # W1-W2 scan demo
├── demo_genetic_xor.py         # Genetic search for XOR gate
└── demo_visualizer.py          # Cluster visualization demo
```
