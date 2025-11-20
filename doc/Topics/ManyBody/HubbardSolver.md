# Hubbard Solver OpenCL Kernels and Buffer Management

This document provides an overview of the OpenCL kernels defined in `cl/hubbard.cl` and their interaction with the buffer management in `pyBall/OCL/HubbardSolver.py`.

## OpenCL Kernels Overview

### 1. `eval_coupling_matrix`
**Description:** Evaluates the full coupling matrix between all sites, including mirror image interactions.
**Arguments:**
*   `nSingle`: Number of single sites.
*   `posE`: Global buffer for site positions and initial energies (`float4`).
*   `Rij`: Global buffer for output distance matrix (`float`).
*   `Wij`: Global buffer for output coupling matrix (`float`).
*   `zMirror`: Z-coordinate of the mirror plane (`float`).
**Buffers Used:** `posE_buff`, `Rij_buff`, `Wij_buff`

### 2. `eval_Oriented_Hopping`
**Description:** Evaluates tunneling amplitudes for each site-tip pair.
**Arguments:**
*   `nSingle`: Number of single sites.
*   `nTips`: Number of tip positions.
*   `posE`: Global buffer for site positions (`float4`).
*   `rots`: Global buffer for site rotation as unit complex number (`float2`).
*   `orbs`: Global buffer for orbital parameters (`float4`).
*   `pTips`: Global buffer for tip positions (`float4`).
*   `Tout`: Global buffer for output tunneling amplitudes (`float`).
**Buffers Used:** `posE_buff`, `rots_buff`, `orbs_buff`, `pTips_buff`, `Tout_buff`

### 3. `solve_minBrute_fly`
**Description:** Finds the minimum energy many-body state using a brute-force search with on-the-fly W_ij calculation.
**Arguments:**
*   `nSingle`: Number of single sites.
*   `posE`: Global buffer for site positions and initial energies (`float4`).
*   `nTips`: Number of tip positions.
*   `pTips`: Global buffer for tip positions (`float4`).
*   `tipDecay`: Decay parameter for tip potential (`float`).
*   `tipRadius`: Radius of the tip (`float`).
*   `zMirror`: Z-coordinate of the mirror plane (`float`).
*   `Wcouple`: Coupling strength parameter (`float`).
*   `Emin`: Global buffer for output minimum energies (`float`).
*   `iMin`: Global buffer for output bitmask of ground-state occupation (`int`).
*   `Itot`: Global buffer for output total currents (`float2`).
**Buffers Used:** `posE_buff`, `pTips_buff`, `Emin_buff`, `iMin_buff`, `Itot_buff`

### 4. `solve_minBrute_boltzmann`
**Description:** Highly efficient single-pass kernel for Boltzmann-weighted averaging of states.
**Arguments:**
*   `nSingle`: Number of single sites.
*   `posE`: Global buffer for site positions and initial energies (`float4`).
*   `nTips`: Number of tip positions.
*   `pTips`: Global buffer for tip positions (`float4`).
*   `tipDecay`: Decay parameter for tip potential (`float`).
*   `tipRadius`: Radius of the tip (`float`).
*   `zMirror`: Z-coordinate of the mirror plane (`float`).
*   `Wcouple`: Coupling strength parameter (`float`).
*   `bBoltzmann`: Inverse temperature (1/kT) (`float`).
*   `Emin`: Global buffer for output minimum energies (`float`).
*   `iMin`: Global buffer for output bitmask of ground-state occupation (`int`).
*   `Itot`: Global buffer for output total currents (`float2`).
**Buffers Used:** `posE_buff`, `pTips_buff`, `Emin_buff`, `iMin_buff`, `Itot_buff`

### 5. `solve_local_updates`
**Description:** Implements an iterative local update algorithm, capable of greedy minimization or simulated annealing (Monte Carlo), to find low-energy configurations of site occupancies. Each work-group processes a single STM tip position, and threads within the group collaborate.

**Parallelization Strategy:**
*   Each work-group (identified by `get_group_id(0)`) is assigned to a unique STM tip position.
*   Within each work-group, threads (`get_local_id(0)`) collaborate to perform iterative updates for that tip.

**Configuration Storage:**
*   Site occupancies (`nSite`) are stored as a bitmask in `__local uchar occ_mask[OCC_BYTES]`. This is highly memory-efficient, with each bit representing the occupancy (0 or 1) of a site.
*   `GET_OCC` and `FLIP_OCC` macros are used for efficient bitwise access and modification of `occ_mask`.

**Local Memory Usage:**
*   `occ_mask`: Stores the current configuration being optimized.
*   `reduction_dE`, `reduction_site`: Used for work-group-wide reduction to find the best proposed local move. Threads store their calculated energy changes (`dE`) and corresponding site indices here.
*   Global buffers (`Esite`, `Tsite`, `W_val`, `W_idx`, `nNeigh`) provide pre-calculated on-site energies, tunneling factors, and interaction parameters.

**Iterative Optimization Loop (`for (int iter = 0; iter < nIter; ++iter)`):**
1.  **Step A: Propose a Move (Parallel):**
    *   Each thread selects a site (`i_site`) to consider flipping its occupancy. Selection can be deterministic (sweeping through sites) or stochastic (random selection).
    *   The energy change (`dE`) resulting from flipping `i_site`'s occupancy is calculated, considering on-site energy and interactions with neighbors.
    *   `dE` and `i_site` are stored in `reduction_dE[lid]` and `reduction_site[lid]`.
    *   A `barrier(CLK_LOCAL_MEM_FENCE)` ensures all proposals are ready.
2.  **Step B: Find Best Move and Decide (Thread 0):**
    *   Only the first thread (`lid == 0`) in the work-group aggregates proposals.
    *   It finds the move with the minimum `dE` (greatest energy decrease).
    *   **Decision Logic:**
        *   If `dEmin < 0.0f` (energy decreases), the move is accepted (greedy minimization).
        *   If `solverMode == 2` (Simulated Annealing) and `dEmin >= 0.0f`, the move is accepted probabilistically based on `exp(-dEmin / kT)`, allowing escape from local minima.
    *   If accepted, `FLIP_OCC` updates `occ_mask`.
    *   A `barrier(CLK_LOCAL_MEM_FENCE)` ensures the updated `occ_mask` is visible for the next iteration.

**Finalization:**
*   After `nIter` iterations, `lid == 0` calculates the final total energy, occupied current, and unoccupied current for the optimized configuration.
*   The optimized `occ_mask`, total energy, and currents are written back to global memory (`occ_out`, `E_out`, `Itot_out`).

### 6. `solve_MC_neigh`
**Description:** Implements a hybrid Monte Carlo (MC) algorithm for finding low-energy configurations, incorporating global mutation strategies (load best, random reset, load neighbor's best) followed by local greedy descent. This kernel is designed for a single phase of the MC process.

**Parallelization Strategy:**
*   Each work-group (identified by `get_group_id(0)`) is assigned to a unique STM tip position (`itip`).
*   Threads (`get_local_id(0)`) within a work-group collaborate on the optimization for that single tip.

**Configuration Storage:**
*   Site occupancies (`nSite`) are stored as a bitmask in `__local uchar occ_mask[OCC_BYTES]`. This is memory-efficient, with each bit representing the occupancy (0 or 1) of a site.

**Local Memory Usage:**
*   `occ_mask`: Stores the current configuration being optimized.
*   `reduction_buffer`, `reduction_site`: Used for work-group-wide reduction to find the best proposed local move.

**Iterative Optimization Loop (`for (int iGlobal = 0; iGlobal < nGlobalSteps; ++iGlobal)`):**
1.  **Step A: Generate a Trial Configuration (Global Mutation):**
    *   Only `lid == 0` selects a global mutation strategy based on `prob_params` (probabilities for loading the best known configuration, random reset, or loading a neighbor's best configuration).
    *   The `occ_mask` is updated with the chosen trial configuration.
    *   A `barrier(CLK_LOCAL_MEM_FENCE)` ensures all threads see the new `occ_mask`.
2.  **Step B: Inner Loop for Local Descent (`for (int iLocal = 0; iLocal < nLocalIter; ++iLocal)`):**
    *   All threads collaborate to perform `nLocalIter` steps of greedy local updates, similar to `solve_local_updates`.
    *   Each thread proposes a site flip, and `lid == 0` executes the best move.
3.  **Step C: Evaluate Relaxed Configuration and Compare with Global Best:**
    *   The total energy of the relaxed configuration (`E_current`) is calculated using `calculate_total_energy_parallel`.
    *   If `E_current` is better than the `E_best_cached` (best energy found so far for this tip), `E_best` and `occ_best` in global memory are updated.
    *   A `barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE)` ensures global writes are visible.

**Finalization:**
*   After all global steps, `lid == 0` reads the final best state from global memory (`occ_best`) and calculates the total tunneling currents (`Itot_out`).

### 7. `solve_MC_2phase`
**Description:** This kernel implements a two-phase Monte Carlo optimization strategy, using a ping-pong buffering approach to alternate between two sets of best-known configurations (`occ_best_A`/`E_best_A` and `occ_best_B`/`E_best_B`). This allows for continuous optimization without overwriting the current best state until a better one is found.

**Parallelization Strategy:**
*   Each work-group (`get_group_id(0)`) is assigned to a unique STM tip position (`itip`).
*   Threads (`get_local_id(0)`) within a work-group collaborate on the optimization for that single tip.

**Configuration Storage:**
*   Site occupancies (`nSite`) are stored as a bitmask in `__local uchar occ_mask[OCC_BYTES]`.

**Local Memory Usage:**
*   `occ_mask`: Stores the current configuration being optimized.
*   `reduction_buffer`, `reduction_site`: Used for work-group-wide reduction to find the best proposed local move.

**Iterative Optimization Loop (`for (int iGlobal = 0; iGlobal < nGlobalSteps; ++iGlobal)`):**
*   This kernel is very similar in structure to `solve_MC_neigh`, but it explicitly manages two sets of global buffers for the best configurations (A and B). In each global step, it reads from one set (e.g., A) as `occ_in`/`E_in` and writes to the other set (e.g., B) as `occ_out`/`E_out`.
*   The `prob_params` are used for global mutation strategies (load best from `occ_in`/`E_in`, load neighbor's best from `occ_in`/`E_in`, or random reset).
*   An inner loop performs `nLocalIter` steps of greedy local descent.
*   The relaxed configuration's energy is compared to the `E_best_current` (read from `E_in`), and if better, the new configuration and energy are written to `occ_out`/`E_out`.

**Finalization:**
*   This kernel does not calculate final currents directly. It is expected that a separate kernel (e.g., `calculate_currents`) will be called after `solve_MC_2phase` completes all its global steps, using the final best configuration from either `occ_best_A` or `occ_best_B`.


---

## An Advanced Monte Carlo Framework for Hubbard Model Ground-State Search

### 1. Motivation: Overcoming the Limitations of Greedy Descent

Solving the Hubbard model for a large number of sites, particularly under the influence of a spatially varying external potential (like an STM tip), is a computationally demanding task. Simple optimization algorithms, such as greedy local descent, are fast but suffer from a critical flaw: they are highly susceptible to getting trapped in local energy minima. For systems with strong inter-site coupling (`Wij`), the energy landscape is rugged, and a simple greedy search is almost certain to find a suboptimal solution.

Furthermore, a purely local solver is non-ergodic; it cannot guarantee that all possible configurations can be reached from an arbitrary starting point. To find the true ground state with high probability, a more robust and explorative algorithm is required.

The framework presented here is a **hybrid metaheuristic** that combines the speed of a greedy local search with a set of powerful global "mutation" strategies. It is designed to efficiently explore the vast configuration space, escape local minima, and leverage the spatial correlations inherent in scanning probe microscopy problems. The algorithm can be classified as a form of **Iterated Local Search (ILS)**.

### 2. The Hybrid Algorithm: A Multi-Strategy Approach

The core of the algorithm is a two-level loop structure: an outer "global" loop for exploration and an inner "local" loop for fast exploitation.

#### 2.1. The Framework: Perturbation Followed by Local Descent

Instead of a simple Metropolis-style acceptance of single moves, our algorithm operates on entire configurations. Each step of the outer loop consists of three phases:

1.  **Perturbation (Global Mutation):** A *trial configuration* is generated using one of several powerful global mutation strategies. This step is designed to move the search to a new region of the configuration space.
2.  **Local Descent (Relaxation):** The newly generated trial configuration is immediately subjected to a fixed number of iterations (`nLocalIter`) of a fast, greedy, local-update optimizer. This quickly relaxes the trial state into a nearby local minimum.
3.  **Selection (Acceptance Criterion):** The energy of the *relaxed* trial configuration is calculated. If this energy is lower than the best energy found so far for the given tip position (stored in global memory), the new configuration and its energy replace the previous best.

The power of this "relax-then-compare" approach is that it allows the algorithm to consider trial states that may initially have very high energy. A randomly generated state, for instance, is likely to have a high energy, but after a few steps of local descent, it may relax into a very deep minimum, potentially even the global one. A simple Metropolis algorithm would have rejected such a high-energy initial move outright.

#### 2.2. Global Mutation Strategies

The choice of the perturbation operator is critical for the algorithm's success. We employ a probabilistic mixture of three distinct strategies, each serving a specific purpose:

**1. Strategy: Reload Best (Exploitation)**
*   **Action:** The current best-known configuration for the given tip position is loaded from global memory.
*   **Rationale:** This is an exploitative move. It allows the algorithm to attempt further refinement of the best solution found so far. If the local descent optimizer finds an even better state starting from the current best, progress is made. This is the most conservative strategy.

**2. Strategy: Random Reset (Ergodicity & Exploration)**
*   **Action:** A completely new configuration is generated by assigning a random occupancy (0 or 1) to every site.
*   **Rationale:** This is the primary mechanism for ensuring ergodicity and global exploration. It can move the search to any point in the configuration space, providing a powerful means of escaping deep, deceptive local minima. An algorithm that periodically tries completely random configurations can, in principle, find the global optimum.

**3. Strategy: Load Neighbor's Best (Spatially-Informed Exploration)**
*   **Action:** The algorithm identifies a random neighboring tip position (pixel) in the 2D scan grid and loads its best-known configuration as the trial state.
*   **Rationale:** This is a novel, problem-specific heuristic that leverages the physical nature of STM simulations. Adjacent tip positions are likely to have similar ground-state configurations. By "borrowing" a successful solution from a neighbor, a work-group can effectively jump-start its search in a highly promising region of the configuration space. This has two major benefits:
    *   **Accelerated Convergence:** It dramatically speeds up the search by sharing information across parallel work-groups.
    *   **Noise Reduction:** It promotes spatial smoothness in the final result, reducing "salt-and-pepper" noise where adjacent pixels would otherwise be stuck in different local minima.

The choice between these three strategies is controlled by a probability vector `(p_load_best, p_random_reset, p_load_neighbor)`.

### 3. OpenCL Kernel Implementation: `solve_MC_neigh` and `solve_MC_2phase`

The advanced Monte Carlo algorithm is now implemented using two primary OpenCL kernels: `solve_MC_neigh` and `solve_MC_2phase`. These kernels work in conjunction to provide a robust and efficient search for ground states, particularly suitable for large-scale simulations and those requiring a multi-phase optimization approach.

`solve_MC_neigh` implements the core hybrid Monte Carlo step, performing global mutations followed by local descent. `solve_MC_2phase` extends this by enabling a ping-pong buffering strategy, allowing for continuous optimization over many iterations while maintaining the best-found configurations.

#### 3.1. `solve_MC_neigh` Kernel Signature & Parameters

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `nSite` | `int` | Number of sites in the Hubbard model. |
| `nTips` | `int` | Total number of tip positions (i.e., independent simulations). |
| `solver_params` | `float4` | Vector of general solver parameters: `{kT, nIter, solverMode, seed}`. `kT` is unused but kept for compatibility. |
| `nLocalIter` | `int` | The number of greedy local descent steps to perform after each global mutation. |
| `prob_params` | `float4` | **Crucial:** Probabilities for global moves: `{p_load_best, p_load_neighbor, p_random_reset, unused}`. Note the order is `p_load_best`, `p_load_neighbor`, `p_random_reset` after cumulative sum. |
| `nx` | `int` | **Crucial:** The width of the 2D simulation grid, required for the "Load Neighbor's Best" strategy. |
| `Esite`, `Tsite`, ... | `__global` | Pointers to global memory buffers containing system parameters (energies, couplings, etc.). |
| `occ_best` | `__global uchar*`| **State:** Global buffer storing the bitmask of the best configuration found for each tip. |
| `E_best` | `__global float*` | **State:** Global buffer storing the energy of the best configuration for each tip. |
| `Itot_out` | `__global float2*`| Output buffer for final calculated currents. |
| `rng_seeds` | `__global const uint*` | Buffer containing initial random number generator seeds for each work-group. |
| `max_neighs` | `int` | Maximum number of neighbors considered for interaction. |
| `initMode` | `int` | Initialization mode for site occupancies. |

#### 3.2. `solve_MC_2phase` Kernel Signature & Parameters

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `nSite` | `int` | Number of sites in the Hubbard model. |
| `nTips` | `int` | Total number of tip positions. |
| `solver_params` | `float4` | Dummy parameter for compatibility. |
| `nLocalIter` | `int` | Number of greedy local descent steps. |
| `nGlobalSteps` | `int` | Number of global mutation steps in this phase. |
| `prob_params` | `float4` | Probabilities for global moves (cumulative sum). |
| `nx` | `int` | Width of the 2D simulation grid. |
| `occ_in`, `E_in` | `__global` | Input buffers for current best configurations and energies (from previous phase). |
| `occ_out`, `E_out` | `__global` | Output buffers for updated best configurations and energies (for next phase). |
| `Esite`, `W_val`, `W_idx`, `nNeigh` | `__global` | System parameters. |
| `rng_seeds` | `__global const uint*` | Random number generator seeds. |
| `max_neighs` | `int` | Maximum number of neighbors. |

#### 3.3. Parallelization and Memory
*   **Work-Group:** Each work-group, identified by `get_group_id(0)`, is assigned to a single, unique tip position (`itip`).
*   **Work-Item:** Threads within a work-group, identified by `get_local_id(0)`, collaborate on the optimization for that single tip.
*   **Local Memory:** To maximize the number of sites (`nSite`) the kernel can handle, local memory usage is minimized. Only a single working copy of the configuration (`__local uchar occ_mask[]`) is stored per work-group. This is a deliberate design choice that prioritizes problem size over minimizing global memory latency, although some operation could be faster if we keep multiple copies of configuration in shared memory (e.g. the best). To minimize local memory footprint, the occupancy of pixels is stored as bitmask (each bit represents the occupancy of a site) and we use `GET_OCC` and `FLIP_OCC` macros for efficient bitwise access and modification of `occ_mask`.

#### 3.4. Execution Flow

Both `solve_MC_neigh` and `solve_MC_2phase` follow a similar iterative structure for each work-group:

1.  **Initialization:** The first thread (`lid == 0`) loads the best-known energy for its assigned tip (`E_best[itip]` or `E_in[itip]`) into a private register for fast access.
2.  **Outer Loop (`nGlobalSteps`):**
    a.  **Global Mutation (`lid == 0`):** A random number is generated to select one of the three mutation strategies based on `prob_params`. The `occ_mask` in local memory is overwritten with the chosen trial configuration (from global best, random, or neighbor's best).
    b.  **Synchronization:** A local memory barrier ensures all threads in the work-group see the new trial `occ_mask`.
    c.  **Inner Loop (Local Descent):** All threads collaborate for `nLocalIter` steps. In each step, every thread proposes flipping a single site, and the `lid == 0` thread identifies and executes the move with the greatest energy decrease (the most "greedy" move).
    d.  **Evaluation:** After local descent, all threads collaborate to calculate the total energy of the relaxed configuration in `occ_mask` using a parallel reduction.
    e.  **Selection (`lid == 0`):** The newly calculated energy is compared to the cached best energy. If it's better, the global `E_best` and `occ_best` buffers (for `solve_MC_neigh`) or `E_out` and `occ_out` (for `solve_MC_2phase`) are updated.
3.  **Finalization:**
    *   For `solve_MC_neigh`, after all global steps are complete, `lid == 0` calculates the final tunnelling currents corresponding to the definitive best configuration stored in `occ_best` and writes them to the output buffer.
    *   For `solve_MC_2phase`, finalization involves writing the best configuration to the output buffers (`occ_out`, `E_out`). The calculation of currents is typically handled by a separate kernel (`calculate_currents`) after the multi-phase optimization is complete.

### 4. Python Routines for Monte Carlo Solvers

The `HubbardSolver.py` module provides Python wrappers and utility functions to manage the OpenCL kernels, handle buffer allocation, and orchestrate the Monte Carlo simulation. These routines simplify the interaction with the low-level OpenCL API.

#### 4.1. Buffer Management: `realloc_mc_buffers`

*   **Function:** `HubbardSolver.realloc_mc_buffers(nSite: int, nTips: int, nMaxNeigh: int)`
*   **Description:** This method is responsible for (re-)allocating the necessary GPU buffers for the Monte Carlo solvers. It sets up buffers for system data (`Esite`, `Tsite`, `W_val`, `W_idx`, `nNeigh`), random number generator seeds (`rng_seeds`), and crucially, the ping-pong state buffers (`E_best_A`, `occ_best_A`, `E_best_B`, `occ_best_B`) for the two-phase approach. It also allocates the final output buffer for currents (`Itot_out`).

#### 4.2. Kernel Setup: `setup_solve_mc_neigh` and `setup_solve_mc_2phase`

*   **Function:** `HubbardSolver.setup_solve_mc_neigh(...)`
    *   **Description:** Prepares the arguments for the `solve_MC_neigh` OpenCL kernel. It takes parameters such as `nSite`, `nTips`, `nLocalIter`, `prob_params` (probabilities for global moves), `nx` (grid width for neighbor logic), and various system buffers. It returns the kernel arguments and the kernel object, ready for execution.

*   **Function:** `HubbardSolver.setup_solve_mc_2phase(...)`
    *   **Description:** Pre-assembles the arguments for the `solve_MC_2phase` OpenCL kernel. This function is designed to set up the ping-pong buffering, taking `nLocalIter`, `nGlobalSteps`, `prob_params`, `nx`, and the `_A` and `_B` buffer pairs as input. It returns the kernel arguments and the kernel object.

#### 4.3. Execution Orchestration: `solve_mc` and `solve_mc_2phase`

*   **Function:** `HubbardSolver.solve_mc(posE, pTips, W_val, W_idx, nNeigh, params={})`
    *   **Description:** This is the high-level Python method that orchestrates the execution of the `solve_MC_neigh` kernel. It handles the initial allocation of buffers using `realloc_mc_buffers`, uploads input data (site positions, tip positions, sparse interaction matrix), sets up the kernel using `setup_solve_mc_neigh`, executes the kernel on the GPU, and downloads the results (`E_best`, `occ_best`, `Itot_out`). This function is suitable for single-phase Monte Carlo runs.

*   **Function:** `HubbardSolver.solve_mc_2phase(posE, pTips, W_val, W_idx, nNeigh, params={})`
    *   **Description:** This method orchestrates the multi-phase Monte Carlo optimization using the `solve_MC_2phase` kernel. It manages the ping-pong buffer swapping, iteratively calling the `solve_MC_2phase` kernel multiple times (controlled by `n_phases` in `params`) to refine the solution. After all phases are complete, it calls the `calculate_currents` kernel to compute the final currents based on the best configuration found across all phases. This function is designed for more complex, long-running optimizations.

#### 4.4. Demo Function: `demo_local_update`

*   **Function:** `run_Hubbard.demo_local_update()`
*   **Description:** Located in `run_Hubbard.py`, this function serves as a comprehensive demonstration of how to use the local update and Monte Carlo solvers. It sets up a system (e.g., a grid of sites), pre-calculates the sparse interaction matrix, and then calls `solver.solve_mc` or `solver.solve_mc_2phase` to run the simulation. It often includes plotting and analysis of the results, showcasing the capabilities of the HubbardSolver. This is the primary entry point for users to experiment with and understand the Monte Carlo framework.


#### 3.2. Parallelization and Memory
*   **Work-Group:** Each work-group, identified by `get_group_id(0)`, is assigned to a single, unique tip position (`itip`).
*   **Work-Item:** Threads within a work-group, identified by `get_local_id(0)`, collaborate on the optimization for that single tip.
*   **Local Memory:** To maximize the number of sites (`nSite`) the kernel can handle, local memory usage is minimized. Only a single working copy of the configuration (`__local uchar occ_mask[]`) is stored per work-group.  This is a deliberate design choice that prioritizes problem size over minimizing global memory latency, although some operation could be faster if we keep multiple copies of configuration in shared memory (e.g. the best). To minimize local memory footprint, the occupancy of pixels is stored as bitmask (each bit represents the occupancy of a site) and we use `GET_OCC` and `FLIP_OCC` macros for efficient bitwise access and modification of `occ_mask`.

#### 3.3. Execution Flow

The kernel executes the hybrid algorithm as follows for each work-group:

1.  **Initialization:** The first thread (`lid == 0`) loads the best-known energy for its assigned tip (`E_best[itip]`) into a private register for fast access.
2.  **Outer Loop (`nGlobalSteps`):**
    a.  **Global Mutation (`lid == 0`):** A random number is generated to select one of the three mutation strategies based on `prob_params`. The `occ_mask` in local memory is overwritten with the chosen trial configuration (from global best, random, or neighbor's best). A helper function, `get_neighbor_tip_id`, handles the 2D-to-1D index conversion.
    b.  **Synchronization:** A local memory barrier ensures all threads in the work-group see the new trial `occ_mask`.
    c.  **Inner Loop (Local Descent):** All threads collaborate for `nLocalIter` steps. In each step, every thread proposes flipping a single site, and the `lid == 0` thread identifies and executes the move with the greatest energy decrease (the most "greedy" move). This mean at the end of loop we have explored `nLocalIter*workgroup_size` configurations.
    d.  **Evaluation:** After local descent, all threads collaborate to calculate the total energy of the relaxed configuration in `occ_mask` using a parallel reduction.
    e.  **Selection (`lid == 0`):** The newly calculated energy is compared to the cached best energy. If it's better, the global `E_best` and `occ_best` buffers are updated.
3.  **Finalization:** After all global steps are complete, `lid == 0` calculates the final tunnelling currents corresponding to the definitive best configuration stored in `occ_best` and writes them to the output buffer.

### 4. Guidelines for Implementation

*   **Parameter Tuning:**
    *   `nLocalIter`: This determines the "depth" of the local search. I controls the ammount of effort we spend on new configurations before we compera them to champion. But it also controls the ratio between the local and global moves. Therefore relation to exploitation and exploration is double edge sword.
    *   `prob_params`: A good starting point is a balance between exploration and exploitation, e.g., `(p_load_best, p_random_reset, p_load_neighbor)=(0.10, 0.45, 0.45)`. This gives a small chance to refine the best, and an equal, large chance to either reset randomly or learn from a neighbor.
*   **Host-Side Logic:**
    *   The Python wrapper must correctly pass the image width `nx` to the kernel.
    *   Before the first run, the `E_best` buffer should be initialized to infinity (`np.inf`) to ensure the first result is always accepted. For subsequent runs (continuation), this buffer should be preserved.