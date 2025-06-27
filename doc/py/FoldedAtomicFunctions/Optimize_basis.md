# Basis Set Optimizer for Folded Atomic Functions

This document outlines the principles, implementation details, and considerations for the basis set optimizer found in `Optimize_basis_new.py`. The goal of this optimizer is to find an optimal set of basis functions for approximating a family of target functions (e.g., Morse potentials) by minimizing a cost function that balances accuracy and computational performance.

## 1. Core Principle: Monte Carlo Optimization

The optimizer employs a **Simulated Annealing-like Monte Carlo** approach. This is a metaheuristic algorithm used to find a good approximation to the global optimum of a given function in a large search space.

Key features:
- **Stochastic Search:** It explores the solution space by making random changes (mutations) to a current solution.
- **Fitness Evaluation:** Each new candidate solution is evaluated using a cost function. Lower cost indicates a better solution.
- **Acceptance Criteria:**
    - Better solutions (lower cost) are always accepted.
    - Worse solutions (higher cost) can be accepted with a certain probability, which depends on the current "temperature" and the difference in cost. This allows the optimizer to escape local minima.
- **Cooling Schedule:** The temperature gradually decreases over iterations, reducing the probability of accepting worse solutions, thus converging towards an optimum.

The generic optimizer logic is encapsulated in `Optimizer.py`.

## 2. Problem Specifics: Basis Set for Morse Potentials

In `Optimize_basis_new.py`, we are specifically optimizing a basis set composed of functions of the form $(1 - (z-z_0)/z_{span})^n$, where $z_{span} = z_{cut} - z_0$. These are derived from `cutoff_poly_basis()` in `basis_utils.py`.

### 2.1. Basis Definition

A basis set is defined as a list of tuples:
`basis_definition = [(z_cut1, [pow1, pow2, ...]), (z_cut2, [pow3, pow4, ...]), ...]`

- `z_cutX`: The cutoff distance for a group of basis functions.
- `[powX, powY, ...]`: A list of integer powers $n$ applied to the base term for that $z_{cut}$.

### 2.2. Cost Function

The cost function `E` is a sum of performance cost `E_performance` and accuracy cost `E_accuracy`.

`E = E_performance + E_accuracy`

**Performance Cost (`E_performance`):**
`E_performance = K_NBASIS \cdot nBasis + K_NZCUT \cdot nZcut + K_NMAX \cdot nMax + K_NSUM \cdot nSum`
- `nBasis`: Total number of individual basis functions.
- `nZcut`: Number of unique `z_cut` values used.
- `nMax`: The highest power $n$ used across all $z_{cut}$ groups.
- `nSum`: The sum of all powers $n$ used.
- `K_*`: Coefficients to weigh the importance of each term.

**Accuracy Cost (`E_accuracy`):**
`E_accuracy = K_logRMSE \cdot \log_{10}(RMSE) + K_RMSE_2 \cdot \max(0, \log_{10}(RMSE/xRMSE_2))`
- `RMSE`: Root Mean Square Error of approximating the target Morse potentials with the current basis set.
- `xRMSE_2`: A threshold RMSE value (e.g., $10^{-2}$).
- `K_logRMSE`, `K_RMSE_2`: Coefficients to weigh the accuracy terms. The second term adds a significant penalty if RMSE exceeds `xRMSE_2`.

The goal is to **minimize** this total cost $E$.

## 3. Mutation Strategies

The optimizer explores the solution space by applying various mutation strategies to the current `basis_definition`. Each mutation has a defined probability of being chosen.

1.  **`mut_add_pow`**: Adds a new power `n` (from the global `N_POWS` list) to a randomly selected existing `z_cut` group, if the power is not already present.
2.  **`mut_rem_pow`**: Removes a randomly selected power `n` from a randomly selected `z_cut` group. If this leaves the `z_cut` group empty, it's removed by `cleanup_bd`.
3.  **`mut_add_zc`**: Adds a new `z_cut` group. The $z_{cut}$ value is chosen from the global `Z_CUTS` list, and it's initialized with a single random power from `N_POWS`.
4.  **`mut_rem_zc`**: Removes an entire randomly selected `z_cut` group.
5.  **`mut_shift_zc`**:
    - Selects a $z_{cut}$ group and shifts its $z_{cut}$ value by $\pm D_{ZCUT}$.
    - The new $z_{cut}$ value is rounded (e.g., to 2 decimal places) and can be continuous (not restricted to `Z_CUTS`).
    - This is useful for fine-tuning $z_{cut}$ values.
6.  **`mut_drop_ld`**:
    - Uses squared cosine similarity matrix to detect linear dependencies between basis functions.
    - For each basis function $f_i$, computes:
      - $LD_{sum} = \sum_{j \neq i} \cos^2(f_i,f_j)$ (total linear dependence)
      - $LD_{max} = \max_{j \neq i} \cos^2(f_i,f_j)$ (strongest pairwise collinearity)
    - If the highest $LD_{max}$ exceeds `LD_DROP_TRASHOLD` (currently 0.7), drops one function from the most collinear pair using a deterministic heuristic:
      1. Prefers higher polynomial power (more complex basis function)
      2. If powers equal, prefers function with fewer powers at same $z_{cut}$ (to encourage shared $z_{cut}$)
      3. If still equal, prefers higher $LD_{sum}$ (more redundant function)
    - Also used in pre-evaluation checks with stricter threshold `LD_DROP_ALWAYS` (0.99) to reject degenerate bases early.

After each mutation, `cleanup_bd` is called to:
- Sort $z_{cut}$ groups by their $z_{cut}$ value.
- Remove any `z_cut` groups that have no powers.
- Ensure powers within each $z_{cut}$ group are sorted and unique.
- Merge $z_{cut}$ groups if their $z_{cut}$ values become identical (e.g., after rounding in `mut_shift_zc`).

## 4. Important Considerations and Insights

### 4.1. Two-Phase Optimization

A common strategy for problems like this involves two phases:
1.  **Coarse Global Optimization:**
    - `mut_shift_zc` might be disabled (probability set to 0) or `D_ZCUT` set to a large value corresponding to steps in the `Z_CUTS` list.
    - Focus is on finding good combinations of `z_cut` groups from the discrete `Z_CUTS` list and the general structure of powers.
    - One might consider using a hash of the basis definition to quickly check if a candidate solution has been evaluated before (not currently implemented, but useful for discrete spaces).
2.  **Fine Local Optimization:**
    - `mut_shift_zc` is enabled with a small `D_ZCUT` (e.g., 0.05-0.25) to allow fine-tuning of `z_cut` values, potentially moving them off the initial `Z_CUTS` grid.
    - Probabilities for other mutations might be adjusted to focus more on refinement.

The `D_ZCUT` parameter and the probabilities in `PROBS` should be configured according to the desired optimization phase.

### 4.2. Linear Dependence Thresholds

Two key thresholds control basis pruning:
1. **`LD_DROP_TRASHOLD` (0.7)**: During mutation, probabilistically drops basis functions when $LD_{max}$ exceeds this value.
2. **`LD_DROP_ALWAYS` (0.99)**: During pre-evaluation, strictly rejects any basis set where $LD_{max}$ exceeds this threshold.

These thresholds balance:
- **Orthogonality**: Lower thresholds produce more orthogonal bases but may be too aggressive.
- **Flexibility**: Higher thresholds allow more similar basis functions when needed for accuracy.

### 4.3. Parameter Pools (`Z_CUTS`, `N_POWS`) and `z0basis`

- `Z_CUTS`: Provides a discrete set of `z_cut` values primarily for the `mut_add_zc` strategy. `mut_shift_zc` can explore values outside this list.
- `N_POWS`: Provides a discrete set of powers `n` for `mut_add_pow` and initializing new `z_cut` groups.

### 4.4. Pre-evaluation Range Checking (`check_ranges`)

- The `check_ranges` function quickly validates a candidate basis definition against constraints defined in `RANGES` (e.g., min/max `nZcut`, min/max `nBasis`).
- This is a "cheap" check performed *before* the more computationally expensive fitness evaluation.
- If a candidate fails this check, it's immediately rejected with infinite cost, saving computation.

### 4.5. `cleanup_bd` Function

- This utility is crucial for maintaining a canonical and valid representation of the `basis_definition` after mutations.
- It handles sorting, removal of empty groups, and merging of identical `z_cut` groups (especially important when `mut_shift_zc` is active and `z_cut` values are rounded).

### 4.6. Weighted Fitting and Masking

- To focus the fitting process on physically relevant regions of the potential energy curves (typically the attractive well and not the highly repulsive core), a weighting scheme is employed.
- In `Optimize_basis_new.py`, this is handled by `WEIGHTS_Y_MASK_GLOBAL`. This mask is generated based on a threshold `V_REPULSIVE_THRESH_GLOBAL`.
- For each sample function in `Y_SAMPLES_GLOBAL`, points where the function's value is *below* `V_REPULSIVE_THRESH_GLOBAL` are given a weight of 1.0, and points at or above this threshold are given a weight of 0.0 (effectively ignoring them).
  ```python
  WEIGHTS_Y_MASK_GLOBAL = (Y_SAMPLES_GLOBAL < V_REPULSIVE_THRESH_GLOBAL).astype(float)
  ```
- This `WEIGHTS_Y_MASK_GLOBAL` is then passed to the `fit_coefficients` function (from `optimize.py`), which performs a weighted least-squares fit. This ensures that the RMSE calculation in `evaluate_fitness_basis` primarily reflects the error in the regions of interest.

### 4.7. `z0basis` Parameter

- The `z0basis` parameter (e.g., `1.0`) is used in `construct_composite_cutoff_basis` and `cutoff_poly_basis`. It defines the origin for the term $(z-z_0)$ in the basis function construction $(1-(z-z_0)/z_{span})^n$.
- Consistent use of `z0basis` is important for the definition and behavior of the basis functions.

## 5. Future Improvements and Directions

- **Adaptive Mutation Probabilities:** Probabilities for different mutation strategies could be adapted during the optimization based on their success rate.
- **More Sophisticated Cooling Schedule:** Explore different cooling schedules for simulated annealing.
- **Population-Based Approaches:** Consider algorithms like Genetic Algorithms or Particle Swarm Optimization, which maintain a population of solutions. This might require a mechanism to check for and handle duplicate solutions in the population, possibly using a hash for discrete basis definitions.
- **Gradient-Based Local Search:** For the fine-tuning phase, if the cost function (or parts of it) were differentiable with respect to `z_cut` values or even power coefficients (if powers were continuous), local gradient-based optimization could be integrated.
- **Multi-Objective Optimization:** Explicitly treat performance and accuracy as separate objectives rather than combining them into a single scalar cost function.
- **Database/Caching:** For coarse optimization phases with discrete `Z_CUTS` and `N_POWS`, caching results for previously evaluated basis definitions could speed up the process.

This documentation provides a snapshot of the current optimizer. Understanding these components is key to effectively using and further developing the tool.
