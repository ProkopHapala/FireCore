from __future__ import annotations
import numpy as np
import random
import copy
import matplotlib.pyplot as plt

from Optimizer import Optimizer
from basis_utils import gen_morse_prms, gen_morse_curves, construct_composite_cutoff_basis
from optimize import fit_coefficients

# --- Configuration for Basis Optimization ---

# Cost function weights (GLOBAL PARAMS K)
K_RMSE = 1.0       # Weight for the RMSE term
K_NBASIS = 0.01    # Cost per basis function
K_NZCUT = 0.1      # Cost per unique z_cut used
K_NMAX = 0.005     # Cost related to the highest power factor used
K_NSUM = 0.001     # Cost related to the sum of all power factors

# Mutation probabilities
PROB_ADD_POWER_TO_EXISTING_ZCUT = 0.30
PROB_REMOVE_POWER_FROM_EXISTING_ZCUT = 0.25
PROB_CHANGE_POWER_EXISTING_ZCUT = 0.15
PROB_ADD_NEW_ZCUT = 0.10
PROB_REMOVE_EMPTY_ZCUT = 0.05 # Or any zcut if not empty
PROB_CHANGE_EXISTING_ZCUT_VAL = 0.10
PROB_CHANGE_DNS_ELEMENT = 0.05 # If dns is part of the state

# Constraints for basis parameters
MIN_ZCUT = 4.0
MAX_ZCUT = 16.0
ZCUT_STEP = 0.5
MIN_POWER_FACTOR = 1
MAX_POWER_FACTOR = 10 # Max 'n' in (zc-z)^(2n)
MAX_POWERS_PER_ZCUT = 6
MAX_TOTAL_ZCUTS = 5

# --- Load or Generate Sample Data ---
z_grid = np.linspace(1.5, 12.0, 250)
n_samples = 20
np.random.seed(0)
prms_list = gen_morse_prms(n_samples, (1.0, 2.0), (2.0, 4.0), 1.0)
samples, _ = gen_morse_curves(z_grid, prms=prms_list)
Y_SAMPLES_GLOBAL = np.vstack(samples)
V_REPULSIVE_THRESH_GLOBAL = 0.5
WEIGHTS_Y_MASK_GLOBAL = (Y_SAMPLES_GLOBAL < V_REPULSIVE_THRESH_GLOBAL).astype(float)

# --- Fitness Evaluation ---
# To store detailed parameters for plotting
OPTIMIZATION_TRAJECTORY_DETAILS = []

def evaluate_fitness_basis(basis_definition: list[tuple[float, list[int]]]) -> tuple[float, str]:
    """
    Evaluates the fitness of a given basis definition. Lower is better.
    basis_definition is like [(z_cut1, [n1, n2]), (z_cut2, [n3, n4, n5]), ...]
    Returns (total_cost, details_string)
    """
    # Initialize details for this evaluation
    current_eval_details = {}

    if not basis_definition: # Handle empty basis definition
        OPTIMIZATION_TRAJECTORY_DETAILS.append({"RMSE": float('inf'), "nBasis": 0, "nZcut":0, "nMax":0, "nSum":0, "cost": float('inf')})
        return float('inf'), "RMSE=inf, nBasis=0"

    phi_composite, _ = construct_composite_cutoff_basis(z_grid, basis_definition)
    nbasis = phi_composite.shape[0]

    if nbasis == 0: # If construction results in no basis functions
        OPTIMIZATION_TRAJECTORY_DETAILS.append({"RMSE": float('inf'), "nBasis": 0, "nZcut":len(basis_definition), "nMax":0, "nSum":0, "cost": float('inf')})
        return float('inf'), f"RMSE=inf, nBasis=0"

    # 1. RMSE
    s_fit = fit_coefficients(Y_SAMPLES_GLOBAL, phi_composite, weights=WEIGHTS_Y_MASK_GLOBAL)
    y_reconstructed = s_fit.T @ phi_composite
    
    total_masked_points = np.sum(WEIGHTS_Y_MASK_GLOBAL)
    if total_masked_points > 0:
        sum_sq_error = np.sum(((Y_SAMPLES_GLOBAL - y_reconstructed) * WEIGHTS_Y_MASK_GLOBAL)**2)
        rmse = np.sqrt(sum_sq_error / total_masked_points)
    else:
        rmse = float('inf') # Should not happen if Y_SAMPLES_GLOBAL is valid

    # 2. nZcut
    n_zcut = len(basis_definition)

    # 3. nmax and nsum
    all_powers = []
    for _, powers in basis_definition:
        all_powers.extend(powers)
    
    n_max = max(all_powers) if all_powers else 0
    n_sum = sum(all_powers) if all_powers else 0

    cost = (K_RMSE * rmse +
            K_NBASIS * nbasis +
            K_NZCUT * n_zcut +
            K_NMAX * n_max +
            K_NSUM * n_sum)
    
    current_eval_details["RMSE"] = rmse
    current_eval_details["nBasis"] = nbasis
    current_eval_details["nZcut"] = n_zcut
    current_eval_details["nMax"] = n_max
    current_eval_details["nSum"] = n_sum
    current_eval_details["cost"] = cost
    OPTIMIZATION_TRAJECTORY_DETAILS.append(current_eval_details)

    details_str = f"RMSE={rmse:.2e}, nBasis={nbasis}, nZcut={n_zcut}, nMax={n_max}, nSum={n_sum}"
    # print(f" Eval: {details_str} -> Cost={cost:.3e}")

    return cost, details_str

# --- Mutation Proposals ---

def propose_mutation_basis(current_basis_def: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    """
    Proposes a mutation to the current basis definition.
    """
    new_basis_def = copy.deepcopy(current_basis_def)
    mutation_type = random.random()

    # Ensure there's at least one z_cut to work with for some mutations
    if not new_basis_def and mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT + PROB_REMOVE_EMPTY_ZCUT + PROB_CHANGE_EXISTING_ZCUT_VAL):
        mutation_type = PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT + PROB_REMOVE_EMPTY_ZCUT + PROB_CHANGE_EXISTING_ZCUT_VAL + 0.01 # Force add new z_cut

    if mutation_type < PROB_ADD_POWER_TO_EXISTING_ZCUT and new_basis_def:
        idx = random.randrange(len(new_basis_def))
        if len(new_basis_def[idx][1]) < MAX_POWERS_PER_ZCUT:
            existing_powers = set(new_basis_def[idx][1])
            possible_new_powers = [p for p in range(MIN_POWER_FACTOR, MAX_POWER_FACTOR + 1) if p not in existing_powers]
            if possible_new_powers:
                new_power = random.choice(possible_new_powers)
                new_basis_def[idx][1].append(new_power)
                new_basis_def[idx][1].sort()
    elif mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT) and new_basis_def:
        idx = random.randrange(len(new_basis_def))
        if new_basis_def[idx][1]:
            new_basis_def[idx][1].pop(random.randrange(len(new_basis_def[idx][1])))
            if not new_basis_def[idx][1] and random.random() < 0.5: # Chance to remove zcut if it becomes empty
                new_basis_def.pop(idx)
    elif mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT) and new_basis_def:
        idx = random.randrange(len(new_basis_def))
        if new_basis_def[idx][1]:
            power_idx_to_change = random.randrange(len(new_basis_def[idx][1]))
            new_power = random.randint(MIN_POWER_FACTOR, MAX_POWER_FACTOR)
            # Avoid duplicate if possible, or just change it
            if new_power not in new_basis_def[idx][1] or len(new_basis_def[idx][1]) == 1:
                 new_basis_def[idx][1][power_idx_to_change] = new_power
                 new_basis_def[idx][1].sort()
    elif mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT + PROB_ADD_NEW_ZCUT):
        if len(new_basis_def) < MAX_TOTAL_ZCUTS:
            new_zc = round(random.uniform(MIN_ZCUT, MAX_ZCUT) / ZCUT_STEP) * ZCUT_STEP
            num_initial_powers = random.randint(1, min(3, MAX_POWERS_PER_ZCUT))
            initial_powers = sorted(random.sample(range(MIN_POWER_FACTOR, MAX_POWER_FACTOR + 1), num_initial_powers))
            # Avoid adding z_cut if it already exists (or merge later)
            if not any(abs(zc_existing - new_zc) < 1e-3 for zc_existing, _ in new_basis_def):
                 new_basis_def.append((new_zc, initial_powers))
    elif mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT + PROB_ADD_NEW_ZCUT + PROB_REMOVE_EMPTY_ZCUT) and new_basis_def:
        idx_to_remove = random.randrange(len(new_basis_def))
        # Prefer removing empty or small ones
        if not new_basis_def[idx_to_remove][1] or random.random() < 0.7:
            new_basis_def.pop(idx_to_remove)
    elif mutation_type < (PROB_ADD_POWER_TO_EXISTING_ZCUT + PROB_REMOVE_POWER_FROM_EXISTING_ZCUT + PROB_CHANGE_POWER_EXISTING_ZCUT + PROB_ADD_NEW_ZCUT + PROB_REMOVE_EMPTY_ZCUT + PROB_CHANGE_EXISTING_ZCUT_VAL) and new_basis_def:
        idx = random.randrange(len(new_basis_def))
        old_zc, powers = new_basis_def.pop(idx)
        new_zc_val = round(random.uniform(MIN_ZCUT, MAX_ZCUT) / ZCUT_STEP) * ZCUT_STEP
        # Avoid collision or merge
        if not any(abs(zc_existing - new_zc_val) < 1e-3 for zc_existing, _ in new_basis_def):
            new_basis_def.append((new_zc_val, powers))
        else: # Put it back if collision and no merge logic
            new_basis_def.append((old_zc, powers))

    # Clean up: remove z_cuts with no powers if not already handled
    new_basis_def = [item for item in new_basis_def if item[1]]
    # Clean up: sort by z_cut and ensure unique z_cuts (simple approach: last one wins if duplicate z_cut)
    if new_basis_def:
        new_basis_def.sort(key=lambda x: x[0])
        unique_def = []
        seen_zcuts = set()
        for zc_val, powers_list in reversed(new_basis_def): # Keep last definition for a z_cut
            zc_rounded = round(zc_val / ZCUT_STEP) * ZCUT_STEP # Group close z_cuts
            if zc_rounded not in seen_zcuts:
                unique_def.append((zc_rounded, powers_list))
                seen_zcuts.add(zc_rounded)
        new_basis_def = sorted(unique_def, key=lambda x: x[0])

    return new_basis_def

# --- Main Optimization ---
if __name__ == "__main__":
    initial_basis_guess = [
        (8.0, [2, 3, 4]),
        (6.0, [1, 2])
    ]

    optimizer = Optimizer(
        initial_solution=initial_basis_guess,
        evaluate_fitness=evaluate_fitness_basis,
        propose_mutation=propose_mutation_basis,
        max_iterations=500, # Increase for better results
        temperature_initial=1.0, # Adjust based on cost function scale
        temperature_decay=0.99,
        verbose=True
    )

    best_basis, best_fitness_val, history_data = optimizer.run()

    print("\n--- Best Basis Found ---")
    for zc, powers in best_basis:
        print(f"z_cut = {zc:.2f}, powers_n = {powers}")
    print(f"Fitness (Cost) = {best_fitness_val:.4e}")

    # Plot optimization history
    iters_hist, fitness_hist, _, temp_hist = zip(*history_data) # _ is for details_string
    
    # Extract detailed parameters from OPTIMIZATION_TRAJECTORY_DETAILS
    # Note: OPTIMIZATION_TRAJECTORY_DETAILS might have more entries than history_data if some mutations are rejected
    # We'll plot based on the length of history_data, assuming one eval per iteration step recorded in history.
    num_recorded_iters = len(iters_hist)
    recorded_details = OPTIMIZATION_TRAJECTORY_DETAILS[:num_recorded_iters]

    rmse_traj = [d.get("RMSE", float('nan')) for d in recorded_details]
    nbasis_traj = [d.get("nBasis", float('nan')) for d in recorded_details]
    nzcut_traj = [d.get("nZcut", float('nan')) for d in recorded_details]

    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # Plot Fitness (Cost)
    color = 'tab:red'
    axs[0].set_ylabel('Fitness (Cost)', color=color)
    axs[0].plot(iters_hist, fitness_hist, color=color, linestyle='-')
    axs[0].tick_params(axis='y', labelcolor=color)
    axs[0].set_yscale('log')
    axs[0].grid(True, linestyle=':')

    ax_temp = axs[0].twinx()
    color = 'tab:blue'
    ax_temp.set_ylabel('Temperature', color=color)
    ax_temp.plot(iters_hist, temp_hist, color=color, linestyle='--')
    ax_temp.tick_params(axis='y', labelcolor=color)
    ax_temp.set_yscale('log')

    # Plot RMSE
    axs[1].plot(iters_hist, rmse_traj, label='RMSE', color='tab:green')
    axs[1].set_ylabel('RMSE')
    axs[1].set_yscale('log')
    axs[1].grid(True, linestyle=':')

    # Plot nBasis and nZcut
    axs[2].plot(iters_hist, nbasis_traj, label='nBasis', color='tab:purple')
    axs[2].set_ylabel('nBasis')
    axs[2].grid(True, linestyle=':')

    axs[3].plot(iters_hist, nzcut_traj, label='nZcut', color='tab:orange')
    axs[3].set_ylabel('nZcut')
    axs[3].set_xlabel('Iteration')
    axs[3].grid(True, linestyle=':')
    
    fig.tight_layout()
    axs[0].set_title('Optimization Progress')
    plt.show()