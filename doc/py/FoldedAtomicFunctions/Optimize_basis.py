from __future__ import annotations
import numpy as np
import random
import copy
import matplotlib.pyplot as plt
import functools

from Optimizer import Optimizer
from basis_utils import gen_morse_prms, gen_morse_curves, construct_composite_cutoff_basis, get_basis_func_map, remove_basis_func_by_flat_index
from optimize import fit_coefficients
import plot_utils as pu # For additional plotting

# --- Configuration for Basis Optimization ---

# This script aims to find optimal basis set for approximating variouse Morse potentials using basis functions of the form (z-z_cut)^n obtained from function cutoff_poly_basis() in basis_utils.py
# We aim to minimize the approximation error and performance cost (number of basis functions, number of unique z_cuts, power of polynominals )

# Cost function:
#  E             = E_performance + E_accuracy 
#  E_performance = K_NBASIS*nBasis     +  K_NZCUT*nZcut     + K_NMAX*nMax    +    K_NSUM*nSum  
#  E_accuracy    = K_logRMSE*log10(RMSE)  +  K_RMSE_2*max(0,log10(RMSE/xRMSE_2)) 

# Cost function parameter:

# Performance cost
K_NBASIS  = 1.0   # Cost per basis function (each basis function needs store coefficient)
K_NMAX    = 1.0   # Cost of highest power over all z_cuts (higher powers need more multiplications)
K_NSUM    = 1.0   # Cost of sum of highest powers over all z_cuts
K_NZCUT   = 3.0   # Cost per unique z_cut (it is cheaper to construct basis functions as powers of same z_cut, then introduce new z_cut)
xRMSE_2   = 1e-2
# Fit Accuracy 
K_logRMSE = 30.0   # Cost for the RMSE term
K_RMSE_2  = 100.0   # Cost if RMSE is above 1e-2 (special penalty if RMSE is above 1e-2)

LD_DROP_TRASHOLD = 0.8
# D_ZCUT: Increment for shifting z_cut in mut_shift_zc.
# For fine-tuning/local optimization, use a small D_ZCUT (e.g., 0.05-0.25).
# For coarse/global optimization with discrete z_cuts, D_ZCUT might be set
# to match the step in Z_CUTS, or mut_shift_zc might be disabled (prob=0).
D_ZCUT = 0.1


# --- Algorithm must respect following ranges, if the trial-solution is outside of these ranges, it is immediately rejected. Some ranges can be evaluated directly in mutation-phase (pre_eval=True), some need to evaluate full fitness function (pre_eval=False)
RANGES = [
# name      range      pre_eval? 
("nZcut",  (1,3),       True  ),   # this is checked before cost function evaluation
("nBasis", (1,8),       True  ),   # this is checked before cost function evaluation
#("RMSE",   (1e-1,1e-9), False ),  # this range must be evaluated after the cost function
]

# --- Mutation probabilities ---
#   here we list all possible mutations and their probabilities, 
#   we pre-calculate the cumulative probabilities and normalize them to sum to 1 for efficient sampling, and store it as simple numpy array, 
#   during sampling we use random number within (0,1) and use binary search to find the mutation event type, this is implemented in Optimizer.py
#   we may use list of callbacsk (functions or lamba function) which are called for each mutation event type, this list of functions is given to Optimizer
PROBS=[
("ADD_POW",     0.35), # Add power to existing z_cut
("REM_POW",     0.35), # Remove power from existing z_cut
("ADD_ZCUT",    0.15 ), # Add new z_cut
("REM_ZCUT",    0.15 ), # Remove z_cut
("DROP_LD",     0.20 ), # drop most linearly dependet basis function
("SHIFT_ZCUT",  0.10), # Adjust probability based on optimization phase
]

# --- In this specific use-case we have pre-defined discrete grid of possible z_cuts and powers, (i.e. we do discrete optimization), this allows us to creat efficient hash-table for searching if some trial-solution is already in the population
Z_CUTS=[4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,14.0,16.0]   # Available z_cuts (discrete values to choose from)
N_POWS=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]               # Available power factors (discrete values to choose from)

# --- Load or Generate Sample Data ---

# --- Load or Generate Sample Data ---
z0basis = 1.0 # z0 for basis construction
z_grid = np.linspace(z0basis, 12.0, 250)
n_samples = 20
np.random.seed(0)
prms_list = gen_morse_prms(n_samples, (1.0, 2.0), (2.0, 4.0), 1.0)
samples, _ = gen_morse_curves(z_grid, prms=prms_list)
Y_SAMPLES_GLOBAL = np.vstack(samples)
V_REPULSIVE_THRESH_GLOBAL = 0.5
WEIGHTS_Y_MASK_GLOBAL = (Y_SAMPLES_GLOBAL < V_REPULSIVE_THRESH_GLOBAL).astype(float)
TOTAL_MASKED_POINTS_GLOBAL = np.sum(WEIGHTS_Y_MASK_GLOBAL)

# --- Fitness Evaluation ---
def evaluate_fitness_basis(bd: list[tuple[float, list[int]]]) -> tuple[float, str]:
    if not bd:
        return float('inf'), "RMSE=inf, nBasis=0, nZcut=0, nMax=0, nSum=0"

    phi_c, _ = construct_composite_cutoff_basis(z_grid, bd, z0=z0basis)
    nb = phi_c.shape[0]

    if nb == 0:
        return float('inf'), f"RMSE=inf, nBasis=0, nZcut={len(bd)}, nMax=0, nSum=0"

    s_fit = fit_coefficients(Y_SAMPLES_GLOBAL, phi_c, weights=WEIGHTS_Y_MASK_GLOBAL)
    y_recon = s_fit.T @ phi_c
    
    rmse = float('inf')
    if TOTAL_MASKED_POINTS_GLOBAL > 0:
        sum_sq_err = np.sum(((Y_SAMPLES_GLOBAL - y_recon) * WEIGHTS_Y_MASK_GLOBAL)**2)
        rmse = np.sqrt(sum_sq_err / TOTAL_MASKED_POINTS_GLOBAL) + 1e-12 
    
    n_zc = len(bd)
    all_pows = [p for _, pows in bd for p in pows]
    n_max = max(all_pows) if all_pows else 0
    n_sum = sum(all_pows) if all_pows else 0

    cost = (K_logRMSE * np.log10(rmse) +
            K_RMSE_2 * max(0, np.log10(rmse / xRMSE_2)) +
            K_NBASIS * nb + K_NZCUT * n_zc + K_NMAX * n_max + K_NSUM * n_sum)

    details = f"RMSE={rmse:.2e}, nBasis={nb}, nZcut={n_zc}, nMax={n_max}, nSum={n_sum}"
    # print(f"eval_fitness_basis: {details} | Cost: {cost:.3e} | Basis: {bd}") # Debug print
    return cost, details

# --- Mutation Callback Functions ---
def cleanup_bd(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    """Sorts by z_cut, removes z_cuts with no powers, and ensures unique z_cuts."""
    cl_bd = [item for item in bd if item[1]]
    if not cl_bd: return []
    cl_bd.sort(key=lambda x: x[0])
    uniq_bd_d = {zc_val: sorted(list(set(pows_list))) for zc_val, pows_list in cl_bd}
    return sorted([(zc, pls) for zc, pls in uniq_bd_d.items()], key=lambda x: x[0])

def mut_add_pow(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    if not new_bd: return cleanup_bd(new_bd)
    idx = random.randrange(len(new_bd))
    zc_val, ex_pows = new_bd[idx]
    poss_pows = [p for p in N_POWS if p not in ex_pows] # N_POWS is global
    if poss_pows:
        new_bd[idx] = (zc_val, sorted(ex_pows + [random.choice(poss_pows)]))
    return cleanup_bd(new_bd)

def mut_rem_pow(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    if not new_bd: return cleanup_bd(new_bd)
    elig_idx = [i for i, (_, pows) in enumerate(new_bd) if pows]
    if not elig_idx: return cleanup_bd(new_bd)
    idx = random.choice(elig_idx)
    zc_val, ex_pows = new_bd[idx]
    pow_rem = random.choice(ex_pows)
    upd_pows = [p for p in ex_pows if p != pow_rem]
    if not upd_pows: new_bd.pop(idx)
    else: new_bd[idx] = (zc_val, sorted(upd_pows))
    return cleanup_bd(new_bd)

def mut_add_zc(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    ex_zcs = {zc for zc, _ in new_bd}
    avail_zcs = [zc for zc in Z_CUTS if zc not in ex_zcs] # Z_CUTS is global
    if avail_zcs and N_POWS: # N_POWS is global
        new_bd.append((random.choice(avail_zcs), [random.choice(N_POWS)]))
    return cleanup_bd(new_bd)

def mut_rem_zc(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    if not new_bd: return cleanup_bd(new_bd)
    new_bd.pop(random.randrange(len(new_bd)))
    return cleanup_bd(new_bd)

def mut_shift_zc(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    if not new_bd: return cleanup_bd(new_bd)
    
    # print(f"mut_shift_zc: START - bd: {new_bd}") # Debug print
    idx = random.randrange(len(new_bd))
    zc_old, pows = new_bd[idx]

    shift = random.choice([-D_ZCUT, D_ZCUT])
    zc_candidate = zc_old + shift
    
    # During fine-tuning, z_cuts can be continuous.
    # During coarse optimization, one might want to snap zc_candidate to the Z_CUTS grid,
    # or disable this mutation / set D_ZCUT to match grid steps.
    # Here, we allow continuous shifts, rounded to a few decimal places.
    # cleanup_bd will merge z_cuts if they become numerically identical after rounding.

    min_allowed_zc = z0basis + 0.1 # Example minimum, adjust as needed
    zc_new = max(min_allowed_zc, round(zc_candidate, 2)) # Round to 2 decimal places to avoid tiny differences

    if not np.isclose(zc_new, zc_old):
        # print(f"mut_shift_zc: Shifting zc {zc_old:.2f} to {zc_new:.2f} for powers {pows}") # Debug print
        new_bd[idx] = (zc_new, pows)
    # else:
        # print(f"mut_shift_zc: No change for zc {zc_old:.2f} (candidate {zc_candidate:.2f} resulted in {zc_new:.2f})") # Debug print

    # print(f"mut_shift_zc: END - bd before cleanup: {new_bd}") # Debug print
    return cleanup_bd(new_bd)

def mut_drop_ld(bd: list[tuple[float, list[int]]]) -> list[tuple[float, list[int]]]:
    new_bd = copy.deepcopy(bd)
    if not new_bd: return cleanup_bd(new_bd)

    phi_c, _ = construct_composite_cutoff_basis(z_grid, new_bd, z0=z0basis)
    n_funcs = phi_c.shape[0]

    if n_funcs < 2: return cleanup_bd(new_bd) # Need at least 2 functions to compare

    # Calculate dot products and norms squared
    dot_products = phi_c @ phi_c.T  # Gram matrix: G_ij = <f_i|f_j>
    norms_sq = np.diag(dot_products) # |f_i|^2

    # Avoid division by zero for zero-norm functions (should ideally not happen with good basis defs)
    norms_sq_safe = np.where(norms_sq < 1e-12, 1e-12, norms_sq)
    
    # Squared cosine similarities: cos_sim_sq_ij = <f_i|f_j>^2 / (|f_i|^2 * |f_j|^2)
    # gram_normalized_sq_ij = G_ij^2 / (norms_sq_i * norms_sq_j)
    gram_normalized_sq = dot_products**2 / (norms_sq_safe[:, np.newaxis] * norms_sq_safe[np.newaxis, :])
    
    # LD_i = sum_j (cos_sim_sq_ij) for j != i
    np.fill_diagonal(gram_normalized_sq, 0) # Exclude self-similarity from sum
    ld_values = np.sum(gram_normalized_sq, axis=1)

    idx_to_drop = np.argmax(ld_values)
    if ld_values[idx_to_drop] > LD_DROP_TRASHOLD:
        func_map = get_basis_func_map(new_bd) # Map flat index to (izcut, ipow_idx)
        new_bd = remove_basis_func_by_flat_index(new_bd, idx_to_drop, func_map)
    return cleanup_bd(new_bd)

# --- Pre-evaluation Range Checker ---
def check_ranges(bd: list[tuple[float, list[int]]]) -> tuple[bool, str]:
    """Checks if basis_definition bd is within pre-defined RANGES."""
    # RANGES is global
    details_parts = []
    n_zc = len(bd)
    details_parts.append(f"nZcut={n_zc}")
    min_zc_cfg, max_zc_cfg = next(item[1] for item in RANGES if item[0] == "nZcut")
    if not (min_zc_cfg <= n_zc <= max_zc_cfg):
        return False, f"Rejected (pre-eval: {', '.join(details_parts)})"

    min_nb_cfg, max_nb_cfg = next(item[1] for item in RANGES if item[0] == "nBasis")
    nb = sum(len(pows) for _, pows in bd) if bd else 0
    details_parts.append(f"nBasis={nb}")
    if not (min_nb_cfg <= nb <= max_nb_cfg):
        # print(f"check_ranges: REJECTED nBasis {nb} not in ({min_nb_cfg}, {max_nb_cfg}) for {bd}") # Debug print
        return False, f"Rejected (pre-eval: {', '.join(details_parts)})"
    return True, f"Pre-eval OK ({', '.join(details_parts)})"

def parse_details_str(s: str) -> dict:
    """ Parses the details string back into a dictionary for plotting. """
    metrics = {"RMSE": float('nan'), "nBasis": float('nan'), "nZcut": float('nan'), 
               "nMax": float('nan'), "nSum": float('nan')}
    if "Rejected" in s or not s: return metrics
    for part in s.split(', '):
        kv = part.split('=')
        if len(kv) == 2:
            key, v_str = kv
            try: metrics[key] = float(v_str) if key == "RMSE" else int(v_str)
            except ValueError: pass
    return metrics

# --- Main Optimization ---
if __name__ == "__main__":
    init_guess = [
        (Z_CUTS[random.randrange(len(Z_CUTS))], random.sample(N_POWS, random.randint(1,2))),
        (Z_CUTS[random.randrange(len(Z_CUTS))], random.sample(N_POWS, random.randint(1,2)))
    ]

    #init_guess = [(5.6, [3, 10, 16]), (9.1, [2, 5, 16]), (13.0, [3, 15])]



    init_guess = cleanup_bd(init_guess)
    
    is_init_valid, _ = check_ranges(init_guess)
    if not is_init_valid:
        print("Warning: Initial guess invalid. Using simplest valid guess.")
        init_guess = [(Z_CUTS[0], [N_POWS[0]])] 

    mut_cb_list = [mut_add_pow, mut_rem_pow, mut_add_zc, mut_rem_zc, mut_drop_ld, mut_shift_zc]
    probs_raw = np.array([p[1] for p in PROBS])
    probs_cum = np.cumsum(probs_raw / np.sum(probs_raw))

    opt = Optimizer(
        initial_solution=init_guess, evaluate_fitness=evaluate_fitness_basis,
        mutation_callbacks=mut_cb_list, mutation_cumulative_probs=probs_cum,
        pre_eval_range_checker=check_ranges, max_iterations=2500, # Restored iterations
        temperature_initial=10.0, temperature_decay=0.997, verbose=True)

    best_bd, best_fit, hist = opt.run()

    print("\n--- Best Basis Found ---")
    for zc, pows in best_bd: print(f"z_cut = {zc:.2f}, powers_n = {pows}")
    print(f"Fitness (Cost) = {best_fit:.4e}")

    iters, fits, details_strs, temps = zip(*hist)
    parsed_details = [parse_details_str(s) for s in details_strs]
    rmses = [d.get("RMSE", float('nan')) for d in parsed_details]
    nbs   = [d.get("nBasis", float('nan')) for d in parsed_details]
    nzcs  = [d.get("nZcut", float('nan')) for d in parsed_details]

    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axs[0].plot(iters, fits, color='tab:red', linestyle='-')
    axs[0].set_ylabel('Fitness (Cost)', color='tab:red')
    axs[0].tick_params(axis='y', labelcolor='tab:red')
    axs[0].grid(True, linestyle=':')
    ax_temp = axs[0].twinx()
    ax_temp.plot(iters, temps, color='tab:blue', linestyle='--')
    ax_temp.set_ylabel('Temperature', color='tab:blue')
    ax_temp.tick_params(axis='y', labelcolor='tab:blue')
    ax_temp.set_yscale('log')

    axs[1].plot(iters, rmses, label='RMSE', color='tab:green')
    axs[1].set_ylabel('RMSE'); axs[1].set_yscale('log'); axs[1].grid(True, linestyle=':')

    axs[2].plot(iters, nbs, label='nBasis', color='tab:purple')
    axs[2].set_ylabel('nBasis'); axs[2].grid(True, linestyle=':')

    axs[3].plot(iters, nzcs, label='nZcut', color='tab:orange')
    axs[3].set_ylabel('nZcut'); axs[3].set_xlabel('Iteration'); axs[3].grid(True, linestyle=':')
    
    fig.tight_layout()
    axs[0].set_title('Optimization Progress')

    if best_bd:
        phi_best, labels_best = construct_composite_cutoff_basis(z_grid, best_bd, z0=z0basis)
        if phi_best.shape[0] > 0:
            pu.plot1D(z_grid, phi_best, title=f"Best Basis Set ({phi_best.shape[0]} functions)", labels=labels_best, ylims=(-0.1, 1.1))
        else:
            print("Best basis resulted in an empty set of functions.")

    if best_bd and Y_SAMPLES_GLOBAL.shape[0] > 0 and phi_best.shape[0] > 0:
        s_best = fit_coefficients(Y_SAMPLES_GLOBAL, phi_best, weights=WEIGHTS_Y_MASK_GLOBAL)
        rng_plot = np.random.default_rng(1)
        n_plot = min(5, Y_SAMPLES_GLOBAL.shape[0])
        idx_plot = rng_plot.choice(Y_SAMPLES_GLOBAL.shape[0], size=n_plot, replace=False)
        plot_pairs = []
        for idx in idx_plot:
            y_orig = Y_SAMPLES_GLOBAL[idx]
            y_recon = s_best[:, idx].T @ phi_best
            plot_pairs.append((y_orig, y_recon))
        pu.plotMultiFunctionApprox(z_grid, plot_pairs, bError=True, errMax=0.05, scMin=1.2, title=f"Sample Approximations with Best Basis ({phi_best.shape[0]} functions)")
    plt.show()