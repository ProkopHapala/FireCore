"""Batch runner for Optimize_basis.

Runs the basis optimiser multiple times starting from random initial bases and
appends the best result of each run (if unique) to a JSON-Lines file.

Usage (defaults shown):
    python Optimize_basis_batch.py \
        --samples Morse1.json \
        --results results.jsonl \
        --runs 10 \
        --seed 0

Each output line is a JSON object of the form
    {"<z0>,<basis_str>": {"rmse": 1.23e-3,
                           "timestamp": "2025-06-18T12:00:00Z",
                           "samples": "Morse1.json"}}
meaning the key itself serves as a unique hash of z0basis and the canonical
basis-definition string.
"""
from __future__ import annotations
import argparse, json, random, pathlib, datetime, sys, os
from typing import List, Tuple
import numpy as np

# Project imports
#sys.path.append(str(pathlib.Path(__file__).parent))  # ensure local import works
#auto = __import__
#ob = auto("Optimize_basis")               # main optimisation module

import Optimize_basis as ob  # main optimisation module
from results_db import ResultsDB
from basis_utils import gen_morse_prms, gen_morse_curves, load_morse_samples_json  # utilities

###############################################################################
# Helper utilities
###############################################################################

def basis_to_string(bd: List[Tuple[float, List[int]]]) -> str:
    """Canonical string representation – sorted by z_cut then powers."""
    bd_sorted = ob.cleanup_bd(bd)
    return str(bd_sorted)

def random_initial_basis(rng: random.Random, max_groups: int = 2) -> List[Tuple[float, List[int]]]:
    n_groups = rng.randint(1, max_groups)
    zcuts = rng.sample(ob.Z_CUTS, k=n_groups)
    bd = []
    for zc in sorted(zcuts):
        n_pows = rng.randint(1, 2)
        powers = rng.sample(ob.N_POWS, k=n_pows)
        bd.append((float(zc), sorted(powers)))
    return ob.cleanup_bd(bd)

def load_samples_json(filename: str) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Load sample data from JSON file using shared utility."""
    z0, z_grid, Y_samples, weights_mask, _ = load_morse_samples_json(filename)
    return z0, z_grid, Y_samples, weights_mask

def append_jsonl(path: str, obj: dict):
    """Append one JSON object as single line (newline terminated)."""
    with open(path, "a", encoding="utf-8") as f:
        json.dump(obj, f, separators=(",", ":"))
        f.write("\n")

def load_existing_keys(path: str):
    keys = set()
    if pathlib.Path(path).is_file():
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                try:
                    d = json.loads(line)
                    keys.update(d.keys())
                except json.JSONDecodeError:
                    continue
    return keys

def gen_sample_json_morse(filename: str, params: dict):
    """Generate and save Morse potential samples to JSON file."""
    from basis_utils import gen_morse_prms, gen_morse_curves
    
    # Generate parameters in dictionary format
    a_vals  = np.random.uniform(params['a_rng'][0], params['a_rng'][1], params['n_samples'])
    r0_vals = np.random.uniform(params['r0_rng'][0], params['r0_rng'][1], params['n_samples'])
    prms = [{'a': a, 'r0': r0, 'D': params['D_val']} for a, r0 in zip(a_vals, r0_vals)]
    
    z_grid = np.linspace(params['z0basis'], params['z_max'], params['nz'])
    curves, _ = gen_morse_curves(z_grid, prms=prms)
    
    data = {
        "z0basis": params['z0basis'],
        "z_grid": {"start": params['z0basis'], "stop": params['z_max'], "num": params['nz']},
        "v_repulsive_thresh": params['V_REPULSIVE_THRESH_GLOBAL'],
        "morse_params": prms
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)

###############################################################################
# Main batch logic
###############################################################################

def run_optimizer(init_bd, max_iter=500, verbose=False):
    # Optimiser parameters mirror defaults in Optimize_basis
    mut_cb_list = [ob.mut_add_pow, ob.mut_rem_pow, ob.mut_add_zc, ob.mut_rem_zc, ob.mut_drop_ld, ob.mut_shift_zc]
    probs_raw = np.array([p[1] for p in ob.PROBS])
    probs_cum = np.cumsum(probs_raw / np.sum(probs_raw))

    opt = ob.Optimizer(
        initial_solution=init_bd,
        evaluate_fitness=ob.evaluate_fitness_basis,
        mutation_callbacks=mut_cb_list,
        mutation_cumulative_probs=probs_cum,
        pre_eval_checker=ob.pre_eval_check,
        max_iterations=max_iter,
        temperature_initial=10.0,
        temperature_decay=0.995,
        verbose=verbose,
    )
    best_bd, best_cost, _ = opt.run()

    # derive rmse from best_cost via parse_details_str -> evaluate again
    phi_best, _ = ob.construct_composite_cutoff_basis(ob.z_grid, best_bd, z0=ob.z0basis)
    s_best = ob.fit_coefficients(ob.Y_SAMPLES_GLOBAL, phi_best,
                                 weights=ob.WEIGHTS_Y_MASK_GLOBAL)
    y_recon = s_best.T @ phi_best
    rmse = ob.weighted_rmse(ob.Y_SAMPLES_GLOBAL, y_recon, weights=ob.WEIGHTS_Y_MASK_GLOBAL)
    return best_bd, rmse

def plot_rmse_vs_nbasis(db_items: Dict[str, dict], *, show=True, ax=None):
    db = ResultsDB(args.results)
    all_dict = dict(db.items())
    # Prepare scatter data
    def xy_from_items(items):
        xs, ys = [], []
        for k, rec in items.items():
            rmse = rec.get("rmse")
            if rmse is None: continue
            _z0, bd = basis_key_to_struct(k)
            nb, *_ = get_basis_metrics(bd)
            xs.append(nb); ys.append(rmse)
        return xs, ys
    x_all, y_all = xy_from_items(all_dict)
    x_new, y_new = xy_from_items({k: all_dict[k] for k in new_keys}) if new_keys else ([], [])

    fig, ax = plt.subplots(figsize=(6,4))
    ax.scatter(x_all, y_all, c='k', marker='.', s=20, label='All runs')
    if x_new:
        ax.scatter(x_new, y_new, c='r', marker='o', s=40, label='This batch')
    ax.set_xlabel('nBasis'); ax.set_ylabel('RMSE'); ax.set_yscale('log'); ax.grid(True, linestyle=':')
    ax.legend()
    plt.tight_layout()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch optimiser runner")
    parser.add_argument("--nsamp", type=int, default=10)
    parser.add_argument("--samples", default="Morse1.json")
    parser.add_argument("--results", default="results.jsonl")
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max_iter", type=int, default=500)
    parser.add_argument("--no_plot", action="store_true", help="Skip plotting at end")
    args = parser.parse_args()

    # Generate samples if file doesn't exist
    if not os.path.exists(args.samples):
        print(f"Sample file {args.samples} not found - generating new samples")
        morse_params = {
            'z0basis': ob.z0basis,
            'V_REPULSIVE_THRESH_GLOBAL': ob.V_REPULSIVE_THRESH_GLOBAL,
            'z_max': ob.z_max,
            'nz': ob.nz,
            'n_samples': args.nsamp,
            'a_rng': ob.a_rng,
            'r0_rng': ob.r0_rng,
            'D_val': ob.D_val
        }
        gen_sample_json_morse(args.samples, params=morse_params)

    # Load sample data
    z0, z_grid, Y_samples, weights_mask = load_samples_json(args.samples)

    # Inject into optimize_basis globals so its evaluate_fitness works
    ob.z0basis = z0
    ob.z_grid = z_grid
    ob.Y_SAMPLES_GLOBAL = Y_samples
    ob.WEIGHTS_Y_MASK_GLOBAL = weights_mask

    # Warm-up RNGs
    rng_py = random.Random(args.seed)
    np.random.seed(args.seed)

    existing_keys = load_existing_keys(args.results)
    print(f"Loaded {len(existing_keys)} existing result keys from {args.results}")
    new_keys = []  # track new records this session

    for run_idx in range(args.runs):
        init_bd = random_initial_basis(rng_py)
        best_bd, rmse = run_optimizer(init_bd, max_iter=args.max_iter, verbose=False)
        key = f"{z0},{basis_to_string(best_bd)}"
        if key in existing_keys:
            print(f"[{run_idx+1}/{args.runs}] Duplicate – skipped: {key}")
            continue
        ts = datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z"
        record = {key: {"rmse": rmse, "timestamp": ts, "samples": args.samples}}
        append_jsonl(args.results, record)
        existing_keys.add(key)
        new_keys.append(key)
        print(f"[{run_idx+1}/{args.runs}] Added result: RMSE={rmse:.2e} | key={key}")

    print("Batch finished. Total unique results now:", len(existing_keys))

    # ------------------------------------------------------ Plot summary
    if not args.no_plot and existing_keys:
        import matplotlib.pyplot as plt
        from results_db import ResultsDB
        from plot_db_utils import basis_key_to_struct
        from basis_utils import get_basis_metrics
        plot_rmse_vs_nbasis(existing_keys, show=True)
        plt.show()
        #plt.savefig('rmse_vs_nbasis_batch.png')
        #print('Plot saved to rmse_vs_nbasis_batch.png')