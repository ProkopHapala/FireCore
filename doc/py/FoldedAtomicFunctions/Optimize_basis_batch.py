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
import argparse, json, random, pathlib, datetime, sys
from typing import List, Tuple

import numpy as np

# Project imports
sys.path.append(str(pathlib.Path(__file__).parent))  # ensure local import works
auto = __import__
ob = auto("Optimize_basis")               # main optimisation module
from basis_utils import gen_morse_curves, gen_morse_prms  # utilities

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

def load_samples_json(path: str):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    z0 = float(data.get("z0basis", 1.0))
    # z_grid can be full list or dict spec
    zg_spec = data["z_grid"]
    if isinstance(zg_spec, dict):
        z_grid = np.linspace(zg_spec["start"], zg_spec["stop"], zg_spec["num"])
    else:
        z_grid = np.array(zg_spec, dtype=float)
    morse_params_raw = data["morse_params"]
    prms_dicts = [{"D": p[0], "a": p[1], "r0": p[2]} for p in morse_params_raw]
    Y_samples, _ = gen_morse_curves(z_grid, prms=prms_dicts)
    Y_samples = np.vstack(Y_samples)
    # masking weights (reuse same heuristic as optimiser)
    v_thresh = float(data.get("v_repulsive_thresh", 0.5))
    weights_mask = (Y_samples < v_thresh).astype(float)
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

###############################################################################
# Main batch logic
###############################################################################

def run_optimizer(init_bd, max_iter=500, verbose=False):
    # Optimiser parameters mirror defaults in Optimize_basis
    mut_cb_list = [ob.mut_add_pow, ob.mut_rem_pow, ob.mut_add_zc,
                   ob.mut_rem_zc, ob.mut_drop_ld, ob.mut_shift_zc]
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

def main():
    parser = argparse.ArgumentParser(description="Batch optimiser runner")
    parser.add_argument("--samples", default="Morse1.json")
    parser.add_argument("--results", default="results.jsonl")
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max_iter", type=int, default=500)
    args = parser.parse_args()

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
        print(f"[{run_idx+1}/{args.runs}] Added result: RMSE={rmse:.2e} | key={key}")

    print("Batch finished. Total unique results now:", len(existing_keys))

if __name__ == "__main__":
    main()
