#!/usr/bin/env python3

import argparse
import json
import math
import os
import re
from typing import Dict, Optional

import numpy as np

kB = 8.617333262145e-5  # eV/K


def _read_data(filename: str) -> np.ndarray:
    data = np.loadtxt(filename, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _detect_columns(filename: str) -> Dict[str, int]:
    col_map = {
        "lambda": 0,
        "dE_dlambda": 1,
        "sigma_dE": 2,
        "cumulative_FE": 3,
        "cumulative_err": 4,
        "distance": 5,
    }
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not (line.startswith("#") and "lambda" in line):
                continue
            parts = line.replace("#", "").split()
            temp_map = {}
            for i, p in enumerate(parts):
                if p == "lambda":
                    temp_map["lambda"] = i
                elif p in ("dE/dlambda", "TI_dE/dlambda"):
                    temp_map["dE_dlambda"] = i
                elif p in ("sigma_dE", "TI_sigma"):
                    temp_map["sigma_dE"] = i
                elif p in ("cumulative_FE", "TI_F"):
                    temp_map["cumulative_FE"] = i
                elif p in ("cumulative_err", "TI_err"):
                    temp_map["cumulative_err"] = i
                elif p == "distance":
                    temp_map["distance"] = i
                elif p == "JE_F":
                    temp_map["JE_F"] = i
                elif p == "JE_F_sigma":
                    temp_map["JE_F_sigma"] = i
                elif p == "JE_W_avg":
                    temp_map["JE_W_avg"] = i
                elif p == "JE_W_sigma":
                    temp_map["JE_W_sigma"] = i
                elif p == "JE_W_skew":
                    temp_map["JE_W_skew"] = i
            if "lambda" in temp_map:
                col_map.update(temp_map)
            break
    return col_map


def _safe_col(data: np.ndarray, col_map: Dict[str, int], name: str) -> Optional[np.ndarray]:
    idx = col_map.get(name)
    if idx is None:
        return None
    if idx < 0 or idx >= data.shape[1]:
        return None
    return data[:, idx]


def _shift_to_zero(curve: Optional[np.ndarray]) -> Optional[np.ndarray]:
    if curve is None:
        return None
    if np.all(np.isnan(curve)):
        return curve
    valid = np.isfinite(curve)
    if not np.any(valid):
        return curve
    return curve - curve[valid][0]


def _r2_score(y: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0.0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def compute_metrics_from_file(
    input_file: str,
    mode: str,
    n_atoms: int,
    temperature: float = 300.0,
    b_segment: float = 1.198,
    lambda_fit_max: float = 0.05,
) -> Dict[str, float]:
    mode = mode.upper()
    data = _read_data(input_file)
    col_map = _detect_columns(input_file)

    # Legacy files may have header column indices that don't match actual data width.
    # Keep behavior aligned with plot_F_interactive.py: if distance index is invalid, use last column.
    if "distance" in col_map and col_map["distance"] >= data.shape[1]:
        col_map["distance"] = data.shape[1] - 1

    lam = _safe_col(data, col_map, "lambda")
    dist = _safe_col(data, col_map, "distance")
    if lam is None:
        raise ValueError(f"Cannot read lambda from {input_file}")
    if dist is None:
        raise ValueError(f"Cannot read distance from {input_file}")

    n_segments = n_atoms - 1
    if n_segments <= 0:
        raise ValueError(f"Invalid n_atoms={n_atoms}; need >=2")

    k_spring = kB * temperature / (n_segments * (b_segment ** 2))
    fe_ref_abs = 0.5 * k_spring * dist * dist
    fe_ref = fe_ref_abs - fe_ref_abs[0]

    ti_curve = _shift_to_zero(_safe_col(data, col_map, "cumulative_FE"))
    je_curve = _shift_to_zero(_safe_col(data, col_map, "JE_F"))

    if mode == "TI":
        curve = ti_curve
    elif mode == "JE":
        curve = je_curve
    else:
        raise ValueError(f"Unsupported mode={mode}, expected TI or JE")

    out: Dict[str, float] = {
        "k_spring_eV_A2": float(k_spring),
        "rmsd_profile_eV": float("nan"),
        "final_abs_error_eV": float("nan"),
        "final_value_eV": float("nan"),
        "final_ref_eV": float("nan"),
        "je_wiggle_slope_eV_per_lambda": float("nan"),
        "je_wiggle_intercept_eV": float("nan"),
        "je_wiggle_r2": float("nan"),
        "je_wiggle_rmse_eV": float("nan"),
    }

    if curve is not None:
        valid = np.isfinite(curve) & np.isfinite(fe_ref)
        if np.any(valid):
            diff = curve[valid] - fe_ref[valid]
            out["rmsd_profile_eV"] = float(np.sqrt(np.mean(diff * diff)))
            idx_last = np.where(valid)[0][-1]
            final_val = float(curve[idx_last])
            final_ref = float(fe_ref[idx_last])
            out["final_value_eV"] = final_val
            out["final_ref_eV"] = final_ref
            out["final_abs_error_eV"] = abs(final_val - final_ref)

    # JE wiggle fit (always from JE curve when available)
    if je_curve is not None:
        mask = (lam >= 0.0) & (lam <= lambda_fit_max) & np.isfinite(je_curve)
        x = lam[mask]
        y = je_curve[mask]
        if len(x) >= 2:
            m, c = np.polyfit(x, y, 1)
            y_pred = m * x + c
            out["je_wiggle_slope_eV_per_lambda"] = float(m)
            out["je_wiggle_intercept_eV"] = float(c)
            out["je_wiggle_rmse_eV"] = float(np.sqrt(np.mean((y - y_pred) ** 2)))
            out["je_wiggle_r2"] = float(_r2_score(y, y_pred))

    return out


def _guess_n_from_filename(path: str) -> Optional[int]:
    m = re.search(r"entropic_spring_(\d+)", os.path.basename(path))
    if m:
        return int(m.group(1))
    return None


def main() -> None:
    ap = argparse.ArgumentParser(description="Analyze FE profile accuracy against entropic spring reference.")
    ap.add_argument("--input", required=True, help="Input *_free_energy.dat file")
    ap.add_argument("--mode", required=True, choices=["TI", "JE"], help="Method to score")
    ap.add_argument("--N", type=int, default=None, help="Number of atoms in entropic spring (10/20/30)")
    ap.add_argument("--T", type=float, default=300.0, help="Temperature [K]")
    ap.add_argument("--b", type=float, default=1.198, help="Segment length [A]")
    ap.add_argument("--lambda-fit-max", type=float, default=0.05, help="Max lambda for JE wiggle linear fit")
    ap.add_argument("--output-json", default=None, help="Optional output JSON path")
    args = ap.parse_args()

    n_atoms = args.N if args.N is not None else _guess_n_from_filename(args.input)
    if n_atoms is None:
        raise ValueError("N could not be guessed from filename. Provide --N explicitly.")

    metrics = compute_metrics_from_file(
        input_file=args.input,
        mode=args.mode,
        n_atoms=n_atoms,
        temperature=args.T,
        b_segment=args.b,
        lambda_fit_max=args.lambda_fit_max,
    )
    payload = {
        "input": args.input,
        "mode": args.mode,
        "N": int(n_atoms),
        "temperature_K": args.T,
        "b_segment_A": args.b,
        "lambda_fit_max": args.lambda_fit_max,
        "metrics": metrics,
    }

    if args.output_json:
        with open(args.output_json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
