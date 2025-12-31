"""
fit_cli.py

Command Line Interface for Rational Power Fitting.
Usage:
    python fit_cli.py --func "exp(-x**2)" --scan_level 2
"""

import sys
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../../")

from pyBall.rational_poly_fit_lib import (
    RationalPowerFitter,
    parse_range,
    generate_configs,
    plot_pareto_and_fits,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit rational power models to 1D data.")
    parser.add_argument("--func",       type=str, default="np.exp(-x**2)", help="Python expression for target function r(x). Use 'x' variable.")
    #parser.add_argument("--func",       type=str, default="np.exp(-x)", help="Python expression for target function r(x). Use 'x' variable.")
    parser.add_argument("--envelope",   type=str, default="(3.0-x)**2", help="Python expression for envelope v(x).")
    parser.add_argument("--range",      type=str, default="0,3", help="x_min,x_max")
    parser.add_argument("--points",     type=int, default=200, help="Number of points")
    parser.add_argument("--scan_level", type=int, default=1, help="1=PolyPowers, 2=Rational, 3=Sparse")
    # Degree controls
    parser.add_argument("--p_deg_min",   type=int, default=3)
    parser.add_argument("--p_deg_max",   type=int, default=7)
    parser.add_argument("--p_deg_step",  type=int, default=1)
    parser.add_argument("--q_deg_min",   type=int, default=1)
    parser.add_argument("--q_deg_max",   type=int, default=4)
    parser.add_argument("--q_deg_step",  type=int, default=1)
    # Term spacing (sparse patterns)
    parser.add_argument("--p_term_min",  type=int, default=0)
    parser.add_argument("--p_term_step", type=int, default=1)
    parser.add_argument("--q_term_min",  type=int, default=0)
    parser.add_argument("--q_term_step", type=int, default=1)
    # Power ranges (min,max,step)
    parser.add_argument("--n_range", type=str, default="1,4,1", help="range for power n: min,max,step (e.g., '0,8,2' for even)")
    parser.add_argument("--m_range", type=str, default="0,4,1", help="range for power m: min,max,step")
    # Basis choice
    parser.add_argument("--basis", type=str, default="cheb", choices=["cheb", "poly"], help="basis for fitting (cheb=Chebyshev, poly=monomial)")
    parser.add_argument("--output",     type=str, default=None, help="File to save coefs")
    parser.add_argument("--verbose",    type=int, default=1, help="Print progress during scan")

    args = parser.parse_args()

    # 1. Prepare Data
    rmin, rmax = map(float, args.range.split(','))
    x = np.linspace(rmin, rmax, args.points)
    
    # Safe eval context
    context = {"x": x, "np": np}
    try:
        r = eval(args.func, context)
        v = eval(args.envelope, context)
        if isinstance(v, float): v = np.ones_like(x)*v
    except Exception as e:
        print(f"Error evaluating function: {e}")
        exit(1)

    # 2. Setup Fitter
    # Simple weighting: focus more on small values? or uniform?
    # Let's use uniform for now.
    fitter = RationalPowerFitter(x, r, v=v)

    # 3. Run Scan
    print("Generating configurations...")
    configs = generate_configs(args)
    print(f"Testing {len(configs)} models...")
    
    pareto, all_res = fitter.pareto_scan(configs, verbose=args.verbose)
    
    print(f"\nFound {len(pareto)} Pareto-optimal models.")
    
    # 4. Display Results
    print(f"{'ncoefs':<8} |  {'Model':<14} |  {'Error (SSE)':<12} | {'Expression'}")
    print("-" * 110)
    
    for res in pareto:
        expr = RationalPowerFitter.format_rational_expression(res, envelope="v(x)", precision=3, fixed_width=True)
        model_label = RationalPowerFitter.compact_pq_label(res)
        print(f"{res['n_coefs_total']:<8} | {model_label:<14} | {res['error']:.4e}   | {expr}")

    # 5. Plotting (delegated to library helper)
    plot_pareto_and_fits(pareto, all_res, x, r, v)
    plt.show()

    # 6. Export Output
    if args.output:
        print(f"\nCoefficients for selected model saved to {args.output}")
        # Simple text dump
        with open(args.output, 'w') as f:
            f.write(f"# Config: {chosen['config']}\n")
            f.write("P_coeffs_std = " + str(list(chosen['poly_p'])) + "\n")
            f.write("Q_coeffs_std = " + str(list(chosen['poly_q'])) + "\n")