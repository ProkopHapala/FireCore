"""
wf_pareto_fit.py

Standalone CLI to run Pareto rational-polynomial fits on a single Fireball wavefunction.

Examples:
- Fit H 1s (shell 1) from ./Fdata and show plots:
    python wf_pareto_fit.py --fdata ../tests/pyFireball/Fdata --nz 1 --shell 1
- Fit C 2p with wider rational search and save figure:
    python wf_pareto_fit.py --fdata ../tests/pyFireball/Fdata --nz 6 --shell 2 \
        --scan_level 2 --p_deg_min 2 --p_deg_max 10 --q_deg_min 1 --q_deg_max 6 \
        --output pareto_C_p.png --no_show
- Use polynomial (monomial) basis and even n,m ranges:
    python wf_pareto_fit.py --fdata ../tests/pyFireball/Fdata --nz 8 --shell 1 \
        --basis poly --n_range 0,8,2 --m_range 0,6,2

Found 6 Pareto-optimal models.
ncoefs   |  Model          |  Error (SSE)  | Expression
--------------------------------------------------------------------------------------------------------------
9        |  P(5,1)/Q(4,2)  |    4.3531e-03 | v(x)*( +2.702e-03 +  -6.792e-03*x +  +8.460e-03*x^2 +  -2.462e-03*x^3 +  +2.287e-04*x^4)^1/( +9.817e-02 +  -1.576e-01*x +  +2.521e-01*x^2 +  -3.002e-02*x^3)^2
10       |  P(7,5)/Q(3,10) |    9.4293e-04 | v(x)*( +9.384e-04 +  -3.481e-03*x +  +7.854e-03*x^2 +  -9.194e-03*x^3 +  +7.050e-03*x^4 +  -1.308e-03*x^5 +  +1.180e-04*x^6)^5/( +3.465e-02 +  -6.466e-02*x +  +8.720e-02*x^2)^10
11       |  P(6,3)/Q(5,2)  |    6.7508e-04 | v(x)*( +2.031e-03 +  -5.849e-03*x +  +9.292e-03*x^2 +  +1.774e-03*x^3 +  -8.662e-04*x^4 +  +1.132e-04*x^5)^3/( +1.716e-04 +  -7.983e-04*x +  +2.286e-03*x^2 +  -3.283e-03*x^3 +  +3.466e-03*x^4)^2
12       |  P(5,1)/Q(7,2)  |    3.9232e-04 | v(x)*( +1.185e-05 +  -2.930e-05*x +  +4.596e-05*x^2 +  -3.526e-05*x^3 +  +4.234e-05*x^4)^1/( +6.368e-03 +  -8.085e-03*x +  +7.600e-03*x^2 +  +1.719e-02*x^3 +  -4.606e-03*x^4 +  +2.940e-03*x^5 +  -4.027e-04*x^6)^2
13       |  P(6,3)/Q(7,4)  |    1.5809e-04 | v(x)*( +1.611e-03 +  -4.839e-03*x +  +8.376e-03*x^2 +  -6.924e-03*x^3 +  +5.610e-03*x^4 +  -7.264e-04*x^5)^3/( +1.095e-02 +  -2.516e-02*x +  +3.685e-02*x^2 +  -1.354e-02*x^3 +  +1.744e-02*x^4 +  -1.722e-03*x^5 +  -1.351e-04*x^6)^4
14       |  P(7,1)/Q(7,2)  |    8.3725e-05 | v(x)*( +9.740e-05 +  -4.204e-04*x +  +9.049e-04*x^2 +  -1.042e-03*x^3 +  +7.189e-04*x^4 +  -1.603e-04*x^5 +  +1.178e-05*x^6)^1/( +1.833e-02 +  -4.195e-02*x +  +5.976e-02*x^2 +  -1.886e-02*x^3 +  +2.440e-02*x^4 +  -4.886e-03*x^5 +  +2.157e-04*x^6)^2
        
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("../../"))

from pyBall.rational_poly_fit_lib import (
    RationalPowerFitter,
    parse_range,
    generate_configs,
    plot_pareto_and_fits,
)
from pyBall.FireballOCL.FdataParser import FdataParser


def load_wf(fdata_dir: str, nz: int, shell: int):
    parser = FdataParser(fdata_dir)
    parser.parse_info()
    wf_files = parser.find_wf(nz)
    if not wf_files or shell < 1 or shell > len(wf_files):
        raise RuntimeError(f"Missing wf shell={shell} for nz={nz} (found {wf_files})")
    rec = parser.read_wf(wf_files[shell - 1])
    r = np.linspace(0.0, rec["rcmax"], rec["mesh"])
    psi = rec["data"]
    return r, psi, rec["rcmax"]


def main():
    ap = argparse.ArgumentParser(description="Pareto rational-poly fits for one Fireball wavefunction")
    ap.add_argument("--fdata",    type=str, default="Fdata", help="Fdata directory")
    ap.add_argument("--nz",       type=int, default=6, help="nuclear charge (element Z)")
    ap.add_argument("--shell",    type=int, default=2, help="wavefunction shell index (s=1, p=2, d=3, ...)")
    #ap.add_argument("--envelope", type=str, default="1", help="envelope expression v(r) in terms of r")
    #ap.add_argument("--envelope", type=str, default="(6.-r)", help="envelope expression v(r) in terms of r")
    ap.add_argument("--envelope", type=str, default="r*(6.-r)**2", help="envelope expression v(r) in terms of r")
    #ap.add_argument("--envelope", type=str, default="(1+r**2)*(6.-r)**2", help="envelope expression v(r) in terms of r")
    ap.add_argument("--output",   type=str, default="pareto_wf.png", help="output figure for Pareto plots")
    ap.add_argument("--no_show",  action="store_true", help="do not call plt.show()")

    # Scan controls (match rational_poly_fit_cli)
    ap.add_argument("--scan_level", type=int, default=2, help="1=PolyPowers, 2=Rational, 3=Sparse")

    ap.add_argument("--p_deg_min",  type=int, default=4)
    ap.add_argument("--p_deg_max",  type=int, default=6)
    ap.add_argument("--p_deg_step", type=int, default=1)

    ap.add_argument("--q_deg_min",  type=int, default=0)
    ap.add_argument("--q_deg_max",  type=int, default=6)
    ap.add_argument("--q_deg_step", type=int, default=1)

    ap.add_argument("--p_term_min", type=int, default=0)
    ap.add_argument("--p_term_step",type=int, default=1)
    ap.add_argument("--q_term_min", type=int, default=0)
    ap.add_argument("--q_term_step",type=int, default=1)

    ap.add_argument("--n_range",    type=str, default="1,16,2",    help="range for power n: min,max,step")
    ap.add_argument("--m_range",    type=str, default="0,16,2",    help="range for power m: min,max,step")
    ap.add_argument("--verbose",    type=int,   default=1,         help="print progress during scan")
    ap.add_argument("--err_max",    type=float, default=1e-2,      help="optional SSE threshold to keep Pareto models")

    ap.add_argument("--basis",      type=str,   default="cheb",    choices=["cheb", "poly"], help="basis for fitting")
    ap.add_argument("--backend",    type=str,   default="process", choices=["thread", "process"], help="parallel backend for fitting")
    ap.add_argument("--workers",    type=int,   default=15,         help="number of worker threads/processes for fitting (>=1)")
    ap.add_argument("--chunksize",  type=int,   default=10,        help="batch size for process backend")

    args = ap.parse_args()

    # Load wavefunction
    r, psi, rcmax = load_wf(args.fdata, args.nz, args.shell)

    # Envelope
    ctx = {"r": r, "np": np}
    try:
        v = eval(args.envelope, ctx)
        if isinstance(v, (int, float)):
            v = np.ones_like(r) * v
    except Exception as e:
        raise RuntimeError(f"Failed to eval envelope '{args.envelope}': {e}")

    fitter = RationalPowerFitter(r, psi, v=v)
    fitter_configs = generate_configs(args)
    pareto, all_res = fitter.pareto_scan(
        fitter_configs,
        verbose=bool(args.verbose),
        workers=args.workers,
        backend=args.backend,
        chunksize=args.chunksize,
    )
    if args.err_max is not None:
        pareto = [res for res in pareto if res.get("error", 0.0) <= args.err_max]

    if len(pareto) == 0:
        print("No Pareto models under the specified error threshold.")
        return

    plt_obj = plot_pareto_and_fits(pareto, all_res, r, psi, v)
    plt_obj.savefig(args.output, dpi=200)
    print(f"Saved Pareto plot to {args.output}")
    print(f"Reference wf cutoff rcmax = {rcmax}")
    print(f"\nFound {len(pareto)} Pareto-optimal models.")
    print(f"{'ncoefs':<8} |  {'Model':<14} |  {'Error (SSE)':<12} | {'Expression'}")
    print("-"*110)
    for res in pareto:
        ncoefs = res.get("n_coefs_total", 0)
        model = RationalPowerFitter.compact_pq_label(res)
        err = res.get("error", 0.0)
        expr = RationalPowerFitter.format_rational_expression(res, envelope="v(x)", precision=3, fixed_width=True)
        print(f"{ncoefs:<8} |  {model:<14} |  {err:>12.4e} | {expr}")
    if not args.no_show:
        plt_obj.show()


if __name__ == "__main__":
    main()




"""

## Example - Results for s-orbital of Carbon iZ=1, shell 1



# ==== Envelop v(x)=1
Found 4 Pareto-optimal models.
ncoefs   |  Model          |  Error (SSE)  | Expression
--------------------------------------------------------------------------------------------------------------
9        |  P(4,7)/Q(5,4)  |    8.4301e-04 | v(x)*( +7.167e-02 +  -7.277e-02*x +  +6.793e-02*x^2 +  -5.321e-03*x^3)^7/( +8.967e-03 +  -1.613e-02*x +  +1.905e-02*x^2 +  -9.792e-03*x^3 +  +4.175e-03*x^4)^4
10       |  P(6,7)/Q(4,14) |    3.5583e-04 | v(x)*( +1.763e-02 +  -3.828e-02*x +  +4.802e-02*x^2 +  -3.157e-02*x^3 +  +1.213e-02*x^4 +  -9.133e-04*x^5)^7/( +1.290e-01 +  -1.405e-01*x +  +8.932e-02*x^2 +  +1.570e-03*x^3)^14
11       |  P(5,3)/Q(6,2)  |    2.3936e-04 | v(x)*( +2.158e-02 +  -1.856e-02*x +  +1.528e-02*x^2 +  +3.799e-03*x^3 +  -7.588e-04*x^4)^3/( +2.585e-03 +  -3.419e-03*x +  +2.037e-03*x^2 +  +1.400e-03*x^3 +  -1.186e-03*x^4 +  +7.055e-04*x^5)^2
12       |  P(6,7)/Q(6,10) |    6.8547e-05 | v(x)*( +4.469e-02 +  -8.728e-02*x +  +9.582e-02*x^2 +  -5.891e-02*x^3 +  +2.385e-02*x^4 +  -2.570e-03*x^5)^7/( +1.090e-01 +  -1.497e-01*x +  +1.084e-01*x^2 +  -2.976e-02*x^3 +  +1.403e-02*x^4 +  -1.608e-03*x^5)^10
prokophapala@carbsisYoga:~/git/FireCore/tests/pyFireball$ 

# ==== Envelop v(x)=(6.-r)
Found 4 Pareto-optimal models.
ncoefs   |  Model           |  Error (SSE)  | Expression
--------------------------------------------------------------------------------------------------------------
9        |  P(4,5 )/Q(5,4)  |    5.8308e-04 | v(x)*( +2.132e-02 +  -2.265e-02*x +  +2.303e-02*x^2 +  +1.319e-03*x^3)^5/( +1.153e-02 +  -1.620e-02*x +  +1.677e-02*x^2 +  -4.205e-03*x^3 +  +3.183e-03*x^4)^4
10       |  P(4,11)/Q(6,6)  |    2.4790e-04 | v(x)*( +1.197e-01 +  -1.292e-01*x +  +1.154e-01*x^2 +  -1.094e-02*x^3)^11/( +2.570e-02 +  -5.204e-02*x +  +6.665e-02*x^2 +  -4.177e-02*x^3 +  +1.866e-02*x^4 +  -1.853e-03*x^5)^6
11       |  P(4,3 )/Q(7,2)  |    2.1999e-04 | v(x)*( +6.382e-02 +  -7.156e-02*x +  +7.091e-02*x^2 +  -8.638e-03*x^3)^3/( +3.223e-02 +  -5.838e-02*x +  +5.971e-02*x^2 +  -2.128e-02*x^3 +  +5.130e-03*x^4 +  +2.426e-03*x^5 +  -4.458e-04*x^6)^2
12       |  P(7,13)/Q(5,12) |    2.3450e-05 | v(x)*( +7.897e-03 +  -1.391e-02*x +  +1.682e-02*x^2 +  -9.750e-03*x^3 +  +6.006e-03*x^4 +  -7.682e-04*x^5 +  +4.623e-05*x^6)^13/( +5.920e-03 +  -1.139e-02*x +  +1.392e-02*x^2 +  -7.987e-03*x^3 +  +3.994e-03*x^4)^12

# ==== Envelop v(x)=(6.-r)**2
Found 5 Pareto-optimal models.
ncoefs   |  Model            |  Error (SSE)  | Expression
--------------------------------------------------------------------------------------------------------------
10       |  P(4,7 )/Q(6,4)   |    2.5425e-04 | v(x)*( +1.112e-01 +  -1.237e-01*x +  +1.147e-01*x^2 +  -1.190e-02*x^3)^7/( +4.742e-02 +  -9.760e-02*x +  +1.237e-01*x^2 +  -7.683e-02*x^3 +  +3.465e-02*x^4 +  -3.933e-03*x^5)^4
11       |  P(5,9 )/Q(6,16)  |    1.3735e-04 | v(x)*( +1.223e-02 +  -2.498e-02*x +  +2.834e-02*x^2 +  -1.641e-02*x^3 +  +5.397e-03*x^4)^9/( +1.024e-01 +  -1.203e-01*x +  +7.686e-02*x^2 +  -6.436e-03*x^3 +  +4.117e-03*x^4 +  -4.099e-04*x^5)^16
12       |  P(7,13)/Q(5,10)  |    4.5483e-05 | v(x)*( +1.737e-02 +  -2.596e-02*x +  +2.768e-02*x^2 +  -1.205e-02*x^3 +  +7.373e-03*x^4 +  -1.304e-03*x^5 +  +9.600e-05*x^6)^13/( +7.075e-03 +  -1.400e-02*x +  +1.714e-02*x^2 +  -9.859e-03*x^3 +  +4.248e-03*x^4)^10
13       |  P(6,3 )/Q(7,2)   |    1.8752e-05 | v(x)*( +9.076e-03 +  +9.913e-04*x +  -7.953e-03*x^2 +  +1.717e-02*x^3 +  -3.668e-03*x^4 +  +2.421e-04*x^5)^3/( +4.228e-03 +  +2.393e-05*x +  -9.948e-03*x^2 +  +2.032e-02*x^3 +  -1.437e-02*x^4 +  +7.437e-03*x^5 +  -8.925e-04*x^6)^2
14       |  P(7,3 )/Q(7,2)   |    1.6104e-05 | v(x)*( +9.170e-03 +  -1.481e-03*x +  -4.409e-03*x^2 +  +1.288e-02*x^3 +  -2.069e-03*x^4 +  -3.539e-05*x^5 +  +2.004e-05*x^6)^3/( +4.294e-03 +  -1.745e-03*x +  -6.880e-03*x^2 +  +1.702e-02*x^3 +  -1.283e-02*x^4 +  +6.597e-03*x^5 +  -7.732e-04*x^6)^2


# ==== Envelop v(x)=(1+r**2)*(6.-r)**2
Found 5 Pareto-optimal models.
ncoefs   |  Model           |  Error (SSE)  | Expression
--------------------------------------------------------------------------------------------------------------
10       |  P(4,5 )/Q(6,4)  |    2.9074e-04 | v(x)*( +4.469e-02 +  -5.583e-02*x +  +5.177e-02*x^2 +  -5.225e-03*x^3)^5/( +4.547e-02 +  -7.589e-02*x +  +8.486e-02*x^2 +  -3.957e-02*x^3 +  +2.139e-02*x^4 +  -2.545e-03*x^5)^4
11       |  P(5,3 )/Q(6,2)  |    1.1023e-04 | v(x)*( +2.671e-02 +  -4.126e-02*x +  +3.718e-02*x^2 +  -8.383e-03*x^3 +  +6.746e-04*x^4)^3/( +2.137e-02 +  -5.383e-02*x +  +7.192e-02*x^2 +  -4.920e-02*x^3 +  +2.048e-02*x^4 +  -1.959e-03*x^5)^2
12       |  P(7,13)/Q(5,10) |    4.1984e-05 | v(x)*( +1.758e-02 +  -2.541e-02*x +  +2.566e-02*x^2 +  -9.404e-03*x^3 +  +5.642e-03*x^4 +  -1.026e-03*x^5 +  +7.602e-05*x^6)^13/( +7.185e-03 +  -1.376e-02*x +  +1.669e-02*x^2 +  -9.395e-03*x^3 +  +4.174e-03*x^4)^10
13       |  P(6,15)/Q(7,14) |    2.3409e-05 | v(x)*( +1.131e-02 +  -1.482e-02*x +  +1.668e-02*x^2 +  -6.836e-03*x^3 +  +7.650e-03*x^4 +  -9.619e-04*x^5)^15/( +1.030e-02 +  -1.471e-02*x +  +1.679e-02*x^2 +  -6.633e-03*x^3 +  +5.316e-03*x^4 +  +1.166e-03*x^5 +  -2.711e-04*x^6)^14
14       |  P(7,5 )/Q(7,4)  |    1.9331e-05 | v(x)*( +6.602e-03 +  -4.298e-03*x +  +1.436e-03*x^2 +  +4.509e-03*x^3 +  +9.943e-04*x^4 +  -4.127e-04*x^5 +  +3.161e-05*x^6)^5/( +4.162e-03 +  -3.738e-03*x +  +9.358e-04*x^2 +  +5.659e-03*x^3 +  -4.318e-03*x^4 +  +3.743e-03*x^5 +  -4.762e-04*x^6)^4


"""