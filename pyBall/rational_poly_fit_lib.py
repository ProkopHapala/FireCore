"""
rational_poly_fit.py

Tools for fitting 1D functions f(x) with the model:
    f(x) = v(x) * ( P(x)^n / Q(x)^m )
using Chebyshev polynomials for numerical stability and converting 
to standard power basis for efficient runtime evaluation.
"""

import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial import chebyshev, polynomial
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

def parse_range(spec, cast=int):
    """
    Parse 'min,max,step' into list of casted values.
    Accepts 'a,b' (step=1) or 'a,b,step'.
    """
    parts = [p.strip() for p in str(spec).split(',') if p.strip() != ""]
    if len(parts) == 1:
        lo = hi = cast(parts[0]); step = 1
    elif len(parts) == 2:
        lo, hi = cast(parts[0]), cast(parts[1]); step = 1
    else:
        lo, hi, step = cast(parts[0]), cast(parts[1]), cast(parts[2])
    vals, v = [], lo
    while v <= hi:
        vals.append(cast(v))
        v += step
    return vals

def generate_configs(args):
    """
    Build list of configs using argparse-style args namespace with fields:
    p_deg_min/max/step, q_deg_min/max/step, p_term_min/step, q_term_min/step,
    n_range, m_range, basis, scan_level.
    """
    configs = []
    p_deg_max_list = range(args.p_deg_min, args.p_deg_max + 1, args.p_deg_step)
    q_deg_max_list = range(args.q_deg_min, args.q_deg_max + 1, args.q_deg_step)

    def build_term_list(max_deg, term_min, term_step):
        return list(range(term_min, max_deg + 1, term_step))

    n_list = parse_range(args.n_range, int)
    m_list = parse_range(args.m_range, int)

    # 1. Polynomial Powers P(x)^n (Q=1)
    if args.scan_level >= 1:
        for n_pow in n_list:
            for dp in p_deg_max_list:
                p_degs = build_term_list(dp, args.p_term_min, args.p_term_step)
                configs.append({
                    'p_degrees': p_degs,
                    'q_degrees': [],
                    'n': n_pow,
                    'm': 0,
                    'basis': args.basis
                })

    # 2. Rational P/Q
    if args.scan_level >= 2:
        for n_pow in n_list:
            for m_pow in m_list:
                for dp in p_deg_max_list:
                    p_degs = build_term_list(dp, args.p_term_min, args.p_term_step)
                    for dq in q_deg_max_list:
                        q_degs = build_term_list(dq, args.q_term_min, args.q_term_step)
                        configs.append({
                            'p_degrees': p_degs,
                            'q_degrees': q_degs,
                            'n': n_pow,
                            'm': m_pow,
                            'basis': args.basis
                        })

    return configs

def plot_pareto_and_fits(pareto, all_res, x, r, v, envelope_label="v(x)"):
    """
    Utility to plot pareto front, fits, and residuals for a set of pareto results.
    Expects entries from RationalPowerFitter.pareto_scan.
    """
    import matplotlib.pyplot as plt
    colors = plt.cm.viridis(np.linspace(0, 1, len(pareto)))

    plt.figure(figsize=(12, 6))

    # Subplot 1: Pareto Front (color-coded by model)
    costs = [r['n_coefs_total'] for r in all_res]
    errs = [r['error'] for r in all_res]
    plt.subplot(1, 3, 1)
    plt.scatter(costs, np.log10(errs), c='gray', alpha=0.4, label='All models')
    for i, model in enumerate(pareto):
        plt.plot(model['n_coefs_total'], np.log10(model['error']), 'o', color=colors[i], label=None)
        label = RationalPowerFitter.compact_pq_label(model)
        plt.annotate(label, (model['n_coefs_total'], np.log10(model['error'])),
                     textcoords="offset points", xytext=(5,5), ha='left', fontsize=7, color=colors[i])
    plt.xlabel('ncoefs (P+Q)')
    plt.ylabel('Log10 SSE Error')
    plt.title('Accuracy vs Efficiency')
    plt.grid(True)

    # Subplot 2: Fits of Pareto models
    plt.subplot(1, 3, 2)
    plt.plot(x, r, 'k:', label='Target', lw=1.5)
    for i, model in enumerate(pareto):
        P_poly = np.polynomial.Polynomial(model['poly_p'])
        Q_poly = np.polynomial.Polynomial(model['poly_q'])
        n_pow = model['config']['n']
        m_pow = model['config']['m']
        y_fit = v * ( (np.abs(P_poly(x)))**n_pow / (np.abs(Q_poly(x)) + 1e-10)**m_pow ) * np.sign(P_poly(x))
        lbl = f"{RationalPowerFitter.compact_pq_label(model)} err={model['error']:.2e}"
        plt.plot(x, y_fit, '-', color=colors[i], label=lbl, lw=0.8)
    plt.axhline(y=0, color='gray', linestyle='--', lw=0.5)
    plt.title('Pareto Fits')
    plt.legend(fontsize=7)

    # Subplot 3: Residuals for Pareto models
    plt.subplot(1, 3, 3)
    for i, model in enumerate(pareto):
        P_poly = np.polynomial.Polynomial(model['poly_p'])
        Q_poly = np.polynomial.Polynomial(model['poly_q'])
        n_pow = model['config']['n']
        m_pow = model['config']['m']
        y_fit = v * ( (np.abs(P_poly(x)))**n_pow / (np.abs(Q_poly(x)) + 1e-10)**m_pow ) * np.sign(P_poly(x))
        plt.plot(x, r - y_fit, '-', color=colors[i], label=None, lw=0.8)
    plt.title('Residuals (Target - Fit)')
    plt.grid(True)

    plt.tight_layout()
    return plt

class RationalPowerFitter:
    def __init__(self, x, y, v=None, w=None):
        """
        x: 1D numpy array of input positions
        y: 1D numpy array of reference values r(x)
        v: 1D numpy array of envelope v(x). Default 1.0
        w: 1D numpy array of weights. Default 1.0
        """
        self.x_orig = np.array(x, dtype=np.float64)
        self.y_orig = np.array(y, dtype=np.float64)
        
        self.v = v if v is not None else np.ones_like(self.x_orig)
        self.w = w if w is not None else np.ones_like(self.x_orig)
        
        # --- Domain Mapping ---
        # Map x from [min, max] to [-1, 1] for Chebyshev stability
        self.x_min = self.x_orig.min()
        self.x_max = self.x_orig.max()
        self.scale = 2.0 / (self.x_max - self.x_min)
        self.offset = -(self.x_max + self.x_min) / (self.x_max - self.x_min)
        
        # u is the variable in [-1, 1]
        self.u = self.x_orig * self.scale + self.offset

    def fit_model(self, p_degrees, q_degrees, n=1, m=1, prune_tol=0.0, basis="cheb"):
        """
        Fits the model.
        
        p_degrees: List of integer powers allowed in P (e.g. [0, 2, 4])
        q_degrees: List of integer powers allowed in Q (e.g. [1, 3])
                   Note: T0 (constant) for Q is always fixed to 1.0 internally.
        n, m:      Powers of the polynomials P and Q.
        prune_tol: If > 0, coefs smaller than this in Chebyshev basis are forced to 0.
        
        Returns: Dictionary with cost, coeffs, and function handles.
        """
        # 1. Setup Active Indices
        # For Q, we remove 0 from indices if present, because we fix Q_c0 = 1.0
        q_degrees_active = [d for d in q_degrees if d != 0]
        
        # 2. Initial Guess
        # Heuristic: Approximate P ~ (y/v)^(1/n). 
        # Calculate mean magnitude to set scale of P's constant term.
        mask = (self.v != 0)
        ratio = np.zeros_like(self.y_orig)
        ratio[mask] = self.y_orig[mask] / self.v[mask]
        
        # Avoid NaNs in root
        avg_ratio = np.mean(np.abs(ratio))
        p_scale = avg_ratio**(1.0/n) if n != 0 else 1.0
        
        # Params vector: [P_coefs..., Q_coefs...]
        # Initialize P constant term to p_scale, others 0
        n_p = len(p_degrees)
        n_q = len(q_degrees_active)
        
        p0 = np.zeros(n_p + n_q)
        if 0 in p_degrees:
            idx_0 = p_degrees.index(0)
            p0[idx_0] = p_scale
            
        # 3. Residual Function
        use_cheb = (basis == "cheb")

        def get_model(params):
            p_coefs_sparse = params[:n_p]
            q_coefs_sparse = params[n_p:]
            
            # Reconstruct full coefficient arrays
            max_p = max(p_degrees) if p_degrees else 0
            max_q = max(q_degrees_active) if q_degrees_active else 0
            
            c_p = np.zeros(max_p + 1)
            c_p[p_degrees] = p_coefs_sparse
            
            c_q = np.zeros(max_q + 1)
            c_q[0] = 1.0 # Fixed scale
            if q_degrees_active:
                c_q[q_degrees_active] = q_coefs_sparse
            
            if use_cheb:
                P_val = chebyshev.chebval(self.u, c_p)
                Q_val = chebyshev.chebval(self.u, c_q)
            else: # basis == "poly"
                P_val = polynomial.polyval(self.x_orig, c_p)
                Q_val = polynomial.polyval(self.x_orig, c_q)
            
            epsilon = 1e-10
            
            # Safe power evaluation (handling signs)
            # P^n
            P_term = np.sign(P_val) * (np.abs(P_val) + epsilon)**n
            # Q^m
            Q_term = (np.abs(Q_val) + epsilon)**m
            
            return self.v * (P_term / Q_term)

        def residuals(params):
            y_est = get_model(params)
            # Weighted squared error
            return np.sqrt(self.w) * (y_est - self.y_orig)

        # 4. Optimize
        res = least_squares(residuals, p0, method='lm', ftol=1e-7, xtol=1e-7)
        
        # 5. Pruning (Optional)
        best_params = res.x
        if prune_tol > 0:
            best_params[np.abs(best_params) < prune_tol] = 0.0
            # Optional: Refit with zero-constraints could be done here, 
            # but usually just zeroing is enough for high-degree pruning.
        
        # 6. Reconstruct Full Coefficients (Chebyshev or Poly)
        p_sparse = best_params[:n_p]
        q_sparse = best_params[n_p:]
        
        if use_cheb:
            c_p_cheb = np.zeros(max(p_degrees) + 1 if p_degrees else 1)
            c_p_cheb[p_degrees] = p_sparse
            
            c_q_cheb = np.zeros(max(q_degrees_active) + 1 if q_degrees_active else 1)
            c_q_cheb[0] = 1.0
            if q_degrees_active:
                c_q_cheb[q_degrees_active] = q_sparse
            
            # 7. Convert to Standard Polynomials on Original Domain
            poly_p_std = self._cheb_to_orig_poly(c_p_cheb)
            poly_q_std = self._cheb_to_orig_poly(c_q_cheb)
            params_cheb_p = c_p_cheb
            params_cheb_q = c_q_cheb
        else:
            c_p_poly = np.zeros(max(p_degrees) + 1 if p_degrees else 1)
            c_p_poly[p_degrees] = p_sparse
            
            c_q_poly = np.zeros(max(q_degrees_active) + 1 if q_degrees_active else 1)
            c_q_poly[0] = 1.0
            if q_degrees_active:
                c_q_poly[q_degrees_active] = q_sparse
            
            poly_p_std = c_p_poly
            poly_q_std = c_q_poly
            params_cheb_p = None
            params_cheb_q = None

        # 8. Calculate Metrics
        final_err = np.sum(residuals(best_params)**2)
        
        # Cost Model
        # Memory: Number of active (non-zero) coefficients
        mem_cost = np.count_nonzero(p_sparse) + np.count_nonzero(q_sparse)
        
        # Ops: Horner steps + Powers + Div
        # P(x) ops: deg_P adds + deg_P mults
        # P^n ops:  ~2 mults (if n small int)
        ops_cost = (len(poly_p_std)-1)*2 + (len(poly_q_std)-1)*2 + 2 # Rough estimate
        if m > 0: ops_cost += 4 # Division cost is higher
        n_coefs_total = len(p_degrees) + len(q_degrees)  # total configured coefficients (including fixed q0 term)
        
        return {
            "error": final_err,
            "params_cheb_p": params_cheb_p,
            "params_cheb_q": params_cheb_q,
            "poly_p": poly_p_std, # [a0, a1, a2...] standard basis
            "poly_q": poly_q_std,
            "mem_cost": mem_cost,
            "ops_cost": ops_cost,
            "n_coefs_total": n_coefs_total,
            "config": {"n": n, "m": m, "p_deg": p_degrees, "q_deg": q_degrees, "basis": basis}
        }

    @staticmethod
    def _poly_to_str(coefs, var="x", precision=3, fixed_width=False):
        """
        Render polynomial coefficients [c0, c1, ...] to string c0 + c1*x + c2*x^2 ...
        """
        if coefs is None or len(coefs) == 0:
            return "0"
        terms = []
        fmt = (lambda c: f"{c:+11.{precision}e}") if fixed_width else (lambda c: f"{c:.{precision}g}")
        for i, c in enumerate(coefs):
            coeff = fmt(c)
            if i == 0:
                terms.append(f"{coeff}")
            elif i == 1:
                terms.append(f"{coeff}*{var}")
            else:
                terms.append(f"{coeff}*{var}^{i}")
        return " + ".join(terms) if terms else "0"

    @staticmethod
    def format_rational_expression(res, envelope="v(x)", precision=3, fixed_width=False):
        """
        Build human-readable expression v(x)*(P(x))^n/(Q(x))^m from fit result dict.
        """
        cfg = res["config"]
        n = cfg.get("n", 1)
        m = cfg.get("m", 0)
        p_str = RationalPowerFitter._poly_to_str(res.get("poly_p", []), precision=precision, fixed_width=fixed_width)
        q_coefs = res.get("poly_q", [])
        if len(q_coefs) == 0:
            q_str = "1"
        else:
            q_str = RationalPowerFitter._poly_to_str(q_coefs, precision=precision, fixed_width=fixed_width)
        expr = f"{envelope}*({p_str})^{n}"
        if m > 0:
            expr += f"/({q_str})^{m}"
        return expr

    def _cheb_to_orig_poly(self, c_cheb):
        """
        Converts Chebyshev coeffs on [-1, 1] to Monomial coeffs on [x_min, x_max].
        Returns array [c0, c1, c2...] for c0 + c1*x + c2*x^2...
        """
        # 1. Cheb [-1,1] -> Poly [-1,1]
        c_std_u = chebyshev.cheb2poly(c_cheb)
        
        # 2. Map variable u -> x
        # u = scale * x + offset
        # We need P(scale*x + offset)
        # Use Polynomial library for composition
        P_u = polynomial.Polynomial(c_std_u)
        # Define line L(x) = offset + scale * x
        L_x = polynomial.Polynomial([self.offset, self.scale])
        
        # Compose
        P_x = P_u(L_x)
        return P_x.coef

    def pareto_scan(self, configs, cost_weights=(1.0, 0.0), verbose=False, workers=1, backend="thread", chunksize=10):
        """
        Brute-force scan over provided configurations to find Pareto front.
        configs: List of dicts like {'p_degrees':[], 'q_degrees':[], 'n':1, 'm':1}
        cost_weights: (w_mem, w_ops) to combine costs for sorting (optional)
        verbose: print progress for large scans
        workers: if >1, parallelize. backend controls executor: 'thread' or 'process'
        chunksize: only for process backend; batches tasks to reduce IPC overhead
        """
        results = []
        n_cfg = len(configs)
        backend = backend.lower()
        use_process = backend == "process"

        if workers <= 1:
            for i, cfg in enumerate(configs):
                try:
                    res = self.fit_model(
                        cfg['p_degrees'], 
                        cfg['q_degrees'], 
                        cfg['n'], 
                        cfg['m'],
                        basis=cfg.get('basis', 'cheb')
                    )
                    results.append(res)
                    if verbose:
                        label = RationalPowerFitter.compact_pq_label(res)
                        print(f"[{i+1}/{n_cfg}] {label} : err={res['error']:.3e} ncoefs={res['n_coefs_total']}")
                except Exception as e:
                    print(f"Fit failed for config {cfg}: {e}")
        else:
            if use_process:
                # ProcessPool for true parallelism
                with ProcessPoolExecutor(max_workers=workers) as ex:
                    # Submit all tasks
                    futures = [ex.submit(self.fit_model, 
                                         cfg['p_degrees'], cfg['q_degrees'], 
                                         cfg['n'], cfg['m'], basis=cfg.get('basis', 'cheb')) 
                               for cfg in configs]
                    
                    for i, fut in enumerate(futures):
                        try:
                            res = fut.result()
                            results.append(res)
                            if verbose:
                                label = RationalPowerFitter.compact_pq_label(res)
                                print(f"[{i+1}/{n_cfg}] {label} : err={res['error']:.3e} ncoefs={res['n_coefs_total']}")
                        except Exception as e:
                            print(f"Fit failed for config {configs[i]}: {e}")
            else:
                with ThreadPoolExecutor(max_workers=workers) as ex:
                    futures = [ex.submit(self.fit_model, 
                                         cfg['p_degrees'], cfg['q_degrees'], 
                                         cfg['n'], cfg['m'], basis=cfg.get('basis', 'cheb')) 
                               for cfg in configs]
                    for i, fut in enumerate(futures):
                        try:
                            res = fut.result()
                            results.append(res)
                            if verbose:
                                label = RationalPowerFitter.compact_pq_label(res)
                                print(f"[{i+1}/{n_cfg}] {label} : err={res['error']:.3e} ncoefs={res['n_coefs_total']}")
                        except Exception as e:
                            print(f"Fit failed for config {configs[i]}: {e}")

        # Sort by Error first
        results.sort(key=lambda x: x['error'])
        
        # Filter Pareto Front (Error vs Memory Cost)
        # A solution is on the front if no other solution has (Lower Error AND Lower Cost)
        pareto = []
        best_cost_so_far = float('inf')
        
        for res in results:
            cost = res['mem_cost']
            # If this model is cheaper than the cheapest model seen so far (which had lower error),
            # then it is a useful trade-off.
            if cost < best_cost_so_far:
                pareto.append(res)
                best_cost_so_far = cost
                
        # Order Pareto by cheapest (n_coefs) then error for consistent display
        pareto = sorted(pareto, key=lambda r: (r['n_coefs_total'], r['error']))
        return pareto, results

    @staticmethod
    def expression_label(res, precision=2):
        """
        Build short label 'ncoefs : error : v(x)*(p0+p1*x+...)^n/(q0+q1*x+...)^m'
        """
        expr = RationalPowerFitter.format_rational_expression(res, envelope="v(x)", precision=precision)
        return f"{res.get('n_coefs_total', 0)}: {res.get('error', 0):.2e}: {expr}"

    @staticmethod
    def compact_pq_label(res):
        """
        Compact label using counts of non-zero coefficients: P(p,n)/Q(q,m)
        p,q are counts of non-zero coefficients before applying powers n,m.
        """
        cfg = res.get("config", {})
        n = cfg.get("n", 1)
        m = cfg.get("m", 0)
        p_cnt = int(np.count_nonzero(res.get("poly_p", [])))
        q_cnt = int(np.count_nonzero(res.get("poly_q", [])))
        if m > 0:
            return f"P({p_cnt},{n})/Q({q_cnt},{m})"
        return f"P({p_cnt},{n})"
