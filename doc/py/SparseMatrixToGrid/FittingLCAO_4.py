import numpy as np
import matplotlib.pyplot as plt
import time
import argparse

def sample_errors(log, sample_iters):
    out = []
    for it in sample_iters:
        idx = it - 1
        if idx < len(log):
            out.append(log[idx])
        else:
            out.append(log[-1])
    return out

def print_error_table(sample_iters, entries):
    header = ["iter"] + [name for name, _ in entries]
    col_w = 14
    print("\nError table (selected iterations):")
    print("".join(f"{h:>{col_w}}" for h in header))
    for it in sample_iters:
        row = [it]
        for _, log in entries:
            row.append(sample_errors(log, [it])[0] if len(log) else float("nan"))
        print("".join(f"{val:>{col_w}.4e}" if i > 0 else f"{val:>{col_w}d}" for i, val in enumerate(row)))

class GPULibrarySimulator:
    def __init__(self, K=200, N=40, noise_level=0.1, width_scale=1.2, seed=None):
        if seed is not None: np.random.seed(seed)
        self.K = K
        self.N = N
        self.x_grid = np.linspace(-5, 5, K)
        self.centers = np.linspace(-4, 4, N)
        self.width = (self.centers[1] - self.centers[0]) * width_scale
        
        # Operator A
        self.A = np.zeros((K, N))
        for j in range(N):
            self.A[:, j] = np.exp(-(self.x_grid - self.centers[j])**2 / (2*self.width**2))
            
        # Ground truth
        self.c_true = np.random.normal(0.0, 1.0, N)
        y_true = self.A @ self.c_true
        self.y_ref = y_true + np.random.normal(0, noise_level, K)
        self.b = self.backward(self.y_ref)

    def forward(self, c):
        return self.A @ c

    def backward(self, r):
        return self.A.T @ r

    def get_residual(self, c):
        y_pred = self.forward(c)
        r = y_pred - self.y_ref
        grad = self.backward(r)
        err = np.linalg.norm(r)
        return grad, err

    # =========================================================================
    # METHOD 1: Conjugate Gradient
    # =========================================================================
    def solve_cg(self, max_iter=100, tol=1e-9):
        c = np.zeros(self.N)
        r_vec = self.b - self.backward(self.forward(c))
        p_vec = r_vec.copy()
        rsold = np.dot(r_vec, r_vec)
        error_log = []
        
        start_time = time.time()
        for i in range(max_iter):
            Hp = self.backward(self.forward(p_vec))
            alpha = rsold / (np.dot(p_vec, Hp) + 1e-16)
            c = c + alpha * p_vec
            r_vec = r_vec - alpha * Hp
            rsnew = np.dot(r_vec, r_vec)
            
            _, real_err = self.get_residual(c)
            error_log.append(real_err)
            
            if np.sqrt(rsnew) < tol: break
            p_vec = r_vec + (rsnew / rsold) * p_vec
            rsold = rsnew
            
        dt = time.time() - start_time
        return c, error_log, dt

    # =========================================================================
    # METHOD 2: Chebyshev Semi-Iteration
    # FIXED: Added 'cond_est' parameter to handle ill-conditioned matrices
    # =========================================================================
    def power_iteration(self, iter=20):
        v = np.random.rand(self.N)
        v /= np.linalg.norm(v)
        lambda_max = 0
        for _ in range(iter):
            Hv = self.backward(self.forward(v))
            vn = np.linalg.norm(Hv)
            if vn < 1e-12: break
            v = Hv / vn
            lambda_max = np.dot(v, Hv)
        return lambda_max

    def solve_chebyshev(self, max_iter=100, cond_est=10000.0):
        c = np.zeros(self.N)
        
        # Tuning: Overestimating lambda_max is safer than underestimating
        eig_max = self.power_iteration(iter=20) * 1.1 
        eig_min = eig_max / cond_est 
        
        d = (eig_max + eig_min) / 2.0
        c_cen = (eig_max - eig_min) / 2.0
        c_prev = np.zeros_like(c)
        
        error_log = []
        start_time = time.time()
        
        for k in range(max_iter):
            grad = self.backward(self.forward(c)) - self.b
            
            if k == 0:
                alpha = 1.0 / d
                c_next = c - alpha * grad
            else:
                beta = (c_cen / 2.0) ** 2
                if k == 1:
                    alpha = 1.0 / (d - beta / d)
                    gamma = 2.0 * d / (2.0 * d * d - c_cen * c_cen)
                else:
                    gamma = 1.0 / (d - (c_cen * c_cen * gamma / 4.0))
                
                # Update
                term1 = c - grad / d
                c_next = gamma * term1 + (1.0 - gamma) * c_prev
                
            c_prev = c.copy()
            c = c_next.copy()
            _, real_err = self.get_residual(c)
            error_log.append(real_err)
            
        dt = time.time() - start_time
        return c, error_log, dt

    # =========================================================================
    # METHOD 3: Anderson Acceleration
    # =========================================================================
    def solve_anderson(self, max_iter=100, m=5, mix_beta=0.1):
        c = np.zeros(self.N)
        X_hist, F_hist = [], []
        error_log = []
        start_time = time.time()
        
        for k in range(max_iter):
            grad = self.backward(self.forward(c)) - self.b
            f_val = -mix_beta * grad 
            
            X_hist.append(c.copy())
            F_hist.append(f_val.copy())
            if len(X_hist) > m: X_hist.pop(0); F_hist.pop(0)
            
            mk = len(X_hist)
            if mk > 1:
                df = [F_hist[-1] - F_hist[i] for i in range(mk-1)]
                dx = [X_hist[-1] - X_hist[i] for i in range(mk-1)]
                DF = np.stack(df, axis=1)
                DX = np.stack(dx, axis=1)
                try:
                    gamma, _, _, _ = np.linalg.lstsq(DF, F_hist[-1], rcond=None)
                except: gamma = np.zeros(mk-1)
                c_new = (X_hist[-1] + F_hist[-1]) - (DX @ gamma + DF @ gamma)
            else:
                c_new = c + f_val

            c = c_new
            _, real_err = self.get_residual(c)
            error_log.append(real_err)
            
        dt = time.time() - start_time
        return c, error_log, dt

    # =========================================================================
    # METHOD 4: Jacobi + Momentum (The "Projective Dynamics" approach)
    # =========================================================================
    def solve_jacobi_momentum(self, max_iter=100, momentum=0.9, taper_last=5, step_scale=0.2):
        c = np.zeros(self.N)
        velocity = np.zeros(self.N)

        # Diagonal Preconditioner (Inverse Mass)
        H_diag = np.sum(self.A * self.A, axis=0)
        H_diag[H_diag < 1e-12] = 1.0
        inv_diag = 1.0 / H_diag

        error_log = []
        start_time = time.time()

        for k in range(max_iter):
            y_pred = self.forward(c)
            residual_sample = y_pred - self.y_ref
            grad = self.backward(residual_sample)

            # Jacobi Step: c* = c - step * D^-1 * grad
            c_star = c - step_scale * inv_diag * grad

            # Momentum Mix
            bmix = momentum
            if taper_last > 0 and k >= max_iter - taper_last: bmix = 0.0
            
            c_new = c_star + bmix * velocity
            velocity = c_new - c
            c = c_new

            _, real_err = self.get_residual(c)
            error_log.append(real_err)

        dt = time.time() - start_time
        return c, error_log, dt

    # =========================================================================
    # METHOD 5: Gradient Descent (Effective Damped Jacobi)
    # =========================================================================
    def solve_gd(self, max_iter=100, step_scale=0.1, precondition=True):
        c = np.zeros(self.N)
        if precondition:
            H_diag = np.sum(self.A * self.A, axis=0)
            H_diag[H_diag < 1e-12] = 1.0
            inv_diag = 1.0 / H_diag
        else:
            inv_diag = 1.0

        error_log = []
        start_time = time.time()

        for k in range(max_iter):
            grad = self.backward(self.forward(c) - self.y_ref)
            
            # If precondition=True, this IS Jacobi
            step = inv_diag * grad
            c = c - step_scale * step

            _, real_err = self.get_residual(c)
            error_log.append(real_err)

        dt = time.time() - start_time
        return c, error_log, dt

    # =========================================================================
    # METHOD 6: Barzilai-Borwein (New Idea)
    # Gradient Descent with adaptive step size based on secant equation.
    # =========================================================================
    def solve_bb(self, max_iter=100, precondition=True):
        c = np.zeros(self.N)
        
        # Preconditioner helps BB significantly too
        if precondition:
            H_diag = np.sum(self.A * self.A, axis=0)
            H_diag[H_diag < 1e-12] = 1.0
            inv_diag = 1.0 / H_diag
        else:
            inv_diag = 1.0

        # Initial Step
        y_pred = self.forward(c)
        grad = self.backward(y_pred - self.y_ref)
        if precondition: grad *= inv_diag
        
        c_prev = c.copy()
        g_prev = grad.copy()
        
        # First step is arbitrary small GD
        c = c - 1e-4 * grad
        
        error_log = []
        start_time = time.time()

        for k in range(max_iter):
            y_pred = self.forward(c)
            grad = self.backward(y_pred - self.y_ref)
            if precondition: grad *= inv_diag

            # Compute Step Size using BB1 rule
            s = c - c_prev
            y = grad - g_prev
            
            # Dot products (Global Reductions)
            sy = np.dot(s, y)
            ss = np.dot(s, s)
            
            if sy > 1e-16:
                alpha = ss / sy 
            else:
                alpha = 0.001 # Fallback
            
            # Safety clamp for alpha to prevent explosions in early steps
            alpha = min(max(alpha, 1e-5), 2.0)

            c_prev = c.copy()
            g_prev = grad.copy()
            
            c = c - alpha * grad
            
            _, real_err = self.get_residual(c)
            error_log.append(real_err)

        dt = time.time() - start_time
        return c, error_log, dt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--K",             type=int,     default=300)
    parser.add_argument("--N",             type=int,     default=60)
    parser.add_argument("--width_scale",   type=float,   default=1.2)
    parser.add_argument("--noise",         type=float,   default=0.0)
    parser.add_argument("--seed",          type=int,     default=42)
    parser.add_argument("--iters",         type=int,     default=128)
    parser.add_argument("--cheb_cond",     type=float,   default=1e5,  help="Condition number est for Chebyshev")
    parser.add_argument("--anderson_m",    type=int,     default=7)
    parser.add_argument("--anderson_beta", type=float,   default=0.1)
    parser.add_argument("--mom_beta",      type=float,   default=0.9)
    parser.add_argument("--mom_step",      type=float,   default=0.2)
    parser.add_argument("--mom_taper",     type=int,     default=5)
    parser.add_argument("--gd_step",       type=float,   default=0.2)
    parser.add_argument("--gd_precond",    type=int,     default=1)
    parser.add_argument("--plot",          type=int,     default=1)
    args = parser.parse_args()

    sample_iters = [1, 2, 4, 8, 16, 32, 64, 128]

    print(f"Initializing Problem (K={args.K}, N={args.N})...")
    sim = GPULibrarySimulator(K=args.K, N=args.N, noise_level=args.noise, width_scale=args.width_scale, seed=args.seed)
    
    print("\nRunning Solvers...")
    
    c_cg, log_cg, time_cg = sim.solve_cg(max_iter=args.iters)
    print(f"CG              : Final Err {log_cg[-1]:.4f} | Time {time_cg*1000:.2f} ms")
    
    c_cheb, log_cheb, time_cheb = sim.solve_chebyshev(max_iter=args.iters, cond_est=args.cheb_cond)
    print(f"Chebyshev       : Final Err {log_cheb[-1]:.4f} | Time {time_cheb*1000:.2f} ms (cond={args.cheb_cond})")
    
    c_and, log_and, time_and = sim.solve_anderson(max_iter=args.iters, m=args.anderson_m, mix_beta=args.anderson_beta)
    print(f"Anderson        : Final Err {log_and[-1]:.4f} | Time {time_and*1000:.2f} ms")

    c_jac, log_jac, time_jac = sim.solve_jacobi_momentum(max_iter=args.iters, momentum=args.mom_beta, taper_last=args.mom_taper, step_scale=args.mom_step)
    print(f"Jacobi+Momentum : Final Err {log_jac[-1]:.4f} | Time {time_jac*1000:.2f} ms")

    c_gd, log_gd, time_gd = sim.solve_gd(max_iter=args.iters, step_scale=args.gd_step, precondition=args.gd_precond)
    print(f"Jacobi (GD+Pre) : Final Err {log_gd[-1]:.4f} | Time {time_gd*1000:.2f} ms")

    c_bb, log_bb, time_bb = sim.solve_bb(max_iter=args.iters, precondition=True)
    print(f"Barzilai-Borwein: Final Err {log_bb[-1]:.4f} | Time {time_bb*1000:.2f} ms")
    
    print_error_table(
        sample_iters,
        [
            ("CG", log_cg),
            ("Chebyshev", log_cheb),
            ("Anderson", log_and),
            ("Jacobi+Mom", log_jac),
            ("Jacobi(GD)", log_gd),
            ("BarzilaiBor", log_bb),
        ],
    )
    
    if args.plot:
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(log_cg,   label='CG', lw=1)
        plt.plot(log_cheb, label='Chebyshev', lw=1)
        plt.plot(log_and,  label='Anderson', lw=1)
        plt.plot(log_jac,  label='Jacobi+Mom', lw=1)
        plt.plot(log_gd,   label='Jacobi (GD)', lw=1)
        plt.plot(log_bb,   label='Barzilai-Borwein', lw=1)
        plt.yscale('log')
        plt.xlabel('Iteration')
        plt.ylabel('L2 Residual')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        plt.plot(sim.x_grid, sim.y_ref, 'k.', alpha=0.3)
        plt.plot(sim.x_grid, sim.forward(c_cg), 'r-', label='Fit')
        plt.title('Fit Result')
        plt.legend()
        plt.tight_layout()
        plt.show()