import numpy as np

from OrderN import _estimate_spectral_bounds, _scale_hamiltonian, _solve_linear


def chebyshev_coeffs_fermi(mu_tilde=0.0, beta=50.0, n_poly=50, n_grid=400):
    """
    Approximate Fermi-Dirac f(x)=1/(1+exp(beta*(x-mu))) on [-1,1] with Chebyshev.
    Returns coefficients c_k such that f(x) ~ sum c_k T_k(x).
    """
    xs = np.cos(np.pi * (np.arange(n_grid) + 0.5) / n_grid)  # Gauss-Cheb nodes
    fxs = 1.0 / (1.0 + np.exp(beta * (xs - mu_tilde)))
    coeffs = np.zeros(n_poly)
    for k in range(n_poly):
        coeffs[k] = (2.0 / n_grid) * np.sum(fxs * np.cos(k * np.arccos(xs)))
    coeffs[0] *= 0.5
    print(f"#DEBUG chebyshev_coeffs_fermi n_poly={n_poly} beta={beta} mu_tilde={mu_tilde}")
    return coeffs


def foe_stochastic_density(
    H,
    S=None,
    n_poly=60,
    n_random=32,
    beta=80.0,
    mu=None,
    jacobi_steps=8,
    solver="jacobi",
    orthogonalize=True,
    use_clenshaw=True,
    ref_rho=None,
    diag_every=0,
    diag_prefix="#INFO FOE_iter",
    calc_mulliken=True,
    beta_is_scaled=False,
):
    """
    Stochastic trace FOE to approximate diagonal of density matrix.
    - H, S: dense arrays
    - n_poly: Chebyshev order
    - n_random: number of Rademacher probe vectors
    - beta: inverse temperature for smoothing (higher=sharper step)
    - mu: chemical potential (defaults to mid-spectrum after scaling)
    - solver: 'jacobi' (default), 'solve', 'cholesky'
    - orthogonalize: if True and S provided, transform to orthonormal basis to stabilize
    Returns rho_est (N,) and auxiliary info dict.
    """
    N = H.shape[0]
    R = n_random

    # Optional orthogonalization to stabilize S handling
    X = None
    H_eff = H
    if S is not None:
        if orthogonalize:
            s, U = np.linalg.eigh(S)
            eps = 1e-10
            s_inv_sqrt = np.diag(1.0 / np.sqrt(np.clip(s, eps, None)))
            X = U @ s_inv_sqrt @ U.T
            H_eff = X @ H @ X
            S_eff = None
            print("#DEBUG foe orthogonalize basis -> S=I")
        else:
            S_eff = S
    else:
        S_eff = None

    lo, hi, mid, span = _estimate_spectral_bounds(H_eff)
    Hs, a, c = _scale_hamiltonian(H_eff, lo, hi)
    mu_val = mid if mu is None else mu
    mu_tilde = (mu_val - c) / a
    beta_eff = beta if beta_is_scaled else beta * a  # scale physical beta (1/eV) by energy scale 'a'
    coeffs = chebyshev_coeffs_fermi(mu_tilde=mu_tilde, beta=beta_eff, n_poly=n_poly)
    eta = np.random.choice([-1.0, 1.0], size=(N, R))

    # if orthogonalized, move probes into orthonormal basis
    if X is not None:
        eta_eff = X.T @ eta
    else:
        eta_eff = eta

    def apply_op(V):
        tmp = Hs @ V
        if S_eff is None:
            return tmp
        return _solve_linear(S_eff, tmp, solver=solver, jacobi_steps=jacobi_steps)

    def to_full(phi_eff):
        return X @ phi_eff if X is not None else phi_eff

    def to_density(phi_eff):
        """Map phi_eff to physical density vector; apply S if Mulliken requested."""
        phi_full = to_full(phi_eff)
        if calc_mulliken and S is not None:
            phi_full = S @ phi_full
        return np.mean(eta * phi_full, axis=1)

    def _print_diag(k, phi_eff):
        if ref_rho is None or diag_every <= 0 or (k % diag_every) != 0:
            return
        rho_k = to_density(phi_eff)
        diff = rho_k - ref_rho
        rmse = np.sqrt(np.mean(diff**2))
        maxdiff = np.max(np.abs(diff))
        print(f"{diag_prefix} k={k:3d} rmse={rmse:.6f} maxdiff={maxdiff:.6f}")

    # If diagnostics requested, force forward recurrence to expose steps
    use_clenshaw = False if (ref_rho is not None and diag_every and diag_every > 0) else use_clenshaw

    if use_clenshaw:
        # Clenshaw evaluation for numerical stability
        b_next = np.zeros_like(eta_eff)
        b_curr = np.zeros_like(eta_eff)
        for k in range(n_poly - 1, 0, -1):
            b_prev = 2.0 * apply_op(b_curr) - b_next + coeffs[k] * eta_eff
            b_next, b_curr = b_curr, b_prev
        phi_eff = coeffs[0] * eta_eff + apply_op(b_curr) - b_next
        print("#DEBUG foe using Clenshaw recurrence")
    else:
        # forward Chebyshev recurrence (less stable)
        v_prev = eta_eff
        v_curr = apply_op(eta_eff)
        phi_eff = coeffs[0] * v_prev + coeffs[1] * v_curr
        _print_diag(1, phi_eff)
        for n in range(2, n_poly):
            v_next = 2.0 * apply_op(v_curr) - v_prev
            phi_eff += coeffs[n] * v_next
            _print_diag(n, phi_eff)
            v_prev, v_curr = v_curr, v_next

    rho_est = to_density(phi_eff)
    
    info = dict(lo=lo, hi=hi, mid=mid, span=span, a=a, c=c, mu=mu_val, 
                mu_tilde=mu_tilde, beta_in=beta, beta_eff=beta_eff, 
                solver=solver, orthogonalize=orthogonalize, 
                clenshaw=use_clenshaw, calc_mulliken=calc_mulliken)
    
    print(f"#DEBUG foe_stochastic_density N={N} R={R} n_poly={n_poly} "
          f"beta_in={beta} beta_eff={beta_eff:.1f} mu={mu_val:.4f} "
          f"solver={solver} mulliken={calc_mulliken}")
          
    return rho_est, info
