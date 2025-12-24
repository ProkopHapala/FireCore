import numpy as np

from OrderN import _estimate_spectral_bounds, _scale_hamiltonian


def get_contour_poles(emin, e_fermi, n_poles=32, margin_e=1.0):
    """
    Semi-circular contour in upper half-plane:
    - highest real point at Fermi (e_fermi)
    - lowest real point below min eigen: emin - margin_e (fixed energy offset)
    """
    e_bottom = emin - margin_e
    center = 0.5 * (e_fermi + e_bottom)
    R = 0.5 * (e_fermi - e_bottom)
    x_gl, w_gl = np.polynomial.legendre.leggauss(n_poles)
    theta = 0.5 * (x_gl + 1.0) * np.pi  # 0..pi upper arc
    dt_dx = 0.5 * np.pi
    z_nodes = center + R * np.exp(1j * theta)
    dz_dt = 1j * R * np.exp(1j * theta)
    complex_weights = w_gl * dt_dx * dz_dt * (-1.0 / np.pi)
    print(f"#DEBUG contour: e_top={e_fermi:.6f} e_bottom={e_bottom:.6f} R={R:.6f} margin_e={margin_e}")
    return z_nodes, complex_weights


def _solve_linear_multi_rhs_dense(A, B):
    return np.linalg.solve(A, B)


def greens_function_probing(
    H,
    S=None,
    mu=0.0,
    emin=None,
    n_poles=32,
    probing_distance=12,
    solver="dense",
    return_probe_idx=None,
    store_probes=False,
    margin_e=1.0,
):
    """
    Deterministic probing Green-function contour integration (no randomness).
    Returns rho (real array).
    """
    N = H.shape[0]
    if emin is None:
        emin = np.min(np.diag(H)) - 5.0
    z_nodes, weights = get_contour_poles(emin, mu, n_poles=n_poles, margin_e=margin_e)
    rho_accum = np.zeros(N, dtype=np.float64)
    collect_probe = (return_probe_idx is not None) or store_probes
    probe_data = None
    all_probes = [] if store_probes else None
    first_z = True

    if collect_probe:
        psel = int(return_probe_idx) % max(1, probing_distance)
        impulses = []
        responses = []
        weights_used = []
        eta_probe = None
        rho_probe = np.zeros(N, dtype=np.float64)

    for z, w in zip(z_nodes, weights):
        A = z * S - H if S is not None else z * np.eye(N) - H
        for p in range(probing_distance):
            sources = np.arange(p, N, probing_distance)
            if sources.size == 0:
                continue
            rhs = np.zeros((N, sources.size), dtype=np.complex128)
            rhs[sources, np.arange(sources.size)] = 1.0
            if store_probes and first_z:
                comb = np.zeros(N, dtype=np.complex128)
                comb[sources] = 1.0
                all_probes.append(comb)
            if solver == "dense":
                X = _solve_linear_multi_rhs_dense(A, rhs)
            else:
                raise ValueError(f"Unsupported solver {solver}")
            # Mulliken: diag( G S ) approximated from each source column
            diag_contrib = np.zeros(sources.size, dtype=np.complex128)
            for col, k in enumerate(sources):
                if S is not None:
                    diag_contrib[col] = np.dot(X[:, col], S[:, k])  # Mulliken column k
                else:
                    diag_contrib[col] = X[k, col]  # diagonal element of G
            rho_accum[sources] -= np.imag(w * diag_contrib)  # sign flip to make density positive

            if collect_probe and p == psel:
                # use whole comb pattern as the probe
                eta_probe = np.zeros(N, dtype=np.complex128)
                eta_probe[sources] = 1.0
                resp = X @ np.ones(sources.size, dtype=np.complex128)
                impulses.append(eta_probe)
                responses.append(resp)
                weights_used.append(w)
                resp_phys = S @ resp if S is not None else resp
                rho_probe -= np.imag(w * (eta_probe * resp_phys))
        first_z = False
    if collect_probe:
        probe_data = dict(
            eta=eta_probe if eta_probe is not None else None,
            responses=responses,
            weights=np.array(weights_used) if responses else np.array([]),
            rho_probe=rho_probe if 'rho_probe' in locals() else None,
            probes=all_probes,
        )
    return rho_accum, probe_data


def greens_function_random(
    H,
    S=None,
    mu=0.0,
    emin=None,
    n_poles=32,
    n_random=32,
    solver="dense",
    return_probe_idx=None,
    store_probes=False,
    margin_e=1.0,
):
    """
    Stochastic (random probe) Green-function contour integration.
    If return_probe_idx is set, returns probe-specific impulse/response and rho_probe.
    """
    N = H.shape[0]
    if emin is None:
        emin = np.min(np.diag(H)) - 5.0
    z_nodes, weights = get_contour_poles(emin, mu, n_poles=n_poles, margin_e=margin_e)
    etas = np.random.choice([-1.0, 1.0], size=(N, n_random))
    if S is not None:
        RHS = S @ etas
    else:
        RHS = etas
    rho_accum = np.zeros(N, dtype=np.float64)

    collect_probe = return_probe_idx is not None
    probe_data = None
    probes_all = [] if store_probes else None
    if collect_probe:
        pidx = int(return_probe_idx) % n_random
        impulses = []
        responses = []
        eta_probe = etas[:, pidx]
        rhs_probe = RHS[:, pidx]

    for z, w in zip(z_nodes, weights):
        A = z * S - H if S is not None else z * np.eye(N) - H
        if solver == "dense":
            X = _solve_linear_multi_rhs_dense(A, RHS)
        else:
            raise ValueError(f"Unsupported solver {solver}")
        diag_est = np.mean(etas * X, axis=1)  # Mulliken if RHS=S@eta
        rho_accum -= np.imag(w * diag_est)
        if collect_probe:
            responses.append(X[:, pidx])
            impulses.append(rhs_probe)
        if store_probes:
            probes_all.append(etas[:, 0].copy())

    if collect_probe:
        rho_probe = np.zeros(N, dtype=np.float64)
        for w, resp in zip(weights, responses):
            resp_phys = S @ resp if S is not None else resp
            rho_probe -= np.imag(w * (eta_probe * resp_phys))
        probe_data = dict(
            eta=eta_probe,
            rhs=rhs_probe,
            responses=responses,
            weights=weights,
            rho_probe=rho_probe,
            probes=probes_all,
        )
    return rho_accum, probe_data


def pole_plot_data(eigs, poles):
    """
    Helper to prepare data for plotting poles vs spectrum.
    Returns dict with eigs (real), pole_re, pole_im.
    """
    return {
        "eigs": np.asarray(eigs),
        "pole_re": np.real(poles),
        "pole_im": np.imag(poles),
    }
