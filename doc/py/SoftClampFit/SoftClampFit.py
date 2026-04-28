
import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# Part 0: Provided Functions and Potentials
# ==============================================================================

def soft_clamp(y, dy, y1, y2):
    """
    Applies a soft clamp to y, smoothly transitioning values above y1 towards y2.
    Also computes the derivative dy accordingly using the chain rule.
    """
    y_new  = y.copy()
    dy_new = dy.copy()
    mask   = y > y1
    if not np.any(mask):
        return y_new, dy_new
    
    y12    = y2 - y1
    invdy  = 1.0 / y12
    # Ensure z is an array for proper indexing
    z = (y[mask] - y1) * invdy
    
    y_new[mask]   = y1 + y12 * (1 - 1 / (1 + z))
    dy_new[mask] *= 1.0 / (1.0 + z)**2
    return y_new, dy_new

def morse_potential(r, R0, E0, alpha):
    """ Reference potential with a soft repulsive wall. """
    return E0 * (np.exp(-2 * alpha * (r - R0)) - 2 * np.exp(-alpha * (r - R0)))

def lj_potential(r, R0, E0):
    """ Model potential with a hard 1/r^12 repulsive wall. """
    ir = 1.0 / r
    u = R0 * ir
    u6 = u**6
    return E0 * (u6**2 - 2 * u6)



# ==============================================================================
# Part 2: Function for Variational Derivatives
# ==============================================================================
print("--- Part 2: Calculating Loss and Variational Derivatives ---")

def calculate_loss_and_grads(params, r_data, E_ref_data, clamp_y1, clamp_y2):
    """
    Calculates the total clamped loss and its analytical gradients
    with respect to the model parameters (R0, E0).
    """
    R0, E0 = params

    # --- 1. FORWARD PASS: Calculate Loss ---
    
    # Model energy and error
    E_model = lj_potential(r_data, R0, E0)
    error = E_model - E_ref_data

    # Initialize containers for loss and its derivative w.r.t model energy
    loss_contributions = np.zeros_like(error)
    d_loss_d_E_model   = np.zeros_like(error)

    # a) Negative error (model is too attractive) -> standard quadratic loss
    neg_mask = error < 0
    loss_contributions[neg_mask] = error[neg_mask]**2
    d_loss_d_E_model[neg_mask]   = 2 * error[neg_mask]

    # b) Positive error (model is too repulsive) -> clamped loss
    pos_mask = error >= 0
    # We apply clamp to the error itself, then square
    # The derivative dy/d(error) is 1.0
    y_in = error[pos_mask]
    dy_in = np.ones_like(y_in)
    
    clamped_err, d_clamped_err_d_err = soft_clamp(y_in, dy_in, clamp_y1, clamp_y2)
    
    loss_contributions[pos_mask] = clamped_err**2

    # --- 2. BACKWARD PASS: Calculate Gradients ---
    
    # Derivative of loss w.r.t. clamped error
    d_loss_d_clamped_err = 2 * clamped_err
    # Chain rule: derivative w.r.t original error
    d_loss_d_err = d_loss_d_clamped_err * d_clamped_err_d_err
    # Chain rule: derivative w.r.t. model energy
    d_loss_d_E_model[pos_mask] = d_loss_d_err

    # Derivatives of the LJ potential w.r.t. its parameters
    ir = 1.0 / r_data
    u = R0 * ir
    u6 = u**6
    
    dE_dE0 = u6**2 - 2 * u6  # dE/dE0 = E/E0
    # From C++ snippet: dE_dR0 = 12.f * (E0/R0) * u6 * (u6 - 1.f);
    dE_dR0 = 12.0 * (E0 / R0) * u6 * (u6 - 1.0)

    # Final chain rule step to get gradient w.r.t. parameters
    # Sum over all data points
    grad_R0 = np.sum(d_loss_d_E_model * dE_dR0)
    grad_E0 = np.sum(d_loss_d_E_model * dE_dE0)

    total_loss = np.sum(loss_contributions)
    
    return total_loss, np.array([grad_R0, grad_E0])

# ==============================================================================
# Part 3: Visual check of numerical vs analytical derivatives along a path
# ==============================================================================
def plot_loss_and_deriv_check(param_samples, r_data, E_ref_data, clamp_y1, clamp_y2):
    """Given an array of parameter samples of shape (N,2) for (R0,E0),
    plot the loss L along the sequence and compare numerical vs analytical dL/dn.
    The analytical derivative uses chain rule: dL/dn = grad . dparams/dn.
    """
    P = np.asarray(param_samples)
    N = P.shape[0]
    L = np.zeros(N)
    G = np.zeros((N,2))
    for i in range(N):
        L[i], G[i] = calculate_loss_and_grads(P[i], r_data, E_ref_data, clamp_y1, clamp_y2)

    # Numerical derivative w.r.t. sample index (unit spacing)
    dL_num = np.zeros(N)
    if N >= 2:
        dL_num[1:-1] = (L[2:] - L[:-2]) * 0.5
        dL_num[0]    = L[1] - L[0]
        dL_num[-1]   = L[-1] - L[-2]

    # Analytical derivative along the discrete path
    dL_an  = np.zeros(N)
    if N >= 2:
        for i in range(1, N-1):
            dp = (P[i+1] - P[i-1]) * 0.5
            dL_an[i] = G[i].dot(dp)
        dL_an[0]  = G[0].dot(P[1] - P[0])
        dL_an[-1] = G[-1].dot(P[-1] - P[-2])

    # RMSE between numerical and analytical derivative along path
    rmse = float(np.sqrt(np.mean((dL_num - dL_an)**2))) if N>1 else float('nan')

    # Plot L and derivatives
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax1.plot(L, 'k-', label='L(params[n])')
    ax1.set_ylabel('Loss L')
    ax1.set_title('Loss along parameter path')
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend()

    ax2.plot(dL_num, ':b', lw=1.5, label='dL/dn (numerical)')
    ax2.plot(dL_an,  '-r', lw=0.5, label='dL/dn (analytical)')
    ax2.set_xlabel('sample index n')
    ax2.set_ylabel('dL/dn')
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend()
    plt.tight_layout()
    plt.show()
    print(f"Derivative check RMSE: {rmse:.3e}")
    return rmse

if __name__ == "__main__":

    import argparse

    # CLI arguments
    parser = argparse.ArgumentParser(description='SoftClampFit demo: clamp on dE and derivative check.')
    parser.add_argument('--R0', type=float, default=3.0, help='Reference/model R0 (Angstrom).')
    parser.add_argument('--E0', type=float, default=1.0, help='Reference/model E0 (kcal/mol).')
    parser.add_argument('--alpha', type=float, default=1.5, help='Morse steepness parameter.')
    parser.add_argument('--rmin', type=float, default=2.0, help='Min distance.')
    parser.add_argument('--rmax', type=float, default=8.0, help='Max distance.')
    parser.add_argument('--nr', type=int, default=200, help='Number of samples.')
    parser.add_argument('--clamp-start', default=0.0, type=float, dest='clamp_start', help='Start clamping when dE>this.')
    parser.add_argument('--clamp-limit', default=1.0, type=float, dest='clamp_limit', help='Asymptotic clamp limit for dE.')
    parser.add_argument('--no-deriv-check', action='store_true', help='Disable derivative check plot.')
    parser.add_argument('--nsamp', type=int, default=201, help='Samples along parameter path for deriv check.')
    parser.add_argument('--dR0', type=float, default=0.2, help='R0 span half-range around center for deriv check.')
    parser.add_argument('--dE0', type=float, default=0.5, help='E0 span half-range around center for deriv check.')
    parser.add_argument('--emin', type=float, default=-1.0, help='Y-min for energy plot (top panel).')
    parser.add_argument('--emax', type=float, default=+3.0, help='Y-max for energy plot (top panel).')
    args = parser.parse_args()

    # ==============================================================================
    # Part 1: Illustrating the Effect of Soft Clamping on the Error
    # ==============================================================================
    print("--- Part 1: Visualizing the Clamping Effect ---")

    # Define shared parameters
    R0 = args.R0  # Angstrom
    E0 = args.E0  # kcal/mol
    alpha = args.alpha # Morse steepness parameter (relatively soft)

    # Generate a range of distances
    r = np.linspace(args.rmin, args.rmax, args.nr)

    # Calculate reference (Morse) and model (LJ) energies
    E_ref = morse_potential(r, R0, E0, alpha)
    E_model = lj_potential(r, R0, E0)

    # Calculate the raw error (model - ref) and squared error
    error = E_model - E_ref
    squared_error = error**2

    # Apply asymmetric soft clamping to the ERROR (not squared)
    # Clamp only where error > 0 (model too repulsive)
    clamp_start_error = args.clamp_start
    clamp_limit_error = args.clamp_limit

    clamped_error = error.copy()
    pos_error_mask = (error > 0)
    if np.any(pos_error_mask):
        y_in  = error[pos_error_mask]
        dy_in = np.ones_like(y_in)
        y_out, _ = soft_clamp(y_in, dy_in, clamp_start_error, clamp_limit_error)
        clamped_error[pos_error_mask] = y_out

    # Squared errors: raw and from clamped dE
    clamped_squared_error = clamped_error**2

    # Print stats
    print(f"Ranges: r=[{args.rmin}, {args.rmax}] nr={args.nr}")
    print(f"Clamp: start={clamp_start_error}, limit={clamp_limit_error}")
    print(f"E_ref:  min={E_ref.min():.3f}, max={E_ref.max():.3f}")
    print(f"E_model:min={E_model.min():.3f}, max={E_model.max():.3f}")
    print(f"dE raw: min={error.min():.3f}, max={error.max():.3f}")
    print(f"dE clp: min={clamped_error.min():.3f}, max={clamped_error.max():.3f}")
    npos = int((error>0).sum()); ncl = int((error>clamp_start_error).sum())
    print(f"pos dE count={npos}/{error.size}, clamped count={ncl}/{error.size}")
    near_limit = clamped_error.max() / clamp_limit_error if clamp_limit_error>0 else np.nan
    print(f"max(clamped dE)/limit = {near_limit:.3f}")


    # Plotting
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Effect of Soft Clamping on Fitting Error", fontsize=16)

    # Top: Energies + Errors (linear scale)
    ax1.plot(r, E_ref, 'g-', label='E_ref (Morse)')
    ax1.plot(r, E_model, 'b--', label='E_model (LJ)')
    ax1.plot(r, error, color='orange', lw=1.5, label='dE = E_model - E_ref')
    ax1.plot(r, clamped_error, 'm-', lw=2, label='clamped dE (pos only)')
    ax1.set_ylabel('Energy / Error (kcal/mol)')
    ax1.axhline(0, color='grey', linestyle=':', lw=1)
    if args.emin is not None or args.emax is not None:
        ymin = ax1.get_ylim()[0] if args.emin is None else args.emin
        ymax = ax1.get_ylim()[1] if args.emax is None else args.emax
        ax1.set_ylim(ymin, ymax)
    ax1.legend(ncol=2)
    ax1.set_title('Potentials with Error Overlay (Clamp on dE)')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Bottom: Squared Error (log scale)
    ax2.plot(r, squared_error, 'r-', label='Raw (dE)^2')
    ax2.plot(r, clamped_squared_error, 'm-', lw=2, label='(clamped dE)^2')
    ax2.set_xlabel('Distance (Angstrom)')
    ax2.set_ylabel('Squared Error')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-4, 1e5)
    ax2.legend()
    ax2.set_title('Contribution to Loss Function (Log Scale)')
    ax2.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    print(f"Max original squared error: {np.max(squared_error):.2e}")
    print(f"Max clamped squared error:  {np.max(clamped_squared_error):.2e}\n")


    # --- Example Usage ---
    # Use the same data as in Part 1
    initial_params = [R0, E0]
    # Use the error thresholds (not squared error) for the function
    clamp_start_val = clamp_start_error
    clamp_limit_val = clamp_limit_error

    total_loss, gradients = calculate_loss_and_grads(
        initial_params, r, E_ref, clamp_start_val, clamp_limit_val
    )

    print(f"Initial Parameters: R0 = {initial_params[0]}, E0 = {initial_params[1]}")
    print(f"Total Clamped Loss: {total_loss:.4f}")
    print(f"Variational Derivatives (Gradients):")
    print(f"  dL/dR0 = {gradients[0]:.4f}")
    print(f"  dL/dE0 = {gradients[1]:.4f}")

    # --- Demo: derivative check along a simple line in (R0,E0) ---
    do_deriv_check = (not args.no_deriv_check)
    if do_deriv_check:
        t = np.linspace(-1.0, 1.0, args.nsamp)
        param_path = np.stack([R0 + t*args.dR0, E0 + t*args.dE0], axis=1)
        _rmse = plot_loss_and_deriv_check(param_path, r, E_ref, clamp_start_val, clamp_limit_val)