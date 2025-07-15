import numpy as np
import matplotlib.pyplot as plt

# --- 1. Core Computational and Plotting Utilities ---

def generate_grid(ns=(300, 300), xlims=(-5., 5.), ylims=(-5., 5.)):
    """Generates a 2D mesh grid for our calculations."""
    x = np.linspace(xlims[0], xlims[1], ns[0])
    y = np.linspace(ylims[0], ylims[1], ns[1])
    return np.meshgrid(x, y)

def aggregate(X, Y, atoms, f_kernel, f_inv):
    """
    A universal aggregation engine.

    This function iterates over a grid and a list of atoms, applies a kernel 
    function `f_kernel` for each atom, sums the results, and then applies an 
    inverse function `f_inv`.

    Args:
        X, Y (np.array): The grid coordinates.
        atoms (list of dict): Each atom must have a 'pos' key.
        f_kernel (function): A function lambda r, atom: ... that returns the value to be summed.
        f_inv (function): The inverse function to apply to the sum.
    """
    sum_f = np.zeros_like(X)
    for atom in atoms:
        # Calculate distance from every grid point to the current atom center
        r = np.sqrt((X - atom['pos'][0])**2 + (Y - atom['pos'][1])**2)
        # Add a small epsilon to avoid singularities (e.g., 1/0 or log(0))
        #r += 1e-9
        sum_f += f_kernel(r, atom)
    #with np.errstate(all='ignore'):
    result = f_inv(sum_f)
    return result

def plt_atoms(atoms, ax=None, alpha=0.2):
    for atom in atoms:
        ax.plot(atom['pos'][0], atom['pos'][1], 'wo', markersize=atom['R0'] * Atom_scale, markeredgecolor='k', alpha=alpha)

def plot_field(X, Y, Z, title=None, bCbar=True, atoms=None, ax=None): 
    """A generic 2D field plotter."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))
    # Use percentile clipping for a robust colormap against outliers
    vmax = np.nanpercentile(Z, 99.5)
    vmin = np.nanpercentile(Z, 0.5)
    extent=[X.min(), X.max(), Y.min(), Y.max()]
    im = ax.imshow(Z, origin='lower', cmap='viridis', extent=extent, vmin=vmin, vmax=vmax, interpolation='nearest')
    if bCbar: plt.colorbar(im, ax=ax)
    if atoms is not None: plt_atoms(atoms, ax=ax)
    # # Plot atom positions
    # for atom in atoms:
    #     ax.plot(atom['pos'][0], atom['pos'][1], 'wo', markersize=8, markeredgecolor='k')
    if title is not None: ax.set_title(title, fontsize=14)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_aspect('equal', 'box')
    plt.tight_layout()


# --- 2. The Main Analysis Orchestrator ---

def analyze_and_plot_coeff(X, Y, atoms, potential_model, dist_model=None,  R_agg=None, ax=None ):
    """
    Calculates and plots the effective coefficient for a given combination
    of a potential and a distance aggregation method.
    """
    # Define the identity function for simple summation
    identity_func = lambda s: s
    
    # 1. Calculate E_ref: The "ground truth" potential summed over all atoms.
    #    We use the aggregate function where f_inv is the identity.
    E_ref = aggregate(X, Y, atoms,  f_kernel=potential_model['ref_kernel'],  f_inv=identity_func)

    # 2. Calculate R_agg: The aggregate distance field.
    #    The kernel function for distance doesn't use atom params, so we wrap it.
    # dist_kernel = lambda r, atom: dist_model['f'](r)
    # R_agg = aggregate(X, Y, atoms,  f_kernel=dist_kernel,  f_inv=dist_model['f_inv'])

    if R_agg is None:
        #dist_kernel = lambda r, atom: (r, atom)
        R_agg = aggregate(X, Y, atoms,  f_kernel=dist_model['f'],  f_inv=dist_model['f_inv'])

    # 3. Calculate E_agg: The potential from the aggregate model.
    E_agg = potential_model['agg_form'](R_agg)
    C_eff = E_ref / E_agg
    
    # 4. Calculate the effective coefficient, C_eff.
    #C_eff = np.divide(E_ref, E_agg, out=np.full_like(X, np.nan), where=np.abs(E_agg) > 1e-9)
    
    # 5. Plot the result using the generic plotter.
    if ax is not None:
        title = (f"Effective '{potential_model['coeff_name']}' for {potential_model['name']}")
        plot_field(X, Y, C_eff, title=title, bCbar=True, atoms=atoms, ax=ax)
    return C_eff




# --- 3. Main Script: Define Models and Run Analyses ---

if __name__ == "__main__":

    # H      3.16876   -0.00000    0.06367   +0.10
    # C      2.07201   -0.00000   -0.03311   +0.20
    # O      1.46158    0.00000    1.01488   -0.3
    # O      1.46212    0.00000   -1.18472   -0.3
    # H      0.50000    0.00000   -1.00000   +0.3

    # H    1.4430    0.00190802059        
    # C    1.9255    0.00455323095           
    # N    1.8300    0.00299212319 
    # O    1.7500    0.00260184625 

    Ri=1.4
    Ei = np.sqrt(0.0019080205)

    eH = np.sqrt(0.00190802059)*Ei
    eC = np.sqrt(0.00455323095)*Ei
    eN = np.sqrt(0.00299212319)*Ei
    eO = np.sqrt(0.00260184625)*Ei
    

    atoms = [
        {'pos': ( 3.16876, 0.06367),  'q':  0.1, 'E0': eH, 'R0': 1.4+Ri},  # H
        {'pos': ( 2.07201, -0.03311), 'q':  0.2, 'E0': eC, 'R0': 1.8+Ri},  # C
        {'pos': ( 1.46158, 1.01488),  'q': -0.3, 'E0': eO, 'R0': 1.6+Ri},  # O
        {'pos': ( 1.46212, -1.18472), 'q': -0.3, 'E0': eO, 'R0': 1.6+Ri},  # O
        {'pos': ( 0.50000, -1.00000), 'q':  0.3, 'E0': eH, 'R0': 1.4+Ri}   # H
    ]

    # atoms = [
    #     {'pos': (-2.0, 0), 'q':  0.5, 'E0': 1.0, 'R0': 1.4},
    #     {'pos': ( 2.0, 0), 'q': -0.5, 'E0': 2.0, 'R0': 2.0}
    # ]

    L = 6.0
    X, Y = generate_grid(xlims=(-L, L), ylims=(-L, L))
    
    # --- Define Aggregate Distance Models ---
    b_dist = 2.0
    agg_lse = {
        'name': f'Log-Sum-Exp(α={b_dist})',
        'f':     lambda r, atom:  np.exp(-b_dist * r),
        'f_inv': lambda s: -np.log(s) / b_dist
    }
    p = -4.0
    agg_lp = {
        'name': f'Lp-norm (p={p})',
        'f':     lambda r, atom: r**p,
        'f_inv': lambda s: s**(1/p)
    }
    dist_models = [agg_lse, agg_lp]

    # --- Analysis 1: Coulomb Potential ---
    b_morse = 1.5
    coulomb_model = {
        'name': 'Coulomb', 'coeff_name': 'Charge (q)',
        'ref_kernel': lambda r, atom: atom['q'] / r,
        'agg_form':   lambda R_agg: 1.0 / R_agg
    }

    # --- Analysis 2: Decomposed Morse Potential ---
    for atom in atoms:
        E0, R0 = atom['E0'], atom['R0']
        atom['A'] =     E0 * np.exp(2 * b_morse * R0)
        atom['B'] = 2 * E0 * np.exp(    b_morse * R0)

    # Define the potential models for the A (repulsive) and B (attractive) terms
    morse_A_model = {
        'name': 'Morse Repulsion', 'coeff_name': 'A',
        'ref_kernel': lambda r, atom: atom['A'] * np.exp(-2 * b_morse * r),
        'agg_form':   lambda R_agg:               np.exp(-2 * b_morse * R_agg)
    }
    morse_B_model = {
        'name': 'Morse Attraction', 'coeff_name': 'B',
        'ref_kernel': lambda r, atom: atom['B'] * np.exp(-b_morse * r),
        'agg_form':   lambda R_agg:               np.exp(-b_morse * R_agg)
    }
    
    Atom_scale = 5.0

    # Run the analysis for both terms and both distance models
    for dist_model in dist_models:
        fig, (ax0,ax1, ax2, ax3) = plt.subplots( 1,4 ,figsize=(20,5))
        R_agg = aggregate(X, Y, atoms,  f_kernel=dist_model['f'],  f_inv=dist_model['f_inv'])
        plot_field(X, Y, R_agg, title=dist_model['name'], bCbar=True, atoms=atoms, ax=ax0)
        Q_eff = analyze_and_plot_coeff(X, Y, atoms, coulomb_model, R_agg=R_agg, ax=ax1)
        A_eff = analyze_and_plot_coeff(X, Y, atoms, morse_A_model, R_agg=R_agg, ax=ax2)
        B_eff = analyze_and_plot_coeff(X, Y, atoms, morse_B_model, R_agg=R_agg, ax=ax3)

        # Q_eff = analyze_and_plot_coeff(X, Y, atoms, coulomb_model, dist_model, ax=ax1)
        # A_eff = analyze_and_plot_coeff(X, Y, atoms, morse_A_model, dist_model, ax=ax2)
        # B_eff = analyze_and_plot_coeff(X, Y, atoms, morse_B_model, dist_model, ax=ax3)
    plt.show()