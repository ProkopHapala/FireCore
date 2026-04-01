import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==========================================
# 1. Physical Parameters & Potential
# ==========================================
K_tip = 4.0      # Tip spring stiffness
F_max = 2.8      # Breaking force limit of the tip-molecule bond
HILL_AMP = 3.0   # Gaussian hill amplitude
HILL_SIG = 1.0   # Gaussian hill width (sigma)
Y_CONF = 0.5     # Harmonic confinement coefficient in y

def V_surf(pos):
    """
    2D Surface Potential: 
    Double well in x, harmonic confinement in y, plus a Gaussian hill at the origin.
    """
    x, y = pos
    well_x = (x**2 - 4)**2 / 8.0
    conf_y = Y_CONF * y**2
    hill   = HILL_AMP * np.exp(-(x**2 + y**2) / (2 * HILL_SIG**2))
    return well_x + conf_y + hill

def dV_surf(pos):
    """ Gradient of the surface potential. """
    x, y = pos
    exp_term = np.exp(-(x**2 + y**2) / (2 * HILL_SIG**2))
    pref = HILL_AMP / (HILL_SIG**2)
    dx = 0.5 * x * (x**2 - 4) - pref * x * exp_term
    dy = 2 * Y_CONF * y - pref * y * exp_term
    return np.array([dx, dy])

def V_tot(pos, xt, yt):
    """ Total energy: Surface + Tip harmonic spring """
    x, y = pos
    return V_surf(pos) + 0.5 * K_tip * ((x - xt)**2 + (y - yt)**2)

def dV_tot(pos, xt, yt):
    """ Gradient of the total energy """
    x, y = pos
    grad_surf = dV_surf(pos)
    dx = grad_surf[0] + K_tip * (x - xt)
    dy = grad_surf[1] + K_tip * (y - yt)
    return np.array([dx, dy])

# ==========================================
# 2. Define Trajectories
# ==========================================
# ==========================================
# 3. Simulate Manipulation (Relaxation)
# ==========================================
def simulate_path(x_path, y_path, init_guess):
    mol_pos = []
    forces = []
    E_surf = []
    E_tot = []
    
    current_guess = init_guess
    for xt, yt in zip(x_path, y_path):
        # Relax molecule in the combined potential
        res = minimize(V_tot, current_guess, args=(xt, yt), jac=dV_tot, method='BFGS')
        opt_pos = res.x
        
        # Calculate force on the tip F = K * (r_mol - r_tip)
        fx = K_tip * (opt_pos[0] - xt)
        fy = K_tip * (opt_pos[1] - yt)
        f_mag = np.sqrt(fx**2 + fy**2)
        
        # Store data
        mol_pos.append(opt_pos)
        forces.append(f_mag)
        E_surf.append(V_surf(opt_pos))
        E_tot.append(V_tot(opt_pos, xt, yt))
        
        current_guess = opt_pos # Update guess for the next step (history)
        
    return np.array(mol_pos), np.array(forces), np.array(E_surf), np.array(E_tot)

def main():
    args = parse_args()

    global K_tip, F_max, HILL_AMP, HILL_SIG, Y_CONF
    K_tip = args.stiffness
    F_max = args.break_force
    HILL_AMP = args.hill_amp
    HILL_SIG = args.hill_sigma
    Y_CONF = args.y_conf

    x_tip_path = np.linspace(-args.x_span, args.x_span, args.steps)
    y_tip_straight = np.zeros_like(x_tip_path)
    y_tip_curved = args.serp_amp * np.cos(np.pi * x_tip_path / (2 * args.x_span))

    start_guess = [-args.x_span, 0.0]

    mol_straight, f_straight, e_surf_straight, e_tot_straight = simulate_path(x_tip_path, y_tip_straight, start_guess)
    mol_curved, f_curved, e_surf_curved, e_tot_curved = simulate_path(x_tip_path, y_tip_curved, start_guess)

    fig = plt.figure(figsize=(14, 7))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1])

    # --- PANEL 1: 2D Top View (Left) ---
    ax_2d = fig.add_subplot(gs[:, 0])

    X, Y = np.meshgrid(np.linspace(-args.x_plot, args.x_plot, args.grid),
                       np.linspace(args.y_plot_min, args.y_plot_max, args.grid))
    Z = np.array([V_surf([x, y]) for x, y in zip(X.ravel(), Y.ravel())]).reshape(X.shape)
    contour = ax_2d.contourf(X, Y, Z, levels=40, cmap='gray', alpha=0.9)
    fig.colorbar(contour, ax=ax_2d, label='Surface Potential Energy')

    # Minima markers
    #ax_2d.plot([-2.05, 2.05], [0, 0], 'y*', markersize=10, label='Minima')

    # Trajectories with markers (points visible)
    ax_2d.plot(x_tip_path, y_tip_straight, 'r-', alpha=0.7,   lw=1.0, label='Tip (Straight)')
    ax_2d.plot(mol_straight[:, 0], mol_straight[:, 1], 'ro-', lw=2.0, ms=10, label='Molecule (Straight)')
    for xt, yt, xm, ym in zip(x_tip_path, y_tip_straight, mol_straight[:,0], mol_straight[:,1]):
        ax_2d.plot([xt, xm], [yt, ym], color='r', lw=0.8, alpha=0.4)

    ax_2d.plot(x_tip_path, y_tip_curved, 'b-', alpha=0.7, lw=1.0, label='Tip (Curved)')
    ax_2d.plot(mol_curved[:, 0], mol_curved[:, 1], 'bo-', lw=2.0, ms=10, label='Molecule (Curved)')
    for xt, yt, xm, ym in zip(x_tip_path, y_tip_curved, mol_curved[:,0], mol_curved[:,1]):
        ax_2d.plot([xt, xm], [yt, ym], color='b', lw=0.8, alpha=0.4)

    ax_2d.set_xlim([-args.x_plot, args.x_plot])
    ax_2d.set_ylim([args.y_plot_min, args.y_plot_max])
    ax_2d.set_title('Top View: 2D Potential and Manipulation Paths')
    ax_2d.set_xlabel('x position')
    ax_2d.set_ylabel('y position')
    ax_2d.legend(loc='upper left')

    # --- PANEL 2: Energy along Path (Top Right) ---
    ax_E = fig.add_subplot(gs[0, 1])
    if args.energy_mode in ('total', 'both'):
        ax_E.plot(x_tip_path, e_tot_straight,  'ro-', lw=2.0,ms=10,  label='Total Energy Straight')
        ax_E.plot(x_tip_path, e_tot_curved,    'bo-', lw=2.0,ms=10,  label='Total Energy Curved')
    if args.energy_mode in ('surface', 'both'):
        ax_E.plot(x_tip_path, e_surf_straight, 'r--', lw=1.0, marker='o', markersize=3, label='Surface Energy Straight')
        ax_E.plot(x_tip_path, e_surf_curved,   'b--', lw=1.0, marker='o', markersize=3, label='Surface Energy Curved')

    ax_E.set_title('Energy along Trajectory')
    ax_E.set_ylabel('Energy')
    ax_E.grid(True, alpha=0.3)
    ax_E.legend(loc='upper left')

    # --- PANEL 3: Force on Tip (Bottom Right) ---
    ax_F = fig.add_subplot(gs[1, 1], sharex=ax_E)
    ax_F.plot(x_tip_path, f_straight, 'ro-', lw=2.0, ms=10,  label='Tip force straight')
    ax_F.plot(x_tip_path, f_curved, 'bo-', lw=2.0, ms=10, label='Tip force curved')

    ax_F.axhline(F_max, color='red', linestyle='--', lw=2, label='Breaking Force Limit')
    ax_F.fill_between(x_tip_path, F_max, max(f_straight.max()+0.5, f_curved.max()+0.5), color='red', alpha=0.15, hatch='//')

    break_idx = np.argmax(f_straight > F_max)
    if break_idx > 0:
        ax_F.annotate('Bond Breaks!',
                      xy=(x_tip_path[break_idx], f_straight[break_idx]),
                      xytext=(x_tip_path[break_idx] - 1.0, f_straight[break_idx] + 0.5),
                      arrowprops=dict(facecolor='red', shrink=0.05),
                      color='red', fontweight='bold')

    ax_F.set_title('Required Force applied by AFM Tip')
    ax_F.set_xlabel('Tip Forward Progress ($x_{tip}$)')
    ax_F.set_ylabel('Force Magnitude $|F_{tip}|$')
    ax_F.set_ylim(0, max(f_straight.max(), f_curved.max()) + 0.7)
    ax_F.grid(True, alpha=0.3)
    ax_F.legend(loc='upper left')

    plt.tight_layout()
    if args.output:
        fname = args.output if args.output.lower().endswith(f'.{args.format}') else f"{args.output}.{args.format}"
        plt.savefig(fname, format=args.format, dpi=300)
    else:
        plt.show()


def parse_args():
    p = argparse.ArgumentParser(description='Serpentine path illustration with force constraints')
    p.add_argument('--stiffness',   type=float, default=3.3, help='tip spring stiffness K_tip')
    p.add_argument('--break-force', type=float, default=4.0, help='maximum sustainable force before bond breaks')
    p.add_argument('--steps',       type=int,   default=100, help='number of tip positions along path')
    p.add_argument('--x-span',      type=float, default=2.5, help='half-span of tip travel in x (sweeps -x to +x)')
    p.add_argument('--serp-amp',    type=float, default=1.5, help='amplitude of curved (serpentine) y-path')
    p.add_argument('--x-plot',      type=float, default=3.0, help='half-range of x shown in top view')
    p.add_argument('--y-plot-min',  type=float, default=-0.5, help='min y shown in top view')
    p.add_argument('--y-plot-max',  type=float, default=2.5, help='max y shown in top view')
    p.add_argument('--grid',        type=int,   default=200, help='grid resolution for surface contour')
    p.add_argument('--hill-amp',    type=float, default=10.0, help='Gaussian hill amplitude')
    p.add_argument('--hill-sigma',  type=float, default=1.0, help='Gaussian hill sigma')
    p.add_argument('--y-conf',      type=float, default=0.5, help='y confinement coefficient (prefactor on y^2)')
    p.add_argument('--energy-mode', choices=['total', 'surface', 'both'], default='total', help='which energy curves to plot')
    p.add_argument('--output',      type=str, default="serpent", help='path to save figure; if not set, show interactively')
    p.add_argument('--format',      type=str, choices=['png', 'svg'], default='svg', help='output format when saving')
    return p.parse_args()


if __name__ == '__main__':
    main()
