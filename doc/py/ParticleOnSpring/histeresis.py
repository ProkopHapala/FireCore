import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# ==========================================
# 1. Physical Parameters
# ==========================================
# A = 1.0         # Amplitude of the surface potential
# k = 1.0         # Wave vector. Creates a hill at x=0, and minima (wells) at x = -pi and x = +pi
# K_tip = 0.1    # Stiffness of the AFM tip spring. Must be < A*k^2 to observe bifurcation/hysteresis

# ==========================================
# 2. Potential Energy Functions
# ==========================================
def V_surf(x):
    """Double-well modeled over the domain by a simple cosine."""
    return A * np.cos(k * x)

def dV_surf(x):
    """Derivative of the surface potential (Force from surface)."""
    return -A * k * np.sin(k * x)

def V_tot(x, x_tip):
    """Total potential: Surface potential + Harmonic tip spring."""
    return V_surf(x) + 0.5 * K_tip * (x - x_tip)**2

def dV_tot(x, x_tip):
    """Derivative of total potential for optimization."""
    return dV_surf(x) + K_tip * (x - x_tip)


def simulate_paths(x_tip_max: float, n_points: int):
    """Compute forward/backward relaxed positions, total energies, and tip forces."""
    x_tips_fwd = np.linspace(-x_tip_max, x_tip_max, n_points)
    x_tips_bwd = np.linspace(x_tip_max, -x_tip_max, n_points)

    x_mols_fwd = []
    x_opt = -np.pi  # Initial guess for forward path (Left well)
    for xt in x_tips_fwd:
        res = minimize(V_tot, x_opt, args=(xt,), jac=dV_tot, method='BFGS')
        x_opt = res.x[0]
        x_mols_fwd.append(x_opt)

    x_mols_bwd = []
    x_opt = np.pi  # Initial guess for backward path (Right well)
    for xt in x_tips_bwd:
        res = minimize(V_tot, x_opt, args=(xt,), jac=dV_tot, method='BFGS')
        x_opt = res.x[0]
        x_mols_bwd.append(x_opt)

    x_mols_fwd = np.array(x_mols_fwd)
    x_mols_bwd = np.array(x_mols_bwd)
    E_tot_fwd = V_tot(x_mols_fwd, x_tips_fwd)
    E_tot_bwd = V_tot(x_mols_bwd, x_tips_bwd)
    F_tip_fwd = -K_tip * (x_mols_fwd - x_tips_fwd)
    F_tip_bwd = -K_tip * (x_mols_bwd - x_tips_bwd)
    return (x_tips_fwd, x_tips_bwd,
            x_mols_fwd, x_mols_bwd,
            E_tot_fwd, E_tot_bwd,
            F_tip_fwd, F_tip_bwd)


def plot_results(step: int, tip_height: float, connector_lw: float, traj_lw: float, force_lw: float,
                 x_tips_fwd, x_tips_bwd,
                 x_mols_fwd, x_mols_bwd,
                 E_tot_fwd, E_tot_bwd,
                 F_tip_fwd, F_tip_bwd,
                 save_path: str | None = None, fmt: str = 'png',
                 marker_tip_size: float = 6.0, marker_mol_size: float = 6.0,
                 traj_style: str = '-'):
    # ==========================================
    # 4. Visualization
    # ==========================================
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), gridspec_kw={'height_ratios':[1.4, 1.0, 1.2]})

    # --- PANEL 1: Energy landscape + molecule energy connectors ---
    ax_energy = axes[0]

    x_plot = np.linspace(-5, 5, 500)
    ax_energy.plot(x_plot, V_surf(x_plot), 'k-', lw=2, label='Surface Potential $V_{surf}(x)$')
    #ax_energy.axvline(0, color='gray', linestyle='--', alpha=0.5, label='Top of the Hill (Barrier)')

    ax_energy.plot(x_mols_fwd, E_tot_fwd, 'b-', lw=traj_lw, alpha=0.7, label='Forward path')
    ax_energy.plot(x_mols_bwd, E_tot_bwd, 'r-', lw=traj_lw, alpha=0.7, label='Backward path')

    first = True
    for xt, xm, Em in zip(x_tips_fwd[::step], x_mols_fwd[::step], E_tot_fwd[::step]):
        ax_energy.plot([xm, xt], [Em, tip_height], 'b-', alpha=0.5, lw=connector_lw)
        ax_energy.plot(xm, Em, 'bo', markersize=marker_mol_size, alpha=0.7)
        ax_energy.plot(xt, tip_height, 'bv', markersize=marker_tip_size, alpha=0.7)

    for i, (xt, xm, Em) in enumerate(zip(x_tips_bwd[::step], x_mols_bwd[::step], E_tot_bwd[::step])):
        ax_energy.plot([xm, xt], [Em, tip_height], 'r-', alpha=0.5, lw=connector_lw)
        ax_energy.plot(xm, Em, 'ro', markersize=marker_mol_size, alpha=0.7)
        ax_energy.plot(xt, tip_height, 'rv', markersize=marker_tip_size, alpha=0.7)

    ax_energy.set_ylabel('Energy (Surface + Spring)')
    ax_energy.set_title('1) Total Energy and Tip/Molecule Connectors')
    ax_energy.grid(True, alpha=0.3)
    ax_energy.set_xlim([-5, 5])
    ax_energy.legend(loc='lower right')

    # --- PANEL 2: Tip force vs tip position ---
    axes[1].plot(x_tips_fwd, F_tip_fwd, 'b'+traj_style, lw=force_lw, alpha=0.8, label='Tip force fwd')
    axes[1].plot(x_tips_bwd, F_tip_bwd, 'r'+traj_style, lw=force_lw, alpha=0.8, label='Tip force bwd')
    axes[1].set_ylabel('Force on tip')
    axes[1].set_title('2) Tip Force vs Tip Position')
    axes[1].grid(True, alpha=0.3)
    axes[1].set_xlim([-5, 5])
    axes[1].legend(loc='lower right')

    # --- PANEL 3: Hysteresis Loop ---
    axes[2].plot(x_tips_fwd, x_mols_fwd, 'b'+traj_style, lw=2.5, label='Forward Trajectory')
    axes[2].plot(x_tips_bwd, x_mols_bwd, 'r'+traj_style, lw=2.5, label='Backward Trajectory')
    axes[2].plot(x_tips_fwd, x_tips_fwd, 'k--', alpha=0.5, label='Ideal Tracking ($x_{mol} = x_{tip}$)')

    axes[2].set_xlabel('Tip Position $x_{tip}$')
    axes[2].set_ylabel('Relaxed Molecule Position $x_{mol}$')
    axes[2].set_title('3) System Hysteresis (History Dependence)')
    axes[2].legend(loc='lower right')
    axes[2].grid(True, alpha=0.3)
    axes[2].set_xlim([-5, 5])

    plt.tight_layout()
    if save_path:
        fname = save_path if save_path.lower().endswith(f'.{fmt}') else f"{save_path}.{fmt}"
        plt.savefig(fname, format=fmt, dpi=300)
    else:
        plt.show()


def parse_args():
    p = argparse.ArgumentParser(description='AFM tip hysteresis visualization')
    p.add_argument('-A', '--amplitude',  type=float, default=1.0,  help='surface potential amplitude A')
    p.add_argument('-k', '--wavevector', type=float, default=1.0,  help='surface potential wave vector k')
    p.add_argument('-K', '--stiffness',  type=float, default=0.5,  help='tip spring stiffness K_tip')
    p.add_argument('-x', '--tip-max',    type=float, default=4.5,  help='maximum absolute tip displacement (sweep from -x to +x)')
    p.add_argument('-n', '--points',     type=int,   default=300,  help='number of sample points per sweep')
    p.add_argument('-m', '--step',       type=int,   default=10,   help='plot every m-th connector point (>=1)')
    p.add_argument('--connector-lw',     type=float, default=0.5,  help='line width for tip-molecule connectors')
    p.add_argument('--traj-lw',          type=float, default=3.0,  help='line width for energy trajectories')
    p.add_argument('--force-lw',         type=float, default=2.0,  help='line width for tip force curves')
    p.add_argument('--tip-height',       type=float, default=None, help='absolute y-level to draw tip markers; default auto = max energy + gap')
    p.add_argument('--tip-gap',          type=float, default=0.5,  help='extra gap above max energy when tip-height not specified')
    p.add_argument('--output',           type=str,   default="histeresis", help='path to save figure; if not set, show interactively')
    p.add_argument('--format',           type=str,   choices=['png', 'svg'], default='svg', help='output format when saving')
    p.add_argument('--marker-tip-size',  type=float, default=8.0,  help='marker size for tip positions (triangles)')
    p.add_argument('--marker-mol-size',  type=float, default=8.0,  help='marker size for molecule positions (circles)')
    p.add_argument('--traj-style',       type=str,   choices=['-', '--'], default='-', help='line style for force and hysteresis trajectories')
    return p.parse_args()


def main():
    args = parse_args()
    step = max(1, args.step)

    global A, k, K_tip
    A = args.amplitude
    k = args.wavevector
    K_tip = args.stiffness

    (x_tips_fwd, x_tips_bwd,
     x_mols_fwd, x_mols_bwd,
     E_tot_fwd, E_tot_bwd,
     F_tip_fwd, F_tip_bwd) = simulate_paths(args.tip_max, args.points)

    tip_height = args.tip_height
    if tip_height is None:
        tip_height = max(E_tot_fwd.max(), E_tot_bwd.max()) + args.tip_gap

    plot_results(step, tip_height, args.connector_lw, args.traj_lw, args.force_lw,
                 x_tips_fwd, x_tips_bwd,
                 x_mols_fwd, x_mols_bwd,
                 E_tot_fwd, E_tot_bwd,
                 F_tip_fwd, F_tip_bwd,
                 save_path=args.output, fmt=args.format,
                 marker_tip_size=args.marker_tip_size,
                 marker_mol_size=args.marker_mol_size,
                 traj_style=args.traj_style)


if __name__ == '__main__':
    main()