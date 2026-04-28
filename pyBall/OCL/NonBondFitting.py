# fitting_driver.py
from operator import iand
import sys
import os
import re
import numpy as np
import pyopencl as cl
import time
from sympy.ntheory import is_abundant
from .OpenCLBase import OpenCLBase
from .FittingDriver import FittingDriver

def run_derivative_test(
    model_name='MODEL_MorseQ_PAIR',
    energy_model_name='ENERGY_MorseQ_PAIR',
    atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
    xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/porcessed/HF-A1_HF-D1.xyz",
    dof_file="/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat",
    eps=1e-4,
    verbose=0,
):
    """Compile templated derivative and energy kernels, then compare analytical
    gradient from get_forces() with finite-difference gradient from evaluate_objective()."""
    print("\n--- Running derivative test ---")
    drv = FittingDriver(verbose=verbose)
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)
    drv.load_dofs(dof_file)
    drv.init_and_upload()

    # Inject both derivative and energy snippets in one build
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_der = extract_macro_block(forces_path, model_name)
    macro_en  = extract_macro_block(forces_path, energy_model_name)
    macros = { 'MODEL_PAIR_ACCUMULATION': macro_der, 'MODEL_PAIR_ENERGY': macro_en }
    drv.compile_with_model(macros=macros, bPrint=False)
    # Bind energy kernel args (derivative kernel args are bound in set_kernel_args() during init/compile)
    drv.set_energy_kernel_args()

    x0 = np.array([d['xstart'] for d in drv.dof_definitions], dtype=np.float32)
    g_an = drv.get_forces(x0)
    g_fd = drv.evaluate_grad_fd(x0, eps=eps)

    diff = g_an - g_fd
    nrm2 = float(np.linalg.norm(diff))
    linf = float(np.max(np.abs(diff))) if diff.size>0 else 0.0
    rel  = nrm2 / (float(np.linalg.norm(g_fd)) + 1e-12)
    print(f"Derivative test: L2={nrm2:.3e} Linf={linf:.3e} Rel={rel:.3e}")
    if verbose>0:
        print("g_an[:8]", g_an[:8])
        print("g_fd[:8]", g_fd[:8])
    return g_an, g_fd

def optimizer_FIRE(driver, initial_dofs, max_steps=1000, dt_start=0.01, fmax=1e-4):
    """
    Performs optimization using the FIRE (Fast Inertial Relaxation Engine) algorithm.
    """
    print("\n--- Starting FIRE Optimization ---")
    # FIRE parameters
    N_min=5; finc=1.1; fdec=0.5; alpha_start=0.1; fa=0.99
    
    dofs = initial_dofs.copy()
    vels = np.zeros_like(dofs)
    alpha = alpha_start
    dt = dt_start
    steps_since_negative_P = 0

    for i_step in range(max_steps):
        forces = driver.get_forces(dofs)
        force_norm = np.linalg.norm(forces)

        if i_step % 10 == 0:
            print(f"FIRE Step {i_step:04d} | Force Norm = {force_norm:.6e} | dt = {dt:.4e} | alpha = {alpha:.4f}")

        if force_norm < fmax:
            print(f"FIRE converged in {i_step} steps. Final force norm: {force_norm:.6e}")
            break

        P = np.dot(forces, vels)

        vels_norm = np.linalg.norm(vels)
        if P > 0:
            steps_since_negative_P += 1
            if steps_since_negative_P > N_min:
                dt = min(dt * finc, 0.1) # Max dt
                alpha *= fa
            # MD step
            vels = (1.0 - alpha) * vels + alpha * forces / force_norm * vels_norm
        else:
            steps_since_negative_P = 0
            dt *= fdec
            alpha = alpha_start
            vels[:] = 0.0 # Stop and reset
        
        dofs += vels * dt
    
    if i_step == max_steps - 1:
        print("FIRE finished due to max steps.")

    return dofs

def setup_driver(
    model_name='ENERGY_MorseQ_PAIR',
    atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
    verbose=0,
):
    # Create driver, load atom types, and compile the energy kernel template.
    # Do NOT upload buffers here; that requires XYZ data loaded first.
    drv = FittingDriver(verbose=verbose)
    drv.load_atom_types(atom_types_file)
    base_path   = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code  = extract_macro_block(forces_path, model_name)
    macros      = { 'MODEL_PAIR_ENERGY': macro_code }
    drv.compile_energy_with_model(macros=macros, bPrint=False)
    return drv


def run_energy_imshow(model_name='ENERGY_MorseQ_PAIR',
                      xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/HBr-D1_HBr-A1.xyz",
                      atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                      kcal=False, sym=True, bColorbar=True,
                      verbose=0, show=True, save=None, cmap='bwr', lines=True,
                      rmax=8.0, drv=None):
    """Evaluate energy-only kernel on a packed XYZ and show 3-panel imshow:
    Reference, Model, Difference. Uses helper functions from tests/tFitREQ/split_scan_imshow_new.py
    to parse headers, detect rows, reshape to grid, and plot.
    """
    import matplotlib.pyplot as plt
    from tests.tFitREQ.split_scan_imshow_new import (
        parse_headers_ra, read_scan_atomicutils, parse_xyz_blocks,
        compute_ra_vec, detect_rows_by_r, reshape_to_grid,
        plot_imshow, compute_shift_from_grid
    )

    # 1) Extract reference energies and geometry (r, a) from headers
    Eh, Rh, Ah = parse_headers_ra(xyz_file)
    if Eh.size == 0:
        raise RuntimeError(f"No reference energies parsed from headers: {xyz_file}")

    title=os.path.basename(xyz_file)

    if drv is None:
        drv = setup_driver(model_name=model_name, atom_types_file=atom_types_file, verbose=verbose)
    # Allow user to call drv.load_data() beforehand to set per-atom hacks
    if not hasattr(drv, 'host_ranges'):
        drv.load_data(xyz_file)

    # 2) Evaluate model energies in the same order
    # Now that data is loaded, allocate/upload minimal buffers and bind kernel args
    drv.init_and_upload_energy_only()
    drv.set_energy_kernel_args()


    Em = drv.evaluate_energies()
    if verbose>0: print(f"Model energies: min={np.nanmin(Em):.6e} max={np.nanmax(Em):.6e}")

    if Em.shape[0] != Eh.shape[0]:
        print(f"Warning: Em(n={Em.shape[0]}) != Eh(n={Eh.shape[0]}). Proceeding with min length.")
        n = min(Em.shape[0], Eh.shape[0])
        Em = Em[:n]
        Eh = Eh[:n]
        Rh = Rh[:n]
        Ah = Ah[:n]

    # 3) Derive geometry r/a from atomic positions with header overrides; then row detection
    Es_geo, Ps = read_scan_atomicutils(xyz_file)
    if Es_geo.size == 0:
        # Fallback to local parser for positions if atomicUtils isn't available
        _, _, Ps = parse_xyz_blocks(xyz_file, natoms=None)
    r, a = compute_ra_vec(Ps, signed=True)
    if Rh.size == r.size and np.any(np.isfinite(Rh)):
        r = np.where(np.isfinite(Rh), Rh, r)
    if Ah.size == a.size and np.any(np.isfinite(Ah)):
        a = np.where(np.isfinite(Ah), Ah, a)
    rows, _ = detect_rows_by_r(r)
    Vref, Rg, Arow, rv = reshape_to_grid(Eh, r, a, rows)
    Vmod, _,   _,   _  = reshape_to_grid(Em, r, a, rows)
    if verbose>0: print(f"Grids: Vref shape={Vref.shape}, Vmod shape={Vmod.shape}")

    # 4) Reference alignment: subtract a common reference (asymptotic min at max r col)
    sref = compute_shift_from_grid(Vref)
    print("sref", sref)
    Vref -= sref    # shift just reference, not model!
    Vdif  = Vmod - Vref

    Vref_min = np.nanmin(Vref)
    Vmod_min = np.nanmin(Vmod)

    if verbose>0:
        print("Vref[:,0]", Vref[:,0])
        print("Vref[0,:]", Vref[0,:])
        print("Vmod[:,0]", Vmod[:,0])
        print("Vmod[0,:]", Vmod[0,:])
        print("Vref.shape , min, max", Vref.shape, np.nanmin(Vref), np.nanmax(Vref))
        print("Vmod.shape , min, max", Vmod.shape, np.nanmin(Vmod), np.nanmax(Vmod))

    # 5) Plot 3 panels via helper
    fig, axs = plot_energy_panels( Vref, Vmod, Vdif, rv, Arow, model_name, sym=sym, cmap=cmap, bColorbar=bColorbar, title=title )

    if save is not None:
        try:
            os.makedirs(os.path.dirname(save), exist_ok=True)
        except Exception:
            pass
        plt.savefig(save, dpi=150)
        if verbose:
            print(f"Saved figure to: {save}")
    # 6) Optional line profiles via helper (each dataset computes its own y-limits etc.)
    if lines:
        # Draw ref and model profiles on the same axes; internals computed per-call
        ax = plot_energy_profile(Vref, rv, Arow, [':', ':', ':', ':'], [1.5, 1.5, 1.5, 1.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, kcal=kcal, sym=sym, name='ref')
        ax = plot_energy_profile(Vmod, rv, Arow, ['-', '-', '-', '-'], [0.5, 0.5, 0.5, 0.5], ['radial', 'radial', 'angular', 'min over r per angle'], marker='.', rmax=rmax, kcal=kcal, sym=sym, name='mod', ax=ax)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return Vref, Vmod, Vdif, rv, Arow

# ----------------------------
# Modular plotting helpers
# ----------------------------

def plot_energy_panels(Vref, Vmod, Vdif, rv, Arow, model_name, title=None, sym=True, cmap='bwr', bColorbar=True):
    """Plot 3-panel imshow (Reference, Model, Difference).
    Returns (emin_ref, vmax_ref, emin_mod, vmax_mod, fig, axs).
    Limits are in displayed units (kcal if requested).
    """
    import matplotlib.pyplot as plt
    from tests.tFitREQ.split_scan_imshow_new import plot_imshow
    fig, axs = plt.subplots(1, 3, figsize=(14, 4))
    plot_imshow(Vref, rv, Arow, title='Reference',           ax=axs[0], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    plot_imshow(Vmod, rv, Arow, title=f'Model {model_name}', ax=axs[1], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    plot_imshow(Vdif, rv, Arow, title='Difference',          ax=axs[2], bColorbar=bColorbar, cmap=cmap, bSym=sym, bByMin=True)
    if title is not None: fig.suptitle(title)
    plt.tight_layout()
    return fig, axs


def plot_energy_profile(
    V, rv, Arow,
    ls, lw, labels,
    marker='.', ms=3.0, rmax=None,
    kcal=False, sym=True,
    row_colors=('k','r'), ang_color='g', min_color='b',
    name='', ylims=None, ax=None
):
    """Plot 4 concise profiles for a single dataset V; computes helpers internally.
    Curves mapping:
      0-> radial @ start angle, 1-> radial @ mid angle,
      2-> angular @ radius of V-min, 3-> per-angle minima.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(7, 4))

    fac = 23.060548 if kcal else 1.0
    # Rows to plot
    ny = V.shape[0] if hasattr(V, 'shape') and len(V.shape) > 0 else 0
    rows = sorted({i for i in (0, ny//2) if 0 <= i < ny})
    # Angle mapping and sorting (0..π)
    theta = np.radians(Arow)
    theta = (theta + (np.pi/2.0)) % np.pi
    srt = np.argsort(theta)
    theta_s = theta[srt]
    # Radius column at global minimum of THIS dataset
    try:
        _, ix_ang = np.unravel_index(np.nanargmin(V), V.shape)
    except Exception:
        ix_ang = None
    # Symmetric y-limits per dataset if requested and not overridden
    if (ylims is None) and sym:
        vmin=np.nanmin(V)
        ylims = (vmin, -vmin)

    rr = rv.astype(float)
    # Radial rows
    for j, irow in enumerate(rows):
        ee = (V[irow, :] * fac).astype(float)
        m = np.isfinite(rr) & np.isfinite(ee)
        if rmax is not None:
            m &= (rr <= float(rmax))
        if not np.any(m):
            continue
        yy = ee[m]
        col = row_colors[j if j < len(row_colors) else -1]
        ax.plot(rr[m], yy, color=col, ls=ls[j], lw=lw[j], marker=marker, ms=ms,
                label=f"{name} {labels[j]} @{Arow[irow]:.0f}°")
    # Angular @ ix_ang
    if ix_ang is not None:
        e_ang = (V[:, int(ix_ang)] * fac)
        yy = e_ang[srt]
        lab_r = rv[int(ix_ang)] if rv is not None else np.nan
        ax.plot(theta_s, yy, color=ang_color, ls=ls[2], lw=lw[2], marker=marker, ms=ms,
                label=f"{name} {labels[2]} @ r≈{lab_r:.2f}")
    # Min over r per angle
    try:
        e_min = np.nanmin(V, axis=1) * fac
        yy = e_min[srt]
        ax.plot(theta_s, yy, color=min_color, ls=ls[3], lw=lw[3], marker=marker,
                label=f"{name} {labels[3]}")
    except Exception:
        pass
    if ylims is not None:
        ax.set_ylim(float(ylims[0]), float(ylims[1]))
    ax.grid(True, ls=':')
    ax.set_xlabel('r [Å] / θ [rad]')
    ax.set_ylabel('E [kcal/mol]' if kcal else 'E [eV]')
    ax.legend(loc='best', fontsize=8)
    return ax


def run_templated_example(model_name='MODEL_MorseQ_PAIR',
                          atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                          xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/input_example.xyz",
                          dof_file="/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat"):
    """Demo: compile and run the templated kernel with a selected model macro.
    Keeps the original baseline runnable separately.
    """
    print("\n--- Running templated fitting demo ---")
    drv = FittingDriver()
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)
    drv.load_dofs(dof_file)
    drv.init_and_upload()

    # Load macro snippet from Forces.cl and compile templated kernel
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code = extract_macro_block(forces_path, model_name)
    macros = { 'MODEL_PAIR_ACCUMULATION': macro_code }
    drv.compile_with_model(macros=macros, bPrint=True)

    x0 = np.array([d['xstart'] for d in drv.dof_definitions])
    x_final = optimizer_FIRE(drv, x0)

    print("\n" + "="*40)
    print(f"Templated model: {model_name}")
    print("Final Optimized DOF values:")
    print(x_final)
    return x_final

def run_energy_example(model_name='ENERGY_MorseQ_PAIR',
                       atom_types_file="/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat",
                       xyz_file="/home/prokop/git/FireCore/tests/tFitREQ/HHalogens/HBr-D1_HBr-A1.xyz"):
    """Compile and run the energy-only templated kernel on provided XYZ file.
    model_name: one of ENERGY_LJQH2_PAIR, ENERGY_LJr8QH2_PAIR, ENERGY_MorseQ_PAIR
    """
    print("\n--- Running energy-only evaluation ---")
    drv = FittingDriver()
    drv.load_atom_types(atom_types_file)
    drv.load_data(xyz_file)

    # Build and upload minimal buffers
    drv.init_and_upload_energy_only()

    # Load energy-only macro and compile energy template
    base_path = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.abspath(os.path.join(base_path, "../../cpp/common_resources/cl/Forces.cl"))
    macro_code = extract_macro_block(forces_path, model_name)
    macros = { 'MODEL_PAIR_ENERGY': macro_code }
    drv.compile_energy_with_model(macros=macros, bPrint=True)
    drv.set_energy_kernel_args()

    Emols = drv.evaluate_energies()
    print("\n" + "="*40)
    print(f"Energy model: {model_name}")
    for i, E in enumerate(Emols):
        ref = drv.host_ErefW[i] if hasattr(drv, 'host_ErefW') and len(drv.host_ErefW)>i else np.nan
        print(f"Sample {i}: E = {E:.6f}  (Eref = {ref:.6f})")
    return Emols

def run_baseline_example():
    """Wraps the existing baseline example into a function."""
    print("\n--- Running baseline example ---")
    driver = FittingDriver()
    driver.load_atom_types("/home/prokop/git/FireCore/cpp/common_resources/AtomTypes.dat")
    driver.load_data("/home/prokop/git/FireCore/tests/tFitREQ/input_example.xyz")
    driver.load_dofs("/home/prokop/git/FireCore/tests/tFitREQ/dofSelection_MorseSR.dat")
    driver.init_and_upload()

    initial_dofs = np.array([d['xstart'] for d in driver.dof_definitions])
    final_dofs_fire = optimizer_FIRE(driver, initial_dofs)
    print("\n" + "="*40)
    print("Final Optimized DOF values:")
    print(final_dofs_fire)
    return final_dofs_fire

# --- Example Usage ---
if __name__ == '__main__':
    # run like
    # python -u -m pyBall.OCL.NonBondFitting

    # Default: baseline example (kept intact)
    #run_baseline_example()

    # To run the templated-kernel demo with a selected model macro, uncomment one:
    #run_templated_example('MODEL_LJQH2_PAIR')
    # run_templated_example('MODEL_LJr8QH2_PAIR')
    #run_templated_example('MODEL_MorseQ_PAIR')

    run_energy_imshow('ENERGY_MorseQ_PAIR')
    #print("="*40)