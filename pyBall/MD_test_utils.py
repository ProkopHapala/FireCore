import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from .OCL.MMFF import MMFF
from .OCL.MolecularDynamics import MolecularDynamics
from . import elements

MONITOR_PROPERTY_TYPES = {
    'F': 'vector',
    'T': 'vector',
    'P': 'vector',
    'L': 'vector',
    'Rcm': 'vector',
    'Vcm': 'vector',
    'Fcm': 'vector',
    'Tcm': 'vector',
    'Ekin': 'scalar',
    'Epot': 'scalar',
    'Etotal': 'scalar',
}

DEFAULT_MONITOR_PROPS = ['F', 'P', 'L', 'Rcm', 'Vcm', 'Ekin', 'Epot', 'Etotal']

_ENERGY_CACHE = {
    'prev_vel': None,
    'prev_pos': None,
}

# __all__ = [
#     'MONITOR_PROPERTY_TYPES',
#     'DEFAULT_MONITOR_PROPS',
#     'parse_monitor_props',
#     'get_atom_masses',
#     'compute_energies',
#     'finalize_monitor_series',
#     'summarize_monitor_series',
#     'plot_monitor_series',
#     'build_mmff_from_mol',
#     'zero_dynamic_buffers',
#     'fetch_arrays',
#     'run_step',
#     'compute_totals',
#     'print_totals',
#     'plot_trajectories',
#     'configure_md',
#     'init_observers',
#     'collect_diagnostics',
#     'dump_buffers',
#     'report_conservation',
#     'finalize_recording',
#     'finalize_monitoring',
#     'overlay_bonds',
# ]

def parse_monitor_props(spec, default_props=None, prop_types=None):
    defaults = list(default_props if default_props is not None else DEFAULT_MONITOR_PROPS)
    registry = prop_types if prop_types is not None else MONITOR_PROPERTY_TYPES
    if spec is None:
        return defaults
    s = spec.strip()
    if s == '':
        return defaults
    key = s.lower()
    if key in {'none', 'off', 'false', '0'}:
        return []
    if key == 'all':
        return list(registry.keys())
    props = [p.strip() for p in s.split(',') if p.strip()]
    invalid = [p for p in props if p not in registry]
    if invalid:
        raise ValueError(f"Unknown monitor properties: {', '.join(invalid)}. Valid options: {sorted(registry.keys())}")
    return props


def get_atom_masses(mol, use_real=False):
    if not use_real:
        return np.ones(len(mol.enames), dtype=np.float32)

    masses = np.zeros(len(mol.enames), dtype=np.float32)
    for i, sym in enumerate(mol.enames):
        entry = elements.ELEMENT_DICT.get(sym)
        if entry is None:
            # Try common suffix patterns (e.g. "C_2", "O1") used in some MOL2 exports
            candidates = []
            if '_' in sym:
                candidates.append(sym.split('_')[0])
            stripped_digits = sym.rstrip('0123456789')
            if stripped_digits != sym:
                candidates.append(stripped_digits)
            alpha_only = ''.join(ch for ch in sym if ch.isalpha())
            if alpha_only and alpha_only != sym:
                candidates.append(alpha_only)
            for cand in candidates:
                if cand in elements.ELEMENT_DICT:
                    entry = elements.ELEMENT_DICT[cand]
                    break
        if entry is None:
            raise ValueError(f"No element data for symbol '{sym}' while building masses")
        masses[i] = float(entry[elements.index_mass])
    return masses


def _estimate_true_velocity(avel_atoms, masses):
    """Return velocity estimate closer to central-difference using cached history."""
    v_now = np.asarray(avel_atoms, dtype=np.float32)
    v_prev = _ENERGY_CACHE['prev_vel']
    if v_prev is not None and v_prev.shape == v_now.shape:
        v_est = 0.5 * (v_now + v_prev)
    else:
        v_est = v_now
    _ENERGY_CACHE['prev_vel'] = v_now.copy()
    return v_est


def compute_energies(avel_atoms, masses, aforce_atoms_full, avel_pi=None, inertia_pi=None):
    v = _estimate_true_velocity(avel_atoms, masses)   # use this with Leap-Frong integrator
    m = np.asarray(masses, dtype=np.float32).reshape(-1, 1)
    kin = 0.5 * m * (v * v).sum(axis=1, keepdims=True)
    Ekin_trans = float(kin.sum())

    Ekin_rot = 0.0
    if avel_pi is not None:
        omega = np.asarray(avel_pi, dtype=np.float32)
        if omega.size:
            if inertia_pi is None:
                inertia = 1.0
            else:
                inertia = np.asarray(inertia_pi, dtype=np.float32)
                if inertia.ndim == 0:
                    inertia = float(inertia)
            omega_sq = (omega * omega).sum(axis=1)
            if isinstance(inertia, float):
                Ekin_rot = 0.5 * inertia * float(np.sum(omega_sq))
            else:
                inertia_vec = inertia.reshape(-1)
                if inertia_vec.shape[0] != omega_sq.shape[0]:
                    raise ValueError("inertia_pi length does not match number of pi orbitals")
                Ekin_rot = 0.5 * float(np.dot(inertia_vec, omega_sq))

    Ekin = Ekin_trans + Ekin_rot

    af = np.asarray(aforce_atoms_full, dtype=np.float32)
    # Each kernel assigns interaction energy to both partners; divide by 2 to avoid double counting.
    Epot = float(af[:, 3].sum())
    return {
        'Ekin_trans': Ekin_trans,
        'Ekin_rot': Ekin_rot,
        'Ekin': Ekin,
        'Epot': Epot,
        'Etotal': Ekin + Epot,
    }


def finalize_monitor_series(monitor_data):
    series = {}
    for key, values in monitor_data.items():
        if not values:
            continue
        arrays = [np.array(v, dtype=np.float32) for v in values]
        if arrays[0].ndim == 0:
            series[key] = np.asarray(arrays, dtype=np.float32).reshape(-1)
        else:
            series[key] = np.stack(arrays, axis=0)
    return series


def summarize_monitor_series(series, props):
    if not props: return
    labels = ['x', 'y', 'z']
    header_printed = False
    for name in props:
        data = series.get(name)
        if data is None:
            continue
        arr = np.asarray(data, dtype=np.float32)
        if not header_printed:
            print('Monitor extrema:')
            header_printed = True
        if arr.ndim <= 1:
            mn = float(np.min(arr))
            mx = float(np.max(arr))
            print(f"  {name}: min={mn:.6e} max={mx:.6e}")
        else:
            stats = []
            ncomp = arr.shape[1]
            for j in range(ncomp):
                mn = float(np.min(arr[:, j]))
                mx = float(np.max(arr[:, j]))
                lab = labels[j] if j < len(labels) else str(j)
                stats.append(f"{lab}=({mn:.6e},{mx:.6e})")
            print(f"  {name}: {' '.join(stats)}")


def plot_monitor_series(series, props, show=True):
    if not props:  return None, None
    fig, axes = plt.subplots(len(props), 1, figsize=(8, 2.6 * len(props)), sharex=True)
    if isinstance(axes, np.ndarray):
        axes_list = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        axes_list = list(axes)
    else:
        axes_list = [axes]
    steps = None
    for ax, name in zip(axes_list, props):
        data = series.get(name)
        if data is None:
            ax.set_visible(False)
            continue
        if steps is None:
            steps = np.arange(data.shape[0])
        if data.ndim == 1:
            ax.plot(steps, data, label=name)
        else:
            comps = ['x', 'y', 'z']
            ncomp = data.shape[1]
            for j in range(ncomp):
                comp_label = comps[j] if j < len(comps) else str(j)
                ax.plot(steps, data[:, j], label=f"{name}_{comp_label}")
            ax.legend(loc='best', fontsize=8)
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
    if axes_list:
        axes_list[-1].set_xlabel('step')
    fig.tight_layout()
    if show:
        plt.show()
    return fig, axes_list


def build_mmff_from_mol(mol2_path):
    from .AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=mol2_path)
    # Ensure neighbors present
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    mm = MMFF(bTorsion=False, verbosity=1)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    return mol, mm


def stats(name, arr):
    A = np.asarray(arr)
    finite = np.isfinite(A)
    n_nan = int(np.isnan(A).sum())
    n_inf = int(np.isinf(A).sum())
    mn = float(np.min(A)) if A.size else 0.0
    mx = float(np.max(A)) if A.size else 0.0
    l2 = float(np.linalg.norm(A)) if A.size else 0.0
    print(f"{name}: shape={A.shape} min={mn:.6e} max={mx:.6e} ||.||_2={l2:.6e} NaN={n_nan} Inf={n_inf} finite={finite.all()}")


def zero_dynamic_buffers(md):
    """Zero avel, aforce, cvf, fneigh for all systems."""
    z_vec = np.zeros(md.nSystems * md.nvecs * 4, dtype=np.float32)
    z_fng = np.zeros(md.nSystems * md.nnode * 8, dtype=np.float32)
    md.toGPU('aforce', z_vec)
    if 'aforce_old' in md.buffer_dict:
        md.toGPU('aforce_old', z_vec)
    md.toGPU('avel',   z_vec)
    md.toGPU('cvf',    z_vec)
    md.toGPU('fneigh', z_fng)
    z_tdrive = np.zeros(md.nSystems * 4, dtype=np.float32)
    md.toGPU('TDrives', z_tdrive)


def fetch_arrays(md, mm):
    """Download common arrays. Returns dict with atoms-only slices for convenience."""
    apos, aforce = md.download_results()
    avel  = md.download_buf('avel').reshape(md.nSystems, md.nvecs, 4)
    fng   = md.download_buf('fneigh').reshape(md.nSystems, md.nnode, 8)
    cvf   = md.download_buf('cvf').reshape(md.nSystems, md.nvecs, 4)
    # Atoms-only views from system 0
    atoms = {
        'apos_atoms':      apos[:mm.natoms, :3],
        'avel_atoms':      avel[0, :mm.natoms, :3],
        'afor_atoms':      aforce[:mm.natoms, :3],
        'afor_atoms_full': aforce[:mm.natoms, :4],
    }
    return {'apos': apos, 'aforce': aforce, 'avel': avel, 'fneigh': fng, 'cvf': cvf, **atoms}


def write_xyz_trajectory(path, trj_atoms, symbols):
    """Write trajectory frames to XYZ.

    Args:
        path: output file path (str or Path)
        trj_atoms: array of shape (steps, natoms, 3)
        symbols: iterable of length natoms with atomic symbols
    """

    path = Path(path)
    nsteps, natoms, _ = trj_atoms.shape
    if len(symbols) != natoms:  raise ValueError(f"symbols length {len(symbols)} does not match natoms {natoms}")

    with path.open('w') as fh:
        for step in range(nsteps):
            fh.write(f"{natoms}\n")
            fh.write(f"step {step}\n")
            coords = trj_atoms[step]
            for sym, (x, y, z) in zip(symbols, coords):
                fh.write(f"{sym} {x:.6f} {y:.6f} {z:.6f}\n")

def compute_totals(apos_atoms, avel_atoms, aforce_atoms, masses=None, origin=None):
    """
    Compute global conservation diagnostics for a single frame.
    Returns dict with keys: F, T, P, L, Rcm, Vcm, Fcm, Tcm.
    """
    r = np.asarray(apos_atoms, dtype=np.float32)
    v = np.asarray(avel_atoms, dtype=np.float32)
    f = np.asarray(aforce_atoms, dtype=np.float32)
    n = r.shape[0]
    if masses is None:
        m = np.ones((n, 1), dtype=np.float32)
    else:
        m = np.asarray(masses, dtype=np.float32).reshape(n, 1)

    mtot = float(m.sum())
    if mtot <= 0.0:
        raise ValueError("Total mass must be positive to compute COM diagnostics")

    rcm = (m * r).sum(axis=0) / mtot
    vcm = (m * v).sum(axis=0) / mtot

    Ftot = f.sum(axis=0)
    P    = (m * v).sum(axis=0)

    r0 = np.zeros(3, dtype=np.float32) if origin is None else np.asarray(origin, dtype=np.float32)
    rrel = r - r0
    Ttot = np.cross(rrel, f).sum(axis=0)
    L    = np.cross(rrel, (m * v)).sum(axis=0)

    Fcm = Ftot.copy()
    Tcm = np.cross(rcm - r0, Ftot)

    return {
        'F': Ftot,
        'T': Ttot,
        'P': P,
        'L': L,
        'Rcm': rcm,
        'Vcm': vcm,
        'Fcm': Fcm,
        'Tcm': Tcm,
    }


def print_totals(label, totals):
    F, T, P, L = totals['F'], totals['T'], totals['P'], totals['L']
    def fmt(x): return f"({x[0]: .3e}, {x[1]: .3e}, {x[2]: .3e})"
    extras = []
    for key in ('Ekin', 'Epot', 'Etotal'):
        if key in totals:
            extras.append(f"{key}={totals[key]:.6e}")
    extra_str = " " + " ".join(extras) if extras else ""
    print(f"[{label}] F={fmt(F)} | T={fmt(T)} | P={fmt(P)} | L={fmt(L)}{extra_str}")


def plot_trajectories(trj_atoms, dim='xy', labels=None, title='Trajectories', save_path=None, show=True):
    """
    Plot atom trajectories from trj_atoms with shape (nSteps, natoms, 3).
    dim: 'xy' | 'xz' | 'yz'
    """
    proj = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(dim, (0,1))
    ii, jj = proj
    nsteps, natoms, _ = trj_atoms.shape
    fig, ax = plt.subplots(figsize=(6,6))
    for a in range(natoms):
        xy = trj_atoms[:, a, :]
        ax.plot(xy[:,ii], xy[:,jj], '-', lw=1, label=(labels[a] if labels is not None else None), zorder=1)
        ax.plot(xy[-1,ii], xy[-1,jj], 'o', ms=3, zorder=2)
    ax.set_xlabel(['x','y','z'][ii])
    ax.set_ylabel(['x','y','z'][jj])
    ax.set_title(title)
    if labels is not None and len(labels) <= 12:
        ax.legend(loc='best', fontsize=8)
    ax.set_aspect('equal', adjustable='box')
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot to {save_path}")
    if show:
        plt.show()
    return fig, ax


def configure_md(mol2_path, dt, damp, flim, subtract_vdw, drive_temp, drive_gamma, drive_seed, builder=build_mmff_from_mol, print_params=False):
    _ENERGY_CACHE['prev_vel'] = None
    _ENERGY_CACHE['prev_pos'] = None
    mol, mm = builder(mol2_path)
    md = MolecularDynamics(nloc=32)
    md.set_pack_system_debug(print_params)
    md.realloc(mm, nSystems=1)
    mm.dt = dt
    mm.damp = damp
    mm.Flimit = flim
    md.upload_all_systems()
    zero_dynamic_buffers(md)
    md.setup_kernels()
    md.kernel_params['bSubtractVdW'] = np.int32(1 if subtract_vdw else 0)
    md.kernel_args_getMMFFf4 = md.generate_kernel_args('getMMFFf4')
    md.kernel_args_runMD     = md.generate_kernel_args('runMD')
    md_params = np.array([dt, damp, damp, 0.0], dtype=np.float32)
    float4_size = 4 * np.float32().itemsize
    for iSys in range(md.nSystems):
        md.toGPU('MDparams', md_params, byte_offset=iSys * float4_size)
    tdrive = np.array([drive_temp, drive_gamma, 0.0, drive_seed], dtype=np.float32)
    for iSys in range(md.nSystems):
        md.toGPU('TDrives', tdrive, byte_offset=iSys * 4 * np.float32().itemsize)
    return mol, mm, md


def init_observers(record, monitor, monitor_props):
    if monitor_props is None:
        monitor_props = list(DEFAULT_MONITOR_PROPS) if monitor else []
    monitor_enabled = monitor and len(monitor_props) > 0
    need_history = record or monitor_enabled
    totals_hist = [] if need_history else None
    monitor_data = {name: [] for name in monitor_props} if monitor_enabled else None
    trj = [] if record else []
    return monitor_props, monitor_enabled, totals_hist, monitor_data, trj


def collect_diagnostics(md, mm, masses, record, monitor_enabled, monitor_props, totals_hist, monitor_data, trj):
    buf = fetch_arrays(md, mm)
    apos_step = buf['apos_atoms'].copy()
    avel_step = buf['avel_atoms'].copy()
    afor_step = buf['afor_atoms'].copy()
    afor_full_step = buf['afor_atoms_full'].copy()
    totals = compute_totals(apos_step, avel_step, afor_step, masses=masses)
    pi_count = max(mm.nvecs - mm.natoms, 0)
    avel_pi = None
    if pi_count > 0:
        avel_pi = buf['avel'][0, mm.natoms:mm.natoms + pi_count, :3]
    totals.update(compute_energies(avel_step, masses, afor_full_step, avel_pi=avel_pi))
    if totals_hist is not None:
        totals_hist.append(totals)
    if record:
        trj.append(apos_step)
    if monitor_enabled and monitor_data is not None:
        for name in monitor_props:
            val = totals.get(name)
            if val is None:
                continue
            monitor_data[name].append(float(val) if np.isscalar(val) else np.array(val, dtype=np.float32))


def dump_buffers(buf, print_stats=False, print_arrays=False):
    if not (print_stats or print_arrays):
        return
    apos, aforce = buf['apos'], buf['aforce']
    avel, fng, cvf = buf['avel'], buf['fneigh'], buf['cvf']
    if print_stats:
        stats('apos', apos)
        stats('avel', avel)
        stats('aforce', aforce)
        stats('fneigh', fng)
        stats('cvf', cvf)
    if print_arrays:
        print('apos:\n', apos)
        print('avel:\n', avel)
        print('aforce:\n', aforce)
        print('fneigh:\n', fng)
        print('cvf:\n', cvf)


def report_conservation(totals_hist, steps, record):
    if not totals_hist:
        return
    F = np.array([x['F'] for x in totals_hist])
    T = np.array([x['T'] for x in totals_hist])
    P = np.array([x['P'] for x in totals_hist])
    L = np.array([x['L'] for x in totals_hist])
    def dmax(A):
        d = A - A[0]
        return float(np.max(np.linalg.norm(d, axis=1)))
    energy_drifts = []
    for key in ('Ekin', 'Epot', 'Etotal'):
        if key in totals_hist[0]:
            arr = np.array([x[key] for x in totals_hist], dtype=np.float64)
            energy_drifts.append(f"|Δ{key}|max={np.max(np.abs(arr - arr[0])):.3e}")
    drift_msg = (
        f"Conservation drift over {steps} steps: "
        f"|ΔF|max={dmax(F):.3e} |ΔT|max={dmax(T):.3e} "
        f"|ΔP|max={dmax(P):.3e} |ΔL|max={dmax(L):.3e}"
    )
    if energy_drifts:
        drift_msg += " " + " ".join(energy_drifts)
    print(drift_msg)
    if record:
        print_totals('totals[0]', totals_hist[0])
        print_totals('totals[-1]', totals_hist[-1])


def overlay_bonds(mol, frame, plot_dim, labels=None):
    from .plotUtils import plotAtoms, plotBonds
    if getattr(mol, 'bonds', None) is None or len(mol.bonds) == 0:
        mol.findBonds(Rcut=3.0, RvdwCut=0.5)
    bonds_arr = np.asarray(mol.bonds, dtype=np.int32)
    if bonds_arr.ndim == 1:
        bonds_arr = bonds_arr.reshape(1, -1)
    links = bonds_arr[:, :2]
    axes = {'xy': (0,1), 'xz': (0,2), 'yz': (1,2)}.get(plot_dim, (0,1))
    plotAtoms(apos=frame, es=getattr(mol, 'enames', None), axes=axes, sizes=60., colors='#404040', marker='o', labels=labels, bNumbers=False)
    plotBonds(links=links, ps=frame, axes=axes, colors='k', lws=1.5)


# def finalize_recording(record, trj, mol, save_trj, plot, plot_dim, plot_label_mode, plot_bonds, plot_title=None):
#     if not (record and trj): return None, None
#     trj_arr = np.stack(trj, axis=0)
#     if save_trj:
#         np.save(save_trj, trj_arr)
#         print(f"Saved trajectory to {save_trj} with shape {trj_arr.shape}")
#     if not plot:  return trj_arr, None
#     labels = None
#     if plot_label_mode == 'number':
#         labels = [str(i) for i in range(trj_arr.shape[1])]
#     elif plot_label_mode == 'type' and getattr(mol, 'enames', None) is not None:
#         labels = mol.enames
#     fig, _ = plot_trajectories(trj_arr, dim=plot_dim, labels=labels, title=f'Trajectories ({plot_dim})', save_path=None, show=False)
#     if fig is not None and plot_title:  fig.suptitle(plot_title, fontsize=12)
#     if plot_bonds:  overlay_bonds(mol, trj_arr[-1], plot_dim, labels=labels)
#     return trj_arr, fig


def finalize_monitoring(monitor_enabled, monitor_data, monitor_props, monitor_save_data, monitor_plot, monitor_save):
    if not (monitor_enabled and monitor_data):  return None
    series = finalize_monitor_series(monitor_data)
    summarize_monitor_series(series, monitor_props)
    if monitor_save_data:
        np.savez(monitor_save_data, **series)
        print(f"Saved monitor data to {monitor_save_data}")
    if monitor_plot or monitor_save:
        fig, _ = plot_monitor_series(series, monitor_props, show=not monitor_save)
        if monitor_save and fig is not None:
            fig.savefig(monitor_save, dpi=150)
            print(f"Saved monitor plot to {monitor_save}")
    return series


def displace_atom_positions(apos, atom_index, displacement):
    apos_new = np.array(apos, copy=True)
    disp = np.zeros_like(apos_new[atom_index])
    disp[:3] = displacement
    apos_new[atom_index] += disp
    return apos_new


def centered_finite_difference(values, dx):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size < 3:
        raise ValueError("Need at least three samples for centered finite difference")
    deriv = (values[2:] - values[:-2]) / (2.0 * dx)
    midpoints = np.arange(1, values.size - 1)
    return deriv, midpoints


def upload_mmff_positions(mm, md, positions):
    positions = np.asarray(positions, dtype=np.float32)
    if positions.shape != mm.apos.shape:
        raise ValueError(f"upload_mmff_positions: expected shape {mm.apos.shape}, got {positions.shape}")
    mm.apos[:, :3] = positions[:, :3]
    if positions.shape[1] > 3:
        mm.apos[:, 3:] = positions[:, 3:]
    md.toGPU('apos', mm.apos.astype(np.float32).ravel())


def evaluate_mmff_gpu(mm, md, positions, do_clean=True, do_nb=False, do_mmff=True, mode='none'):
    upload_mmff_positions(mm, md, positions)
    run_step(md, do_clean=do_clean, do_nb=do_nb, do_mmff=do_mmff, mode=mode)
    if mode == 'none':
        orig_mdparams = md.download_buf('MDparams').copy().reshape(md.nSystems, 4)
        zero_mdparams = np.zeros(4, dtype=np.float32)
        float4_size = 4 * np.float32().itemsize
        for iSys in range(md.nSystems):
            md.toGPU('MDparams', zero_mdparams, byte_offset=iSys * float4_size)
        md.run_updateAtomsMMFFf4()
        for iSys in range(md.nSystems):
            md.toGPU('MDparams', orig_mdparams[iSys], byte_offset=iSys * float4_size)
    buf = fetch_arrays(md, mm)
    forces = buf['afor_atoms'].astype(np.float64)
    energy = float(buf['afor_atoms_full'][:, 3].sum())
    return energy, forces


def scan_energy_force(mm, md, atom_index=0, axis=0, dx=1e-2, nsamp=100, restore=True, evaluator=None, do_nb=False):
    if nsamp < 3 or nsamp % 2 == 0:
        raise ValueError("nsamp must be odd and >= 3 for centered differences")
    if evaluator is None:
        evaluator = lambda apos: evaluate_mmff_gpu(mm, md, apos, do_clean=True, do_nb=do_nb, do_mmff=True, mode='none')
    base_apos = np.array(mm.apos, copy=True)
    natoms = mm.natoms
    energies = np.zeros(nsamp, dtype=np.float64)
    forces = np.zeros((nsamp, natoms, 3), dtype=np.float64)
    offsets = (np.arange(nsamp) - nsamp // 2) * dx
    disp = np.zeros(3, dtype=np.float64)
    for i, offset in enumerate(offsets):
        disp[:] = 0.0
        disp[axis] = offset
        apos_shifted = displace_atom_positions(base_apos, atom_index, disp)
        energy, force = evaluator(apos_shifted)
        energies[i] = energy
        if force.shape != (natoms, 3):
            raise ValueError(f"Evaluator returned forces with shape {force.shape}, expected {(natoms, 3)}")
        forces[i] = force
    if restore:
        upload_mmff_positions(mm, md, base_apos)
    analytic_force = -forces[:, atom_index, axis]
    numeric_force, idx = centered_finite_difference(energies, dx)
    mid_offsets = offsets[idx]
    diff = analytic_force[idx] - numeric_force
    stats = {
        'energy_min': float(energies.min()),
        'energy_max': float(energies.max()),
        'force_min': float(analytic_force.min()),
        'force_max': float(analytic_force.max()),
        'diff_min': float(diff.min()),
        'diff_max': float(diff.max()),
        'diff_rms': float(np.sqrt(np.mean(diff ** 2))),
    }
    return {
        'offsets': offsets,
        'energies': energies,
        'forces': forces,
        'analytic_force': analytic_force,
        'numeric_force': numeric_force,
        'force_mid_offsets': mid_offsets,
        'diff': diff,
        'diff_stats': stats,
    }


def plot_energy_force_scan(scan, axis_label='displacement [$\AA$]', show=True, save_path=None):
    offsets = scan['offsets']
    energies = scan['energies']
    mid_offsets = scan['force_mid_offsets']
    analytic_mid = scan['analytic_force'][1:-1]
    numeric = scan['numeric_force']
    diff = scan['diff']
    stats = scan['diff_stats']
    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    axE, axF, axD = axes
    axE.plot(offsets, energies, '-', lw=0.5)
    axE.set_ylabel('Energy')
    axE.grid(alpha=0.3)
    axF.plot(mid_offsets, analytic_mid, '-', lw=0.5, label='Analytic')
    axF.plot(mid_offsets, numeric,      ':', lw=1.5, label='Numeric')
    axF.set_ylabel('Force component')
    axF.grid(alpha=0.3)
    axF.legend(loc='best')
    axD.plot(mid_offsets, diff, '-', lw=0.5, label='Analytic - Numeric')
    axD.axhline(0.0, color='k', lw=0.8)
    axD.set_xlabel(axis_label)
    axD.set_ylabel('Force diff')
    axD.grid(alpha=0.3)
    axD.text(0.05, 0.95,
             f"min={stats['diff_min']:.3e}\nmax={stats['diff_max']:.3e}\nrms={stats['diff_rms']:.3e}",
             transform=axD.transAxes, va='top', ha='left', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"Saved scan plot to {save_path}")
    if show:
        plt.show()
    return fig, axes