import numpy as np
import matplotlib
import os
if os.environ.get('DISPLAY', '') == '':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
 

import pyopencl as cl

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
    'Ekin_rot': 'scalar',
    'Epot': 'scalar',
    'Etotal': 'scalar',
}

DEFAULT_MONITOR_PROPS = ['F', 'P', 'L', 'Rcm', 'Vcm', 'Ekin', 'Ekin_rot', 'Epot', 'Etotal']

_ENERGY_CACHE = {
    'prev_vel': None,
    'prev_pos': None,
}


def clean_element_symbol(s):
    if s is None:
        return ''
    ss = str(s)
    a = ''.join(ch for ch in ss if ch.isalpha())
    if a:
        if len(a) >= 2 and a[0].isupper() and a[1].isupper():
            a = a[0] + a[1].lower() + a[2:]
        return a
    if ss:
        return ss[0]
    return ''


def clean_element_symbols_list(enames):
    return [clean_element_symbol(s) for s in list(enames)]


def apply_random_displacement(pos, dmax, seed=None):
    if dmax is None:
        return pos
    dmax = float(dmax)
    if dmax <= 0.0:
        return pos
    rng = np.random.default_rng(seed)
    d = rng.uniform(-dmax, dmax, size=pos.shape).astype(np.float32)
    pos[:] = (np.asarray(pos, dtype=np.float32) + d)
    return pos


def max_norm3(A):
    A = np.asarray(A)
    return float(np.sqrt(np.max(np.sum(A*A, axis=1)))) if A.size else 0.0


def write_xyz_frame(fout, enames, pos3, comment=""):
    pos3 = np.asarray(pos3)
    nat = pos3.shape[0]
    fout.write(f"{nat}\n")
    fout.write(f"{comment}\n")
    for i in range(nat):
        x, y, z = float(pos3[i, 0]), float(pos3[i, 1]), float(pos3[i, 2])
        fout.write(f"{enames[i]} {x: .10f} {y: .10f} {z: .10f}\n")


def iter_should_write(istep, write_stride, steps_total):
    if istep == 0:
        return True
    if istep == steps_total:
        return True
    w = int(write_stride)
    if w <= 0:
        return False
    return (istep % w) == 0


def plot_fmax_history_png(png_path, steps, fmax_hist, fconv=None, title=None):
    png_path = str(png_path)
    if not png_path:
        return
    steps = np.asarray(steps, dtype=np.int32)
    fmax_hist = np.asarray(fmax_hist, dtype=np.float64)
    if steps.size == 0 or fmax_hist.size == 0:
        return
    if not np.isfinite(fmax_hist).all():
        return
    if os.environ.get('DISPLAY', '') == '':
        try:
            import matplotlib
            matplotlib.use('Agg')
        except Exception:
            pass
    fig = plt.figure(figsize=(6.0, 3.6))
    ax = fig.add_subplot(1, 1, 1)
    ax.semilogy(steps, fmax_hist, '-', lw=1.5, label='max|F|')
    if fconv is not None:
        ax.axhline(float(fconv), color='k', lw=1.0, ls='--', label='fconv')
    ax.set_xlabel('step')
    ax.set_ylabel('max|F|')
    if title:
        ax.set_title(title)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)


def conv_plot_path_from_xyz(out_xyz_path):
    s = str(out_xyz_path)
    if s.lower().endswith('.xyz'):
        return s[:-4] + '.conv.png'
    return s + '.conv.png'


def format_relax_summary(ff_name, mol_path, converged, nstep, fmax, out_xyz, out_png, wall_s, abort_reason=None):
    m = os.path.basename(str(mol_path))
    c = 1 if converged else 0
    msg = f"RELAX_SUMMARY ff={ff_name} mol={m} converged={c} nstep={int(nstep)} fmax={float(fmax):.6e} wall={float(wall_s):.3f}s xyz={out_xyz} png={out_png}"
    if abort_reason:
        msg += f" abort={abort_reason}"
    return msg


def effective_verbosity(quiet=0, verbosity=1):
    v = int(verbosity)
    if int(quiet):
        v = 0
    return v

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


def compute_energies(avel_all, natoms, masses, aforce_atoms_full):
    vel_atoms = np.asarray(avel_all[:natoms, :3], dtype=np.float32)
    v = _estimate_true_velocity(vel_atoms, masses)   # use this with Leap-Frong integrator
    #v = np.asarray(avel_atoms, dtype=np.float32)       # use this with Verlet integrator
    m = np.asarray(masses, dtype=np.float32).reshape(-1, 1)
    kin = 0.5 * m * (v * v).sum(axis=1, keepdims=True)
    Ekin = float(kin.sum())
    omega = np.asarray(avel_all[natoms:, :3], dtype=np.float32)
    if omega.size:
        omega_sq = (omega * omega).sum(axis=1)
        # if np.allclose(omega_sq, 0.0):
        #     print(f"DEBUG Ekin_rot zero: max |omega|={float(np.max(np.abs(omega))):.3e} n_pi={omega.shape[0]}")
        # else:
        #     print(f"DEBUG Ekin_rot stats: max |omega|={float(np.max(np.sqrt(omega_sq))):.3e} min |omega|={float(np.min(np.sqrt(omega_sq))):.3e} sum |omega|^2={float(np.sum(omega_sq)):.3e}")
        Ekin_rot = 0.5 * float(np.sum(omega_sq))
        Ekin += Ekin_rot
    else:
        Ekin_rot = 0.0
    af = np.asarray(aforce_atoms_full, dtype=np.float32)
    # Each kernel assigns interaction energy to both partners; divide by 2 to avoid double counting.
    Epot = float(af[:, 3].sum())
    return {
        'Ekin': Ekin,
        'Ekin_rot': Ekin_rot,
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
    if not props:
        return None, None

    display_props = []
    for name in props:
        if name == 'Ekin_rot' and 'Ekin' in props:
            continue
        display_props.append(name)

    fig, axes = plt.subplots(len(display_props), 1, figsize=(8, 2.6 * len(display_props)), sharex=True)
    if isinstance(axes, np.ndarray):
        axes_list = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        axes_list = list(axes)
    else:
        axes_list = [axes]

    axes_map = {name: ax for name, ax in zip(display_props, axes_list)}
    if 'Ekin' in axes_map and 'Ekin_rot' in props:
        axes_map['Ekin_rot'] = axes_map['Ekin']

    steps = None
    used_axes = set()
    for name in props:
        ax = axes_map.get(name)
        if ax is None:
            continue
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
        if name == 'Ekin_rot' and 'Ekin' in props:
            ax.plot(steps, data, label='Ekin_rot', linestyle='--')
        if ax not in used_axes:
            ax.set_ylabel(name)
            ax.grid(True, alpha=0.3)
            used_axes.add(ax)
        if name in {'Ekin', 'Ekin_rot'}:
            ax.legend(loc='best', fontsize=8)

    if axes_list:
        axes_list[-1].set_xlabel('step')
    fig.tight_layout()
    if show:
        plt.show()
    return fig, axes


def build_mmff_from_mol(mol2_path):
    from .AtomicSystem import AtomicSystem
    mol = AtomicSystem(fname=mol2_path)
    # Ensure neighbors present
    if getattr(mol, 'ngs', None) is None:
        mol.neighs()
    mm = MMFF(bTorsion=False, verbosity=1)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    return mol, mm


def build_2node_system(params_dict=None):
    """Create minimal 2-atom system with pi orbitals matching pipi_dynamics.py.
    
    Args:
        params_dict: optional dict or str (ignored if str to match mol2_path signature).
                     If dict, may contain 'bL0', 'ha', 'hb' keys.
    
    Returns:
        mol: AtomicSystem with 2 carbon atoms
        mm: MMFF object with force field populated
    """
    from .AtomicSystem import AtomicSystem
    
    if isinstance(params_dict, str):
        params_dict = None
    params = params_dict if params_dict is not None else {}
    
    bL0 = params.get('bL0', 1.5)
    kBL = np.float32(params.get('kBL', 10.0))
    kSP = np.float32(params.get('kSP', 1.0))
    kPP = np.float32(params.get('kPP', 1.0))
    ha_init = params.get('ha', None)
    hb_init = params.get('hb', None)
    
    if ha_init is None:
        ha_init = np.array([0.1, 1.0, 0.15], dtype=np.float32)
        ha_init /= np.linalg.norm(ha_init)
    if hb_init is None:
        hb_init = np.array([0.2, 1.0, -0.10], dtype=np.float32)
        hb_init /= np.linalg.norm(hb_init)
    
    apos = np.array([[0.0, 0.0, 0.0], [bL0, 0.0, 0.0]], dtype=np.float64)
    enames = ['C', 'C']
    atypes = np.array([6, 6], dtype=np.int32)
    qs = np.array([0.0, 0.0], dtype=np.float64)
    bonds = [(0, 1, 1.0)]
    
    mol = AtomicSystem(apos=apos, atypes=atypes, enames=enames, qs=qs, bonds=bonds)
    mol.neighs()
    
    mm = MMFF(bTorsion=False, verbosity=1)
    mm.toMMFFsp3_loc(mol=mol, atom_types=mm.atom_types, bRealloc=True, bEPairs=False, bUFF=False)
    
    if mm.nnode == 2:
        mm.pipos[0, :] = ha_init.astype(np.float32)
        mm.pipos[1, :] = hb_init.astype(np.float32)
        mm.apos[mm.natoms, :3] = ha_init.astype(np.float32)
        mm.apos[mm.natoms+1, :3] = hb_init.astype(np.float32)
    
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
    avel  = md.download_buf('avel')  .reshape(md.nSystems, md.nvecs, 4)
    fng   = md.download_buf('fneigh', dtype=np.float32).reshape(md.nSystems, md.nnode, 8)
    cvf   = md.download_buf('cvf')   .reshape(md.nSystems, md.nvecs, 4)
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

def compute_totals(apos_atoms, avel_atoms, aforce_atoms, masses=None, origin=None, pi_omega=None, pi_inertia=None):
    """
    Compute global conservation diagnostics for a single frame.
    Returns dict with keys: F, T, P, L, Rcm, Vcm, Fcm, Tcm.

    `pi_omega` can be provided to include rotational angular momentum from
    pi-orbital degrees of freedom. If `pi_inertia` is None a unit inertia is
    assumed, otherwise it can be scalar or array broadcastable to `pi_omega`.
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
    L_atoms = np.cross(rrel, (m * v)).sum(axis=0)

    L_rot = np.zeros(3, dtype=np.float32)
    if pi_omega is not None:
        omega = np.asarray(pi_omega, dtype=np.float32)
        if omega.ndim != 2 or omega.shape[1] != 3:
            raise ValueError("pi_omega must have shape (n_pi, 3)")
        if pi_inertia is None:
            L_rot = omega.sum(axis=0)
        else:
            I = np.asarray(pi_inertia, dtype=np.float32)
            if I.ndim == 0:
                L_rot = (I * omega).sum(axis=0)
            else:
                I = I.reshape(-1, 1)
                if I.shape[0] != omega.shape[0]:
                    raise ValueError("pi_inertia length does not match pi_omega")
                L_rot = (I * omega).sum(axis=0)
    L = L_atoms + L_rot

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
    avel_all = buf['avel'][0, :, :]
    avel_step = avel_all[:mm.natoms, :3].copy()
    afor_step = buf['afor_atoms'].copy()
    afor_full_step = buf['afor_atoms_full'].copy()
    pi_slice = slice(mm.natoms, mm.natoms + mm.nnode)
    pi_omega = None
    if mm.nnode > 0:
        pi_omega = avel_all[pi_slice, :3].copy()
    pi_inertia = getattr(mm, 'pi_inertia', None)
    totals = compute_totals(apos_step, avel_step, afor_step, masses=masses, pi_omega=pi_omega, pi_inertia=pi_inertia)
    totals.update(compute_energies(avel_all, mm.natoms, masses, afor_full_step))
    if totals_hist is not None:
        totals_hist.append(totals)
    if record:
        if mm.nnode > 0:
            host_nodes = apos_step[:mm.nnode, :]
            pi_vecs = buf['apos'][mm.natoms:mm.natoms + mm.nnode, :3].copy()
            norms = np.linalg.norm(pi_vecs, axis=1, keepdims=True)
            norms[norms == 0.0] = 1.0
            pi_offsets = pi_vecs / norms
            pi_positions = host_nodes + pi_offsets
            frame = np.concatenate((apos_step, pi_positions.astype(np.float32)), axis=0)
        else:
            frame = apos_step
        trj.append(frame)
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


def run_step(md, do_clean=True, do_nb=False, do_mmff=True, mode='none', use_rot_force=False):
    if do_clean:
        md.run_cleanForceMMFFf4()
    if do_nb:
        md.run_getNonBond()
    if do_mmff:
        if use_rot_force:
            md.run_getMMFFf4_rot()
        else:
            md.run_getMMFFf4()
    if mode == 'basic':
        md.run_updateAtomsMMFFf4()
    elif mode == 'rot':
        md.run_updateAtomsMMFFf4_rot()
    elif mode == 'none':
        pass
    else:
        raise ValueError(f"run_step: unknown mode '{mode}'")


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


def mmff_cpp_init_from_atoms(apos, enames, data_dir, nPBC=(0,0,0), b141=True, bSimple=False, bConj=True, bCumulene=True, tmp_xyz="/tmp/mmff_parity_tmp.xyz"):
    from . import MMFF as mmff_cpp
    apos = np.asarray(apos, dtype=np.float64)
    if apos.ndim != 2 or apos.shape[1] < 3:
        raise ValueError(f"mmff_cpp_init_from_atoms: expected apos (natoms,>=3), got {apos.shape}")
    natoms = apos.shape[0]
    if len(enames) != natoms:
        raise ValueError(f"mmff_cpp_init_from_atoms: len(enames)={len(enames)} != natoms={natoms}")
    # write temporary xyz for C++ loader
    with open(tmp_xyz, 'w') as f:
        f.write(f"{natoms}\n")
        f.write("mmff_parity\n")
        for i in range(natoms):
            x,y,z = apos[i,:3]
            f.write(f"{enames[i]} {x:.16e} {y:.16e} {z:.16e}\n")
    sElementTypes  = os.path.join(data_dir, "ElementTypes.dat")
    sAtomTypes     = os.path.join(data_dir, "AtomTypes.dat")
    sBondTypes     = os.path.join(data_dir, "BondTypes.dat")
    sAngleTypes    = os.path.join(data_dir, "AngleTypes.dat")
    sDihedralTypes = os.path.join(data_dir, "DihedralTypes.dat")
    mmff_cpp.init(
        xyz_name=tmp_xyz,
        surf_name=None,
        smile_name=None,
        sElementTypes=sElementTypes,
        sAtomTypes=sAtomTypes,
        sBondTypes=sBondTypes,
        sAngleTypes=sAngleTypes,
        sDihedralTypes=sDihedralTypes,
        bMMFF=True,
        bEpairs=False,
        nPBC=nPBC,
        gridStep=0.1,
        bUFF=False,
        b141=b141,
        bSimple=bSimple,
        bConj=bConj,
        bCumulene=bCumulene,
    )
    mmff_cpp.getBuffs()
    return {
        'tmp_xyz': tmp_xyz,
        'data_dir': data_dir,
    }


def mmff_cpp_set_switches(do_nonbond=False, do_angles=True, do_pisigma=True, do_pipii=True, do_pbc=False, do_check_invariants=False):
    from . import MMFF as mmff_cpp
    def sw(b):
        return 1 if b else -1
    # explicitly force (0 means keep)
    mmff_cpp.setSwitches(
        CheckInvariants=sw(do_check_invariants),
        PBC=sw(do_pbc),
        NonBonded=sw(do_nonbond),
        NonBondNeighs=sw(do_nonbond),
        SurfAtoms=-1,
        GridFF=-1,
        MMFF=1,
        Angles=sw(do_angles),
        PiSigma=sw(do_pisigma),
        PiPiI=sw(do_pipii),
    )


def evaluate_mmff_cpp(positions, do_nonbond=False, do_angles=True, do_pisigma=True, do_pipii=True, do_pbc=False):
    from . import MMFF as mmff_cpp
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] < 3:
        raise ValueError(f"evaluate_mmff_cpp: expected positions (natoms,>=3), got {positions.shape}")
    # set switches first (avoid 0=keep trap)
    mmff_cpp_set_switches(do_nonbond=do_nonbond, do_angles=do_angles, do_pisigma=do_pisigma, do_pipii=do_pipii, do_pbc=do_pbc, do_check_invariants=False)
    # push geometry into C++ buffers
    mmff_cpp.apos[:,:3] = positions[:,:3]
    E = float(mmff_cpp.eval())
    F = np.array(mmff_cpp.fapos[:,:3], copy=True, dtype=np.float64)
    return E, F


def force_parity_report(F_ref, F_tst, tol=5e-5, label_ref="REF", label_tst="TST"):
    F_ref = np.asarray(F_ref, dtype=np.float64)
    F_tst = np.asarray(F_tst, dtype=np.float64)
    if F_ref.shape != F_tst.shape:
        raise ValueError(f"force_parity_report: shape mismatch {F_ref.shape} vs {F_tst.shape}")
    dF = F_tst - F_ref
    max_abs = float(np.max(np.abs(dF))) if dF.size else 0.0
    rms = float(np.sqrt(np.mean(dF*dF))) if dF.size else 0.0
    idx = np.unravel_index(int(np.argmax(np.abs(dF))), dF.shape) if dF.size else (0,0)
    ok = max_abs <= tol
    msg = {
        'ok': ok,
        'max_abs': max_abs,
        'rms': rms,
        'idx': idx,
        'ref': F_ref[idx[0]],
        'tst': F_tst[idx[0]],
        'diff': dF[idx[0]],
        'tol': tol,
        'label_ref': label_ref,
        'label_tst': label_tst,
    }
    return msg


def evaluate_uff_gpu(uff_cl, positions, bClearForce=True):
    positions = np.asarray(positions, dtype=np.float32)
    if positions.ndim != 2 or positions.shape[1] < 3:
        raise ValueError(f"evaluate_uff_gpu: expected positions shape (natoms,>=3), got {positions.shape}")
    if positions.shape[0] != uff_cl.natoms:
        raise ValueError(f"evaluate_uff_gpu: expected natoms={uff_cl.natoms}, got {positions.shape[0]}")
    uff_cl.upload_positions(positions)
    uff_cl.run_eval_step(bClearForce=bClearForce)
    f4 = np.zeros(uff_cl.natoms * 4, dtype=np.float32)
    cl.enqueue_copy(uff_cl.queue, f4, uff_cl.buffer_dict['fapos'])
    f4 = f4.reshape(uff_cl.natoms, 4)
    forces = f4[:, :3].astype(np.float64)
    energy = float(f4[:, 3].sum())
    return energy, forces


def scan_energy_force_uff(uff_cl, atom_index=0, axis=0, dx=1e-2, nsamp=100, restore=True, evaluator=None, base_apos=None):
    if nsamp < 3 or nsamp % 2 == 0:
        raise ValueError("nsamp must be odd and >= 3 for centered differences")
    if evaluator is None:
        evaluator = lambda apos: evaluate_uff_gpu(uff_cl, apos, bClearForce=True)
    if base_apos is None:
        raise ValueError("scan_energy_force_uff: base_apos must be provided (natoms,>=3)")
    base_apos = np.asarray(base_apos, dtype=np.float32)
    natoms = base_apos.shape[0]
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
        uff_cl.upload_positions(base_apos)
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