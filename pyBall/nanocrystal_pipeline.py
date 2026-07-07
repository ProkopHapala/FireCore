"""nanocrystal_pipeline.py — NPZ-stage orchestration: MMFF relax, topology-linear Hessian, eigh spectrum.

Design: per-crystal cache files 01→05 (see doc/Topics/FTIR_Nanocrystals/NPZ_Crystal_Schema.md).
Topology/FF params in 03_topology.npz are fixed at export; relax updates coordinates only.
Hessian default (--source topology) assembles exported springs at 02_relaxed.pos.

Open issues:
- Relax uses C++ MMFF; Hessian uses MMFFL linear sticks — full FF parity at relax pending LFF path.
- Eigenvectors not written to 05_spectrum.npz (v1.2); viewers use solve_normal_modes_from_hessian_npz.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]

_EV_TO_J = 1.602176634e-19
_AMU_TO_KG = 1.66053906660e-27
_ANG_TO_M = 1e-10
_C_CM_S = 2.99792458e10


def omega_internal_to_cm1(omega_internal: np.ndarray) -> np.ndarray:
    """sqrt(eV/amu)/Å (MMFF internal) → wavenumber cm⁻¹."""
    om = np.asarray(omega_internal, dtype=np.float64)
    omega2 = om * om
    omega_si = np.sqrt(np.maximum(omega2, 0.0) * _EV_TO_J / (_AMU_TO_KG * _ANG_TO_M ** 2))
    return omega_si / _C_CM_S


def filter_vibrational_modes(omegas: np.ndarray, rigid_shift: float = 1e6, vib_floor: float = 1e-4) -> np.ndarray:
    """Drop rigid-mode shift artifacts (ω ≈ sqrt(rigid_shift)) and near-zero modes."""
    om = np.asarray(omegas, dtype=np.float64)
    om = om[np.isfinite(om) & (om > vib_floor)]
    if om.size == 0:
        raise ValueError("filter_vibrational_modes: no modes left")
    cap = 0.5 * np.sqrt(float(rigid_shift))
    phys = om[om < cap]
    if phys.size == 0:
        raise ValueError(f"filter_vibrational_modes: all modes above cap {cap:.3g} (rigid_shift={rigid_shift})")
    return phys


def ensure_mmff_data() -> Path:
    tmmff = REPO / "tests/tMMFF"
    data = tmmff / "data"
    if not data.exists() and not data.is_symlink():
        os.symlink("../../cpp/common_resources", data, target_is_directory=True)
    return tmmff


def write_xyz(path: Path, pos: np.ndarray, Z: np.ndarray) -> None:
    sym_map = {1: "H", 6: "C", 14: "Si"}
    n = len(Z)
    lines = [str(n), "nanocrystal pipeline"]
    for i in range(n):
        s = sym_map.get(int(Z[i]), "C")
        lines.append(f"{s} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def init_path_from_crystal_npz(npz_path: Path, work_dir: Path) -> tuple[Path, dict]:
    """Load 01_crystal.npz; write temp mol2 (if bonds_ij) or xyz for MMFF.init."""
    from pyBall.io.crystal_npz import load_crystal_npz, enames_from_Z
    from pyBall.atomicUtils import save_mol2

    npz_path = Path(npz_path).resolve()
    crystal = load_crystal_npz(npz_path)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    bonds_ij = crystal.get("bonds_ij")
    if bonds_ij is not None and len(bonds_ij) > 0:
        mol2 = work_dir / "_init_from_npz.mol2"
        save_mol2(str(mol2), enames_from_Z(crystal["Z"]), crystal["pos"], bonds_ij, comment=f"from {npz_path.name}")
        return mol2, crystal
    xyz = work_dir / "_init_from_npz.xyz"
    write_xyz(xyz, crystal["pos"], crystal["Z"])
    return xyz, crystal


def write_status(out_npz: Path, status: dict) -> None:
    out_npz.with_suffix(".status.json").write_text(json.dumps(status))


def atom_Z_from_mmff(MMFF, natoms: int) -> np.ndarray:
    from pyBall import elements

    Z = np.zeros(natoms, dtype=np.int32)
    for i in range(natoms):
        tn = MMFF.getTypeName(i)
        sym = tn.split("_")[0] or "C"
        Z[i] = int(elements.ELEMENT_DICT[sym][0]) if sym in elements.ELEMENT_DICT else 6
    return Z


def relax(
    init_path: Path,
    out_npz: Path,
    out_xyz: Path | None = None,
    nstep_max: int = 2000,
    fconv: float = 1e-4,
    allow_unconverged: bool = False,
    bonds_ij: np.ndarray | None = None,
    Z_init: np.ndarray | None = None,
) -> dict:
    init_path = Path(init_path).resolve()
    if not init_path.is_file():
        raise FileNotFoundError(init_path)
    tmmff = ensure_mmff_data()
    orig = os.getcwd()
    t0 = time.perf_counter()
    natoms = 0
    pos = None
    Z = None
    fmax = float("nan")
    n_steps = 0
    converged = False
    try:
        os.chdir(tmmff)
        from pyBall import MMFF

        MMFF.setVerbosity(0, 0)
        MMFF.init(xyz_name=str(init_path), bEpairs=False, bMMFF=True, nPBC=(0, 0, 0))
        MMFF.getBuffs()
        natoms = MMFF.natoms
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
        MMFF.set_opt(dt_max=0.05, dt_min=0.01, damp_max=0.1, cvf_min=-0.1, cvf_max=+0.1)
        n_steps = int(MMFF.run(nstepMax=nstep_max, dt=0.05, Fconv=fconv, damping=0.1, ialg=2, omp=False))
        pos = MMFF.apos[:natoms].copy()
        fapos = MMFF.fapos[:natoms].copy()
        fmax = float(np.max(np.abs(fapos)))
        converged = fmax < fconv
        if not converged and not allow_unconverged:
            raise RuntimeError(f"MMFF relax did not converge: fmax={fmax:.3e} after {n_steps} steps")
        Z = np.asarray(Z_init, dtype=np.int32).reshape(-1) if Z_init is not None else atom_Z_from_mmff(MMFF, natoms)
        if Z.shape[0] != natoms:
            raise ValueError(f"relax: Z_init length {Z.shape[0]} != natoms {natoms}")
    finally:
        os.chdir(orig)
    if pos is None:
        raise RuntimeError("MMFF relax failed before producing coordinates")
    timing_ms = (time.perf_counter() - t0) * 1000.0
    out_npz = Path(out_npz)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    status = {"status": "ok", "natoms": int(natoms), "fmax": fmax, "converged": converged, "n_steps": n_steps, "timing_ms": timing_ms}
    save_kw = {"pos": pos, "Z": Z, "fmax": fmax, "n_steps": n_steps, "converged": converged, "timing_ms": timing_ms, "natoms": natoms}
    if bonds_ij is not None and len(bonds_ij) > 0:
        save_kw["bonds_ij"] = np.asarray(bonds_ij, dtype=np.int32)
    np.savez(out_npz, **save_kw)
    if out_xyz:
        write_xyz(Path(out_xyz), pos, Z)
    write_status(out_npz, status)
    return status


def _mmff_hessian_fd(init_path: Path, dx: float = 1e-4) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """MMFF finite-difference Hessian at init_path (mol2/xyz); same switches as relax."""
    init_path = Path(init_path).resolve()
    if not init_path.is_file():
        raise FileNotFoundError(init_path)
    tmmff = ensure_mmff_data()
    orig = os.getcwd()
    try:
        os.chdir(tmmff)
        from pyBall import MMFF, FTIR

        MMFF.setVerbosity(0, 0)
        MMFF.init(xyz_name=str(init_path), bEpairs=False, bMMFF=True, nPBC=(0, 0, 0))
        MMFF.getBuffs()
        natoms = MMFF.natoms
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
        inds = np.arange(natoms, dtype=np.int32)
        H = MMFF.getHessian3Nx3N(inds, dx=dx)
        H = 0.5 * (H + H.T)
        if np.isnan(H).any() or np.isinf(H).any():
            raise ValueError("NaN/Inf in MMFF Hessian")
        M = FTIR.get_mass_matrix(MMFF, natoms)
        pos = MMFF.apos[:natoms].copy()
        Z = atom_Z_from_mmff(MMFF, natoms)
        return H, M, pos, Z, natoms
    finally:
        os.chdir(orig)


def _topology_surface_arrays(topology: dict) -> dict:
    """Surface / nonbond group arrays from 03_topology.npz for cache in 04_hessian.npz."""
    out = {}
    for k in ('group_bbox_min', 'group_bbox_max', 'icolGroup', 'icol', 'radius', 'group_atoms', 'group_nAtoms', 'n_groups', 'group_cap', 'excl_icol', 'excl_count'):
        if k in topology:
            v = topology[k]
            out[k] = np.asarray(v) if k != 'defects_json' else v
    return out


def build_hessian_bundle(
    relaxed_npz: Path,
    topology_npz: Path,
    out_npz: Path,
    work_dir: Path | None = None,
    dx: float = 1e-4,
    rigid_shift: float = 1e6,
    source: str = 'topology',
    compare_mmff: bool = False,
) -> dict:
    """Hessian: exported 03_topology springs/params at relaxed coordinates (02_relaxed.npz pos).

    Topology (neigh_idx, bond_k, bond_l0, stick_class, surface AABBs) is fixed at export; only
    coordinates change during relax. Hessian evaluates the same springs at relaxed pos.

    source:
      topology — analytical K from 03_topology.npz (default; same params as export)
      mmff     — C++ MMFF FD diagnostic (different angle model; not the export FF)
    compare_mmff — if True, also run MMFF FD and record parity_max_rel_diff in status
    """
    from pyBall import FTIR
    from pyBall.io.crystal_npz import load_crystal_npz, load_topology_npz, topology_has_mmffl, validate_topology_crystal_parity

    relaxed_npz = Path(relaxed_npz).resolve()
    topology_npz = Path(topology_npz).resolve()
    crystal = load_crystal_npz(relaxed_npz)
    topology = load_topology_npz(topology_npz, full=True)
    if not topology_has_mmffl(topology):
        raise ValueError(f"build_hessian_bundle: {topology_npz} missing MMFFL keys (neigh_idx, bond_k, bond_l0)")
    init_npz = relaxed_npz.parent / '01_crystal.npz'
    Z_ref = np.asarray(load_crystal_npz(init_npz)['Z'], dtype=np.int32).reshape(-1) if init_npz.is_file() else np.asarray(crystal['Z'], dtype=np.int32).reshape(-1)
    bonds_ref = load_crystal_npz(init_npz).get('bonds_ij') if init_npz.is_file() else crystal.get('bonds_ij')
    crystal_chk = {'pos': crystal['pos'], 'Z': Z_ref, 'bonds_ij': bonds_ref}
    parity = validate_topology_crystal_parity(crystal_chk, topology)
    pos = np.asarray(crystal['pos'], dtype=np.float64)
    Z = Z_ref
    natoms = Z.shape[0]
    t0 = time.perf_counter()
    H_topo = FTIR.build_hessian_from_linear_topology(pos, topology['neigh_idx'], topology['bond_k'], topology['bond_l0'], stick_class=topology.get('stick_class'), neigh_count=topology.get('neigh_count'))
    M = FTIR.mass_matrix_from_Z(Z)
    parity_max_rel = float('nan')
    H_mmff = M_mmff = None
    if compare_mmff or source == 'mmff':
        work = Path(work_dir) if work_dir else Path(out_npz).parent / 'work_hessian'
        init_path, _ = init_path_from_crystal_npz(relaxed_npz, work)
        H_mmff, M_mmff, pos_mmff, _, _ = _mmff_hessian_fd(init_path, dx=dx)
        if pos_mmff.shape != pos.shape or not np.allclose(pos_mmff, pos, rtol=1e-5, atol=1e-4):
            raise ValueError("build_hessian_bundle: MMFF init pos != relaxed NPZ pos (check mol2/Z)")
        parity_max_rel = float(np.max(np.abs(H_topo - H_mmff)) / (np.max(np.abs(H_mmff)) + 1e-12))
    if source == 'mmff':
        if H_mmff is None:
            raise RuntimeError("build_hessian_bundle: mmff source failed to produce Hessian")
        H, M, hess_source = H_mmff, M_mmff, 'mmff_fd'
    elif source == 'topology':
        H, hess_source = H_topo, 'topology_linear'
    else:
        raise ValueError(f"build_hessian_bundle: unknown source={source!r} (use topology or mmff)")
    H_proj = FTIR.project_rigid_modes(H, M, pos, shift=rigid_shift)
    timing_ms = (time.perf_counter() - t0) * 1000.0
    out = Path(out_npz)
    out.parent.mkdir(parents=True, exist_ok=True)
    status = {
        'status': 'ok', 'natoms': int(natoms), 'ndof': int(3 * natoms), 'timing_ms': timing_ms,
        'hessian_source': hess_source, 'parity_max_rel_diff': parity_max_rel, 'n_sticks': parity.get('n_sticks'),
        'topology_npz': str(topology_npz), 'relaxed_npz': str(relaxed_npz),
    }
    save_kw = {
        'K': H, 'K_projected': H_proj, 'M': M, 'pos': pos, 'Z': Z,
        'source': np.array(hess_source), 'timing_ms': timing_ms, 'natoms': natoms,
        'parity_max_rel_diff': parity_max_rel, 'topology_path': np.array(str(topology_npz)),
    }
    if bonds_ref is not None and len(bonds_ref) > 0:
        save_kw['bonds_ij'] = np.asarray(bonds_ref, dtype=np.int32)
    save_kw.update(_topology_surface_arrays(topology))
    np.savez(out, **save_kw)
    write_status(out, status)
    return status


def build_hessian_mmff(relaxed_xyz: Path, out_npz: Path, dx: float = 1e-4, rigid_shift: float = 1e6) -> dict:
    """Legacy: MMFF FD from bare xyz (re-infers bonds). Prefer build_hessian_bundle."""
    from pyBall import FTIR

    H, M, pos, Z, natoms = _mmff_hessian_fd(relaxed_xyz, dx=dx)
    H_proj = FTIR.project_rigid_modes(H, M, pos, shift=rigid_shift)
    timing_ms = 0.0
    out = Path(out_npz)
    out.parent.mkdir(parents=True, exist_ok=True)
    status = {'status': 'ok', 'natoms': int(natoms), 'ndof': int(3 * natoms), 'timing_ms': timing_ms, 'hessian_source': 'mmff_fd_xyz', 'warning': 'bonds re-inferred from xyz'}
    np.savez(out, K=H, K_projected=H_proj, M=M, pos=pos, Z=Z, source=np.array('mmff_fd_xyz'), timing_ms=timing_ms, natoms=natoms)
    write_status(out, status)
    return status


def gaussian_blur_hist(hist: np.ndarray, sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return hist.copy()
    radius = int(max(1, round(3 * sigma_bins)))
    k = np.arange(-radius, radius + 1, dtype=np.float64)
    w = np.exp(-0.5 * (k / sigma_bins) ** 2)
    w /= w.sum()
    return np.convolve(hist.astype(np.float64), w, mode="same")


def modes_to_histogram(omegas, margin_frac=0.02, min_bins=64, sigma_bins=0.0, rigid_shift=1e6, weights=None):
    om = filter_vibrational_modes(omegas, rigid_shift=rigid_shift)
    om_cm = omega_internal_to_cm1(om)
    if weights is not None:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape[0] != np.asarray(omegas, dtype=np.float64).shape[0]:
            raise ValueError(f"modes_to_histogram: weights length {w.shape[0]} != omegas {np.asarray(omegas).shape[0]}")
        cap = 0.5 * np.sqrt(float(rigid_shift))
        mask = np.isfinite(omegas) & (omegas > 1e-4) & (omegas < cap)
        w = w[mask]
        if w.size != om.size:
            raise ValueError("modes_to_histogram: weights mask mismatch")
    else:
        w = None
    om_min, om_max = float(om_cm.min()), float(om_cm.max())
    span = max(om_max - om_min, 1e-8)
    margin = margin_frac * span
    lo = max(0.0, om_min - margin)
    hi = om_max + margin
    n_bins = max(min_bins, int(np.ceil(om_cm.size * 0.5)))
    edges = np.linspace(lo, hi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(om_cm, bins=edges, weights=w)
    hist = hist.astype(np.float64)
    if sigma_bins > 0:
        hist = gaussian_blur_hist(hist, sigma_bins)
    kind = "weighted_probe" if w is not None else "mode_count"
    grid_meta = json.dumps({"omega_min_cm1": lo, "omega_max_cm1": hi, "margin_frac": margin_frac, "n_bins": n_bins, "n_modes": int(om_cm.size), "units": "cm-1", "spectrum_kind": kind})
    return centers, hist, om, om_cm, grid_meta


def vibrational_mode_mask(omegas: np.ndarray, rigid_shift: float = 1e6) -> np.ndarray:
    om = np.asarray(omegas, dtype=np.float64)
    cap = 0.5 * np.sqrt(float(rigid_shift))
    return np.isfinite(om) & (om > 1e-4) & (om < cap)


def load_spectrum_npz(path) -> dict:
    """Load 05_spectrum.npz arrays required for plotting / mode picking."""
    d = np.load(path, allow_pickle=True)
    _require = ('omega_centers', 'hist', 'omegas_modes')
    missing = [k for k in _require if k not in d.files]
    if missing:
        raise KeyError(f"load_spectrum_npz: missing keys {missing} in {path}")
    out = {
        'path': str(path),
        'omega_centers': np.asarray(d['omega_centers'], dtype=np.float64),
        'hist': np.asarray(d['hist'], dtype=np.float64),
        'omegas_modes': np.asarray(d['omegas_modes'], dtype=np.float64),
    }
    for k in ('omegas_modes_vib', 'probe_weight_x', 'probe_weight_y', 'probe_weight_z', 'sigma_bins', 'timing_ms'):
        if k in d.files:
            out[k] = np.asarray(d[k])
    if 'units' in d.files:
        out['units'] = str(np.asarray(d['units']).item()) if np.asarray(d['units']).shape == () else str(d['units'])
    if 'grid_meta' in d.files:
        out['grid_meta'] = str(np.asarray(d['grid_meta']).item()) if np.asarray(d['grid_meta']).shape == () else str(d['grid_meta'])
    return out


def solve_normal_modes_from_hessian_npz(hessian_npz, rigid_shift: float = 1e6) -> dict:
    """Mass-weighted eigh of K_projected — eigenvectors for mode visualization (not stored in 05 v1.2)."""
    from pyBall import FTIR

    data = np.load(hessian_npz, allow_pickle=True)
    K = np.asarray(data["K_projected"] if "K_projected" in data else data["K"], dtype=np.float64)
    M = np.asarray(data["M"], dtype=np.float64)
    pos = np.asarray(data["pos"], dtype=np.float64)
    natoms = int(data["natoms"]) if "natoms" in data else pos.shape[0]
    if "K_projected" not in data:
        K = FTIR.project_rigid_modes(K, M, pos, shift=rigid_shift)
    m = np.diag(M).reshape(natoms, 3)[:, 0]
    m_sqrt = np.repeat(np.sqrt(m), 3)
    m_inv_sqrt = 1.0 / m_sqrt
    K_mw = (m_inv_sqrt[:, None] * K) * m_inv_sqrt[None, :]
    w, v = np.linalg.eigh(K_mw)
    omegas_modes = np.sqrt(np.maximum(w, 0.0))
    omegas_cm = omega_internal_to_cm1(omegas_modes)
    modes_cart = v * m_inv_sqrt[:, None]
    mask = vibrational_mode_mask(omegas_modes, rigid_shift=rigid_shift)
    Z = np.asarray(data["Z"], dtype=np.int32).reshape(-1) if "Z" in data else None
    return {
        'path': str(hessian_npz),
        'pos': pos,
        'Z': Z,
        'natoms': natoms,
        'omegas_modes': omegas_modes,
        'omegas_cm': omegas_cm,
        'modes_cart': modes_cart,
        'vib_mask': mask,
        'vib_indices': np.flatnonzero(mask),
        'omegas_modes_vib': omegas_modes[mask],
        'omegas_cm_vib': omegas_cm[mask],
        'rigid_shift': float(rigid_shift),
    }


def mode_probe_weights(v: np.ndarray, m_inv_sqrt: np.ndarray, natoms: int, pol_dir, charges=None) -> np.ndarray:
    """Squared dipole coupling |e·u_mode|² per mode (FTIR vibration_spectrum_from_modes convention)."""
    pol = np.asarray(pol_dir, dtype=np.float64).reshape(3)
    pn = float(np.linalg.norm(pol))
    if pn <= 0:
        raise ValueError("mode_probe_weights: zero polarization")
    pol /= pn
    if charges is None:
        charges = np.ones(natoms, dtype=np.float64)
    else:
        charges = np.asarray(charges, dtype=np.float64)
    f = np.zeros(v.shape[0], dtype=np.float64)
    for n in range(natoms):
        f[n * 3 : n * 3 + 3] = charges[n] * pol
    f_tilde = v.T @ (f * m_inv_sqrt)
    return (f_tilde ** 2).astype(np.float64)


POLARIZATION_PRESETS = {
    "hist": None,
    "x": np.array([1.0, 0.0, 0.0]),
    "y": np.array([0.0, 1.0, 0.0]),
    "z": np.array([0.0, 0.0, 1.0]),
}


def plot_spectrum_line(path, omega_centers, hist, title="spectrum", xlabel="Wavenumber (cm$^{-1}$)", ylabel="Mode count"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.asarray(omega_centers, dtype=np.float64)
    y = np.asarray(hist, dtype=np.float64)
    if x.size == 0:
        raise ValueError("plot_spectrum_line: empty spectrum")
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=150)
    ax.plot(x, y, color="#1d4ed8", lw=1.6)
    ax.fill_between(x, 0, y, color="#3b82f6", alpha=0.18)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(float(x[0]), float(x[-1]))
    ax.set_ylim(0, float(np.max(y)) * 1.08 if np.max(y) > 0 else 1.0)
    ax.grid(True, alpha=0.25, ls=":")
    fig.tight_layout()
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def solve_spectrum(
    hessian_npz: Path,
    out_npz: Path,
    margin_frac: float = 0.02,
    min_bins: int = 64,
    sigma_bins: float = 1.5,
    out_plot: Path | None = None,
    rigid_shift: float = 1e6,
) -> dict:
    t0 = time.perf_counter()
    nm = solve_normal_modes_from_hessian_npz(hessian_npz, rigid_shift=rigid_shift)
    omegas_modes = nm['omegas_modes']
    natoms = nm['natoms']
    pos = nm['pos']
    om_vib = filter_vibrational_modes(omegas_modes, rigid_shift=rigid_shift)
    m = np.diag(np.load(hessian_npz, allow_pickle=True)["M"]).reshape(natoms, 3)[:, 0]
    m_inv_sqrt = 1.0 / np.repeat(np.sqrt(m), 3)
    v = nm['modes_cart'] / m_inv_sqrt[:, None]
    probe_w = {}
    for pname, pdir in POLARIZATION_PRESETS.items():
        if pname == "hist":
            continue
        probe_w[pname] = mode_probe_weights(v, m_inv_sqrt, natoms, pdir)
    centers, hist, _om_raw, om_cm, grid_meta = modes_to_histogram(omegas_modes, margin_frac=margin_frac, min_bins=min_bins, sigma_bins=sigma_bins, rigid_shift=rigid_shift)
    timing_ms = (time.perf_counter() - t0) * 1000.0
    out = Path(out_npz)
    out.parent.mkdir(parents=True, exist_ok=True)
    status = {"status": "ok", "n_modes": int(om_vib.size), "omega_min_cm1": float(om_cm.min()), "omega_max_cm1": float(om_cm.max()), "timing_ms": timing_ms, "spectrum_kind_hist": "mode_count"}
    np.savez(
        out,
        omega_centers=centers,
        hist=hist,
        omegas_modes=omegas_modes,
        omegas_modes_vib=om_vib,
        probe_weight_x=probe_w["x"],
        probe_weight_y=probe_w["y"],
        probe_weight_z=probe_w["z"],
        grid_meta=np.array(grid_meta),
        timing_ms=timing_ms,
        sigma_bins=sigma_bins,
        units=np.array("cm-1"),
    )
    write_status(out, status)
    if out_plot:
        plot_spectrum_line(out_plot, centers, hist, title=Path(out_plot).stem)
    return status


def accumulate_spectra(
    input_paths: list[Path],
    out_dir: Path,
    rigid_shift: float = 1e6,
    margin_frac: float = 0.02,
    min_bins: int = 128,
    sigma_bins: float = 1.5,
) -> dict:
    """Naive mode-count histogram + probe-weighted spectra (pol x/y/z)."""
    om_list, wx_list, wy_list, wz_list = [], [], [], []
    for p in input_paths:
        d = np.load(p)
        if "omegas_modes" not in d:
            raise ValueError(f"accumulate_spectra: missing omegas_modes in {p}")
        om = np.asarray(d["omegas_modes"], dtype=np.float64)
        mask = vibrational_mode_mask(om, rigid_shift=rigid_shift)
        om_list.append(om[mask])
        for key, lst in [("probe_weight_x", wx_list), ("probe_weight_y", wy_list), ("probe_weight_z", wz_list)]:
            if key in d:
                lst.append(np.asarray(d[key], dtype=np.float64)[mask])
            else:
                lst.append(np.ones(int(mask.sum()), dtype=np.float64))
    if not om_list:
        raise ValueError("no modes in inputs")
    om_all = np.concatenate(om_list)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = out_dir.parent / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    results = {"status": "ok", "n_crystals": len(input_paths), "n_modes_total": int(om_all.size), "spectra": {}}
    variants = [
        ("histogram", None, "Mode count (naive histogram)", "spectrum_histogram.png"),
        ("pol_x", "x", "Probe polarization ‖x̂ (|e·u|² weight)", "spectrum_pol_x.png"),
        ("pol_y", "y", "Probe polarization ‖ŷ", "spectrum_pol_y.png"),
        ("pol_z", "z", "Probe polarization ‖ẑ", "spectrum_pol_z.png"),
    ]
    weight_map = {"x": wx_list, "y": wy_list, "z": wz_list}
    for name, wkey, title, plot_name in variants:
        if wkey is None:
            centers, hist, _, _om_cm, grid_meta = modes_to_histogram(om_all, margin_frac=margin_frac, min_bins=min_bins, sigma_bins=sigma_bins, rigid_shift=rigid_shift)
            w_all = None
        else:
            w_all = np.concatenate(weight_map[wkey])
            centers, hist, _, _om_cm, grid_meta = modes_to_histogram(om_all, margin_frac=margin_frac, min_bins=min_bins, sigma_bins=sigma_bins, rigid_shift=rigid_shift, weights=w_all)
        npz_path = out_dir / f"spectrum_{name}.npz"
        np.savez(npz_path, omega_centers=centers, hist=hist, n_crystals=len(input_paths), n_modes_total=int(om_all.size), units=np.array("cm-1"), grid_meta=np.array(grid_meta), spectrum_kind=np.array(name))
        plot_path = plot_dir / plot_name
        ylabel = "Weighted mode sum" if name != "histogram" else "Mode count"
        plot_spectrum_line(plot_path, centers, hist, title=f"{title} (n={len(input_paths)} crystals)", ylabel=ylabel)
        results["spectra"][name] = {"npz": str(npz_path), "plot": str(plot_path), "omega_min_cm1": float(centers[0]), "omega_max_cm1": float(centers[-1])}
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description="Nanocrystal pipeline Python stages")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_relax = sub.add_parser("relax", help="MMFF relax init → 02_relaxed.npz")
    p_relax.add_argument("--init-xyz", default=None)
    p_relax.add_argument("--init-mol2", default=None)
    p_relax.add_argument("--init-npz", default=None, help="01_crystal.npz (writes temp mol2 with bonds_ij when present)")
    p_relax.add_argument("--work-dir", default=None, help="scratch for NPZ→mol2/xyz conversion")
    p_relax.add_argument("--out-npz", required=True)
    p_relax.add_argument("--out-xyz", default=None)
    p_relax.add_argument("--nstep-max", type=int, default=2000)
    p_relax.add_argument("--fconv", type=float, default=1e-4)
    p_relax.add_argument("--allow-unconverged", action="store_true")

    p_hess = sub.add_parser("hessian", help="Hessian → 04_hessian.npz (requires relaxed + topology NPZ bundle)")
    p_hess.add_argument("--relaxed-npz", default=None, help="02_relaxed.npz (pos, Z, bonds_ij)")
    p_hess.add_argument("--topology-npz", default=None, help="03_topology.npz (MMFFL springs + surface AABBs)")
    p_hess.add_argument("--relaxed-xyz", default=None, help="legacy: bare xyz (re-infers bonds — avoid for defects)")
    p_hess.add_argument("--out-npz", required=True)
    p_hess.add_argument("--work-dir", default=None)
    p_hess.add_argument("--dx", type=float, default=1e-4)
    p_hess.add_argument("--rigid-shift", type=float, default=1e6)
    p_hess.add_argument("--source", choices=('topology', 'mmff'), default='topology', help="topology=exported springs (default); mmff=C++ FD diagnostic only")
    p_hess.add_argument("--compare-mmff", action='store_true', help="also run MMFF FD and record parity_max_rel_diff vs topology K")

    p_spec = sub.add_parser("spectrum", help="eigh + histogram → 05_spectrum.npz")
    p_spec.add_argument("--hessian-npz", required=True)
    p_spec.add_argument("--out-npz", required=True)
    p_spec.add_argument("--margin-frac", type=float, default=0.02)
    p_spec.add_argument("--min-bins", type=int, default=64)
    p_spec.add_argument("--sigma-bins", type=float, default=1.5)
    p_spec.add_argument("--out-plot", default=None, help="matplotlib PNG line spectrum")
    p_spec.add_argument("--rigid-shift", type=float, default=1e6)

    p_acc = sub.add_parser("accumulate", help="merge per-crystal spectra (histogram + pol x/y/z)")
    p_acc.add_argument("--inputs", nargs="+", required=True)
    p_acc.add_argument("--out-dir", required=True, help="directory for spectrum_*.npz (plots go to ../plots/)")

    args = ap.parse_args(argv)
    if args.cmd == "relax":
        bonds_ij = None
        if args.init_npz:
            work = Path(args.work_dir) if args.work_dir else Path(args.out_npz).parent / "work_relax"
            init, crystal = init_path_from_crystal_npz(args.init_npz, work)
            bonds_ij = crystal.get("bonds_ij")
            Z_init = crystal.get("Z")
        else:
            init = args.init_xyz or args.init_mol2
            Z_init = None
            if not init:
                raise ValueError("need --init-npz, --init-xyz, or --init-mol2")
        relax(init, args.out_npz, out_xyz=args.out_xyz, nstep_max=args.nstep_max, fconv=args.fconv, allow_unconverged=args.allow_unconverged, bonds_ij=bonds_ij, Z_init=Z_init)
    elif args.cmd == "hessian":
        if args.relaxed_npz and args.topology_npz:
            build_hessian_bundle(args.relaxed_npz, args.topology_npz, args.out_npz, work_dir=args.work_dir, dx=args.dx, rigid_shift=args.rigid_shift, source=args.source, compare_mmff=args.compare_mmff)
        elif args.relaxed_xyz:
            build_hessian_mmff(args.relaxed_xyz, args.out_npz, dx=args.dx, rigid_shift=args.rigid_shift)
        else:
            raise ValueError("hessian requires --relaxed-npz + --topology-npz (preferred) or legacy --relaxed-xyz")
    elif args.cmd == "spectrum":
        solve_spectrum(args.hessian_npz, args.out_npz, margin_frac=args.margin_frac, min_bins=args.min_bins, sigma_bins=args.sigma_bins, out_plot=args.out_plot, rigid_shift=args.rigid_shift)
    elif args.cmd == "accumulate":
        status = accumulate_spectra(args.inputs, args.out_dir)
        print(json.dumps(status))


if __name__ == "__main__":
    main()
