import os
import numpy as np
import matplotlib.pyplot as plt

from .STM import DEFAULT_EXPORT_DIR, project_pdos_to_grid
import pyBall.FireCore as fc
import pyBall.atomicUtils as pu


def set_export_dir(export_dir=None):
    """Create export directory if needed, return absolute path."""
    if export_dir is None:
        export_dir = DEFAULT_EXPORT_DIR
    export_dir = os.path.abspath(export_dir)
    os.makedirs(export_dir, exist_ok=True)
    return export_dir


def save_plot(fig, name, export_dir=None, dpi=130):
    """Save figure to export_dir (default: DEFAULT_EXPORT_DIR)."""
    export_dir = set_export_dir(export_dir)
    path = os.path.join(export_dir, f"{name}.png")
    fig.savefig(path, dpi=dpi)
    print(f"  → {path}")
    return path


def plot_2d_map(ax, data_xy, x_arr, y_arr, title, cmap='viridis', vmin=None, vmax=None):
    """Plot 2D data as heatmap with proper axes.

    Theory: I(x,y) is a scalar STM signal evaluated at constant energy and height, e.g.
      - Caroli transmission: T(E,x,y) = Tr[Γ_L G^r Γ_R G^a]
      - Response probe: resp(E,x,y) = x† Γ_measure x, where x solves [(E+iη)S - H - Σ]x=b
      - MO overlap intensity: |t_n(x,y)|^2 with t_n = v(x,y)·C_n
    """
    im = ax.imshow(np.asarray(data_xy).T, origin='lower',
                   extent=[x_arr[0], x_arr[-1], y_arr[0], y_arr[-1]],
                   cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_title(title)
    return im


def plot_1d_scan(ax, s_arr, y_arr, title, ylabel):
    """Plot 1D scan data."""
    ax.plot(s_arr, y_arr, 'o-', lw=1.2)
    ax.set_xlabel('s (Å)')
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def plot_orbital_on_plane(atomPos, atomTypes, mo_idx, z=1.0, size=20.0, nx=200, ny=200,
                        export_dir=None, filename=None, overlay_atoms=True, label=None):
    """Plot molecular orbital on a plane above the molecular plane.

    Args:
        atomPos: atomic positions (natoms, 3)
        atomTypes: atomic types (natoms,)
        mo_idx: MO index (1-based)
        z: height above molecular plane (Å)
        size: grid size (Å)
        nx, ny: grid resolution
        export_dir: export directory path
        filename: output filename (without .png)
        overlay_atoms: whether to overlay atoms
        label: orbital label (e.g., "HOMO", "LUMO+1") for title

    Returns:
        mo_2d: 2D orbital values (ny, nx)
    """
    export_dir = set_export_dir(export_dir)

    # Center grid on atom center of mass for orbital evaluation
    atom_center = atomPos.mean(axis=0)
    origin = atom_center
    normal = np.array([0.0, 0.0, 1.0])

    # Create orthonormal in-plane basis
    n = normal / np.linalg.norm(normal)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial, n)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    x_axis = trial - np.dot(trial, n) * n
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(n, x_axis)
    y_axis /= np.linalg.norm(y_axis)

    # Build grid
    xs = np.linspace(-size/2, size/2, nx)
    ys = np.linspace(-size/2, size/2, ny)
    X_local, Y_local = np.meshgrid(xs, ys)

    # Generate 3D points at height z above plane
    plane_origin = origin + z * n
    points = plane_origin[None, :] + X_local[..., None] * x_axis[None, None, :] + Y_local[..., None] * y_axis[None, None, :]
    points = points.reshape(-1, 3)

    # Evaluate orbital
    # Fortran expects points as vector from evaluation point to origin (negative of absolute coords)
    mo_values = fc.orb2points(-points, iMO=mo_idx, ikpoint=1)
    mo_2d = mo_values.reshape(ny, nx)

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    extent = [X_local.min(), X_local.max(), Y_local.min(), Y_local.max()]
    vmax = np.max(np.abs(mo_2d))
    vmin = -vmax
    im = ax.imshow(mo_2d, origin="lower", extent=extent, cmap="bwr", aspect="equal", vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=ax)

    # Overlay atoms
    if overlay_atoms:
        # Project atoms relative to grid origin (atom center of mass)
        apos_plane = []
        for pos in atomPos:
            vec = pos - origin
            x = np.dot(vec, x_axis)
            y = np.dot(vec, y_axis)
            apos_plane.append([x, y])
        apos_plane = np.array(apos_plane)
        ax.scatter(apos_plane[:, 0], apos_plane[:, 1], c='green', s=50, edgecolors='black', linewidths=0.5, zorder=10)

    ax.set_xlabel("x in plane [Å]")
    ax.set_ylabel("y in plane [Å]")
    title_label = f" ({label})" if label else ""
    ax.set_title(f"MO {mo_idx}{title_label} at z={z} Å")
    plt.tight_layout()

    if filename is None:
        filename = f"mo{mo_idx:04d}_z{z:.2f}"
    path = os.path.join(export_dir, f"{filename}.png")
    plt.savefig(path, dpi=200)
    plt.close()
    print(f"Saved to {path}")

    return mo_2d


def plot_atoms(ax, atomPos, atomTypes, color='green', ms=3, label=None):
    """Plot atoms as distinguishable dots on 2D map.

    Theory: Atomic positions provide structural context for STM images.
    Green dots are visible with 'bwr' colormap (blue-white-red).

    Args:
        ax: matplotlib axes
        atomPos: (natoms, 3) array of positions in Å
        atomTypes: (natoms,) array of atomic numbers
        color: dot color (default='green' for visibility with bwr)
        ms: marker size
        label: legend label (optional)
    """
    ax.plot(atomPos[:, 0], atomPos[:, 1], '.', color=color, ms=ms, label=label)
    if label:
        ax.legend(fontsize=7)


def plot_dos(ax, E_grid, A, eigen, title, E_F):
    """Plot PDOS with eigenvalues as vertical lines.

    Theory: PDOS A(E) = (1/π) Im[G^r(E)] shows spectral weight distribution.
    Eigenvalues ε_n appear as peaks (broadened by coupling to leads).

    Args:
        ax: matplotlib axes
        E_grid: (nE,) energy grid in eV
        A: (nE, norb) spectral function
        eigen: (norb,) orbital eigenvalues
        title: plot title
        E_F: Fermi level
    """
    A_total = A.sum(axis=1)
    ax.plot(E_grid, A_total, 'k-', lw=1.2)
    ax.axvline(E_F, color='r', ls='--', lw=1.0, label='E_F')
    for e in eigen:
        ax.axvline(e, color='gray', lw=0.5, alpha=0.5)
    ax.set_xlabel('E (eV)'); ax.set_ylabel('PDOS (arb)')
    ax.set_title(title); ax.legend(fontsize=8)


def project_orbital_to_grid(C, mo_idx, atomTypes, atomPos, orb2atom, norb_per,
                           fdata_basis_dir, step=0.2, margin=3.5, grid_spec=None):
    """Project a single MO (signed and squared) onto 3D grid.

    Theory: Molecular orbital ψ_n(r) = Σ_μ C_μn φ_μ(r), where φ_μ are atomic orbitals.
    The squared density |ψ_n(r)|² is the real-space PDOS for that orbital.
    The signed wavefunction ψ_n(r) reveals nodal planes (sign flips).

    IMPORTANT: Fireball uses Fortran/Ortega order [s, py, pz, px] for p-orbitals,
    but OpenCL GridProjector expects [s, px, py, pz]. We must remap coefficients.

    Args:
        C: (norb, nmo) MO eigenvectors (in Fortran order)
        mo_idx: index of MO to project
        atomTypes: (natoms,) atomic numbers
        atomPos: (natoms, 3) positions in Å
        orb2atom: (norb,) orbital-to-atom mapping
        norb_per: (natoms,) orbitals per atom
        fdata_basis_dir: path to Fdata/basis for GridProjector
        step: grid spacing in Å
        margin: margin around atoms in Å
        grid_spec: optional dict with grid parameters (origin, dA, dB, dC, ngrid)

    Returns:
        psi_signed: (nx, ny, nz) signed wavefunction
        psi_sq: (nx, ny, nz) squared density
        grid_spec: dict with grid parameters
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    from pyBall.FireballOCL.FdataParser import FdataParser

    fparser = FdataParser(os.path.dirname(fdata_basis_dir))
    fparser.parse_info()
    z_list = sorted(set(int(z) for z in atomTypes))
    norb_for = {nz: sum(2*l+1 for l in fparser.species_info[nz]['lssh']) for nz in z_list}
    norb_per_check = [norb_for[int(z)] for z in atomTypes]
    starts = np.zeros(len(atomTypes)+1, dtype=int)
    for ia, no in enumerate(norb_per_check):
        starts[ia+1] = starts[ia] + no
    norb_total = starts[-1]

    # Extract MO coefficients for this orbital
    mo_coeffs = C[:, mo_idx] if C.shape[1] > mo_idx else C[mo_idx, :]
    if len(mo_coeffs) != norb_total:
        raise ValueError(f"MO coefficients size mismatch: {len(mo_coeffs)} vs {norb_total}")

    # Remap coefficients from Fortran order [s, py, pz, px] to OpenCL order [s, px, py, pz]
    # Permutation: [0, 3, 1, 2] maps: 0->0 (s), 1->3 (py->pz), 2->1 (pz->px), 3->2 (px->py)
    _ORT_SPP_TO_STD = np.array([0, 3, 1, 2], dtype=int)
    mo_coeffs_remapped = mo_coeffs.copy()
    for ia in range(len(atomTypes)):
        no = norb_per_check[ia]
        i0 = starts[ia]
        if no == 4:  # sp atom: apply remapping
            # Fortran: [s, py, pz, px] -> OpenCL: [s, px, py, pz]
            coeffs_fortran = mo_coeffs[i0:i0+4]
            coeffs_opencl = coeffs_fortran[_ORT_SPP_TO_STD]
            mo_coeffs_remapped[i0:i0+4] = coeffs_opencl

    # Build density matrix for single MO (signed) - still diagonal for orbital projection
    # Note: This is a workaround - we're using density projection for orbital projection
    # The diagonal elements are the orbital coefficients
    numorb_max = max(norb_per_check)
    rho_signed = np.zeros((len(atomTypes), 1, numorb_max, numorb_max))
    rho_sq = np.zeros((len(atomTypes), 1, numorb_max, numorb_max))
    neigh_j = np.ones((len(atomTypes), 1), dtype=np.int32)
    for ia in range(len(atomTypes)):
        neigh_j[ia, 0] = ia + 1
        no = norb_per_check[ia]
        i0 = starts[ia]
        for io in range(no):
            idx = i0 + io
            if idx < len(mo_coeffs_remapped):
                c = mo_coeffs_remapped[idx]
                rho_signed[ia, 0, io, io] = float(c)
                rho_sq[ia, 0, io, io] = float(c * c)
    neighn = np.ones(len(atomTypes), dtype=np.int32)

    class Neighs:
        pass
    neighs = Neighs()
    neighs.neigh_j = neigh_j
    neighs.neighn = neighn
    neighs.neigh_max = 1
    neighs.iatyp = atomTypes.copy()
    num_orb_arr = np.zeros(20, dtype=np.int32)
    for nz in z_list:
        num_orb_arr[nz] = norb_for[nz]
    neighs.num_orb = num_orb_arr
    neighs.nzx = np.array(z_list, dtype=np.int32)

    if grid_spec is None:
        pmin = atomPos.min(0) - margin
        pmax = atomPos.max(0) + margin
        ngrid = ((np.ceil((pmax - pmin) / step).astype(int) + 7) // 8) * 8
        grid_spec = {'origin': pmin, 'dA': [step, 0, 0], 'dB': [0, step, 0], 'dC': [0, 0, step], 'ngrid': ngrid}

    projector = ocl_grid.GridProjector(fdata_basis_dir)
    projector.load_basis(z_list)
    atoms_dict = {'pos': atomPos, 'Rcut': np.array([3.0 for z in atomTypes]), 'type': atomTypes}

    psi_signed = projector.project(rho_signed, neighs, atoms_dict, grid_spec, nMaxAtom=64)
    psi_sq = projector.project(rho_sq, neighs, atoms_dict, grid_spec, nMaxAtom=64)

    return psi_signed, psi_sq, grid_spec


def project_orbital_to_grid_v2(C, mo_idx, atomTypes, atomPos, orb2atom, norb_per,
                              fdata_basis_dir, grid_spec):
    """
    Project a single MO onto 3D grid using the orbital projection kernel (not density kernel).

    Computes ψ(r) = Σ_i C_i φ_i(r) (signed wavefunction, not squared density).

    IMPORTANT: Fireball uses Fortran/Ortega order [s, py, pz, px] for p-orbitals,
    but OpenCL orbital projection kernel expects [s, px, py, pz]. We must remap coefficients.

    Args:
        C: (norb, nmo) MO eigenvectors (in Fortran order)
        mo_idx: index of MO to project
        atomTypes: (natoms,) atomic numbers
        atomPos: (natoms, 3) positions in Å
        orb2atom: (norb,) orbital-to-atom mapping
        norb_per: (natoms,) orbitals per atom
        fdata_basis_dir: path to Fdata/basis for GridProjector
        grid_spec: dict with grid parameters (origin, dA, dB, dC, ngrid)

    Returns:
        psi: (nx, ny, nz) signed wavefunction
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    from pyBall.FireballOCL.FdataParser import FdataParser

    fparser = FdataParser(os.path.dirname(fdata_basis_dir))
    fparser.parse_info()
    z_list = sorted(set(int(z) for z in atomTypes))
    ABOHR = 0.529177
    rcut_map = {}
    for z in z_list:
        wfs = fparser.find_wf(int(z))
        if len(wfs) == 0:
            raise RuntimeError(f"project_orbital_to_grid_v2(): no .wf files found for species Z={z}")
        rcut_map[int(z)] = max(ABOHR * fparser.read_wf(f)['rcutoff'] for f in wfs)

    mo_coeffs = C[:, mo_idx].copy() if C.shape[0] == len(orb2atom) else C[mo_idx, :].copy()

    # Build per-atom coefficient arrays
    natoms = len(atomTypes)
    norb_per = np.asarray(norb_per, dtype=np.int32)
    starts = np.zeros(natoms + 1, dtype=np.int32)
    starts[1:] = np.cumsum(norb_per)
    norb_total = int(starts[-1])
    if len(mo_coeffs) != norb_total:
        raise ValueError(f"project_orbital_to_grid_v2(): MO coefficient size mismatch {len(mo_coeffs)} != {norb_total}")
    numorb_max = int(norb_per.max())
    coeffs_atoms = np.zeros((natoms, numorb_max), dtype=np.float32)

    # Remap from Fortran [s, py, pz, px] to OpenCL [s, px, py, pz] to match density kernel
    _ORT_SPP_TO_STD = np.array([0, 3, 1, 2], dtype=int)
    for ia in range(natoms):
        no = int(norb_per[ia])
        i0 = int(starts[ia])
        if no == 4:
            coeffs_atoms[ia, :4] = mo_coeffs[i0:i0 + 4][_ORT_SPP_TO_STD]
        else:
            coeffs_atoms[ia, :no] = mo_coeffs[i0:i0 + no]

    # Prepare atoms dict for GridProjector
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([rcut_map[int(z)] for z in atomTypes], dtype=np.float32),
        'type': atomTypes
    }

    # Use new orbital projection kernel
    projector = ocl_grid.GridProjector(fdata_basis_dir)
    projector.load_basis(z_list)
    psi = projector.project_orbital(coeffs_atoms, norb_per, atoms_dict, grid_spec, nMaxAtom=64)

    return psi

