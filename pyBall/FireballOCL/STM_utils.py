import os
import numpy as np
import matplotlib.pyplot as plt

from .STM import DEFAULT_EXPORT_DIR, project_pdos_to_grid
import pyBall.FireCore as fc
import pyBall.atomicUtils as pu


# =============================================================================
# Orbital Mapping Conventions (Critical for Fortran/OpenCL parity)
# =============================================================================

# Fortran Fireball order: [s, py, pz, px] (spherical harmonics m=-1,0,+1)
# OpenCL Grid kernel order: [px, py, pz, s] (Cartesian x,y,z)
# OpenCL Hamiltonian order: [s, px, py, pz] (intermediate convention)

_PERM_FORT_TO_GRID = np.array([3, 1, 2, 0], dtype=np.int32)   # [s,py,pz,px] -> [px,py,pz,s]
_PERM_FORT_TO_HAM  = np.array([0, 3, 1, 2], dtype=np.int32)   # [s,py,pz,px] -> [s,px,py,pz]


def remap_coeffs_fortran_to_grid(coeffs_fortran, norb_per):
    """
    Remap MO coefficients from Fortran order to OpenCL Grid kernel order.

    CRITICAL: This is the correct mapping for Grid.cl kernels (project_orbital,
    project_orbital_points, etc.) which expect [px, py, pz, s] order.

    Fortran Fireball uses spherical harmonic order [s, py, pz, px] where:
      - py corresponds to m=-1 (Y_{1,-1})
      - pz corresponds to m=0 (Y_{1,0})
      - px corresponds to m=+1 (Y_{1,+1})

    OpenCL Grid kernels use Cartesian order [px, py, pz, s] for vectorized compute.

    Args:
        coeffs_fortran: (natoms, max_norb) array in Fortran order, or (norb_total,) flat array
        norb_per: (natoms,) orbitals per atom (1 for H, 4 for C/N/O, etc.)

    Returns:
        coeffs_opencl: (natoms, 4) array in Grid kernel order [px,py,pz,s]
                      Padded orbitals are set to 0.

    Example:
        >>> coeffs_f = np.array([[0.5, 0.0, 0.0, 0.0],   # H: only s
        ...                      [0.1, 0.2, 0.3, 0.4]]) # O: s,py,pz,px
        >>> remap_coeffs_fortran_to_grid(coeffs_f, [1, 4])
        array([[0.   , 0.   , 0.   , 0.5  ],   # H: [0,0,0,s]
               [0.4  , 0.2  , 0.3  , 0.1  ]])  # O: [px,py,pz,s]
    """
    coeffs_fortran = np.asarray(coeffs_fortran)
    norb_per = np.asarray(norb_per, dtype=np.int32)
    natoms = len(norb_per)

    # Handle both (natoms, max_norb) and flat (norb_total,) input
    if coeffs_fortran.ndim == 1:
        # Flat array - need to reshape per atom
        starts = np.zeros(natoms + 1, dtype=np.int32)
        starts[1:] = np.cumsum(norb_per)
        coeffs_2d = np.zeros((natoms, 4), dtype=coeffs_fortran.dtype)
        for ia in range(natoms):
            i0 = starts[ia]
            no = norb_per[ia]
            coeffs_2d[ia, :no] = coeffs_fortran[i0:i0+no]
        coeffs_fortran = coeffs_2d

    coeffs_opencl = np.zeros((natoms, 4), dtype=np.float32)

    for ia in range(natoms):
        no = int(norb_per[ia])
        if no == 1:
            # H or similar: only s orbital -> pack as [0,0,0,s]
            coeffs_opencl[ia, 3] = float(coeffs_fortran[ia, 0])
        elif no == 4:
            # sp atom: remap [s,py,pz,px] -> [px,py,pz,s]
            coeffs_opencl[ia, :] = coeffs_fortran[ia, :4][_PERM_FORT_TO_GRID]
        else:
            # Other cases: just copy what we have, pad with zeros
            coeffs_opencl[ia, :no] = coeffs_fortran[ia, :no]

    return coeffs_opencl


def remap_coeffs_fortran_to_hamiltonian(coeffs_fortran, norb_per):
    """
    Remap MO coefficients from Fortran order to OpenCL Hamiltonian order.

    For OCL_Hamiltonian kernels which use [s, px, py, pz] order.

    Args:
        coeffs_fortran: (natoms, max_norb) array in Fortran order
        norb_per: (natoms,) orbitals per atom

    Returns:
        coeffs_opencl: (natoms, 4) array in Hamiltonian order [s,px,py,pz]
    """
    coeffs_fortran = np.asarray(coeffs_fortran)
    norb_per = np.asarray(norb_per, dtype=np.int32)
    natoms = len(norb_per)

    coeffs_opencl = np.zeros((natoms, 4), dtype=np.float32)

    for ia in range(natoms):
        no = int(norb_per[ia])
        if no == 4:
            coeffs_opencl[ia, :4] = coeffs_fortran[ia, :4][_PERM_FORT_TO_HAM]
        else:
            coeffs_opencl[ia, :no] = coeffs_fortran[ia, :no]

    return coeffs_opencl


def build_grid_spec(atom_center, Lx, Ly, Lz, nx, ny, nz, voxel_center_sampling=True):
    """
    Build grid specification for OpenCL projection.

    IMPORTANT: Fortran orb2points() samples at voxel centers (0.5*dx offsets from origin).
    OpenCL grid kernels sample at voxel corners (origin + i*dA + j*dB + k*dC).
    For parity, we shift origin by half-step when voxel_center_sampling=True.

    Args:
        atom_center: (3,) center of molecule
        Lx, Ly, Lz: grid extents in Å
        nx, ny, nz: grid resolution
        voxel_center_sampling: If True, shift origin so samples are at voxel centers
                              (matches Fortran orb2points behavior)

    Returns:
        grid_spec: dict with 'origin', 'dA', 'dB', 'dC', 'ngrid'
        extent: [xmin, xmax, ymin, ymax] for matplotlib plotting
    """
    grid_origin = atom_center - np.array([Lx/2.0, Ly/2.0, Lz/2.0])
    grid_spec = {
        'origin': grid_origin.copy(),
        'dA': np.array([Lx/nx, 0.0, 0.0]),
        'dB': np.array([0.0, Ly/ny, 0.0]),
        'dC': np.array([0.0, 0.0, Lz/nz]),
        'ngrid': np.array([nx, ny, nz], dtype=np.int32)
    }

    if voxel_center_sampling:
        # Shift to voxel centers to match Fortran sampling
        grid_spec['origin'] = grid_spec['origin'] + 0.5 * (grid_spec['dA'] + grid_spec['dB'] + grid_spec['dC'])

    extent = [grid_origin[0], grid_origin[0] + Lx, grid_origin[1], grid_origin[1] + Ly]
    return grid_spec, extent


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

    CRITICAL: Uses correct coefficient remapping from Fortran [s,py,pz,px]
    to OpenCL Grid kernel order [px,py,pz,s] for parity with orb2points.

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

    # Extract MO coefficients - use row indexing like test_h2o_orbital_comparison.py
    # C is typically (nmo, norb) from fc.get_wfcoef()
    norb_total = sum(norb_per)
    mo_coeffs = C[mo_idx, :].copy() if C.shape[0] > C.shape[1] else C[mo_idx, :].copy()
    
    # Verify size
    if len(mo_coeffs) != norb_total:
        # Try column indexing as fallback
        mo_coeffs = C[:, mo_idx].copy()
        if len(mo_coeffs) != norb_total:
            raise ValueError(f"project_orbital_to_grid_v2(): MO coefficient size mismatch {len(mo_coeffs)} != {norb_total}")

    # Use the correct remapping function [px,py,pz,s] for Grid kernels
    coeffs_atoms = remap_coeffs_fortran_to_grid(mo_coeffs, norb_per)

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


def project_orbital_to_points_exp(C, mo_idx, atomTypes, atomPos, orb2atom, norb_per,
                                 fdata_basis_dir, points, beta=1.0, r0=3.0, rcut=20.0):
    """Evaluate a single MO at arbitrary points using OpenCL exponential radial kernel.

    This mirrors project_orbital_to_points(), but replaces Fireball basis radial functions
    with exp(-beta*(r-r0)) and uses a larger cutoff (rcut) for vacuum tunneling.

    Args:
        C: (norb, nmo) or (nmo, norb) MO eigenvectors
        mo_idx: index of MO to project (0-based)
        atomTypes: (natoms,) atomic numbers
        atomPos: (natoms, 3) positions in Å
        orb2atom: unused here (kept for API symmetry)
        norb_per: (natoms,) orbitals per atom
        fdata_basis_dir: path to Fdata/basis for GridProjector (used for species mapping)
        points: (n_points, 3) points to evaluate at
        beta, r0: exponential decay parameters
        rcut: cutoff radius in Å for the exponential kernel
    Returns:
        psi: (n_points,) signed wavefunction values
    """
    from pyBall.FireballOCL import Grid as ocl_grid
    from pyBall.FireballOCL.FdataParser import FdataParser

    fparser = FdataParser(os.path.dirname(fdata_basis_dir))
    fparser.parse_info()
    z_list = sorted(set(int(z) for z in atomTypes))
    rcut_map = {int(z): float(rcut) for z in z_list}

    norb_total = sum(norb_per)
    mo_coeffs = C[mo_idx, :].copy() if C.shape[0] > C.shape[1] else C[mo_idx, :].copy()
    if len(mo_coeffs) != norb_total:
        mo_coeffs = C[:, mo_idx].copy()
        if len(mo_coeffs) != norb_total:
            raise ValueError(f"project_orbital_to_points_exp(): MO coefficient size mismatch {len(mo_coeffs)} != {norb_total}")

    norb_per_arr = np.asarray(norb_per, dtype=np.int32)
    coeffs_atoms = remap_coeffs_fortran_to_grid(mo_coeffs, norb_per_arr)

    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([rcut_map[int(z)] for z in atomTypes], dtype=np.float32),
        'type': atomTypes
    }

    projector = ocl_grid.GridProjector(fdata_basis_dir)
    projector.load_basis(z_list)
    psi = projector.project_orbital_points_exp(
        points=points,
        coeffs=coeffs_atoms,
        norb_per=norb_per_arr,
        atoms_dict=atoms_dict,
        beta=float(beta),
        r0=float(r0)
    )
    return psi


def project_orbital_to_points(C, mo_idx, atomTypes, atomPos, orb2atom, norb_per,
                              fdata_basis_dir, points):
    """
    Evaluate a single MO at arbitrary points using OpenCL project_orbital_points kernel.

    This is the most rigorous method for parity testing with Fortran orb2points()
    because it avoids any grid sampling / slicing ambiguity.

    Args:
        C: (norb, nmo) MO eigenvectors (in Fortran order)
        mo_idx: index of MO to project (0-based)
        atomTypes: (natoms,) atomic numbers
        atomPos: (natoms, 3) positions in Å
        orb2atom: (norb,) orbital-to-atom mapping
        norb_per: (natoms,) orbitals per atom
        fdata_basis_dir: path to Fdata/basis for GridProjector
        points: (n_points, 3) array of points to evaluate at

    Returns:
        psi: (n_points,) signed wavefunction values
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
            raise RuntimeError(f"project_orbital_to_points(): no .wf files found for species Z={z}")
        rcut_map[int(z)] = max(ABOHR * fparser.read_wf(f)['rcutoff'] for f in wfs)

    # Extract MO coefficients - use row indexing like test_h2o_orbital_comparison.py
    # C is typically (nmo, norb) from fc.get_wfcoef()
    norb_total = sum(norb_per)
    mo_coeffs = C[mo_idx, :].copy() if C.shape[0] > C.shape[1] else C[mo_idx, :].copy()
    
    # Verify size
    if len(mo_coeffs) != norb_total:
        # Try column indexing as fallback
        mo_coeffs = C[:, mo_idx].copy()
        if len(mo_coeffs) != norb_total:
            raise ValueError(f"project_orbital_to_points(): MO coefficient size mismatch {len(mo_coeffs)} != {norb_total}")

    # Use the correct remapping function [px,py,pz,s] for Grid kernels
    norb_per_arr = np.asarray(norb_per, dtype=np.int32)
    coeffs_atoms = remap_coeffs_fortran_to_grid(mo_coeffs, norb_per_arr)

    # Prepare atoms dict for GridProjector
    atoms_dict = {
        'pos': atomPos,
        'Rcut': np.array([rcut_map[int(z)] for z in atomTypes], dtype=np.float32),
        'type': atomTypes
    }

    # Use pointwise orbital projection kernel
    projector = ocl_grid.GridProjector(fdata_basis_dir)
    projector.load_basis(z_list)
    psi = projector.project_orbital_points(
        points=points,
        coeffs=coeffs_atoms,
        norb_per=norb_per_arr,
        atoms_dict=atoms_dict
    )

    return psi


# =============================================================================
# Orbital Comparison Testing Utilities (from test_h2o_orbital_comparison.py)
# =============================================================================

def parse_orbital_indices(arg_str, max_mo):
    """Parse orbital indices from string like '0,1,2,3,4' or '1-5' or 'all'.

    Args:
        arg_str: String specification of orbitals (e.g., "0,1,2", "1-5", "all")
        max_mo: Maximum number of MOs available

    Returns:
        List of integer MO indices (0-based)

    Examples:
        >>> parse_orbital_indices("0,1,2", 10)
        [0, 1, 2]
        >>> parse_orbital_indices("1-3", 10)
        [1, 2, 3]
        >>> parse_orbital_indices("all", 5)
        [0, 1, 2, 3, 4]
    """
    if arg_str.lower() == 'all':
        return list(range(max_mo))

    indices = []
    parts = arg_str.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    return sorted(set(indices))  # Remove duplicates and sort


def compute_correlation_stats(psi_ref, psi_test):
    """Compute correlation and scaling factors between two orbital arrays.

    Args:
        psi_ref: Reference array (e.g., Fortran orb2points result)
        psi_test: Test array (e.g., OpenCL projection result)

    Returns:
        dict with correlation coefficient and various scaling estimates
    """
    flat_ref = np.asarray(psi_ref).flatten()
    flat_test = np.asarray(psi_test).flatten()

    # Correlation coefficient
    corr = np.corrcoef(flat_ref, flat_test)[0, 1]

    # Linear regression scale: psi_ref ≈ scale * psi_test
    scale_linreg = np.sum(flat_ref * flat_test) / np.sum(flat_test * flat_test)

    # Alternative scaling estimates
    scale_ratio_mean = np.mean(np.abs(flat_ref)) / np.mean(np.abs(flat_test))
    scale_ratio_std = np.std(flat_ref) / np.std(flat_test)
    scale_ratio_max = np.max(np.abs(flat_ref)) / np.max(np.abs(flat_test))

    return {
        'corr': corr,
        'scale_linreg': scale_linreg,
        'scale_ratio_mean': scale_ratio_mean,
        'scale_ratio_std': scale_ratio_std,
        'scale_ratio_max': scale_ratio_max,
        'ref_min': np.min(flat_ref),
        'ref_max': np.max(flat_ref),
        'test_min': np.min(flat_test),
        'test_max': np.max(flat_test)
    }


def write_scaling_summary(results, export_dir, filename='scaling_summary.txt'):
    """Write scaling analysis summary to text file.

    Args:
        results: List of dicts from compute_correlation_stats() for each orbital
        export_dir: Output directory
        filename: Output filename
    """
    filepath = os.path.join(export_dir, filename)

    with open(filepath, 'w') as fsum:
        fsum.write("=" * 80 + "\n")
        fsum.write("SCALING ANALYSIS: Reference vs Test\n")
        fsum.write("=" * 80 + "\n\n")

        fsum.write(f"{'MO':<8} {'Label':<12} {'Corr':>10} {'Scale_LR':>12} {'Scale_Mean':>12} {'Scale_Std':>12} {'Scale_Max':>12}\n")
        fsum.write("-" * 80 + "\n")

        for r in results:
            label = r.get('mo_label', f"MO{r.get('mo_idx', '?')}")
            mo_idx = r.get('mo_idx', -1)
            fsum.write(f"{mo_idx:<8} {label:<12} {r['corr']:>10.4f} {r['scale_linreg']:>12.4f} "
                      f"{r['scale_ratio_mean']:>12.4f} {r['scale_ratio_std']:>12.4f} {r['scale_ratio_max']:>12.4f}\n")

        fsum.write("\n" + "=" * 80 + "\n")
        fsum.write("SCALING FACTOR INTERPRETATION:\n")
        fsum.write("-" * 80 + "\n")
        fsum.write("Scale_LR:    Linear regression scale (psi_ref ≈ scale * psi_test)\n")
        fsum.write("Scale_Mean:  Ratio of mean absolute values\n")
        fsum.write("Scale_Std:   Ratio of standard deviations\n")
        fsum.write("Scale_Max:   Ratio of max absolute values\n")
        fsum.write("\nExpected: If agreement is perfect, all scales ≈ 1.0\n")
        fsum.write("If scales differ consistently, there's a systematic scaling factor.\n")

        # Compute average scale across all orbitals
        avg_scale_linreg = np.mean([r['scale_linreg'] for r in results])
        std_scale_linreg = np.std([r['scale_linreg'] for r in results])
        avg_corr = np.mean([r['corr'] for r in results])

        fsum.write("\n" + "=" * 80 + "\n")
        fsum.write("SUMMARY STATISTICS:\n")
        fsum.write("-" * 80 + "\n")
        fsum.write(f"Average correlation: {avg_corr:.4f}\n")
        fsum.write(f"Average scale (linear reg): {avg_scale_linreg:.4f} ± {std_scale_linreg:.4f}\n")

        # Check for s vs p orbital scaling differences
        s_orbitals = [r for r in results if '_s' in r.get('mo_label', '')]
        p_orbitals = [r for r in results if any(x in r.get('mo_label', '') for x in ['px', 'py', 'pz'])]

        if s_orbitals and p_orbitals:
            avg_s = np.mean([r['scale_linreg'] for r in s_orbitals])
            avg_p = np.mean([r['scale_linreg'] for r in p_orbitals])
            fsum.write(f"\nAverage scale for s-orbitals: {avg_s:.4f}\n")
            fsum.write(f"Average scale for p-orbitals: {avg_p:.4f}\n")
            fsum.write(f"Ratio s/p: {avg_s/avg_p:.4f}\n")

        fsum.write("\n" + "=" * 80 + "\n")

    return filepath


def plot_orbital_comparison(ax, data_xy, atomPos, atomTypes, extent, title, coeffs_text=None, enames=None):
    """Plot orbital with atoms overlaid and symmetric color range.

    IMPORTANT: Both atomPos and extent must be in the SAME coordinate system.
    If extent is absolute coordinates (e.g., from grid_origin), atomPos must be absolute.
    If extent is relative to molecular center, atomPos must also be relative.

    Args:
        ax: matplotlib axes
        data_xy: (nx, ny) 2D orbital data
        atomPos: (natoms, 3) atom positions (absolute or relative, must match extent)
        atomTypes: (natoms,) atom types
        extent: [xmin, xmax, ymin, ymax] for imshow (must match atomPos coordinate system)
        title: Plot title
        coeffs_text: Optional text to display on plot (coefficient info)
        enames: Optional atom names for labels
    """
    vmax = np.max(np.abs(data_xy))
    if vmax < 1e-10:
        vmax = 1e-3
    vmin = -vmax

    # With indexing='ij', data is [x,y]. imshow expects [y,x] (row=y, col=x), so we transpose
    im = ax.imshow(data_xy.T, origin='lower', extent=extent, cmap='bwr', aspect='equal', vmin=vmin, vmax=vmax)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Overlay atoms - atomPos must be in same coordinate system as extent
    ax.scatter(atomPos[:, 0], atomPos[:, 1], c='green', s=100, edgecolors='black', linewidths=1, zorder=10, marker='o')

    # Add atom labels if provided
    if enames is not None:
        for i, (pos, name) in enumerate(zip(atomPos, enames)):
            ax.annotate(f'{name}{i}', (pos[0], pos[1]),
                       textcoords="offset points", xytext=(5, 5), fontsize=8)

    # Add coefficients text if provided
    if coeffs_text is not None:
        ax.text(0.02, 0.98, coeffs_text, transform=ax.transAxes, fontsize=8,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    ax.set_title(title, fontsize=10)


def generate_mo_labels(mo_indices, eigenvalues, nelec, mol_name=None):
    """Generate MO labels (HOMO, LUMO, etc.) based on energy and index.

    Args:
        mo_indices: List of MO indices (0-based)
        eigenvalues: Array of eigenvalues in eV
        nelec: Number of electrons (to determine HOMO)
        mol_name: Optional molecule name for special labels

    Returns:
        List of MO labels
    """
    homo_idx = nelec // 2 - 1  # 0-based HOMO index
    lumo_idx = homo_idx + 1
    nmo = len(eigenvalues)

    labels = []
    for idx in mo_indices:
        # Special labels for core orbitals (only for small molecules)
        if nmo <= 10 and idx == 0 and mol_name == 'H2O':
            labels.append('O1s')
        elif nmo <= 10 and idx == 1 and mol_name == 'H2O':
            labels.append('O2s')
        elif idx == homo_idx:
            labels.append('HOMO')
        elif idx == lumo_idx:
            labels.append('LUMO')
        elif idx < homo_idx:
            labels.append(f'HOMO{idx - homo_idx:+d}')
        elif idx > lumo_idx:
            labels.append(f'LUMO+{idx - lumo_idx}')
        else:
            labels.append(f'MO{idx+1}')

    return labels

