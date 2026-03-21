"""Substrate PLQ field generation for JAX GridFF workflows."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as sp_ndimage

from pyBall import atomicUtils as au
from pyBall import elements

from .array_api import HAS_JAX, as_numpy
from .config import FeatureToggles, GridConfig
from .interfaces import DensityData
from .utils import fractional_to_cartesian
from .hybrid_energy.reactive import DensityReactiveModel

if HAS_JAX:
    import jax
    import jax.numpy as jnp


COULOMB_CONST = 14.3996448915


def _is_nice(n: int, allowed_factors=(2, 3, 5)) -> bool:
    value = int(max(n, 1))
    for factor in allowed_factors:
        while value % factor == 0 and value > 1:
            value //= factor
    return value == 1


def poisson_potential_from_density(rho_zyx, voxel):
    rho_zyx = np.asarray(rho_zyx, dtype=float)
    dx, dy, dz = np.asarray(voxel, dtype=float)
    nz, ny, nx = rho_zyx.shape
    dvol = dx * dy * dz
    charge_density = -rho_zyx / max(dvol, 1.0e-12)
    rho_fft = np.fft.fftn(charge_density)
    kz = np.fft.fftfreq(nz, d=dz) * 2.0 * np.pi
    ky = np.fft.fftfreq(ny, d=dy) * 2.0 * np.pi
    kx = np.fft.fftfreq(nx, d=dx) * 2.0 * np.pi
    kzv, kyv, kxv = np.meshgrid(kz, ky, kx, indexing="ij")
    kernel = 4.0 * np.pi / np.where(kxv * kxv + kyv * kyv + kzv * kzv > 0.0, kxv * kxv + kyv * kyv + kzv * kzv, 1.0)
    kernel[0, 0, 0] = 0.0
    potential = np.fft.ifftn(rho_fft * kernel).real
    return COULOMB_CONST * potential


def bspline_prefilter_numpy(values, wrap_z: bool = False):
    values = np.asarray(values, dtype=float)
    if values.ndim == 3:
        coeff = sp_ndimage.spline_filter1d(values, order=3, axis=2, mode="wrap")
        coeff = sp_ndimage.spline_filter1d(coeff, order=3, axis=1, mode="wrap")
        coeff = sp_ndimage.spline_filter1d(coeff, order=3, axis=0, mode="wrap" if wrap_z else "nearest")
        return coeff
    channels = [bspline_prefilter_numpy(values[..., index], wrap_z=wrap_z) for index in range(values.shape[-1])]
    return np.stack(channels, axis=-1)


def _bspline_basis5_scalar(t: float) -> np.ndarray:
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    t5 = t3 * t2
    return np.array(
        [
            -0.008333333333333333 * t5 + 0.041666666666666664 * t4 - 0.08333333333333333 * t3 + 0.08333333333333333 * t2 - 0.041666666666666664 * t + 0.008333333333333333,
            0.041666666666666664 * t5 - 0.16666666666666666 * t4 + 0.16666666666666666 * t3 + 0.16666666666666666 * t2 - 0.4166666666666667 * t + 0.21666666666666667,
            -0.08333333333333333 * t5 + 0.25 * t4 - 0.5 * t2 + 0.55,
            0.08333333333333333 * t5 - 0.16666666666666666 * t4 - 0.16666666666666666 * t3 + 0.16666666666666666 * t2 + 0.4166666666666667 * t + 0.21666666666666667,
            -0.041666666666666664 * t5 + 0.041666666666666664 * t4 + 0.08333333333333333 * t3 + 0.08333333333333333 * t2 + 0.041666666666666664 * t + 0.008333333333333333,
            0.008333333333333333 * t5,
        ],
        dtype=float,
    )


def _project_atoms_on_grid_quintic_pbc_zyx(
    positions: np.ndarray,
    charges: np.ndarray,
    origin: np.ndarray,
    voxel: np.ndarray,
    grid_shape_xyz: tuple[int, int, int],
) -> np.ndarray:
    nx, ny, nz = (int(v) for v in grid_shape_xyz)
    qgrid = np.zeros((nz, ny, nx), dtype=np.float64)
    shift = np.asarray(voxel, dtype=float)
    for atom_pos, charge in zip(np.asarray(positions, dtype=float), np.asarray(charges, dtype=float)):
        g = (atom_pos + shift - origin) / voxel
        gi = np.floor(g).astype(int)
        t = g - gi.astype(float)
        bx = _bspline_basis5_scalar(float(t[0]))
        by = _bspline_basis5_scalar(float(t[1]))
        bz = _bspline_basis5_scalar(float(t[2]))
        xq = (gi[0] - 2 + np.arange(6, dtype=int)) % nx
        yq = (gi[1] - 2 + np.arange(6, dtype=int)) % ny
        zq = (gi[2] - 2 + np.arange(6, dtype=int)) % nz
        for iz_local in range(6):
            gz = int(zq[iz_local])
            qbz = float(charge) * float(bz[iz_local])
            for iy_local in range(6):
                gy = int(yq[iy_local])
                qbyz = qbz * float(by[iy_local])
                qgrid[gz, gy, xq] += qbyz * bx
    return qgrid


def _poisson_opencl_style(qgrid_zyx: np.ndarray, voxel: np.ndarray) -> np.ndarray:
    nz, ny, nx = qgrid_zyx.shape
    dx, dy, dz = np.asarray(voxel, dtype=float)
    dvol = dx * dy * dz
    nxyz = float(nx * ny * nz)
    rho_k = np.fft.fftn(np.asarray(qgrid_zyx, dtype=np.float64))
    kx = np.fft.fftfreq(nx, d=dx) * (2.0 * np.pi)
    ky = np.fft.fftfreq(ny, d=dy) * (2.0 * np.pi)
    kz = np.fft.fftfreq(nz, d=dz) * (2.0 * np.pi)
    kzv, kyv, kxv = np.meshgrid(kz, ky, kx, indexing="ij")
    k2 = kxv * kxv + kyv * kyv + kzv * kzv
    # FFTW backward transforms used in the current GridFF path are unnormalized.
    # NumPy's ifftn includes a 1/N factor, so the reciprocal-space kernel has to
    # omit that normalization to reproduce the same real-space amplitude.
    scale = COULOMB_CONST * 4.0 * np.pi / max(dvol, 1.0e-30)
    kernel = np.zeros_like(k2, dtype=np.float64)
    mask = k2 > 1.0e-32
    kernel[mask] = scale / k2[mask]
    return np.fft.ifftn(rho_k * kernel).real


def _laplace_real_pbc_numpy(values_zyx: np.ndarray, velocity_zyx: np.ndarray | None, c_sor: float, c_v: float):
    values = np.asarray(values_zyx, dtype=np.float64)
    neighbors = (
        np.roll(values, 1, axis=2)
        + np.roll(values, -1, axis=2)
        + np.roll(values, 1, axis=1)
        + np.roll(values, -1, axis=1)
        + np.roll(values, 1, axis=0)
        + np.roll(values, -1, axis=0)
    ) / 6.0
    updated = neighbors + (neighbors - values) * float(c_sor)
    if velocity_zyx is None:
        return updated, None
    velocity_old = np.asarray(velocity_zyx, dtype=np.float64)
    velocity_new = (updated - values) * float(c_v) + velocity_old * (1.0 - float(c_v))
    return values + velocity_new, velocity_new


def _laplace_real_loop_inert_opencl_style(values_zyx: np.ndarray, niter: int = 2, c_sor: float = 0.0, c_v: float = 0.8) -> np.ndarray:
    v1 = np.asarray(values_zyx, dtype=np.float64)
    velocity = np.zeros_like(v1)
    v2, _ = _laplace_real_pbc_numpy(v1, None, c_sor=c_sor, c_v=0.0)
    last = v2
    for iteration in range(int(niter)):
        if (iteration % 2) == 0:
            v1, velocity = _laplace_real_pbc_numpy(v2, velocity, c_sor=c_sor, c_v=c_v)
            last = v1
        else:
            v2, velocity = _laplace_real_pbc_numpy(v1, velocity, c_sor=c_sor, c_v=c_v)
            last = v2
    return last


def _slab_crop_opencl_style_zyx(
    extended_zyx: np.ndarray,
    original_shape_xyz: tuple[int, int, int],
    voxel: np.ndarray,
    lz_slab: float,
    dipole: float = 0.0,
) -> np.ndarray:
    nx, ny, nz = (int(v) for v in original_shape_xyz)
    dz = float(voxel[2])
    lz = float(nz) * dz
    lz_total = lz + float(lz_slab)
    volume = float(nx) * float(ny) * float(nz) * float(voxel[0]) * float(voxel[1]) * dz
    volume = volume * lz_total / max(lz, 1.0e-30)
    dvcor = 4.0 * np.pi * COULOMB_CONST * float(dipole) / max(volume, 1.0e-30)
    vcor0 = -dvcor * lz_total * 0.5
    ext_nz = int(extended_zyx.shape[0])
    out = np.zeros((nz, ny, nx), dtype=np.float64)
    for iz in range(nz):
        vcor_z = vcor0 + dvcor * (float(iz) * dz)
        jz = ext_nz - iz - 1
        out[iz] = extended_zyx[jz, ::-1, ::-1] + vcor_z
    return out


def coulomb_ewald_slab_opencl_numpy(
    density: DensityData,
    lz_slab: float = 20.0,
    niter: int = 2,
    dipole: float = 0.0,
) -> np.ndarray:
    if density.charges is None:
        return np.zeros(density.grid_shape_zyx, dtype=np.float64)
    nx, ny, nz = density.grid_shape_xyz
    dz = float(density.voxel[2])
    extra_nz = int(np.ceil(float(lz_slab) / max(dz, 1.0e-12)))
    ext_nz = extra_nz + int(nz)
    ext_nz = int(ext_nz)
    while not _is_nice(ext_nz):
        ext_nz += 1
    qgrid_ext = _project_atoms_on_grid_quintic_pbc_zyx(
        positions=density.positions,
        charges=density.charges,
        origin=density.origin,
        voxel=density.voxel,
        grid_shape_xyz=(nx, ny, ext_nz),
    )
    vgrid_ext = _poisson_opencl_style(qgrid_ext, density.voxel)
    vin_ext = _laplace_real_loop_inert_opencl_style(vgrid_ext, niter=niter, c_sor=0.0, c_v=0.8)
    return _slab_crop_opencl_style_zyx(
        extended_zyx=vin_ext,
        original_shape_xyz=(nx, ny, nz),
        voxel=density.voxel,
        lz_slab=float(lz_slab),
        dipole=float(dipole),
    )


def _symbol_to_atomic_number(symbol: str) -> int:
    if symbol not in elements.ELEMENT_DICT:
        raise KeyError(f"Unknown element symbol '{symbol}'")
    return int(elements.ELEMENT_DICT[symbol][0])


def build_req_array(
    symbols: list[str],
    charges: np.ndarray | None,
    element_types_path: str,
    radius_scale: dict[str, float] | None = None,
    energy_scale: dict[str, float] | None = None,
):
    atomic_numbers = np.asarray([_symbol_to_atomic_number(symbol) for symbol in symbols], dtype=int)
    re_vdw = np.asarray(au.getVdWparams(atomic_numbers, fname=element_types_path), dtype=float)
    reqs = np.zeros((len(symbols), 4), dtype=float)
    reqs[:, 0] = re_vdw[:, 0]
    reqs[:, 1] = np.sqrt(re_vdw[:, 1])
    radius_scale = radius_scale or {}
    energy_scale = energy_scale or {}
    for index, symbol in enumerate(symbols):
        if symbol in radius_scale:
            reqs[index, 0] *= float(radius_scale[symbol])
        if symbol in energy_scale:
            reqs[index, 1] *= np.sqrt(float(energy_scale[symbol]))
    if charges is not None:
        reqs[:, 2] = np.asarray(charges, dtype=float)
    return reqs


def auto_pbc_counts(cell: np.ndarray, rcut: float, mask=(1, 1, 0)) -> tuple[int, int, int]:
    counts = [0, 0, 0]
    cell = np.asarray(cell, dtype=float)
    for index in range(3):
        if mask[index]:
            length = float(np.linalg.norm(cell[index]))
            counts[index] = int(rcut / max(length, 1.0e-12)) + 1
    return tuple(counts)


def make_pbc_shifts(cell: np.ndarray, counts: tuple[int, int, int]) -> np.ndarray:
    cell = np.asarray(cell, dtype=float)
    shifts = []
    for ia in range(-counts[0], counts[0] + 1):
        for ib in range(-counts[1], counts[1] + 1):
            for ic in range(-counts[2], counts[2] + 1):
                shifts.append((ia * cell[0]) + (ib * cell[1]) + (ic * cell[2]))
    return np.asarray(shifts, dtype=float)


def iter_grid_point_chunks(density: DensityData, chunk_size: int):
    nx, ny, nz = density.grid_shape_xyz
    total = nx * ny * nz
    xy = nx * ny
    for start in range(0, total, chunk_size):
        stop = min(total, start + chunk_size)
        linear = np.arange(start, stop, dtype=np.int64)
        iz = linear // xy
        rem = linear - iz * xy
        iy = rem // nx
        ix = rem - iy * nx
        frac = np.stack(
            [
                ix.astype(float) / max(nx, 1),
                iy.astype(float) / max(ny, 1),
                iz.astype(float) / max(nz, 1),
            ],
            axis=1,
        )
        points = density.origin[None, :] + fractional_to_cartesian(frac, density.cell)
        yield start, stop, iz.astype(np.int32), iy.astype(np.int32), ix.astype(np.int32), points


def _pairwise_chunk_numpy(points, atom_positions, reqs, alpha_morse: float):
    dp = points[:, None, :] - atom_positions[None, :, :]
    r2 = np.sum(dp * dp, axis=-1) + 1.0e-32
    r = np.sqrt(r2)
    e = np.exp((-alpha_morse) * (r - reqs[None, :, 0]))
    e_m = e * reqs[None, :, 1]
    pauli = np.sum(e_m * e, axis=1)
    london = np.sum(e_m * -2.0, axis=1)
    coulomb = np.sum(COULOMB_CONST * reqs[None, :, 2] / r, axis=1)
    return pauli, london, coulomb


if HAS_JAX:

    @jax.jit
    def _pairwise_chunk_jax(points, atom_positions, reqs, alpha_morse):
        dp = points[:, None, :] - atom_positions[None, :, :]
        r2 = jnp.sum(dp * dp, axis=-1) + 1.0e-32
        r = jnp.sqrt(r2)
        e = jnp.exp((-alpha_morse) * (r - reqs[None, :, 0]))
        e_m = e * reqs[None, :, 1]
        pauli = jnp.sum(e_m * e, axis=1)
        london = jnp.sum(e_m * -2.0, axis=1)
        coulomb = jnp.sum(COULOMB_CONST * reqs[None, :, 2] / r, axis=1)
        return pauli, london, coulomb


def accumulate_pairwise_fields(
    density: DensityData,
    reqs: np.ndarray,
    alpha_morse: float,
    point_chunk: int,
    atom_chunk: int,
    morse_counts: tuple[int, int, int],
    coulomb_counts: tuple[int, int, int] | None,
    prefer_jax: bool,
):
    nx, ny, nz = density.grid_shape_xyz
    pauli = np.zeros((nz, ny, nx), dtype=float)
    london = np.zeros((nz, ny, nx), dtype=float)
    coulomb = np.zeros((nz, ny, nx), dtype=float)

    morse_shifts = make_pbc_shifts(density.cell, morse_counts)
    morse_positions = (density.positions[None, :, :] + morse_shifts[:, None, :]).reshape(-1, 3)
    morse_reqs = np.repeat(reqs[None, :, :], morse_shifts.shape[0], axis=0).reshape(-1, reqs.shape[1])

    if coulomb_counts is None:
        coulomb_positions = None
        coulomb_reqs = None
    else:
        coulomb_shifts = make_pbc_shifts(density.cell, coulomb_counts)
        coulomb_positions = (density.positions[None, :, :] + coulomb_shifts[:, None, :]).reshape(-1, 3)
        coulomb_reqs = np.repeat(reqs[None, :, :], coulomb_shifts.shape[0], axis=0).reshape(-1, reqs.shape[1])

    use_jax = bool(prefer_jax and HAS_JAX)
    alpha_value = float(alpha_morse)
    for _, _, iz, iy, ix, points in iter_grid_point_chunks(density, point_chunk):
        chunk_pauli = np.zeros(points.shape[0], dtype=float)
        chunk_london = np.zeros(points.shape[0], dtype=float)
        chunk_coulomb = np.zeros(points.shape[0], dtype=float)

        for start in range(0, morse_positions.shape[0], atom_chunk):
            stop = min(start + atom_chunk, morse_positions.shape[0])
            if use_jax:
                p_val, l_val, _ = _pairwise_chunk_jax(
                    jnp.asarray(points, dtype=jnp.float64),
                    jnp.asarray(morse_positions[start:stop], dtype=jnp.float64),
                    jnp.asarray(morse_reqs[start:stop], dtype=jnp.float64),
                    jnp.asarray(alpha_value, dtype=jnp.float64),
                )
                chunk_pauli += as_numpy(p_val)
                chunk_london += as_numpy(l_val)
            else:
                p_val, l_val, _ = _pairwise_chunk_numpy(points, morse_positions[start:stop], morse_reqs[start:stop], alpha_value)
                chunk_pauli += p_val
                chunk_london += l_val

        if coulomb_positions is not None:
            for start in range(0, coulomb_positions.shape[0], atom_chunk):
                stop = min(start + atom_chunk, coulomb_positions.shape[0])
                if use_jax:
                    _, _, c_val = _pairwise_chunk_jax(
                        jnp.asarray(points, dtype=jnp.float64),
                        jnp.asarray(coulomb_positions[start:stop], dtype=jnp.float64),
                        jnp.asarray(coulomb_reqs[start:stop], dtype=jnp.float64),
                        jnp.asarray(alpha_value, dtype=jnp.float64),
                    )
                    chunk_coulomb += as_numpy(c_val)
                else:
                    _, _, c_val = _pairwise_chunk_numpy(points, coulomb_positions[start:stop], coulomb_reqs[start:stop], alpha_value)
                    chunk_coulomb += c_val

        pauli[iz, iy, ix] = chunk_pauli
        london[iz, iy, ix] = chunk_london
        coulomb[iz, iy, ix] = chunk_coulomb
    return pauli, london, coulomb


def _reactive_channels(density: DensityData, reactive_elements: tuple[str, ...], grid_config: GridConfig, enabled: bool):
    nz, ny, nx = density.grid_shape_zyx
    if (not enabled) or (density.rho_zyx is None) or (len(reactive_elements) == 0):
        return np.zeros((nz, ny, nx, len(reactive_elements)), dtype=float)
    z_ref = float(np.max(density.positions[:, 2]))
    reactive_model = DensityReactiveModel(
        z_ref=z_ref,
        voxel=density.voxel,
        power=grid_config.reactive_power,
        sigma_z=grid_config.reactive_sigma_z,
    )
    return reactive_model.build_channels(density, reactive_elements)


def _surrogate_density_fields(density: DensityData, grid_config: GridConfig):
    if density.rho_zyx is None:
        raise ValueError("surrogate_density builder requires rho_zyx")
    rho = np.asarray(density.rho_zyx, dtype=float)
    rho_norm = rho / max(float(np.max(rho)), 1.0e-12)
    pauli = np.power(np.clip(rho_norm, 0.0, None), grid_config.pauli_power)
    london = -np.power(np.clip(rho_norm, 0.0, None), grid_config.london_power)
    if density.v_loc_zyx is not None:
        coulomb = np.asarray(density.v_loc_zyx, dtype=float)
        source = "locpot"
    else:
        coulomb = poisson_potential_from_density(rho, density.voxel)
        source = "poisson_from_rho"
    return pauli, london, coulomb, source


def _metal_density_fields(density: DensityData, grid_config: GridConfig):
    """Volumetric P/L from CHGCAR + LOCPOT for Q. Physically correct for metals."""
    if density.rho_zyx is None:
        raise ValueError("metal_density_plq builder requires rho_zyx from CHGCAR")
    rho = np.clip(np.asarray(density.rho_zyx, dtype=float), 0.0, None)
    if grid_config.metal_density_rho_smoothing_sigma > 0.0:
        rho = sp_ndimage.gaussian_filter(rho, grid_config.metal_density_rho_smoothing_sigma)
    bulk_rho = float(grid_config.metal_bulk_electron_density)
    rho_scale = max(bulk_rho, float(np.max(rho)), 1.0e-30)
    rho_norm = np.clip(rho / rho_scale, 0.0, 1.0)
    pauli = np.power(rho_norm, grid_config.metal_density_pauli_power)
    london = -np.power(rho_norm, grid_config.metal_density_london_power)
    if grid_config.london_damping_d0 > 0.0:
        z_surface = float(np.max(density.positions[:, 2]))
        nz = london.shape[0]
        voxel_z = float(density.voxel[2, 2]) if density.voxel.ndim == 2 else float(density.voxel[2])
        origin_z = float(density.origin[2])
        z_grid = np.arange(nz) * voxel_z + origin_z
        d_from_surface = z_grid - z_surface
        w = max(float(grid_config.london_damping_width), 1.0e-6)
        fermi = 1.0 / (1.0 + np.exp((d_from_surface - grid_config.london_damping_d0) / w))
        london *= fermi[:, np.newaxis, np.newaxis]
    if density.v_loc_zyx is not None:
        coulomb = np.asarray(density.v_loc_zyx, dtype=float)
        source = "locpot"
    else:
        coulomb = poisson_potential_from_density(rho, density.voxel)
        source = "poisson_from_rho"
    return pauli, london, coulomb, source


def build_substrate_grids(
    density: DensityData,
    reactive_elements: tuple[str, ...],
    grid_config: GridConfig | None = None,
    toggles: FeatureToggles | None = None,
    prefer_jax: bool = True,
):
    grid_config = grid_config or GridConfig()
    toggles = toggles or FeatureToggles()

    if grid_config.builder_mode == "surrogate_density":
        pauli, london, coulomb, coulomb_source = _surrogate_density_fields(density, grid_config)
    elif grid_config.builder_mode == "metal_density_plq":
        pauli, london, coulomb, coulomb_source = _metal_density_fields(density, grid_config)
    else:
        reqs = build_req_array(
            symbols=density.symbols,
            charges=density.charges,
            element_types_path=grid_config.element_types_path,
            radius_scale=grid_config.req_scale_radius,
            energy_scale=grid_config.req_scale_energy,
        )
        morse_counts = auto_pbc_counts(density.cell, grid_config.morse_rcut, mask=(1, 1, 0))
        use_pairwise_coulomb = not (grid_config.builder_mode == "metal_dft_plq" and toggles.use_locpot and density.v_loc_zyx is not None)
        use_opencl_style_ewald = bool(grid_config.builder_mode == "parity_core" and density.charges is not None)
        if use_opencl_style_ewald:
            use_pairwise_coulomb = False
        coulomb_counts = auto_pbc_counts(density.cell, grid_config.coulomb_rcut, mask=(1, 1, 0)) if use_pairwise_coulomb else None
        pauli, london, pairwise_coulomb = accumulate_pairwise_fields(
            density=density,
            reqs=reqs,
            alpha_morse=grid_config.alpha_morse,
            point_chunk=int(grid_config.pairwise_point_chunk),
            atom_chunk=int(grid_config.pairwise_atom_chunk),
            morse_counts=morse_counts,
            coulomb_counts=coulomb_counts,
            prefer_jax=prefer_jax,
        )
        if grid_config.builder_mode == "metal_dft_plq" and toggles.use_locpot and density.v_loc_zyx is not None:
            coulomb = np.asarray(density.v_loc_zyx, dtype=float)
            coulomb_source = "locpot"
        elif use_opencl_style_ewald:
            coulomb = coulomb_ewald_slab_opencl_numpy(density, lz_slab=20.0, niter=2, dipole=0.0)
            coulomb_source = "ewald_slab_quintic_numpy"
        elif np.any(np.abs(reqs[:, 2]) > 0.0):
            coulomb = pairwise_coulomb
            coulomb_source = "pairwise_direct"
        elif density.rho_zyx is not None:
            coulomb = poisson_potential_from_density(density.rho_zyx, density.voxel)
            coulomb_source = "poisson_from_rho"
        else:
            coulomb = np.zeros_like(pauli)
            coulomb_source = "zeros"

    reactive = _reactive_channels(
        density=density,
        reactive_elements=reactive_elements,
        grid_config=grid_config,
        enabled=bool(toggles.use_reactive_grid),
    )
    wrap_z = bool(grid_config.builder_mode == "parity_core")
    pauli_coeff = bspline_prefilter_numpy(pauli, wrap_z=wrap_z)
    london_coeff = bspline_prefilter_numpy(london, wrap_z=wrap_z)
    coulomb_coeff = bspline_prefilter_numpy(coulomb, wrap_z=wrap_z)
    reactive_coeff = bspline_prefilter_numpy(reactive, wrap_z=wrap_z)
    z_ref = float(np.max(density.positions[:, 2]))
    metadata = {
        "builder_mode": grid_config.builder_mode,
        "coulomb_source": coulomb_source,
        "reactive_elements": list(reactive_elements),
        "z_ref": z_ref,
        "interpolation_order": int(grid_config.interpolation_order),
        "alpha_morse": float(grid_config.alpha_morse),
        "morse_rcut": float(grid_config.morse_rcut),
        "coulomb_rcut": float(grid_config.coulomb_rcut),
        "req_scale_radius": dict(grid_config.req_scale_radius),
        "req_scale_energy": dict(grid_config.req_scale_energy),
    }
    if grid_config.builder_mode == "metal_density_plq":
        metadata["metal_density_pauli_power"] = float(grid_config.metal_density_pauli_power)
        metadata["metal_density_london_power"] = float(grid_config.metal_density_london_power)
        metadata["metal_bulk_electron_density"] = float(grid_config.metal_bulk_electron_density)
        metadata["london_damping_d0"] = float(grid_config.london_damping_d0)
        metadata["london_damping_width"] = float(grid_config.london_damping_width)
    return {
        "pauli_zyx": pauli,
        "london_zyx": london,
        "coulomb_zyx": coulomb,
        "reactive_zyxc": reactive,
        "pauli_coeff_zyx": pauli_coeff,
        "london_coeff_zyx": london_coeff,
        "coulomb_coeff_zyx": coulomb_coeff,
        "reactive_coeff_zyxc": reactive_coeff,
        "metadata": metadata,
    }
