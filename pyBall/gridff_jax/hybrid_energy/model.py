"""Hybrid GridFF student model with JAX-first evaluation and NumPy fallback."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from pyBall import atomicUtils as au
from pyBall import elements

from ..array_api import HAS_JAX, as_numpy, finite_difference_force
from ..config import FeatureToggles, GridConfig, HybridModelConfig
from ..interfaces import DensityData, ModelEvaluation, PoseBatch
from ..substrate_fields import build_substrate_grids, poisson_potential_from_density
from ..utils import cartesian_to_fractional
from .qeq import solve_qeq_with_reservoir

if HAS_JAX:
    import jax
    import jax.numpy as jnp

    from .qeq import solve_qeq_with_reservoir_jax


DEFAULT_QEQ = {
    "H": (-4.528, 13.8904),
    "C": (-5.343, 10.126),
    "N": (-6.899, 11.760),
    "O": (-8.741, 13.364),
    "Cu": (-4.48, 8.50),
    "Ag": (-4.44, 7.58),
    "Au": (-5.10, 8.95),
}


@dataclass
class HybridParameters:
    pauli: dict[str, float] = field(default_factory=dict)
    london: dict[str, float] = field(default_factory=dict)
    reactive: dict[str, float] = field(default_factory=dict)
    static_charge: dict[str, float] = field(default_factory=dict)
    req_radius_offset: dict[str, float] = field(default_factory=dict)
    req_energy_scale: dict[str, float] = field(default_factory=dict)
    chi: dict[str, float] = field(default_factory=dict)
    hardness: dict[str, float] = field(default_factory=dict)
    image_scale: float = 1.0
    image_plane: float = 0.0
    image_damping: float = 0.5
    sample_shift_z: float = 0.0
    coulomb_shift_z: float = 0.0
    reservoir_chi: float = -4.5
    reservoir_hardness: float = 7.0
    total_charge: float = 0.0


def default_hybrid_parameters(elements, z_image: float):
    params = HybridParameters(image_plane=float(z_image))
    for symbol in elements:
        chi, hardness = DEFAULT_QEQ.get(symbol, (-4.0, 10.0))
        params.pauli[symbol] = 1.0
        params.london[symbol] = 1.0
        params.reactive[symbol] = 1.0
        params.static_charge[symbol] = 0.0
        params.req_radius_offset[symbol] = 0.0
        params.req_energy_scale[symbol] = 1.0
        params.chi[symbol] = chi
        params.hardness[symbol] = hardness
    return params


def _symbol_to_atomic_number(symbol: str) -> int:
    if symbol not in elements.ELEMENT_DICT:
        raise KeyError(f"Unknown element symbol '{symbol}'")
    return int(elements.ELEMENT_DICT[symbol][0])


def _base_req_for_elements(symbols, element_types_path: str):
    atomic_numbers = np.asarray([_symbol_to_atomic_number(symbol) for symbol in symbols], dtype=int)
    re_vdw = np.asarray(au.getVdWparams(atomic_numbers, fname=element_types_path), dtype=float)
    radius = np.asarray(re_vdw[:, 0], dtype=float)
    sqrt_energy = np.sqrt(np.clip(np.asarray(re_vdw[:, 1], dtype=float), 0.0, None))
    return radius, sqrt_energy


def build_physics_grids(
    density: DensityData,
    reactive_elements: tuple[str, ...],
    grid_config: GridConfig | None = None,
    toggles: FeatureToggles | None = None,
    prefer_jax: bool = True,
):
    return build_substrate_grids(
        density=density,
        reactive_elements=reactive_elements,
        grid_config=grid_config,
        toggles=toggles,
        prefer_jax=prefer_jax,
    )


def _bspline_basis_numpy(u):
    inv6 = 1.0 / 6.0
    u2 = u * u
    t = 1.0 - u
    return np.stack(
        [
            inv6 * t * t * t,
            inv6 * (3.0 * u2 * (u - 2.0) + 4.0),
            inv6 * (3.0 * u * (1.0 + u - u2) + 1.0),
            inv6 * u2 * u,
        ],
        axis=-1,
    )


def _sample_trilinear_numpy(values, density: DensityData, positions):
    values = np.asarray(values, dtype=float)
    positions = np.asarray(positions, dtype=float)
    frac = cartesian_to_fractional(positions - density.origin[None, :], density.cell)
    frac[:, 0:2] = frac[:, 0:2] - np.floor(frac[:, 0:2])
    nz, ny, nx = values.shape[:3]
    u = frac[:, 0] * nx
    v = frac[:, 1] * ny
    w = np.clip(frac[:, 2] * nz, 0.0, nz - 1.000001)

    ix0 = np.floor(u).astype(np.int32) % nx
    iy0 = np.floor(v).astype(np.int32) % ny
    iz0 = np.floor(w).astype(np.int32)
    ix1 = (ix0 + 1) % nx
    iy1 = (iy0 + 1) % ny
    iz1 = np.minimum(iz0 + 1, nz - 1)
    tx = (u - np.floor(u))[:, None]
    ty = (v - np.floor(v))[:, None]
    tz = (w - np.floor(w))[:, None]

    if values.ndim == 3:
        values = values[..., None]
    c000 = values[iz0, iy0, ix0]
    c001 = values[iz0, iy0, ix1]
    c010 = values[iz0, iy1, ix0]
    c011 = values[iz0, iy1, ix1]
    c100 = values[iz1, iy0, ix0]
    c101 = values[iz1, iy0, ix1]
    c110 = values[iz1, iy1, ix0]
    c111 = values[iz1, iy1, ix1]
    c00 = c000 * (1.0 - tx) + c001 * tx
    c01 = c010 * (1.0 - tx) + c011 * tx
    c10 = c100 * (1.0 - tx) + c101 * tx
    c11 = c110 * (1.0 - tx) + c111 * tx
    c0 = c00 * (1.0 - ty) + c01 * ty
    c1 = c10 * (1.0 - ty) + c11 * ty
    out = c0 * (1.0 - tz) + c1 * tz
    if out.shape[1] == 1:
        return out[:, 0]
    return out


def _sample_bspline_numpy(values, density: DensityData, positions):
    values = np.asarray(values, dtype=float)
    positions = np.asarray(positions, dtype=float)
    frac = cartesian_to_fractional(positions - density.origin[None, :], density.cell)
    frac[:, 0:2] = frac[:, 0:2] - np.floor(frac[:, 0:2])
    if values.ndim == 3:
        values = values[..., None]
    nz, ny, nx = values.shape[:3]
    ux = frac[:, 0] * nx
    uy = frac[:, 1] * ny
    uz = frac[:, 2] * nz
    ix = np.floor(ux).astype(np.int32)
    iy = np.floor(uy).astype(np.int32)
    iz = np.floor(uz).astype(np.int32)
    tx = ux - ix
    ty = uy - iy
    tz = uz - iz
    inside = (iz >= 1) & (iz < (nz - 2))
    bix = _bspline_basis_numpy(tx)
    biy = _bspline_basis_numpy(ty)
    biz = _bspline_basis_numpy(tz)
    xinds = np.stack([ix - 1, ix, ix + 1, ix + 2], axis=1) % nx
    yinds = np.stack([iy - 1, iy, iy + 1, iy + 2], axis=1) % ny
    zinds = np.clip(np.stack([iz - 1, iz, iz + 1, iz + 2], axis=1), 0, nz - 1)
    out = np.zeros((positions.shape[0], values.shape[-1]), dtype=float)
    for ia in range(4):
        for ib in range(4):
            for ic in range(4):
                weight = (bix[:, ia] * biy[:, ib] * biz[:, ic])[:, None]
                out += weight * values[zinds[:, ic], yinds[:, ib], xinds[:, ia]]
    out[~inside] = 0.0
    if out.shape[1] == 1:
        return out[:, 0]
    return out


if HAS_JAX:

    def _bspline_basis_jax(u):
        inv6 = jnp.asarray(1.0 / 6.0, dtype=jnp.float64)
        u2 = u * u
        t = 1.0 - u
        return jnp.stack(
            [
                inv6 * t * t * t,
                inv6 * (3.0 * u2 * (u - 2.0) + 4.0),
                inv6 * (3.0 * u * (1.0 + u - u2) + 1.0),
                inv6 * u2 * u,
            ],
            axis=-1,
        )

    def _sample_trilinear_jax(values, positions, origin, cell_inv):
        values = jnp.asarray(values, dtype=jnp.float64)
        positions = jnp.asarray(positions, dtype=jnp.float64)
        frac = (positions - origin[None, :]) @ cell_inv
        frac = frac.at[:, 0:2].set(frac[:, 0:2] - jnp.floor(frac[:, 0:2]))
        nz, ny, nx = values.shape[:3]
        u = frac[:, 0] * nx
        v = frac[:, 1] * ny
        w = jnp.clip(frac[:, 2] * nz, 0.0, nz - 1.000001)

        ix0 = jnp.floor(u).astype(jnp.int32) % nx
        iy0 = jnp.floor(v).astype(jnp.int32) % ny
        iz0 = jnp.floor(w).astype(jnp.int32)
        ix1 = (ix0 + 1) % nx
        iy1 = (iy0 + 1) % ny
        iz1 = jnp.minimum(iz0 + 1, nz - 1)
        tx = (u - jnp.floor(u))[:, None]
        ty = (v - jnp.floor(v))[:, None]
        tz = (w - jnp.floor(w))[:, None]

        if values.ndim == 3:
            values = values[..., None]
        c000 = values[iz0, iy0, ix0]
        c001 = values[iz0, iy0, ix1]
        c010 = values[iz0, iy1, ix0]
        c011 = values[iz0, iy1, ix1]
        c100 = values[iz1, iy0, ix0]
        c101 = values[iz1, iy0, ix1]
        c110 = values[iz1, iy1, ix0]
        c111 = values[iz1, iy1, ix1]
        c00 = c000 * (1.0 - tx) + c001 * tx
        c01 = c010 * (1.0 - tx) + c011 * tx
        c10 = c100 * (1.0 - tx) + c101 * tx
        c11 = c110 * (1.0 - tx) + c111 * tx
        c0 = c00 * (1.0 - ty) + c01 * ty
        c1 = c10 * (1.0 - ty) + c11 * ty
        out = c0 * (1.0 - tz) + c1 * tz
        if out.shape[-1] == 1:
            return out[:, 0]
        return out

    def _sample_bspline_jax(values, positions, origin, cell_inv):
        values = jnp.asarray(values, dtype=jnp.float64)
        positions = jnp.asarray(positions, dtype=jnp.float64)
        frac = (positions - origin[None, :]) @ cell_inv
        frac = frac.at[:, 0:2].set(frac[:, 0:2] - jnp.floor(frac[:, 0:2]))
        if values.ndim == 3:
            values = values[..., None]
        nz, ny, nx = values.shape[:3]
        ux = frac[:, 0] * nx
        uy = frac[:, 1] * ny
        uz = frac[:, 2] * nz
        ix = jnp.floor(ux).astype(jnp.int32)
        iy = jnp.floor(uy).astype(jnp.int32)
        iz = jnp.floor(uz).astype(jnp.int32)
        tx = ux - ix
        ty = uy - iy
        tz = uz - iz
        inside = (iz >= 1) & (iz < (nz - 2))
        bix = _bspline_basis_jax(tx)
        biy = _bspline_basis_jax(ty)
        biz = _bspline_basis_jax(tz)
        xinds = jnp.mod(jnp.stack([ix - 1, ix, ix + 1, ix + 2], axis=1), nx)
        yinds = jnp.mod(jnp.stack([iy - 1, iy, iy + 1, iy + 2], axis=1), ny)
        zinds = jnp.clip(jnp.stack([iz - 1, iz, iz + 1, iz + 2], axis=1), 0, nz - 1)
        out = jnp.zeros((positions.shape[0], values.shape[-1]), dtype=jnp.float64)
        for ia in range(4):
            for ib in range(4):
                for ic in range(4):
                    weight = (bix[:, ia] * biy[:, ib] * biz[:, ic])[:, None]
                    out = out + weight * values[zinds[:, ic], yinds[:, ib], xinds[:, ia]]
        out = jnp.where(inside[:, None], out, 0.0)
        if out.shape[-1] == 1:
            return out[:, 0]
        return out


class HybridGridFFModel:
    def __init__(
        self,
        density: DensityData,
        reactive_elements: tuple[str, ...],
        toggles: FeatureToggles | None = None,
        grid_config: GridConfig | None = None,
        hybrid_config: HybridModelConfig | None = None,
        prefer_jax: bool = True,
    ):
        self.density = density
        self.reactive_elements = tuple(reactive_elements)
        self.parameter_elements = tuple(reactive_elements)
        self.toggles = toggles or FeatureToggles()
        self.grid_config = grid_config or GridConfig()
        self.hybrid_config = hybrid_config or HybridModelConfig(reactive_channels=self.reactive_elements)
        self.use_jax = bool(prefer_jax and HAS_JAX)
        self.grids = build_physics_grids(
            density=density,
            reactive_elements=self.reactive_elements,
            grid_config=self.grid_config,
            toggles=self.toggles,
            prefer_jax=self.use_jax,
        )
        if self.grid_config.interpolation_order >= 3:
            self.interpolation_kind = "bspline_cubic"
            self.sample_grids = {
                "pauli_zyx": self.grids["pauli_coeff_zyx"],
                "london_zyx": self.grids["london_coeff_zyx"],
                "coulomb_zyx": self.grids["coulomb_coeff_zyx"],
                "reactive_zyxc": self.grids["reactive_coeff_zyxc"],
            }
        else:
            self.interpolation_kind = "trilinear"
        self.sample_grids = {
            "pauli_zyx": self.grids["pauli_zyx"],
            "london_zyx": self.grids["london_zyx"],
            "coulomb_zyx": self.grids["coulomb_zyx"],
            "reactive_zyxc": self.grids["reactive_zyxc"],
        }
        self.element_to_index = {element: i for i, element in enumerate(self.parameter_elements)}
        self.base_req_radius, self.base_req_sqrt_energy = _base_req_for_elements(
            self.parameter_elements,
            self.grid_config.element_types_path,
        )
        self._jax_predictor_cache = {}

        if HAS_JAX:
            self._jax_origin = jnp.asarray(np.asarray(density.origin, dtype=float), dtype=jnp.float64)
            self._jax_cell_inv = jnp.asarray(np.linalg.inv(np.asarray(density.cell, dtype=float)), dtype=jnp.float64)
            self._jax_base_req_radius = jnp.asarray(self.base_req_radius, dtype=jnp.float64)
            self._jax_base_req_sqrt_energy = jnp.asarray(self.base_req_sqrt_energy, dtype=jnp.float64)
            self._jax_sample_grids = {
                "pauli_zyx": jnp.asarray(np.asarray(self.sample_grids["pauli_zyx"], dtype=float), dtype=jnp.float64),
                "london_zyx": jnp.asarray(np.asarray(self.sample_grids["london_zyx"], dtype=float), dtype=jnp.float64),
                "coulomb_zyx": jnp.asarray(np.asarray(self.sample_grids["coulomb_zyx"], dtype=float), dtype=jnp.float64),
                "reactive_zyxc": jnp.asarray(np.asarray(self.sample_grids["reactive_zyxc"], dtype=float), dtype=jnp.float64),
            }

    def make_species_indices(self, symbols):
        return np.asarray([self.element_to_index[symbol] for symbol in symbols], dtype=np.int32)

    def parameter_arrays(self, params: HybridParameters, use_jax: bool | None = None):
        use_jax = self.use_jax if use_jax is None else bool(use_jax and HAS_JAX)
        xp = jnp if use_jax else np
        arrays = {
            "pauli": xp.asarray([params.pauli.get(symbol, 1.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "london": xp.asarray([params.london.get(symbol, 1.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "reactive": xp.asarray([params.reactive.get(symbol, 0.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "static_charge": xp.asarray([params.static_charge.get(symbol, 0.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "req_radius_offset": xp.asarray([params.req_radius_offset.get(symbol, 0.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "req_energy_scale": xp.asarray([params.req_energy_scale.get(symbol, 1.0) for symbol in self.parameter_elements], dtype=xp.float64),
            "chi": xp.asarray([params.chi.get(symbol, DEFAULT_QEQ.get(symbol, (-4.0, 10.0))[0]) for symbol in self.parameter_elements], dtype=xp.float64),
            "hardness": xp.asarray([params.hardness.get(symbol, DEFAULT_QEQ.get(symbol, (-4.0, 10.0))[1]) for symbol in self.parameter_elements], dtype=xp.float64),
            "image_scale": xp.asarray(params.image_scale, dtype=xp.float64),
            "image_plane": xp.asarray(params.image_plane, dtype=xp.float64),
            "image_damping": xp.asarray(params.image_damping, dtype=xp.float64),
            "sample_shift_z": xp.asarray(params.sample_shift_z, dtype=xp.float64),
            "coulomb_shift_z": xp.asarray(params.coulomb_shift_z, dtype=xp.float64),
            "reservoir_chi": xp.asarray(params.reservoir_chi, dtype=xp.float64),
            "reservoir_hardness": xp.asarray(params.reservoir_hardness, dtype=xp.float64),
            "total_charge": xp.asarray(params.total_charge, dtype=xp.float64),
        }
        return arrays

    def _atom_params(self, params: HybridParameters, symbols):
        return {
            "pauli": np.array([params.pauli.get(symbol, 1.0) for symbol in symbols], dtype=float),
            "london": np.array([params.london.get(symbol, 1.0) for symbol in symbols], dtype=float),
            "reactive": np.array([params.reactive.get(symbol, 0.0) for symbol in symbols], dtype=float),
            "static_charge": np.array([params.static_charge.get(symbol, 0.0) for symbol in symbols], dtype=float),
            "req_radius_offset": np.array([params.req_radius_offset.get(symbol, 0.0) for symbol in symbols], dtype=float),
            "req_energy_scale": np.array([params.req_energy_scale.get(symbol, 1.0) for symbol in symbols], dtype=float),
            "chi": np.array([params.chi.get(symbol, DEFAULT_QEQ.get(symbol, (-4.0, 10.0))[0]) for symbol in symbols], dtype=float),
            "hardness": np.array([params.hardness.get(symbol, DEFAULT_QEQ.get(symbol, (-4.0, 10.0))[1]) for symbol in symbols], dtype=float),
            "species_idx": np.array([self.element_to_index[symbol] for symbol in symbols], dtype=np.int32),
        }

    def _sample_field_numpy(self, values, positions):
        if self.grid_config.interpolation_order >= 3:
            return _sample_bspline_numpy(values, self.density, positions)
        return _sample_trilinear_numpy(values, self.density, positions)

    def _shift_positions_numpy(self, positions, z_shift: float):
        shifted = np.asarray(positions, dtype=float).copy()
        shifted[:, 2] += float(z_shift)
        return shifted

    def _static_charges_numpy(self, params: HybridParameters, symbols, fixed_charges):
        scale = np.array([params.static_charge.get(symbol, 0.0) for symbol in symbols], dtype=float)
        if fixed_charges is not None and np.max(np.abs(fixed_charges)) > 0.0:
            return np.asarray(fixed_charges, dtype=float) * (1.0 + scale)
        return scale

    def _molecule_plq_numpy(self, atom_params):
        if not self.hybrid_config.use_req_plq:
            return atom_params["pauli"], atom_params["london"]
        species_idx = np.asarray(atom_params["species_idx"], dtype=np.int32)
        radius = self.base_req_radius[species_idx] + atom_params["req_radius_offset"]
        sqrt_energy = self.base_req_sqrt_energy[species_idx] * atom_params["req_energy_scale"]
        expo = np.exp(float(self.grid_config.alpha_morse) * radius)
        london = expo * sqrt_energy
        pauli = expo * london
        return pauli, london

    def energy_single(self, positions, symbols, params: HybridParameters, fixed_charges=None):
        atom_params = self._atom_params(params, symbols)
        pauli_coeff, london_coeff = self._molecule_plq_numpy(atom_params)
        positions_pl = self._shift_positions_numpy(positions, params.sample_shift_z)
        positions_q = self._shift_positions_numpy(positions_pl, params.coulomb_shift_z)
        pauli_sample = self._sample_field_numpy(self.sample_grids["pauli_zyx"], positions_pl)
        london_sample = self._sample_field_numpy(self.sample_grids["london_zyx"], positions_pl)
        coulomb_sample = self._sample_field_numpy(self.sample_grids["coulomb_zyx"], positions_q)
        reactive_sample = self._sample_field_numpy(self.sample_grids["reactive_zyxc"], positions_pl)
        if np.ndim(reactive_sample) == 1:
            reactive_selected = reactive_sample
        else:
            reactive_selected = reactive_sample[np.arange(len(symbols)), atom_params["species_idx"]]

        coulomb_charges = self._static_charges_numpy(params, symbols, fixed_charges)
        e_pauli = float(np.dot(pauli_coeff, pauli_sample))
        e_london = float(np.dot(london_coeff, london_sample))
        e_coul = 0.0
        charges = np.asarray(coulomb_charges, dtype=float)
        e_qeq = 0.0
        e_image = 0.0
        if self.toggles.use_density_charge:
            e_coul = float(np.dot(coulomb_charges, coulomb_sample))
        if self.toggles.use_ct_qeq:
            qeq = solve_qeq_with_reservoir(
                positions=positions,
                chi=atom_params["chi"],
                hardness=atom_params["hardness"],
                phi_external=coulomb_sample,
                image_plane=params.image_plane,
                image_damping=params.image_damping,
                image_scale=params.image_scale,
                reservoir_chi=params.reservoir_chi,
                reservoir_hardness=params.reservoir_hardness,
                total_charge=params.total_charge,
                include_image=self.toggles.use_image_charge,
            )
            charges = qeq["charges"]
            e_qeq = float(qeq["energy"])
            e_image = float(qeq["image_energy"])

        if self.toggles.use_image_charge_fixed and not self.toggles.use_ct_qeq:
            dz = positions_q[:, 2] - float(params.image_plane)
            dz_safe = np.where(dz > float(params.image_damping), dz, float(params.image_damping))
            e_image = -float(params.image_scale) * float(np.sum(coulomb_charges ** 2 / (4.0 * dz_safe)))

        e_reactive = 0.0
        if self.toggles.use_reactive_grid:
            e_reactive = float(np.dot(atom_params["reactive"], reactive_selected))

        total = e_pauli + e_london + e_coul + e_qeq + e_reactive + e_image
        return total, {
            "pauli": e_pauli,
            "london": e_london,
            "coulomb": e_coul,
            "qeq": e_qeq,
            "image": e_image,
            "reactive": e_reactive,
            "charges": charges,
        }

    def _force_single(self, positions, symbols, params: HybridParameters, fixed_charges=None):
        energy_fn = lambda pos: self.energy_single(pos, symbols, params, fixed_charges=fixed_charges)[0]
        return finite_difference_force(energy_fn, positions)

    if HAS_JAX:

        def _molecule_plq_jax(self, species_idx, params_tree):
            if not self.hybrid_config.use_req_plq:
                return params_tree["pauli"][species_idx], params_tree["london"][species_idx]
            radius = self._jax_base_req_radius[species_idx] + params_tree["req_radius_offset"][species_idx]
            sqrt_energy = self._jax_base_req_sqrt_energy[species_idx] * params_tree["req_energy_scale"][species_idx]
            expo = jnp.exp(jnp.asarray(self.grid_config.alpha_morse, dtype=jnp.float64) * radius)
            london = expo * sqrt_energy
            pauli = expo * london
            return pauli, london

        def _sample_field_jax(self, values, positions):
            if self.grid_config.interpolation_order >= 3:
                return _sample_bspline_jax(values, positions, self._jax_origin, self._jax_cell_inv)
            return _sample_trilinear_jax(values, positions, self._jax_origin, self._jax_cell_inv)

        def _energy_single_jax(self, positions, species_idx, fixed_charges, charge_mode, params_tree):
            positions_pl = positions.at[:, 2].add(params_tree["sample_shift_z"])
            positions_q = positions_pl.at[:, 2].add(params_tree["coulomb_shift_z"])
            pauli_sample = self._sample_field_jax(
                self._jax_sample_grids["pauli_zyx"], positions_pl
            )
            london_sample = self._sample_field_jax(
                self._jax_sample_grids["london_zyx"], positions_pl
            )
            coulomb_sample = self._sample_field_jax(
                self._jax_sample_grids["coulomb_zyx"], positions_q
            )
            reactive_sample = self._sample_field_jax(
                self._jax_sample_grids["reactive_zyxc"], positions_pl
            )
            if reactive_sample.ndim == 1:
                reactive_selected = reactive_sample
            else:
                reactive_selected = reactive_sample[jnp.arange(species_idx.shape[0]), species_idx]

            pauli_coeff, london_coeff = self._molecule_plq_jax(species_idx, params_tree)
            reactive_coeff = params_tree["reactive"][species_idx]
            static_charge = params_tree["static_charge"][species_idx]
            chi = params_tree["chi"][species_idx]
            hardness = params_tree["hardness"][species_idx]
            coulomb_charges = (fixed_charges * (1.0 + static_charge) * charge_mode) + (static_charge * (1.0 - charge_mode))

            e_pauli = jnp.dot(pauli_coeff, pauli_sample)
            e_london = jnp.dot(london_coeff, london_sample)
            e_coul = jnp.array(0.0, dtype=jnp.float64)
            charges = coulomb_charges
            e_qeq = jnp.array(0.0, dtype=jnp.float64)
            e_image = jnp.array(0.0, dtype=jnp.float64)
            if self.toggles.use_density_charge:
                e_coul = jnp.dot(coulomb_charges, coulomb_sample)
            if self.toggles.use_ct_qeq:
                qeq = solve_qeq_with_reservoir_jax(
                    positions=positions,
                    chi=chi,
                    hardness=hardness,
                    phi_external=coulomb_sample,
                    image_plane=params_tree["image_plane"],
                    image_damping=params_tree["image_damping"],
                    image_scale=params_tree["image_scale"],
                    reservoir_chi=params_tree["reservoir_chi"],
                    reservoir_hardness=params_tree["reservoir_hardness"],
                    total_charge=params_tree["total_charge"],
                    include_image=self.toggles.use_image_charge,
                )
                charges = qeq["charges"]
                e_qeq = qeq["energy"]
                e_image = qeq["image_energy"]

            if self.toggles.use_image_charge_fixed and not self.toggles.use_ct_qeq:
                dz = positions_q[:, 2] - params_tree["image_plane"]
                dz_safe = jnp.where(dz > params_tree["image_damping"], dz, params_tree["image_damping"])
                e_image = -params_tree["image_scale"] * jnp.sum(coulomb_charges ** 2 / (4.0 * dz_safe))

            e_reactive = jnp.array(0.0, dtype=jnp.float64)
            if self.toggles.use_reactive_grid:
                e_reactive = jnp.dot(reactive_coeff, reactive_selected)

            total = e_pauli + e_london + e_coul + e_qeq + e_reactive + e_image
            return total, {
                "pauli": e_pauli,
                "london": e_london,
                "coulomb": e_coul,
                "qeq": e_qeq,
                "image": e_image,
                "reactive": e_reactive,
                "charges": charges,
            }

        def _get_jax_batch_predictor(self, natoms: int, compute_forces: bool):
            key = (natoms, bool(compute_forces))
            if key in self._jax_predictor_cache:
                return self._jax_predictor_cache[key]

            if compute_forces:

                def per_pose(positions, species_idx, fixed_charges, charge_mode, params_tree):
                    (energy, aux), grad = jax.value_and_grad(self._energy_single_jax, argnums=0, has_aux=True)(
                        positions, species_idx, fixed_charges, charge_mode, params_tree
                    )
                    return energy, -grad, aux

                fn = jax.jit(jax.vmap(per_pose, in_axes=(0, None, None, None, None)))
            else:

                def per_pose(positions, species_idx, fixed_charges, charge_mode, params_tree):
                    return self._energy_single_jax(positions, species_idx, fixed_charges, charge_mode, params_tree)

                fn = jax.jit(jax.vmap(per_pose, in_axes=(0, None, None, None, None)))
            self._jax_predictor_cache[key] = fn
            return fn

        def loss_for_batch_jax(self, positions_batch, species_idx, fixed_charges, charge_mode, params_tree, energies_true, forces_true, force_weight, energy_weight):
            positions_batch = jnp.asarray(positions_batch, dtype=jnp.float64)
            species_idx = jnp.asarray(species_idx, dtype=jnp.int32)
            fixed_charges = jnp.asarray(fixed_charges, dtype=jnp.float64)
            charge_mode = jnp.asarray(charge_mode, dtype=jnp.float64)
            energies_true = jnp.asarray(energies_true, dtype=jnp.float64)
            forces_true = jnp.asarray(forces_true, dtype=jnp.float64)
            if force_weight > 0.0:
                predictor = self._get_jax_batch_predictor(species_idx.shape[0], True)
                energies_pred, forces_pred, _ = predictor(positions_batch, species_idx, fixed_charges, charge_mode, params_tree)
            else:
                predictor = self._get_jax_batch_predictor(species_idx.shape[0], False)
                energies_pred, _ = predictor(positions_batch, species_idx, fixed_charges, charge_mode, params_tree)
                forces_pred = jnp.zeros_like(forces_true)
            energy_loss = jnp.mean((energies_pred - energies_true) ** 2)
            force_loss = jnp.mean((forces_pred - forces_true) ** 2) if force_weight > 0.0 else jnp.array(0.0, dtype=jnp.float64)
            return energy_weight * energy_loss + force_weight * force_loss

    def evaluate_batch(self, poses: PoseBatch, params: HybridParameters | None = None, compute_forces: bool = True):
        params = params or default_hybrid_parameters(self.reactive_elements, poses.metadata["z_image"])
        fixed_charges = poses.adsorbate.charges
        charge_mode = 1.0 if fixed_charges is not None and np.max(np.abs(fixed_charges)) > 0.0 else 0.0
        if self.use_jax:
            params_tree = self.parameter_arrays(params, use_jax=True)
            species_idx = jnp.asarray(self.make_species_indices(poses.adsorbate.symbols), dtype=jnp.int32)
            fixed_charge_array = jnp.asarray(
                np.asarray(fixed_charges, dtype=float) if fixed_charges is not None else np.zeros(poses.adsorbate.natoms, dtype=float),
                dtype=jnp.float64,
            )
            charge_mode_array = jnp.asarray(charge_mode, dtype=jnp.float64)
            positions_batch = jnp.asarray(np.asarray(poses.positions, dtype=float), dtype=jnp.float64)
            if compute_forces:
                predictor = self._get_jax_batch_predictor(poses.adsorbate.natoms, True)
                energies, forces, aux = predictor(positions_batch, species_idx, fixed_charge_array, charge_mode_array, params_tree)
            else:
                predictor = self._get_jax_batch_predictor(poses.adsorbate.natoms, False)
                energies, aux = predictor(positions_batch, species_idx, fixed_charge_array, charge_mode_array, params_tree)
                forces = jnp.zeros_like(positions_batch)
            return ModelEvaluation(
                energies=as_numpy(energies),
                forces=as_numpy(forces),
                components={key: as_numpy(value) for key, value in aux.items() if key != "charges"},
                charges=as_numpy(aux["charges"]),
            )

        energies = []
        forces = []
        components = {name: [] for name in ["pauli", "london", "coulomb", "qeq", "image", "reactive"]}
        charge_list = []
        for positions in poses.positions:
            energy, detail = self.energy_single(positions, poses.adsorbate.symbols, params, fixed_charges=fixed_charges)
            energies.append(energy)
            if compute_forces:
                forces.append(self._force_single(positions, poses.adsorbate.symbols, params, fixed_charges=fixed_charges))
            for key in components:
                components[key].append(detail[key])
            charge_list.append(detail["charges"])
        return ModelEvaluation(
            energies=np.asarray(energies, dtype=float),
            forces=np.asarray(forces, dtype=float) if compute_forces else np.zeros_like(poses.positions),
            components={key: np.asarray(value, dtype=float) for key, value in components.items()},
            charges=np.asarray(charge_list, dtype=float),
        )
