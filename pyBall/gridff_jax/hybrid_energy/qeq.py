"""Charge-equilibration with metal reservoir and image plane response."""

from __future__ import annotations

import numpy as np

from ..array_api import HAS_JAX

if HAS_JAX:
    import jax.numpy as jnp

COULOMB_CONST = 14.3996448915


def build_image_matrix(positions, image_plane: float, damping: float = 0.5, scale: float = 1.0):
    positions = np.asarray(positions, dtype=float)
    n_atoms = positions.shape[0]
    matrix = np.zeros((n_atoms, n_atoms), dtype=float)
    xy = positions[:, :2]
    z = positions[:, 2]
    for i in range(n_atoms):
        dz = max(z[i] - image_plane, 1.0e-4)
        r_self = np.sqrt((2.0 * dz) ** 2 + damping * damping)
        matrix[i, i] = -scale * COULOMB_CONST / r_self
        for j in range(i + 1, n_atoms):
            dxy = np.linalg.norm(xy[i] - xy[j])
            dz_pair = z[i] + z[j] - 2.0 * image_plane
            rij = np.sqrt(dxy * dxy + dz_pair * dz_pair + damping * damping)
            value = -scale * COULOMB_CONST / rij
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix


def build_coulomb_matrix(positions, damping: float = 0.5):
    positions = np.asarray(positions, dtype=float)
    n_atoms = positions.shape[0]
    matrix = np.zeros((n_atoms, n_atoms), dtype=float)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            distance = np.linalg.norm(positions[i] - positions[j])
            value = COULOMB_CONST / np.sqrt(distance * distance + damping * damping)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix


def solve_qeq_with_reservoir(
    positions,
    chi,
    hardness,
    phi_external,
    image_plane: float,
    image_damping: float = 0.5,
    image_scale: float = 1.0,
    reservoir_chi: float = -4.5,
    reservoir_hardness: float = 7.0,
    total_charge: float = 0.0,
    include_image: bool = True,
):
    positions = np.asarray(positions, dtype=float)
    chi = np.asarray(chi, dtype=float)
    hardness = np.asarray(hardness, dtype=float)
    phi_external = np.asarray(phi_external, dtype=float)
    n_atoms = positions.shape[0]

    interaction = build_coulomb_matrix(positions, damping=image_damping)
    if include_image:
        interaction += build_image_matrix(
            positions,
            image_plane=image_plane,
            damping=image_damping,
            scale=image_scale,
        )
    interaction[np.diag_indices(n_atoms)] += hardness

    matrix = np.zeros((n_atoms + 2, n_atoms + 2), dtype=float)
    rhs = np.zeros(n_atoms + 2, dtype=float)
    matrix[:n_atoms, :n_atoms] = interaction
    matrix[n_atoms, n_atoms] = reservoir_hardness
    matrix[:n_atoms, n_atoms + 1] = 1.0
    matrix[n_atoms, n_atoms + 1] = 1.0
    matrix[n_atoms + 1, :n_atoms] = 1.0
    matrix[n_atoms + 1, n_atoms] = 1.0
    rhs[:n_atoms] = -(chi + phi_external)
    rhs[n_atoms] = -reservoir_chi
    rhs[n_atoms + 1] = total_charge

    solution = np.linalg.solve(matrix, rhs)
    charges = solution[:n_atoms]
    reservoir_charge = solution[n_atoms]
    energy = (
        np.dot(charges, chi + phi_external)
        + 0.5 * np.dot(charges, interaction @ charges)
        + reservoir_chi * reservoir_charge
        + 0.5 * reservoir_hardness * reservoir_charge * reservoir_charge
    )
    image_energy = 0.0
    if include_image:
        image_matrix = build_image_matrix(
            positions,
            image_plane=image_plane,
            damping=image_damping,
            scale=image_scale,
        )
        image_energy = 0.5 * np.dot(charges, image_matrix @ charges)
    return {
        "charges": charges,
        "reservoir_charge": reservoir_charge,
        "energy": energy,
        "image_energy": image_energy,
        "interaction": interaction,
    }


if HAS_JAX:

    def build_image_matrix_jax(positions, image_plane: float, damping: float = 0.5, scale: float = 1.0):
        positions = jnp.asarray(positions, dtype=jnp.float64)
        n_atoms = positions.shape[0]
        eye = jnp.eye(n_atoms, dtype=jnp.float64)
        xy = positions[:, :2]
        z = positions[:, 2]
        dxy2 = jnp.sum((xy[:, None, :] - xy[None, :, :]) ** 2, axis=-1)
        dz_pair = z[:, None] + z[None, :] - 2.0 * image_plane
        off_diag = -scale * COULOMB_CONST / jnp.sqrt(dxy2 + dz_pair * dz_pair + damping * damping)
        dz = jnp.maximum(z - image_plane, 1.0e-4)
        diag = -scale * COULOMB_CONST / jnp.sqrt((2.0 * dz) ** 2 + damping * damping)
        return off_diag * (1.0 - eye) + jnp.diag(diag)


    def build_coulomb_matrix_jax(positions, damping: float = 0.5):
        positions = jnp.asarray(positions, dtype=jnp.float64)
        diff = positions[:, None, :] - positions[None, :, :]
        dist2 = jnp.sum(diff * diff, axis=-1)
        n_atoms = positions.shape[0]
        eye = jnp.eye(n_atoms, dtype=jnp.float64)
        matrix = COULOMB_CONST / jnp.sqrt(dist2 + damping * damping)
        return matrix * (1.0 - eye)


    def solve_qeq_with_reservoir_jax(
        positions,
        chi,
        hardness,
        phi_external,
        image_plane: float,
        image_damping: float = 0.5,
        image_scale: float = 1.0,
        reservoir_chi: float = -4.5,
        reservoir_hardness: float = 7.0,
        total_charge: float = 0.0,
        include_image: bool = True,
    ):
        positions = jnp.asarray(positions, dtype=jnp.float64)
        chi = jnp.asarray(chi, dtype=jnp.float64)
        hardness = jnp.asarray(hardness, dtype=jnp.float64)
        phi_external = jnp.asarray(phi_external, dtype=jnp.float64)
        n_atoms = positions.shape[0]
        eye = jnp.eye(n_atoms, dtype=jnp.float64)

        interaction = build_coulomb_matrix_jax(positions, damping=image_damping)
        image_matrix = jnp.zeros_like(interaction)
        if include_image:
            image_matrix = build_image_matrix_jax(
                positions,
                image_plane=image_plane,
                damping=image_damping,
                scale=image_scale,
            )
            interaction = interaction + image_matrix
        interaction = interaction + jnp.diag(hardness)

        matrix = jnp.zeros((n_atoms + 2, n_atoms + 2), dtype=jnp.float64)
        rhs = jnp.zeros((n_atoms + 2,), dtype=jnp.float64)
        matrix = matrix.at[:n_atoms, :n_atoms].set(interaction)
        matrix = matrix.at[n_atoms, n_atoms].set(reservoir_hardness)
        matrix = matrix.at[:n_atoms, n_atoms + 1].set(1.0)
        matrix = matrix.at[n_atoms, n_atoms + 1].set(1.0)
        matrix = matrix.at[n_atoms + 1, :n_atoms].set(1.0)
        matrix = matrix.at[n_atoms + 1, n_atoms].set(1.0)
        rhs = rhs.at[:n_atoms].set(-(chi + phi_external))
        rhs = rhs.at[n_atoms].set(-reservoir_chi)
        rhs = rhs.at[n_atoms + 1].set(total_charge)

        solution = jnp.linalg.solve(matrix, rhs)
        charges = solution[:n_atoms]
        reservoir_charge = solution[n_atoms]
        energy = (
            jnp.dot(charges, chi + phi_external)
            + 0.5 * jnp.dot(charges, interaction @ charges)
            + reservoir_chi * reservoir_charge
            + 0.5 * reservoir_hardness * reservoir_charge * reservoir_charge
        )
        image_energy = 0.5 * jnp.dot(charges, image_matrix @ charges) if include_image else jnp.array(0.0, dtype=jnp.float64)
        return {
            "charges": charges,
            "reservoir_charge": reservoir_charge,
            "energy": energy,
            "image_energy": image_energy,
            "interaction": interaction,
        }
