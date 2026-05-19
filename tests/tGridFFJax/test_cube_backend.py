#!/usr/bin/env python3
"""Minimal unit test for CubeDensityBackend.

Constructs a tiny synthetic Gaussian-cube file with hand-set values,
parses it through the backend, and verifies that the values are
recovered correctly (modulo Bohr→Å unit conversion and the
x-slowest→zyx layout transpose).

Run:
    JAX_PLATFORMS=cpu python3 tests/tGridFFJax/test_cube_backend.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(ROOT))

from pyBall.gridff_jax.config import DensityBackendConfig, GridConfig
from pyBall.gridff_jax.density_backends import make_density_backend
from pyBall.gridff_jax.density_backends.cube import BOHR_TO_ANG, BOHR3_TO_ANG3


def make_cube_text(nx: int, ny: int, nz: int, values_xyz: np.ndarray,
                   origin_bohr=(0.0, 0.0, 0.0), spacing_bohr=0.5,
                   atoms=()) -> str:
    """Build a Gaussian-cube file as a string. atoms = list of (Z, x, y, z) in Bohr."""
    lines = [
        "test cube",
        "synthetic",
        f"{len(atoms)} {origin_bohr[0]:.6f} {origin_bohr[1]:.6f} {origin_bohr[2]:.6f}",
        f"{nx} {spacing_bohr:.6f} 0.0 0.0",
        f"{ny} 0.0 {spacing_bohr:.6f} 0.0",
        f"{nz} 0.0 0.0 {spacing_bohr:.6f}",
    ]
    for z, x, y, w in atoms:
        lines.append(f"{z} 0.0 {x:.6f} {y:.6f} {w:.6f}")
    # x slowest, z fastest
    chunks = []
    for ix in range(nx):
        for iy in range(ny):
            row = " ".join(f"{values_xyz[ix, iy, iz]:.10e}" for iz in range(nz))
            chunks.append(row)
    lines.extend(chunks)
    return "\n".join(lines) + "\n"


def test_cube_roundtrip():
    nx, ny, nz = 4, 3, 5
    rng = np.random.default_rng(0)
    values_xyz = rng.standard_normal((nx, ny, nz)).astype(float)
    atoms = [(11, 0.5, 0.5, 0.5), (17, 1.0, 1.0, 1.0)]  # Na + Cl

    with tempfile.TemporaryDirectory() as td:
        cube_path = Path(td) / "test.cube"
        cube_path.write_text(make_cube_text(nx, ny, nz, values_xyz,
                                            atoms=atoms, spacing_bohr=0.5))

        cfg = DensityBackendConfig(kind="cube", chgcar_path=str(cube_path))
        density = make_density_backend(cfg, grid=GridConfig()).load()

    # Shape check (zyx order)
    assert density.rho_zyx is not None
    assert density.rho_zyx.shape == (nz, ny, nx), (
        f"shape mismatch: {density.rho_zyx.shape} vs expected ({nz},{ny},{nx})"
    )
    # Value check: cube was in e/Bohr^3, backend converts to e/Å^3
    # values_xyz (nx,ny,nz) -> rho_zyx (nz,ny,nx) via transpose(2,1,0)
    expected_zyx = values_xyz.transpose(2, 1, 0) / BOHR3_TO_ANG3
    np.testing.assert_allclose(density.rho_zyx, expected_zyx, rtol=1e-10, atol=1e-12)

    # Symbols + atomic positions
    assert density.symbols == ["Na", "Cl"]
    expected_positions = np.array([
        [0.5, 0.5, 0.5],
        [1.0, 1.0, 1.0],
    ]) * BOHR_TO_ANG
    np.testing.assert_allclose(density.positions, expected_positions, rtol=1e-12)

    # Cell: spacing * n in Bohr -> Å
    expected_cell_diag = np.array([nx, ny, nz]) * 0.5 * BOHR_TO_ANG
    cell_diag = np.array([density.cell[i, i] for i in range(3)])
    np.testing.assert_allclose(cell_diag, expected_cell_diag, rtol=1e-12)

    print("test_cube_roundtrip: PASSED")


def test_cube_with_potential():
    """Potential cubes must be converted Hartree→eV, NOT divided by Bohr^3."""
    nx, ny, nz = 2, 2, 2
    rng = np.random.default_rng(1)
    rho = rng.standard_normal((nx, ny, nz)).astype(float)
    pot = rng.standard_normal((nx, ny, nz)).astype(float)

    with tempfile.TemporaryDirectory() as td:
        rho_path = Path(td) / "rho.cube"
        pot_path = Path(td) / "pot.cube"
        rho_path.write_text(make_cube_text(nx, ny, nz, rho))
        pot_path.write_text(make_cube_text(nx, ny, nz, pot))

        cfg = DensityBackendConfig(kind="cube",
                                   chgcar_path=str(rho_path),
                                   locpot_path=str(pot_path))
        density = make_density_backend(cfg, grid=GridConfig()).load()

    HARTREE_TO_EV = 27.211386245988
    # Density: e/Bohr^3 → e/Å^3 (divide by Bohr^3)
    expected_rho_zyx = rho.transpose(2, 1, 0) / BOHR3_TO_ANG3
    np.testing.assert_allclose(density.rho_zyx, expected_rho_zyx, rtol=1e-10, atol=1e-12)
    # Potential: Hartree → eV (multiply by 27.211)
    assert density.v_loc_zyx is not None
    expected_pot_zyx = pot.transpose(2, 1, 0) * HARTREE_TO_EV
    np.testing.assert_allclose(density.v_loc_zyx, expected_pot_zyx, rtol=1e-10, atol=1e-12)
    print("test_cube_with_potential: PASSED (density /Bohr^3, potential *Hartree)")


if __name__ == "__main__":
    test_cube_roundtrip()
    test_cube_with_potential()
    print("\nAll tests PASSED")
