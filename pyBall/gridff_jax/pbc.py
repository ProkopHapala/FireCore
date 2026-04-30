"""Periodic Boundary Condition (PBC) utilities for variable-size DFT slabs.

When the DFT slab has a specific supercell size (e.g., 3×3, 4×4, 12×12),
the PLQ grid is built from that cell's CHGCAR/LOCPOT. Molecules placed near
cell boundaries may have atoms extending beyond the cell. This module provides:

1. **Position wrapping**: Wrap Cartesian atom positions into the cell via PBC
2. **Cell-size validation**: Check that the DFT cell is large enough for the molecule
3. **Supercell detection**: Detect how many primitive cells fit in the DFT supercell
4. **Tiling computation**: Minimum teacher slab tiling for molecule+PBC correctness
5. **PBCInfo metadata**: Propagate PBC context through sampler → teacher → student

The grid interpolation functions in ``hybrid_energy/model.py`` already wrap
fractional (u,v) coordinates via ``frac[:, 0:2] - floor(frac[:, 0:2])`` and
use modulo indexing ``% nx``, ``% ny``.  This module ensures the *sampler*
and *diagnostics* also handle PBC consistently.

Usage
-----
    from pyBall.gridff_jax.pbc import (
        PBCInfo, wrap_positions_into_cell, check_molecule_cell_compatibility,
        detect_supercell_multiplicity, compute_minimum_teacher_tiling,
    )

    pbc_info = PBCInfo.from_density(density)
    wrapped = wrap_positions_into_cell(positions, density.cell, density.origin)
    ok, msg = check_molecule_cell_compatibility(adsorbate, density)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class PBCInfo:
    """PBC metadata container propagated through the pipeline.

    Attributes
    ----------
    cell : (3, 3) array
        Unit cell vectors of the DFT slab.
    cell_lengths : (3,) array
        Lengths of the three cell vectors (Å).
    cell_angles_deg : (3,) array
        Angles between cell vectors (degrees): [alpha, beta, gamma].
    lateral_area : float
        Area of the lateral (a × b) cell face (Å²).
    supercell_multiplicity : tuple[int, int]
        Estimated (nx, ny) primitive cells in the DFT supercell.
        (1, 1) if no supercell detected or detection not applicable.
    primitive_cell : (3, 3) array or None
        Estimated primitive cell vectors, if supercell detected.
    min_teacher_tile : tuple[int, int]
        Minimum (nx, ny) tiling the teacher needs for the largest molecule.
    warnings : list[str]
        Any PBC-related warnings (e.g., cell too small).
    """
    cell: np.ndarray = field(default_factory=lambda: np.eye(3))
    cell_lengths: np.ndarray = field(default_factory=lambda: np.ones(3))
    cell_angles_deg: np.ndarray = field(default_factory=lambda: np.array([90.0, 90.0, 90.0]))
    lateral_area: float = 0.0
    supercell_multiplicity: tuple[int, int] = (1, 1)
    primitive_cell: np.ndarray | None = None
    min_teacher_tile: tuple[int, int] = (1, 1)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_density(cls, density, adsorbate=None) -> PBCInfo:
        """Build PBCInfo from a DensityData object and optional adsorbate.

        Parameters
        ----------
        density : DensityData
            Substrate data with cell, positions, symbols.
        adsorbate : AdsorbateDefinition, optional
            If provided, compute minimum tiling and compatibility check.
        """
        cell = np.asarray(density.cell, dtype=float)
        lengths = np.array([np.linalg.norm(cell[i]) for i in range(3)])
        angles = _cell_angles(cell)
        area = float(np.linalg.norm(np.cross(cell[0], cell[1])))

        nx, ny = detect_supercell_multiplicity(density)
        prim_cell = None
        if nx > 1 or ny > 1:
            prim_cell = cell.copy()
            prim_cell[0] = cell[0] / nx
            prim_cell[1] = cell[1] / ny

        info = cls(
            cell=cell,
            cell_lengths=lengths,
            cell_angles_deg=angles,
            lateral_area=area,
            supercell_multiplicity=(nx, ny),
            primitive_cell=prim_cell,
        )

        if adsorbate is not None:
            tile = compute_minimum_teacher_tiling(adsorbate, cell)
            info.min_teacher_tile = tile
            ok, msg = check_molecule_cell_compatibility(adsorbate, density)
            if not ok:
                info.warnings.append(msg)

        info.metadata = {
            "cell_lengths_A": lengths.tolist(),
            "cell_angles_deg": angles.tolist(),
            "lateral_area_A2": area,
            "supercell_nx": nx,
            "supercell_ny": ny,
        }
        return info


# ---------------------------------------------------------------------------
#  Position wrapping
# ---------------------------------------------------------------------------

def wrap_positions_into_cell(
    positions: np.ndarray,
    cell: np.ndarray,
    origin: np.ndarray | None = None,
    wrap_xy_only: bool = True,
) -> np.ndarray:
    """Wrap Cartesian atom positions into the periodic cell.

    For molecule-on-surface calculations, only the lateral (x, y) directions
    are periodic.  The z-direction is NOT wrapped by default because the
    molecule sits above the slab and z is not periodic in a slab calculation.

    Parameters
    ----------
    positions : (N, 3) or (n_poses, N, 3) array
        Cartesian coordinates.
    cell : (3, 3) array
        Unit cell vectors.
    origin : (3,) array, optional
        Cell origin. Default: [0, 0, 0].
    wrap_xy_only : bool
        If True (default), only wrap x and y; leave z unchanged.

    Returns
    -------
    wrapped : same shape as positions
        Positions wrapped into [origin, origin+cell).
    """
    positions = np.asarray(positions, dtype=float).copy()
    cell = np.asarray(cell, dtype=float)
    if origin is None:
        origin = np.zeros(3, dtype=float)
    else:
        origin = np.asarray(origin, dtype=float)

    original_shape = positions.shape
    if positions.ndim == 2:
        positions = positions[np.newaxis, :, :]

    cell_inv = np.linalg.inv(cell)

    for i in range(positions.shape[0]):
        shifted = positions[i] - origin[None, :]
        frac = shifted @ cell_inv.T
        if wrap_xy_only:
            frac[:, 0] = frac[:, 0] - np.floor(frac[:, 0])
            frac[:, 1] = frac[:, 1] - np.floor(frac[:, 1])
        else:
            frac = frac - np.floor(frac)
        positions[i] = frac @ cell.T + origin[None, :]

    return positions.reshape(original_shape)


def wrap_anchor_into_cell(
    positions: np.ndarray,
    anchor_index: int,
    cell: np.ndarray,
    origin: np.ndarray | None = None,
) -> np.ndarray:
    """Wrap molecule so that the anchor atom is inside the cell.

    Unlike ``wrap_positions_into_cell`` (which wraps each atom independently),
    this function translates the ENTIRE molecule by a lattice vector so the
    anchor atom lies inside the cell.  This preserves the molecule's internal
    geometry while ensuring the anchor is in-cell for grid interpolation.

    Parameters
    ----------
    positions : (N, 3) array
        Cartesian coordinates of one pose's atoms.
    anchor_index : int
        Index of the anchor atom.
    cell : (3, 3) array
        Unit cell vectors.
    origin : (3,) array, optional
        Cell origin.

    Returns
    -------
    wrapped : (N, 3) array
    """
    positions = np.asarray(positions, dtype=float).copy()
    cell = np.asarray(cell, dtype=float)
    if origin is None:
        origin = np.zeros(3, dtype=float)
    else:
        origin = np.asarray(origin, dtype=float)

    cell_inv = np.linalg.inv(cell)
    anchor_shifted = positions[anchor_index] - origin
    anchor_frac = anchor_shifted @ cell_inv.T
    shift_frac = np.floor(anchor_frac[:2])
    shift_cart = np.zeros(3, dtype=float)
    shift_cart += shift_frac[0] * cell[0]
    shift_cart += shift_frac[1] * cell[1]
    positions -= shift_cart[None, :]
    return positions


# ---------------------------------------------------------------------------
#  Cell-size validation
# ---------------------------------------------------------------------------

def molecule_extent_in_cell(
    adsorbate,
    cell: np.ndarray,
) -> dict[str, float]:
    """Compute the molecule extent along each lateral cell axis.

    Returns
    -------
    dict with keys: 'extent_a', 'extent_b', 'cell_a', 'cell_b',
        'ratio_a', 'ratio_b'.
    """
    cell = np.asarray(cell, dtype=float)
    pos = np.asarray(adsorbate.positions, dtype=float)
    pos_centered = pos - pos[adsorbate.anchor_index]

    cell_a = float(np.linalg.norm(cell[0]))
    cell_b = float(np.linalg.norm(cell[1]))

    extent_a = float(pos_centered[:, 0].max() - pos_centered[:, 0].min()) if len(pos) > 1 else 0.0
    extent_b = float(pos_centered[:, 1].max() - pos_centered[:, 1].min()) if len(pos) > 1 else 0.0

    # Also compute extent along cell vector directions (for non-orthogonal cells)
    a_hat = cell[0] / max(cell_a, 1e-12)
    b_hat = cell[1] / max(cell_b, 1e-12)
    proj_a = pos_centered @ a_hat
    proj_b = pos_centered @ b_hat
    extent_a_proj = float(proj_a.max() - proj_a.min()) if len(pos) > 1 else 0.0
    extent_b_proj = float(proj_b.max() - proj_b.min()) if len(pos) > 1 else 0.0

    # Use the larger of Cartesian and projected extents
    extent_a = max(extent_a, extent_a_proj)
    extent_b = max(extent_b, extent_b_proj)

    return {
        "extent_a": extent_a,
        "extent_b": extent_b,
        "cell_a": cell_a,
        "cell_b": cell_b,
        "ratio_a": extent_a / max(cell_a, 1e-12),
        "ratio_b": extent_b / max(cell_b, 1e-12),
    }


def check_molecule_cell_compatibility(
    adsorbate,
    density,
    max_ratio: float = 0.45,
) -> tuple[bool, str]:
    """Check if the DFT cell is large enough for the molecule.

    A molecule is "compatible" if its lateral extent is < max_ratio × cell length
    in both a and b directions.  If the extent is too large, periodic images
    in the PLQ grid will overlap, producing unphysical energies.

    Parameters
    ----------
    adsorbate : AdsorbateDefinition
    density : DensityData
    max_ratio : float
        Maximum allowed molecule-extent / cell-length ratio.

    Returns
    -------
    (ok, message) : tuple[bool, str]
        ok is True if compatible, False otherwise. message explains.
    """
    info = molecule_extent_in_cell(adsorbate, density.cell)

    issues = []
    if info["ratio_a"] > max_ratio:
        issues.append(
            f"Molecule '{adsorbate.name}' extent along a ({info['extent_a']:.2f} Å) "
            f"is {info['ratio_a']:.0%} of cell a ({info['cell_a']:.2f} Å). "
            f"Max allowed: {max_ratio:.0%}. Consider using a larger DFT supercell."
        )
    if info["ratio_b"] > max_ratio:
        issues.append(
            f"Molecule '{adsorbate.name}' extent along b ({info['extent_b']:.2f} Å) "
            f"is {info['ratio_b']:.0%} of cell b ({info['cell_b']:.2f} Å). "
            f"Max allowed: {max_ratio:.0%}. Consider using a larger DFT supercell."
        )

    if issues:
        return False, " | ".join(issues)
    return True, f"OK: molecule fits within cell (a: {info['ratio_a']:.1%}, b: {info['ratio_b']:.1%})"


# ---------------------------------------------------------------------------
#  Supercell detection
# ---------------------------------------------------------------------------

def detect_supercell_multiplicity(
    density,
    nn_tol: float = 0.3,
) -> tuple[int, int]:
    """Detect how many primitive cells fit in the DFT supercell.

    Uses nearest-neighbor distances in the surface layer to estimate the
    primitive cell lattice parameter, then divides the DFT cell lengths by it.

    Parameters
    ----------
    density : DensityData
    nn_tol : float
        Tolerance factor for nearest-neighbor matching.

    Returns
    -------
    (nx, ny) : tuple[int, int]
        Number of primitive cells along a and b.  (1, 1) if detection fails.
    """
    cell = np.asarray(density.cell, dtype=float)
    positions = np.asarray(density.positions, dtype=float)
    symbols = list(density.symbols)

    # Work with the topmost layer (surface atoms)
    zs = positions[:, 2]
    z_top = zs.max()
    layer_mask = np.abs(zs - z_top) < 0.5
    surface_pos = positions[layer_mask]

    if len(surface_pos) < 2:
        return (1, 1)

    # Find nearest-neighbor distance in the surface layer
    nn_dist = _nearest_neighbor_distance(surface_pos, cell)
    if nn_dist < 0.5:
        return (1, 1)

    cell_a = float(np.linalg.norm(cell[0]))
    cell_b = float(np.linalg.norm(cell[1]))

    # Number of atoms in the surface layer gives us the multiplicity
    n_surface = len(surface_pos)

    # For a simple metal surface (one atom per primitive cell in the layer),
    # n_surface = nx * ny. For fcc(111), each layer has one atom per primitive cell.
    # Use nearest-neighbor distance to estimate primitive cell parameter.
    nx = max(1, round(cell_a / nn_dist))
    ny = max(1, round(cell_b / nn_dist))

    # Validate: n_surface should be approximately nx * ny (for one-atom-per-cell surfaces)
    # or nx * ny * basis_size for multi-atom surfaces
    if n_surface > 0:
        basis_size = round(n_surface / max(nx * ny, 1))
        if basis_size < 1:
            basis_size = 1
        expected = nx * ny * basis_size
        if abs(expected - n_surface) > nn_tol * n_surface:
            # Detection unreliable, fall back
            return (1, 1)

    return (nx, ny)


def _nearest_neighbor_distance(positions: np.ndarray, cell: np.ndarray) -> float:
    """Find the minimum nearest-neighbor distance including PBC images."""
    cell = np.asarray(cell, dtype=float)
    n = len(positions)
    if n < 2:
        return 0.0

    min_dist = float("inf")
    for ia in range(-1, 2):
        for ib in range(-1, 2):
            shift = ia * cell[0] + ib * cell[1]
            shifted = positions + shift[None, :]
            for i in range(n):
                for j in range(n):
                    if ia == 0 and ib == 0 and i == j:
                        continue
                    d = float(np.linalg.norm(positions[i] - shifted[j]))
                    if d < min_dist and d > 0.5:
                        min_dist = d
    return min_dist


# ---------------------------------------------------------------------------
#  Teacher tiling computation
# ---------------------------------------------------------------------------

def compute_minimum_teacher_tiling(
    adsorbate,
    cell: np.ndarray,
    safety_margin: float = 1.0,
) -> tuple[int, int]:
    """Compute the minimum slab tiling for teacher evaluation.

    The teacher (MAD-SURF/MACE) evaluates the molecule + slab as an ASE Atoms
    object with PBC.  If the slab cell is small relative to the molecule, the
    molecule interacts with its own periodic image.  This function computes the
    minimum (nx, ny) tiling to avoid self-interaction.

    Parameters
    ----------
    adsorbate : AdsorbateDefinition
    cell : (3, 3) array
        DFT slab cell vectors.
    safety_margin : float
        Extra margin in Angstrom beyond the molecule extent.

    Returns
    -------
    (nx, ny) : tuple[int, int]
    """
    cell = np.asarray(cell, dtype=float)
    pos = np.asarray(adsorbate.positions, dtype=float)
    if len(pos) <= 1:
        return (1, 1)

    pos_centered = pos - pos.mean(axis=0, keepdims=True)
    tiles = []
    for axis in range(2):
        cell_length = float(np.linalg.norm(cell[axis]))
        mol_extent = float(pos_centered[:, axis].max() - pos_centered[:, axis].min())
        n = math.ceil((mol_extent + safety_margin) / max(cell_length, 1e-12))
        tiles.append(max(n, 1))
    return (tiles[0], tiles[1])


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _cell_angles(cell: np.ndarray) -> np.ndarray:
    """Compute angles (alpha, beta, gamma) between cell vectors in degrees."""
    cell = np.asarray(cell, dtype=float)
    a, b, c = cell[0], cell[1], cell[2]
    la, lb, lc = [max(np.linalg.norm(v), 1e-12) for v in (a, b, c)]
    alpha = math.degrees(math.acos(np.clip(np.dot(b, c) / (lb * lc), -1, 1)))
    beta = math.degrees(math.acos(np.clip(np.dot(a, c) / (la * lc), -1, 1)))
    gamma = math.degrees(math.acos(np.clip(np.dot(a, b) / (la * lb), -1, 1)))
    return np.array([alpha, beta, gamma], dtype=float)


def fractional_lateral_to_primitive(
    uv_frac: np.ndarray,
    supercell_nx: int,
    supercell_ny: int,
) -> np.ndarray:
    """Map fractional (u,v) in the DFT supercell to the primitive cell index.

    Parameters
    ----------
    uv_frac : (N, 2) array
        Fractional coordinates in the DFT supercell [0, 1).
    supercell_nx, supercell_ny : int
        Number of primitive cells along a and b.

    Returns
    -------
    (N, 2) array
        Fractional coordinates within the first primitive cell [0, 1/nx) × [0, 1/ny).
    """
    uv = np.asarray(uv_frac, dtype=float).copy()
    uv[:, 0] = (uv[:, 0] * supercell_nx) % 1.0 / supercell_nx
    uv[:, 1] = (uv[:, 1] * supercell_ny) % 1.0 / supercell_ny
    return uv


def primitive_uv_grid(
    n_u: int,
    n_v: int,
    supercell_nx: int = 1,
    supercell_ny: int = 1,
) -> np.ndarray:
    """Build a (u,v) grid sampling only the first primitive cell.

    Instead of sampling u ∈ [0, 1) over the full DFT supercell, sample
    u ∈ [0, 1/nx) where nx is the supercell multiplicity.  This avoids
    redundant sampling over translationally equivalent positions.

    Parameters
    ----------
    n_u, n_v : int
        Number of grid points per direction.
    supercell_nx, supercell_ny : int
        Primitive cell counts.

    Returns
    -------
    (n_u * n_v, 2) array
        Fractional coordinates in the DFT supercell frame.
    """
    u_range = 1.0 / max(supercell_nx, 1)
    v_range = 1.0 / max(supercell_ny, 1)
    u_vals = np.linspace(0.0, u_range, n_u, endpoint=False)
    v_vals = np.linspace(0.0, v_range, n_v, endpoint=False)
    uu, vv = np.meshgrid(u_vals, v_vals, indexing="ij")
    return np.column_stack([uu.ravel(), vv.ravel()])
