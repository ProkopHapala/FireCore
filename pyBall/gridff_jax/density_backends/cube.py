"""Gaussian-cube density backend.

Reads volumetric density (and optionally a second cube for electrostatic
potential) in the Gaussian-cube format produced by Psi4, Gaussian, NWChem,
ORCA, Q-Chem and most other quantum-chemistry codes.

Cube format (units: atomic units, Bohr, e/Bohr^3)
-------------------------------------------------
  line 1: title (free text)
  line 2: title (free text)
  line 3: natoms  origin_x origin_y origin_z
  line 4: nx  dx_x dx_y dx_z
  line 5: ny  dy_x dy_y dy_z
  line 6: nz  dz_x dz_y dz_z
  next natoms lines: Z  charge  x y z   (atomic positions; charge often 0)
  data: nx*ny*nz floats, x-slowest order (i.e., for ix in range(nx): for iy
        in range(ny): for iz in range(nz): write value(ix,iy,iz))

If natoms is negative on line 3, an additional "Dset Ids" record appears
before the data (this backend treats |natoms| atoms and ignores the Dset
Ids record — sufficient for single-volume densities).

This backend converts:
  - Lengths Bohr -> Angstrom (factor 0.529177210903)
  - Densities e/Bohr^3 -> e/Angstrom^3 (factor 1/Bohr^3)

Configuration
-------------
Reuses ``DensityBackendConfig`` fields with this mapping:

  config.kind          = "cube"
  config.chgcar_path   = path to density cube   (required)
  config.locpot_path   = path to potential cube (optional)
  config.grid_shape    = optional resample (Nx, Ny, Nz)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..interfaces import DensityBackend, DensityData
from ..utils import voxel_spacing_from_cell

BOHR_TO_ANG = 0.529177210903
BOHR3_TO_ANG3 = BOHR_TO_ANG ** 3


@dataclass
class ParsedCube:
    cell: np.ndarray            # (3, 3) Å
    origin: np.ndarray          # (3,) Å
    symbols: list[str]
    positions: np.ndarray       # (natoms, 3) Å
    values_xyz: np.ndarray      # (nx, ny, nz) — natural cube ordering, raw cube units
    voxel_volume_bohr3: float   # for unit-aware conversion downstream


def _z_to_symbol(z: int) -> str:
    """Map atomic number → symbol via ASE's full periodic table."""
    from ase.data import chemical_symbols
    z = int(z)
    if z < 1 or z >= len(chemical_symbols):
        raise ValueError(f"Atomic number {z} is outside ase.data.chemical_symbols range")
    return chemical_symbols[z]


def read_cube(path: str | Path) -> ParsedCube:
    """Parse a Gaussian cube file. Returns a ParsedCube with Å + e/Å^3 units."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        lines = fh.readlines()
    if len(lines) < 6:
        raise ValueError(f"Cube file {path} has fewer than 6 header lines")

    parts3 = lines[2].split()
    natoms_signed = int(parts3[0])
    natoms = abs(natoms_signed)
    origin_bohr = np.array([float(x) for x in parts3[1:4]], dtype=float)

    def _voxel_line(line: str) -> tuple[int, np.ndarray]:
        ps = line.split()
        n = int(ps[0])
        v = np.array([float(x) for x in ps[1:4]], dtype=float)
        return n, v

    nx, vx = _voxel_line(lines[3])
    ny, vy = _voxel_line(lines[4])
    nz, vz = _voxel_line(lines[5])

    # Voxel vectors are in Bohr if any n>0, in Angstrom if all n<0 (rare).
    # We treat the common case: n>0 → units are Bohr.
    if nx < 0 or ny < 0 or nz < 0:
        raise NotImplementedError(
            "Cube files with negative voxel counts (Å-unit indicator) are not "
            "yet supported by this backend; please convert to standard Bohr units."
        )

    symbols: list[str] = []
    positions_bohr = np.zeros((natoms, 3), dtype=float)
    for i in range(natoms):
        ps = lines[6 + i].split()
        z = int(ps[0])
        # ps[1] is "charge" — usually 0.0 in cubes from densities
        positions_bohr[i] = [float(ps[2]), float(ps[3]), float(ps[4])]
        symbols.append(_z_to_symbol(z))

    # Skip the Dset Ids record if natoms_signed < 0
    data_start = 6 + natoms
    if natoms_signed < 0:
        data_start += 1

    # Read the data block via np.fromiter on whitespace-split tokens.
    # (Previously used np.fromstring, which is deprecated since NumPy 1.14
    # and emits a DeprecationWarning that will become a hard error.)
    import re as _re
    _tokens = (t for line in lines[data_start:] for t in _re.split(r"\s+", line.strip()) if t)
    flat = np.fromiter((float(t) for t in _tokens), dtype=float)
    expected = nx * ny * nz
    if flat.size != expected:
        raise ValueError(
            f"Cube file {path} data block has {flat.size} entries, expected "
            f"{expected} (= {nx}*{ny}*{nz})"
        )
    # Cube ordering: x slowest, z fastest → reshape (nx, ny, nz)
    values_xyz = flat.reshape(nx, ny, nz)

    # Convert lengths Bohr → Å. The cube DATA block is returned RAW (in cube's
    # own units, typically Bohr^-3 for densities or Hartree for potentials) —
    # the caller must apply the unit conversion appropriate for the quantity
    # it is reading. We also return the voxel volume in Bohr^3 to make the
    # conversion mechanical.
    origin = origin_bohr * BOHR_TO_ANG
    vx_ang = vx * BOHR_TO_ANG
    vy_ang = vy * BOHR_TO_ANG
    vz_ang = vz * BOHR_TO_ANG
    cell = np.stack([vx_ang * nx, vy_ang * ny, vz_ang * nz], axis=0)
    positions = positions_bohr * BOHR_TO_ANG
    voxel_volume_bohr3 = float(np.abs(np.linalg.det(np.stack([vx, vy, vz], axis=0))))

    return ParsedCube(
        cell=cell, origin=origin, symbols=symbols,
        positions=positions,
        values_xyz=values_xyz.copy(),
        voxel_volume_bohr3=voxel_volume_bohr3,
    )


def _resample_zyx(values_zyx: np.ndarray, target_xyz) -> np.ndarray:
    """Nearest-neighbour resample (NOT trilinear).

    WARNING: this aliases the field. For force-field fitting on potential
    grids, aliasing creates artificial corrugation. Do not use this to
    coarsen a fine input — supply the cube at the resolution you want.
    Kept for compatibility with the VASP backend, which has the same
    limitation."""
    nx, ny, nz = (int(v) for v in target_xyz)
    src_nz, src_ny, src_nx = values_zyx.shape
    ix = (np.linspace(0, src_nx - 1, nx)).astype(int)
    iy = (np.linspace(0, src_ny - 1, ny)).astype(int)
    iz = (np.linspace(0, src_nz - 1, nz)).astype(int)
    return values_zyx[np.ix_(iz, iy, ix)]


class CubeDensityBackend(DensityBackend):
    """Read substrate density (and optional potential) from cube file(s).

    Configuration fields used:
      - ``config.chgcar_path`` -> path to the density cube (required)
      - ``config.locpot_path`` -> path to a second cube for V(r) (optional)
      - ``config.grid_shape``  -> resample to (Nx, Ny, Nz) if set
    """

    def __init__(self, config, surface=None, grid=None):
        self.config = config
        self.surface = surface
        self.grid = grid

    def load(self) -> DensityData:
        if not self.config.chgcar_path:
            raise ValueError("cube backend requires chgcar_path (density cube file)")
        rho_cube = read_cube(self.config.chgcar_path)
        # values_xyz (nx, ny, nz) → convention used elsewhere is zyx.
        # DENSITY cube convention: values are e/Bohr^3 → divide by Bohr^3
        # to convert to e/Å^3.
        rho_zyx = rho_cube.values_xyz.transpose(2, 1, 0).copy() / BOHR3_TO_ANG3

        v_loc_zyx = None
        if self.config.locpot_path:
            pot_cube = read_cube(self.config.locpot_path)
            if pot_cube.values_xyz.shape != rho_cube.values_xyz.shape:
                raise ValueError(
                    f"Potential cube grid {pot_cube.values_xyz.shape} does not "
                    f"match density cube grid {rho_cube.values_xyz.shape}"
                )
            # POTENTIAL cube convention: values are in Hartree (energy per
            # electron). The framework consumes potentials in eV, so apply
            # ONLY the energy-unit conversion — do NOT divide by Bohr^3.
            HARTREE_TO_EV = 27.211386245988
            v_loc_zyx = pot_cube.values_xyz.transpose(2, 1, 0).copy() * HARTREE_TO_EV

        grid_shape_xyz = tuple(int(v) for v in rho_cube.values_xyz.shape)
        if self.config.grid_shape is not None:
            grid_shape_xyz = tuple(int(v) for v in self.config.grid_shape)
            rho_zyx = _resample_zyx(rho_zyx, grid_shape_xyz)
            if v_loc_zyx is not None:
                v_loc_zyx = _resample_zyx(v_loc_zyx, grid_shape_xyz)

        voxel = voxel_spacing_from_cell(rho_cube.cell, grid_shape_xyz)
        metadata = {
            "backend": "cube",
            "source_files": {
                "density": str(self.config.chgcar_path),
                "potential": str(self.config.locpot_path) if self.config.locpot_path else None,
            },
            "grid_shape_xyz": grid_shape_xyz,
            "grid_shape_zyx": tuple(reversed(grid_shape_xyz)),
            "grid_order": "zyx",
            "n_atoms": len(rho_cube.symbols),
            "symbol_counts": {s: rho_cube.symbols.count(s) for s in sorted(set(rho_cube.symbols))},
            "units_note": "lengths Å, density e/Å^3 (converted from Bohr/e per Bohr^3)",
        }
        return DensityData(
            cell=rho_cube.cell,
            origin=rho_cube.origin,
            voxel=voxel,
            symbols=rho_cube.symbols,
            positions=rho_cube.positions,
            charges=None,
            rho_zyx=rho_zyx,
            v_loc_zyx=v_loc_zyx,
            grid_shape_xyz_hint=(int(grid_shape_xyz[0]), int(grid_shape_xyz[1]), int(grid_shape_xyz[2])),  # type: ignore[arg-type]
            metadata=metadata,
        )
