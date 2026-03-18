"""Built-in rigid benchmark adsorbates for Ag(111) validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pyBall import atomicUtils as au

from ..interfaces import AdsorbateDefinition


def benchmark_adsorbates():
    return {
        "H": AdsorbateDefinition(
            name="H",
            symbols=["H"],
            positions=np.array([[0.0, 0.0, 0.0]], dtype=float),
            charges=np.array([0.0], dtype=float),
            anchor_index=0,
        ),
        "C": AdsorbateDefinition(
            name="C",
            symbols=["C"],
            positions=np.array([[0.0, 0.0, 0.0]], dtype=float),
            charges=np.array([0.0], dtype=float),
            anchor_index=0,
        ),
        "O": AdsorbateDefinition(
            name="O",
            symbols=["O"],
            positions=np.array([[0.0, 0.0, 0.0]], dtype=float),
            charges=np.array([0.0], dtype=float),
            anchor_index=0,
        ),
        "N": AdsorbateDefinition(
            name="N",
            symbols=["N"],
            positions=np.array([[0.0, 0.0, 0.0]], dtype=float),
            charges=np.array([0.0], dtype=float),
            anchor_index=0,
        ),
        "CO": AdsorbateDefinition(
            name="CO",
            symbols=["C", "O"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.128],
                ],
                dtype=float,
            ),
            charges=np.array([-0.0207, 0.0207], dtype=float),
            anchor_index=0,
        ),
        "H2O": AdsorbateDefinition(
            name="H2O",
            symbols=["O", "H", "H"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [0.7586, 0.0, 0.5043],
                    [-0.7586, 0.0, 0.5043],
                ],
                dtype=float,
            ),
            charges=np.array([-0.834, 0.417, 0.417], dtype=float),
            anchor_index=0,
        ),
        "CHONH2": AdsorbateDefinition(
            name="CHONH2",
            symbols=["C", "O", "N", "H", "H", "H"],
            # Formamide (H2N-CHO) planar geometry, B3LYP/6-311G** quality bond lengths:
            # C=O 1.200 Å along +x, C-N 1.352 Å, molecule in xy plane; C at origin (anchor)
            positions=np.array(
                [
                    [ 0.000,  0.000,  0.000],  # C  (anchor)
                    [ 1.200,  0.000,  0.000],  # O  (C=O along +x)
                    [-0.689,  1.134,  0.000],  # N
                    [-0.548, -0.948,  0.000],  # H on C
                    [-1.693,  1.101,  0.000],  # H on N (1st)
                    [-0.220,  2.038,  0.000],  # H on N (2nd)
                ],
                dtype=float,
            ),
            # AMBER GAFF2 RESP partial charges (sum = 0.000 e)
            charges=np.array([0.616, -0.571, -0.862, 0.065, 0.376, 0.376], dtype=float),
            anchor_index=0,
        ),
    }


def load_adsorbate_from_xyz(
    xyz_path: str | Path,
    name: str | None = None,
    anchor_index: int = 0,
    use_input_charges: bool = False,
) -> AdsorbateDefinition:
    atoms = au.AtomicSystem(fname=str(xyz_path))
    if atoms.apos is None or atoms.enames is None:
        raise ValueError(f"Failed to load adsorbate geometry from '{xyz_path}'")
    symbols = [str(symbol) for symbol in atoms.enames]
    positions = np.asarray(atoms.apos, dtype=float)
    charges = None
    if use_input_charges and atoms.qs is not None:
        charges = np.asarray(atoms.qs, dtype=float)
    return AdsorbateDefinition(
        name=name or Path(xyz_path).stem,
        symbols=symbols,
        positions=positions,
        charges=charges,
        anchor_index=int(anchor_index),
        metadata={
            "source_file": str(xyz_path),
            "used_input_charges": bool(use_input_charges),
            "input_charge_sum": float(np.sum(atoms.qs)) if atoms.qs is not None else 0.0,
        },
    )


def get_adsorbate(
    name: str,
    xyz_path: str | Path | None = None,
    anchor_index: int = 0,
    use_input_charges: bool = False,
) -> AdsorbateDefinition:
    if xyz_path:
        return load_adsorbate_from_xyz(
            xyz_path=xyz_path,
            name=name,
            anchor_index=anchor_index,
            use_input_charges=use_input_charges,
        )
    adsorbates = benchmark_adsorbates()
    if name not in adsorbates:
        raise KeyError(f"Unknown built-in adsorbate '{name}'")
    return adsorbates[name]
