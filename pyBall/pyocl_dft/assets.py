"""Asset discovery helpers for :mod:`pyBall.pyocl_dft`.

These routines centralise the logic described in
@/doc/Topics/AFM/AFM.md#115-123 so that both density projection and potential
assembly can fail fast when required data are missing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

BASIS_SUBDIR = Path("Fdata") / "basis"
BASIS_SUFFIXES = (".wf1", ".wf2")


@dataclass(frozen=True)
class BasisInfo:
    root: Path
    elements: Iterable[str]

    def files_for(self, element: str) -> list[Path]:
        return [self.root / f"{element}{suffix}" for suffix in BASIS_SUFFIXES]


def find_basis_dir(start: Optional[Path] = None) -> Path:
    """Locate the Fireball basis directory relative to ``start`` (default: repo)."""

    if start is None:
        start = Path(__file__).resolve().parents[2]
    candidate = start / BASIS_SUBDIR
    if candidate.is_dir():
        return candidate
    raise FileNotFoundError(
        f"Cannot locate basis data under {candidate}. Set FIRECORE_BASIS_DIR "
        "or pass an explicit path to load_wf_basis()."
    )


def verify_basis(root: Path, elements: Iterable[str]) -> BasisInfo:
    """Ensure the required basis files exist for every element symbol."""

    missing: list[str] = []
    for elem in elements:
        for suffix in BASIS_SUFFIXES:
            path = root / f"{elem}{suffix}"
            if not path.is_file():
                missing.append(path.name)
    if missing:
        raise FileNotFoundError(
            "Missing basis files: " + ", ".join(missing) + " in " + str(root)
        )
    return BasisInfo(root=root, elements=tuple(elements))


def ensure_acum_coef(default: tuple[float, float] = (0.0, 2.0)) -> np.ndarray:
    """Return an accumulator coefficient vector as ``float2``."""

    arr = np.asarray(default, dtype=np.float32)
    if arr.shape != (2,):  # Should match usage in projectAtomsDens
        raise ValueError("acumCoef must contain exactly two coefficients")
    return arr
