"""Atom type registry for transferable per-element-pair parameters.

Phase 1 (current): Per-element parameters (H, C, N, O, ...) that are
    fitted jointly across multiple molecules on a given substrate.
    Parameters are stored as element-substrate pairs: H-Ag, C-Ag, etc.

Phase 2 (future): Per-atom-type parameters (C_R, C_3, H_alkyl, H_hydroxyl)
    using connectivity-based typing similar to UFF.  The data structures
    below are designed to support this extension without refactoring.

The registry loads/saves parameter sets as JSON for portability.

Usage
-----
    from pyBall.gridff_jax.atom_types import (
        ElementPairParams,
        AtomTypeRegistry,
        load_element_pair_params,
        save_element_pair_params,
    )

    registry = AtomTypeRegistry(substrate="Ag")
    registry.set_pair("H", "Ag", radius_offset=0.12, energy_scale=1.5, c6_coeff=23.0)
    save_element_pair_params(registry, "params_Ag.json")
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class ElementPairParams:
    """Parameters for one element-substrate pair (e.g., H-Ag).

    These are the transferable coefficients analogous to UFF parameters,
    fitted from multi-molecule 6D data on a specific substrate.
    """
    element: str = ""
    substrate: str = ""
    radius_offset: float = 0.0
    energy_scale: float = 1.0
    c6_coeff: float = 0.0
    pauli_scale: float = 1.0
    london_scale: float = 1.0
    static_charge: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class AtomTypeParams:
    """Parameters for one atom-type-substrate pair (e.g., C_R-Ag).

    FUTURE: Not used in Phase 1.  The data structure is in place so that
    adding connectivity-based atom typing later does not require refactoring.
    """
    atom_type: str = ""
    base_element: str = ""
    substrate: str = ""
    hybridization: str = ""
    radius_offset: float = 0.0
    energy_scale: float = 1.0
    c6_coeff: float = 0.0
    pauli_scale: float = 1.0
    london_scale: float = 1.0
    static_charge: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
#  UFF-inspired atom type definitions (for future Phase 2)
# ---------------------------------------------------------------------------

UFF_ATOM_TYPES = {
    "H_": {"element": "H", "hybridization": "any", "description": "Generic hydrogen"},
    "C_3": {"element": "C", "hybridization": "sp3", "description": "sp3 carbon"},
    "C_2": {"element": "C", "hybridization": "sp2", "description": "sp2 carbon"},
    "C_R": {"element": "C", "hybridization": "aromatic", "description": "Aromatic carbon"},
    "C_1": {"element": "C", "hybridization": "sp", "description": "sp carbon"},
    "N_3": {"element": "N", "hybridization": "sp3", "description": "sp3 nitrogen"},
    "N_2": {"element": "N", "hybridization": "sp2", "description": "sp2 nitrogen"},
    "N_R": {"element": "N", "hybridization": "aromatic", "description": "Aromatic nitrogen"},
    "N_1": {"element": "N", "hybridization": "sp", "description": "sp nitrogen"},
    "O_3": {"element": "O", "hybridization": "sp3", "description": "sp3 oxygen"},
    "O_2": {"element": "O", "hybridization": "sp2", "description": "sp2 oxygen"},
    "O_R": {"element": "O", "hybridization": "aromatic", "description": "Aromatic oxygen"},
    "O_1": {"element": "O", "hybridization": "sp", "description": "sp oxygen"},
    "F_": {"element": "F", "hybridization": "any", "description": "Fluorine"},
    "S_3+2": {"element": "S", "hybridization": "sp3", "description": "Divalent sulfur"},
    "S_R": {"element": "S", "hybridization": "aromatic", "description": "Aromatic sulfur"},
    "Cl": {"element": "Cl", "hybridization": "any", "description": "Chlorine"},
    "Br": {"element": "Br", "hybridization": "any", "description": "Bromine"},
    "P_3+3": {"element": "P", "hybridization": "sp3", "description": "Trivalent phosphorus"},
    "Si3": {"element": "Si", "hybridization": "sp3", "description": "sp3 silicon"},
}

COINAGE_METALS = {"Cu", "Ag", "Au", "Pt", "Pd", "Ni"}


class AtomTypeRegistry:
    """Registry of element-pair (Phase 1) and atom-type-pair (Phase 2) parameters.

    Parameters
    ----------
    substrate : str
        Substrate element (e.g., "Ag", "Cu", "Au").
    """

    def __init__(self, substrate: str = "Ag"):
        self.substrate = substrate
        self._element_pairs: dict[str, ElementPairParams] = {}
        self._atom_type_pairs: dict[str, AtomTypeParams] = {}

    def pair_key(self, element: str) -> str:
        return f"{element}-{self.substrate}"

    def set_pair(
        self,
        element: str,
        substrate: str | None = None,
        **kwargs: Any,
    ) -> None:
        """Set or update parameters for an element-substrate pair."""
        sub = substrate or self.substrate
        key = f"{element}-{sub}"
        if key in self._element_pairs:
            for k, v in kwargs.items():
                if hasattr(self._element_pairs[key], k):
                    setattr(self._element_pairs[key], k, v)
        else:
            self._element_pairs[key] = ElementPairParams(
                element=element, substrate=sub, **kwargs
            )

    def get_pair(self, element: str, substrate: str | None = None) -> ElementPairParams:
        """Get parameters for an element-substrate pair."""
        sub = substrate or self.substrate
        key = f"{element}-{sub}"
        if key not in self._element_pairs:
            return ElementPairParams(element=element, substrate=sub)
        return self._element_pairs[key]

    def all_pairs(self) -> dict[str, ElementPairParams]:
        return dict(self._element_pairs)

    def elements(self) -> list[str]:
        return sorted(set(p.element for p in self._element_pairs.values()))

    def to_hybrid_params_dict(self) -> dict[str, dict[str, float]]:
        """Convert to the dict format used by HybridParameters."""
        result: dict[str, dict[str, float]] = {
            "req_radius_offset": {},
            "req_energy_scale": {},
            "c6_coeff": {},
            "pauli": {},
            "london": {},
            "static_charge": {},
        }
        for key, params in self._element_pairs.items():
            el = params.element
            result["req_radius_offset"][el] = params.radius_offset
            result["req_energy_scale"][el] = params.energy_scale
            result["c6_coeff"][el] = params.c6_coeff
            result["pauli"][el] = params.pauli_scale
            result["london"][el] = params.london_scale
            result["static_charge"][el] = params.static_charge
        return result

    def from_hybrid_params(
        self,
        params: Any,
        elements: list[str] | None = None,
    ) -> None:
        """Import parameters from a HybridParameters object."""
        if elements is None:
            elements = sorted(
                set(list(getattr(params, "pauli", {}).keys())
                    + list(getattr(params, "req_radius_offset", {}).keys()))
            )
        for el in elements:
            self.set_pair(
                el,
                radius_offset=getattr(params, "req_radius_offset", {}).get(el, 0.0),
                energy_scale=getattr(params, "req_energy_scale", {}).get(el, 1.0),
                c6_coeff=getattr(params, "c6_coeff", {}).get(el, 0.0),
                pauli_scale=getattr(params, "pauli", {}).get(el, 1.0),
                london_scale=getattr(params, "london", {}).get(el, 1.0),
                static_charge=getattr(params, "static_charge", {}).get(el, 0.0),
            )


def save_element_pair_params(
    registry: AtomTypeRegistry,
    path: str | Path,
) -> None:
    """Save element-pair parameters to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "substrate": registry.substrate,
        "pairs": {k: asdict(v) for k, v in registry.all_pairs().items()},
    }
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)


def load_element_pair_params(path: str | Path) -> AtomTypeRegistry:
    """Load element-pair parameters from JSON."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    registry = AtomTypeRegistry(substrate=payload.get("substrate", "Ag"))
    for key, pair_dict in payload.get("pairs", {}).items():
        pair_dict.pop("metadata", None)
        registry._element_pairs[key] = ElementPairParams(**pair_dict)
    return registry
