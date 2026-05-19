"""Substrate-class factory for the GridFF-JAX framework.

A *substrate class* groups every config choice that must change together
when you move from a metallic substrate (Ag, Au, Cu) to an ionic one
(NaCl, KBr) or any future class. Picking a class deterministically sets:

  - ``GridConfig.builder_mode``        — which substrate-field builder runs
  - ``FeatureToggles``                 — which channels are assembled
  - ``DensityBackendConfig.kind``      — default density source
  - ``TeacherBackendConfig.teacher_tile`` — periodic-image tiling default
  - ``HybridModelConfig.use_image``    — image-charge plane (metals only)

Why this exists: previously, ``benchmark_ag_6d.py`` hardcoded
``builder_mode='metal_density_plq'`` and assumed LOCPOT-derived work-function
physics. That's fine for coin metals, wrong for ionic salts. NaCl needs the
``parity_core`` Ewald-slab path. Mixing those by hand at the CLI level is
error-prone, so this module owns the mapping.

Usage
-----
    from pyBall.gridff_jax.substrate_classes import apply_substrate_class

    cfg = RunConfig()
    apply_substrate_class(cfg, "ionic")    # mutates cfg in place
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from typing import Any

from .config import FeatureToggles, RunConfig


# ---------------------------------------------------------------------------
#  Per-class defaults
# ---------------------------------------------------------------------------
#  Channel semantics (see hybrid_energy/model.py and substrate_fields.py):
#    use_density_charge       — derive Pauli/London magnitudes from density
#    use_locpot               — use VASP LOCPOT for image-plane Coulomb
#                               (only meaningful for conductors)
#    use_ct_qeq               — charge-transfer QEQ
#    use_image_charge         — image-charge plane (metal screening)
#    use_image_charge_fixed   — image plane with fixed z (metal)
#    use_reactive_grid        — element-wise reactive shells
#    use_teacher_residual     — keep teacher residual channel for QC

SUBSTRATE_CLASS_DEFAULTS: dict[str, dict[str, Any]] = {
    # ---- metals (Ag, Au, Cu, ...) ---------------------------------------
    "metal": {
        "description": "Conductive substrate (Ag, Au, Cu, ...). "
                       "Uses CHGCAR density-derived PLQ + LOCPOT image plane.",
        "builder_mode": "metal_density_plq",
        "toggles": FeatureToggles(
            use_density_charge=True,
            use_locpot=True,
            use_ct_qeq=False,
            use_image_charge=False,
            use_image_charge_fixed=True,
            use_reactive_grid=False,
            use_teacher_residual=True,
        ),
        "density_kind_default": "vasp_volumetric",
        "teacher_tile_default": (0, 0),      # 'auto'
        # Hybrid model flags (NOTE: at runtime the model reads FeatureToggles,
        # not these. These exist for use_req_plq and informational provenance.)
        "use_qeq": False,
        "use_image": True,
        "use_reactive_grid": False,
        # ---- Physics knobs that belong with the class, not the driver ----
        # Metals benefit from London Fermi damping near the surface (avoids
        # the divergent r→0 1/r⁶ tail), pairwise C₆ refinement, and Stage 2
        # raw 1/r⁶ refit on top of REQ-PLQ.
        "use_pairwise_c6": True,
        "london_damping_d0": 3.70,
        "london_damping_width": 0.35,
        "use_raw_r6_stage2": True,
        "fit_c6_coeff": True,
        "two_stage_c6": True,
        "fit_req_radius_offset": True,
        "fit_req_energy_scale": True,
    },
    # ---- ionic salts (NaCl, KBr, CsCl, ...) ------------------------------
    "ionic": {
        "description": "Ionic substrate (NaCl, KBr, CsCl, ...). "
                       "Uses point-charge Ewald-slab Coulomb (no work function).",
        "builder_mode": "parity_core",
        "toggles": FeatureToggles(
            use_density_charge=True,
            use_locpot=False,
            use_ct_qeq=False,
            use_image_charge=False,
            use_image_charge_fixed=False,
            use_reactive_grid=False,
            use_teacher_residual=True,
        ),
        "density_kind_default": "surface_xyz",
        "teacher_tile_default": (1, 1),
        "use_qeq": False,
        "use_image": False,
        "use_reactive_grid": False,
        # ---- Ionic physics: no Fermi damping (no metallic screening),
        # no pairwise C₆ refinement (dispersion is small vs the ~eV ionic
        # potential), no raw 1/r⁶ Stage 2 (would just absorb Coulomb noise).
        # REQ-PLQ fitting stays on for the short-range part.
        "use_pairwise_c6": False,
        "london_damping_d0": 0.0,
        "london_damping_width": 0.5,
        "use_raw_r6_stage2": False,
        "fit_c6_coeff": False,
        "two_stage_c6": False,
        "fit_req_radius_offset": True,
        "fit_req_energy_scale": True,
    },
}


# ---------------------------------------------------------------------------
#  Public API
# ---------------------------------------------------------------------------
def available_substrate_classes() -> list[str]:
    return sorted(SUBSTRATE_CLASS_DEFAULTS.keys())


def apply_substrate_class(cfg: RunConfig, name: str,
                          override_density_kind: bool = True,
                          override_teacher_tile: bool = True,
                          verbose: bool = True) -> RunConfig:
    """Mutate ``cfg`` to match the named substrate class.

    Parameters
    ----------
    cfg : RunConfig
        The config to mutate.
    name : str
        Substrate class name (``metal``, ``ionic``).
    override_density_kind : bool
        If True (default), set ``cfg.density_backend.kind`` to the class
        default. Set to False if the user has explicitly chosen a density
        backend via CLI (so their choice wins).
    override_teacher_tile : bool
        Same logic for the teacher tile.
    verbose : bool
        Print what was changed.

    Returns
    -------
    RunConfig
        The same object, mutated.

    Raises
    ------
    KeyError if ``name`` is not a known substrate class.
    """
    if name not in SUBSTRATE_CLASS_DEFAULTS:
        raise KeyError(
            f"Unknown substrate class '{name}'. "
            f"Available: {available_substrate_classes()}"
        )
    spec = SUBSTRATE_CLASS_DEFAULTS[name]

    if verbose:
        print(f"[substrate-class] applying '{name}': {spec['description']}", flush=True)

    cfg.grid.builder_mode = spec["builder_mode"]
    cfg.toggles = deepcopy(spec["toggles"])
    # NOTE: hybrid_model flags are informational; FeatureToggles drives
    # runtime energy assembly (see hybrid_energy/model.py). We still set them
    # for provenance / sidecar export consistency.
    cfg.hybrid_model.use_qeq = spec["use_qeq"]
    cfg.hybrid_model.use_image = spec["use_image"]
    cfg.hybrid_model.use_reactive_grid = spec["use_reactive_grid"]

    # ---- Physics knobs owned by the class ----
    cfg.grid.use_pairwise_c6 = bool(spec["use_pairwise_c6"])
    cfg.grid.london_damping_d0 = float(spec["london_damping_d0"])
    cfg.grid.london_damping_width = float(spec["london_damping_width"])
    cfg.training.fit_c6_coeff = bool(spec["fit_c6_coeff"])
    cfg.training.two_stage_c6 = bool(spec["two_stage_c6"])
    cfg.training.fit_req_radius_offset = bool(spec["fit_req_radius_offset"])
    cfg.training.fit_req_energy_scale = bool(spec["fit_req_energy_scale"])

    if override_density_kind:
        cfg.density_backend.kind = spec["density_kind_default"]
    if override_teacher_tile:
        cfg.teacher_backend.teacher_tile = spec["teacher_tile_default"]

    if verbose:
        print(f"  builder_mode      = {cfg.grid.builder_mode}", flush=True)
        print(f"  toggles           = {asdict(cfg.toggles)}", flush=True)
        print(f"  density.kind      = {cfg.density_backend.kind}", flush=True)
        print(f"  teacher_tile      = {cfg.teacher_backend.teacher_tile}", flush=True)
        print(f"  use_pairwise_c6   = {cfg.grid.use_pairwise_c6}", flush=True)
        print(f"  london_damp_d0    = {cfg.grid.london_damping_d0}", flush=True)
        print(f"  two_stage_c6      = {cfg.training.two_stage_c6}", flush=True)
        print(f"  fit_c6_coeff      = {cfg.training.fit_c6_coeff}", flush=True)

    return cfg


def get_substrate_class_recipe(name: str) -> dict[str, Any]:
    """Return the (deepcopy of) the recipe dict for callers that need to query
    fields not exposed via RunConfig — e.g. ``use_raw_r6_stage2``."""
    if name not in SUBSTRATE_CLASS_DEFAULTS:
        raise KeyError(f"Unknown substrate class '{name}'")
    return deepcopy(SUBSTRATE_CLASS_DEFAULTS[name])
