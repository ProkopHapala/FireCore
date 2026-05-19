"""Density backend registry."""

from .cube import CubeDensityBackend
from .ml_density import MLDensityBackend
from .pseudo_density import PseudoDensityBackend
from .surface_xyz import SurfaceXYZBackend
from .vasp_volumetric import VaspVolumetricBackend

BACKENDS = {
    "vasp_volumetric": VaspVolumetricBackend,
    "cube": CubeDensityBackend,
    "pseudo_density": PseudoDensityBackend,
    "ml_density": MLDensityBackend,
    "surface_xyz": SurfaceXYZBackend,
}


def make_density_backend(config, surface=None, grid=None):
    if config.kind not in BACKENDS:
        raise KeyError(f"Unsupported density backend '{config.kind}'")
    backend_cls = BACKENDS[config.kind]
    return backend_cls(config=config, surface=surface, grid=grid)
