"""Teacher backend registry."""

from .generic_calc import GenericASECalcBackend
from .madsurf import MADSurfTeacherBackend
from .precomputed import PrecomputedTeacherBackend
from .synthetic import SyntheticTeacherBackend

BACKENDS = {
    "madsurf": MADSurfTeacherBackend,
    "generic_calc": GenericASECalcBackend,
    "precomputed": PrecomputedTeacherBackend,
    "synthetic": SyntheticTeacherBackend,
}


def make_teacher_backend(config, model_factory=None):
    if config.kind not in BACKENDS:
        raise KeyError(f"Unsupported teacher backend '{config.kind}'")
    backend_cls = BACKENDS[config.kind]
    return backend_cls(config=config, model_factory=model_factory)
