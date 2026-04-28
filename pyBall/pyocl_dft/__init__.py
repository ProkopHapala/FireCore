"""pyocl_dft
=================

Experimental Python/PyOpenCL host layer that mirrors the public API of
``pyBall.DFT.oclfft`` without relying on the legacy ``libOCL_GridFF.so``
shared library. The implementation is intentionally modular so that the
existing C/C++ path can continue to operate in parallel while we port the
host orchestration logic to Python, as described in
``doc/Topics/AFM/AFM_migration_plan.md``.

Only a subset of the original functions is functional at this stage. The
module raises ``NotImplementedError`` for unfinished features to ensure that
missing pieces fail loudly rather than silently degrading behaviour (see the
project coding guidelines).
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np

from .context import (
    OCLContext,
    context_singleton,
    ensure_default_context,
    release_default_context,
)
from . import assets

__all__ = [
    "OCLContext",
    "assets",
    "init",
    "release",
    "get_context",
    "load_kernels",
    "require_program",
]


def init(
    cl_src_dir: Optional[str] = None,
    platform: Optional[int] = None,
    device: Optional[int] = None,
    build_options: Optional[Sequence[str]] = None,
) -> OCLContext:
    """Initialise (or reconfigure) the shared :class:`OCLContext`.

    Parameters
    ----------
    cl_src_dir
        Directory that contains ``myprog.cl``, ``GridFF.cl`` and ``relax.cl``.
        Defaults to the project path used by the C++ host.
    platform, device
        Optional indices for selecting a specific OpenCL platform/device.
    build_options
        Additional options passed to ``cl.Program.build``.
    """

    ctx = context_singleton()
    ctx.configure(
        cl_src_dir=cl_src_dir,
        platform_index=platform,
        device_index=device,
        build_options=build_options,
    )
    return ctx


def release() -> None:
    """Dispose of the shared context, releasing OpenCL resources."""

    release_default_context()


def get_context() -> OCLContext:
    """Return the lazily created default context."""

    return ensure_default_context()


def load_kernels(names: Iterable[str]) -> None:
    """Pre-build one or more OpenCL programs."""

    ctx = ensure_default_context()
    for name in names:
        ctx.build_program(name)


def require_program(name: str) -> "pyopencl.Program":  # type: ignore[name-defined]
    """Ensure *name* is compiled and return the corresponding program."""

    ctx = ensure_default_context()
    return ctx.build_program(name)


# === Unimplemented API surface =================================================


def _raise_unimplemented(func: str) -> None:
    raise NotImplementedError(
        f"pyocl_dft.{func}() is not implemented yet; continue using "
        "pyBall.DFT.oclfft for this call or contribute the missing port."
    )


def initFFT(*args, **kwargs):  # noqa: N802 (legacy API signature)
    _raise_unimplemented("initFFT")


def convolve(*args, **kwargs):
    _raise_unimplemented("convolve")


def poisson(*args, **kwargs):
    _raise_unimplemented("poisson")


def gradient(*args, **kwargs):
    _raise_unimplemented("gradient")


def projectAtomsDens(*args, **kwargs):
    _raise_unimplemented("projectAtomsDens")


def projectDenmat(*args, **kwargs):
    _raise_unimplemented("projectDenmat")


def evalLJC_QZs(*args, **kwargs):
    _raise_unimplemented("evalLJC_QZs")


def relaxStrokesTilted(*args, **kwargs):
    _raise_unimplemented("relaxStrokesTilted")
