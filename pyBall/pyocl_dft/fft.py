"""FFT helpers built on top of :mod:`pyopencl` and :mod:`pyclfft`.

This module provides a thin abstraction similar to the legacy C++ helper
functions ``initFFT`` and ``runFFT`` from ``OCL_GridFF.cpp``. For now only the
plan cache and allocation utilities are implemented; kernel orchestration will
be added once the density/potential pipelines are ported.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import numpy as np

try:
    import pyopencl as cl
    import pyopencl.array as cl_array
    import pyclfft
except ModuleNotFoundError as exc:  # pragma: no cover
    raise ImportError(
        "pyclfft (and pyopencl) are required for pyocl_dft. Install them or "
        "use the legacy pyBall.DFT.oclfft bindings."
    ) from exc

from .context import OCLContext, ensure_default_context


@dataclass
class FFTPlan:
    shape: Tuple[int, ...]
    precision: str = "single"

    def make_forward_inverse(self, queue: cl.CommandQueue) -> Tuple[pyclfft.Plan, pyclfft.Plan]:
        forward = pyclfft.Plan(
            queue.context,
            queue=queue,
            shape=self.shape,
            inplace=True,
            dtype=np.complex64 if self.precision == "single" else np.complex128,
        )
        inverse = pyclfft.Plan(
            queue.context,
            queue=queue,
            shape=self.shape,
            inplace=True,
            dtype=np.complex64 if self.precision == "single" else np.complex128,
            direction=pyclfft.Direction.BACKWARD,
        )
        return forward, inverse


class FFTManager:
    """Cache of clFFT/pyclfft plans keyed by grid shape."""

    def __init__(self, context: OCLContext) -> None:
        self.context = context
        self._plans: Dict[Tuple[int, ...], Tuple[pyclfft.Plan, pyclfft.Plan]] = {}

    def ensure_plan(self, shape: Iterable[int], *, precision: str = "single") -> Tuple[pyclfft.Plan, pyclfft.Plan]:
        key = tuple(int(v) for v in shape)
        if key not in self._plans:
            plan = FFTPlan(key, precision)
            self._plans[key] = plan.make_forward_inverse(self.context.queue)
        return self._plans[key]

    def forward(self, buffer: cl_array.Array) -> None:
        plan, _ = self.ensure_plan(buffer.shape)
        plan.enqueue(field=buffer.data)

    def inverse(self, buffer: cl_array.Array) -> None:
        _, plan = self.ensure_plan(buffer.shape)
        plan.enqueue(field=buffer.data)


_fft_manager: FFTManager | None = None


def manager() -> FFTManager:
    global _fft_manager
    if _fft_manager is None:
        _fft_manager = FFTManager(ensure_default_context())
    return _fft_manager
