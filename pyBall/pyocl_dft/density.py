"""Density projection helpers for :mod:`pyBall.pyocl_dft`.

At this stage the module only defines function skeletons that wire up kernel
arguments; numerical validation and chunked projections will be implemented in
subsequent iterations. The API mirrors the ``pyBall.DFT.oclfft`` exports so
existing scripts can switch hosts by changing the import path.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import pyopencl as cl

from ..OCL.OpenCLBase import OpenCLBase
from .context import OCLContext, ensure_default_context
from . import assets


class OCLDensityProjector(OpenCLBase):
    """Thin OpenCL host wrapping the ``projectOrbDenToGrid_texture`` kernel.

    The class derives from :class:`pyBall.OCL.OpenCLBase` to reuse the
    project-wide infrastructure for context management, kernel compilation, and
    buffer bookkeeping.  It loads the legacy ``myprog.cl`` kernel bundle and
    exposes minimal helper methods for allocating/occupying the buffers needed
    by the AFM density projection workflow.
    """

    def __init__(self, cl_src_dir: Optional[Path] = None, *, nloc: int = 32, device_index: int = 0):
        super().__init__(nloc=nloc, device_index=device_index)
        self._kernel_path = self._resolve_kernel_path(cl_src_dir)
        self.load_program(kernel_path=str(self._kernel_path), bMakeHeaders=True)
        self.kernel_params = {}
        self._init_kernel_metadata()

    # ------------------------------------------------------------------ setup
    @staticmethod
    def _resolve_kernel_path(cl_src_dir: Optional[Path]) -> Path:
        if cl_src_dir is None:
            cl_src_dir = Path(__file__).resolve().parents[2] / "cpp" / "common_resources" / "cl"
        kernel_path = cl_src_dir / "myprog.cl"
        if not kernel_path.is_file():  # pragma: no cover - fail fast when misconfigured
            raise FileNotFoundError(f"Cannot find myprog.cl under {kernel_path.parent}")
        return kernel_path

    def _init_kernel_metadata(self) -> None:
        if "projectOrbDenToGrid_texture" not in self.kernelheaders:
            raise RuntimeError("Kernel headers missing projectOrbDenToGrid_texture; rebuild failed?")
        header = self.kernelheaders["projectOrbDenToGrid_texture"]
        params = self.parse_kernel_header(header)
        self.kernel_params["projectOrbDenToGrid_texture"] = params

    # ---------------------------------------------------------------- buffers
    def ensure_buffers(
        self,
        *,
        atoms: np.ndarray,
        coefs: np.ndarray,
        grid_shape: Tuple[int, int, int],
        buffer_prefix: str = "dens",
    ) -> ProjectionBuffers:
        atoms4 = _ensure_float4(atoms, "atoms")
        coefs4 = _ensure_float4(coefs, "coefs")
        natoms = atoms4.shape[0]
        ngrid = int(np.prod(grid_shape))

        atom_bytes = atoms4.nbytes
        coef_bytes = coefs4.nbytes
        out_bytes = ngrid * 2 * np.dtype(np.float32).itemsize

        self.check_buf(f"{buffer_prefix}_atoms", atom_bytes)
        self.check_buf(f"{buffer_prefix}_coefs", coef_bytes)
        self.check_buf(f"{buffer_prefix}_out", out_bytes)

        self.toGPU(f"{buffer_prefix}_atoms", atoms4)
        self.toGPU(f"{buffer_prefix}_coefs", coefs4)

        return ProjectionBuffers(
            atoms=self.buffer_dict[f"{buffer_prefix}_atoms"],
            coefs=self.buffer_dict[f"{buffer_prefix}_coefs"],
            output=self.buffer_dict[f"{buffer_prefix}_out"],
            grid_shape=grid_shape,
        )

    # ---------------------------------------------------------------- dispatch
    def project_density(
        self,
        bufs: ProjectionBuffers,
        *,
        basis_image: cl.Image,
        grid_origin: Sequence[float],
        grid_vectors: Sequence[Sequence[float]],
        acum_coef: np.ndarray,
    ) -> np.ndarray:
        queue = self.queue
        kernel = self.prg.projectOrbDenToGrid_texture

        n_atoms = np.int32(bufs.atoms.size // (4 * np.dtype(np.float32).itemsize))
        nx, ny, nz = (np.int32(v) for v in bufs.grid_shape)
        grid_origin4 = np.asarray(list(grid_origin) + [0.0], dtype=np.float32)
        dA, dB, dC = grid_vectors
        grid_dA = np.asarray(list(dA) + [0.0], dtype=np.float32)
        grid_dB = np.asarray(list(dB) + [0.0], dtype=np.float32)
        grid_dC = np.asarray(list(dC) + [0.0], dtype=np.float32)
        ngrid = cl.array.vec.make_int4(nx, ny, nz, nx * ny * nz)

        global_size = (int(self.roundUpGlobalSize(int(nx * ny * nz))),)
        local_size = (self.nloc,)

        kernel(
            queue,
            global_size,
            local_size,
            n_atoms,
            bufs.atoms,
            bufs.coefs,
            bufs.output,
            basis_image,
            ngrid,
            grid_origin4,
            grid_dA,
            grid_dB,
            grid_dC,
            acum_coef,
        )

        out_host = np.empty(bufs.grid_shape + (2,), dtype=np.float32)
        cl.enqueue_copy(queue, out_host, bufs.output)
        return out_host


def projector_singleton() -> OCLDensityProjector:
    global _PROJECTOR
    try:
        return _PROJECTOR
    except NameError:
        _PROJECTOR = OCLDensityProjector()
        return _PROJECTOR


@dataclass
class ProjectionBuffers:
    atoms: cl.Buffer
    coefs: cl.Buffer
    output: cl.Buffer
    grid_shape: tuple[int, int, int]


def _ensure_float4(array: np.ndarray, name: str) -> np.ndarray:
    arr = np.asarray(array, dtype=np.float32)
    if arr.ndim != 2 or arr.shape[1] != 4:
        raise ValueError(f"{name} must be an (N,4) array")
    return arr


def prepare_buffers(
    ctx: Optional[OCLContext],
    atoms: np.ndarray,
    coefs: np.ndarray,
    out_shape: Sequence[int],
    *,
    buffer_prefix: str = "dens",
) -> ProjectionBuffers:
    ctx = ctx or ensure_default_context()
    natoms = atoms.shape[0]
    atoms4 = _ensure_float4(atoms, "atoms")
    coefs4 = _ensure_float4(coefs, "coefs")

    atoms_buf = ctx.new_buffer(
        f"{buffer_prefix}_atoms",
        spec=ctx.BufferSpec(name="atoms", shape=atoms4.shape, dtype=np.float32),
    )
    coefs_buf = ctx.new_buffer(
        f"{buffer_prefix}_coefs",
        spec=ctx.BufferSpec(name="coefs", shape=coefs4.shape, dtype=np.float32),
    )
    out_buf = ctx.new_buffer(
        f"{buffer_prefix}_out",
        spec=ctx.BufferSpec(name="output", shape=tuple(out_shape) + (2,), dtype=np.float32),
    )
    ctx.upload(f"{buffer_prefix}_atoms", atoms4)
    ctx.upload(f"{buffer_prefix}_coefs", coefs4)

    return ProjectionBuffers(atoms=atoms_buf, coefs=coefs_buf, output=out_buf, grid_shape=tuple(out_shape))


def project_atoms_dens(
    atoms: np.ndarray,
    coefs: np.ndarray,
    *,
    grid_shape: Sequence[int],
    acum_coef: Optional[Sequence[float]] = None,
    ctx: Optional[OCLContext] = None,
) -> np.ndarray:
    """Project MO coefficients onto a density grid using ``projectOrbDenToGrid_texture``.

    Returns the resulting complex grid as an ``np.complex64`` array. The
    implementation is currently a prototype that allocates temporary buffers
    for each invocation; future revisions will cache buffers and expose chunked
    projections mirroring ``atoms2box`` from the C++ host.
    """

    ctx = ctx or ensure_default_context()
    acum = assets.ensure_acum_coef(tuple(acum_coef) if acum_coef is not None else (0.0, 2.0))

    bufs = prepare_buffers(ctx, atoms, coefs, grid_shape)
    out_array = np.zeros(bufs.grid_shape + (2,), dtype=np.float32)

    program = ctx.build_program("myprog.cl")
    kernel = program.projectOrbDenToGrid_texture

    kernel.set_scalar_arg_dtypes([np.int32, None, None, None, np.int32, np.int32, np.int32, np.float32, np.float32])

    # TODO: wire actual arguments once texture handling and grid descriptors are in place.
    raise NotImplementedError(
        "project_atoms_dens kernel dispatch not implemented yet; groundwork prepared"
    )
