import os
from dataclasses import dataclass

import numpy as np
import pyopencl as cl

from .OpenCLBase import OpenCLBase

DEFAULT_WORKGROUP_SIZE = 32
MAX_NEIGHBORS          = 8


def _ensure_float4(arr: np.ndarray, *, w_default: float = 0.0) -> np.ndarray:
    data = np.asarray(arr, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array for float4 conversion, got shape {data.shape}")
    if data.shape[1] == 3:
        pad = np.full((data.shape[0], 1), np.float32(w_default), dtype=np.float32)
        data = np.hstack((data, pad))
    if data.shape[1] != 4:
        raise ValueError(f"Expected array with 4 columns after padding, got shape {data.shape}")
    return np.ascontiguousarray(data, dtype=np.float32)


def _ensure_int32(arr: np.ndarray) -> np.ndarray:
    data = np.asarray(arr, dtype=np.int32)
    if data.ndim != 1:
        raise ValueError(f"Expected 1D array for int32 conversion, got shape {data.shape}")
    return np.ascontiguousarray(data, dtype=np.int32)


def _pack_float3(vec) -> np.ndarray:
    data = np.asarray(vec, dtype=np.float32)
    if data.shape != (3,):
        raise ValueError(f"Expected length-3 vector for float3 packing, got shape {data.shape}")
    out = np.zeros(4, dtype=np.float32)
    out[:3] = data
    return out


@dataclass
class LFFParams:
    dt: float = 1e-2
    n_outer: int = 1
    n_inner: int = 8
    efield: tuple = (0.0, 0.0, 0.0)
    bmix: float = 0.0


class LFFSolver(OpenCLBase):
    """pyOpenCL wrapper around the bonded (lff_jacobi) and non-bonded (lff_nb_jacobi) kernels."""

    def __init__(self, *, workgroup_size: int = DEFAULT_WORKGROUP_SIZE, max_neighbors: int = MAX_NEIGHBORS):
        if workgroup_size <= 0 or (workgroup_size & (workgroup_size - 1)) != 0:
            raise ValueError("workgroup_size must be a positive power of two")
        super().__init__(nloc=workgroup_size, device_index=0)

        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "../../cpp/common_resources/cl/LFF.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False, bMakeHeaders=False):
            raise RuntimeError("Failed to load LFF.cl kernel")

        self.max_neighbors = int(max_neighbors)
        if self.max_neighbors != MAX_NEIGHBORS:
            raise ValueError(f"Kernel compiled with MAX_NEIGHBORS={MAX_NEIGHBORS}, got {self.max_neighbors}")

        self.kernelheaders = {
            "lff_jacobi": """__kernel
void lff_jacobi(
    __global const int*    mols,
    __global float4*       apos,
    __global float4*       avel,
    __global const int*    neighs,
    __global const float2* KLs,
    __global const int*    fixed_mask,
    const float3           Efield,
    const float            dt,
    const int              nOuter,
    const int              nInner,
    const float            bMix
)""",
            "lff_nb_jacobi": """__kernel
void lff_nb_jacobi(
    __global const int*    mols,
    __global float4*       apos,
    __global float4*       avel,
    __global float4*       REQHs,
    __global const int*    neighs,
    __global const float2* KLs,
    __global const int*    fixed_mask,
    const float3           Efield,
    const float            dt,
    const int              nOuter,
    const int              nInner,
    const float            bMix
)""",
        }

        self.n_mols = 0
        self.nAtomTot = 0
        self.n_systems = 1
        self.mol_offsets = None
        self.params = LFFParams()
        self.kernel_params = {
            'Efield': _pack_float3(self.params.efield),
            'dt':     np.float32(self.params.dt),
            'nOuter': np.int32  (self.params.n_outer),
            'nInner': np.int32  (self.params.n_inner),
            'bMix':   np.float32(self.params.bmix),
        }
        self.kernel_args: dict[str, tuple] = {}
        self.has_reqhs = False

    # ------------------------------------------------------------------
    # Allocation helpers
    # ------------------------------------------------------------------
    def realloc(self, atom_counts):
        mol_sizes = _ensure_int32(atom_counts)
        if mol_sizes.size == 0: raise ValueError("atom_counts cannot be empty")

        self.n_mols      = mol_sizes.size
        self.mol_offsets = np.zeros(self.n_mols + 1, dtype=np.int32)
        np.cumsum(mol_sizes, out=self.mol_offsets[1:])
        self.nAtomTot    = int(self.mol_offsets[-1])

        float_size   = np.float32().itemsize
        int_size     = np.int32().itemsize
        neigh_stride = self.max_neighbors

        mf = cl.mem_flags

        # Allocate buffers
        self.create_buffer('mols', (self.n_mols + 1) * int_size, mf.READ_ONLY)
        self.create_buffer('apos',       self.nAtomTot * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('avel',       self.nAtomTot * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('REQHs',      self.nAtomTot * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('neighs',     self.nAtomTot * neigh_stride * int_size, mf.READ_ONLY)
        self.create_buffer('KLs',        self.nAtomTot * neigh_stride * 2 * float_size, mf.READ_ONLY)
        self.create_buffer('fixed_mask', self.nAtomTot * int_size, mf.READ_ONLY)

        self.kernel_args.clear()
        self.toGPU('mols', self.mol_offsets)
        self.has_reqhs = False

    # ------------------------------------------------------------------
    # Upload helpers
    # ------------------------------------------------------------------
    def upload_state(self, *, pos, vel, neighs, KLs, fixed_mask=None, REQHs=None):
        if self.nAtomTot == 0: raise RuntimeError("Call realloc(atom_counts) before uploading state")

        pos_dev = _ensure_float4(pos)
        vel_dev = _ensure_float4(vel)

        if pos_dev.shape[0] != self.nAtomTot: raise ValueError(f"pos length {pos_dev.shape[0]} does not match total atoms {self.nAtomTot}")
        if vel_dev.shape[0] != self.nAtomTot: raise ValueError(f"vel length {vel_dev.shape[0]} does not match total atoms {self.nAtomTot}")

        neigh_arr = np.ascontiguousarray(neighs, dtype=np.int32)
        KLs_arr   = np.ascontiguousarray(KLs, dtype=np.float32)
        expected_neigh_shape = (self.nAtomTot, self.max_neighbors)
        if neigh_arr.shape != expected_neigh_shape:
            raise ValueError(f"neighs expected shape {expected_neigh_shape}, got {neigh_arr.shape}")
        expected_KLs_shape = expected_neigh_shape + (2,)
        if KLs_arr.shape != expected_KLs_shape:
            raise ValueError(f"KLs expected shape {expected_KLs_shape}, got {KLs_arr.shape}")

        if fixed_mask is None:
            fixed_arr = np.zeros(self.nAtomTot, dtype=np.int32)
        else:
            fixed_arr = np.ascontiguousarray(fixed_mask, dtype=np.int32)
            if fixed_arr.size != self.nAtomTot: raise ValueError(f"fixed_mask length {fixed_arr.size} does not match total atoms {self.nAtomTot}")

        self.toGPU('apos',   pos_dev)
        self.toGPU('avel',   vel_dev)
        self.toGPU('neighs', neigh_arr)
        self.toGPU('KLs',    KLs_arr)
        self.toGPU('fixed_mask', fixed_arr)
        if REQHs is not None:
            req_dev = _ensure_float4(REQHs)
            if req_dev.shape[0] != self.nAtomTot:
                raise ValueError(f"REQHs length {req_dev.shape[0]} does not match total atoms {self.nAtomTot}")
            self.toGPU('REQHs', req_dev)
            self.has_reqhs = True
        self.queue.finish()

    # ------------------------------------------------------------------
    # Parameter configuration
    # ------------------------------------------------------------------
    def set_params(self, *, dt=None, n_outer=None, n_inner=None, efield=None, bmix=None):
        if dt      is not None: self.params.dt = float(dt)
        if n_outer is not None: self.params.n_outer = int(n_outer)
        if n_inner is not None: self.params.n_inner = int(n_inner)
        if efield  is not None: self.params.efield = tuple(float(x) for x in efield)
        if bmix    is not None: self.params.bmix = float(bmix)

        self.kernel_params['Efield'] = _pack_float3(self.params.efield)
        self.kernel_params['dt']     = np.float32(self.params.dt)
        self.kernel_params['nOuter'] = np.int32  (self.params.n_outer)
        self.kernel_params['nInner'] = np.int32  (self.params.n_inner)
        self.kernel_params['bMix']   = np.float32(self.params.bmix)

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def _build_common_overrides(self) -> dict[str, np.ndarray | np.float32 | np.int32]:
        params = self.params
        return {
            'Efield': _pack_float3(params.efield),
            'dt':     np.float32(params.dt),
            'nOuter': np.int32  (params.n_outer),
            'nInner': np.int32  (params.n_inner),
            'bMix':   np.float32(params.bmix),
        }

    def run_jacobi(self):
        if self.n_mols == 0: raise RuntimeError("Call realloc(atom_counts) before running the kernel")
        overrides = self._build_common_overrides()
        self.kernel_args['lff_jacobi'] = self.generate_kernel_args("lff_jacobi", overrides=overrides)
        global_size = (self.roundUpGlobalSize(self.n_mols * self.nloc),)
        local_size = (self.nloc,)
        self.prg.lff_jacobi(self.queue, global_size, local_size, *self.kernel_args['lff_jacobi'])
        self.queue.finish()

    def run_nb_jacobi(self):
        if self.n_mols == 0: raise RuntimeError("Call realloc(atom_counts) before running the kernel")
        if not self.has_reqhs:
            raise RuntimeError("Upload REQHs parameters (via upload_state(..., REQHs=...)) before calling run_nb_jacobi()")
        overrides = self._build_common_overrides()
        self.kernel_args['lff_nb_jacobi'] = self.generate_kernel_args("lff_nb_jacobi", overrides=overrides)
        global_size = (self.roundUpGlobalSize(self.n_mols * self.nloc),)
        local_size = (self.nloc,)
        self.prg.lff_nb_jacobi(self.queue, global_size, local_size, *self.kernel_args['lff_nb_jacobi'])
        self.queue.finish()

    def run(self, *, include_nonbonded: bool = False):
        if include_nonbonded:
            self.run_nb_jacobi()
        else:
            self.run_jacobi()

    # ------------------------------------------------------------------
    # Download helpers
    # ------------------------------------------------------------------
    def download_state(self):
        #if self.nAtomTot == 0: return {'pos': np.empty((0, 4), dtype=np.float32), 'vel': np.empty((0, 4), dtype=np.float32)}
        pos_host = np.empty((self.nAtomTot, 4), dtype=np.float32)
        vel_host = np.empty((self.nAtomTot, 4), dtype=np.float32)
        self.fromGPU('apos', pos_host)
        self.fromGPU('avel', vel_host)
        self.queue.finish()
        return {'pos': pos_host, 'vel': vel_host}

