"""
RndNormal.py  –  pyOpenCL wrapper for RndNormal.cl normal-distribution generators.

Usage:
    rng = RndNormal(N=2**20)
    results, timings = rng.run_all(seed=42, n_warmup=2, n_bench=10)
"""

import os
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array

# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_CL_FILE = os.path.join(_HERE, "RndNormal.cl")

KERNEL_NAMES = [
    "box_muller_32",
    "box_muller_64",
    "inverse_cdf",
    "marsaglia_polar",
    "global_buffer",
    "shared_buffer",
]

# Buffer kernels that require a pre-generated normal pool
_BUFFER_KERNELS = {"global_buffer", "shared_buffer"}

def _select_gpu(prefer_nvidia=True):
    """Return (ctx, queue) for the best available GPU."""
    platforms = cl.get_platforms()
    if prefer_nvidia:
        for p in platforms:
            if 'NVIDIA' in p.name:
                devs = p.get_devices(cl.device_type.GPU)
                if devs:
                    ctx = cl.Context([devs[0]])
                    print(f"[RndNormal] Using GPU: {devs[0].name}  ({p.name})")
                    return ctx
    for p in platforms:
        devs = p.get_devices(cl.device_type.GPU)
        if devs:
            ctx = cl.Context([devs[0]])
            print(f"[RndNormal] Using GPU: {devs[0].name}  ({p.name})")
            return ctx
    print("[RndNormal] No GPU found, falling back to CPU/default")
    return cl.create_some_context(interactive=False)


class RndNormal:
    """
    Manages OpenCL context, kernel compilation, normal-buffer, and kernel dispatch.

    Parameters
    ----------
    N           : number of samples to generate per call
    buffer_size : size of pre-generated normal pool (should be power-of-2 for fast % )
    lds_size    : local work-group size for shared_buffer kernel (power-of-2)
    prefer_nvidia : prefer NVIDIA platform if available
    build_opts  : extra build flags (list of str)
    """

    def __init__(self, N=2**20, buffer_size=2**20, lds_size=256,
                 prefer_nvidia=True, build_opts=None):
        self.N           = N
        self.buffer_size = buffer_size
        self.lds_size    = lds_size

        self.ctx   = _select_gpu(prefer_nvidia)
        self.queue = cl.CommandQueue(self.ctx,
                     properties=cl.command_queue_properties.PROFILING_ENABLE)

        self._build(build_opts or ["-cl-fast-relaxed-math", "-cl-mad-enable"])
        self._upload_normal_buffer()

    # ------------------------------------------------------------------
    def _build(self, extra_opts):
        src = open(_CL_FILE).read()
        defines = (
            f"#define BUFFER_SIZE {self.buffer_size}\n"
            f"#define LDS_SIZE    {self.lds_size}\n"
            f"#define INV_SQRT2   0.7071067811865475f\n"
            f"#define M_PI_F      3.14159265358979323846f\n"
        )
        full_src = defines + src
        opts = " ".join(extra_opts)
        print("[RndNormal] Building OpenCL program...")
        self.program = cl.Program(self.ctx, full_src).build(options=opts)
        print("[RndNormal] Build OK")

    def _upload_normal_buffer(self):
        """Pre-generate a high-quality normal pool on CPU and upload to GPU."""
        cpu_buf = np.random.normal(0.0, 1.0, self.buffer_size).astype(np.float32)
        self.normal_buf = cl_array.to_device(self.queue, cpu_buf)
        print(f"[RndNormal] Uploaded normal buffer ({self.buffer_size} floats, "
              f"{cpu_buf.nbytes/1e6:.1f} MB)")

    # ------------------------------------------------------------------
    def _kernel_args(self, name, seed):
        """Return (extra_args, local_size) for a given kernel name."""
        s = np.uint32(seed)
        if name in _BUFFER_KERNELS:
            extra = [self.normal_buf.data, s]
            ls    = (self.lds_size,) if name == "shared_buffer" else None
        else:
            extra = [s]
            ls    = None
        return extra, ls

    def run_kernel(self, name, seed=42, n_warmup=1, n_bench=1):
        """
        Run a single kernel, return (samples_float32, elapsed_ms).
        elapsed_ms is the GPU event time of the *last* bench iteration.
        """
        kern       = getattr(self.program, name)
        extra, ls  = self._kernel_args(name, seed)
        N          = self.N
        gs         = N if ls is None else ((N + ls[0] - 1) // ls[0]) * ls[0]

        out = cl_array.empty(self.queue, (gs,), dtype=np.float32)

        for _ in range(n_warmup):
            kern(self.queue, (gs,), ls, out.data, *extra)
        self.queue.finish()

        evt = kern(self.queue, (gs,), ls, out.data, *extra)
        evt.wait()
        elapsed = (evt.profile.end - evt.profile.start) * 1e-6  # ms

        result = out.get()[:N]
        return result, elapsed

    def run_all(self, seed=42, n_warmup=2, n_bench=1, kernels=None):
        """
        Run all (or selected) kernels.
        Returns:
            results : dict  name -> np.ndarray float32 (N,)
            timings : dict  name -> float  (ms, GPU event time)
        """
        names   = kernels or KERNEL_NAMES
        results = {}
        timings = {}
        print(f"\n[RndNormal] Running {len(names)} kernels  N={self.N}")
        for name in names:
            data, t = self.run_kernel(name, seed=seed, n_warmup=n_warmup)
            results[name] = data
            timings[name] = t
            print(f"  {name:22s}  {t:7.3f} ms   "
                  f"{self.N / t * 1e-6:.0f} M samples/s")
        return results, timings

    # ------------------------------------------------------------------
    def bench(self, name, seed=42, n_warmup=3, n_iters=50):
        """
        Wall-clock benchmark over n_iters kernel launches.
        Returns mean elapsed time in ms.
        """
        import time
        kern       = getattr(self.program, name)
        extra, ls  = self._kernel_args(name, seed)
        N          = self.N
        gs         = N if ls is None else ((N + ls[0] - 1) // ls[0]) * ls[0]
        out        = cl_array.empty(self.queue, (gs,), dtype=np.float32)

        for _ in range(n_warmup):
            kern(self.queue, (gs,), ls, out.data, *extra)
        self.queue.finish()

        t0 = time.perf_counter()
        for i in range(n_iters):
            kern(self.queue, (gs,), ls, out.data, *[extra[0]] + [np.uint32(seed + i)] if name in _BUFFER_KERNELS else [np.uint32(seed + i)])
        self.queue.finish()
        return (time.perf_counter() - t0) / n_iters * 1e3  # ms


# ---------------------------------------------------------------------------
# Quick CLI test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    rng = RndNormal(N=2**20, buffer_size=2**20, lds_size=256)
    results, timings = rng.run_all(seed=42, n_warmup=3)

    print("\n[RndNormal] Sample stats (mean, std):")
    for name, data in results.items():
        print(f"  {name:22s}  mean={np.mean(data):+.4f}  std={np.std(data):.4f}")
