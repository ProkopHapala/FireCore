import pyopencl as cl
import pyopencl.array
import numpy as np

class Conv3D:
    def __init__(self, ctx, queue, tile_size=(8, 8, 8)):
        """
        :param ctx: PyOpenCL context
        :param queue: PyOpenCL queue
        :param tile_size: Tuple (lx, ly, lz) for local workgroup size
        """
        self.ctx = ctx
        self.queue = queue
        self.ls = tile_size

        from pathlib import Path
        kernel_path = Path(__file__).resolve().parent / "cl" / "covolve.cl"
        KERNEL_SOURCE = kernel_path.read_text()
        
        # Compile with dynamic tile sizes
        build_options = f"-D LS_X={self.ls[0]} -D LS_Y={self.ls[1]} -D LS_Z={self.ls[2]}"
        self.prg = cl.Program(ctx, KERNEL_SOURCE).build(options=build_options)
        self.kernel = self.prg.Convolution3D_General

    def run(self, input_arr, output_arr=None, G0_arr=None,  coef_mul=1.0, coef_add=1.0,  pbc=(True, True, True)):
        """
        :param input_arr: numpy array or cl.Buffer (3D, float32)
        :param pbc: tuple of bools for Periodic Boundary Conditions (x, y, z)
        """
        # Ensure data types
        ns = np.int32(input_arr.shape[::-1]) # {nx, ny, nz} reversed for OpenCL (x is fast)
        nx, ny, nz = input_arr.shape[2], input_arr.shape[1], input_arr.shape[0] # Numpy is Z, Y, X
        
        # Buffers
        if isinstance(input_arr, np.ndarray):
            gs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=input_arr.astype(np.float32))
        else:
            gs_buf = input_arr

        if output_arr is None:
            out_buf = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, input_arr.nbytes)
        elif isinstance(output_arr, np.ndarray):
             out_buf = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, output_arr.nbytes)
        else:
            out_buf = output_arr

        g0_ptr = None
        if G0_arr is not None:
            if isinstance(G0_arr, np.ndarray):
                g0_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=G0_arr.astype(np.float32))
                g0_ptr = g0_buf
            else:
                g0_ptr = G0_arr

        # B-Spline Weights (B0=2/3, B1=1/6)
        # Note: In 3D, weight is product of 1D weights.
        # w_center = B0*B0*B0
        # w_face   = B0*B0*B1 (1 neighbor dim)
        # w_edge   = B0*B1*B1 (2 neighbor dims)
        # w_corner = B1*B1*B1 (3 neighbor dims)
        
        b0 = 2.0/3.0
        b1 = 1.0/6.0
        
        w_center = b0**3
        w_face   = (b0**2) * b1
        w_edge   = b0 * (b1**2)
        w_corner = b1**3
        
        weights = np.array([w_center, w_face, w_edge, w_corner], dtype=np.float32)
        
        # Boundary params
        pbc_int = np.array([1 if p else 0 for p in pbc] + [0], dtype=np.int32)
        coefs = np.array([coef_mul, coef_add], dtype=np.float32)
        dims = np.array([nx, ny, nz, 0], dtype=np.int32) # x, y, z

        # Calculate padded Global Size
        # Global size must be multiple of local size
        def align(n, block):
            return ((n + block - 1) // block) * block

        g_dim_x = align(nx, self.ls[0])
        g_dim_y = align(ny, self.ls[1])
        g_dim_z = align(nz, self.ls[2])
        
        global_size = (g_dim_x, g_dim_y, g_dim_z)
        local_size  = self.ls

        # Enqueue
        self.kernel(self.queue, global_size, local_size, dims, gs_buf, g0_ptr, out_buf, weights, pbc_int, coefs)

        # Read back if needed
        if isinstance(output_arr, np.ndarray):
            cl.enqueue_copy(self.queue, output_arr, out_buf)
            return output_arr
        elif output_arr is None:
            res = np.empty_like(input_arr, dtype=np.float32)
            cl.enqueue_copy(self.queue, res, out_buf)
            return res
        else:
            return out_buf