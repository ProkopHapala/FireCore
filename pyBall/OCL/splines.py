import sys
import os
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import matplotlib.pyplot as plt

def make_inds_pbc(n): 
    return np.array([
    [ 0, 1,   2,   3   ],
    [ 0, 1,   2,   3-n ],
    [ 0, 1,   2-n, 3-n ],
    [ 0, 1-n, 2-n, 3-n ]], dtype=np.int32 );

def roundup_global_size(global_size, local_size):
    remainder = global_size % local_size
    if remainder == 0: return global_size
    return global_size + (local_size - remainder)


class OCLSplines:
    def __init__(self, nloc = 32 ):
        self.nloc= nloc
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        
        try:
            with open('../../cpp/common_resources/cl/splines.cl', 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print( "OCLSplines() called from path=", os.getcwd() )
            print(f"Error compiling OpenCL program: {e}")

    def sample3D_comb3(self, g0, dg, ng, Eg, ps, C):
        """
        __kernel void sample3D_comb3(
            const float3 g0,
            const float3 dg,
            const int3 ng,
            __global const float3* Eg,
            const int n,
            __global const float3* ps,
            __global float4* fes,
            const float3 C,
            __global int4* xqs,
            __global int4* yqs
        ) {
        """
        n = len(ps)
        g0 = np.array(g0, dtype=np.float32)
        dg = np.array(dg, dtype=np.float32)
        ng = np.array(ng, dtype=np.int32)
        C  = np.array(C,  dtype=np.float32)
        
        Eg_buf  = cl_array.to_device(self.queue, Eg.astype(np.float32))
        ps_buf  = cl_array.to_device(self.queue, ps.astype(np.float32))
        fes_buf = cl_array.empty(    self.queue, (n, 4), dtype=np.float32)
        xqs = make_inds_pbc(ng[0])
        yqs = make_inds_pbc(ng[1])

        nG = roundup_global_size( n, self.nloc )
        
        self.prg.sample3D_comb3(self.queue, (nG,), (self.nloc,),      g0, dg, ng, Eg_buf.data, np.int32(n),  ps_buf.data, fes_buf.data, C, xqs, yqs )
        
        return fes_buf.get()

    def sample1D_pbc(self, g0, dg, ng, Gs, ps):
        """
        __kernel void sample1D_pbc(
            const float g0,
            const float dg,
            const int ng,
            __global const float* Gs,
            const int n,
            __global const float* ps,
            __global float2* fes,
            __global int4* xqs
        ) {
        """
        n   = len(ps)
        g0  = np.float32(g0)
        dg  = np.float32(dg)
        ng  = np.int32(ng)
        #xqs = make_inds_pbc(ng)
        
        Gs_buf  = cl_array.to_device(self.queue, Gs.astype(np.float32))
        ps_buf  = cl_array.to_device(self.queue, ps.astype(np.float32))
        fes_buf = cl_array.empty(self.queue, (n, 2), dtype=np.float32)
        fes_buf = cl_array.empty(self.queue, (n, 2), dtype=np.float32)
        
        nG = roundup_global_size( n, self.nloc )

        self.prg.sample1D_pbc(self.queue, (nG,), (self.nloc,),        g0, dg, ng, Gs_buf.data, np.int32(n),   ps_buf.data, fes_buf.data )
        
        return fes_buf.get()

# Example usage
if __name__ == "__main__":
    ocl_splines = OCLSplines()

    # 3D example
    # g0_3d = [0, 0, 0]
    # dg_3d = [0.1, 0.1, 0.1]
    # ng_3d = [64, 64, 64]
    # Eg_3d = np.random.rand(64, 64, 64, 3).astype(np.float32)
    # ps_3d = np.random.rand(1000, 3).astype(np.float32) * 6.3
    # C_3d = [1, 1, 1]

    # result_3d = ocl_splines.sample3D_comb3(g0_3d, dg_3d, ng_3d, Eg_3d, ps_3d, C_3d)

    # plt.figure(figsize=(10, 8))
    # plt.scatter(ps_3d[:, 0], ps_3d[:, 1], c=result_3d[:, 3], cmap='viridis')
    # plt.colorbar(label='Energy')
    # plt.title('3D Spline Interpolation')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.show()

    # 1D example
    g0_1d = 0
    dg_1d = 0.1
    ng_1d = 100
    Gs_1d = np.sin(np.linspace(0, 2*np.pi, ng_1d)).astype(np.float32)
    ps_1d = np.linspace(0, 10, 1000).astype(np.float32)

    result_1d = ocl_splines.sample1D_pbc(g0_1d, dg_1d, ng_1d, Gs_1d, ps_1d)

    plt.figure(figsize=(10, 6))
    plt.plot(ps_1d, result_1d[:, 0], label='Interpolated')
    plt.plot(np.linspace(0, 10, ng_1d), Gs_1d, 'ro', label='Original')
    plt.title('1D Spline Interpolation with PBC')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()
