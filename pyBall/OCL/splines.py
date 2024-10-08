import sys
import os
import numpy as np
import ctypes
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

bytePerFloat = 4

# cl_Mat3_dtype = np.dtype([
#     ('a', cltypes.float4),
#     ('b', cltypes.float4),
#     ('c', cltypes.float4)
# ], align=True)

# cl_Mat3_dtype = np.dtype([
#     ('a', np.float32, 4),
#     ('b', np.float32, 4),
#     ('c', np.float32, 4)
# ], align=True)

# class cl_Mat3(ctypes.Structure):
#     _fields_ = [
#         ("a", cltypes.float4),
#         ("b", cltypes.float4),
#         ("c", cltypes.float4)
#     ]

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

def roundup_global_size_3d( global_size,  local_size):
    return (
         roundup_global_size(global_size[0], local_size[0]),
         roundup_global_size(global_size[1], local_size[1]),
         roundup_global_size(global_size[2], local_size[2])
     )
def get_cl_info( device ):
    
    '''
    also use:
        nvidia-smi -q
    '''

    # Example for NVIDIA GPUs:
    # Pascal: 128 CUDA cores per SM
    # Turing: 64 CUDA cores per SM
    # Ampere: 128 CUDA cores per SM  10496 cores / 82 compute units

    print(f"Device Name: {device.name}")
    print(f"Max Compute Units: {device.max_compute_units}")
    print(f"Max Work Group Size: {device.max_work_group_size}")
    print(f"Global Memory Size: {device.global_mem_size / (1024*1024)} MB")
    print(f"Local Memory Size: {device.local_mem_size / 1024} KB")
    print(f"Max Clock Frequency: {device.max_clock_frequency} MHz")
              
    # Get local memory characteristics
    granularity      = cl.characterize.local_memory_access_granularity(device)
    bank_count       = cl.characterize.local_memory_bank_count(device)
    usable_local_mem = cl.characterize.usable_local_mem_size(device)
    # Print results
    print(f"Local Memory Access Granularity: {granularity} bytes")
    print(f"Number of Local Memory Banks: {bank_count}")
    print(f"Usable Local Memory Size: {usable_local_mem} bytes")

    # Retrieve various characteristics
    fast_math_options = cl.characterize.get_fast_inaccurate_build_options(device)
    simd_group_size   = cl.characterize.get_simd_group_size(device, 4)  # Assuming float4
    double_support    = cl.characterize.has_amd_double_support(device)
    #src_cache_support = cl.characterize.has_src_build_cache(device)

    # Print results
    print(f"Fast Math Options: {fast_math_options}")
    print(f"SIMD Group Size: {simd_group_size}")
    print(f"Double Precision Support: {double_support}")
    #print(f"Source Build Cache Support: {src_cache_support}")

def local_memory_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return device.local_mem_size /( sp_per_cu/local_size )

def local_memory_float_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/4

def local_memory_float4_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/(4*4)

class OCLSplines:

    def __init__(self, nloc=32 ):
        self.nloc  = nloc
        self.ctx   = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)

        get_cl_info( self.ctx.devices[0] )

        local_size = 64
        #print( " local_memory_per_workgroup() size=", local_size, " __local []  ", local_memory_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " Byte " );
        print( " local_memory_per_workgroup() size=", local_size, " __local []  ", local_memory_float_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " float32 " );

        try:
            with open('../../cpp/common_resources/cl/splines.cl', 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print( "OCLSplines() called from path=", os.getcwd() )
            print(f"Error compiling OpenCL program: {e}")
            exit(0)


        self.atoms_buff  = None
        self.REQs_buff       = None
        self.V_Paul_buff     = None
        self.V_Lond_buff     = None

    def prepare_sample3D(self, g0, dg, ng, Eg ):
        g0 = np.array(g0+(0,), dtype=np.float32)
        dg = np.array(dg+(0,), dtype=np.float32)
        ng = np.array(ng+(0,), dtype=np.int32)
        self.grid3D_shape = (g0, dg, ng)
        self.E3D_buf      = cl_array.to_device(self.queue, Eg.astype(np.float32) )

    def sample3D_comb(self, ps, C):
        """
        __kernel void sample3D_comb(
            const float4 g0,
            const float4 dg,
            const int4 ng,
            __global const float4* Eg,
            const int n,
            __global const float4* ps,
            __global float4* fes,
            const float4 C
        ) {
        """
        n = len(ps)

        #self.prepare_sample3D( g0, dg, ng, Eg )

        C       = np.array(C , dtype=np.float32)
        ps_buf  = cl_array.to_device(self.queue, ps.astype(np.float32) )
        fes_buf = cl_array.empty(    self.queue, (n, 4), dtype=np.float32)

        nG = roundup_global_size( n, self.nloc )
        (g0, dg, ng) = self.grid3D_shape

        # print("g0", g0)
        # print("dg", dg)
        # print("ng", ng)
        self.prg.sample3D_comb(self.queue, (nG,), (self.nloc,),  g0, dg, ng, self.E3D_buf.data, np.int32(n),  ps_buf.data, fes_buf.data, C )
        fe = fes_buf.get()
        return fe

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
    

    # ============ Fitting


    def allocate_cl_buffers(self, ns=None, E0=None, nByte=4):

        print(f"OCLSplines::allocate_cl_buffers().1 Queue: {self.queue}, Context: {self.ctx}")

        if ns is None:
            ns = E0.shape
        ntot = np.int32(ns[0] * ns[1] * ns[2])

        self.G0_buf  = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=ntot *nByte)   # Reference
        self.Gs_buf  = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)
        self.fGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # variational derivative
        self.dGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # fitting error
        self.vGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # velocity of Gs update

        #self.E3D_tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, cl.ImageFormat(cl.channel_order.R, cl.channel_type.FLOAT),  shape=ns  )
        
        if E0 is not None:
            cl.enqueue_copy(self.queue, self.G0_buf,  E0.astype(np.float32))
            cl.enqueue_copy(self.queue, self.Gs_buf,  E0.astype(np.float32))
            #cl.enqueue_copy(self.queue, self.E3D_tex, E0.astype(np.float32))
        
        self.ns = np.array(ns, dtype=np.int32)

        print(f"OCLSplines::allocate_cl_buffers().2 Queue: {self.queue}, Context: {self.ctx}")

    def prepare_fit_buffers(self, g0, dg, ng, Eg):
        img_format = cl.ImageFormat(cl.channel_order.R, cl.channel_type.FLOAT)
        #self.E3D_buf = cl_array.to_device(self.queue, Eg.astype(np.float32))
        self.E3D_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE | cl.mem_flags.COPY_HOST_PTR  )
        self.E3D_tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY   | cl.mem_flags.COPY_HOST_PTR,  img_format ,  shape=ng, hostbuf=Eg.astype(np.float32))


    def BsplineConv3D(self, Gs, G0):
        ns = np.array(Gs.shape, dtype=np.int32)
        Gs_buf  = cl_array.to_device(self.queue, Gs.astype(np.float32))
        G0_buf  = cl_array.to_device(self.queue, G0.astype(np.float32))
        out_buf = cl_array.empty(self.queue, Gs.shape, dtype=np.float32)

        self.prg.BsplineConv3D(self.queue, Gs.shape, None,     np.int32(ns), Gs_buf.data, G0_buf.data, out_buf.data)
        
        return out_buf.get()

    def BsplineConv3D_tex(self, G0):
        ns      = np.array(self.E3D_buf.shape, dtype=np.int32)
        G0_buf  = cl_array.to_device(self.queue, G0.astype(np.float32))
        out_buf = cl_array.empty(self.queue, self.E3D_buf.shape, dtype=np.float32)

        self.prg.BsplineConv3D_tex(self.queue, self.E3D_buf.shape, None,   np.int32(ns), self.E3D_tex, G0_buf.data, out_buf.data)
        
        return out_buf.get()

    def move(self, p, v, f, par):
        ntot = len(p)
        p_buf = cl_array.to_device(self.queue, p.astype(np.float32))
        v_buf = cl_array.to_device(self.queue, v.astype(np.float32))
        f_buf = cl_array.to_device(self.queue, f.astype(np.float32))

        self.prg.move(self.queue, (ntot,), None,  np.int32(ntot), p_buf.data, v_buf.data, f_buf.data, np.float32(par))
        
        return p_buf.get(), v_buf.get()

    def fit3D(self, ns=None, E_ref=None, nmaxiter=100, dt=0.3, Ftol=1e-9, cdamp=0.98, bAlloc=True, nPerStep=10, bConvTrj=False ):
        
        print(f"OCLSplines::fit3D().1 Queue: {self.queue}, Context: {self.ctx}")

        if bAlloc:
            ns  = E_ref.shape
            ns_bak = ns
            self.allocate_cl_buffers(E0=E_ref)

        print(f"OCLSplines::fit3D().2 Queue: {self.queue}, Context: {self.ctx}")

        ns = np.array( ns+(0,), dtype=np.int32 )
        #ns = np.array( ns, dtype=np.int32 )
        lsh = (4,4,4)


        cs_Err = np.array( [-1.0,1.0], dtype=np.float32 )
        cs_F   = np.array( [ 1.0,0.0], dtype=np.float32 )

        MDpar = [dt, dt, cdamp, 0]
        MDpar = np.array( MDpar, dtype=np.float32 )

        gsh  = roundup_global_size_3d( ns, lsh )
        ntot = np.int32( ns[0]*ns[1]*ns[2] )
        nL   = 32
        nG   = roundup_global_size( ntot, nL )
        nStepMax = nmaxiter // nPerStep

        print( "gsh ", gsh, " lsh ", lsh, " ns ", ns )

        out = np.zeros( ns_bak, dtype=np.float32)

        ConvTrj=None
        if bConvTrj:
            ConvTrj = np.zeros( (nStepMax,3) ) + np.nan
            
        print(f"OCLSplines::fit3D().3 Queue: {self.queue}, Context: {self.ctx}")

        for i in range(nStepMax):
            for j in range(nPerStep):
                self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns,   self.Gs_buf,  self.G0_buf,  self.dGs_buf, cs_Err )   # self.queue.finish() 
                self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns,   self.dGs_buf, None       ,  self.fGs_buf, cs_F   )   # self.queue.finish() 
                self.prg.move          (self.queue, (nG,), (nL,), ntot, self.Gs_buf,  self.vGs_buf, self.fGs_buf, MDpar  )   # self.queue.finish() 
                self.queue.finish() 
            cl.enqueue_copy(self.queue, out, self.fGs_buf); Ftot = np.max(np.abs(out))
            if bConvTrj:
                cl.enqueue_copy(self.queue, out, self.dGs_buf); Etot = np.max(np.abs(out))
                ConvTrj[i,:] = (0.0+i*nPerStep,Ftot,Etot)

            print( f"Ftot[{i*nPerStep}] {Ftot}" )
            if Ftot<Ftol:  break

        # cl.enqueue_copy(self.queue, out, self.fGs_buf); Ftot = out.norm()
        return out, ConvTrj
    
    def prepare_Morse_buffers(self, na, nxyz, bytePerFloat=4 ):
        print( "OCLSplines::prepare_Morse_buffers() " )
        #ntot=ng[0]*ng[1]*ng[2]
        self.atoms_buff  = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        self.REQs_buff   = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        self.V_Paul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat)
        self.V_Lond_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat)
        
    def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=(0.1, 0.1, 0.1), ng=None,            lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],                     grid_p0=(0.0, 0.0, 0.0), GFFParams=(0.1, 1.5, 0.0, 0.0)):

        T00 = time.perf_counter()

        na = len(atoms)

        if ng is None:
            ng = (int(lvec[0][0] / dg[0]), int(lvec[1][1] / dg[1]), int(lvec[2][2] / dg[2]))
        nxyz = ng[0]*ng[1]*ng[2]

        if self.atoms_buff is None:
            self.prepare_Morse_buffers(na, nxyz)

        atoms_np = np.asarray(atoms, dtype=np.float32)
        reqs_np = np.asarray(REQs, dtype=np.float32)

        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)
        cl.enqueue_copy(self.queue, self.REQs_buff, reqs_np)

        na_cl   = np.int32(na)
        nPBC_cl = np.array(nPBC+(0,),      dtype=np.int32   )
        ng_cl   = np.array(ng+(0,),        dtype=np.int32   )

        lvec_a = np.array( [*lvec[0],0.],   dtype=np.float32   )
        lvec_b = np.array( [*lvec[1],0.],   dtype=np.float32   )
        lvec_c = np.array( [*lvec[2],0.],   dtype=np.float32   )
        grid_p0_cl   = np.array(grid_p0+(0.0,),   dtype=np.float32 )
        GFFParams_cl = np.array(GFFParams, dtype=np.float32 )

        print( "lvec_a ", lvec_a )
        print( "lvec_b ", lvec_b )
        print( "lvec_c ", lvec_c )
        print( "grid_p0_cl ",   grid_p0_cl )
        print( "GFFParams_cl ", GFFParams_cl )
        ng_ = ng[::-1]; print("ng_ ", ng_)

        V_Paul = np.zeros( ng_, dtype=np.float32)
        V_Lond = np.zeros( ng_, dtype=np.float32)

        nL = self.nloc
        nG = roundup_global_size(nxyz, nL)

        #measure exact execution time
        
        T0 = time.perf_counter()

        self.prg.make_MorseFF(self.queue, (nG,), (nL,),
                            na_cl, self.atoms_buff, self.REQs_buff, self.V_Paul_buff, self.V_Lond_buff,
                            nPBC_cl, ng_cl, lvec_a,lvec_b,lvec_c, grid_p0_cl, GFFParams_cl)
    
        # Copy results back to host
        cl.enqueue_copy(self.queue, V_Paul, self.V_Paul_buff)
        cl.enqueue_copy(self.queue, V_Lond, self.V_Lond_buff)
        print( "make_MorseFF() time[s] ", time.perf_counter()-T0,   "   preparation  ", T0-T00 )
        # print( "V_Lond.min,max ", V_Lond.min(), V_Lond.max() )
        # print( "V_Paul.min,max ", V_Paul.min(), V_Paul.max() )
        return V_Paul, V_Lond



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
