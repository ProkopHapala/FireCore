import sys
import os
import numpy as np
import ctypes
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

#from gpyfft.fft import FFT

FFT = None

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

def try_load_clFFT():
    global FFT
    if FFT is None:
        from gpyfft.fft import FFT as FFT_
        FFT = FFT_

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

class GridShape:
    def __init__(self, ns=None, dg=(0.0,0.0,0.0), lvec=None, Ls=None, g0=(0.0,0.0,0.0) ):
        if lvec is None: lvec = np.array( [ [Ls[0],0.,0.], [0.,Ls[1],0.], [0.,0.,Ls[2]] ] )
        if Ls   is None: Ls = (lvec[0][0],lvec[1][1],lvec[2][2])
        if ns   is None: ns = (int((Ls[0]/dg[0])+0.5),int((Ls[1]/dg[1])+0.5),int((Ls[2]/dg[2])+0.5))
        self.ns   = ns
        self.Ls   = Ls
        self.g0 = g0
        self.dg   = dg
        self.lvec = lvec
        self.dV   = dg[0]*dg[1]*dg[2] 

    def __str__(self):
        return f"GridShape(ns={self.ns}, Ls={self.Ls}, pos0={self.pos0}, dg={self.dg}, lvec={self.lvec})"

class GridCL:
    def __init__(self, gsh : GridShape ):
        self.nxyz = np.int32( gsh.ns[0]*gsh.ns[1]*gsh.ns[2] )
        self.ns = np.array( gsh.ns+(self.nxyz,), dtype=np.int32   )
        self.g0 = np.array( gsh.g0+(0.,),   dtype=np.float32 )
        self.dg = np.array( gsh.dg+(0.,),   dtype=np.float32 )
        self.a  = np.array( ( *gsh.lvec[0], 0.), dtype=np.float32 )
        self.b  = np.array( ( *gsh.lvec[1], 0.), dtype=np.float32 )
        self.c  = np.array( ( *gsh.lvec[2], 0.), dtype=np.float32 )
        self.da = np.array( ( *(gsh.lvec[0]/gsh.ns[0]), 0.), dtype=np.float32 )
        self.db = np.array( ( *(gsh.lvec[1]/gsh.ns[1]), 0.), dtype=np.float32 )
        self.dc = np.array( ( *(gsh.lvec[2]/gsh.ns[2]), 0.), dtype=np.float32 )
        
    def __str__(self):
        return f"GridShape(ng={self.ng}, g0={self.g0}, dg={self.dg}\n   La={self.La},Lb={self.Lb},Lc={self.Lc}\n    da={self.da},db={self.db},dc={self.dc})"


class OCLSplines:

    flag_atoms  = 0b000001  # 1
    flag_REQs   = 0b000010  # 2
    flag_qs     = 0b000100  # 4
    flag_VMorse = 0b001000  # 8
    flag_VCoul  = 0b010000  # 16
    flag_Qgrid  = 0b100000  # 32

    def __init__(self, nloc=32 ):
        self.nloc  = nloc
        self.ctx   = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        self.grid  = None   # instance of GridShape, if initialized
        self.gcl   = None   # instance of GridCL, if initialized
 
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

    def set_grid(self, gsh : GridShape ):
        self.gsh = gsh
        self.gcl = GridCL(gsh)

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
    
    def try_buff(self, name, names, sz ):
        if name in names:
            buff_name = name+"_buff"
            if hasattr(self,buff_name): 
                buff = getattr(self, buff_name)
                if not ( buff is None or buff.size != sz ):
                    return
            print( "try_buff(",buff_name,") reallocate to [bytes]: ", sz )
            buff = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, sz )
            setattr(self, buff_name, buff )

    # def try_make_buffs(self, names, na, nxyz, bytePerFloat=4):
    #     self.try_buff("atoms",  names, na*4*bytePerFloat)
    #     self.try_buff("REQs",   names, na*4*bytePerFloat)
    #     self.try_buff("V_Paul", names, nxyz*bytePerFloat)
    #     self.try_buff("V_Lond", names, nxyz*bytePerFloat)
    #     self.try_buff("Qgrid",  names, nxyz*2*bytePerFloat)
    #     self.try_buff("Vgrid",  names, nxyz*2*bytePerFloat)
    #     self.try_buff("V_Coul", names, nxyz*bytePerFloat)
    
    def try_make_buffs(self, names, na, nxyz, bytePerFloat=4):
        self.try_buff("atoms",  names, na*4*bytePerFloat)
        self.try_buff("REQs",   names, na*4*bytePerFloat)
        self.try_buff("V_Paul", names, nxyz*bytePerFloat)
        self.try_buff("V_Lond", names, nxyz*bytePerFloat)
        self.try_buff("Qgrid",  names, nxyz*2*bytePerFloat)
        self.try_buff("Vgrid",  names, nxyz*2*bytePerFloat)
        self.try_buff("V_Coul", names, nxyz*bytePerFloat)
        
        # New buffers for laplace_real_loop_inert
        self.try_buff("V1", names, nxyz*bytePerFloat)
        self.try_buff("V2", names, nxyz*bytePerFloat)
        self.try_buff("vV", names, nxyz*bytePerFloat)
        
        # New buffers for laplace_real_pbc and slabPotential
        # self.try_buff("Vin", names, nxyz*bytePerFloat)
        # self.try_buff("Vout", names, nxyz*bytePerFloat)


    def prepare_Morse_buffers(self, na, nxyz, bytePerFloat=4 ):
        #print( "OCLSplines::prepare_Morse_buffers() " )
        #ntot=ng[0]*ng[1]*ng[2]
        # -- used_by: make_MorseFF, project_atoms_on_grid
        self.atoms_buff  = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        # -- used_by: make_MorseFF
        self.REQs_buff   = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        self.V_Paul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat)
        self.V_Lond_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat)
        # -- used_by: project_atoms_on_grid, 
        self.Qgrid_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )
        # -- used_by: project_atoms_on_grid, poisson
        self.V_Coul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )

    def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=(0.1, 0.1, 0.1), ng=None,            lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],                     g0=(0.0, 0.0, 0.0), GFFParams=(0.1, 1.5, 0.0, 0.0), bReturn=True ):

        T00 = time.perf_counter()

        grid = GridShape( ns=ng, dg=dg, lvec=lvec, g0=g0 )
        self.set_grid( grid )

        na   = len(atoms)
        nxyz = self.gcl.nxyz

        buff_names={'atoms','REQs','V_Paul','V_Lond'}
        self.try_make_buffs(buff_names, na, nxyz)

        atoms_np = np.asarray(atoms, dtype=np.float32)
        reqs_np  = np.asarray(REQs, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)
        cl.enqueue_copy(self.queue, self.REQs_buff, reqs_np)

        na_cl        = np.int32(na)
        nPBC_cl      = np.array(nPBC+(0,), dtype=np.int32   )
        GFFParams_cl = np.array(GFFParams, dtype=np.float32 )

        nL = self.nloc
        nG = roundup_global_size(nxyz, nL)

        T0 = time.perf_counter()
        self.prg.make_MorseFF(self.queue, (nG,), (nL,),
                            na_cl, self.atoms_buff, self.REQs_buff, self.V_Paul_buff, self.V_Lond_buff,
                            nPBC_cl, self.gcl.ns, self.gcl.a,self.gcl.b,self.gcl.c, self.gcl.g0, GFFParams_cl)
    
        if bReturn:
            npbc = (nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
            sh = self.gsh.ns[::-1]
            V_Paul = np.zeros( sh, dtype=np.float32)
            V_Lond = np.zeros( sh, dtype=np.float32)
            cl.enqueue_copy(self.queue, V_Paul, self.V_Paul_buff)
            cl.enqueue_copy(self.queue, V_Lond, self.V_Lond_buff)
            nops = na_cl * nxyz * npbc
            dT=time.perf_counter()-T0
            print( "make_MorseFF() time[s] ", dT,   "  preparation[s]  ", T0-T00, "[s] nGOPs: ", nops*1e-9," speed[GOPs/s]: ", (nops*1e-9)/dT , " V_Paul.nbytes: ", V_Paul.nbytes*1e-6, "[Mbyte] na,nxyz,npbc ", na_cl, nxyz, npbc  )
            # print( "V_Lond.min,max ", V_Lond.min(), V_Lond.max() )
            # print( "V_Paul.min,max ", V_Paul.min(), V_Paul.max() )
            return V_Paul, V_Lond
    

    def _project_atoms_on_grid_quintic_pbc(self, sz_glob, sz_loc, na ):
        """
        __kernel void project_atoms_on_grid_quintic_pbc(
            const int num_atoms,            // 1 number of atoms
            __global const float4* atoms,   // 2 Atom positions and charges
            __global       float*  Qgrid,   // 3 Output grid
            const int4   ng,                // 4 Grid size
            const float4 g0,                // 5 Grid origin
            const float4 dg                 // 6 Grid dimensions
        ) {
        """
        self.prg.project_atoms_on_grid_quintic_pbc( self.queue, sz_glob, sz_loc,
            na, self.atoms_buff, self.Qgrid_buff,
            self.gcl.ns,  self.gcl.g0, self.gcl.dg
        )

    def project_atoms_on_grid_quintic_pbc(self, atoms, ng=None,   dg=(0.1, 0.1, 0.1),    lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]], g0=(0.0, 0.0, 0.0), bReturn=True ):

        grid = GridShape( ns=ng, dg=dg, lvec=lvec, g0=g0 )
        self.set_grid( grid )

        na = len(atoms)
        buff_names={'atoms','Qgrid'}
        self.try_make_buffs(buff_names, na, self.gcl.nxyz )

        atoms_np = np.array(atoms, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)
        
        nL = self.nloc
        nG = roundup_global_size( self.gcl.nxyz, nL)
        self._project_atoms_on_grid_quintic_pbc( (nG,), (nL,), np.int32(na) )
        
        if bReturn:
            sh = self.gsh.ns[::-1]+(2,)
            Qgrid = np.zeros( sh, dtype=np.float32)
            print( "Qgrid.shape", Qgrid.shape )
            cl.enqueue_copy(self.queue, Qgrid, self.Qgrid_buff)
            Qgrid = Qgrid[:,:,:,0].copy()   # take just the real part
            #print( Qgrid[:,0,0], Qgrid[0,:,0], Qgrid[0,0,:], )
            print( "project_atoms_on_grid_quintic_pbc() DONE Qtot=", Qgrid.sum()," Qabs=", np.abs(Qgrid).sum()," Qmin=", Qgrid.min()," Qmax=", Qgrid.max() )
            #self.Qgrid=Qgrid
            return Qgrid

    def _poisson(self, sz_glob, sz_loc ):
        event, = self.transform.enqueue()
        event.wait()
        coefs_cl = np.array( (0.0,0.0,0.0,4*np.pi*self.gsh.dV ), dtype=np.float32 )
        self.prg.poissonW( self.queue, sz_glob, sz_loc,
            self.gcl.ns, self.Qgrid_buff, self.Vgrid_buff, coefs_cl )
        event, = self.inverse_transform.enqueue()
        event.wait()

    def prepare_poisson(self):
        sh=self.gsh.ns[::-1] 
        self.Qgrid_cla = cl_array.Array(self.queue, shape=sh, dtype=np.complex64, data=self.Qgrid_buff )
        self.Vgrid_cla = cl_array.Array(self.queue, shape=sh, dtype=np.complex64, data=self.Vgrid_buff )
        self.transform         = FFT(self.ctx, self.queue, self.Qgrid_cla, axes=(0, 1, 2))
        self.inverse_transform = FFT(self.ctx, self.queue, self.Vgrid_cla, axes=(0, 1, 2))

    def poisson(self, bReturn=True):

        try_load_clFFT()

        buff_names={'Qgrid','Vgrid'}
        self.try_make_buffs(buff_names, 0, self.gcl.nxyz)
        self.prepare_poisson()

        nL = self.nloc
        nG = roundup_global_size( self.gcl.nxyz, nL)
        self._poisson( (nG,), (nL,) )        

        if bReturn:
            sh = self.gsh.ns[::-1]+(2,)
            print( "sh ", sh )
            Vgrid = np.zeros( sh, dtype=np.float32   )  
            cl.enqueue_copy ( self.queue, Vgrid, self.Vgrid_buff         )
            return Vgrid

    def poissonW(self, sz_glob, sz_loc, ng , dV=0.001 ):
        """
        __kernel void poissonW(
            const int4   ns,
            __global float2* A,
            __global float2* out,
            const float4 dCell     // 
        ){ 
        """
        nxyz = self.gcl.nxyz
        
        buff_names={'Qgrid','V_Coul'}
        self.try_make_buffs(buff_names, 0, nxyz)
                
        coefs    =(0.0,0.0,0.0,4*np.pi*dV)
        coefs_cl = np.array(coefs, dtype=np.float32)
        ns_cl    = np.array( ng+(nxyz,),               dtype=np.int32   )
        
        self.prg.poissonW( self.queue, sz_glob, sz_loc,
            ns_cl, self.Qgrid_buff, self.V_Coul_buff,
            coefs_cl
        )
        
        out = np.zeros(N, dtype=np.complex64)
        cl.enqueue_copy(self.queue, out, self.out_buff)
    
        return out

    def slabPotential(self, Vin, nz_slab, dz, Vol, dVcor, Vcor0):
        nxyz = Vin.shape[0] * Vin.shape[1] * Vin.shape[2]
        buff_names = {'Vin', 'Vout'}
        self.try_make_buffs(buff_names, 0, nxyz)

        ng = np.array(Vin.shape + (0,), dtype=np.int32)
        params = np.array([dz, Vol, dVcor, Vcor0], dtype=np.float32)
        
        cl.enqueue_copy(self.queue, self.Vin_buff, Vin.astype(np.float32))
        
        nL = self.nloc
        nG = roundup_global_size( self.gcl.nxyz, nL)
        self.prg.slabPotential(   self.queue, (nG,), (nL,),
            np.int32(nz_slab), self.Vin_buff, self.Vout_buff, params, ng
        )

        Vout = np.empty_like(Vin)
        cl.enqueue_copy(self.queue, Vout, self.Vout_buff)
        return Vout
    
    def _laplace_real_pbc(self, sz_glob, sz_loc, cSOR):

        self.prg.laplace_real_pbc(  self.queue, sz_glob, sz_loc,
            self.gcl.ng, self.Vin_buff, self.Vout_buff, cSOR
        )
    
    def laplace_real_pbc(self, Vin, cSOR=0.0):
        nxyz = Vin.shape[0] * Vin.shape[1] * Vin.shape[2]
        buff_names = {'Vin', 'Vout'}
        self.try_make_buffs(buff_names, 0, nxyz)        
        cl.enqueue_copy(self.queue, self.Vin_buff, Vin.astype(np.float32))
        
        sz_loc = (4,4,4)
        sz_glob = roundup_global_size_3d( self.gsh.ns, sz_loc )
        self._laplace_real_pbc( sz_loc, sz_glob, cSOR, np.float32(cSOR) )
       
        Vout = np.empty_like(Vin)
        cl.enqueue_copy(self.queue, Vout, self.Vout_buff)
        return Vout


    def laplace_real_loop_inert(self, V, nmaxiter=1000, tol=1e-6, bPBC=True, cSOR=0.0, cV=0.5, ):
        print( "OCLSplines::laplace_real_loop_inert() " )
        ntot = V.size
        ng = np.array(V.shape + (0,), dtype=np.int32)

        buff_names = {'V1', 'V2', 'vV'}
        self.try_make_buffs(buff_names, 0, ntot)

        cl.enqueue_copy(self.queue, self.V_buff, V.astype(np.float32))

        sz_loc = (4,4,4)
        sz_glob = roundup_global_size_3d( self.gsh.ns, sz_loc )

        cV_cl = np.float32(cV)
        cSOR_cl = np.float32(cSOR)

        self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,
            ng, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl,cV_cl)
        for iter in range(nmaxiter):
            self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,
                ng, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl,cV_cl )

        if iter % 2 == 1:
            cl.enqueue_copy(self.queue, self.V__buff, self.V_buff)
            self.V_buff, self.V__buff = self.V__buff, self.V_buff

        result = np.empty_like(V)
        cl.enqueue_copy(self.queue, result, self.V_buff)
        return result, iter
    
    def makeCoulombEwald(self, atoms ):
        print( "OCLSplines::makeCoulombEwald() " )
        try_load_clFFT()
        if self.gcl is None: 
            print("ERROR in OCLSplines::makeCoulombEwald() gcl is None, => please call set_grid() first " )
            exit()

        na = len(atoms)
        buff_names={'atoms','Qgrid','Vgrid'}

        sz_loc  = (self.nloc,)
        sz_glob = ( roundup_global_size( self.gcl.nxyz, self.nloc), )

        self.try_make_buffs(buff_names, na, self.gcl.nxyz )
        atoms_np = np.array(atoms, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)

        self._project_atoms_on_grid_quintic_pbc( sz_glob, sz_loc, np.int32(na) )

        #self._poisson( sz_glob, sz_loc, )
        # Vcoul = np.zeros( self.gsh.ns+(2,), dtype=np.float32   )  
        # cl.enqueue_copy(self.queue, Vcoul, self.Vcoul_buff )
        # return Vcoul
    
        Vgrid = self.poisson()
        Vgrid = Vgrid[:,:,:,0].copy()

        print( "Vgrid min,max", Vgrid.min(), Vgrid.max() )
        return Vgrid



