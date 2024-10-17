import sys
import os
import numpy as np
import ctypes
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .clUtils import GridShape,GridCL

#from gpyfft.fft import FFT

COULOMB_CONST  =    14.3996448915 

bytePerFloat = 4

class GridFF_cl:

    def __init__(self, nloc=32 ):
        self.nloc  = nloc
        #self.ctx   = cl.create_some_context()
        #self.queue = cl.CommandQueue(self.ctx)
        self.grid  = None   # instance of GridShape, if initialized
        self.gcl   = None   # instance of GridCL, if initialized

        self.cl = cl
 
        # # Get a list of available OpenCL platforms
        # platforms = cl.get_platforms()
        # # Print information about each platform
        # print("Available OpenCL Platforms:")
        # for i, platform in enumerate(platforms):
        #     print(f"{i}. {platform.name}")
        # # Prompt the user to select a platform
        # platform_choice = int(input("Enter the number of the platform you want to use: "))
        # platform = platforms[platform_choice]
        #clu.get_cl_info( self.ctx.devices[0] )

        self.ctx,self.queue = clu.get_nvidia_device( what="nvidia")

        local_size = 64
        #print( " local_memory_per_workgroup() size=", local_size, " __local []  ", clu.local_memory_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " Byte " );
        print( " local_memory_per_workgroup() size=", local_size, " __local []  ", clu.local_memory_float_per_workgroup( self.ctx.devices[0], local_size=32, sp_per_cu=128 ), " float32 " );

        try:
            with open('../../cpp/common_resources/cl/GridFF.cl', 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print( "GridFF_cl() called from path=", os.getcwd() )
            print(f"Error compiling OpenCL program: {e}")
            exit(0)

        # self.atoms_buff    = None
        # self.REQs_buff     = None
        # self.V_Paul_buff   = None
        # self.V_Lond_buff   = None

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

        nG = clu.roundup_global_size( n, self.nloc )
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
        
        nG = clu.roundup_global_size( n, self.nloc )

        self.prg.sample1D_pbc(self.queue, (nG,), (self.nloc,),        g0, dg, ng, Gs_buf.data, np.int32(n),   ps_buf.data, fes_buf.data )
        
        return fes_buf.get()
    

    # ============ Fitting

    def allocate_cl_buffers(self, ns=None, E0=None, nByte=4):

        print(f"GridFF_cl::allocate_cl_buffers().1 Queue: {self.queue}, Context: {self.ctx}")

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

        print(f"GridFF_cl::allocate_cl_buffers().2 Queue: {self.queue}, Context: {self.ctx}")

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
    
    def fit3D(self, Ref_buff, nmaxiter=300, dt=0.5, Ftol=1e-16, damp=0.15, nPerStep=50, bConvTrj=False, bReturn=True, bPrint=False, bTime=True, bDebug=True ):
        # NOTE / ToDo : It is a bit strange than GridFF.h::makeGridFF_Bspline_d() the fit is fastes with damp=0.0 but here damp=0.15 performs better
        #print(f"GridFF_cl::fit3D().1 Queue: {self.queue}, Context: {self.ctx}")
        print(f"GridFF_cl::fit3D() dt={dt}, damp={damp} nmaxiter={nmaxiter} Ftol{Ftol}")
        T00=time.perf_counter()
        cdamp=1.-damp

        ns   = self.gsh.ns #[::-1] 
        nxyz = self.gcl.nxyz

        buff_names={'Gs','dGs','fGs','vGs'}
        self.try_make_buffs( buff_names, 0, nxyz )

        ns_cl = np.array( ns+(0,), dtype=np.int32 )
        #ns = np.array( ns, dtype=np.int32 )
        cs_Err   = np.array( [-1.0,1.0], dtype=np.float32 )
        cs_F     = np.array( [ 1.0,0.0], dtype=np.float32 )
        MDpar_cl = np.array( [dt, dt, cdamp, 0], dtype=np.float32 )

        lsh = (4,4,4)
        gsh  = clu.roundup_global_size_3d( ns, lsh )

        nL   = 32
        nG   = clu.roundup_global_size( nxyz, nL )
        nStepMax = nmaxiter // nPerStep
        #print( "GridFF_cl::::fit3D() gsh ", gsh, " lsh ", lsh, " ns ", ns )
        out = np.zeros( ns[::-1], dtype=np.float32)

        ConvTrj=None
        if bConvTrj:
            ConvTrj = np.zeros( (nStepMax,3) ) + np.nan
            
        #print(f"GridFF_cl::fit3D() Queue: {self.queue}, Context: {self.ctx}")
        """
        __kernel void setMul(
            const int  ntot,
            __global float* v,
            __global float* out,  
            float c
        ) {
        """

        T0=time.perf_counter()

        # if bDebug:
        #     cl.enqueue_copy(self.queue, out, Ref_buff);
        #     print( "GridFF_Cl::fit3D() debug Ref_buff.min,max: ", out.min(), out.max() )

        self.prg.setMul(self.queue, (nG,), (nL,), nxyz,  Ref_buff,  self.Gs_buff,  np.float32(1.0) ) # setup force
        self.prg.setMul(self.queue, (nG,), (nL,), nxyz,  Ref_buff,  self.vGs_buff, np.float32(0.0) ) # setup velocity
        
        # if bDebug:
        #     cl.enqueue_copy(self.queue, self.Gs_buff, Ref_buff);
        #     print( "GridFF_Cl::fit3D() debug Gs_buff.min,max: ", out.min(), out.max() )

        # --- debug run
        # for i in range(1):
        #     plt.figure(); cl.enqueue_copy(self.queue, out, self.Gs_buff); self.queue.finish(); plt.imshow( out[1,:,:] ); plt.title("Gs_buff"); #plt.show()
        #     self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns_cl, self.Gs_buff,  Ref_buff,      self.dGs_buff, cs_Err    )  
        #     plt.figure(); cl.enqueue_copy(self.queue, out, self.dGs_buff); self.queue.finish(); plt.imshow( out[1,:,:] ); plt.title("dGs_buff"); #plt.show()
        #     self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns_cl, self.dGs_buff, None       ,   self.fGs_buff, cs_F      ) 
        #     plt.figure(); cl.enqueue_copy(self.queue, out, self.fGs_buff); self.queue.finish(); plt.imshow( out[1,:,:] ); plt.title("fGs_buff"); 
        #     plt.show()
        #     exit()

        nstepdone=0
        for i in range(nStepMax):
            for j in range(nPerStep):
                self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns_cl, self.Gs_buff,  Ref_buff,      self.dGs_buff, cs_Err    )   # self.queue.finish() 
                self.prg.BsplineConv3D (self.queue, gsh,   lsh,   ns_cl, self.dGs_buff, None       ,   self.fGs_buff, cs_F      )   # self.queue.finish() 
                self.prg.move          (self.queue, (nG,), (nL,), nxyz,  self.Gs_buff,  self.vGs_buff, self.fGs_buff, MDpar_cl  )   # self.queue.finish() 
                nstepdone+=1
            cl.enqueue_copy(self.queue, out, self.fGs_buff); self.queue.finish(); Ftot = np.max(np.abs(out))
            if bConvTrj:
                cl.enqueue_copy(self.queue, out, self.dGs_buff); self.queue.finish(); Etot = np.max(np.abs(out))
                ConvTrj[i,:] = (0.0+i*nPerStep,Ftot,Etot)
                if bPrint: print( f"GridFF::fit3D()[{i*nPerStep}] |F|={Ftot} |E|={Etot}" )
            else:
                if bPrint: print( f"GridFF::fit3D()[{i*nPerStep}] |F|={Ftot} " )
            if Ftot<Ftol:  break

        if bTime:
            self.queue.finish()
            dT=time.perf_counter()-T0
            nops = nstepdone*nxyz
            self.queue.finish()
            print( "GridFF_cl::fit3D() time[s] ", dT, "  preparation[s]  ", T0-T00, "[s] nGOPs: ", nops*1e-9," speed[GOPs/s]: ", (nops*1e-9)/dT , " nstep, nxyz ", nstepdone, nxyz  )

        if bReturn:
            cl.enqueue_copy(self.queue, out, self.Gs_buff);
            # finish opencl
            self.queue.finish()
            if bDebug:
                print( "GridFF_cl::fit3D() DONE Gs_buff.min,max: ", out.min(), out.max() )
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


    def try_make_buffs(self, names, na, nxyz, bytePerFloat=4):
        self.try_buff("atoms",  names, na*4*bytePerFloat)
        self.try_buff("REQs",   names, na*4*bytePerFloat)
        self.try_buff("V_Paul", names, nxyz*bytePerFloat)
        self.try_buff("V_Lond", names, nxyz*bytePerFloat)
        self.try_buff("Qgrid",  names, nxyz*2*bytePerFloat)
        self.try_buff("Vgrid",  names, nxyz*2*bytePerFloat)
        self.try_buff("V_Coul", names, nxyz*bytePerFloat)

        # fitting buffers      
        self.try_buff("Gs",  names, nxyz*bytePerFloat )
        self.try_buff("fGs", names, nxyz*bytePerFloat )
        self.try_buff("dGs", names, nxyz*bytePerFloat )
        self.try_buff("vGs", names, nxyz*bytePerFloat )

        # New buffers for laplace_real_loop_inert
        #self.try_buff("V1", names, nxyz*bytePerFloat)
        #self.try_buff("V2", names, nxyz*bytePerFloat)
        #self.try_buff("vV", names, nxyz*bytePerFloat)
        
        # New buffers for laplace_real_pbc and slabPotential
        # self.try_buff("Vin", names, nxyz*bytePerFloat)
        # self.try_buff("Vout", names, nxyz*bytePerFloat)


    def prepare_Morse_buffers(self, na, nxyz, bytePerFloat=4 ):
        #print( "GridFF_cl::prepare_Morse_buffers() " )
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

    def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=(0.1, 0.1, 0.1), ng=None,            lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],                     g0=(0.0, 0.0, 0.0), GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True ):

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
        nG = clu.roundup_global_size(nxyz, nL)

        T0 = time.perf_counter()
        self.prg.make_MorseFF(self.queue, (nG,), (nL,),
                            na_cl, self.atoms_buff, self.REQs_buff, self.V_Paul_buff, self.V_Lond_buff,
                            nPBC_cl, self.gcl.ns, self.gcl.a,self.gcl.b,self.gcl.c, self.gcl.g0, GFFParams_cl)
    
        if bTime:
            self.queue.finish()
            dT=time.perf_counter()-T0
            npbc = (nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
            nops = na_cl * nxyz * npbc
            print( "GridFF_cl::make_MorseFF() time[s] ", dT,   "  preparation[s]  ", T0-T00, "[s] nGOPs: ", nops*1e-9," speed[GOPs/s]: ", (nops*1e-9)/dT , " na,nxyz,npbc ", na_cl, nxyz, npbc  )

        if bReturn:
            sh = self.gsh.ns[::-1]
            V_Paul = np.zeros( sh, dtype=np.float32)
            V_Lond = np.zeros( sh, dtype=np.float32)
            cl.enqueue_copy(self.queue, V_Paul, self.V_Paul_buff)
            cl.enqueue_copy(self.queue, V_Lond, self.V_Lond_buff)
            # print( "V_Lond.min,max ", V_Lond.min(), V_Lond.max() )
            # print( "V_Paul.min,max ", V_Paul.min(), V_Paul.max() )
            return V_Paul, V_Lond
    

    def _project_atoms_on_grid_quintic_pbc(self, sz_glob, sz_loc, na, ns ):
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
        nxyz2 = np.int32( ns[0]*ns[1]*ns[2] * 2 )
        self.prg.set( self.queue, sz_glob, sz_loc, nxyz2, self.Qgrid_buff, np.float32(0.0) )
        self.prg.project_atoms_on_grid_quintic_pbc( self.queue, sz_glob, sz_loc,
            na, self.atoms_buff, self.Qgrid_buff,
            ns, self.gcl.g0, self.gcl.dg
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
        nG = clu.roundup_global_size( self.gcl.nxyz, nL)
        self._project_atoms_on_grid_quintic_pbc( (nG,), (nL,), np.int32(na), self.gcl.ns )
        
        if bReturn:
            sh = self.gsh.ns[::-1]+(2,)
            Qgrid = np.zeros( sh, dtype=np.float32)
            #print( "Qgrid.shape", Qgrid.shape )
            cl.enqueue_copy(self.queue, Qgrid, self.Qgrid_buff)
            Qgrid = Qgrid[:,:,:,0].copy()   # take just the real part
            #print( Qgrid[:,0,0], Qgrid[0,:,0], Qgrid[0,0,:], )
            print( "GridFF_cl::project_atoms_on_grid_quintic_pbc() DONE Qtot=", Qgrid.sum()," Qabs=", np.abs(Qgrid).sum()," Qmin=", Qgrid.min()," Qmax=", Qgrid.max() )
            #self.Qgrid=Qgrid
            return Qgrid

    def _poisson(self, sz_glob, sz_loc, ns, coefs_cl ):
        event, = self.transform.enqueue() 
        event.wait()
        self.prg.poissonW( self.queue, sz_glob, sz_loc,       ns, self.Qgrid_buff, self.Vgrid_buff, coefs_cl )
        event, = self.inverse_transform.enqueue()
        event.wait()

    def prepare_poisson(self, sh=None):
        if sh is None: sh=self.gsh.ns[::-1] 
        print("GridFF_cl::prepare_poisson() sh=", sh )
        self.Qgrid_cla = cl_array.Array(self.queue, shape=sh, dtype=np.complex64, data=self.Qgrid_buff )
        self.Vgrid_cla = cl_array.Array(self.queue, shape=sh, dtype=np.complex64, data=self.Vgrid_buff )
        self.transform         = clu.FFT(self.ctx, self.queue, self.Qgrid_cla, axes=(0, 1, 2))
        self.inverse_transform = clu.FFT(self.ctx, self.queue, self.Vgrid_cla, axes=(0, 1, 2))

    def poisson(self, bReturn=True, sh=None, dV=None ):

        clu.try_load_clFFT()
        if sh is None: sh=self.gsh.ns[::-1] 
        nxyz = np.int32( sh[0]*sh[1]*sh[2] )

        if dV is None: dV = self.gsh.dV

        buff_names={'Qgrid','Vgrid'}
        self.try_make_buffs(buff_names, 0, nxyz)
        self.prepare_poisson( sh=sh )

        ns_cl    = np.array( (*sh[::-1],nxyz), dtype=np.int32 )
        sc_ewald = 4*np.pi*dV*COULOMB_CONST
        coefs_cl = np.array( (0.0,0.0,0.0, sc_ewald ), dtype=np.float32 )

        nL = self.nloc
        nG = clu.roundup_global_size( nxyz, nL)
        self._poisson( (nG,), (nL,), ns_cl, coefs_cl )        

        if bReturn:
            sh_ = (*sh,2)
            print( "sh ", sh_ )
            Vgrid = np.zeros( sh_, dtype=np.float32   )  
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
        
        buff_names={'Qgrid','Vgrid'}
        self.try_make_buffs(buff_names, 0, nxyz)
                
        coefs    =(0.0,0.0,0.0,4*np.pi*dV)
        coefs_cl = np.array( coefs,        dtype=np.float32 )
        ns_cl    = np.array( ng+(nxyz,),   dtype=np.int32   )
        
        self.prg.poissonW( self.queue, sz_glob, sz_loc,
            ns_cl, self.Qgrid_buff, self.Vgrid_buff,
            coefs_cl
        )
        
        out = np.zeros(N, dtype=np.complex64)
        cl.enqueue_copy(self.queue, out, self.out_buff)
    
        return out

    def slabPotential(self, nz_slab, dz, Vol, dVcor, Vcor0, bDownload=True):
        buff_names = {'Vgrid', 'V_coul'}
        self.try_make_buffs(buff_names, 0, self.gcl.nxyz )
        ns_cl  = np.array( self.gsh.ns+(nz_slab,),  dtype=np.int32   )
        params = np.array( [dz, Vol, dVcor, Vcor0], dtype=np.float32 )
        sz_loc  = (4,4,4,)
        sz_glob = clu.roundup_global_size_3d( self.gsh.ns, sz_loc)
        self.prg.slabPotential( self.queue, sz_glob, sz_loc,
            ns_cl, self.Vgrid_buff, self.V_coul_buff, params )
        if bDownload:
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
        sz_glob = clu.roundup_global_size_3d( self.gsh.ns, sz_loc )
        self._laplace_real_pbc( sz_loc, sz_glob, cSOR, np.float32(cSOR) )
       
        Vout = np.empty_like(Vin)
        cl.enqueue_copy(self.queue, Vout, self.Vout_buff)
        return Vout


    def laplace_real_loop_inert(self, nmaxiter=1000, tol=1e-6, bPBC=True, cSOR=0.0, cV=0.5, bReturn=False, sh=None ):
        print( "GridFF_cl::laplace_real_loop_inert() " )

        if sh is None: sh=self.gsh.ns[::-1]
        nxyz = np.int32( sh[0]*sh[1]*sh[2] )

        buff_names = {'V1', 'V2', 'vV'}
        self.try_make_buffs(buff_names, 0, nxyz )

        sz_loc = (4,4,4)
        sz_glob = clu.roundup_global_size_3d( self.gsh.ns, sz_loc )

        ns_cl   = np.array( (*sh,0), dtype=np.int32 )
        C_cl    = np.array( (1.0,0.0), dtype=np.float32 )
        cV_cl   = np.float32(cV)
        cSOR_cl = np.float32(cSOR)

        szl = (self.nloc,)
        szg = ( clu.roundup_global_size( self.gcl.nxyz, self.nloc), )
        self.prg.setCMul(self.queue, szg, szl, nxyz,  self.Vgrid_buff,  self.V1_buff, C_cl ) # copy real part from complex Vgrid to scalar V1

        self.prg.laplace_real_pbc(     self.queue, sz_glob, sz_loc,       ns_cl, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl, np.float32(0.0) )
        for iter in range(nmaxiter):
            if iter % 2 == 0:
                self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,    ns_cl, self.V2_buff, self.V1_buff, self.vV_buff, cSOR_cl,cV_cl )
            else:
                self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,    ns_cl, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl,cV_cl )
        #if iter % 2 == 1:
        #    cl.enqueue_copy(self.queue, self.V__buff, self.V_buff)    self.V_buff, self.V__buff = self.V__buff, self.V_buff
        if bReturn:
            result = np.empty_like(V)
            cl.enqueue_copy(self.queue, result, self.V_buff)
        return result, iter
    
    def makeCoulombEwald(self, atoms ):
        print( "GridFF_cl::makeCoulombEwald() " )
        clu.try_load_clFFT()
        if self.gcl is None: 
            print("ERROR in GridFF_cl::makeCoulombEwald() gcl is None, => please call set_grid() first " )
            exit()

        na = len(atoms)
        buff_names={'atoms','Qgrid','Vgrid'}

        sz_loc  = (self.nloc,)
        sz_glob = ( clu.roundup_global_size( self.gcl.nxyz, self.nloc), )

        self.try_make_buffs(buff_names, na, self.gcl.nxyz )
        atoms_np = np.array(atoms, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)

        self._project_atoms_on_grid_quintic_pbc( sz_glob, sz_loc, np.int32(na), self.gcl.ns )

        #self._poisson( sz_glob, sz_loc, )
        # Vcoul = np.zeros( self.gsh.ns+(2,), dtype=np.float32   )  
        # cl.enqueue_copy(self.queue, Vcoul, self.Vcoul_buff )
        # return Vcoul
    
        Vgrid = self.poisson()
        Vgrid = Vgrid[:,:,:,0].copy()

        #print( "Vgrid min,max", Vgrid.min(), Vgrid.max() )
        return Vgrid


    def makeCoulombEwald_slab(self, atoms, Lz_slab=20.0, bDipoleCoorection=False, bReturn=True ):
        print( "GridFF_cl::makeCoulombEwald() " )
        clu.try_load_clFFT()
        if self.gcl is None: 
            print("ERROR in GridFF_cl::makeCoulombEwald() gcl is None, => please call set_grid() first " )
            exit()

        nz_slab   = Lz_slab/self.gsh.dg[2]; print("nz_slab", nz_slab, " ns[2] ", self.gcl.ns[2] )
        ns_cl     = np.array( (self.gcl.ns[0],self.gcl.ns[1],self.gcl.ns[2]+nz_slab,0), dtype=np.int32 ); print("ns_cl", ns_cl)
        nxyz_slab = ns_cl[0]*ns_cl[1]*ns_cl[2]

        #if sh is None: sh=self.gsh.ns[::-1] 
        sh=ns_cl[:3][::-1] 

        na = len(atoms)
        buff_names={'atoms','Qgrid','Vgrid'}

        sz_loc  = (self.nloc,)
        sz_glob = ( clu.roundup_global_size( nxyz_slab, self.nloc), )

        self.try_make_buffs(buff_names, na, nxyz_slab )
        atoms_np = np.array(atoms, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)

        print("GridFF_cl::makeCoulombEwald_slab()._project_atoms_on_grid_quintic_pbc")
        self._project_atoms_on_grid_quintic_pbc( sz_glob, sz_loc, np.int32(na), ns_cl  )    

        # if bReturn:
        #     sh_ = (*sh,2)
        #     Qgrid = np.zeros( sh_, dtype=np.float32)
        #     #print( "Qgrid.shape", Qgrid.shape )
        #     cl.enqueue_copy(self.queue, Qgrid, self.Qgrid_buff)
        #     Qgrid = Qgrid[:,:,:,0].copy()   # take just the real part
        #     #print( Qgrid[:,0,0], Qgrid[0,:,0], Qgrid[0,0,:], )
        #     print( "GridFF_cl::makeCoulombEwald_slab().project DONE Qgrid.shape=", Qgrid.shape,"  Qtot= ", Qgrid.sum()," Qabs=", np.abs(Qgrid).sum()," Qmin=", Qgrid.min()," Qmax=", Qgrid.max() )
        #     #self.Qgrid=Qgrid
        #     return Qgrid
        
        if bDipoleCoorection:
            self.poisson( bReturn=False, sh=sh )
            #Vgrid = Vgrid[:,:,:,0].copy()
            #print( "Vgrid min,max", Vgrid.min(), Vgrid.max() )
            # ToDo: We should distinguis between z-dimension of the slab and of the extended grid !!!
            dipole_z = np.sum( atoms[:,2]*atoms[:,3] )
            Vol      = self.gsh.getVolume()
            dVcor    = 4.0 *np.pi * COULOMB_CONST * dipole_z/Vol;
            Vcor0    = -dVcor * self.gsh.Ls[2]/2;
            Vcoul = self.slabPotential( nz_slab, dz, Vol, dVcor, Vcor0, bDownload=True)
        else:
            print("GridFF_cl::makeCoulombEwald_slab().poisson")
            Vgrid = self.poisson( bReturn=True, sh=sh )
            Vcoul = Vgrid[:,:,:,0].copy()
            print("GridFF_cl::makeCoulombEwald_slab() DONE Vgrid.shape ", Vgrid.shape, " Vcoul.shape ", Vcoul.shape )

        return Vcoul


