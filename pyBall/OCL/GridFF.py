import sys
import os
import numpy as np
import ctypes
import pyopencl as cl

import time
import json
import os
import matplotlib.pyplot as plt

import pyopencl.tools
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .clUtils import GridShape,GridCL

# Dictionary to store kernel execution times
kernel_times = {}

def profile_kernel(name, queue, kernel, global_size, local_size, *args, **kwargs):
    """Profile a kernel execution"""
    # Record start time
    start_time = time.time()
    
    # Call the kernel
    event = kernel(queue, global_size, local_size, *args, **kwargs)
    
    # Wait for the kernel to complete
    event.wait()
    
    # Record end time
    end_time = time.time()
    
    # Calculate duration
    duration = (end_time - start_time) * 1000  # Convert to ms
    
    # Add to kernel times
    if name not in kernel_times:
        kernel_times[name] = []
    kernel_times[name].append(duration)
    
    # Print profiling information
    print(f"Kernel: {name}, Duration: {duration:.3f} ms")
    
    return event

def save_profiling_results(output_dir="firecore_profile_results"):
    """Save profiling results to a file"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create results dictionary
    results = {
        "kernel_count": len(kernel_times),
        "kernels": {}
    }
    
    # Add kernel information
    for kernel_name, durations in kernel_times.items():
        if len(durations) > 0:
            results["kernels"][kernel_name] = {
                "executions": len(durations),
                "total_duration": sum(durations),
                "average_duration": sum(durations) / len(durations),
                "min_duration": min(durations),
                "max_duration": max(durations),
                "durations": durations
            }
    
    # Save results to a file
    results_file = os.path.join(output_dir, "firecore_profiling_results.json")
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"Profiling results saved to {results_file}")
    
    # Create a bar chart of kernel durations if there are any kernels
    if len(results["kernels"]) > 0:
        plot_kernel_durations(results, os.path.join(output_dir, "firecore_kernel_durations.png"))
    else:
        print("No kernel executions to plot")
    
    return results

def plot_kernel_durations(results, filename):
    """Plot kernel durations as a bar chart"""
    # Extract kernel names and total durations
    names = []
    durations = []
    
    for kernel_name, kernel_info in results["kernels"].items():
        names.append(kernel_name)
        durations.append(kernel_info["total_duration"])
    
    # Sort by duration
    sorted_data = sorted(zip(names, durations), key=lambda x: x[1], reverse=True)
    names = [x[0] for x in sorted_data]
    durations = [x[1] for x in sorted_data]
    
    # Create a bar chart
    plt.figure(figsize=(12, 6))
    bars = plt.bar(range(len(names)), durations)
    
    # Add labels and title
    plt.xlabel("Kernel")
    plt.ylabel("Total Duration (ms)")
    plt.title("OpenCL Kernel Execution Times")
    plt.xticks(range(len(names)), names, rotation=45, ha="right")
    
    # Add duration values on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f"{height:.2f}", ha="center", va="bottom")
    
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Kernel duration plot saved to {filename}")



#from gpyfft.fft import FFT

COULOMB_CONST  =    14.3996448915 

bytePerFloat = 4


def prepare_grid3d(g0, dg, ng):
    g0 = np.array(g0+(0,), dtype=np.float32)
    dg = np.array(dg+(0,), dtype=np.float32)
    ng = np.array(ng+(0,), dtype=np.int32)
    return (g0, dg, ng)


class GridFF_cl:
    def save_profiling_results(self, output_dir="firecore_profile_results"):
        """Save profiling results to a file"""
        return save_profiling_results(output_dir)


    def __init__(self, nloc=32 ): #, desired_voxel=0.15, allowed_factors={2, 3, 5}):  # <--- Add parameters here
        self.nloc  = nloc
        # self.desired_voxel = desired_voxel       # <--- Store the parameter
        # self.allowed_factors = allowed_factors   # <--- Store the parame

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
            # Try to find the kernel file using an absolute path
            script_dir = os.path.dirname(os.path.abspath(__file__))
            firecore_root = os.path.abspath(os.path.join(script_dir, "../.."))
            kernel_path = os.path.join(firecore_root, "cpp/common_resources/cl/GridFF.cl")
            print(f"Looking for kernel at: {kernel_path}")
            
            with open(kernel_path, 'r') as f:
                self.prg = cl.Program(self.ctx, f.read()).build()
        except Exception as e:
            print( "GridFF_cl() called from path=", os.getcwd() )
            print(f"Error compiling OpenCL program: {e}")
            exit(0)

        # self.atoms_buff    = None
        # self.REQs_buff     = None
        # self.V_Paul_buff   = None
        # self.V_Lond_buff   = None

    def try_make_buff( self, buff_name, sz):
        if hasattr(self,buff_name): 
            buff = getattr(self, buff_name)
            # If the buffer exists but its size is different, release it first to free GPU memory
            if buff is not None and buff.size != sz:
                print(f"Releasing outdated buffer {buff_name} (size {buff.size}) -> realloc {sz}")
                buff.release()
            # If after the potential release the buffer already has the right size, we can reuse it
            if buff is not None and buff.size == sz:
                return
        print( "try_make_buff(",buff_name,") reallocate to [bytes]: ", sz )
        buff = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, sz )
        setattr(self, buff_name, buff )

    def try_buff(self, name, names, sz ):
        if name in names: self.try_make_buff( name+"_buff" , sz)

    def try_make_buffs(self, names, na, nps, bytePerFloat=4):
        # print("try_make_buffs na=", na," nxyz=",  np, " nxyz*bytePerFloat=", nps*bytePerFloat ) 
        print("try_make_buffs na=", na," nps=",  nps, " nps*bytePerFloat=", nps*bytePerFloat )
        self.try_buff("atoms",  names, na*4*bytePerFloat)
        self.try_buff("REQs",   names, na*4*bytePerFloat)
        self.try_buff("V_Paul", names, nps*bytePerFloat)
        self.try_buff("V_Lond", names, nps*bytePerFloat)
        self.try_buff("V_Coul", names, nps*bytePerFloat)
        self.try_buff("Qgrid",  names, nps*2*bytePerFloat)
        self.try_buff("Vgrid",  names, nps*2*bytePerFloat)

        # fitting buffers      
        self.try_buff("Gs",  names, nps*bytePerFloat )
        self.try_buff("fGs", names, nps*bytePerFloat )
        self.try_buff("dGs", names, nps*bytePerFloat )
        self.try_buff("vGs", names, nps*bytePerFloat )

        # New buffers for laplace_real_loop_inert
        self.try_buff("V1", names, nps*bytePerFloat)
        self.try_buff("V2", names, nps*bytePerFloat)
        self.try_buff("vV", names, nps*bytePerFloat)
        
        self.try_buff("FEps", names, nps*4*bytePerFloat)
        self.try_buff("ps",   names, nps*4*bytePerFloat)

        self.try_buff("FE_Paul", names, nps*4*bytePerFloat)
        self.try_buff("FE_Lond", names, nps*4*bytePerFloat)
        self.try_buff("FE_Coul", names, nps*4*bytePerFloat)

        # New buffers for laplace_real_pbc and slabPotential
        # self.try_buff("Vin", names, nxyz*bytePerFloat)
        # self.try_buff("Vout", names, nxyz*bytePerFloat)

    def set_grid(self, gsh : GridShape ):
        self.gsh = gsh
        self.gcl = GridCL(gsh)

    # def prepare_sample3D(self, g0, dg, ng, Eg ):
    #     g0 = np.array(g0+(0,), dtype=np.float32)
    #     dg = np.array(dg+(0,), dtype=np.float32)
    #     ng = np.array(ng+(0,), dtype=np.int32)
    #     #self.grid3D_shape = (g0, dg, ng)
    #     #self.E3D_buf      = cl_array.to_device(self.queue, Eg.astype(np.float32) )

    # def set_grid(self, gsh: GridShape):
    #     # Here we assume that gsh.lvec contains the lattice vectors (each a 3d array).
    #     # We calculate the physical length in each direction.norms.
    #     Ls = [np.linalg.norm(v) for v in gsh.lvec]
    #     ns = []
    #     dgs = []
    #     desired_voxel = 0.15  # Or pass this as a parameter if you wish
    #     for L in Ls:
    #         N, new_voxel = adjust_grid_size(L, desired_voxel)
    #         ns.append(N)
    #         dgs.append(new_voxel)
    #         print(f"For L = {L:.4f}, chosen grid points: {N} with voxel size: {new_voxel:.4f}")
    #     # Update the grid shape values with the FFT-friendly dimensions.
    #     gsh.ns = ns
    #     gsh.dg = dgs
    #     self.gsh = gsh
    #     self.gcl = GridCL(gsh)



    def sample3D_comb(self, ps, C):
        """
        wrapper for 
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

        #self.try_make_buff( FE_buff, n*   )

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
    

    def sample3D(self, ps, buff):
        '''
         wrapper for 
        __kernel void sample3D(
            const float4 g0,
            const float4 dg,
            const int4 ng,
            __global const float* Eg,
            const int n,
            __global const float4* ps,
            __global float4* fes
        '''
        n = len(ps)

        #self.prepare_sample3D( g0, dg, ng, Eg )

        ps_buf  = cl_array.to_device(self.queue, ps.astype(np.float32) )
        fes_buf = cl_array.empty(    self.queue, (n, 4), dtype=np.float32)

        nG = clu.roundup_global_size( n, self.nloc )
        #(g0, dg, ng) = self.grid3D_shape
        

        # print("g0", g0)
        # print("dg", dg)
        # print("ng", ng)
        self.prg.sample3D(self.queue, (nG,), (self.nloc,),  self.gcl.g0, self.gcl.dg, self.gcl.ns, buff, np.int32(n),  ps_buf.data, fes_buf.data )
        fe = fes_buf.get()
        return fe

    def sample3D_grid(self, V_buff, samp_grid, bReturn=True ):
        '''
         wrapper for 
        __kernel void sample3D_grid(
            const float4 g0,
            const float4 dg,
            const int4   ng,
            __global const float* Eg,
            const float4 samp_g0,
            const float4 samp_dg,
            const int4   samp_ng,
            __global float4* fes
        '''

        #self.prepare_sample3D( g0, dg, ng, Eg )
        (samp_g0, samp_dg, samp_ng) = samp_grid
        #ntot = samp_ng[0]*samp_ng[1]*samp_ng[2]
        ntot = samp_ng[3]

        self.try_make_buff( "FE_buff", ntot*4*bytePerFloat )
        #fes_buf = cl_array.empty( self.queue, (ntot, 4), dtype=np.float32)

        nG = clu.roundup_global_size( ntot, self.nloc )
        #(g0, dg, ng) = self.grid3D_shape
        
        self.prg.sample3D_grid(self.queue, (nG,), (self.nloc,),  self.gcl.g0, self.gcl.dg, self.gcl.ns, V_buff, samp_g0, samp_dg, samp_ng, self.FE_buff )
        #fe = fes_buf.get()
        if bReturn:
            sh = samp_ng[:3][::-1]
            FE = np.zeros( (*sh,4), dtype=np.float32)
            cl.enqueue_copy(self.queue, FE, self.FE_buff )

            print( "GridFF_cl::sample3D_grid() FE min max ",  FE.min(), FE.max() )
            return FE


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

    # def allocate_cl_buffers(self, ns=None, E0=None, nByte=4):
    #     print(f"GridFF_cl::allocate_cl_buffers().1 Queue: {self.queue}, Context: {self.ctx}")
    #     if ns is None:
    #         ns = E0.shape
    #     ntot = np.int32(ns[0] * ns[1] * ns[2])
    #     self.G0_buf  = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=ntot *nByte)   # Reference
    #     self.Gs_buf  = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)
    #     self.fGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # variational derivative
    #     self.dGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # fitting error
    #     self.vGs_buf = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=ntot *nByte)   # velocity of Gs update
    #     #self.E3D_tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, cl.ImageFormat(cl.channel_order.R, cl.channel_type.FLOAT),  shape=ns  )
    #     if E0 is not None:
    #         cl.enqueue_copy(self.queue, self.G0_buf,  E0.astype(np.float32))
    #         cl.enqueue_copy(self.queue, self.Gs_buf,  E0.astype(np.float32))
    #         #cl.enqueue_copy(self.queue, self.E3D_tex, E0.astype(np.float32))
    #     self.ns = np.array(ns, dtype=np.int32)
    #     print(f"GridFF_cl::allocate_cl_buffers().2 Queue: {self.queue}, Context: {self.ctx}")

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
    
    def fit3D_with_buffer(self, buffer, nPerStep=10, nmaxiter=300, dt=0.5, Ftol=1e-16, damp=0.15, bConvTrj=False, bReturn=True, bPrint=False, bTime=True, bDebug=True):
        """
        A wrapper around fit3D that handles buffer size mismatches
        """
        # Get the actual size of the input buffer in bytes
        buffer_size_bytes = buffer.size
        
        # Calculate how many float32 elements that corresponds to
        num_elements = buffer_size_bytes // 4  # 4 bytes per float32
        
        # Check if this matches our expected grid size
        expected_elements = self.gcl.nxyz
        
        if num_elements != expected_elements:
            print(f"Warning: Buffer size ({num_elements}) doesn't match expected grid size ({expected_elements})")
            print(f"Creating a temporary buffer of the correct size...")
            
            # Create a temporary buffer of the correct size
            temp_buff = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=expected_elements * 4)
            
            # Determine how many elements we can safely copy
            copy_elements = min(num_elements, expected_elements)
            
            # Use a kernel to copy the data (up to the smaller size)
            nL = self.nloc
            nG = clu.roundup_global_size(copy_elements, nL)
            
            # Create a simple copy kernel if needed
            if not hasattr(self.prg, 'copyPartialBuffer'):
                copy_kernel = """
                __kernel void copyPartialBuffer(
                    const int n,
                    __global const float* src,
                    __global float* dst
                ) {
                    int i = get_global_id(0);
                    if (i < n) {
                        dst[i] = src[i];
                    }
                }
                """
                self.prg = cl.Program(self.ctx, self.prg.get_info(cl.program_info.SOURCE) + copy_kernel).build()
            
            # Copy the data
            self.prg.copyPartialBuffer(self.queue, (nG,), (nL,), np.int32(copy_elements), buffer, temp_buff)
            self.queue.finish()
            
            # Use the temporary buffer for fitting
            result, trj = self.fit3D(temp_buff, nPerStep=nPerStep, nmaxiter=nmaxiter, dt=dt, 
                                    Ftol=Ftol, damp=damp, bConvTrj=bConvTrj, 
                                    bReturn=bReturn, bPrint=bPrint, bTime=bTime, bDebug=bDebug)
            
            return result, trj
        else:
            # If sizes match, just call the original fit3D
            return self.fit3D(buffer, nPerStep=nPerStep, nmaxiter=nmaxiter, dt=dt, 
                            Ftol=Ftol, damp=damp, bConvTrj=bConvTrj, 
                            bReturn=bReturn, bPrint=bPrint, bTime=bTime, bDebug=bDebug)



    def fit3D(self, Ref_buff, nmaxiter=300, dt=0.5, Ftol=1e-16, damp=0.15, nPerStep=50, bConvTrj=False, bReturn=True, bPrint=False, bTime=True, bDebug=True ):
        # NOTE / TODO : It is a bit strange than GridFF.h::makeGridFF_Bspline_d() the fit is fastes with damp=0.0 but here damp=0.15 performs better
        #print(f"GridFF_cl::fit3D().1 Queue: {self.queue}, Context: {self.ctx}")
        print(f"GridFF_cl::fit3D() dt={dt}, damp={damp} nmaxiter={nmaxiter} Ftol{Ftol}")
        T00=time.perf_counter()
        cdamp=1.-damp

        ns   = self.gsh.ns #[::-1] 
        nxyz = self.gcl.nxyz

        buff_names={'Gs','dGs','fGs','vGs'}
        self.try_make_buffs( buff_names, 0, nxyz )

        ns_cl = np.array( ns+[0], dtype=np.int32 )
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

        profile_kernel("setMul", self.queue, self.prg.setMul, (nG,), (nL,), nxyz,  Ref_buff,  self.Gs_buff,  np.float32(1.0) ) # setup force
        profile_kernel("setMul", self.queue, self.prg.setMul, (nG,), (nL,), nxyz,  Ref_buff,  self.vGs_buff, np.float32(0.0) ) # setup velocity
        
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
        # self.V_Coul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )
    
    def prepare_Coulomb_buffers(self, na, nxyz, bytePerFloat=4 ):
        #print( "GridFF_cl::prepare_Coulomb_buffers() " )
        self.atoms_buff  = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        self.REQs_buff   = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY,  size=na*4*bytePerFloat)
        self.Qgrid_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )
        self.V_Coul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )


    # def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=None, ng=None, lvec=None, g0=None, GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True ):

    #     T00 = time.perf_counter()

    #     # Only create a new grid if one doesn't exist or if parameters are explicitly provided
    #     # if self.gsh is None or dg is not None or ng is not None or lvec is not None or g0 is not None:
    #     #     # Use existing grid parameters if not explicitly provided
    #     #     if dg is None and hasattr(self, 'gsh'):
    #     #         dg = self.gsh.dg
    #     #         print(f"Using existing grid step size: {dg}")
    #     #     if lvec is None and hasattr(self, 'gsh') and self.gsh.lvec is not None:
    #     #         lvec = self.gsh.lvec
    #     #         print(f"Using existing lattice vectors: {lvec}")
    #     #     if g0 is None and hasattr(self, 'gsh'):
    #     #         g0 = self.gsh.g0
    #     #         print(f"Using existing grid origin: {g0}")
            
    #     #     grid = GridShape(ns=ng, dg=dg, lvec=lvec, g0=g0)
    #     #     print(f"Creating grid with: ns={grid.ns}, dg={grid.dg}, lvec={lvec}")
    #     #     self.set_grid(grid)
        
    #     # atoms.apos = xyzq[:,:3]
    #     # print("New_Atoms:",atoms.shape)
    #     # print("xyzq[:,:3] =", atoms[:,:3])
    #     print("GridFF_cl::make_MorseFF() dg = ",dg)


    #     na   = len(atoms)
    #     nxyz = self.gcl.nxyz

    #     buff_names={'atoms','REQs','V_Paul','V_Lond'}
    #     self.try_make_buffs(buff_names, na, nxyz)

    #     atoms_np = np.asarray(atoms, dtype=np.float32)
    #     reqs_np  = np.asarray(REQs, dtype=np.float32)
    #     cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)
    #     cl.enqueue_copy(self.queue, self.REQs_buff, reqs_np)

    #     na_cl        = np.int32(na)
    #     nPBC_cl      = np.array(nPBC+(0,), dtype=np.int32   )
    #     GFFParams_cl = np.array(GFFParams, dtype=np.float32 )

    #     nL = self.nloc
    #     nG = clu.roundup_global_size(nxyz, nL)

    #     T0 = time.perf_counter()
    #     self.prg.make_MorseFF(self.queue, (nG,), (nL,),
    #                         na_cl, self.atoms_buff, self.REQs_buff, self.V_Paul_buff, self.V_Lond_buff,
    #                         nPBC_cl, self.gcl.ns, self.gcl.a,self.gcl.b,self.gcl.c, self.gcl.g0, GFFParams_cl)
    
    #     if bTime:
    #         self.queue.finish()
    #         dT=time.perf_counter()-T0
    #         npbc = float((nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1))
    #         nops = float(na_cl) * float(nxyz) * npbc
    #         print( "GridFF_cl::make_MorseFF() time[s] ", dT,   "  preparation[s]  ", T0-T00, "[s] nGOPs: ", nops*1e-9," speed[GOPs/s]: ", (nops*1e-9)/dT , " na,nxyz,npbc ", na_cl, nxyz, npbc  )

    #     if bReturn:
    #         sh = self.gsh.ns[::-1]
    #         V_Paul = np.zeros( sh, dtype=np.float32)
    #         V_Lond = np.zeros( sh, dtype=np.float32)
    #         cl.enqueue_copy(self.queue, V_Paul, self.V_Paul_buff)
    #         cl.enqueue_copy(self.queue, V_Lond, self.V_Lond_buff)
    #         # print( "V_Lond.min,max ", V_Lond.min(), V_Lond.max() )
    #         # print( "V_Paul.min,max ", V_Paul.min(), V_Paul.max() )
    #         return V_Paul, V_Lond
    
    def make_MorseFF_f4(self, atoms, REQs, nPBC=(4, 4, 0), dg=(0.1, 0.1, 0.1), ng=None,            lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]],                     g0=(0.0, 0.0, 0.0), GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True ):

        T00 = time.perf_counter()

        grid = GridShape( ns=ng, dg=dg, lvec=lvec, g0=g0 )
        self.set_grid( grid )

        na   = len(atoms)
        nxyz = self.gcl.nxyz

        buff_names={'atoms','REQs','FE_Paul','FE_Lond'}
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
        self.prg.make_MorseFF_f4(self.queue, (nG,), (nL,),
                            na_cl, self.atoms_buff, self.REQs_buff, self.FE_Paul_buff, self.FE_Lond_buff,
                            nPBC_cl, self.gcl.ns, self.gcl.a,self.gcl.b,self.gcl.c, self.gcl.g0, GFFParams_cl)
    
        if bTime:
            self.queue.finish()
            dT=time.perf_counter()-T0
            npbc = (nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
            nops = na_cl * nxyz * npbc
            print( "GridFF_cl::make_MorseFF() time[s] ", dT,   "  preparation[s]  ", T0-T00, "[s] nGOPs: ", nops*1e-9," speed[GOPs/s]: ", (nops*1e-9)/dT , " na,nxyz,npbc ", na_cl, nxyz, npbc  )

        if bReturn:
            sh = self.gsh.ns[::-1]
            FE_Paul = np.zeros( sh+(4,), dtype=np.float32)
            FE_Lond = np.zeros( sh+(4,), dtype=np.float32)
            cl.enqueue_copy(self.queue, FE_Paul, self.FE_Paul_buff)
            cl.enqueue_copy(self.queue, FE_Lond, self.FE_Lond_buff)
            # print( "V_Lond.min,max ", V_Lond.min(), V_Lond.max() )
            # print( "V_Paul.min,max ", V_Paul.min(), V_Paul.max() )
            return FE_Paul, FE_Lond

    def make_Coulomb_points(self, atoms, ps, nPBC=(25, 25, 0), lvec=None, Ls=None, GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True):
        if Ls   is None: Ls   = self.gsh.Ls 
        if lvec is None: lvec = [[Ls[0],0.0,0.0],[0.0,Ls[1],0.0],[0.0,0.0,Ls[2]]]
        T00 = time.perf_counter()

        na  = len(atoms)
        nps = len(ps)

        # Buffer creation: atoms, ps, and FE_Coul buffers
        buff_names = {'atoms', 'ps', 'FEps'}
        self.try_make_buffs(buff_names, na, nps )

        # nxyz = self.gcl.nxyz
        # self.V_Coul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=nxyz*bytePerFloat )

        # Convert input data to NumPy arrays and enqueue them to the device buffers
        atoms_np = np.asarray(atoms, dtype=np.float32 )
        ps_np    = np.asarray(ps,    dtype=np.float32 )
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)
        cl.enqueue_copy(self.queue, self.ps_buff,    ps_np)

        #nPBC=(20,20,0)
        # Set up kernel parameters
        na_cl = np.int32(na)
        np_cl = np.int32(nps)
        nPBC_cl = np.array(nPBC + (0,), dtype=np.int32)
        #lvec_a = np.array( lvec, dtype=np.float32)
        #lvec_b = np.array( lvec, dtype=np.float32)
        lvec_c = np.array( [lvec[2][0],lvec[2][1],lvec[2][2],0.0], dtype=np.float32)
        GFFParams_cl = np.array(GFFParams, dtype=np.float32)

        # Work-group sizes
        nL = self.nloc
        nG = clu.roundup_global_size(nps, nL)

        # Kernel invocation
        T0 = time.perf_counter()
        self.prg.make_Coulomb_points(self.queue, (nG,), (nL,),
                                na_cl, np_cl, self.atoms_buff, self.ps_buff, self.FEps_buff,
                                nPBC_cl, self.gcl.a, self.gcl.b, lvec_c, GFFParams_cl)

        if bTime:
            self.queue.finish()
            dT = time.perf_counter() - T0
            npbc = (nPBC[0] * 2 + 1) * (nPBC[1] * 2 + 1) * (nPBC[2] * 2 + 1)
            nops = na_cl * nps * npbc
            print("GridFF_cl::make_Coulomb_points() time[s]:", dT, " preparation[s]:", T0 - T00, "[s] nGOPs:", nops * 1e-9, " speed[GOPs/s]:", (nops * 1e-9) / dT, " na, nps, npbc:", na_cl, nps, npbc)

        # Return results if needed
        if bReturn:
            FE_Coul = np.zeros( (nps,4), dtype=np.float32)
            cl.enqueue_copy(self.queue, FE_Coul, self.FEps_buff)
            return FE_Coul


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
        print("Debug Project Atoms on Grid Quintic PBC [%f, %f, %f, %f]" % (self.gcl.g0[0], self.gcl.g0[1], self.gcl.g0[2], self.gcl.g0[3]))
        nxyz2 = np.int32( ns[0]*ns[1]*ns[2] * 2 )
        self.prg.set( self.queue, sz_glob, sz_loc, nxyz2, self.Qgrid_buff, np.float32(0.0) )
        self.prg.project_atoms_on_grid_quintic_pbc( self.queue, sz_glob, sz_loc,
            na, self.atoms_buff, self.Qgrid_buff,
            ns, self.gcl.g0, self.gcl.dg
        )

    def project_atoms_on_grid_quintic_pbc(self, atoms, ng=None,   dg=(0.1, 0.1, 0.1),    lvec=[[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]], g0=(0.0, 0.0, 0.0), bReturn=True ):

        g0=self.gcl.g0
        grid = GridShape( ns=ng, dg=dg, lvec=lvec, g0=g0 )
        self.set_grid( grid )

        print("Debug Project Atoms on Grid Quintic PBC %f %f %f" % g0 )


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

    def _poisson_old(self, sz_glob, sz_loc, ns, coefs_cl ):
        event, = self.transform.enqueue() 
        event.wait()
        self.prg.poissonW_old( self.queue, sz_glob, sz_loc,       ns, self.Qgrid_buff, self.Vgrid_buff, coefs_cl )
        event, = self.inverse_transform.enqueue()
        event.wait()

    def prepare_poisson(self, sh=None):
        if sh is None: sh=self.gsh.ns[::-1] 
        print("GridFF_cl::prepare_poisson() sh=", sh )
        self.Qgrid_cla = cl_array.Array(self.queue, shape=tuple(sh), dtype=np.complex64, data=self.Qgrid_buff )
        print("FFT is about to be created with grid shape:", self.Qgrid_cla.shape)
        self.Vgrid_cla = cl_array.Array(self.queue, shape=tuple(sh), dtype=np.complex64, data=self.Vgrid_buff )
        self.transform         = clu.FFT(self.ctx, self.queue, self.Qgrid_cla, axes=(0, 1, 2))
        self.inverse_transform = clu.FFT(self.ctx, self.queue, self.Vgrid_cla, axes=(0, 1, 2))

    def poisson_old(self, bReturn=True, sh=None, dV=None ):

        clu.try_load_clFFT()
        if sh is None: sh=self.gsh.ns[::-1]
        nxyz = np.int32( sh[0]*sh[1]*sh[2] )

        if dV is None: dV = self.gsh.dV

        buff_names={'Qgrid','Vgrid'}
        self.try_make_buffs(buff_names, 0, nxyz)
        self.prepare_poisson( sh=sh )

        ns_cl    = np.array( (*sh[::-1],nxyz), dtype=np.int32 )
        #sc_ewald = 4*np.pi*dV*COULOMB_CONST
        sc_ewald = 4*np.pi*dV*dV*COULOMB_CONST   # This Seems wrong

        #sc_ewald = ( 4.0*np.pi*dV*COULOMB_CONST ) / ( nxyz )    # This normalization makes physical sense, and the result is independnet on number of grid points, however, the poential is too small (much smaller than real space reference)
        # See discussion of normalization constant here: https://chatgpt.com/share/672d095b-9304-8012-ac50-de5cbe5cba18
        # also chack how normalization constant is done in /home/prokop/git/FireCore/cpp/common/molecular/EwaldGrid.h   EwaldGrid::laplace_reciprocal_kernel()
        # see https://github.com/ProkopHapala/FireCore/blob/079f9244f486d365362e1766b837e3aa4df6d368/cpp/common/molecular/EwaldGrid.h#L502

        coefs_cl = np.array( (0.0,0.0,0.0, sc_ewald ), dtype=np.float32 )

        nL = self.nloc
        nG = clu.roundup_global_size( nxyz, nL)
        self._poisson_old( (nG,), (nL,), ns_cl, coefs_cl )        

        if bReturn:
            sh_ = (*sh,2)
            print( "sh ", sh_ )
            Vgrid = np.zeros( sh_, dtype=np.float32   )  
            cl.enqueue_copy ( self.queue, Vgrid, self.Vgrid_buff )
            print( "Vgrid min,max ", Vgrid.min(), Vgrid.max(), " sc_ewald ", sc_ewald, " dV ", dV, " nxyz ", nxyz )
            return Vgrid

    def poisson(self, bReturn=True, sh=None, dV=None):
        if sh is None:
            sh = self.gsh.ns[::-1]
        nxyz = np.int32(sh[0] * sh[1] * sh[2])

        if dV is None:
            dV = self.gsh.dV

        buff_names = {'Qgrid', 'Vgrid'}
        self.try_make_buffs(buff_names, 0, nxyz)
        self.prepare_poisson(sh=sh)

        ns_cl = np.array((*sh[::-1], nxyz), dtype=np.int32)

        # Compute frequencies
        Lx = sh[2] * self.gsh.dg[0]
        Ly = sh[1] * self.gsh.dg[1]
        Lz = sh[0] * self.gsh.dg[2]

        freq_x = 2.0 * np.pi / Lx
        freq_y = 2.0 * np.pi / Ly
        freq_z = 2.0 * np.pi / Lz

        scEwald = COULOMB_CONST * 4.0 * np.pi / (nxyz * dV)

        #coefs_cl = np.array((freq_x, freq_y, freq_z, 0.0), dtype=np.float32)
        coefs_cl = np.array((freq_x, freq_y, freq_z, scEwald), dtype=np.float32)

        nL = self.nloc
        nG = clu.roundup_global_size(nxyz, nL)
        self._poisson((nG,), (nL,), ns_cl, coefs_cl)

        if bReturn:
            sh_ = (*sh, 2)
            Vgrid = np.zeros(sh_, dtype=np.float32)
            cl.enqueue_copy(self.queue, Vgrid, self.Vgrid_buff)
            # Apply scaling factor after inverse FFT
            scEwald = COULOMB_CONST * 4.0 * np.pi / (nxyz * dV)
            Vgrid *= scEwald
            print("Vgrid min,max ", Vgrid.min(), Vgrid.max(), " scEwald ", scEwald, " dV ", dV, " nxyz ", nxyz)
            return Vgrid



    def slabPotential_old(self, nz_slab, dz, Vol, dVcor, Vcor0, bDownload=True):
        buff_names = {'Vgrid', 'V_Coul'}
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
        
    def slabPotential(self,  Vin_buff, nz_slab, dipol=0.0,  bDownload=True, bTranspose=False):
        dz  = self.gsh.dg[2]
        Lz  = self.gsh.Ls[2] 
        dL_slab = nz_slab * dz
        Lz_slab = Lz + dL_slab
        Vol = self.gsh.V * Lz_slab/Lz
        dVcor  = 4.0 * np.pi * COULOMB_CONST * dipol/Vol;
        Vcor0 = -dVcor * Lz_slab/2;
        print( "dz, Lz, dL_slab, Lz_slab, Vol, dVcor, Vcor0 ,g0 ,g0", dz, Lz, dL_slab, Lz_slab, Vol, dVcor, Vcor0, self.gcl.g0,self.gsh.g0 )
        buff_names = {'V_Coul'}
        self.try_make_buffs(buff_names, 0, self.gcl.nxyz )
        # self.V_Coul_buff = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, size=self.gcl.nxyz*bytePerFloat )

        
        raw_nz = self.gsh.ns[2] + nz_slab
        raw_nz_int = int(np.ceil(raw_nz))
        adj_nz = clu.next_nice(raw_nz_int, allowed_factors={2, 3, 5})
        ns_cl = np.array([self.gsh.ns[0], self.gsh.ns[1], self.gsh.ns[2], adj_nz-self.gsh.ns[2]], dtype=np.int32)
        #ns_cl = np.array([self.gsh.ns[0], self.gsh.ns[1], adj_nz, 0], dtype=np.int32)
        # ns_cl = np.array([self.gsh.ns[0], self.gsh.ns[1], adj_nz-nz_slab, (nz_slab) ], dtype=np.int32)
        # ns_cl = np.array((self.gsh.ns[0], self.gsh.ns[1], self.gsh.ns[2], nz_slab), dtype=np.int32)
        # ns_cl  = np.array( self.gsh.ns+(nz_slab,),  dtype=np.int32   )
        # Debug dimensions
        print(f"Buffer dimensions: {self.gsh.ns}")
        print(f"Kernel parameters: ns_cl={ns_cl}")

        params = np.array( [dz, Vol, dVcor, Vcor0], dtype=np.float32 )
        sz_loc  = (4,4,4,)
        sz_glob = clu.roundup_global_size_3d( self.gsh.ns, sz_loc)
        print("!!!!!! slabPotential()  bTranspose  ", bTranspose )
        if bTranspose:
            print( "GridFF_cl::slabPotential() Transpose => slabPotential_zyx" )
            event=self.prg.slabPotential_zyx( self.queue, sz_glob, sz_loc, ns_cl, Vin_buff, self.V_Coul_buff, params )
        else:
            event=self.prg.slabPotential( self.queue, sz_glob, sz_loc, ns_cl, Vin_buff, self.V_Coul_buff, params )
        event.wait()
        self.queue.finish()
        slabPotential_zyx = np.empty(self.gsh.ns[::-1], dtype=np.float32)
        cl.enqueue_copy(self.queue, slabPotential_zyx, self.V_Coul_buff).wait()
        
        # # Plot three slices through the middle
        # plt.figure(figsize=(15, 5))
        # mid_x = slabPotential_zyx.shape[2]//2
        # mid_y = slabPotential_zyx.shape[1]//2
        # mid_z = slabPotential_zyx.shape[0]//2
        # # mid_x= 0
        # # mid_y= 40
        # # mid_z= 2
        # # XY plane
        # plt.subplot(131)
        # plt.imshow(slabPotential_zyx[mid_z, :, :], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XY plane (z={mid_z})')
        # # XZ plane
        # plt.subplot(132)
        # plt.imshow(slabPotential_zyx[:, mid_y, :], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XZ plane (y={mid_y})')
        # # YZ plane
        # plt.subplot(133)
        # plt.imshow(slabPotential_zyx[:, :, mid_x], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'YZ plane (x={mid_x})')
        # plt.suptitle('slabPotential_zyx')
        # plt.tight_layout()
        # plt.show()
        
     


        # ===== DEBUG: Output Buffer Check =====
        debug_check = np.empty(self.gsh.ns[::-1], dtype=np.float32)
        cl.enqueue_copy(self.queue, debug_check, self.V_Coul_buff)
        print(f"V_Coul_buff state in slabPotential: range [{debug_check.min():.3f}, {debug_check.max():.3f}]")
        if bDownload:
            if bTranspose:
                V_Coul = np.empty( self.gsh.ns, dtype=np.float32  )
            else:
                V_Coul = np.empty( self.gsh.ns[::-1], dtype=np.float32  )

            # Debug before copy
            debug_before = np.empty_like(V_Coul)
            cl.enqueue_copy(self.queue, debug_before, self.V_Coul_buff)
            self.queue.finish()
            print("\nDebug V_Coul_buff before final copy:")
            print(f"Shape: {debug_before.shape}")
            print(f"Range: {debug_before.min():.3f} to {debug_before.max():.3f}")

            cl.enqueue_copy(self.queue, V_Coul, self.V_Coul_buff)
            # cl.enqueue_copy(self.queue, self.V_Coul_buff, V_Coul)
            self.queue.finish()

            # Debug after copy
            debug_after = np.empty_like(V_Coul)
            cl.enqueue_copy(self.queue, debug_after, self.V_Coul_buff)
            self.queue.finish()
            print("\nDebug V_Coul_buff after final copy:")
            print(f"Shape: {debug_after.shape}")
            print(f"Range: {debug_after.min():.3f} to {debug_after.max():.3f}")
            
            return V_Coul
        
        # if bDownload:
        #     result = np.empty(self.gsh.ns[::-1] if bTranspose else self.gsh.ns, dtype=np.float32)
        #     evt = cl.enqueue_copy(self.queue, result, self.V_Coul_buff)
        #     evt.wait()
        #     return result
        # return None

    def laplace_real_loop_inert(self, niter=16, cSOR=0.0, cV=0.8, bReturn=False, sh=None ):
        print( "GridFF_cl::laplace_real_loop_inert() " )

        if sh is None: sh=self.gsh.ns[::-1]
        nxyz = np.int32( sh[0]*sh[1]*sh[2] )
        print(f"Grid shape (ZYX): {sh}, total elements: {nxyz}")

        buff_names = {'V1', 'V2', 'vV'}
        self.try_make_buffs(buff_names, 0, nxyz )

        sz_loc = (4,4,4)
        sz_glob = clu.roundup_global_size_3d( sh[::-1], sz_loc )

        ns_cl   = np.array( (*sh[::-1],0), dtype=np.int32 )
        C_cl    = np.array( (1.0,0.0), dtype=np.float32 )
        cV_cl   = np.float32(cV)
        cSOR_cl = np.float32(cSOR)

        szl = (self.nloc,)
        szg = ( clu.roundup_global_size( nxyz, self.nloc), )

        self.prg.setCMul( self.queue, szg,  szl, nxyz, self.Vgrid_buff,  self.V1_buff, C_cl ) # copy real part from complex Vgrid to scalar V1
        self.prg.set    ( self.queue, szg,  szl, nxyz, self.vV_buff, np.float32(0.0) )        # initialize velocity to zero
        last_buff=1
        #self.prg.laplace_real_pbc(     self.queue, sz_glob, sz_loc,       ns_cl, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl, np.float32(0.0) )
        self.prg.laplace_real_pbc(     self.queue, sz_glob, sz_loc,       ns_cl, self.V1_buff, self.V2_buff, None,         cSOR_cl, np.float32(0.0) )
        last_buff=2

        for iter in range(niter):
            if iter % 2 == 0:
                self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,    ns_cl, self.V2_buff, self.V1_buff, self.vV_buff, cSOR_cl,cV_cl )
                last_buff=1
            else:
                self.prg.laplace_real_pbc( self.queue, sz_glob, sz_loc,    ns_cl, self.V1_buff, self.V2_buff, self.vV_buff, cSOR_cl,cV_cl )
                last_buff=2

        if bReturn:
            V = np.empty(sh, dtype=np.float32)
            nxyz_ = V.shape[0] * V.shape[1] * V.shape[2]
            print(f"V shape: {V.shape}|{nxyz_},{nxyz}, V2_buff size: {self.V2_buff.size}")
            if last_buff == 2:
                cl.enqueue_copy(self.queue, V, self.V2_buff)
            elif last_buff == 1:
                cl.enqueue_copy(self.queue, V, self.V1_buff)
                
            # DEBUG: Analyze Laplace output to find max value position
            max_idx = np.unravel_index(np.argmax(np.abs(V)), V.shape)
            max_val = V[max_idx]
            print(f"\n=== DEBUG: Laplace Output Analysis ====")
            print(f"Laplace output shape: {V.shape} (ZYX order)")
            print(f"Laplace output max at (z,y,x): {max_idx}, value: {max_val:.6f}")
            return V
        else:
            if last_buff == 2:
                return self.V2_buff
            elif last_buff == 1:
                return self.V1_buff
    
    def makeCoulombEwald(self, atoms, bOld=False ):
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
    
        if bOld:
            Vgrid = self.poisson_old()
        else:
            Vgrid = self.poisson()

        Vgrid = Vgrid[:,:,:,0].copy()
        #print( "Vgrid min,max", Vgrid.min(), Vgrid.max() )
        return Vgrid


    def makeCoulombEwald_slab(self, atoms, Lz_slab=20.0, dipol=0.0, niter=16, bDipoleCoorection=False, bReturn=True, bTranspose=False, bSaveQgrid=False, bCheckVin=False, bCheckPoisson=False ):
        print( f"GridFF_cl::makeCoulombEwald_slab()  Lz_slab {Lz_slab}, dipol {dipol} niter {niter} bDipoleCoorectio {bDipoleCoorection} bReturn {bReturn} bTranspose {bTranspose} bSaveQgrid {bSaveQgrid} bCheckVin{bCheckVin} bCheckPoisson {bCheckPoisson} " )
        
        clu.try_load_clFFT()
        if self.gcl is None: 
            print("ERROR in GridFF_cl::makeCoulombEwald() gcl is None, => please call set_grid() first " )
            exit()

        nz_slab   = Lz_slab/self.gsh.dg[2]; print("nz_slab", nz_slab, " ns[2] ", self.gcl.ns[2] )
        dz = self.gsh.dg[2]; print("Using grid step size dz =", dz)        
        raw_nz = self.gsh.ns[2] + nz_slab
        raw_nz_int = int(np.ceil(raw_nz))
        adj_nz = clu.next_nice(raw_nz_int, allowed_factors={2, 3, 5})
        ns_cl = np.array([self.gsh.ns[0], self.gsh.ns[1], adj_nz, adj_nz-self.gsh.ns[2]], dtype=np.int32); print("ns_cl", ns_cl)

        nxyz_slab = ns_cl[0]*ns_cl[1]*ns_cl[2]

        #if sh is None: sh=self.gsh.ns[::-1] 
        sh=ns_cl[:3][::-1] 

        na = len(atoms)
        buff_names={'atoms','Qgrid','Vgrid'}

        sz_loc  = (self.nloc,)
        sz_glob = ( clu.roundup_global_size( nxyz_slab, self.nloc), )

        self.try_make_buffs(buff_names, na, nxyz_slab )
        atoms_np = np.array(atoms, dtype=np.float32)
        print("Checking GridFF_cl::makeCoulombEwald_slab() dg", self.gsh.dg)
        print("Grid dimensions (ns):", self.gsh.ns)
        print("Grid lattice vectors (Ls):", self.gsh.Ls)
        # nxyz = self.gcl.nxyz
        # na = len(atoms)
        # self.prepare_Coulomb_buffers(na, nxyz)
        # self.try_make_buffs('V_Coul', na, nxyz_slab )

        # NOTE / TODO : This is strange, not sure why we need to shift the coordinates
        # atoms_np[:,0] += self.gcl.dg[0] * 1 
        # atoms_np[:,1] += self.gcl.dg[1] * 1
        # #atoms_np[:,2] += self.gcl.dg[2] * (-1 + 4 - 0.075)   # NOTE / TODO : shift -1 is because we do the same shift on CPU (in GridFF::tryLoadGridFF_potentials() ), this make sure that Qgrid_gpu and Qgrid_cpu are the same.   But align V_Coul with CPU we need to shift by 4 (see below)

        # atoms_np[:,2] += self.gcl.dg[2] * 1      # NOTE / TODO : This works with new normalization (tested by compare_npy_z.py with respect to C++ EwaldGrid.h for NaCl_1x1_L3 ) ( see GridFF_cl::poisson()  vs GridFF_cl::poisson_old() )
        #atoms_np[:,2] += self.gcl.dg[2] * (-1 )
        #atoms_np[:,2] += self.gcl.dg[2] * 1
        
        atoms_np[:,0] += self.gcl.dg[0] * 1
        atoms_np[:,1] += self.gcl.dg[1] * 1
        atoms_np[:,2] += self.gcl.dg[2] * 1
        
        cl.enqueue_copy(self.queue, self.atoms_buff, atoms_np)

        print("GridFF_cl::makeCoulombEwald_slab()._project_atoms_on_grid_quintic_pbc")
        self._project_atoms_on_grid_quintic_pbc(sz_glob, sz_loc, np.int32(na), ns_cl)

        if bSaveQgrid: 
            #sh    = self.gsh.ns[::-1]
            Qgrid = np.zeros( (*sh,2,), dtype=np.float32 )
            cl.enqueue_copy(self.queue, Qgrid, self.Qgrid_buff )
            # print("Qgrid min,max ", Qgrid[:,:,:,0].min(), Qgrid[:,:,:,0].max() )
            # print("\n--- Specific Z-slices from Full Qgrid (GPU data) ---")
            # for z_in in [0, 200, 599, 799]:
            #     print(f"Qgrid[Z_in={z_in}] range: [{Qgrid[z_in,:,:,0].min():.6f}, {Qgrid[z_in,:,:,0].max():.6f}]")
            # print("------------------------------------------------------")
            # fname = "NaCl_1x1_L1/"
            # fname = "Na_0.9_Cl_-0.9/"
            # path = "./data/" + os.path.splitext( fname )[0] + "Qgrid_gpu.npy"
            # print(f"Saving Qgrid to {path}")
            # np.save( path, Qgrid[:,:,:,0] )
            ####### Plot three slices through the middle
            # plt.figure(figsize=(15, 5))
            # # mid_x= 0
            # # mid_y= 40
            # # mid_z= 2
            # mid_x = Qgrid.shape[2] // 2
            # mid_y = Qgrid.shape[1] // 2
            # mid_z = Qgrid.shape[0] // 2
            
            # # XY plane
            # plt.subplot(131)
            # plt.imshow(Qgrid[mid_z, :, :, 0], origin='lower')
            # plt.colorbar(label='Charge density')
            # plt.title(f'XY plane (z={mid_z})')
            
            # # XZ plane
            # plt.subplot(132)
            # plt.imshow(Qgrid[:, mid_y, :, 0], origin='lower')
            # plt.colorbar(label='Charge density')
            # plt.title(f'XZ plane (y={mid_y})')
            
            # # YZ plane
            # plt.subplot(133)
            # plt.imshow(Qgrid[:, :, mid_x, 0], origin='lower')
            # plt.colorbar(label='Charge density')
            # plt.title(f'YZ plane (x={mid_x})')
            
            # plt.suptitle('Atom Projection on Grid (Charge Density)')
            # plt.tight_layout()
            # # plt.show()
            
        
        self.poisson( bReturn=bCheckPoisson, sh=sh )
        Vin_buff = self.laplace_real_loop_inert( bReturn=False, niter=niter, sh=sh )

        # Visualize V1_buff (extended grid with nz_slab)
        V_after_poisson_laplace = np.empty(sh[::-1], dtype=np.float32)  # Note: sh includes extended z-dimension
        cl.enqueue_copy(self.queue, V_after_poisson_laplace, Vin_buff)
        
        ####### Plot three slices through the middle
        # plt.figure(figsize=(15, 5))
        # # mid_x= 0
        # # mid_y= 40
        # # mid_z= 7
        # mid_x = V_after_poisson_laplace.shape[2] // 2
        # mid_y = V_after_poisson_laplace.shape[1] // 2
        # mid_z = V_after_poisson_laplace.shape[0] // 2
        
        
        # # XY plane
        # plt.subplot(131)
        # plt.imshow(V_after_poisson_laplace[mid_z, :, :], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XY plane (z={mid_z})')
        
        # # XZ plane
        # plt.subplot(132)
        # plt.imshow(V_after_poisson_laplace[:, mid_y, :], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XZ plane (y={mid_y})')
        
        # # YZ plane
        # plt.subplot(133)
        # plt.imshow(V_after_poisson_laplace[ :, :,mid_x], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'YZ plane (x={mid_x})')
        
        # plt.suptitle('Potential After Poisson and Laplace')
        # plt.tight_layout()
        # # plt.show()
            
        # Before calling slabPotential
        print("\nDebug before slabPotential:")
        if hasattr(self, 'V_Coul_buff'):
            debug_before = np.empty(self.gsh.ns[::-1], dtype=np.float32)
            cl.enqueue_copy(self.queue, debug_before, self.V_Coul_buff)
            print(f"V_Coul_buff exists with shape: {debug_before.shape}")
            print(f"Range: {debug_before.min():.3f} to {debug_before.max():.3f}")
        else:
            print("V_Coul_buff does not exist yet")

        V_after_slab = self.slabPotential( Vin_buff, nz_slab, dipol=dipol, bDownload=bReturn, bTranspose=bTranspose )
        # V_after_slab = np.empty(sh[::-1], dtype=np.float32)  # Note: sh includes extended z-dimension
        # cl.enqueue_copy(self.queue, V_after_slab, V_Coul_test)
        # verify_vcoul = np.empty(self.gsh.ns[::-1], dtype=np.float32)
        # cl.enqueue_copy(self.queue, verify_vcoul, self.V_Coul_buff)
        # self.queue.finish()

        # Get the size of the buffer in bytes
        buffer_size_bytes = self.V_Coul_buff.size

        # Calculate the shape based on the known dimensions
        # For example, if you know the first two dimensions:
        shape_xy = self.gsh.ns[0:2][::-1]  # First two dimensions in reverse order
        elements_xy = shape_xy[0] * shape_xy[1]
        shape_z = buffer_size_bytes // (4 * elements_xy)  # 4 bytes per float32

        # Create an array with the correct shape
        verify_shape = (*shape_xy, shape_z)
        verify_vcoul = np.empty(verify_shape, dtype=np.float32)

        # Now copy should work
        cl.enqueue_copy(self.queue, verify_vcoul, self.V_Coul_buff)
        self.queue.finish()

        print(f"Verification buffer shape: {verify_vcoul.shape}")
        print(f"Verification buffer range: [{verify_vcoul.min():.3f}, {verify_vcoul.max():.3f}]")

        ####### Plot three slices through the middle
        # plt.figure(figsize=(15, 5))
        # mid_x= 0
        # mid_y= 40
        # mid_z= 7
        
        # # XY plane
        # plt.subplot(131)
        # plt.imshow(V_after_slab[:, :, mid_z], origin='lower')
        # # plt.imshow(verify_vcoul[:, :, mid_z], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XY plane (z={mid_z})')

        # # XZ plane
        # plt.subplot(132)
        # plt.imshow(V_after_slab[:, mid_y, :], origin='lower')
        # # plt.imshow(verify_vcoul[:, mid_y, :], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'XZ plane (y={mid_y})')

        # # YZ plane
        # plt.subplot(133)
        # plt.imshow(V_after_slab[ mid_x, :,:], origin='lower')
        # # plt.imshow(verify_vcoul[ mid_x, :,:], origin='lower')
        # plt.colorbar(label='Potential')
        # plt.title(f'YZ plane (x={mid_x})')
        # plt.suptitle('Potential After Slab')
        # plt.tight_layout()
        # plt.show()
            



        if bCheckVin:
            sh=self.gsh.ns[::-1]
            Vin = np.empty(sh, dtype=np.float32)
            cl.enqueue_copy(self.queue, Vin, Vin_buff)
            print("Vin min,max ", Vin.min(), Vin.max() )

        # return self.slabPotential( Vin_buff, nz_slab, dipol=dipol, bDownload=bReturn, bTranspose=bTranspose )
        return V_after_slab



    def release_unused_buffs(self, keep_names=set()):
        """Release all OpenCL buffers except for those whose *base* names are in
        *keep_names* (a set of strings like {'atoms', 'REQs'}). This is useful
        to free GPU memory between logically separate operations (e.g. after
        finishing Coulomb grids but before building Morse grids)."""
        for attr, val in list(self.__dict__.items()):
            if not attr.endswith('_buff'):
                continue
            base = attr[:-5]  # strip '_buff'
            if base in keep_names:
                continue
            if val is not None:
                print(f"Releasing buffer {attr} to free GPU memory")
                try:
                    val.release()
                except Exception as e:
                    print(f"Warning: failed to release {attr}: {e}")
            setattr(self, attr, None)

    def make_MorseFF(self, atoms, REQs, nPBC=(4, 4, 0), dg=None, ng=None, lvec=None, g0=None, GFFParams=(0.1, 1.5, 0.0, 0.0), bTime=True, bReturn=True ):

        T00 = time.perf_counter()

        # Only create a new grid if one doesn't exist or if parameters are explicitly provided
        # if self.gsh is None or dg is not None or ng is not None or lvec is not None or g0 is not None:
        #     # Use existing grid parameters if not explicitly provided
        #     if dg is None and hasattr(self, 'gsh'):
        #         dg = self.gsh.dg
        #         print(f"Using existing grid step size: {dg}")
        #     if lvec is None and hasattr(self, 'gsh') and self.gsh.lvec is not None:
        #         lvec = self.gsh.lvec
        #         print(f"Using existing lattice vectors: {lvec}")
        #     if g0 is None and hasattr(self, 'gsh'):
        #         g0 = self.gsh.g0
        #         print(f"Using existing grid origin: {g0}")
            
        #     grid = GridShape(ns=ng, dg=dg, lvec=lvec, g0=g0)
        #     print(f"Creating grid with: ns={grid.ns}, dg={grid.dg}, lvec={lvec}")
        #     self.set_grid(grid)
        
        # atoms.apos = xyzq[:,:3]
        # print("New_Atoms:",atoms.shape)
        # print("xyzq[:,:3] =", atoms[:,:3])
        print("GridFF_cl::make_MorseFF() dg = ",dg)


        na   = len(atoms)
        nxyz = self.gcl.nxyz

        buff_names={'atoms','REQs','V_Paul','V_Lond'}
        self.try_make_buffs(buff_names, na, nxyz)

        # Free large buffers that are no longer needed to prevent GPU out-of-memory
        self.release_unused_buffs({"atoms", "REQs", "V_Paul", "V_Lond"})

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
            npbc = float((nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1))
            nops = float(na_cl) * float(nxyz) * npbc
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
